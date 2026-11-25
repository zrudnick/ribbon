#include "predictions.h"
#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>
#include <cmath>
#include "ribbon.h"
#include "coverage.h"

// Decides which of the predictions of a finished bundle are kept: merges the
// copies found in different cells (merge_cell_copies), drops the ones that are
// overlapped by better covered predictions (filter_overlapping_predictions) or
// too weakly covered (drop_weak_predictions), and the single exons that are
// not transcripts of their own (drop_unsupported_single_exons). It also takes
// cells off transcripts that a cell's own reads contradict
// (drop_cells_with_competing_junctions). With guides, the predictions that
// have a guide's intron chain are labeled with it first
// (label_predictions_with_guides): the labeled ones get a lower coverage floor
// and skip the overlap, weak-junction and UMI-start filters. out.cpp only
// prints the result.

// Import globals from main program:

extern int32_t min_transcript_len;      // minimum length for a transcript to be printed
extern uint32_t splice_site_error;      // maximum allowed error at splice sites
extern float min_read_coverage;         // minimum read coverage per bundle bp for a multi-exon transcript
extern float min_single_exon_coverage;  // minimum read coverage per bundle bp for a single-exon transcript
extern float isoform_fraction;           // minimum coverage of a prediction as a fraction of the best one of its gene
extern bool smartseq3;                   // --sq3: read 1 with a UMI starts at the 5' end of its molecule
extern bool require_umi_starts;          // --sq3 rule, off with --keep-paths-without-umi-starts
extern bool require_junction_evidence;   // off with --keep-weak-junction-paths

namespace ribbon {

  // Minimum pooled coverage of a reference (guide) prediction; novel ones need
  // min_read_coverage (-c). Most weak reference predictions are right, so
  // they get a lower floor (see drop_weak_predictions).
  const float min_reference_coverage = 1;

  // A novel single-exon prediction must stand out from its surroundings: more
  // cells must have a read starting inside it than read starts at
  // single_exon_min_fold times the density from single_exon_background_gap to
  // single_exon_background_reach bp beyond its ends explain, at a p-value
  // below single_exon_alpha (a Poisson-binomial test over the cells, see
  // single_exon_stands_out). The fold is an effect size: with many cells the
  // test alone keeps intronic hotspots that are only a few times denser than
  // the intron (ribbon_data/benchmark/small_n/se_effect/README.md); 50 is the
  // enrichment the earlier ratio rule asked for. The gap keeps most fragments
  // of the exon itself out of the background. Most single exons that fail
  // this are unspliced pre-mRNA inside introns (see
  // drop_unsupported_single_exons).
  const double single_exon_alpha = 0.001;
  const double single_exon_min_fold = 50;
  const int32_t single_exon_background_gap = 200;
  const int32_t single_exon_background_reach = 2000;
  // With --sq3 it also needs min_single_exon_umi_starts UMI molecules whose
  // read 1 starts within single_exon_start_window bp of its 5' end; with
  // --tso-filter, molecules marked as TSO strand invasion do not count
  // (ribbon_data/benchmark/umi/tss/README.md, lead 3).
  const float min_single_exon_umi_starts = 3;
  const int32_t single_exon_start_window = 25;

  // With --sq3 a novel multi-exon prediction needs min_umi_start_molecules
  // UMI molecules whose read 1 starts within umi_start_window bp of its 5'
  // end (see drop_predictions_without_umi_starts).
  const float min_umi_start_molecules = 1;
  const int32_t umi_start_window = 100;

  // With --sq3, on the mitochondrial chromosome (and on every chromosome with
  // --tso-filter), a transcript starts where at least start_peak_min_molecules
  // UMI molecules start within
  // start_peak_merge_len bp, coverage over start_peak_window bp after the pile
  // is at least start_peak_min_step times the coverage before it, and the pile
  // is at least start_peak_min_share of that coverage after (a new population
  // of molecules, not a few starts inside a transcript, such as TSO strand
  // invasion). On real chrMT this finds 13 of the 15 mRNA / rRNA starts;
  // genome-wide it costs about 0.35 precision points without --tso-filter and
  // almost nothing with it (ribbon_data/benchmark/missed_isoforms/README.md);
  // see find_start_peaks.
  const double start_peak_min_molecules = 10;
  const int32_t start_peak_merge_len = 5;
  const int32_t start_peak_window = 50;
  const double start_peak_min_step = 2;
  const double start_peak_min_share = 0.05;

  // A novel multi-exon prediction needs min_junction_evidence distinct
  // pieces of evidence, pooled over cells, for each of its junctions (see
  // drop_paths_with_weak_junctions).
  const int32_t min_junction_evidence = 2;

  // The 5' ends of the reads that are read 1 with a UMI (--sq3), on one
  // strand, sorted, with their weights (reads or molecules).
  struct ReadStarts {
    std::vector<std::pair<int32_t, float>> by_position;

    // Total weight of the starts within window bp of position.
    float count_near(int32_t position, int32_t window) const {
      auto lo = std::lower_bound(by_position.begin(), by_position.end(),
                                 std::make_pair(position - window, -1.0f));
      float total = 0;
      for (auto it = lo; it != by_position.end()
           && it->first <= position + window; ++it)
        total += it->second;
      return total;
    }
  };

  // The start positions of all reads of a bundle, sorted, with their cell and
  // weight (read_count: with --sq3 the reads of a UMI molecule weigh 1 in
  // total), and each cell's total weight in the bundle (its depth).
  struct CellReadStarts {
    struct Start {
      int32_t position;
      int32_t cell;
      float weight;
      bool operator<(const Start& o) const { return position < o.position; }
    };
    std::vector<Start> by_position;
    std::vector<double> cell_depth; // per cell
    double total_depth = 0; // all cells

    // Total weight of the starts in [lo, hi].
    double weight_in(int32_t lo, int32_t hi) const {
      auto it = std::lower_bound(by_position.begin(), by_position.end(),
                                 Start{lo, 0, 0});
      double total = 0;
      for (; it != by_position.end() && it->position <= hi; ++it)
        total += it->weight;
      return total;
    }

    // Number of cells with at least one start in [lo, hi].
    int32_t cells_in(int32_t lo, int32_t hi) const {
      std::vector<bool> seen(cell_depth.size(), false);
      int32_t cells = 0;
      auto it = std::lower_bound(by_position.begin(), by_position.end(),
                                 Start{lo, 0, 0});
      for (; it != by_position.end() && it->position <= hi; ++it)
        if (!seen[it->cell]) {
          seen[it->cell] = true;
          cells++;
        }
      return cells;
    }
  };

  struct CoverageRank { // predinfo
    int32_t p;
    float cov;
    float max_cov;
    CoverageRank(int32_t p = 0, float c = 0.0, float m = 0.0)
        : p(p), cov(c), max_cov(m) {}
  };

  // Comparator establishing the print/priority order among predictions:
  // leftmost start, then highest coverage, then fewest exons, then
  // lexicographically by exon coordinates.
  int32_t compare_prediction_by_position(const void *p1, const void *p2) {
    CPrediction *a = (CPrediction *)p1;
    CPrediction *b = (CPrediction *)p2;

    // Order by start, then by higher coverage.
    if (a->start < b->start) return -1; // order based on start
    if (a->start > b->start) return 1;

    if (a->cov > b->cov) return -1; // the one with higher coverage comes first
    if (a->cov < b->cov) return 1;

    // Order based on number of exons
    if (count_of(a->exons) < count_of(b->exons)) return -1;
    if (count_of(a->exons) > count_of(b->exons)) return 1;

    // Same number of exons
    int32_t i = 0;
    uint32_t a1 = 0;
    uint32_t b1 = 0;
    while (i < count_of(a->exons)) {
      if (a->exons[i].start < b->exons[i].start) {
        a1 = a->exons[i].start;
        b1 = b->exons[i].start;
        break;
      }
      if (a->exons[i].end < b->exons[i].end) {
        a1 = a->exons[i].end;
        b1 = b->exons[i].end;
        break;
      }
      i++;
    }
    if (a1) {
      if (a1 < b1) return -1; // the one with the first exon comes first
      if (a1 > b1) return 1;
    }
    return 0;
  }

  // Returns whether predictions n1 and n2 overlap and have the same intron
  // chain: the same exon count and the same internal exon boundaries (only
  // the first start and the last end may differ).
  bool have_same_intron_chain(PredictionList& predictions, int32_t n1, int32_t n2) {
    // Genes don't overlap
    if ((predictions[n1]->end < predictions[n2]->start) || (predictions[n2]->end < predictions[n1]->start))
      return (false);

    // Compare the internal exon boundaries.
    int32_t nex = count_of(predictions[n1]->exons);
    if (nex != count_of(predictions[n2]->exons)) return (false);
    for (int32_t i = 0; i < nex; i++) {
      if (i > 0 && (predictions[n1]->exons[i].start != predictions[n2]->exons[i].start))
        return (false);
      if (i < nex - 1 && (predictions[n1]->exons[i].end != predictions[n2]->exons[i].end))
        return (false);
    }

    return (true);
  }

  // Comparator for CoverageRank entries: sorts from highest to lowest coverage.
  int32_t compare_coverage_rank(const void *p1,
                                const void *p2) {
    CoverageRank *a = (CoverageRank *)p1;
    CoverageRank *b = (CoverageRank *)p2;
    if (a->cov < b->cov) return 1;
    if (a->cov > b->cov) return -1;
    return 0;
  }

  // Walks the exons of a from index a_idx and the exons of b from index b_idx,
  // both in order. Returns 2 if an exon of a overlaps an exon of b, otherwise
  // 1 (the exons of one run out before any overlap is found).
  int32_t find_exon_overlap(std::vector<GSeg>& a_exons,
                            std::vector<GSeg>& b_exons, int32_t a_idx,
                            int32_t b_idx) {
    for (; a_idx < count_of(a_exons); a_idx++) {
      // Advance the exons of b if before exon a_idx of a
      while (b_idx < count_of(b_exons)
             && b_exons[b_idx].end < a_exons[a_idx].start) b_idx++;
      if (b_idx == count_of(b_exons)) return (1); // no exon overlap found
      if (b_exons[b_idx].overlap(a_exons[a_idx])) return (2);
    }
    return (1);
  }

  // Classifies how predictions n1 and n2 overlap: 0 = none, 1 = locus only,
  // 2 = exon overlap, 3 = n1 included in n2, 4 = n2 included in n1. The
  // overhangs get the distances or how far an inclusion sticks out.
  int32_t classify_prediction_overlap(PredictionList& predictions, int32_t n1,
                                      int32_t n2, int32_t& left_overhang,
                                      int32_t& right_overhang) {
    // Genes don't overlap
    if ((predictions[n1]->end < predictions[n2]->start)
        || (predictions[n2]->end < predictions[n1]->start)) {
      // Distance from n2 to n1 -> positive if n1 is at the left of n2
      left_overhang = predictions[n2]->start - predictions[n1]->end;
      // Distance from n1 to n2 -> positive if n1 is at the right of n2
      right_overhang = predictions[n1]->start - predictions[n2]->end;
      return (0); // genes don't overlap
    }

    // Predictions are not equal as we checked this before -> check for overlap
    if (count_of(predictions[n1]->exons) == count_of(predictions[n2]->exons))
      return find_exon_overlap(predictions[n1]->exons, predictions[n2]->exons, 0, 0);

    // Different number of exons -> one prediction might be included in the
    // other. Find the prediction with more exons (big) and with fewer (small).
    int32_t small = n1;
    int32_t big = n2;
    if (count_of(predictions[n1]->exons) > count_of(predictions[n2]->exons)) {
      big = n1;
      small = n2;
    }
    std::vector<GSeg>& small_exons = predictions[small]->exons;
    std::vector<GSeg>& big_exons = predictions[big]->exons;

    // Check if small is included into big
    int32_t e1 = 0; // small
    int32_t e2 = 0; // big
    // Advance exons in predictions[big]
    while (e2 < count_of(big_exons)
           && big_exons[e2].end < small_exons[e1].start) e2++;

    // Not enough exons for big to contain small -> check for overlap instead
    if (count_of(big_exons) - e2 < count_of(small_exons))
      return find_exon_overlap(small_exons, big_exons, 0, e2);

    // big might contain small
    // Leftoverhang -> positive if small not included
    left_overhang = big_exons[e2].start - small_exons[e1].start;

    bool overlap = false;
    while (e1 < count_of(small_exons) - 1) {
      // Intron is different -> revert to overlap check
      if (small_exons[e1].end != big_exons[e2].end
          || small_exons[e1 + 1].start != big_exons[e2 + 1].start) {
        // No need to check anymore if there was overlap before
        if (overlap) return (2);
        return find_exon_overlap(small_exons, big_exons, e1, e2);
      }
      overlap = true;
      e1++;
      e2++;
    }

    // Check that e1 overlaps e2
    // small is included in big if we made it so far
    if (overlap || big_exons[e2].overlap(small_exons[e1])) {
      // Rightoverhang -> positive if small not included
      right_overhang = small_exons[e1].end - big_exons[e2].end;
      if (small == n1) return (3);
      else return (4);
    }

    return (1); // no exon overlap found
  }

  // For two multi-exon predictions on opposite strands, returns whether one's
  // first exon reaches past the other's last exon start. Single-exon
  // predictions never count.
  bool is_antisense_overlap(PredictionList& predictions, int32_t n1, int32_t n2) {
    // Only different strands can be considered non-overlapping
    if (predictions[n1]->strand == predictions[n2]->strand) return (false);
    // Consider only intronic genes for this
    if (count_of(predictions[n1]->exons) == 1 || count_of(predictions[n2]->exons) == 1)
      return (false);

    // Check if the first exon of one reaches past the last one of the other.
    if (predictions[n1]->exons[0].end > predictions[n2]->exons.back().start
        || predictions[n2]->exons[0].end > predictions[n1]->exons.back().start) return (true);

    return (false);
  }

  // Merges copies of the same prediction found in different cells: adds up
  // their coverages, keeps the coverage of each cell (cell_cov) and drops the
  // non-reference predictions that are too weak or too short.
  void merge_cell_copies(PredictionList& predictions, int32_t cell_count) {
    int32_t prediction_count = count_of(predictions);

    for (int32_t p = 0; p < prediction_count; p++)
      if (predictions[p] && predictions[p]->flag) {
        predictions[p]->cell_cov.resize(cell_count, 0.0);
        float best_coverage = predictions[p]->cov;
        // Coverage, start and end of the best copy so far (a reference copy's
        // start and end take precedence).
        uint32_t best_start = predictions[p]->start;
        uint32_t best_end = predictions[p]->end;

        // gene_no holds the number of cells the prediction is expressed in
        // until out.cpp renumbers it as the gene number.
        predictions[p]->gene_no = 1;
        // A guide prediction has cell -1 and arrives with the cell_cov of the
        // cells credited with it (see credit_guide_cell).
        int32_t guide_cells = 0;
        for (float guide_cell_coverage : predictions[p]->cell_cov) {
          if (guide_cell_coverage) guide_cells++;
        }
        if (guide_cells) predictions[p]->gene_no = guide_cells;
        // Predictions of the cell-agnostic fallback path (core.cpp) have cell
        // -1: they have no per-cell coverage to record.
        if (predictions[p]->cell >= 0)
          predictions[p]->cell_cov[predictions[p]->cell] += predictions[p]->cov;

        // Merge the following predictions that have the same intron chain.
        int32_t m = p + 1;
        // Removed predictions and predictions that start inside the first exon
        // of p are the ones to look at.
        auto is_removed_or_inside_first_exon = [&](int32_t idx) {
          return !predictions[idx]
                 || predictions[idx]->start <= predictions[p]->exons[0].end;
        };
        while (m < prediction_count && is_removed_or_inside_first_exon(m)) {
          if (predictions[m] && predictions[m]->flag) {
            if (predictions[p]->strand == predictions[m]->strand) {
              // Predictions n and m are equal
              if (have_same_intron_chain(predictions, p, m)) {
                if (predictions[p]->t_eq && predictions[m]->t_eq
                    && predictions[p]->t_eq != predictions[m]->t_eq) {
                  m++;
                  continue;
                } // both are equal but represent different transcripts

                if (predictions[m]->t_eq) { // reference prediction takes precedence
                  predictions[p]->t_eq = predictions[m]->t_eq;
                  best_start = predictions[m]->start;
                  best_end = predictions[m]->end;
                } else if (!predictions[p]->t_eq) {
                  // Use start/end of maximally covered prediction
                  if (predictions[m]->cov > best_coverage) {
                    best_coverage = predictions[m]->cov;
                    best_start = predictions[m]->start;
                    best_end = predictions[m]->end;
                  }
                }

                // See note above on predictions[n]->cell >= 0
                // A guide prediction m (cell -1) carries its cells' cell_cov.
                for (int32_t c = 0; c < count_of(predictions[m]->cell_cov); c++) {
                  if (!predictions[m]->cell_cov[c]) continue;
                  if (!predictions[p]->cell_cov[c]) predictions[p]->gene_no++;
                  predictions[p]->cell_cov[c] += predictions[m]->cell_cov[c];
                }
                if (predictions[m]->cell >= 0) {
                  if (!predictions[p]->cell_cov[predictions[m]->cell])
                    predictions[p]->gene_no++;
                  predictions[p]->cell_cov[predictions[m]->cell] += predictions[m]->cov;
                }
                // The merged raw ratio: the mean of both, weighted by coverage
                double merged_cov = predictions[p]->cov + predictions[m]->cov;
                if (merged_cov > 0)
                  predictions[p]->raw_ratio
                      = (predictions[p]->cov * predictions[p]->raw_ratio
                         + predictions[m]->cov * predictions[m]->raw_ratio)
                        / merged_cov;
                predictions[p]->cov += predictions[m]->cov;
                for (int32_t k = 0; k < count_of(predictions[p]->exons); k++) {
                  predictions[p]->exon_cov[k] += predictions[m]->exon_cov[k];
                }
                predictions[m].reset();
              }
            }
          }
          m++;
        }

        // Drop non-reference predictions that are weakly covered or too short
        // and extend the others to their best start and end.
        if (!predictions[p]->t_eq) {
          bool is_weakly_covered = predictions[p]->cov < min_read_coverage;
          bool is_weak_single_exon = count_of(predictions[p]->exons) == 1
              && predictions[p]->cov < min_single_exon_coverage;
          if (is_weakly_covered || is_weak_single_exon) {
            predictions[p]->flag = false;
            continue;
          }
          if (best_start != predictions[p]->start) {
            predictions[p]->tlen += predictions[p]->start - best_start;
            predictions[p]->start = best_start;
            predictions[p]->exons[0].start = best_start;
          }
          if (best_end != predictions[p]->end) {
            predictions[p]->tlen += best_end - predictions[p]->end;
            predictions[p]->end = best_end;
            predictions[p]->exons.back().end = best_end;
          }
          if (predictions[p]->tlen < min_transcript_len) {
            predictions[p]->flag = false;
            continue;
          }
        }
      }

    // Remove the predictions that were merged into others.
    predictions.erase(std::remove_if(predictions.begin(), predictions.end(),
                              [](const std::unique_ptr<CPrediction>& p) {
                                return !p;
                              }),
               predictions.end());
  }

  // Lists every prediction with its coverage and its best single-cell coverage,
  // sorted from the most to the least covered.
  std::vector<CoverageRank> build_coverage_order(PredictionList& predictions) {
    std::vector<CoverageRank> coverage_order;
    for (int32_t i = 0; i < count_of(predictions); i++) {
      // A prediction that was unflagged before the merge has no cell_cov
      float max_cell_coverage = predictions[i]->cell_cov.empty() ? 0
                                                         : predictions[i]->cell_cov[0];
      for (int32_t j = 1; j < count_of(predictions[i]->cell_cov); j++) {
        if (predictions[i]->cell_cov[j] > max_cell_coverage)
          max_cell_coverage = predictions[i]->cell_cov[j];
      }

      CoverageRank p(i, predictions[i]->cov, max_cell_coverage);
      coverage_order.push_back(p);
    }
    sort_by(coverage_order, compare_coverage_rank);
    return coverage_order;
  }

  // Goes through the predictions from the most covered to the least and clears
  // the flag of the ones that are included in or lose to a better covered
  // overlapping prediction.
  void filter_overlapping_predictions(PredictionList& predictions,
                                      const std::vector<CoverageRank>& coverage_order,
                                      int32_t cell_count) {
    for (int32_t i = 0; i < count_of(coverage_order) - 1; i++) {
      // If predictions[n1] was not eliminated before
      if (predictions[coverage_order[i].p]->flag) {
        int32_t n1 = coverage_order[i].p;
        for (int32_t j = i + 1; j < count_of(coverage_order); j++) {
          if (predictions[coverage_order[j].p]->flag) {
            // Always predord.cov(n1)>=predord.cov(n2)
            int32_t n2 = coverage_order[j].p;
            int32_t left_overhang = 0;
            int32_t right_overhang = 0;
            // Overlap code: 0 = no overlap; 1 = locus overlap; 2 = exon
            // overlap; 3 = n1 included in n2; 4 = n2 included in n1
            int32_t ovl_code = classify_prediction_overlap(predictions, n1, n2,
                left_overhang, right_overhang);
            if (ovl_code) { // genes overlap
              // Both single exon -> keep the most covered one -> they are on
              // different strands because otherwise they are considered equal
              if (count_of(predictions[n1]->exons) == 1
                  && count_of(predictions[n2]->exons) == 1) {
                predictions[n2]->flag = false;
                continue; // n2 is less covered here
              }
              // Only one prediction is single exon now
              bool has_single_exon_prediction
                  = count_of(predictions[n1]->exons) == 1
                    || count_of(predictions[n2]->exons) == 1;
              bool is_same_strand
                  = predictions[n1]->strand == predictions[n2]->strand;
              if (has_single_exon_prediction && is_same_strand) {
                int32_t big = -1;
                int32_t small = -1;
                // The predictions[n1] is a single exon included in
                // predictions[n2] -> no need to check because one has to be
                // single exon
                if (ovl_code == 3) {
                  big = n2;
                  small = n1;
                // The predictions[n2] is a single exon included in predictions[n1]
                } else if (ovl_code == 4) {
                  big = n1;
                  small = n2;
                }
                if (big > -1) { // this should be a polyA-seq feature only
                  int32_t e = 0;
                  int32_t overlap_len = 0;
                  if (predictions[big]->strand == '+') { // last exon
                    overlap_len
                        = predictions[big]->exons.back().overlapLen(predictions[small].get());
                    e = count_of(predictions[big]->exons) - 1;
                  } else if (predictions[big]->strand == '-') { // first exon only
                    overlap_len
                        = predictions[big]->exons[0].overlapLen(predictions[small].get());
                  }

                  if (overlap_len > 0.7 * predictions[small]->len()) { // v17

                    // Add cell coverages from predictions[small] to
                    // predictions[big]
                    predictions[small]->flag = false;

                    float added_coverage = predictions[small]->cov * overlap_len;
                    predictions[big]->cov += added_coverage / predictions[big]->tlen;
                    predictions[big]->exon_cov[e] += added_coverage
                        / predictions[big]->exons[e].len();
                    for (int32_t c = 0; c < cell_count; c++)
                      // If cell is non zero in predictions[small]
                      if (predictions[small]->cell_cov[c]) {
                        // Record cell coverage in predictions[big]
                        if (!predictions[big]->cell_cov[c]) predictions[big]->gene_no++;
                        predictions[big]->cell_cov[c] += predictions[small]->cell_cov[c]
                            * overlap_len / predictions[big]->tlen;
                      }

                    if (!predictions[n1]->flag) break; // either n1 or n2 are false
                    else continue;
                  }
                }
              }

              // Do not consider genes overlapping in this case
              if (ovl_code == 2 && is_antisense_overlap(predictions, n1, n2)) continue;

              // All genes overlap on the same strand now

              // A novel single exon inside a better covered prediction is
              // dropped; one with a guide's chain is never dropped here
              if (!predictions[n2]->t_eq
                  && count_of(predictions[n2]->exons) == 1) {
                predictions[n2]->flag = false;
                continue;
              }

              if (!predictions[n1]->t_eq && count_of(predictions[n1]->exons) == 1) {
                if (predictions[n1]->cov < predictions[n2]->cov + min_single_exon_coverage) {
                  predictions[n1]->flag = false;
                  break;
                }
              }

              // Don't allow overhangs that are too big
              if (abs(left_overhang) < (int32_t)splice_site_error
                  && abs(right_overhang) < (int32_t)splice_site_error) {
                // 3 = n1 included in n2;
                if (ovl_code == 3
                    && coverage_order[i].max_cov < coverage_order[j].max_cov) {
                  predictions[n1]->flag = false;
                  break;
                }

                //  4 = n2 included in n1
                if (ovl_code == 4
                    && coverage_order[i].max_cov < coverage_order[j].max_cov) {
                  predictions[n2]->flag = false;
                  continue;
                }
              }
            }
          }
        }
      }
    }
  }

  // Drops the predictions covered below min_read_coverage (novel ones) or
  // min_reference_coverage (reference ones).
  void drop_weak_predictions(PredictionList& predictions) {
    for (auto& p : predictions) {
      float min_coverage = p->t_eq ? min_reference_coverage : min_read_coverage;
      if (p->flag && p->cov < min_coverage) p->flag = false;
    }
  }

  // Drops the predictions that no cell is credited with: a transcript is only
  // reported with the cells that have it.
  void drop_predictions_without_cells(PredictionList& predictions) {
    for (auto& p : predictions) {
      if (!p->flag) continue;
      bool has_cell = false;
      for (float coverage : p->cell_cov)
        if (coverage > 0) has_cell = true;
      if (!has_cell) p->flag = false;
    }
  }

  // Returns the root of the gene that prediction idx belongs to (gene_parent is
  // the union-find parent array of the predictions), and points idx at the
  // root it walked to.
  int32_t find_gene_root(std::vector<int32_t>& gene_parent, int32_t idx) {
    int32_t root = gene_parent[idx];
    while (root != gene_parent[root]) {
      root = gene_parent[root];
      gene_parent[idx] = root;
    }
    return root;
  }

  // Groups the unflagged predictions that share an exon on the same strand
  // into genes (union-find). Returns, for each prediction, the index of the
  // root of its gene.
  std::vector<int32_t> group_predictions_into_genes(PredictionList& predictions) {
    int32_t prediction_count = count_of(predictions);
    std::vector<int32_t> gene_parent;
    for (int32_t i = 0; i < prediction_count; i++) gene_parent.push_back(i);

    for (int32_t i = 0; i < prediction_count - 1; i++)
      if (predictions[i]->flag) {
        int32_t root_i = find_gene_root(gene_parent, i);
        int32_t j = i + 1;
        while (j < prediction_count && predictions[i]->end >= predictions[j]->start) {
          if (predictions[j]->flag && predictions[i]->strand == predictions[j]->strand) {
            // Check overlap betweeg genes
            bool overlap = false;
            int32_t m = 0;
            for (int32_t k = 0; k < count_of(predictions[i]->exons); k++) {
              while (m < count_of(predictions[j]->exons)
                     && predictions[j]->exons[m].end < predictions[i]->exons[k].start) m++;
              if (m == count_of(predictions[j]->exons)) break;
              if (predictions[j]->exons[m].overlap(predictions[i]->exons[k])) {
                overlap = true;
                break;
              }
            }
            if (overlap) {
              int32_t root_j = find_gene_root(gene_parent, j);
              if (root_j < root_i) {
                gene_parent[root_i] = root_j;
                gene_parent[i] = root_j;
                // i now belongs to root_j: the next merge must link root_j, not
                // the old root_i, or it would overwrite this link
                root_i = root_j;
              } else if (root_i < root_j) {
                gene_parent[root_j] = root_i;
                gene_parent[j] = root_i;
              }
            }
          }
          j++;
        }
      }

    return gene_parent;
  }

  // Drops (flag = false) the predictions whose coverage, pooled over the cells,
  // is below isoform_fraction times that of the most covered prediction of
  // their gene (as group_predictions_into_genes defines it). This includes the
  // reference ones. The predictions are sorted by position.
  void drop_weak_isoforms(PredictionList& predictions) {
    std::vector<int32_t> gene_parent = group_predictions_into_genes(predictions);
    int32_t prediction_count = count_of(predictions);

    // The best coverage of each gene, by the index of its root
    std::vector<double> best_coverage(prediction_count, 0);
    for (int32_t p = 0; p < prediction_count; p++)
      if (predictions[p]->flag) {
        int32_t root = find_gene_root(gene_parent, p);
        best_coverage[root] = std::max(best_coverage[root], predictions[p]->cov);
      }

    for (int32_t p = 0; p < prediction_count; p++)
      if (predictions[p]->flag
          && predictions[p]->cov < isoform_fraction
                 * best_coverage[find_gene_root(gene_parent, p)])
        predictions[p]->flag = false;
  }

  // The junctions that have at least one read in one cell of a bundle, and the
  // donors (junction starts) and acceptors (junction ends) of those junctions.
  struct CellJunctionReads {
    std::unordered_set<uint64_t> junctions;
    std::unordered_set<uint32_t> donors;
    std::unordered_set<uint32_t> acceptors;
  };

  uint64_t junction_key(uint32_t start, uint32_t end) {
    return ((uint64_t)start << 32) | end;
  }

  // Lists, for each cell of the bundle, the junctions its reads are split
  // across.
  std::vector<CellJunctionReads> collect_cell_junction_reads(BundleData *bundle,
                                                             int32_t cell_count) {
    std::vector<CellJunctionReads> cell_reads(cell_count);
    for (auto& read : bundle->reads) {
      if (read->bc < 0 || read->bc >= cell_count) continue;
      for (CJunction *junction : read->juncs) {
        // The first junction of the bundle is a placeholder
        if (!junction || !junction->start) continue;
        CellJunctionReads& cell = cell_reads[read->bc];
        cell.junctions.insert(junction_key(junction->start, junction->end));
        cell.donors.insert(junction->start);
        cell.acceptors.insert(junction->end);
      }
    }
    return cell_reads;
  }

  // Returns whether, in a cell, one intron of the prediction has no read but
  // another junction with the same donor or the same acceptor has reads: the
  // cell then shows a different way of splicing at that site.
  bool has_competing_junction(CPrediction& prediction,
                              const CellJunctionReads& cell) {
    for (int32_t i = 0; i < count_of(prediction.exons) - 1; i++) {
      uint32_t donor = prediction.exons[i].end;
      uint32_t acceptor = prediction.exons[i + 1].start;
      if (cell.junctions.count(junction_key(donor, acceptor))) continue;
      // The intron has no read, so any junction found at its sites is another one
      if (cell.donors.count(donor) || cell.acceptors.count(acceptor)) return true;
    }
    return false;
  }

  // Takes a cell off each novel multi-exon prediction that has a competing
  // junction in that cell (see has_competing_junction). The coverage of the
  // cell leaves the prediction, which is dropped when no cell is left. This
  // holds for predictions with a guide's chain too: the cell's own reads
  // contradict them.
  void drop_cells_with_competing_junctions(BundleData *bundle,
                                           int32_t cell_count) {
    std::vector<CellJunctionReads> cell_reads
        = collect_cell_junction_reads(bundle, cell_count);

    for (auto& prediction : bundle->predictions) {
      if (!prediction->flag || count_of(prediction->exons) < 2) continue;

      float removed_coverage = 0;
      bool has_cell_left = false;
      for (int32_t c = 0; c < count_of(prediction->cell_cov); c++) {
        if (prediction->cell_cov[c] <= 0) continue;
        if (has_competing_junction(*prediction, cell_reads[c])) {
          removed_coverage += prediction->cell_cov[c];
          prediction->cell_cov[c] = 0;
          prediction->gene_no--;
        } else {
          has_cell_left = true;
        }
      }
      if (!removed_coverage) continue;

      if (!has_cell_left) {
        prediction->flag = false;
        continue;
      }

      // The exon coverages are pooled over the cells: scale them down with the
      // total coverage
      float old_coverage = prediction->cov;
      prediction->cov = std::max(0.0f, old_coverage - removed_coverage);
      for (float& exon_coverage : prediction->exon_cov)
        exon_coverage *= prediction->cov / old_coverage;
    }
  }

  // True if the genome right upstream of a UMI read's 5' end (with
  // --tso-filter) shows TSO strand invasion: GGG within the 4 bases before it
  // on the transcript's strand. The Smart-seq3 TSO ends in GGG; where it
  // anneals inside a cDNA instead of at the RNA's 5' end, the genome there
  // matches it (Tang et al. 2013; Hagemann-Jensen et al. 2022).
  bool is_tso_strand_invasion(BundleData *bundle, char strand,
                              int32_t five_prime) {
    if (!bundle->tso_gseq) return false;
    int32_t len = bundle->end - bundle->start + 1;
    // The 4 bases upstream on the strand: before the 5' end on '+', after it
    // on '-' (where GGG reads CCC on the forward strand)
    int32_t first = strand == '+' ? five_prime - 4 : five_prime + 1;
    int32_t from = first - bundle->start;
    if (from < 0 || from + 4 > len) return false;
    char g = strand == '+' ? 'G' : 'C';
    for (int32_t i = from; i + 3 <= from + 4; i++)
      if (bundle->tso_gseq[i] == g && bundle->tso_gseq[i + 1] == g
          && bundle->tso_gseq[i + 2] == g)
        return true;
    return false;
  }

  // The 5' ends of the bundle's UMI read 1s (--sq3), per strand. Weighted by
  // read count, or with by_molecule by molecule: each (cell, UMI) counts once
  // per 5' end, 1 / nh for a multi-mapped one. The strand is the one the UMI
  // gave the read (umi_strand). With skip_invasion, molecules marked as TSO
  // strand invasion are left out (is_tso_strand_invasion).
  PerStrand<ReadStarts> collect_umi_starts(BundleData *bundle,
                                           bool by_molecule,
                                           bool skip_invasion = false) {
    PerStrand<ReadStarts> umi_starts;
    // (5' end, cell, UMI) -> weight, for by_molecule
    PerStrand<std::map<std::tuple<int32_t, int32_t, std::string>, float>>
        molecules;
    for (const auto& read : bundle->reads) {
      if (read->ub.empty() || !read->is_read1) continue;
      char strand = read->umi_strand;
      if (strand != '+' && strand != '-') continue;
      int32_t five_prime = strand == '+' ? read->start : read->end;
      if (skip_invasion && is_tso_strand_invasion(bundle, strand, five_prime))
        continue;
      if (!by_molecule) {
        umi_starts[strand].by_position.push_back(
            {five_prime, read->read_count});
        continue;
      }
      float& weight = molecules[strand][{five_prime, read->bc, read->ub}];
      weight = std::max(weight, 1.0f / std::max<int16_t>(read->nh, 1));
    }

    // The molecules, then sort each strand by position
    for (char strand : {'-', '+'}) {
      for (const auto& [key, weight] : molecules[strand])
        umi_starts[strand].by_position.push_back({std::get<0>(key), weight});
      std::sort(umi_starts[strand].by_position.begin(),
                umi_starts[strand].by_position.end());
    }
    return umi_starts;
  }

  // Collects the read starts of the bundle with their cells and weights, and
  // each cell's depth (see CellReadStarts).
  CellReadStarts collect_cell_read_starts(BundleData *bundle,
                                          int32_t cell_count) {
    CellReadStarts starts;
    starts.cell_depth.assign(cell_count, 0);
    for (const auto& read : bundle->reads) {
      if (read->read_count <= 0 || read->bc < 0 || read->bc >= cell_count)
        continue;
      starts.by_position.push_back(
          {(int32_t)read->start, read->bc, read->read_count});
      starts.cell_depth[read->bc] += read->read_count;
      starts.total_depth += read->read_count;
    }
    std::sort(starts.by_position.begin(), starts.by_position.end());
    return starts;
  }

  // P(K >= k), where K is the number of successes in independent trials with
  // the given success probabilities (a Poisson-binomial distribution),
  // computed exactly by dynamic programming over the trials.
  double poisson_binomial_tail(const std::vector<double>& probabilities,
                               int32_t k) {
    if (k <= 0) return 1;

    // dist[j]: probability of exactly j successes so far (j < k); dist[k]:
    // of k or more, which stays k or more
    std::vector<double> dist(k + 1, 0.0);
    dist[0] = 1;
    for (double p : probabilities) {
      if (p <= 0) continue; // cannot succeed, changes nothing
      // From the top down, so each value is updated from the previous trial
      dist[k] += dist[k - 1] * p;
      for (int32_t j = k - 1; j > 0; j--)
        dist[j] = dist[j] * (1 - p) + dist[j - 1] * p;
      dist[0] *= 1 - p;
    }
    return dist[k];
  }

  // True if more cells have a read starting inside the single-exon prediction
  // p than its surroundings explain (p-value below single_exon_alpha). Under
  // the null, reads start inside at single_exon_min_fold times the density of
  // the two background windows beyond its ends, shared over the cells by their
  // depth in the bundle.
  bool single_exon_stands_out(const CellReadStarts& starts,
                              const CPrediction& p) {
    int32_t start = p.exons[0].start, end = p.exons[0].end;
    int32_t cells_inside = starts.cells_in(start, end);
    if (cells_inside == 0 || starts.total_depth <= 0) return false;

    // Density of read starts per bp in the two windows beyond the ends (no
    // reads outside the bundle)
    int32_t gap = single_exon_background_gap;
    int32_t reach = single_exon_background_reach;
    double around = starts.weight_in(start - reach, start - gap)
        + starts.weight_in(end + gap, end + reach);
    double density = around / (2 * (reach - gap + 1));

    // Each cell's chance of a read starting inside under the null: Poisson
    // with mean fold x density x length x the cell's share of the bundle's
    // reads
    double length = end - start + 1;
    double null_density = single_exon_min_fold * density;
    std::vector<double> chance;
    chance.reserve(starts.cell_depth.size());
    for (double depth : starts.cell_depth)
      chance.push_back(
          1 - std::exp(-null_density * length * depth / starts.total_depth));
    return poisson_binomial_tail(chance, cells_inside) < single_exon_alpha;
  }

  // Drops the single-exon predictions that do not look like transcripts of
  // their own: those seen in no more cells than the reads around them explain
  // and, with --sq3, those with too few UMI molecules starting at their 5'
  // end. Single exons with a guide are tested too: many annotated single exons
  // are pseudogenes or retained introns (ribbon_data/benchmark/guides). On the
  // real data most dropped ones are unspliced pre-mRNA inside introns
  // (ribbon_data/benchmark/README.md).
  void drop_unsupported_single_exons(BundleData *bundle) {
    PredictionList& predictions = bundle->predictions;
    bool any = false;
    for (const auto& p : predictions)
      if (p->flag && count_of(p->exons) == 1) any = true;
    if (!any) return;

    // With --sq3: the 5' ends of the UMI read 1s of each strand
    PerStrand<ReadStarts> umi_starts;
    if (smartseq3) umi_starts = collect_umi_starts(bundle, true, true);

    // The read starts of every cell, for the background test
    CellReadStarts starts
        = collect_cell_read_starts(bundle, count_of(bundle->cell_name));

    // The background test is against unspliced pre-mRNA in introns; the
    // mitochondrial genome has none, and its genes lie next to each other, so
    // there it would compare a gene with its neighbours
    bool mitochondrial = is_mitochondrial(bundle->ref_seq);

    for (auto& p : predictions) {
      if (!p->flag || count_of(p->exons) != 1) continue;
      if (!mitochondrial && !single_exon_stands_out(starts, *p)) {
        p->flag = false;
        continue;
      }
      if (smartseq3 && (p->strand == '+' || p->strand == '-')) {
        int32_t five_prime = p->strand == '+' ? p->exons[0].start
                                              : p->exons[0].end;
        if (umi_starts[p->strand].count_near(five_prime,
                single_exon_start_window) < min_single_exon_umi_starts)
          p->flag = false;
      }
    }
  }

  // With --sq3, drops the novel multi-exon predictions with too few UMI
  // molecules starting at their 5' end (min_umi_start_molecules). On the real
  // data most of them are fragments that start inside a later exon or far from
  // any annotated start (ribbon_data/benchmark/umi/tss/README.md).
  void drop_predictions_without_umi_starts(BundleData *bundle) {
    if (!smartseq3 || !require_umi_starts) return;
    PredictionList& predictions = bundle->predictions;
    bool any = false;
    for (const auto& p : predictions)
      if (p->flag && !p->t_eq && count_of(p->exons) > 1) any = true;
    if (!any) return;

    PerStrand<ReadStarts> umi_starts = collect_umi_starts(bundle, true);
    for (auto& p : predictions) {
      if (!p->flag || p->t_eq || count_of(p->exons) < 2) continue;
      if (p->strand != '+' && p->strand != '-') continue;
      int32_t five_prime = p->strand == '+' ? p->exons.front().start
                                            : p->exons.back().end;
      if (umi_starts[p->strand].count_near(five_prime, umi_start_window)
          < min_umi_start_molecules)
        p->flag = false;
    }
  }

  // Drops the novel multi-exon predictions with a junction that has fewer
  // than min_junction_evidence distinct pieces of evidence over all cells: a
  // distinct read alignment (cell, start, end; identical copies are already
  // one read), and with --sq3 a UMI molecule (cell, UMI, strand) instead of
  // its reads. Most junctions seen once are noise on the real data
  // (ribbon_data/benchmark/umi/README.md, lead 4). Reference ones are kept.
  void drop_paths_with_weak_junctions(BundleData *bundle) {
    if (!require_junction_evidence) return;
    PredictionList& predictions = bundle->predictions;

    // The junctions (exon end, next exon start) of the novel predictions
    using Junction = std::pair<int32_t, int32_t>;
    using Evidence = std::tuple<int32_t, int32_t, int32_t, std::string>;
    std::map<Junction, std::set<Evidence>> evidence;
    for (const auto& p : predictions) {
      if (!p->flag || p->t_eq || count_of(p->exons) < 2) continue;
      for (int32_t i = 0; i + 1 < count_of(p->exons); i++)
        evidence[{p->exons[i].end, p->exons[i + 1].start}];
    }
    if (evidence.empty()) return;

    // The evidence of each read for those junctions (the gaps between its
    // aligned blocks)
    for (const auto& read : bundle->reads) {
      if (!read->nh || count_of(read->segs) < 2) continue;
      bool is_molecule = smartseq3 && !read->ub.empty();
      Evidence piece = is_molecule
          ? Evidence(read->bc, -1, read->umi_strand, read->ub)
          : Evidence(read->bc, read->start, read->end, "");
      for (int32_t i = 0; i + 1 < count_of(read->segs); i++) {
        auto it = evidence.find({read->segs[i].end, read->segs[i + 1].start});
        if (it != evidence.end()) it->second.insert(piece);
      }
    }

    // Drop the predictions with a weak junction
    for (auto& p : predictions) {
      if (!p->flag || p->t_eq || count_of(p->exons) < 2) continue;
      for (int32_t i = 0; i + 1 < count_of(p->exons); i++)
        if (count_of(evidence[{p->exons[i].end, p->exons[i + 1].start}])
            < min_junction_evidence) {
          p->flag = false;
          break;
        }
    }
  }

  // A UMI molecule of the bundle (--sq3): its cell, the strand its UMI gives,
  // the 5' end of its read 1, its weight (1 / nh) and its reads (read 1, read
  // 2 and the PCR copies).
  struct UmiMolecule {
    int32_t cell = -1;
    char strand = '.';
    int32_t five_prime = -1;
    float weight = 0;
    std::vector<const CReadAln *> reads;
  };

  // The UMI molecules of the bundle, one per (cell, UMI, strand), as in
  // weight_reads_by_molecule. A molecule's 5' end is the most common 5' end
  // of its read 1s (copies usually share it); one without a read 1 here has
  // five_prime -1.
  std::vector<UmiMolecule> collect_umi_molecules(BundleData *bundle) {
    std::map<std::tuple<int32_t, std::string, char>, UmiMolecule> by_key;
    std::map<std::tuple<int32_t, std::string, char>,
             std::map<int32_t, float>> read1_ends;
    for (const auto& read : bundle->reads) {
      if (read->ub.empty()) continue;
      if (read->umi_strand != '+' && read->umi_strand != '-') continue;
      auto key = std::make_tuple(read->bc, read->ub, read->umi_strand);
      UmiMolecule& molecule = by_key[key];
      molecule.cell = read->bc;
      molecule.strand = read->umi_strand;
      molecule.weight = std::max(molecule.weight,
                                 1.0f / std::max<int16_t>(read->nh, 1));
      molecule.reads.push_back(read.get());
      if (read->is_read1) {
        int32_t five_prime = read->umi_strand == '+' ? read->start : read->end;
        read1_ends[key][five_prime] += read->read_count;
      }
    }

    // The most common read 1 end of each molecule is its 5' end
    std::vector<UmiMolecule> molecules;
    for (auto& [key, molecule] : by_key) {
      float best = 0;
      for (const auto& [five_prime, count] : read1_ends[key])
        if (count > best) {
          best = count;
          molecule.five_prime = five_prime;
        }
      molecules.push_back(std::move(molecule));
    }
    return molecules;
  }

  // True if the read lies on the exons of p: each aligned block inside one
  // exon, and each gap between blocks exactly an intron of p (between two
  // consecutive exons). The read's first block may start before p's first
  // exon and its last block may end after p's last exon (the 5' read can
  // start a little upstream of the transcript).
  bool read_fits_prediction(const CReadAln& read, const CPrediction& p) {
    int32_t exon_count = count_of(p.exons);
    int32_t block_count = count_of(read.segs);
    int32_t prev_exon = -1;
    for (int32_t b = 0; b < block_count; b++) {
      const GSeg& block = read.segs[b];
      // The exon that the block overlaps
      int32_t e = 0;
      while (e < exon_count && p.exons[e].end < block.start) e++;
      if (e == exon_count || p.exons[e].start > block.end) return false;

      // Inside the exon, except outward at the transcript's two ends
      bool open_left = b == 0 && e == 0;
      bool open_right = b == block_count - 1 && e == exon_count - 1;
      if (block.start < p.exons[e].start && !open_left) return false;
      if (block.end > p.exons[e].end && !open_right) return false;

      // The gap before the block is the intron before exon e
      if (b > 0) {
        if (e != prev_exon + 1) return false;
        if (read.segs[b - 1].end != p.exons[e - 1].end
            || block.start != p.exons[e].start)
          return false;
      }
      prev_exon = e;
    }
    return block_count > 0;
  }

  // The window where a molecule's 5' end may start for a prediction: from
  // umi_start_window bp upstream of its 5' end to the end of its first exon
  // (the start of its last exon on '-'), so that a molecule counts even where
  // the prediction's 5' end is placed upstream of the true start.
  struct StartWindow {
    int32_t lo;
    int32_t hi;
    int32_t p;
    bool operator<(const StartWindow& o) const { return lo < o.lo; }
  };

  // With --quant: in each cell and graph (QuantGroup, see save_quant_paths),
  // re-splits the coverage that the flow gave the cell's reported transcripts
  // at the nodes they share, in proportion to how strongly each is expressed:
  // its own reads (read transfrags that fit only it among the reported ones)
  // per position where an own read can start (its unique nodes' length plus one
  // fragment per run of unique nodes or edges). The greedy flow gives shared
  // nodes mostly to the first path. Paths that the filters dropped keep their
  // share (it is not reported). Sets cell_quant_cov; assembly is not changed.
  void resplit_quant_coverage(BundleData *bundle) {
    if (bundle->quant_groups.empty()) return;
    PredictionList& predictions = bundle->predictions;
    int32_t cell_count = count_of(bundle->cell_name);

    // The kept spliced predictions by strand and intron chain
    std::map<std::pair<char, std::vector<int32_t>>, CPrediction *> by_chain;
    for (auto& p : predictions) {
      if (!p->flag || count_of(p->exons) < 2) continue;
      std::vector<int32_t> chain;
      for (int32_t i = 0; i + 1 < count_of(p->exons); i++) {
        chain.push_back(p->exons[i].end);
        chain.push_back(p->exons[i + 1].start);
      }
      by_chain[{p->strand, chain}] = p.get();
    }
    // The kept prediction a path became in the cell (NULL if none)
    auto reported_as = [&](const QuantPath& path,
                           int32_t cell) -> CPrediction * {
      CPrediction *found = NULL;
      if (!path.intron_chain.empty()) {
        auto it = by_chain.find({path.strand, path.intron_chain});
        if (it != by_chain.end()) found = it->second;
      } else {
        for (auto& p : predictions)
          if (p->flag && count_of(p->exons) == 1 && p->strand == path.strand
              && p->start <= path.end && path.start <= p->end) {
            found = p.get();
            break;
          }
      }
      if (found && cell < count_of(found->cell_cov)
          && found->cell_cov[cell] > 0)
        return found;
      return NULL;
    };

    for (const QuantGroup& group : bundle->quant_groups) {
      int32_t cell = group.cell;
      std::vector<int32_t> kept; // the group's paths that were reported
      std::vector<CPrediction *> prediction_of;
      for (int32_t k = 0; k < count_of(group.paths); k++) {
        CPrediction *p = reported_as(group.paths[k], cell);
        prediction_of.push_back(p);
        if (p) kept.push_back(k);
      }
      if (count_of(kept) < 2) continue; // nothing to share
      std::vector<bool> is_kept(group.paths.size(), false);
      for (int32_t k : kept) is_kept[k] = true;

      // Own reads among the reported paths
      std::vector<double> own(group.paths.size(), 0);
      for (const auto& [abundance, fitting] : group.transfrags) {
        int32_t only = -1, n = 0;
        for (int32_t k : fitting)
          if (is_kept[k]) { only = k; n++; }
        if (n == 1) own[only] += abundance;
      }

      // The reported paths through each node and each edge
      std::map<int32_t, std::vector<int32_t>> paths_of_node;
      std::map<std::pair<int32_t, int32_t>, int32_t> paths_on_edge;
      for (int32_t k : kept) {
        const std::vector<int32_t>& nodes = group.paths[k].nodes;
        for (int32_t j = 0; j < count_of(nodes); j++) {
          if (nodes[j] > 0 && nodes[j] < group.sink)
            paths_of_node[nodes[j]].push_back(k);
          if (j) paths_on_edge[{nodes[j - 1], nodes[j]}]++;
        }
      }

      // Own reads per position where an own read can start
      std::vector<double> density(group.paths.size(), 0);
      for (int32_t k : kept) {
        const std::vector<int32_t>& nodes = group.paths[k].nodes;
        double positions = 0;
        bool in_run = false;
        for (int32_t j = 0; j < count_of(nodes); j++) {
          int32_t node = nodes[j];
          bool unique_node = node > 0 && node < group.sink
              && count_of(paths_of_node[node]) == 1;
          bool unique_edge = j && nodes[j - 1] > 0 && node < group.sink
              && paths_on_edge[{nodes[j - 1], node}] == 1;
          if (unique_node) positions += group.node_len.at(node);
          if (unique_node || unique_edge) {
            if (!in_run) positions += group.fragment_len;
            in_run = true;
          } else in_run = false;
        }
        if (positions > 0) density[k] = own[k] / positions;
      }

      // Split what the flow gave the reported paths at each node by density;
      // where one of them has no own reads, keep the flow's split there
      std::vector<double> bases(group.paths.size(), 0);
      std::vector<double> length(group.paths.size(), 0);
      for (const auto& [node, ks] : paths_of_node) {
        double density_total = 0, given_out = 0;
        bool all_have_own = true;
        std::vector<double> taken;
        for (int32_t k : ks) {
          if (density[k] <= 0) all_have_own = false;
          density_total += density[k];
          const QuantPath& path = group.paths[k];
          int32_t j = std::find(path.nodes.begin(), path.nodes.end(), node)
              - path.nodes.begin();
          taken.push_back(path.taken[j]);
          given_out += taken.back();
        }
        for (int32_t i = 0; i < count_of(ks); i++) {
          bases[ks[i]] += all_have_own && density_total > 0
              ? given_out * density[ks[i]] / density_total : taken[i];
          length[ks[i]] += group.node_len.at(node);
        }
      }

      // Each path's coverage; copies of one transcript in the cell (paths with
      // the same intron chain) add up, as in its cell_cov
      std::map<CPrediction *, double> quant_cov;
      for (int32_t k : kept)
        if (length[k] > 0) quant_cov[prediction_of[k]] += bases[k] / length[k];
      for (auto& [p, coverage] : quant_cov) {
        p->cell_quant_cov.resize(cell_count, 0.0);
        p->cell_quant_cov[cell] = coverage;
      }
    }
  }

  // A prediction's coverage in a cell for --quant: the re-split coverage
  // (cell_quant_cov) if it has one there, else its cell coverage.
  double quant_coverage(const CPrediction& prediction, int32_t cell) {
    if (cell < count_of(prediction.cell_quant_cov)
        && prediction.cell_quant_cov[cell] > 0)
      return prediction.cell_quant_cov[cell];
    return prediction.cell_cov[cell];
  }

  // With --sq3, counts for each kept prediction and cell the UMI molecules
  // whose 5' end lies in its StartWindow and whose reads all fit its exons
  // (cell_molecules). A molecule that fits several predictions of its cell is
  // split in proportion to their coverage in that cell; one that fits none is
  // not counted. With --tso-filter, molecules marked as TSO strand invasion
  // are not counted.
  void count_molecules_of_predictions(BundleData *bundle) {
    if (!smartseq3) return;
    PredictionList& predictions = bundle->predictions;
    int32_t cell_count = count_of(bundle->cell_name);

    // The windows of the kept predictions of each strand, sorted by their
    // start, and the longest window (to bound the search)
    PerStrand<std::vector<StartWindow>> windows;
    PerStrand<int32_t> longest(0);
    for (int32_t p = 0; p < count_of(predictions); p++) {
      CPrediction& prediction = *predictions[p];
      if (!prediction.flag || prediction.exons.empty()) continue;
      char strand = prediction.strand;
      if (strand != '+' && strand != '-') continue;
      StartWindow window;
      if (strand == '+') {
        window.lo = prediction.exons.front().start - umi_start_window;
        window.hi = prediction.exons.front().end;
      } else {
        window.lo = prediction.exons.back().start;
        window.hi = prediction.exons.back().end + umi_start_window;
      }
      window.p = p;
      windows[strand].push_back(window);
      longest[strand] = std::max(longest[strand], window.hi - window.lo);
      prediction.cell_molecules.assign(cell_count, 0);
    }
    for (char strand : {'-', '+'})
      std::sort(windows[strand].begin(), windows[strand].end());

    for (const UmiMolecule& molecule : collect_umi_molecules(bundle)) {
      if (molecule.five_prime < 0) continue;
      if (is_tso_strand_invasion(bundle, molecule.strand, molecule.five_prime))
        continue;
      const auto& strand_windows = windows[molecule.strand];
      StartWindow first = {molecule.five_prime - longest[molecule.strand], 0,
                           -1};
      auto it = std::lower_bound(strand_windows.begin(), strand_windows.end(),
                                 first);

      // The predictions of the molecule's cell that it fits
      std::vector<int32_t> fitting;
      double fitting_cov = 0;
      for (; it != strand_windows.end() && it->lo <= molecule.five_prime;
           ++it) {
        if (it->hi < molecule.five_prime) continue;
        CPrediction& prediction = *predictions[it->p];
        if (molecule.cell >= count_of(prediction.cell_cov)
            || prediction.cell_cov[molecule.cell] <= 0)
          continue;
        bool fits = true;
        for (const CReadAln *read : molecule.reads)
          if (!read_fits_prediction(*read, prediction)) {
            fits = false;
            break;
          }
        if (!fits) continue;
        fitting.push_back(it->p);
        fitting_cov += quant_coverage(prediction, molecule.cell);
      }

      // Split the molecule by the predictions' coverage in the cell (--quant's
      // coverage, see resplit_quant_coverage)
      for (int32_t p : fitting)
        predictions[p]->cell_molecules[molecule.cell] += molecule.weight
            * quant_coverage(*predictions[p], molecule.cell) / fitting_cov;
    }
  }

  // True if the prediction has the intron chain of guide (same strand and
  // number of exons); a single exon instead needs a single-exon guide that
  // overlaps at least half of each.
  bool has_guide_chain(const CPrediction& prediction, GffObj *guide) {
    int32_t exon_count = count_of(prediction.exons);
    if (guide->strand != prediction.strand
        || guide->exons.Count() != exon_count)
      return false;
    if (exon_count == 1) {
      int32_t overlap = std::min<int32_t>(guide->end, prediction.end)
          - std::max<int32_t>(guide->start, prediction.start) + 1;
      return 2 * overlap >= (int32_t)(guide->end - guide->start + 1)
             && 2 * overlap >= (int32_t)(prediction.end - prediction.start + 1);
    }
    for (int32_t e = 0; e < exon_count; e++) {
      if (e > 0 && prediction.exons[e].start != guide->exons[e]->start)
        return false;
      if (e + 1 < exon_count && prediction.exons[e].end != guide->exons[e]->end)
        return false;
    }
    return true;
  }

  // Labels each prediction with the first guide (t_eq) that has its intron
  // chain (has_guide_chain). The reads built the prediction; the guide only
  // names it and lets it skip some filters (see finalize_predictions).
  void label_predictions_with_guides(BundleData *bundle) {
    std::vector<GffObj *>& guides = bundle->guides;
    if (guides.empty()) return;
    for (auto& p : bundle->predictions) {
      if (!p || !p->flag || p->t_eq) continue;
      // The guides are sorted by start
      for (GffObj *guide : guides) {
        if (guide->start > p->end) break;
        if (guide->end < p->start) continue;
        if (has_guide_chain(*p, guide)) {
          p->t_eq = guide;
          break;
        }
      }
    }
  }

  // Gives each kept prediction with a guide the guide's start and end: the
  // reads rarely cover a transcript's ends fully.
  void snap_ends_to_guides(PredictionList& predictions) {
    for (auto& p : predictions) {
      if (!p || !p->flag || !p->t_eq
          || count_of(p->exons) != p->t_eq->exons.Count())
        continue;
      p->tlen += (int32_t)p->start - (int32_t)p->t_eq->start
          + (int32_t)p->t_eq->end - (int32_t)p->end;
      p->start = p->exons.front().start = p->t_eq->start;
      p->end = p->exons.back().end = p->t_eq->end;
    }
  }

  // Decides which predictions of the bundle are kept: sorts them by position,
  // merges the copies found in different cells, and clears the flag of the
  // ones that are dropped. The flagged predictions remain sorted by position
  // and carry their per-cell coverage (cell_cov).
  void finalize_predictions(BundleData *bundle) {
    PredictionList& predictions = bundle->predictions;
    sort_by(predictions, compare_prediction_by_position);

    // Name the predictions that have a guide's intron chain, so that the
    // filters below can treat them as known transcripts.
    label_predictions_with_guides(bundle);

    // Cells are numbered in the order they were first seen.
    int32_t cell_count = count_of(bundle->cell_name);
    merge_cell_copies(predictions, cell_count);

    // A cell's own reads can contradict a transcript the pooled graph gave it.
    drop_cells_with_competing_junctions(bundle, cell_count);

    // Every kept transcript has at least one cell.
    drop_predictions_without_cells(predictions);

    // Filter overlapping predictions, from the most covered to the least.
    std::vector<CoverageRank> coverage_order = build_coverage_order(predictions);
    filter_overlapping_predictions(predictions, coverage_order, cell_count);

    drop_weak_predictions(predictions);
    drop_weak_isoforms(predictions);

    // Single exons need coverage of their own (and 5' reads with --sq3).
    drop_unsupported_single_exons(bundle);

    // With --sq3, novel spliced transcripts need UMI molecules at their 5' end.
    drop_predictions_without_umi_starts(bundle);

    // Novel spliced transcripts need distinct evidence for each junction.
    drop_paths_with_weak_junctions(bundle);

    // The kept transcripts with a guide take its ends.
    snap_ends_to_guides(predictions);

    // With --quant, shared coverage re-split between each cell's transcripts.
    resplit_quant_coverage(bundle);

    // With --sq3, the UMI molecules of each kept transcript and cell.
    count_molecules_of_predictions(bundle);
  }

  // With --sq3, on the mitochondrial chromosome, and on every chromosome with
  // --tso-filter (bundle->tso_gseq; without it strand invasion piles pass):
  // fills bundle->start_peaks with the positions where a transcript starts by
  // its 5' UMI molecules (see start_peak_min_molecules). Mitochondrial mRNAs
  // are cut out of one long precursor, so coverage alone rarely separates them.
  // The source (or sink) link of a peak gets the coverage jump / DROP, like a
  // coverage trim.
  void find_start_peaks(BundleData *bundle) {
    if (!smartseq3
        || (!is_mitochondrial(bundle->ref_seq) && !bundle->tso_gseq))
      return;
    PerStrand<ReadStarts> umi_starts = collect_umi_starts(bundle, true, true);

    // Mean coverage per bp of a genomic window (no reads outside the bundle)
    auto mean_coverage = [&](int32_t lo, int32_t hi) {
      lo = std::max(lo, (int32_t)bundle->start);
      hi = std::min(hi, (int32_t)bundle->end);
      if (hi < lo) return 0.0;
      return get_region_coverage('.', lo - bundle->start, hi - bundle->start,
                                 bundle->bp_coverage)
             / start_peak_window;
    };

    for (char strand : {'-', '+'}) {
      // The molecules' 5' ends, sorted, with their weights
      const auto& starts = umi_starts[strand].by_position;
      size_t i = 0;
      while (i < starts.size()) {
        // A pile: 5' ends at most start_peak_merge_len apart; its peak is the
        // position with the most molecules
        double molecules = 0, peak_molecules = -1;
        int32_t peak = starts[i].first, last = starts[i].first;
        while (i < starts.size()
               && starts[i].first - last <= start_peak_merge_len) {
          int32_t position = starts[i].first;
          double here = 0;
          while (i < starts.size() && starts[i].first == position)
            here += starts[i++].second;
          molecules += here;
          if (here > peak_molecules) {
            peak_molecules = here;
            peak = position;
          }
          last = position;
        }
        if (molecules < start_peak_min_molecules) continue;

        // The coverage step across the peak, in the transcript's direction
        double before, after;
        if (strand == '+') {
          before = mean_coverage(peak - start_peak_window, peak - 1);
          after = mean_coverage(peak, peak + start_peak_window - 1);
        } else {
          before = mean_coverage(peak + 1, peak + start_peak_window);
          after = mean_coverage(peak - start_peak_window + 1, peak);
        }
        if (after < start_peak_min_step * before
            || molecules < start_peak_min_share * after)
          continue;
        bundle->start_peaks[strand].push_back(
            {(uint32_t)peak, (float)((after - before) / DROP)});
      }
    }
  }

} // namespace ribbon
