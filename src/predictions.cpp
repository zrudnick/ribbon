#include "predictions.h"
#include <algorithm>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>
#include "ribbon.h"
#include "coverage.h"

// Decides which of the predictions of a finished bundle are kept: merges the
// copies found in different cells (merge_cell_copies), drops the ones that are
// overlapped by better covered predictions (filter_overlapping_predictions) or
// too weakly covered (drop_weak_predictions), and the single exons that are
// not transcripts of their own (drop_unsupported_single_exons). It also takes
// cells off transcripts that a cell's own reads contradict
// (drop_cells_with_competing_junctions). out.cpp only prints the result.

// Import globals from main program:

extern int32_t min_transcript_len;      // minimum length for a transcript to be printed
extern uint32_t splice_site_error;      // maximum allowed error at splice sites
extern float min_read_coverage;         // minimum read coverage per bundle bp for a multi-exon transcript
extern float min_single_exon_coverage;  // minimum read coverage per bundle bp for a single-exon transcript
extern float isoform_fraction;           // minimum coverage of a prediction as a fraction of the best one of its gene
extern bool smartseq3;                   // --sq3: read 1 with a UMI starts at the 5' end of its molecule

namespace ribbon {

  // Minimum pooled coverage of a reference (guide) prediction; novel ones need
  // min_read_coverage (-c). Most weak reference predictions are right, so
  // they get a lower floor (see drop_weak_predictions).
  const float min_reference_coverage = 1;

  // A novel single-exon prediction must stand out from the coverage around
  // it: the mean coverage from single_exon_background_gap to
  // single_exon_background_reach bp beyond its ends may be at most
  // max_single_exon_background times its own. Most single exons that fail
  // this are unspliced pre-mRNA inside introns (see
  // drop_unsupported_single_exons).
  const float max_single_exon_background = 0.02;
  const int32_t single_exon_background_gap = 200;
  const int32_t single_exon_background_reach = 2000;
  // With --sq3 it also needs min_single_exon_umi_starts UMI reads (read 1)
  // that start within single_exon_start_window bp of its 5' end.
  const float min_single_exon_umi_starts = 3;
  const int32_t single_exon_start_window = 50;

  // The 5' ends of the reads that are read 1 with a UMI (--sq3), on one
  // strand, sorted, with their read counts.
  struct ReadStarts {
    std::vector<std::pair<int32_t, float>> by_position;

    // Total read count of the starts within window bp of position.
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
    for (int32_t i = 0; i < count_of(coverage_order) - 1; i++)
      // If predictions[n1] was not eliminated before
      if (predictions[coverage_order[i].p]->flag) {
        int32_t n1 = coverage_order[i].p;
        for (int32_t j = i + 1; j < count_of(coverage_order); j++)
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

              // This is not a known gene -> only then we can eliminate it
              if (!predictions[n2]->t_eq) {
                // More strict about novel predictions if annotation is
                // available
                if (predictions[n1]->t_eq
                    && predictions[n2]->cov < ERROR_PERC * predictions[n1]->cov) {
                  predictions[n2]->flag = false;
                  continue;
                } else if (count_of(predictions[n2]->exons) == 1) {
                  predictions[n2]->flag = false;
                  continue;
                }
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
        if (predictions[n1]->flag && count_of(predictions[n1]->exons) == 1) { // v22
          // A single-exon prediction survives only if it is seen in enough
          // cells (SCELL_NO) and well covered in at least one
          if (predictions[n1]->gene_no < SCELL_NO || coverage_order[i].max_cov < 1.3)
            predictions[n1]->flag = false;
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

  // Gives each guide that no cell was credited with (de novo found it in no
  // cell) the one cell where the flow gave it the most coverage, so that the
  // guide keeps a cell.
  void credit_best_cell_to_guides_without_cells(PredictionList& predictions) {
    for (auto& p : predictions) {
      if (!p->flag || p->flow_cov.empty()) continue;
      bool has_cell = false;
      for (float coverage : p->cell_cov)
        if (coverage > 0) has_cell = true;
      if (has_cell) continue;

      // The cell with the most flow coverage
      int32_t best_cell = 0;
      for (int32_t c = 1; c < count_of(p->flow_cov); c++)
        if (p->flow_cov[c] > p->flow_cov[best_cell]) best_cell = c;
      if (p->flow_cov[best_cell] <= 0) continue;

      p->cell_cov.resize(count_of(p->flow_cov), 0.0);
      p->cell_cov[best_cell] = p->flow_cov[best_cell];
      p->gene_no = 1;
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
  std::vector<int32_t> group_predictions_into_genes(
      PredictionList& predictions) {
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
  std::vector<CellJunctionReads> collect_cell_junction_reads(
      BundleData *bundle, int32_t cell_count) {
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
  // cell leaves the prediction, which is dropped when no cell is left.
  // Reference predictions are not touched.
  void drop_cells_with_competing_junctions(BundleData *bundle,
                                           int32_t cell_count) {
    std::vector<CellJunctionReads> cell_reads
        = collect_cell_junction_reads(bundle, cell_count);

    for (auto& prediction : bundle->predictions) {
      if (!prediction->flag || prediction->t_eq
          || count_of(prediction->exons) < 2)
        continue;

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

  // Coverage (all reads, per bp summed) of the genomic region [start, end],
  // clipped to the bundle: outside the bundle there are no reads.
  double bundle_region_coverage(BundleData *bundle, int64_t start,
                                int64_t end) {
    int64_t lo = std::max<int64_t>(start, bundle->start);
    int64_t hi = std::min<int64_t>(end, bundle->end);
    if (hi < lo) return 0;
    return get_region_coverage('.', (uint32_t)(lo - bundle->start),
                               (uint32_t)(hi - bundle->start),
                               bundle->bp_coverage);
  }

  // True if the single-exon prediction p does not stand out from the
  // coverage around it (max_single_exon_background).
  bool is_single_exon_in_background(BundleData *bundle, const CPrediction& p) {
    int64_t start = p.exons[0].start, end = p.exons[0].end;
    double inside = bundle_region_coverage(bundle, start, end)
        / (end - start + 1);
    if (inside <= 0) return true;

    // The two windows beyond the ends, the same length on each side
    int64_t gap = single_exon_background_gap;
    int64_t reach = single_exon_background_reach;
    double around = bundle_region_coverage(bundle, start - reach, start - gap)
        + bundle_region_coverage(bundle, end + gap, end + reach);
    double background = around / (2 * (reach - gap + 1));
    return background > max_single_exon_background * inside;
  }

  // Drops the novel single-exon predictions that do not look like transcripts
  // of their own: those that do not stand out from the coverage around them
  // and, with --sq3, those with too few UMI reads starting at their 5' end.
  // Reference ones are kept. On the real data most dropped ones are unspliced
  // pre-mRNA inside introns (ribbon_data/benchmark/README.md).
  void drop_unsupported_single_exons(BundleData *bundle) {
    PredictionList& predictions = bundle->predictions;
    bool any = false;
    for (const auto& p : predictions)
      if (p->flag && !p->t_eq && count_of(p->exons) == 1) any = true;
    if (!any) return;

    // With --sq3: the 5' ends of the UMI read 1s of each strand
    PerStrand<ReadStarts> umi_starts;
    if (smartseq3) {
      for (const auto& read : bundle->reads) {
        if (read->ub.empty() || !read->is_read1) continue;
        if (read->strand != '+' && read->strand != '-') continue;
        int32_t five_prime = read->strand == '+' ? read->start : read->end;
        umi_starts[read->strand].by_position.push_back(
            {five_prime, read->read_count});
      }
      for (char strand : {'-', '+'})
        std::sort(umi_starts[strand].by_position.begin(),
                  umi_starts[strand].by_position.end());
    }

    for (auto& p : predictions) {
      if (!p->flag || p->t_eq || count_of(p->exons) != 1) continue;
      if (is_single_exon_in_background(bundle, *p)) {
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

  // Decides which predictions of the bundle are kept: sorts them by position,
  // merges the copies found in different cells, and clears the flag of the
  // ones that are dropped. The flagged predictions remain sorted by position
  // and carry their per-cell coverage (cell_cov).
  void finalize_predictions(BundleData *bundle) {
    PredictionList& predictions = bundle->predictions;
    sort_by(predictions, compare_prediction_by_position);

    // Cells are numbered in the order they were first seen.
    int32_t cell_count = count_of(bundle->cell_name);
    merge_cell_copies(predictions, cell_count);

    // A cell's own reads can contradict a transcript the pooled graph gave it.
    drop_cells_with_competing_junctions(bundle, cell_count);

    // Every kept transcript has at least one cell.
    credit_best_cell_to_guides_without_cells(predictions);
    drop_predictions_without_cells(predictions);

    // Filter overlapping predictions, from the most covered to the least.
    std::vector<CoverageRank> coverage_order = build_coverage_order(predictions);
    filter_overlapping_predictions(predictions, coverage_order, cell_count);

    drop_weak_predictions(predictions);
    drop_weak_isoforms(predictions);

    // Single exons need coverage of their own (and 5' reads with --sq3).
    drop_unsupported_single_exons(bundle);
  }

} // namespace ribbon
