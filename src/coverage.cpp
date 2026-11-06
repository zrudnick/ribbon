#include "ribbon.h"
#include "GThreads.h"

// Per-base read coverage (bp_coverage) of a bundle: building it from
// the reads, answering region queries, and finding the coverage-drop trim
// points of the splice graph nodes.

namespace ribbon {

  // Records a coverage change of delta over [start,end) as +delta/-delta at the
  // interval's ends in the difference array bp_coverage[strand] (a prefix
  // sum gives per-base coverage). Stranded reads are also added to '.'.
  void add_coverage_diff(BpCoverage& bp_coverage, char strand,
                         int32_t start, int32_t end, float delta) {
    bool fold_into_total = false;
    // Stranded (-/+) coverage is also added to the all-reads total in
    // bp_coverage['.']
    if (strand != '.') fold_into_total = true;
    // If strand is '.' then we add delta to it here
    bp_coverage[strand][start + 1] += delta;
    bp_coverage[strand][end + 1] -= delta;
    // If stranded, it gets added to bp_coverage['.'] here too, so that
    // bp_coverage['.'] is all coverage
    if (fold_into_total) {
      bp_coverage['.'][start + 1] += delta;
      bp_coverage['.'][end + 1] -= delta;
    }
  }

  // Adds the coverage of a read and of the mates that come before it to
  // bp_coverage, one segment at a time; overlapping segments of read and mate
  // are merged first so that shared bases are not counted twice.
  void accumulate_read_coverage(ReadList& reads, int32_t read_idx,
                                BpCoverage& bp_coverage,
                                BundleData *bundle) {
    char strand = reads[read_idx]->strand;

    // Start with the whole read count as unpaired, then subtract the pairs.
    float single_count = reads[read_idx]->read_count;
    for (int32_t i = 0; i < count_of(reads[read_idx]->pair_idx); i++) {
      single_count -= reads[read_idx]->pair_count[i];
      if (reads[read_idx]->pair_idx[i] <= read_idx) { // pair comes before
        int32_t mate_idx = reads[read_idx]->pair_idx[i];
        int32_t mate_seg_count = count_of(reads[mate_idx]->segs);
        int32_t read_seg_count = count_of(reads[read_idx]->segs);
        char mate_strand = reads[mate_idx]->strand;
        int32_t mate_seg_idx = 0;
        int32_t read_seg_idx = 0;
        int32_t segs_remaining = mate_seg_count + read_seg_count;
        while (segs_remaining) {
          int32_t start;
          int32_t end;
          if (read_seg_idx < read_seg_count) {
            // No more pair segments or pair segment comes after read segment
            if (mate_seg_idx == mate_seg_count
              || reads[mate_idx]->segs[mate_seg_idx].start
              > reads[read_idx]->segs[read_seg_idx].end) {

              start = reads[read_idx]->segs[read_seg_idx].start;
              end = reads[read_idx]->segs[read_seg_idx].end;
              read_seg_idx++;
              segs_remaining--;

            } else if (reads[mate_idx]->segs[mate_seg_idx].end
              < reads[read_idx]->segs[read_seg_idx].start) {

              start = reads[mate_idx]->segs[mate_seg_idx].start;
              end = reads[mate_idx]->segs[mate_seg_idx].end;
              mate_seg_idx++;
              segs_remaining--;

            } else { // segments overlap

              start = reads[mate_idx]->segs[mate_seg_idx].start;
              end = reads[mate_idx]->segs[mate_seg_idx].end;
              if ((int32_t)reads[read_idx]->segs[read_seg_idx].start < start)
                start = reads[read_idx]->segs[read_seg_idx].start;
              mate_seg_idx++;
              segs_remaining--;
              bool merging = true;

              while (merging) {
                merging = false;
                if (read_seg_idx < read_seg_count
                    && (int32_t)(reads[read_idx]->segs[read_seg_idx].start)
                    <= end) {
                  if ((int32_t)reads[read_idx]->segs[read_seg_idx].end > end)
                    end = reads[read_idx]->segs[read_seg_idx].end;
                  read_seg_idx++;
                  segs_remaining--;
                  merging = true;
                }
                if (mate_seg_idx < mate_seg_count
                    && (int32_t)(reads[mate_idx]->segs[mate_seg_idx].start)
                    <= end) {
                  if ((int32_t)reads[mate_idx]->segs[mate_seg_idx].end > end)
                    end = reads[mate_idx]->segs[mate_seg_idx].end;
                  mate_seg_idx++;
                  segs_remaining--;
                  merging = true;
                }
              }
            }
          } else {
            start = reads[mate_idx]->segs[mate_seg_idx].start;
            end = reads[mate_idx]->segs[mate_seg_idx].end;
            mate_seg_idx++;
            segs_remaining--;
          }
          char merged_strand = strand;
          if (strand != mate_strand) {
            if (strand == '.') merged_strand = mate_strand;
            else if (mate_strand != '.') merged_strand = '.';
          }
          // Both the read and it's pair are added here; if pairs overlap then
          // bp_coverage only counts one coverage not both
          add_coverage_diff(bp_coverage, merged_strand, start
                            - bundle->start,
                            end - bundle->start + 1,
                            reads[read_idx]->pair_count[i]);
        }
      }
    }

    // Add the part of the read count that is not in any pair.
    if (single_count > epsilon) {
      for (int32_t i = 0; i < count_of(reads[read_idx]->segs); i++) {
        add_coverage_diff(bp_coverage, strand,
          reads[read_idx]->segs[i].start - bundle->start,
          reads[read_idx]->segs[i].end - bundle->start + 1, 
          single_count);
      }
    }
  }

  // Total coverage over [start,end]. bp_coverage[strand] ('.' for all reads)
  // holds a running sum that restarts every BSIZE bp, so we add the whole
  // blocks and the two partial blocks.
  static inline double compute_region_coverage(char strand,
                                               uint32_t start,
                                               uint32_t end,
                                               BpCoverage& bp_coverage,
                                               bool stranded) {
    int32_t end_block = int32_t((end + 1) / BSIZE);
    int32_t start_block = int32_t(start / BSIZE);
    double cov = 0;
    if (stranded) {
      const std::vector<float>& all_reads = bp_coverage['.'];
      const std::vector<float>& opposite = bp_coverage[opposite_strand(strand)];
      for (int32_t i = start_block + 1; i <= end_block; i++)
        cov += all_reads[i * BSIZE - 1] - opposite[i * BSIZE - 1];
      cov += all_reads[end + 1] - all_reads[start]
          - opposite[end + 1] + opposite[start];
    } else {
      const std::vector<float>& coverage = bp_coverage[strand];
      for (int32_t i = start_block + 1; i <= end_block; i++)
        cov += coverage[i * BSIZE - 1];
      cov += coverage[end + 1] - coverage[start];
    }
    if (cov < 0) cov = 0;
    return (cov);
  }

  // Coverage of [start,end] of the given strand in bp_coverage ('.' = all
  // reads, regardless of strand).
  double get_region_coverage(char strand, uint32_t start, uint32_t end,
                             BpCoverage& bp_coverage) {
    return compute_region_coverage(strand, start, end, bp_coverage, false);
  }

  // Like get_region_coverage, but returns the coverage attributable to
  // strand alone (all-reads total bp_coverage['.'] minus the coverage of the
  // opposite strand).
  double get_stranded_coverage(char strand, uint32_t start, uint32_t end,
                               BpCoverage& bp_coverage) {
    return compute_region_coverage(strand, start, end, bp_coverage, true);
  }

  // Records or updates a source/sink trim candidate at position i from the
  // drop ratios; abundance_divisor scales the recorded abundance (DROP
  // normally, DROP * window_len where the coverages are raw sums).
  static void update_trim_point_from_drop(std::vector<CTrimPoint>& trim_point,
                                          uint32_t i, double coverage_left,
                                          double coverage_right,
                                          double source_ratio,
                                          double sink_ratio, double local_drop,
                                          double abundance_divisor,
                                          double& last_drop) {
    if (coverage_left < coverage_right) { // possible source trimming
      if (source_ratio < local_drop) { // this drop passes threshold
        // Add this point to trimpoints
        // (the first point, or a source point far enough after a sink point)
        bool starts_new_trim_point = !count_of(trim_point)
          || (!trim_point.back().start
              && i + 1 - (int32_t)trim_point.back().pos > CHI_THR);

        if (starts_new_trim_point) {
          CTrimPoint t(i + 1,
            (coverage_right - coverage_left) / abundance_divisor,
            true);

          trim_point.push_back(t);
          last_drop = source_ratio;
        // Smaller drop than before --> replace drop; it also has to be smallest
        // within a 50bp prevwindow
        } else if (source_ratio < last_drop) {
          last_drop = source_ratio;
          trim_point.back().pos = i + 1;
          trim_point.back().abundance = (coverage_right - coverage_left)
            / abundance_divisor;
          trim_point.back().start = true;
        }
      }
    // Possible sink trimming: coverage_right<coverage_left here
    } else if (coverage_right != coverage_left) {
      if (sink_ratio < local_drop) { // this drop passes threshold
        // Add this point to trimpoints
        // (the first point, or a sink point far enough after a source point)
        bool starts_new_trim_point = !count_of(trim_point)
          || (trim_point.back().start
              && i - (int32_t)trim_point.back().pos > CHI_THR);
        if (starts_new_trim_point) {
          CTrimPoint t(i, (coverage_left - coverage_right) / abundance_divisor,
            false);
          trim_point.push_back(t);
          last_drop = sink_ratio;
        // Smaller drop than before --> replace drop
        } else if (sink_ratio < last_drop) {
          last_drop = sink_ratio;
          trim_point.back().pos = i;
          trim_point.back().abundance = (coverage_left - coverage_right)
            / abundance_divisor;
          trim_point.back().start = false;
        }
      }
    }
  }

  // Finds every coverage-drop point of the exon [start,end] that looks like a
  // source or sink trim and appends it to trim_point (the scan window depends
  // on the exon length); finally weak candidates are dropped on re-check.
  void find_all_trim_points(int32_t ref_start, char strand,
                            uint32_t start,
                            uint32_t end, BpCoverage& bp_coverage,
                            std::vector<CTrimPoint>& trim_point) {
    // Length of the region where we look for trims
    int32_t region_len = end - start + 1;
    // Very short exon -> do not check. The scan below skips long_intron_anchor
    // bp at each end, so an exon this short has no positions to check anyway
    // (it gives the same output with or without this return).
    if (region_len < CHI_THR) return;

    double local_drop = ERROR_PERC / DROP;
    if (region_len < 2 * (CHI_WIN + CHI_THR) + 1) {
      // Check very tightly, and only accept one trim
      if (region_len < CHI_WIN) {
        local_drop = ERROR_PERC / (10 * DROP); // very sharp drop
      }
      double last_drop = local_drop;
      uint32_t loose_start = start + long_intron_anchor;
      uint32_t loose_end = end - long_intron_anchor;
      for (uint32_t i = loose_start; i < loose_end; i++) {
        // Shouldn't this go to i?
        double coverage_left = get_stranded_coverage(strand,
            start - ref_start, i - 1 - ref_start, bp_coverage) / (i - start);
        double coverage_right = get_stranded_coverage(strand,
            i - ref_start, end - ref_start, bp_coverage) / (end - i + 1);

        // Make sure we add one read to mitigate gaps in coverage
        double source_ratio = coverage_left / coverage_right;
        double sink_ratio = (coverage_right + 1) / (coverage_left + 1);
        update_trim_point_from_drop(trim_point, i, coverage_left,
                                    coverage_right, source_ratio, sink_ratio,
                                    local_drop, DROP, last_drop);
      }
      return;
    }

    // Longer exon: scan the beginning, the middle and the end with their own
    // windows.
    int32_t window_len = CHI_WIN + CHI_THR;
    double last_drop = local_drop;

    for (uint32_t i = start + CHI_THR - 1; i < start + window_len - 1; i++) {
      double coverage_left = get_stranded_coverage(strand,
          start - ref_start, i - ref_start, bp_coverage) / (i - start);
      double coverage_right = get_stranded_coverage(strand,
          i + 1 - ref_start, i + window_len - ref_start,
          bp_coverage) / window_len;

      double source_ratio = coverage_left / coverage_right;
      double sink_ratio = coverage_right / coverage_left;
      update_trim_point_from_drop(trim_point, i, coverage_left, coverage_right,
                                  source_ratio, sink_ratio, local_drop, DROP,
                                  last_drop);
    }

    // We might want to use different parameters here
    local_drop = ERROR_PERC / DROP;
    if (!count_of(trim_point)) last_drop = local_drop;
    region_len += start - window_len;
    for (uint32_t i = start + window_len - 1; i < (uint32_t)region_len; i++) {
      // We add 1 bp coverage to all so we can avoid 0 coverages
      double coverage_left = get_stranded_coverage(strand,
          i - window_len + 1 - ref_start, i - ref_start, bp_coverage);
      double coverage_right = get_stranded_coverage(strand,
          i + 1 - ref_start, i + window_len - ref_start, bp_coverage);
      double source_ratio = coverage_left / coverage_right;
      double sink_ratio = coverage_right / coverage_left;
      update_trim_point_from_drop(trim_point, i, coverage_left, coverage_right,
                                  source_ratio, sink_ratio, local_drop,
                                  DROP * window_len, last_drop);
    }

    local_drop = ERROR_PERC;

    for (uint32_t i = region_len; i < end - CHI_THR - 1; i++) {
      double coverage_left = get_stranded_coverage(strand,
          i - window_len + 1 - ref_start, i - ref_start,
          bp_coverage) / window_len;
      double coverage_right = get_stranded_coverage(strand,
          i + 1 - ref_start, end - ref_start, bp_coverage) / (end - i);
      double source_ratio = coverage_left / coverage_right;
      double sink_ratio = coverage_right / coverage_left;
      update_trim_point_from_drop(trim_point, i, coverage_left, coverage_right,
                                  source_ratio, sink_ratio, local_drop, DROP,
                                  last_drop);
    }

    // Re-check each candidate against its neighbors (up to 3 * CHI_WIN away)
    // and discard the ones that do not hold up.
    local_drop = ERROR_PERC / DROP;
    // We might adjust local_drop again here
    uint32_t last_start = start;
    for (int32_t i = 0; i < count_of(trim_point); i++) {
      uint32_t mid_pos = trim_point[i].pos;
      if (trim_point[i].start) mid_pos--;
      // Do not consider a two big of a window
      if (mid_pos - last_start > 3 * CHI_WIN) { last_start = mid_pos - 3
                                               * CHI_WIN; }
      // We add 1 bp coverage to all so we can avoid 0 coverages
      double coverage_left = get_stranded_coverage(strand,
          last_start - ref_start, mid_pos - ref_start,
          bp_coverage) / (mid_pos - last_start + 1);
      uint32_t end_pos = end;
      if (i < count_of(trim_point) - 1) end_pos = trim_point[i + 1].pos;
      // Do not consider a two big of a window
      if (end_pos - mid_pos > 3 * CHI_WIN) { end_pos = mid_pos + 3 * CHI_WIN; }
      double coverage_right = get_stranded_coverage(strand,
          mid_pos + 1 - ref_start, end_pos - ref_start,
          bp_coverage) / (end_pos - mid_pos);

      // Source|sink above threshold
      bool is_source_drop = trim_point[i].start
                            && coverage_left < local_drop * coverage_right;
      bool is_sink_drop = coverage_right < local_drop * coverage_left;
      if (is_source_drop || is_sink_drop) {
        last_start = trim_point[i].pos + 1;
      // The trim_point is not above threshold --> do not change last_start, but
      // we need to reevaluate previous trimpoints
      } else {
        trim_point[i].pos = 0;
        int32_t review_idx = i - 1;
        while (review_idx >= 0) {
          while (review_idx >= 0 && !trim_point[review_idx].pos) review_idx--;
          if (review_idx >= 0) { // found a point above threshold
            mid_pos = trim_point[review_idx].pos;
            if (trim_point[review_idx].start) mid_pos--;
            uint32_t left_start = start;
            int32_t j = review_idx - 1;
            while (j >= 0) {
              if (trim_point[j].pos) {
                left_start = trim_point[j].pos;
                break;
              }
              j--;
            }
            // Do not consider a two big of a window
            if (mid_pos - left_start > 3 * CHI_WIN) { left_start = mid_pos - 3
                * CHI_WIN; }
            // We add 1 bp coverage to all so we can avoid 0 coverages
            coverage_left = get_stranded_coverage(strand,
                                                  left_start - ref_start,
                                                  mid_pos - ref_start,
                                                  bp_coverage)
                                                  / (mid_pos - left_start + 1);
            // Do not consider a two big of a window
            if (end_pos - mid_pos > 3 * CHI_WIN) { end_pos = mid_pos + 3
                                                  * CHI_WIN; }
            coverage_right = get_stranded_coverage(strand,
                mid_pos + 1 - ref_start, end_pos - ref_start,
                bp_coverage) / (end_pos - mid_pos);
            // Source|sink above threshold
            bool is_earlier_source_drop
                = trim_point[review_idx].start
                  && coverage_left < local_drop * coverage_right;
            bool is_earlier_sink_drop
                = coverage_right < local_drop * coverage_left;
            if (is_earlier_source_drop || is_earlier_sink_drop) {
              break; // finished looking for review_idx
            } else {
              trim_point[review_idx].pos = 0;
            }
            review_idx--;
          }
        }
      }
    }
  }

  // Splits graph_node at its coverage-drop trim points (find_all_trim_points),
  // left to right: a source trim links the new right piece from the source, a
  // sink trim links the left piece to the sink. Returns the rightmost piece.
  CGraphnode *trim_node_at_all_points(const GraphBuildContext& context,
      uint32_t new_end, CGraphnode *graph_node,
      PendingTransfragList& pending_transfrags, int32_t& node_count,
      PreGraphSegment *pre_graph_segment, int32_t& edge_count) {
    // Find the trim points of the node.
    std::vector<CTrimPoint> trim_point;
    find_all_trim_points(context.ref_start, context.strand,
                         graph_node->start, new_end,
                         *context.bp_coverage, trim_point);

    // Split the node at each trim point and record its source or sink link.
    for (int32_t trim_i = 0; trim_i < count_of(trim_point); trim_i++)
      if (trim_point[trim_i].pos) {
        if (trim_point[trim_i].start) { // source trim
          graph_node->end = trim_point[trim_i].pos - 1;

          CGraphnode *prev_node = graph_node;
          graph_node = add_graph_node(context, trim_point[trim_i].pos, new_end,
                                      node_count, pre_graph_segment);
          node_count++;
          // This node is the child of source
          context.source->child.push_back(graph_node->node_id);
          // This node has source as parent
          graph_node->parent.push_back(context.source->node_id);
          // This node is the child of previous node
          prev_node->child.push_back(graph_node->node_id);
          // This node has as parent the previous node
          graph_node->parent.push_back(prev_node->node_id);
          pending_transfrags.push_back({0.0f, float(node_count - 1),
                                        float(trim_point[trim_i].abundance)});
          // The link between the two halves has abundance 0: reads that span
          // the split make their own transfrag. StringTie gives it 1.0 (and 1.0
          // more to the source/sink link); we tested it, no difference
          pending_transfrags.push_back({float(prev_node->node_id),
                                        float(graph_node->node_id), 0.0f});
          // COUNT 2 EDGES HERE
          edge_count += 2;
        } else { // this is sink
          graph_node->end = trim_point[trim_i].pos;
          CGraphnode *prev_node = graph_node;
          graph_node = add_graph_node(context, trim_point[trim_i].pos + 1,
                                      new_end, node_count, pre_graph_segment);
          node_count++;
          // This node is the child of previous node
          prev_node->child.push_back(graph_node->node_id);
          // This node has as parent the previous node
          graph_node->parent.push_back(prev_node->node_id);
          // The prev_node is the parent of sink
          context.sink->parent.push_back(prev_node->node_id);
          // Remember to create transfrag as well -> we don't know the gno yet,
          // so we can not create it here
          pending_transfrags.push_back({float(node_count - 2), -1.0f,
                                        float(trim_point[trim_i].abundance)});
          // The link between the two halves has abundance 0: reads that span
          // the split make their own transfrag. StringTie gives it 1.0 (and 1.0
          // more to the source/sink link); we tested it, no difference
          pending_transfrags.push_back({float(prev_node->node_id),
                                        float(graph_node->node_id), 0.0f});
          // COUNT 2 EDGES HERE
          edge_count += 2;
        }
      }

    return (graph_node);
  }

} // namespace ribbon
