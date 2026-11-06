#pragma once

#include "ribbon.h"

// Per-base read coverage (bp_coverage) of a bundle: building it from
// the reads, answering region queries, and finding the coverage-drop trim
// points of the splice graph nodes.

namespace ribbon {
  // Block size (bp) of the coverage index in bp_coverage (see
  // compute_region_coverage)
  const int32_t BSIZE = 10000;

  void add_coverage_diff(BpCoverage& bp_coverage, char strand,
                         int32_t start, int32_t end, float delta);

  void accumulate_read_coverage(ReadList& reads, int32_t read_idx,
                                BpCoverage& bp_coverage,
                                BundleData *bundle);

  double get_region_coverage(char strand, uint32_t start, uint32_t end,
                             BpCoverage& bp_coverage);

  double get_stranded_coverage(char strand, uint32_t start, uint32_t end,
                               BpCoverage& bp_coverage);

  void find_all_trim_points(int32_t ref_start, char strand,
                            uint32_t start,
                            uint32_t end, BpCoverage& bp_coverage,
                            std::vector<CTrimPoint>& trim_point);

  CGraphnode *trim_node_at_all_points(const GraphBuildContext& context,
      uint32_t new_end, CGraphnode *graph_node,
      PendingTransfragList& pending_transfrags, int32_t& node_count,
      PreGraphSegment *pre_graph_segment, int32_t& edge_count);
} // namespace ribbon
