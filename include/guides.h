#pragma once

#include "ribbon.h"

// Guide-transcript (reference annotation) matching against the splice graph;
// see guides.cpp.

namespace ribbon {
  std::shared_ptr<CTransfrag> match_guide_to_graph(
      GffObj *guide, SpliceGraph& graph, int32_t cell_count);

  int32_t compare_guide_by_path_bits_size(const void *p1, const void *p2);

  void process_reference_guides(SpliceGraph& graph, char strand,
                                std::vector<CGuide>& matched_guides,
                                BundleData *bundle);

  bool is_within_guide_exon(std::vector<GffObj *>& guides, char strand,
                            uint32_t start, uint32_t end);

  void set_guide_status(GffObj *guide, GuideBundleStatus status);

  void compute_guide_boundaries(
      std::vector<GffObj *>& guides,
          std::vector<GuideBoundary>& guide_boundaries,
      std::vector<int32_t>& guide_prediction);

  int32_t count_guides_overlapping_pre_graph(
      std::vector<GffObj *>& guides, char strand_char, uint32_t pre_graph_start,
      uint32_t pre_graph_end, int32_t& guide_cursor);
} // namespace ribbon
