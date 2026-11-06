#pragma once

#include "ribbon.h"

// Joins read groups into pre-graphs (chains of PreGraphSegments), one per
// stranded locus, that each become a splice graph.

namespace ribbon {
  int32_t start_new_pre_graph(PreGraphList& pre_graph_list, CGroup *group,
                       PreGraphSegmentList& pre_graph_segments);

  void add_group_to_pre_graph(CGroup *group, PreGraph *pre_graph,
                                PreGraphSegmentList& pre_graph_segments,
                                uint32_t merge_distance);

  void propagate_strand_link_id(CGroup *prev_group, CGroup *group,
                                   int32_t base_link_id,
                                   std::vector<int32_t>& strand_link_id_root,
                                   std::vector<int32_t>& link_id_root);

  char next_group_by_start(PerStrand<CGroup *>& curr_group);

  void form_pre_graphs_from_groups(AssemblyState& state);
} // namespace ribbon
