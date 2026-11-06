#pragma once

#include "ribbon.h"

// Pruning of a splice graph that has too many nodes: removes the weakest
// junctions, merges the nodes left singly linked and renumbers the graph.

namespace ribbon {
  int32_t compare_graph_junction_by_coverage(const void *p1,
                                             const void *p2);

  void remove_graph_edge(std::vector<CGNode>& node_state, int32_t parent_node,
                         int32_t child_node, int32_t pending_idx,
                             GraphNodeList& nodes,
                         int32_t& edge_count, int32_t& estimated_node_count,
                         CGraphnode *sink,
                         PendingTransfragList& pending_transfrags);

  int32_t count_nodes_after_pruning(std::vector<CGNode>& node_state,
                                    int32_t node_count,
                                GraphNodeList& nodes);

  int32_t prune_graph_to_node_limit(const GraphBuildContext& context,
                                int32_t node_count,
                                    int32_t pre_graph_segment_count,
                                JunctionList& junction, int32_t& edge_count,
                                PendingTransfragList& pending_transfrags);
} // namespace ribbon
