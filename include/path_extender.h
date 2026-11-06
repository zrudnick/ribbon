#pragma once

#include "ribbon.h"

// Greedy extension of a seed node into a source-to-sink path, choosing at
// each step the best supported neighbor.

namespace ribbon {
  // Extends a path (ending at node i) one node at a time toward the source
  // (toward_source) or the sink, picking at each step the neighbor that the
  // compatible transfrags support best (cell abundance^alpha * pooled
  // abundance^beta) and avoiding low-coverage ones. When going toward the
  // source the nodes are added at the end of path (the caller reverses it).
  class PathExtender {
    const FlowContext& context;
    int32_t i; // the node where the path ends now
    bool toward_source;
    std::vector<int32_t>& path;
    BitVec& current_path_bits;
    std::vector<float>& node_coverage;
    bool& path_empty;
    // Set when a spliced step was a weak borrowing: chosen by pooled support
    // among several supported neighbors, with less than min_borrowed_share
    // of their support (the caller then does not store the path)
    bool& weak_borrow;
    BitVec& low_cov_flags;

    // What differs between the two directions: the node we walk to, and the
    // step to the neighbor next to i.
    int32_t end_node;
    int32_t step;

    // The choice of the next node at the current node: the best neighbor by
    // support in this cell, or else by pooled support. If the adjacent node was
    // excluded as a coverage drop, it is assessed again when no neighbor has
    // support.
    int32_t best_neighbor = -1;
    // If there are missing spliced reads between node and neighbor we still
    // want to consider them
    int32_t best_pooled_neighbor = -1;
    float best_support = 0;
    float best_pooled_support = 0;
    // The pooled support summed over the neighbors, and how many have some
    float pooled_support_sum = 0;
    int32_t pooled_supported_count = 0;
    bool adjacent_excluded = false;
    std::vector<std::pair<int32_t, int32_t *>> scored_neighbors;

    // How one step ended: the walk goes on, it reached the source (or sink), or
    // there is no way.
    enum class StepResult { Continue, ReachedEnd, NoWay };

    int32_t *find_edge_bit(int32_t toward, int32_t away);

    void report_missing_edge(int32_t toward, int32_t source_code,
                             int32_t sink_code);

    void set_path_bits(int32_t node, int32_t *edge_pos, bool value);

    void add_transfrag_support(int32_t neighbor, float& support,
                               float& pooled_support, bool want_support,
                               bool want_pooled_support);

    bool is_excluded_adjacent_neighbor(int32_t neighbor,
                                       CGraphnode *current_node,
                                       CGraphnode *neighbor_node,
                                       bool adjacent_to_neighbor);

    void score_neighbors_by_cell_support(CGraphnode *current_node,
                                         std::vector<int32_t>& neighbors);

    void score_neighbors_by_pooled_support();
    void check_weak_borrow();
    StepResult choose_neighbor_without_cell_support();

    void move_to_best_neighbor();
    StepResult walk_one_node();

  public:
    PathExtender(const FlowContext& context_in, int32_t i_in,
                 bool toward_source_in, std::vector<int32_t>& path_in,
                 BitVec& current_path_bits_in,
                 std::vector<float>& node_coverage_in,
                 bool& path_empty_in, bool& weak_borrow_in,
                 BitVec& low_cov_flags_in);

    bool walk();
  };

  bool is_transfrag_on_path(const FlowContext& context,
                            BitVec& transfrag_path_bits,
                            std::vector<int32_t>& transfrag_nodes,
                            BitVec& current_path_bits, int32_t first_node,
                            int32_t last_node);

  bool extend_path(const FlowContext& context, int32_t i, bool toward_source,
                   std::vector<int32_t>& path, BitVec& current_path_bits,
                   std::vector<float>& node_coverage, bool& path_empty,
                   bool& weak_borrow, BitVec& low_cov_flags);
} // namespace ribbon
