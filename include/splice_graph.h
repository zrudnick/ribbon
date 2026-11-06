#pragma once

#include "ribbon.h"

// Builds the splice graph of one pre-graph: nodes between junctions, the
// source and sink links, the guide boundaries and the parent/child bit vectors.

namespace ribbon {
  // Builds the splice graph of one pre-graph on a strand. The pre-graph's
  // segments are walked left to right while two cursors move over the junction
  // starts and the junction ends; each stretch between junction events becomes
  // a node.
  class SpliceGraphBuilder {
    BundleData *bundle;
    AssemblyState& state;
    SpliceGraph& graph;
    char strand;
    int32_t graph_idx;
    int32_t& edge_count;
    std::vector<GuideBoundary>& guide_boundaries;

    int32_t ref_start;
    JunctionList& junction;
    int32_t n_junctions;
    PreGraph *pre_graph;
    PreGraphSegmentList& pre_graph_segments;
    uint32_t pre_graph_end; // end of the last segment of the pre-graph
    GraphNodeList& nodes; // this graph's nodes, indexed by node id
    BpCoverage& bp_coverage;
    int32_t cell_count; // number of cells (barcodes) in the bundle
    CGraphnode *source; // node 0
    // The sink gets its real id (node_count) and joins `nodes` only at the end.
    std::unique_ptr<CGraphnode> sink_owner;
    CGraphnode *sink;
    GraphBuildContext context; // what the node-building helpers need

    int32_t node_count = 1; // number of nodes in the graph
    int32_t junc_start_idx = 0; // index of sorted junction starts
    int32_t junc_end_idx = 0; // index of sorted junction ends
    int32_t guide_boundary_idx = 0;
    bool pre_graph_has_guide_boundary = false;
    // Junction end position -> ids of the nodes that have a junction starting
    // toward that position. Filled when a junction start closes a node; read
    // when the walk reaches the junction's end.
    std::unordered_map<int32_t, std::vector<int32_t>> junction_origins;
    // Links to the source and sink made by guide, trim and coverage splits;
    // they are built once node_count is final. Triples (from, to, abundance):
    // to = -1 is the sink, from = 0 the source; from = -1 marks a removed one.
    PendingTransfragList pending_transfrags;

    // The segment being walked, its first node start and end, and the node
    // that is still open.
    PreGraphSegment *pre_graph_segment = NULL;
    uint32_t curr_start = 0;
    uint32_t pre_graph_segment_end = 0;
    CGraphnode *graph_node = NULL;
    // No node is left open in the segment
    bool pre_graph_segment_closed = false;

    void find_guide_boundary_in_pre_graph();
    bool next_junction_end_is_by(uint32_t pos);
    bool next_junction_start_is_by(uint32_t pos);
    bool read_junction_ends_up_to(uint32_t position);

    void skip_ends_on_other_strands();
    void skip_starts_on_other_strands();
    bool link_from_junction_origins(uint32_t pos, CGraphnode *node);

    void open_next_node(uint32_t start);
    void close_open_node_at(uint32_t end_pos);

    bool skip_hanging_piece_at_segment_start();
    void open_first_node(bool junction_ends_here);

    int32_t next_junction_event();
    bool tail_has_coverage_drop(uint32_t pos);

    void process_junction_start();
    void process_junction_end();

    void process_segment();
    int32_t finish_graph();

  public:
    SpliceGraphBuilder(BundleData *bundle_in, AssemblyState& state_in,
                       GraphTables& tables_in, char strand_in,
                       int32_t graph_idx_in,
                       std::vector<GuideBoundary>& guide_boundaries_in);

    int32_t build();
  };

  CGraphnode *add_graph_node(const GraphBuildContext& context, uint32_t start,
                             uint32_t end, int32_t node_id,
                             PreGraphSegment *pre_graph_segment);

  CGraphnode *link_source_to_guide_start(const GraphBuildContext& context,
      uint32_t new_start, uint32_t new_end, CGraphnode *graph_node,
      PendingTransfragList& pending_transfrags, int32_t& node_count,
      PreGraphSegment *pre_graph_segment, int32_t& edge_count);

  CGraphnode *link_guide_end_to_sink(const GraphBuildContext& context,
      uint32_t new_start, uint32_t new_end, CGraphnode *graph_node,
      PendingTransfragList& pending_transfrags, int32_t& node_count,
      PreGraphSegment *pre_graph_segment, int32_t& edge_count);

  BitVec build_reachable_bits_depth_first(
      const GraphBuildContext& context,
      CGraphnode *node, BitVec parents,
      int32_t node_count,
      std::vector<bool>& visited,
      int32_t& edge_count);

  int32_t build_splice_graph(BundleData *bundle, AssemblyState& state,
                         GraphTables& tables, char strand, int32_t graph_idx,
                         std::vector<GuideBoundary>& guide_boundaries);
} // namespace ribbon
