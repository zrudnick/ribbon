#pragma once

#include "ribbon.h"

// The types the graph-building stages share (groups, pre-graphs, graph nodes
// and tables); the stages are coverage.cpp, junctions.cpp, grouping.cpp,
// pre_graphs.cpp, splice_graph.cpp and graph_pruning.cpp.
// Conventions: a strand is '-' (negative), '+' (positive) or '.' (unknown),
// and the graphs are indexed [strand][graph_idx] (GraphTables::splice_graphs).
// Node id 0 is the source, 1..node_count-2 are the real nodes in genomic order
// and node_count-1 is the sink. Every BitVec path_bits has one bit per node
// followed by one bit per edge; path_bit_layout gives the bit of each edge.

namespace ribbon {
  // For each pre-graph segment id of a strand: the graph nodes that cover it.
  using SegmentNodesList = std::vector<std::vector<CGraphinfo>>;

  // A transfrag that build_splice_graph cannot create yet, because the final
  // node ids are unknown while the graph is built: a link between two nodes,
  // from the source (from == 0) or to the sink (to == -1). A negative from
  // marks a cancelled one.
  struct PendingTransfrag {
    float from;
    float to;
    float abundance;
    // True for the link a guide's start or end made to the source or sink.
    bool made_by_guide = false;
  };
  using PendingTransfragList = std::vector<PendingTransfrag>;

  struct SpliceGraph;

  // What the node-building helpers of build_splice_graph need to know about
  // the graph being built (one strand of one bundle).
  struct GraphBuildContext {
    char strand;
    int32_t graph_idx;
    int32_t ref_start;
    int32_t cell_count; // number of cells (barcodes) in the bundle
    BpCoverage *bp_coverage;
    SpliceGraph *graph; // the graph being built
    // pre_graph_segment_to_graph_nodes of the strand
    SegmentNodesList *segment_nodes;
    CGraphnode *source;
    CGraphnode *sink;
  };

  // This can work as a guide keeper too, where pos is the guideidx, abundance
  // is the flow, and start is the included status
  struct CTrimPoint {
    uint32_t pos;
    float abundance;
    bool start : 1;
    CTrimPoint(uint32_t _pos = 0, float abund = 0.0, bool _start = true)
        : pos(_pos), abundance(abund), start(_start) {}
  };

  struct CGNode {
    int32_t id; // initial id in graphno
    bool last; // if this is last node (to be linked to sink later)
    bool keep; // if we keep it in the final count (true by default)
    bool merge; // if this node needs to be merged to its adjacent node
    bool future;
    CGNode(int32_t _id = 0, bool _last = false, bool _keep = true,
           bool _merge = false, bool _future = false)
        : id(_id), last(_last), keep(_keep), merge(_merge), future(_future) {}
  };

  // The PreGraphs of one strand (the list owns them). A stranded pre-graph
  // that is worth processing gets one splice graph, with the pre-graph's index
  // as its graph index.
  using PreGraphList = std::vector<std::unique_ptr<PreGraph>>;

  // The PreGraphSegment chains of one strand, indexed by pre-graph segment id.
  // The list owns the segments.
  using PreGraphSegmentList = std::vector<std::unique_ptr<PreGraphSegment>>;

  // Shared per-bundle assembly state threaded through infer_transcripts's
  // stages (junction support propagation, read grouping, bundle formation,
  // graph construction). Grows as that function's monolithic body is decomposed
  // into named stage functions.
  struct AssemblyState {
    // This bundle's junctions, indexed/sorted by end coordinate
    std::vector<CJunction *> end_sorted_junctions;

    // Same-strand read groups formed by assign_read_to_groups: the group
    // list itself, the working/first CGroup per strand ('-', '.' unknown, '+'),
    // the next link_id to assign, each read's group membership, and the
    // union-find-style link_id merge map.
    GroupList group;
    PerStrand<CGroup *> curr_group;
    PerStrand<CGroup *> start_group;
    int32_t link_id = 0; // next unused link id
    std::vector<int32_t> merge_map; // remembers merged groups
    // Union-find root array: link_id_root[c] is the canonical link_id
    // that c has been merged into
    std::vector<int32_t> link_id_root;
    // One std::vector<int> per read, listing the groups it belongs to
    std::vector<std::vector<int32_t>> read_group;

    // Junction coordinates seen as a guide-derived group boundary while
    // assigning reads to groups; consulted (and cleared) when merging
    // adjacent/overlapping groups afterward, so a boundary that coincides with
    // a guide boundary is never merged across.
    std::unordered_set<int32_t> junction_start_positions;
    std::unordered_set<int32_t> junction_end_positions;

    // Per-strand ('-', '.' unknown, '+') PreGraphs formed from the merged
    // groups, their PreGraphSegment chains, and the group-id ->
    // pre-graph-segment-id lookup used to retrace which pre-graph segment a
    // given read's group ended up in.
    PerStrand<PreGraphList> pre_graph_list;
    PerStrand<PreGraphSegmentList> pre_graph_segment_list;
    PerStrand<std::vector<int32_t>> group_to_pre_graph_segment;

    AssemblyState(JunctionList& junction, int32_t n_reads)
        : read_group(n_reads) {
      for (auto& j : junction) end_sorted_junctions.push_back(j.get());
    }
  };

  // One splice graph of a strand: the pre-graph with the same index it is built
  // from. A pre-graph that is not worth processing keeps a node count of 0.
  struct SpliceGraph {
    int32_t node_count = 0;
    int32_t edge_count = 0;
    // The nodes of the graph; node i is nodes[i]
    GraphNodeList nodes;
    // The abundance and the path_bits of each transfrag of the graph
    TransfragList transfrag;
    // The bit of each node and edge of the graph in the path_bits of a
    // transfrag
    PathBitLayout path_bit_layout;
    // Finds the transfrag of the graph with a given path without scanning the
    // transfrag list
    TransfragIndex transfrag_index;
  };

  // The graphs of a bundle: splice_graphs is indexed by strand ('-' or '+')
  // and then by graph g (the index of the stranded PreGraph the graph is built
  // from).
  struct GraphTables {
    // Maps a pre-graph segment to the graph nodes that cover it
    PerStrand<SegmentNodesList> pre_graph_segment_to_graph_nodes;
    // The graphs of a strand, one slot per pre-graph (built or not)
    PerStrand<std::vector<SpliceGraph>> splice_graphs;
  };

  // One read, or one pair of reads, to place in the groups of a strand ('-',
  // '.' unknown or '+') with the coverage of this pair; the
  // link id is the one shared by the whole fragment.
  struct ReadPlacement {
    int32_t read_idx;
    int32_t mate_idx; // -1 if the read is placed alone
    int32_t pair_slot; // -1 if the read is placed alone
    float read_coverage;
    char strand;
    int32_t read_link_id;
  };

  int32_t predict_unstranded_pre_graphs(BundleData *bundle,
      AssemblyState& state);

  void get_splice_graphs(BundleData *bundle, AssemblyState& state,
                    GraphTables& tables,
                    std::vector<GuideBoundary>& guide_boundaries);

  void free_pre_graphs(AssemblyState& state);

  void fill_transfrags(BundleData *bundle, AssemblyState& state,
                       GraphTables& tables);

  void free_read_groups(AssemblyState& state);

  int32_t extract_transcripts_from_graphs(
      BundleData *bundle, GraphTables& tables,
      std::vector<int32_t>& guide_prediction, int32_t gene_count);

} // namespace ribbon
