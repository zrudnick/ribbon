#pragma once

#include "ribbon.h"

// Transfrag path index and transfrag simplification/stitching; see
// transfrag.cpp.

namespace ribbon {
  std::vector<int32_t> make_path_key(PathBitLayout& layout,
                                     std::vector<int32_t>& nodes,
                                     BitVec& path_bits);

  void remove_low_abundance_transfrags(SpliceGraph& graph, float threshold);

  int32_t compare_transfrag_by_reach(const void *p1, const void *p2);

  bool insert_transfrag_into_node_list(int32_t t,
                                       std::vector<int32_t>& transfrag,
                                       int32_t first,
                                       int32_t last);

  bool assign_incomplete_transfrag_to_nodes(int32_t t, int32_t from_node,
                                            int32_t to_node,
                                            GraphNodeList& graph_nodes);

  bool resolve_transfrag_paths(int32_t t, SpliceGraph& graph);

  void process_transfrags(SpliceGraph& graph,
                          std::vector<CGuide>& matched_guides,
                          BundleData *bundle);

  void get_read_graph_nodes(
      char strand, float read_coverage, std::vector<int32_t>& read_graph_idxs,
      float strand_prop, std::vector<std::vector<int32_t>>& read_nodes,
      ReadList& reads, int32_t read_idx,
      std::vector<std::vector<int32_t>>& read_group,
      std::vector<int32_t>& merge_map,
      PerStrand<std::vector<int32_t>>& group_to_pre_graph_segment,
      PerStrand<SegmentNodesList>& segment_nodes,
      PerStrand<std::vector<SpliceGraph>>& splice_graphs);

  CTransfrag *add_read_abundance_to_transfrag(
      SpliceGraph& graph, std::vector<int32_t>& nodes, BitVec& path_bits,
      float abundance, float raw_scale, uint32_t read_start, uint32_t read_end,
      int32_t cell_count, int32_t cell, uint32_t gap_start, uint32_t gap_end);

  void assign_fragment_to_transfrags(BundleData *bundle,
                                     AssemblyState& state, GraphTables& tables,
                                     int32_t read_idx, int32_t mate_idx,
                                     float read_coverage, float raw_scale);

  TransfragIndex build_transfrag_index(SpliceGraph& graph);
} // namespace ribbon
