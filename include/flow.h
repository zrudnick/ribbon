#pragma once

#include "ribbon.h"

// Transcript extraction from the splice graph: greedy paths, each with a max
// flow through its transfrags; see flow.cpp.

namespace ribbon {
  // What the transcript extraction functions of one graph and one cell share:
  // the graph (its nodes, transfrags and path bit layout) and the cell
  // whose transfrag abundances we use.
  struct FlowContext {
    SpliceGraph& graph;
    char strand;
    BundleData *bundle;
    int32_t cell;
    // For each node, the transfrags that have abundance in this cell (in the
    // order of the node's transfrags); set once the cell's abundances are known
    const std::vector<std::vector<int32_t>> *cell_transfrags = nullptr;
    // True during the de novo pass that credits guides to cells
    // (extract_cell_denovo_predictions), where paths are not checked with
    // is_path_contradicted: there it only took guides away from their cells
    bool credit_pass = false;
    // For each node, the cell's abundance of transfrags that start there over
    // that of transfrags that end there (0 if either is 0); set once per cell,
    // see compute_cell_node_rates
    const std::vector<float> *node_rate = nullptr;
  };

  void extract_transcripts_greedily(const FlowContext& context,
                                    int32_t seed_node,
                                    int32_t& gene_count,
                                    bool is_first_transcript,
                                    PredictionList& predictions,
                                    std::vector<float>& node_coverage,
                                    BitVec& transfrag_on_path,
                                    BitVec& used_node, float prev_best_coverage,
                                    BitVec& covered_path_bits,
                                    BitVec& low_cov_flags,
                                    int32_t prev_prediction_count);

  void save_quant_paths(const FlowContext& context,
                        const std::vector<std::vector<float>>& original_cell_abundance,
                        int32_t first_prediction, double fragment_len);

  void drop_shadowed_paths(const FlowContext& context,
                           const std::vector<std::vector<float>>& original_cell_abundance,
                           int32_t first_prediction);

  int32_t find_transcripts(const FlowContext& shared_context,
                           int32_t gene_count,
                           std::vector<CGuide>& matched_guides,
                           std::vector<GffObj *>& guides,
                           std::vector<int32_t>& guide_prediction,
                           BitVec& low_cov_flags);

} // namespace ribbon
