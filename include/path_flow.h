#pragma once

#include "ribbon.h"

// The flow along a path of the splice graph: how much a path carries (a max
// flow through the cell's transfrags), depleting the transfrags that support
// it, and storing the path as a transcript prediction.

namespace ribbon {
  bool deplete_transfrag_abundance(CTransfrag *t, int32_t cell, double amount,
                                   bool zero_if_pooled_depleted);

  const std::vector<int32_t>& transfrags_of_cell(const FlowContext& context,
                                                 int32_t node);

  std::vector<float> compute_cell_node_rates(const FlowContext& context);

  float max_flow_and_deplete_transfrags(const FlowContext& context,
                                        std::vector<int32_t>& path,
                                        BitVec& transfrag_on_path,
                                        std::vector<float>& node_flux,
                                        BitVec& current_path_bits);

  // True (and merges the coverage) when the transcript is the same as an
  // already stored prediction from the same graph, used by store_transcript to
  // avoid re-storing near-duplicate per-cell predictions.
  bool merge_if_equal_prediction(CPrediction *prediction,
                                 std::vector<GSeg>& exons,
                                 std::vector<float>& exon_cov);

  double store_transcript(const FlowContext& context,
                          int32_t prev_prediction_count,
                          PredictionList& predictions,
                          std::vector<int32_t>& path,
                          std::vector<float>& node_flux,
                          std::vector<float>& node_coverage,
                          int32_t& gene_count,
                          bool& is_first_transcript, bool& included,
                          BitVec& covered_path_bits);
} // namespace ribbon
