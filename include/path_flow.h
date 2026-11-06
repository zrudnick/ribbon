#pragma once

#include "ribbon.h"

// Pushing flow along a path of the splice graph: how much a path carries,
// depleting the transfrags that support it, and storing the path as a
// transcript prediction.

namespace ribbon {
  void apply_sparse_coverage_fallback(int32_t i,
      std::vector<float>& path_in_abundance,
      std::vector<float>& total_in_abundance,
      std::vector<float>& path_out_abundance,
      std::vector<float>& total_out_abundance);

  float compute_path_flow(int32_t path_len,
                                 std::vector<float>& path_in_abundance,
                                 std::vector<float>& total_in_abundance,
                                 std::vector<float>& path_out_abundance,
                                 std::vector<float>& total_out_abundance);

  bool deplete_transfrag_abundance(CTransfrag *t, int32_t cell,
                                          double amount,
                                          bool zero_if_pooled_depleted);

  const std::vector<int32_t>& transfrags_of_cell(
      const FlowContext& context, int32_t node);

  float push_flow_and_deplete_transfrags(const FlowContext& context,
                                         std::vector<int32_t>& path,
                                         BitVec& transfrag_on_path,
                                         std::vector<float>& node_flux,
                                         BitVec& current_path_bits);

  // True (and merges the coverage) when the transcript is the same as an
  // already stored prediction from the same graph, used by store_transcript to
  // avoid re-storing near-duplicate per-cell predictions.
  bool merge_if_equal_prediction(CPrediction *prediction, std::vector<GSeg>& exons,
                                 std::vector<float>& exon_cov);

  double store_transcript(const FlowContext& context, int32_t prev_prediction_count,
                          PredictionList& predictions, std::vector<int32_t>& path,
                          std::vector<float>& node_flux,
                          std::vector<float>& node_coverage,
                              int32_t& gene_count,
                          bool& is_first_transcript, bool& included,
                          BitVec& covered_path_bits,
                              BundleData *bundle = NULL,
                          GffObj *guide = NULL);
} // namespace ribbon
