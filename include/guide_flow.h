#pragma once

#include "ribbon.h"

// Flow for reference (guide) transcripts: guides claim their evidence first,
// their ends can move to better supported nodes, and their cells are credited.

namespace ribbon {
  void update_guide_prediction(const FlowContext& context,
                         PredictionList& predictions,
                         int32_t p, std::vector<int32_t>& path,
                         std::vector<float>& node_flux,
                         std::vector<float>& node_coverage);

  float push_guide_flow_and_deplete_transfrags(const FlowContext& context,
      int32_t g, std::vector<CGuide>& matched_guides, BitVec& transfrag_on_path,
      std::vector<float>& node_flux);

  float try_alternative_guide_ends(const FlowContext& context, int32_t g,
                                   std::vector<CGuide>& matched_guides,
                                   BitVec& transfrag_on_path,
                                   std::vector<float>& node_flux);

  float estimate_guide_flow(const FlowContext& context,
                            std::vector<int32_t>& path,
                            const GuideFlowCache& cache);

  int32_t compare_guide_by_abundance(const void *p1, const void *p2);

  void assign_flow_to_guides(const FlowContext& context,
                             std::vector<CGuide>& matched_guides,
                             int32_t& gene_count, PredictionList& predictions,
                             std::vector<float>& node_coverage,
                             BitVec& transfrag_on_path,
                             BitVec& current_path_bits,
                             bool& is_first_transcript,
                             std::vector<GffObj *>& guides,
                             std::vector<int32_t>& guide_prediction,
                             std::vector<int32_t>& printed_guides);

  int32_t assign_flow_to_guides_and_pick_seed(const FlowContext& context,
                                          std::vector<CGuide>& matched_guides,
                                          int32_t& gene_count,
                                              PredictionList& predictions,
                                          std::vector<float>& node_coverage,
                                          BitVec& transfrag_on_path,
                                          BitVec& current_path_bits,
                                          bool& is_first_transcript,
                                          std::vector<GffObj *>& guides,
                                          std::vector<int32_t>& guide_prediction);
} // namespace ribbon
