#pragma once

#include "ribbon.h"

// Selection of the predictions of a finished bundle (finalize_predictions);
// see predictions.cpp.

namespace ribbon {
  void finalize_predictions(BundleData *bundle);

  // Transcript starts from 5' UMI molecules (--sq3; the mitochondrial
  // chromosome, or every chromosome with --tso-filter); see predictions.cpp.
  void find_start_peaks(BundleData *bundle);

  // Genes of the predictions (union-find over the flagged predictions that
  // share an exon on the same strand); see predictions.cpp.
  int32_t find_gene_root(std::vector<int32_t>& gene_parent, int32_t idx);
  std::vector<int32_t> group_predictions_into_genes(PredictionList& predictions);
} // namespace ribbon
