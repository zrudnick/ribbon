#pragma once

#include "ribbon.h"

// Marks the nodes and edges of a graph whose coverage is too low to follow.

namespace ribbon {
  // Marks nodes/edges whose pooled coverage is negligible (< ERROR_PERC of
  // their neighbors) in low_cov_flags, so the per-cell traversal functions
  // know which ones to route around.
  void mark_low_coverage_nodes_and_edges(const FlowContext& context,
                                         BitVec& low_cov_flags);
} // namespace ribbon
