#pragma once

#include "ribbon.h"

// GTF output for a finished bundle: turns the kept predictions into printed
// gene/transcript records (print_results) and finishes the output file
// (write_final_output); see out.cpp.

namespace ribbon {
  void write_final_output();

  int32_t print_results(BundleData *bundle, int32_t gene_count,
                    const std::string& ref_name);
} // namespace ribbon
