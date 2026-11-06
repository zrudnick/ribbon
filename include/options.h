#pragma once

#include "ribbon.h"

// Command-line parsing: process_options reads the options into the global
// parameters of the run.

#define VERSION "3.0.3"

namespace ribbon {
  void process_options(int32_t argc, char *argv[],
                       std::vector<std::string>& in_files);
} // namespace ribbon
