#pragma once

#include "ribbon.h"

// Turns BAM records into reads (CReadAln): strand from the tags, exons and
// junctions from the CIGAR, junction checks against the MD tag and pairing of
// mates.

namespace ribbon {
  char get_strand(GSamRecord *brec);

  void process_read_in(int32_t curr_start, int32_t curr_end, BundleData& bundle,
                       std::unordered_map<std::string, int32_t>& hash_read,
                       GReadAlnData& aln_data);
} // namespace ribbon
