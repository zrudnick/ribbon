#pragma once

#include "ribbon.h"

// Junction support: counts the reads that support each junction, rescues or
// redirects weak junctions and filters the ones that stay unsupported.

namespace ribbon {
  void compute_junction_support(BundleData *bundle);

  void propagate_junction_support(BundleData *bundle,
                                  AssemblyState& state);

  // Looks around index idx of a junction list sorted by start (or by end if
  // by_end) for the junction with this start, end and strand: first backwards
  // from idx, then forwards, among the junctions with the same sorted
  // coordinate. Returns NULL if there is none.
  template <class JunctionPtrList>
  CJunction *find_junction_around(JunctionPtrList& list, int32_t idx,
                                  bool by_end, uint32_t start, uint32_t end,
                                  char strand) {
    auto has_same_coordinate = [&](const CJunction *junction) {
      return by_end ? junction->end == end : junction->start == start;
    };
    auto is_wanted = [&](const CJunction *junction) {
      return junction->strand == strand && junction->start == start
          && junction->end == end;
    };
    for (int32_t k = idx; k > 0 && has_same_coordinate(&*list[k]); k--) {
      if (is_wanted(&*list[k])) return &*list[k];
    }
    for (int32_t k = idx + 1;
         k < count_of(list) && has_same_coordinate(&*list[k]); k++) {
      if (is_wanted(&*list[k])) return &*list[k];
    }
    return NULL;
  }

  CJunction *find_added_junction(JunctionList& junction,
                                 int32_t original_count, uint32_t start,
                                 uint32_t end, char strand);

  bool is_junction_well_supported(CJunction& junction, int32_t ref_start,
                                  BpCoverage& bp_coverage);

  int32_t compare_junction_by_end(const void *p1, const void *p2);
} // namespace ribbon
