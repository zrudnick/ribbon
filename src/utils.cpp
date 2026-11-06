#include "utils.h"

// General helpers used by many files. The templates are in utils.h.

namespace ribbon {

  // Returns the root of x in a union-find array (parent[x] == x for a root).
  int32_t find_root(const std::vector<int32_t>& parent, int32_t x) {
    while (parent[x] != x) x = parent[x];
    return x;
  }

  // Erases the first occurrence of value from v; false if there is none.
  bool erase_value(std::vector<int32_t>& v, int32_t value) {
    auto it = std::find(v.begin(), v.end(), value);
    if (it == v.end()) return false;
    v.erase(it);
    return true;
  }

  // Inserts value into the sorted vector v, unless it is already there.
  void insert_sorted_unique(std::vector<int32_t>& v, int32_t value) {
    auto it = std::lower_bound(v.begin(), v.end(), value);
    if (it == v.end() || *it != value) v.insert(it, value);
  }

} // namespace ribbon
