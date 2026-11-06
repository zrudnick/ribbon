#pragma once

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>

// General helpers used by many files. They do not depend on the ribbon types.

// Number of elements of a container as an int. (size() is unsigned, so
// `v.size() - 1` on an empty vector would wrap around.)
template <class Container>
inline int32_t count_of(const Container& c) {
  return static_cast<int32_t>(c.size());
}

// Index of the element equal to item (neither is less than the other) in a
// vector kept sorted by operator<, or -1 if there is none.
template <class T>
inline int32_t find_sorted(const std::vector<T>& v, const T& item) {
  auto it = std::lower_bound(v.begin(), v.end(), item);
  if (it == v.end() || item < *it) return -1;
  return static_cast<int32_t>(it - v.begin());
}

// Inserts item at its place in a vector kept sorted by operator<.
template <class T>
inline void insert_sorted(std::vector<T>& v, const T& item) {
  v.insert(std::lower_bound(v.begin(), v.end(), item), item);
}

// Pointer to an element for the old three-way compare functions
// (int compare(const void*, const void*)), which take plain pointers.
template <class T>
inline T *compare_arg(T& x) { return &x; }
template <class T>
inline T *compare_arg(std::unique_ptr<T>& x) { return x.get(); }
template <class T>
inline T *compare_arg(T *&x) { return x; }

// Sorts v with a three-way compare function (negative: first goes before
// second). Elements the function considers equal keep their original order.
template <class T>
inline void sort_by(std::vector<T>& v,
                    int32_t (*compare)(const void *, const void *)) {
  std::stable_sort(v.begin(), v.end(), [compare](const T& a, const T& b) {
    return compare(compare_arg(const_cast<T&>(a)),
                   compare_arg(const_cast<T&>(b))) < 0;
  });
}

namespace ribbon {
  int32_t find_root(const std::vector<int32_t>& parent, int32_t x);

  bool erase_value(std::vector<int32_t>& v, int32_t value);

  void insert_sorted_unique(std::vector<int32_t>& v, int32_t value);
} // namespace ribbon
