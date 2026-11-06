#pragma once

#include "ribbon.h"

// Which bit of a transfrag's path_bits stands for which node or edge of a
// graph.

namespace ribbon {
  class PathBitLayout {
    int32_t node_count = 0;
    int32_t next_edge_bit = 0;
    std::unordered_map<int32_t, int32_t> bit_of_edge;

    int32_t edge_key(int32_t from, int32_t to) const;

  public:
    void init(int32_t node_count_in);
    int32_t *find_edge_bit(int32_t from, int32_t to);
    bool has_edge_bit(const BitVec& path_bits, int32_t from, int32_t to);
    int32_t get_or_create_edge_bit(int32_t from, int32_t to);
    void set_edge_bit(BitVec& path_bits, int32_t from, int32_t to,
                      bool value = true);
  };
} // namespace ribbon
