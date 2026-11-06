#include "ribbon.h"
#include "GThreads.h"

// Which bit of a transfrag's path_bits stands for which node or edge of a graph.

namespace ribbon {

  // Maps an edge (from < to, node_count nodes with the source) to a unique
  // integer (triangular-matrix packing): the key into bit_of_edge.
  int32_t PathBitLayout::edge_key(int32_t from, int32_t to) const {
    // This includes source to node edges
    return ((node_count - 1) * (from + 1) - from * (from - 1) / 2
            + to - from);
  }

  // Starts the layout of a graph with node_count nodes: the nodes take the
  // first node_count bits, so the first free bit for an edge is node_count.
  void PathBitLayout::init(int32_t node_count_in) {
    node_count = node_count_in;
    next_edge_bit = node_count_in;
  }

  // Bit of the edge from -> to in the path bit vectors, or NULL if the edge
  // has no bit yet.
  int32_t *PathBitLayout::find_edge_bit(int32_t from, int32_t to) {
    auto it = bit_of_edge.find(edge_key(from, to));
    return it == bit_of_edge.end() ? NULL : &it->second;
  }

  // True if the edge from -> to has a bit and it is set in the path bit
  // vector.
  bool PathBitLayout::has_edge_bit(const BitVec& path_bits, int32_t from,
                                   int32_t to) {
    int32_t *edge_bit = find_edge_bit(from, to);
    return edge_bit && path_bits[*edge_bit];
  }

  // Bit of the edge from -> to, taking the next free bit (and advancing it)
  // if the edge has none yet.
  int32_t PathBitLayout::get_or_create_edge_bit(int32_t from, int32_t to) {
    int32_t *edge_bit = find_edge_bit(from, to);
    if (edge_bit) return *edge_bit;
    bit_of_edge[edge_key(from, to)] = next_edge_bit;
    return next_edge_bit++;
  }

  // Sets (or clears) the bit of the edge from -> to in a path bit vector. Does
  // nothing if the edge has no bit yet.
  void PathBitLayout::set_edge_bit(BitVec& path_bits, int32_t from, int32_t to,
                                   bool value) {
    int32_t *edge_bit = find_edge_bit(from, to);
    if (edge_bit) path_bits[*edge_bit] = value;
  }

} // namespace ribbon
