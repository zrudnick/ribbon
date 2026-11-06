#include "ribbon.h"
#include "GThreads.h"

// Marks the nodes and edges of a graph whose coverage is too low to follow.

namespace ribbon {

  struct CPred {
    int32_t n;
    float cov;
    CPred(int32_t n = 0, float c = 0) : n(n), cov(c) {}
  };

  // Comparator for CPred entries (n = node id, cov = coverage): most
  // covered node first.
  int32_t compare_node_by_coverage(const void *p1, const void *p2) {
    CPred *a = (CPred *)p1;
    CPred *b = (CPred *)p2;
    if (a->cov < b->cov) return 1; // most abundant comes first
    if (a->cov > b->cov) return -1;
    return 0;
  }

  // Width of the region on each side of a boundary whose coverage we compare.
  const int32_t adjacent_region_len = 1;

  // Coverage per bp and per cell under which a node that touches no other node
  // is flagged as low coverage. The node coverage is pooled over the cells that
  // have reads in the bundle, so the cutoff is this value times their number.
  // On the 192 simulated cells this gave a small gain in precision for a
  // smaller loss in sensitivity (about 1.5 reads per bp over all 192 cells).
  const double lonely_node_min_coverage_per_cell = 0.0078;

  // Looks at the edges that join node i to its parents (at its start) or to
  // its children (at its end). If the coverage outside the boundary is high
  // next to the coverage inside, the edges to the nodes that are not next to i
  // are flagged as low coverage (set in low_cov_flags). Otherwise returns the
  // edges to the nodes other than the source or sink, to be compared with each
  // other by flag_weak_edges.
  std::vector<CPred> flag_or_collect_edges(const FlowContext& context,
                                           int32_t i, bool parents,
                                           BitVec& low_cov_flags) {
    int32_t ref_start = context.bundle->start;
    BpCoverage& bp_coverage = context.bundle->bp_coverage;
    CGraphnode *node = context.graph.nodes[i].get();
    std::vector<int32_t>& neighbors = parents ? node->parent : node->child;

    std::vector<CPred> edges;
    if (count_of(neighbors) <= 1) return edges;

    float outside_coverage;
    float inside_coverage;
    if (parents) {
      outside_coverage = get_stranded_coverage(context.strand,
          node->start - adjacent_region_len - ref_start,
          node->start - 1 - ref_start, bp_coverage);
      inside_coverage = get_stranded_coverage(context.strand,
          node->start - ref_start,
          node->start + adjacent_region_len - 1 - ref_start, bp_coverage);
    } else {
      inside_coverage = get_stranded_coverage(context.strand,
          node->end - adjacent_region_len + 1 - ref_start,
          node->end - ref_start, bp_coverage);
      outside_coverage = get_stranded_coverage(context.strand,
          node->end + 1 - ref_start,
          node->end + adjacent_region_len - ref_start, bp_coverage);
    }

    if (outside_coverage > (1 - ERROR_PERC) * inside_coverage) {
      // For all neighbors that are not the source (or sink) and are not
      // adjacent -> the connection to the source should not be deleted
      for (int32_t neighbor : neighbors) {
        if (parents ? (neighbor && neighbor < i - 1)
                    : (neighbor < context.graph.node_count - 1 && neighbor > i
                        + 1)) {
          context.graph.path_bit_layout.set_edge_bit(low_cov_flags, parents
              ? neighbor : i,
                                      parents ? i : neighbor);
        }
      }
    } else {
      // For all neighbors that are not the source (or sink)
      for (int32_t neighbor : neighbors) {
        if (parents ? neighbor != 0 : neighbor < context.graph.node_count - 1) {
          int32_t *edge_pos = parents
              ? context.graph.path_bit_layout.find_edge_bit(neighbor, i)
              : context.graph.path_bit_layout.find_edge_bit(i, neighbor);
          if (edge_pos) edges.push_back(CPred(*edge_pos, 0));
        }
      }
    }
    return edges;
  }

  // Adds up the abundance of the transfrags of node i that cross each of the
  // edges (to the parents if left_side, else to the children) and flags the
  // edges that carry less than low_iso_frac of the strongest one. Does nothing
  // for less than two edges.
  void flag_weak_edges(const FlowContext& context, int32_t i,
                       std::vector<CPred>& edges, bool left_side,
                       BitVec& low_cov_flags) {
    if (count_of(edges) <= 1) return;

    CGraphnode *node = context.graph.nodes[i].get();
    float max_abundance = 0;
    for (int32_t j = 0; j < count_of(node->transfrag_ids); j++) {
      int32_t t = node->transfrag_ids[j];

      // The transfrag starts before node i (left) or ends after it (right)
      if (left_side ? context.graph.transfrag[t]->nodes[0] < i
                    : context.graph.transfrag[t]->nodes.back() > i) {
        for (CPred& edge_entry : edges) {
          // Transfrag connects the nodes with an actual junction
          if (context.graph.transfrag[t]->path_bits[edge_entry.n]) {
            edge_entry.cov += context.graph.transfrag[t]->abundance;
            if (edge_entry.cov > max_abundance) max_abundance = edge_entry.cov;
            break; // transfrag can only connect two nodes by an edge
          }
        }
      }
    }

    // If there is a maximal junction -> check junctions for lowcov
    if (max_abundance) {
      for (CPred& edge_entry : edges) {
        if (edge_entry.cov < low_iso_frac * max_abundance) {
          // Assume a 1% error for splice junctions
          low_cov_flags[edge_entry.n] = 1;
        }
      }
    }
  }

  // Walks from node seed to the left (or right) over the nodes that touch each
  // other, as long as they have less coverage per bp than the seed, and flags
  // them as low coverage unless they have a good connection elsewhere. Does
  // nothing if the seed has only one parent (or child), besides the source (or
  // sink).
  void flag_low_coverage_adjacent_nodes(const FlowContext& context,
                                        int32_t seed, float seed_coverage,
                                        bool leftwards,
                                        BitVec& low_cov_flags) {
    int32_t ref_start = context.bundle->start;
    BpCoverage& bp_coverage = context.bundle->bp_coverage;
    GraphNodeList& nodes = context.graph.nodes;
    int32_t step = leftwards ? -1 : 1;
    // The source when going left, the sink when going right
    int32_t terminal = leftwards ? 0 : context.graph.node_count - 1;

    // Check for low coverages of node only if there is more than one
    // parent/child
    std::vector<int32_t>& seed_neighbors
        = leftwards ? nodes[seed]->parent : nodes[seed]->child;
    if (count_of(seed_neighbors) == 1
        || (count_of(seed_neighbors) == 2
        && (seed_neighbors[0] == terminal || seed_neighbors[1] == terminal)))
      return;

    // Is there a good (not low coverage) edge between these nodes?
    auto is_good_edge = [&](int32_t from, int32_t to) {
      int32_t *edge_pos = context.graph.path_bit_layout.find_edge_bit(from, to);
      return edge_pos && !low_cov_flags[*edge_pos];
    };

    int32_t n = seed;
    while (true) {
      int32_t adjacent = n + step;
      if (leftwards ? adjacent <= 0 : adjacent >= context.graph.node_count - 1)
        break;
      CGraphnode *lower = nodes[leftwards ? adjacent : n].get();
      CGraphnode *upper = nodes[leftwards ? n : adjacent].get();
      CGraphnode *adjacent_node = nodes[adjacent].get();
      bool nodes_touch = lower->end == upper->start - 1;
      bool has_less_coverage
          = seed_coverage > adjacent_node->cov / adjacent_node->len();
      if (!nodes_touch || !has_less_coverage) break;

      // First check if nearby region is above threshold too --> maybe big node
      // gets above threshold but nearby it is not
      float left_coverage = get_stranded_coverage(context.strand,
          upper->start - adjacent_region_len - ref_start,
          lower->end - ref_start, bp_coverage);
      float right_coverage = get_stranded_coverage(context.strand,
          upper->start - ref_start,
          lower->end + adjacent_region_len - ref_start, bp_coverage);
      // V20 <-- adjacent_region_len
      if (leftwards ? left_coverage > SCELL_TOL * right_coverage
                    : right_coverage > SCELL_TOL * left_coverage) {
        break; // found adjacent node over threshold
      }

      // The adjacent node is low coverage unless it has a good connection
      // elsewhere: to another node on the side that faces n, or to a node on
      // the far side other than the one next to it and the source (or sink).
      bool mark_low = true;
      std::vector<int32_t>& facing
          = leftwards ? adjacent_node->child : adjacent_node->parent;
      std::vector<int32_t>& far_side
          = leftwards ? adjacent_node->parent : adjacent_node->child;
      if (count_of(facing) > 1) {
        for (int32_t neighbor : facing)
          if (neighbor != n
              && (leftwards ? is_good_edge(adjacent, neighbor)
                            : is_good_edge(neighbor, adjacent))) {
            mark_low = false;
            break;
          }
      }
      for (int32_t neighbor : far_side)
        if (neighbor != terminal && neighbor != adjacent + step
            && (leftwards ? is_good_edge(neighbor, adjacent)
                          : is_good_edge(adjacent, neighbor))) {
          mark_low = false;
          break;
        }

      if (mark_low) {
        low_cov_flags[adjacent] = 1;
      }
      n = adjacent;
    }
  }

  // Flags (in low_cov_flags, by node id and edge bit position) the nodes and
  // edges whose pooled coverage is negligible next to their neighbors, so that
  // the per-cell path extension can avoid them. Computed once per graph.
  void mark_low_coverage_nodes_and_edges(const FlowContext& context,
                                         BitVec& low_cov_flags) {
    // The per-bp coverage of every node (n = node id, cov = coverage).
    std::vector<CPred> node_coverages;
    for (int32_t i = 1; i < context.graph.node_count - 1; i++) {
      // First mark junctions that are less than the error noise threshold
      std::vector<CPred> left_edges
          = flag_or_collect_edges(context, i, true, low_cov_flags);
      std::vector<CPred> right_edges
          = flag_or_collect_edges(context, i, false, low_cov_flags);

      // Only if more than one junction on either side leave the node do we
      // check for maximal junction
      flag_weak_edges(context, i, left_edges, true, low_cov_flags);
      flag_weak_edges(context, i, right_edges, false, low_cov_flags);

      // Remember the node and its coverage per bp.
      CPred p(i, context.graph.nodes[i]->cov / context.graph.nodes[i]->len());
      node_coverages.push_back(p);
    }
    // Sort the nodes, most covered first, so that the weaker nodes next to a
    // node are judged against the strongest ones first.
    sort_by(node_coverages, compare_node_by_coverage);

    // Flag the less covered nodes next to each node, unless they are well
    // connected elsewhere.
    for (int32_t i = 0; i < count_of(node_coverages); i++) {
      // This is not a low coverage node --> set lowcov nodes and junctions
      if (!low_cov_flags[node_coverages[i].n]) {
        flag_low_coverage_adjacent_nodes(context, node_coverages[i].n,
            node_coverages[i].cov, true, low_cov_flags);
        flag_low_coverage_adjacent_nodes(context, node_coverages[i].n,
            node_coverages[i].cov, false, low_cov_flags);
      }
    }

    // Flag the lonely nodes: a node that touches no other node is never
    // compared with a neighbor above, so we judge it by its own coverage.
    GraphNodeList& nodes = context.graph.nodes;
    double lonely_node_min_coverage = lonely_node_min_coverage_per_cell
                                      * count_of(context.bundle->cell_name);
    for (int32_t i = 1; i < context.graph.node_count - 1; i++) {
      bool touches_previous = i > 1 && nodes[i - 1]->end + 1 == nodes[i]->start;
      bool touches_next = i + 1 < context.graph.node_count - 1
                          && nodes[i]->end + 1 == nodes[i + 1]->start;
      if (!touches_previous && !touches_next && !low_cov_flags[i]
          && nodes[i]->cov / nodes[i]->len() < lonely_node_min_coverage) {
        low_cov_flags[i] = 1;
      }
    }
  }

} // namespace ribbon
