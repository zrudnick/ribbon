#include "ribbon.h"
#include "GThreads.h"

// Weight of the per-cell abundance in the path-selection score
extern double alpha;
// Weight of the pooled (all-cells) abundance in the path-selection score
extern double beta;
// --min-borrowed-share: see PathExtender::check_weak_borrow (0: off)
extern float min_borrowed_share;

// Greedy extension of a seed node into a source-to-sink path, choosing at
// each step the best supported neighbor.

namespace ribbon {

  // Returns whether a transfrag fits the path segment from first_node to
  // last_node: inside it the transfrag nodes and edges must be on the path
  // (current_path_bits), and outside it the nodes must be reachable from the
  // path.
  bool is_transfrag_on_path(const FlowContext& context,
                            BitVec& transfrag_path_bits,
                            std::vector<int32_t>& transfrag_nodes,
                            BitVec& current_path_bits, int32_t first_node,
                            int32_t last_node) {
    // Graph and cell data shared through the context.

    // Check that the parts of the transfrag outside the segment connect to it.
    // Mini can be reached through transcript
    if (transfrag_nodes[0] < first_node)
      if (!context.graph.nodes[first_node]
               ->parent_reachable_bits[transfrag_nodes[0]])
          return false;

    // From last_node we can reach the end of the transcript
    if (transfrag_nodes.back() > last_node)
      if (!context.graph.nodes[last_node]
               ->child_reachable_bits[transfrag_nodes.back()])
        return false;

    // Check the nodes and edges inside the segment against the path.
    int32_t is_first_in_range = 1;

    for (int32_t i = 0; i < count_of(transfrag_nodes); i++) {
      if (transfrag_nodes[i] >= first_node && transfrag_nodes[i] <= last_node) {
        if (!current_path_bits[transfrag_nodes[i]]) return false;
        int32_t *edge_pos = NULL;
        if (i)
          edge_pos = context.graph.path_bit_layout.find_edge_bit(
              transfrag_nodes[i - 1], transfrag_nodes[i]);
        if (is_first_in_range) {
          is_first_in_range = 0;
          if (i && transfrag_nodes[i] > first_node && edge_pos
              && transfrag_path_bits[*edge_pos])
            return false;
          // We can not reach mini from previous node in transfrag
          if (i
              && !context.graph.nodes[first_node]
                      ->parent_reachable_bits[transfrag_nodes[i - 1]])
            return false;
        } else if (i && edge_pos && transfrag_path_bits[*edge_pos]
                   && !current_path_bits[*edge_pos]) return false;
      }
      if (transfrag_nodes[i] > last_node) {
        int32_t *edge_pos = context.graph.path_bit_layout.find_edge_bit(
            transfrag_nodes[i - 1], transfrag_nodes[i]);
        if (i && transfrag_nodes[i - 1] < last_node && edge_pos
            && transfrag_path_bits[*edge_pos])
          return false;
        // We can not reach this node from last_node
        if (!context.graph.nodes[last_node]
                 ->child_reachable_bits[transfrag_nodes[i]])
            return false;
        break;
      }
    }

    return true;
  }

  // The bit of the edge between two nodes of a graph that extend_path walks:
  // toward is the node closer to the end of the walk (the source going toward
  // the source, else the sink) and away the other one.
  int32_t *find_walk_edge_bit(const FlowContext& context, bool toward_source,
                              int32_t toward, int32_t away) {
    return toward_source
        ? context.graph.path_bit_layout.find_edge_bit(toward, away)
                         : context.graph.path_bit_layout.find_edge_bit(away,
                             toward);
  }

  // True if node a has more coverage per bp than node b (compared without
  // dividing).
  bool has_higher_coverage(const FlowContext& context, int32_t a, int32_t b) {
    return context.graph.nodes[a]->cov * context.graph.nodes[b]->len()
        > context.graph.nodes[b]->cov * context.graph.nodes[a]->len();
  }

  // Walks from the node adjacent to node i over the nodes that touch each
  // other, by topology alone, toward the source (toward_source) or the sink.
  // If it gets there, adds the nodes it passed to the path and returns true.
  bool follow_touching_nodes_to_end(const FlowContext& context, int32_t i,
                                    int32_t adjacent, bool toward_source,
                                    std::vector<int32_t>& path,
                                    BitVec& current_path_bits) {
    int32_t end_node = toward_source ? 0 : context.graph.node_count - 1;
    std::vector<int32_t> walk_to_end;
    walk_to_end.push_back(i);
    int32_t node = adjacent;
    while (toward_source ? node > 0 : node < context.graph.node_count - 1) {
      walk_to_end.push_back(node);
      int32_t next_node = toward_source ? -1 : context.graph.node_count;
      std::vector<int32_t>& next_candidates = toward_source
          ? context.graph.nodes[node]->parent
              : context.graph.nodes[node]->child;
      for (int32_t j = 0; j < count_of(next_candidates); j++) {
        int32_t candidate = next_candidates[j];
        if (candidate == end_node) { // candidate is the source (or sink)
          next_node = end_node;
          break;
        }
        bool touches = toward_source
            ? context.graph.nodes[candidate]->end
                == context.graph.nodes[node]->start - 1
            : context.graph.nodes[candidate]->start
                == context.graph.nodes[node]->end + 1;
        if (touches) next_node = candidate;
      }
      node = next_node;
    }
    if (node != end_node) return false;

    // Reached the source (or sink)
    for (int32_t j = 1; j < count_of(walk_to_end); j++) {
      path.push_back(walk_to_end[j]);
      current_path_bits[walk_to_end[j]] = 1;
      int32_t *edge_pos = find_walk_edge_bit(context, toward_source,
                                             walk_to_end[j],
                                             walk_to_end[j - 1]);
      if (edge_pos) current_path_bits[*edge_pos] = 1;
    }
    return true;
  }

  int32_t *PathExtender::find_edge_bit(int32_t toward, int32_t away) {
    return find_walk_edge_bit(context, toward_source, toward, away);
  }

  void PathExtender::report_missing_edge(int32_t toward, int32_t source_code,
                           int32_t sink_code) {
    fprintf(stderr, "%d Found parent-child %d-%d not linked by edge\n",
            toward_source ? source_code : sink_code,
            toward_source ? toward : i, toward_source ? i : toward);
    exit(1);
  }

  // Sets or clears the bits of a node and of its edge to i in the path
  // path_bits.
  void PathExtender::set_path_bits(int32_t node, int32_t *edge_pos,
                                   bool value) {
    current_path_bits[node] = value;
    if (edge_pos) current_path_bits[*edge_pos] = value;
  }

  // Adds the support for the path through the neighbor node of i that the
  // transfrags going through it give: with this cell's abundance, or only
  // the pooled abundance if the cell has none.
  // The transfrags with abundance in this cell give support, the others
  // pooled_support. A caller that needs only one of them says so, and the
  // transfrags of the other kind are not checked.
  void PathExtender::add_transfrag_support(int32_t neighbor, float& support,
                             float& pooled_support, bool want_support,
                             bool want_pooled_support) {
    CGraphnode *neighbor_node = context.graph.nodes[neighbor].get();
    // The pooled support comes from the transfrags that the cell lacks
    const std::vector<int32_t>& node_transfrags = want_pooled_support
        ? neighbor_node->transfrag_ids : transfrags_of_cell(context, neighbor);
    for (int32_t j = 0; j < count_of(node_transfrags); j++) {
      int32_t t = node_transfrags[j];
      bool in_cell = !(context.graph.transfrag[t]->cell_abundance[context.cell]
                       < epsilon);
      if (in_cell ? !want_support : !want_pooled_support) continue;
      // The transfrag goes between i and the neighbor and is compatible with
      // the path.
      int32_t left_node = toward_source ? neighbor : i;
      int32_t right_node = toward_source ? i : neighbor;
      bool spans_i_and_neighbor
          = context.graph.transfrag[t]->nodes[0] <= left_node
          && context.graph.transfrag[t]->nodes.back() >= right_node;
      if (spans_i_and_neighbor
          && is_transfrag_on_path(context,
              context.graph.transfrag[t]->path_bits,
                                  context.graph.transfrag[t]->nodes,
                                  current_path_bits,
                                  toward_source ? neighbor : path[0],
                                  toward_source ? path[0] : neighbor)) {
        if (!in_cell) {
          // If cell abundance is zero than use transfrag abundances to guide
          // most likely path
          pooled_support += context.graph.transfrag[t]->abundance;
        } else {
          // Not sure why this wouldn't be the case
          if (context.graph.transfrag[t]->abundance > epsilon) {
            // Take global abundances into account as well as they are more
            // informative
            support += exp(
                alpha * log(context.graph.transfrag[t]
                                ->cell_abundance[context.cell])
                + beta * log(context.graph.transfrag[t]->abundance));
          } else {
            support += context.graph.transfrag[t]->cell_abundance[context.cell];
          }
          path_empty = false;
        }
      }
    }
  }

  // True if the neighbor is the next node, touching the current node, and
  // the coverage drops from the current node to it. The pooled coverage drops
  // from the current node to the neighbor. The drop counts for this cell only
  // if the cell sees it too (or has no reads on the current node): in the
  // pool, a longer isoform that one cell type makes looks like a drop.
  bool PathExtender::is_excluded_adjacent_neighbor(int32_t neighbor,
                                     CGraphnode *current_node,
                                     CGraphnode *neighbor_node,
                                     bool adjacent_to_neighbor) {
    bool pooled_drop = neighbor_node->cov / neighbor_node->len() < 1000
        && current_node->cov * DROP * neighbor_node->len()
        > neighbor_node->cov * current_node->len();
    float cell_current_coverage = current_node->cell_cov[context.cell];
    bool cell_drop = cell_current_coverage <= 0
        || cell_current_coverage * DROP * neighbor_node->len()
        > neighbor_node->cell_cov[context.cell] * current_node->len();

    bool inside_graph = toward_source ? i > 1 : i < context.graph.node_count
        - 2;
    bool is_next_touching_node = neighbor == i + step && inside_graph
        && adjacent_to_neighbor;
    return is_next_touching_node
           && (low_cov_flags[i + step] || (pooled_drop && cell_drop));
  }

  // Scores each neighbor by the support of the compatible transfrags that
  // cross between it and i. The pooled support is only used when no neighbor
  // has support in this cell, so it is computed after this loop.
  void PathExtender::score_neighbors_by_cell_support(CGraphnode *current_node,
                                       std::vector<int32_t>& neighbors) {
    for (int32_t n = 0; n < count_of(neighbors); n++) {
      int32_t neighbor = neighbors[n];
      CGraphnode *neighbor_node = context.graph.nodes[neighbor].get();

      // Do not consider neighbors that are lowcoverage in universal splice
      // graph unless they are the only option v18 (or have good abundance
      // locally, A3)
      bool gap_to_neighbor = toward_source
          ? neighbor_node->end < current_node->start - 1
          : current_node->end < neighbor_node->start - 1;
      if (low_cov_flags[neighbor] && gap_to_neighbor) continue;

      int32_t *edge_pos = find_edge_bit(neighbor, i);
      // Do not follow low coverage edges
      if (edge_pos && low_cov_flags[*edge_pos]) continue;

      float neighbor_support = 0;
      float unused_pooled_support = 0;

      // Adjacent to the neighbor
      // v18
      bool adjacent_to_neighbor = toward_source
          ? current_node->start == neighbor_node->end + 1
          : current_node->end + 1 == neighbor_node->start;

      if (is_excluded_adjacent_neighbor(neighbor, current_node,
                                        neighbor_node, adjacent_to_neighbor)) {
        adjacent_excluded = true;
      } else {
        set_path_bits(neighbor, edge_pos, 1);
        // This shouldn't happen (we have more tolerance for a link to sink)
        if (!edge_pos && (toward_source || neighbor < context.graph.node_count
            - 1))
          report_missing_edge(neighbor, 4, 1);

        add_transfrag_support(neighbor, neighbor_support,
                              unused_pooled_support, true, false);
        scored_neighbors.push_back({neighbor, edge_pos});

        if (neighbor_support > best_support) {
          best_support = neighbor_support;
          best_neighbor = neighbor;
        } else if (neighbor_support == best_support && best_neighbor != -1) {
          // Prefer node with higher coverage in the global graph; has
          // problem with source (coverage 0 by default) and with preferring
          // sink
          if (has_higher_coverage(context, neighbor, best_neighbor)) {
            best_neighbor = neighbor;
          }
        }

        set_path_bits(neighbor, edge_pos, 0);
      }
    }
  }

  // No neighbor has support in this cell: scores the neighbors by their
  // pooled support.
  void PathExtender::score_neighbors_by_pooled_support() {
    for (const auto& scored : scored_neighbors) {
      int32_t neighbor = scored.first;
      int32_t *edge_pos = scored.second;
      float unused_support = 0;
      float neighbor_pooled_support = 0;

      set_path_bits(neighbor, edge_pos, 1);
      add_transfrag_support(neighbor, unused_support,
                            neighbor_pooled_support, false, true);
      set_path_bits(neighbor, edge_pos, 0);
      pooled_support_sum += neighbor_pooled_support;
      if (neighbor_pooled_support > 0) pooled_supported_count++;

      if (neighbor_pooled_support > best_pooled_support) {
        best_pooled_support = neighbor_pooled_support;
        best_pooled_neighbor = neighbor;
      // If both neighbors have the same empty coverage then choose the one
      // with the best node coverage
      } else if (neighbor_pooled_support == best_pooled_support
                 && best_pooled_neighbor != -1) {
        // Prefer node with higher coverage in the global graph
        if (has_higher_coverage(context, neighbor, best_pooled_neighbor)) {
          best_pooled_neighbor = neighbor;
        }
      }
    }
  }

  // Called when the best pooled neighbor is taken: flags the path if this is a
  // spliced step won with less than min_borrowed_share of the pooled support
  // of several supported neighbors (such paths are mostly wrong).
  void PathExtender::check_weak_borrow() {
    if (min_borrowed_share <= 0 || pooled_supported_count < 2) return;
    if (best_pooled_support >= min_borrowed_share * pooled_support_sum) return;

    // Only a splice between two real nodes counts, not a step to a touching
    // node or to the source (or sink).
    if (best_pooled_neighbor == end_node) return;
    CGraphnode *current_node = context.graph.nodes[i].get();
    CGraphnode *neighbor_node = context.graph.nodes[best_pooled_neighbor].get();
    bool touching = toward_source
        ? neighbor_node->end + 1 >= current_node->start
        : current_node->end + 1 >= neighbor_node->start;
    if (!touching) weak_borrow = true;
  }

  // No neighbor has cell support: falls back to the adjacent neighbor, then
  // to the best pooled neighbor, or follows the unspliced nodes to the source
  // (or sink).
  PathExtender::StepResult
  PathExtender::choose_neighbor_without_cell_support() {
    int32_t adjacent = i + step;
    // Assess the nearby node as well v18
    if (adjacent_excluded && !low_cov_flags[adjacent]) {
      float adjacent_support = 0;
      float adjacent_pooled_support = 0;
      int32_t *edge_pos = find_edge_bit(adjacent, i);
      set_path_bits(adjacent, edge_pos, 1);
      if (!edge_pos && (toward_source || i))
        report_missing_edge(adjacent, 5, 2);

      add_transfrag_support(adjacent, adjacent_support,
                            adjacent_pooled_support, true, true);
      set_path_bits(adjacent, edge_pos, 0);

      if (adjacent_support) {
        best_neighbor = adjacent;
      // Found max before (the same test in both directions; toward the sink
      // it used to require pooled support of the adjacent node too)
      } else if (best_pooled_support && best_pooled_neighbor != -1) {
        if (adjacent_pooled_support > best_pooled_support
            && node_coverage[adjacent]
            > node_coverage[best_pooled_neighbor]) best_neighbor = adjacent;
        else {
          best_neighbor = best_pooled_neighbor;
          check_weak_borrow();
        }
      } else if (adjacent_pooled_support) best_neighbor = adjacent;
      else {
        return StepResult::NoWay;
      }
    } else if (best_pooled_neighbor != -1) {
      best_neighbor = best_pooled_neighbor;
      check_weak_borrow();
    } else if (adjacent_excluded && low_cov_flags[adjacent]) { // v18
      // Walk the adjacent nodes to the source (or sink) by topology alone.
      return follow_touching_nodes_to_end(context, i, adjacent, toward_source,
                                          path, current_path_bits)
                 ? StepResult::ReachedEnd : StepResult::NoWay;
    } else {
      return StepResult::NoWay;
    }
    return StepResult::Continue;
  }

  // Adds the chosen neighbor to the path and moves to it. The source is not
  // added to the path, only to the path's path_bits.
  void PathExtender::move_to_best_neighbor() {
    if (!toward_source || best_neighbor) {
      path.push_back(best_neighbor);
    }

    int32_t *edge_pos = find_edge_bit(best_neighbor, i);
    set_path_bits(best_neighbor, edge_pos, 1);
    // Allow a missing source->i (or i->sink) edge mapping: some graphs might
    // not have the explicit edge stored in the layout.
    if (!edge_pos && best_neighbor != end_node)
      report_missing_edge(best_neighbor, 6, 3);

    i = best_neighbor;
  }

  // Finds the neighbor to go to from node i and moves there.
  PathExtender::StepResult PathExtender::walk_one_node() {
    // Find all neighbors -> if there are none then this node is the source
    // (or the sink)
    CGraphnode *current_node = context.graph.nodes[i].get();
    std::vector<int32_t>& neighbors
        = toward_source ? current_node->parent : current_node->child;
    if (!count_of(neighbors)) return StepResult::ReachedEnd;

    best_neighbor = -1;
    best_pooled_neighbor = -1;
    best_support = 0;
    best_pooled_support = 0;
    pooled_support_sum = 0;
    pooled_supported_count = 0;
    adjacent_excluded = false;
    scored_neighbors.clear();

    score_neighbors_by_cell_support(current_node, neighbors);

    if (best_neighbor == -1) score_neighbors_by_pooled_support();

    if (best_neighbor == -1) {
      StepResult result = choose_neighbor_without_cell_support();
      if (result != StepResult::Continue) return result;
    }

    move_to_best_neighbor();
    return StepResult::Continue;
  }

  PathExtender::PathExtender(const FlowContext& context_in, int32_t i_in,
               bool toward_source_in, std::vector<int32_t>& path_in,
               BitVec& current_path_bits_in,
               std::vector<float>& node_coverage_in,
               bool& path_empty_in, bool& weak_borrow_in,
               BitVec& low_cov_flags_in)
      : context(context_in), i(i_in), toward_source(toward_source_in),
        path(path_in), current_path_bits(current_path_bits_in),
        node_coverage(node_coverage_in), path_empty(path_empty_in),
        weak_borrow(weak_borrow_in), low_cov_flags(low_cov_flags_in),
        end_node(toward_source_in ? 0 : context_in.graph.node_count - 1),
        step(toward_source_in ? -1 : 1) {}

  // Walks to the end one node at a time. Returns false if there is no way.
  bool PathExtender::walk() {
    while (true) {
      StepResult result = walk_one_node();
      if (result != StepResult::Continue)
        return result == StepResult::ReachedEnd;
    }
  }

  // Extends path (ending at node i) toward the source (toward_source) or the
  // sink; returns false if there is no way. Sets weak_borrow (never clears it)
  // if a spliced step was a weak borrowing.
  bool extend_path(const FlowContext& context, int32_t i, bool toward_source,
                   std::vector<int32_t>& path, BitVec& current_path_bits,
                   std::vector<float>& node_coverage, bool& path_empty,
                   bool& weak_borrow, BitVec& low_cov_flags) {
    PathExtender extender(context, i, toward_source, path, current_path_bits,
                          node_coverage, path_empty, weak_borrow,
                          low_cov_flags);
    return extender.walk();
  }

} // namespace ribbon
