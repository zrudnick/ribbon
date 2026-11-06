#include "ribbon.h"
#include "GThreads.h"

extern int32_t allowed_nodes;

// Pruning of a splice graph that has too many nodes: removes the weakest
// junctions, merges the nodes left singly linked and renumbers the graph.

namespace ribbon {

  struct CGJunc {
    int32_t left_node;
    int32_t right_node;
    double cov; // ngood
    double good_cov; // ngood_reads
    CGJunc(int32_t n1 = 0, int32_t n2 = 0, double _cov = 0, double _good_cov = 0)
        : left_node(n1), right_node(n2), cov(_cov), good_cov(_good_cov) {}
  };

  // Comparator: orders candidate-for-removal CGJunc edges by their "good"
  // coverage then raw coverage, ascending, so the weakest edges sort first and
  // get pruned before stronger ones (see prune_graph_to_node_limit).
  int32_t compare_graph_junction_by_coverage(const void *p1,
                                             const void *p2) {
    CGJunc *a = (CGJunc *)p1;
    CGJunc *b = (CGJunc *)p2;

    if (a->good_cov > b->good_cov) return 1; // less coverage come first
    if (a->good_cov < b->good_cov) return -1;

    if (a->cov > b->cov) return 1; // less coverage come first
    if (a->cov < b->cov) return -1;

    return 0;
  }

  // True if node has the previous node as its only parent and no gap between
  // them, so it can be merged into the previous node.
  bool can_merge_into_previous(GraphNodeList& nodes, int32_t node) {
    return count_of(nodes[node]->parent) == 1
        && nodes[node]->parent[0] == node - 1
        && nodes[node - 1]->end == nodes[node]->start - 1;
  }

  // True if node has the next node as its only child and no gap between them,
  // so it can be merged into the next node.
  bool can_merge_into_next(GraphNodeList& nodes, int32_t node) {
    return count_of(nodes[node]->child) == 1
        && nodes[node]->child[0] == node + 1
        && nodes[node]->end == nodes[node + 1]->start - 1;
  }

  // Cancels the pending source or sink link at idx and the split transfrag
  // that came with it, if any: the link of a node that was split off, from the
  // previous part of the node. The split is the entry next to the link, before
  // or after it depending on where the link was created.
  void cancel_pending_link_and_split(PendingTransfragList& pending,
                                     int32_t idx) {
    const PendingTransfrag link = pending[idx];
    bool from_source = link.from == 0;
    auto is_split = [&](const PendingTransfrag& other) {
      return from_source ? (other.from > 0 && other.to == link.to)
                         : (other.to > 0 && other.from == link.from);
    };

    pending[idx].from = -1;
    if (idx + 1 < count_of(pending) && is_split(pending[idx + 1])) {
      pending[idx + 1].from = -1;
    } else if (idx > 0 && is_split(pending[idx - 1])) {
      pending[idx - 1].from = -1;
    }
  }

  // Deletes one graph edge while pruning to allowed_nodes: from the source if
  // parent_node == 0, to the sink if child_node < 0. Cancels its pending
  // transfrag, updates edge_count and the node_state flags and node estimate.
  void remove_graph_edge(std::vector<CGNode>& node_state, int32_t parent_node,
                         int32_t child_node, int32_t pending_idx,
                         GraphNodeList& nodes,
                         int32_t& edge_count, int32_t& estimated_node_count,
                         CGraphnode *sink,
                         PendingTransfragList& pending_transfrags) {
    if (!parent_node) { // this is a connection to source

      // Cancel the pending source transfrag and unlink the node.
      cancel_pending_link_and_split(pending_transfrags, pending_idx);
      edge_count--;
      // The parent_node is no longer linked to source
      node_state[child_node].future = false;

      // Eliminate child_node from source connection and vice-versa
      erase_value(nodes[child_node]->parent, 0);
      erase_value(nodes[0]->child, child_node);

      // The node can be merged into the previous one if the source link was its
      // only other parent.
      if (can_merge_into_previous(nodes, child_node)) estimated_node_count--;
    } else if (child_node < 0) { // this is a connection to sink
      // Cancel the pending sink transfrag and unlink the node.
      cancel_pending_link_and_split(pending_transfrags, pending_idx);
      edge_count--;
      // The parent_node is no longer linked to sink
      node_state[parent_node].last = false;

      // Eliminate parent_node from sink parents
      erase_value(sink->parent, parent_node);

      // The node can be merged into the next one if the sink link was its only
      // other child.
      if (can_merge_into_next(nodes, parent_node)) estimated_node_count--;
    } else {
      // Remove the child from the parent, then update the state of the parent.
      if (erase_value(nodes[parent_node]->child, child_node)) {
        if (!count_of(nodes[parent_node]->child)
            && !count_of(nodes[parent_node]->parent)) {
          // The parent_node has no children or parents
          node_state[parent_node].keep = false;
          estimated_node_count--;
        } else if (!count_of(nodes[parent_node]->child)
                   && !node_state[parent_node].last) {
          // The parent_node is childless now: link it to the sink
          node_state[parent_node].last = true;
          insert_sorted_unique(sink->parent, parent_node);
          // Add parent_node to pending_transfrags
          pending_transfrags.push_back({float(parent_node), -1.0f, 1.0f});
          edge_count++;
        } else if (can_merge_into_next(nodes, parent_node)) {
          // The node_state parent_node can be merged with next one because
          // there is no other junction from parent_node
          estimated_node_count--;
        }
      }

      // Remove the parent from the child, then update the state of the child.
      if (erase_value(nodes[child_node]->parent, parent_node)) {
        if (!count_of(nodes[child_node]->child)
            && !count_of(nodes[child_node]->parent)) {
          // The child_node has no children or parents
          node_state[child_node].keep = false;
          estimated_node_count--;
        } else if (!count_of(nodes[child_node]->parent)
                   && !node_state[child_node].future) {
          // The child_node is parentless now: link it to the source
          node_state[child_node].future = true;
          insert_sorted_unique(nodes[0]->child, child_node);
          edge_count++;
          // Add source to child_node
          nodes[child_node]->parent.push_back(0);
          pending_transfrags.push_back({0.0f, float(child_node), 1.0f});
          edge_count++;
        } else if (can_merge_into_previous(nodes, child_node)) {
          // The node_state child_node can be merged with previous one because
          // there is no other junction to child_node
          estimated_node_count--;
        }
      }

      edge_count--;
    }
  }

  // Counts the nodes that survive pruning (node_state[].keep). A kept node that
  // collapses into its adjacent predecessor (single parent and child, no gap)
  // is marked as merge and not counted.
  int32_t count_nodes_after_pruning(std::vector<CGNode>& node_state,
                                    int32_t node_count,
                                GraphNodeList& nodes) {
    // First include source and first node_state because this is always kept
    int32_t kept_node_count = 1;

    for (int32_t node_idx = 1; node_idx < node_count; node_idx++) {
      if (node_state[node_idx].keep) {
        bool update_count = true; // by default we update the count
        // Only in the case that previous node_state is adjacent we can collapse
        // nodes; never merge a future node_state
        if (!node_state[node_idx].future && node_state[node_idx - 1].keep
            && !node_state[node_idx - 1].last
            && nodes[node_idx - 1]->end + 1 == nodes[node_idx]->start) {
          // The adjacent parent is the only parent, and this node_state is the
          // only child of it
          if (count_of(nodes[node_idx - 1]->child) == 1
              && nodes[node_idx - 1]->child[0] == node_idx
              && count_of(nodes[node_idx]->parent) == 1) {
            update_count = false;
          }
        }

        if (update_count) {
          kept_node_count++;
          node_state[node_idx].merge = false;
        } else node_state[node_idx].merge = true;
      }
    }
    return (kept_node_count);
  }

  // Gives the neighbors that are kept their new id, and removes the others
  // together with their edge.
  void renumber_neighbors(std::vector<int32_t>& neighbors,
                          std::vector<CGNode>& node_state,
                          int32_t& edge_count) {
    for (int32_t k = count_of(neighbors) - 1; k >= 0; k--) {
      // If we keep it we have to give it the new id
      if (node_state[neighbors[k]].keep) {
        neighbors[k] = node_state[neighbors[k]].id;
      } else { // we need to delete the neighbor
        neighbors.erase(neighbors.begin() + k);
        edge_count--; // delete edge btwn node_state and neighbor
      }
    }
  }

  // Lists the junctions that connect two nodes of the graph and that we can
  // potentially delete.
  std::vector<CGJunc> find_prunable_junctions(const GraphBuildContext& context,
                                              GraphNodeList& nodes,
                                              int32_t node_count,
                                              JunctionList& junction) {
    std::vector<CGJunc> prunable_junctions;

    // Add the junctions that connect two nodes of this graph as candidates.

    int32_t scan_node = 1; // first node_state in graph
    for (int32_t junction_idx = 0; junction_idx < count_of(junction);
         junction_idx++) {
      // All future junctions are out of the graph
      if (junction[junction_idx]->start > nodes[node_count - 1]->start) break;
      // Junction in graph
      if (!junction[junction_idx]->guide_match
          && junction[junction_idx]->strand == context.strand
          && junction[junction_idx]->start >= nodes[1]->end
          && junction[junction_idx]->end <= nodes[node_count - 1]->start) {
        // Advance node_state scan_node until junction start is >= than its
        // start
        while (scan_node < node_count
               && nodes[scan_node]->end
               < junction[junction_idx]->start) scan_node++;
        if (scan_node == node_count) break; // this should never happen
        // End of node_state equal to junction start
        if (nodes[scan_node]->end == junction[junction_idx]->start) {
          int32_t child_node = 0;
          for (int32_t child_pos = 0;
               child_pos < count_of(nodes[scan_node]->child); child_pos++) {
            if (nodes[nodes[scan_node]->child[child_pos]]->start
                == junction[junction_idx]->end) {
              // We found a node_state that is the child of scan_node that
              // junction junction_idx connects
              child_node = nodes[scan_node]->child[child_pos];
              break;
            }
          }
          if (child_node) {
            // Add junction between scan_node and child_node
            CGJunc jc(scan_node, child_node,
                      junction[junction_idx]->effective_n_reads(),
                      junction[junction_idx]->effective_n_reads_good());
            prunable_junctions.push_back(jc);
          }
        }
      }
    }

    return prunable_junctions;
  }

  // Marks the nodes that will be linked to the source or sink (in node_state)
  // and adds those links as candidates for removal.
  void add_source_sink_links_as_prunable(
      GraphNodeList& nodes, PendingTransfragList& pending_transfrags,
      std::vector<CGNode>& node_state,
      std::vector<CGJunc>& prunable_junctions) {
    for (int32_t entry_idx = 0; entry_idx < count_of(pending_transfrags);
         entry_idx++) {
      int32_t to_node = int32_t(pending_transfrags[entry_idx].to);

      // The node_state linked to source
      if (!pending_transfrags[entry_idx].from && to_node > 0) {
        node_state[to_node].future = true;
        // Add source link
        // only for computed links
        if (nodes[to_node - 1]->end + 1 == nodes[to_node]->start) {
          double pos = -entry_idx;
          // Add junction between n and to_node
          CGJunc jc(0, to_node, pos,
            double(pending_transfrags[entry_idx].abundance));
          prunable_junctions.push_back(jc);
        }
      } else if (to_node < 0) { // node_state linked to sink
        int32_t from_node = int32_t(pending_transfrags[entry_idx].from);
        // When to_node==-1 there is a link to sink from from_node*/
        node_state[from_node].last = true;
        // Only for computed links
        if (from_node < count_of(nodes) - 1
            && nodes[from_node]->end + 1 == nodes[from_node + 1]->start) {
          double pos = -entry_idx;
          // Add junction between n and to_node
          CGJunc jc(from_node, -1, pos,
                    double(pending_transfrags[entry_idx].abundance));
          prunable_junctions.push_back(jc);
        }
      }
    }

  }

  // Gives the nodes that are kept their new ids, from 1 to kept_node_count,
  // merges the nodes marked to be merged, and links a node without parents to
  // the source.
  void renumber_kept_nodes(GraphNodeList& nodes,
                           std::vector<CGNode>& node_state, int32_t node_count,
                           int32_t& edge_count) {
    int32_t next_new_id = 1;
    for (int32_t node_idx = 1; node_idx < node_count; node_idx++) {
      if (node_state[node_idx].keep) {
        CGraphnode *curr_node = nodes[node_idx].get();
        if (node_state[node_idx].merge) { // this node_state needs to be merged
          int32_t merge_root_idx = node_idx - 1;
          while (node_state[merge_root_idx].merge) merge_root_idx--;
          CGraphnode *merge_target = nodes[merge_root_idx].get();

          merge_target->end = curr_node->end;
          merge_target->child = curr_node->child;
          node_state[node_idx].keep = false;
          curr_node = merge_target;

        } else {
          node_state[node_idx].id = next_new_id++;

          // Also link a node without parents to the source: the only link that
          // build_reachable_bits_depth_first does not recreate.
          if (!count_of(curr_node->parent)) {
            CGraphnode *source = nodes[0].get();
            source->child.push_back(curr_node->node_id);
            curr_node->parent.push_back(0);
            edge_count++;
          }
          if (!count_of(curr_node->child)) {
            edge_count++;
          }
        }
      }
    }

  }

  // Rewrites the child and parent lists of the nodes, the source and the sink
  // with the new node ids (and edge_count), and removes the nodes that are not
  // kept.
  void rewrite_node_references(const GraphBuildContext& context,
                               GraphNodeList& nodes,
                               std::vector<CGNode>& node_state,
                               int32_t node_count, int32_t& edge_count) {
    // We go from the last node to the first so that the deletions stay valid.
    for (int32_t i = node_count - 1; i > 0; i--) {
      CGraphnode *graph_node = nodes[i].get();

      if (node_state[i].keep) {
        renumber_neighbors(graph_node->child, node_state, edge_count);
        for (int32_t p = count_of(graph_node->parent) - 1; p >= 0; p--) {

          if (graph_node->parent[p]) { // if the parent is not the source
            // If we keep it we have to give it the new id
            if (node_state[graph_node->parent[p]].keep) {
              graph_node->parent[p] = node_state[graph_node->parent[p]].id;
            } else if (node_state[graph_node->parent[p]].merge) {
              int32_t merge_root_idx = graph_node->parent[p] - 1;
              while (node_state[merge_root_idx].merge){
                merge_root_idx--;
              } 
              graph_node->parent[p] = node_state[merge_root_idx].id;
            } else { // we need to delete the parent
              graph_node->parent.erase(graph_node->parent.begin() + p);
              edge_count--;
            }
          }
        }
        graph_node->node_id = node_state[i].id;
      } else { // we need to delete the node_state
        // Need to update edge count first
        for (int32_t c = 0; c < count_of(graph_node->child); c++) {
          // Link to a deleted node_state needs to be removed
          if (!node_state[graph_node->child[c]].keep) edge_count--;
        }
        // This deletes the actual element too if gpvec was not created with
        // free element =false
        nodes.erase(nodes.begin() + i);
      }
    }

    // Deal with the source
    CGraphnode *source = nodes[0].get();
    // (A child to delete should never happen, as we should never delete a link
    // to source.)
    renumber_neighbors(source->child, node_state, edge_count);

    // Deal with the sink.
    renumber_neighbors(context.sink->parent, node_state, edge_count);

  }

  // Updates the pending transfrags to the new numbering and cancels the ones
  // of removed nodes.
  void renumber_pending_transfrags(PendingTransfragList& pending_transfrags,
                                   std::vector<CGNode>& node_state,
                                   int32_t& edge_count) {
    for (PendingTransfrag& pending : pending_transfrags)
      // This also asumes future transcripts to be between source/sink and a
      // node_state -> NOT TRUE
      if (pending.from >= 0) {
        int32_t from_node = int32_t(pending.from);
        int32_t to_node = int32_t(pending.to);

        if (from_node > 0) {
          if (node_state[from_node].keep || node_state[from_node].merge) {
            while (node_state[from_node].merge) from_node--;
            if (to_node > 0) {
              if (node_state[to_node].keep) {
                pending.from = (float)node_state[from_node].id;
                pending.to = (float)node_state[to_node].id;
              } else {
                pending.from = -1;
                edge_count--;
              }
            } else pending.from = (float)node_state[from_node].id;
          } else {
            pending.from = -1;
            edge_count--;
          } // no longer interested in this pending transfrag
        } else { // from_node is the source
          if (node_state[to_node].keep) {
            pending.to = (float)node_state[to_node].id;
          } else {
            pending.from = -1;
            edge_count--;
          }
        }
      }

  }

  // Updates the map from the pre-graph segments to the graph nodes.
  void renumber_segment_nodes(const GraphBuildContext& context,
                              int32_t pre_graph_segment_count,
                              std::vector<CGNode>& node_state) {
    for (int32_t i = 0; i < pre_graph_segment_count; i++) {
      std::vector<CGraphinfo>& graph_infos
          = (*context.segment_nodes)[i];

      for (int32_t j = count_of(graph_infos) - 1; j >= 0; j--) {
        if (graph_infos[j].n_graph == context.graph_idx) {
          if (node_state[graph_infos[j].node_no].keep) {
            graph_infos[j].node_no = node_state[graph_infos[j].node_no].id;
          } else {
            graph_infos.erase(graph_infos.begin() + j);
          }
        }
      }
    }
  }

  // If the graph has more than allowed_nodes nodes, removes its weakest edges
  // first until it fits (remove_graph_edge), then renumbers and merges the kept
  // nodes and rewrites the lists and maps. Returns the new node count.
  int32_t prune_graph_to_node_limit(const GraphBuildContext& context,
                                    int32_t node_count,
                                    int32_t pre_graph_segment_count,
                                    JunctionList& junction, int32_t& edge_count,
                                    PendingTransfragList& pending_transfrags) {
    // This graph's nodes, indexed by node id
    GraphNodeList& nodes
      = context.graph->nodes;

    // Scan all junctions and keep the ones that are in the graph and we can
    // potentially delete
    std::vector<CGJunc> prunable_junctions
      = find_prunable_junctions(context, nodes, node_count, junction);

    // Keep the state of every node: kept, merged, or linked to source or sink.
    std::vector<CGNode> node_state;
    for (int32_t node_id = 0; node_id < node_count; node_id++) { // create nodes
      CGNode state(node_id);
      node_state.push_back(state); // node_state[0] has id 1
    }

    // Mark the nodes that will be linked to the source or sink and add those
    // links as candidates for removal.
    add_source_sink_links_as_prunable(nodes, pending_transfrags, 
      node_state, prunable_junctions);

    // Remove the weakest candidates until the estimated node count fits.
    sort_by(prunable_junctions, compare_graph_junction_by_coverage);

    // Estimate the new node count while we remove the edges.
    int32_t j = 0;
    int32_t estimated_node_count = node_count;
    while (estimated_node_count > allowed_nodes
      && j < count_of(prunable_junctions)) {

      remove_graph_edge(node_state, prunable_junctions[j].left_node,
        prunable_junctions[j].right_node,
        int32_t(-prunable_junctions[j].cov), nodes,
        edge_count, estimated_node_count, context.sink, pending_transfrags);
      j++;
    }

    // Count the nodes that survive pruning.
    int32_t kept_node_count = count_nodes_after_pruning(node_state, 
      node_count, nodes);

    // Re-number nodes that are kept from 1 to kept_node_count instead
    renumber_kept_nodes(nodes, node_state, node_count, edge_count);

    // Rewrite the graph structures with the new node numbers: the child and
    // parent lists (and edge_count), graph_nodes and
    // pre_graph_segment_to_graph_nodes.
    rewrite_node_references(context, nodes, node_state, node_count,
      edge_count);
    renumber_pending_transfrags(pending_transfrags, node_state, edge_count);
    renumber_segment_nodes(context, pre_graph_segment_count, node_state);
    return (kept_node_count);
  }

} // namespace ribbon
