#include "ribbon.h"
#include "GThreads.h"

// Turns the reads of a bundle into "transfrags" (read-supported paths through
// one splice graph) and prepares them for assembly. A transfrag has its ordered
// node ids (nodes), a bit vector (path_bits) over the nodes and edges of the
// graph (see PathBitLayout) and an abundance, pooled and per cell. A group is a
// run of overlapping same-strand reads (built in grouping.cpp):
// group_to_pre_graph_segment maps it to a pre-graph segment and
// pre_graph_segment_to_graph_nodes maps a pre-graph segment to the graph nodes
// that cover it. The transfrag index (TransfragIndex) finds the transfrag of a
// graph with a given path of nodes and edges. assign_fragment_to_transfrags
// adds each read or fragment, then process_transfrags cleans up the set.

namespace ribbon {

  // Key of a path in a TransfragIndex: its nodes in order, each node id doubled
  // and plus 1 if a graph edge (not a gap) leads to it from the previous node.
  // Does not handle paths that include the source node.
  std::vector<int32_t> make_path_key(PathBitLayout& layout,
                                     std::vector<int32_t>& nodes,
                                     BitVec& path_bits) {
    std::vector<int32_t> key;
    key.reserve(nodes.size());
    for (int32_t node_idx = 0; node_idx < count_of(nodes); node_idx++) {
      bool by_edge = false;
      if (node_idx) { // not the first node in path_bits
        // There is an edge between node[n-1] and node[n]
        by_edge = layout.has_edge_bit(path_bits, nodes[node_idx - 1],
                                      nodes[node_idx]);
      }
      key.push_back(2 * nodes[node_idx] + by_edge);
    }
    return key;
  }

  // True if the transfrag is not a guide, has less abundance than threshold,
  // and neither starts at the source nor ends at the sink.
  static bool is_weak_internal_transfrag(const CTransfrag& transfrag,
                                         int32_t node_count, float threshold) {
    return !transfrag.guide && transfrag.abundance < threshold
           && transfrag.nodes[0]
           && transfrag.nodes.back() < node_count - 1;
  }

  // Discards internal (not linked to the source or sink) non-guide transfrags
  // with an abundance below threshold, and from the index too. While more than
  // max_transfrag_number remain, raises the threshold by 1 and repeats.
  void remove_low_abundance_transfrags(SpliceGraph& graph, float threshold) {
    // Remove the weak internal transfrags.
    for (int32_t t = count_of(graph.transfrag) - 1; t >= 0; t--) {
      // Need to delete transfrag that doesn't come from source or ends at sink
      if (is_weak_internal_transfrag(*graph.transfrag[t], graph.node_count,
          threshold)) {
        // This should be eliminated if we want to store transcripts from 0 node
        graph.transfrag_index.erase(make_path_key(graph.path_bit_layout,
            graph.transfrag[t]->nodes, graph.transfrag[t]->path_bits));
        std::swap(graph.transfrag[t], graph.transfrag.back());
        graph.transfrag.pop_back();
      }
    }

    // Raise the threshold until at most max_transfrag_number transfrags remain.
    while (count_of(graph.transfrag) > max_transfrag_number) {
      threshold++;
      for (int32_t t = count_of(graph.transfrag) - 1; t >= 0; t--)
        // Need to delete transfrag that doesn't come from source
        if (is_weak_internal_transfrag(*graph.transfrag[t], graph.node_count,
            threshold)) {
          // This should be eliminated if we want to store transcripts from 0
          // node
          graph.transfrag_index.erase(make_path_key(graph.path_bit_layout,
              graph.transfrag[t]->nodes, graph.transfrag[t]->path_bits));
          std::swap(graph.transfrag[t], graph.transfrag.back());
          graph.transfrag.pop_back();
        }
    }
  }

  // Comparator for CTransfrag: sorts transfrags by decreasing reach (last
  // node - first node), then by decreasing node count, decreasing path_bits
  // size (nodes + edges set) and decreasing abundance.
  int32_t compare_transfrag_by_reach(const void *p1, const void *p2) {
    CTransfrag *a = (CTransfrag *)p1;
    CTransfrag *b = (CTransfrag *)p2;

    // By decreasing genomic span (last node - first node)
    // transfrag with further reach comes first
    if (a->nodes.back() - a->nodes[0] < b->nodes.back() - b->nodes[0]) return 1;
    if (a->nodes.back() - a->nodes[0] > b->nodes.back() - b->nodes[0])
      return -1;

    // By decreasing node count
    // transfrag with more nodes comes first
    if (count_of(a->nodes) < count_of(b->nodes)) return 1;
    if (count_of(a->nodes) > count_of(b->nodes)) return -1;

    // By decreasing path_bits completeness
    if (a->path_bits.count() < b->path_bits.count()) return 1;
    if (a->path_bits.count() > b->path_bits.count()) return -1;

    // By decreasing abundance
    if (a->abundance < b->abundance) return 1;
    if (a->abundance > b->abundance) return -1;
    return 0;
  }

  // Inserts transfrag index t into transfrag (a node's sorted list of
  // transfrags) by binary search over [first,last] unless t is already there.
  // Returns true if t was inserted.
  bool insert_transfrag_into_node_list(int32_t t,
                                       std::vector<int32_t>& transfrag,
                                       int32_t first,
                                       int32_t last) {
    while (first <= last) {
      int32_t mid = (first + last) / 2; // compute mid point.
      if (t > transfrag[mid]) first = mid + 1;
      else if (t < transfrag[mid]) last = mid - 1;
      else return false; // transcript already inserted
    }

    int32_t insert_pos = first;

    transfrag.insert(transfrag.begin() + insert_pos, t);
    return true;
  }

  // Registers the gapped transfrag t, whose path from from_node to to_node is
  // unknown, with every node that is a child of from_node and a parent of
  // to_node. Returns true if at least one node got the transfrag.
  bool assign_incomplete_transfrag_to_nodes(int32_t t, int32_t from_node,
                                            int32_t to_node,
                                            GraphNodeList& graph_nodes) {
    // A better version of this function could take into account also the
    // fragments inner length

    bool assigned = false;

    // Try every node between them that could lie on the missing path.
    for (int32_t i = from_node + 1; i < to_node; i++) {
      // Node i is a child of from_node, and to_node is a child of node i
      if (graph_nodes[from_node]->child_reachable_bits[i]
          && graph_nodes[i]->child_reachable_bits[to_node]) {
        assigned = insert_transfrag_into_node_list(t,
            graph_nodes[i]->transfrag_ids, 0,
            count_of(graph_nodes[i]->transfrag_ids) - 1) or assigned;
      }
    }

    return assigned;
  }

  // Returns whether transfrag t is real (unambiguous). At each gap in its
  // path_bits we follow the graph to the next node; if several children lead
  // there, the transfrag branches and is not real.
  bool resolve_transfrag_paths(int32_t t, SpliceGraph& graph) {
    // Index of the current node (the one we are looking at) in transfrag
    int32_t i = 0;
    // Index of last node in transfrag
    int32_t last_idx = count_of(graph.transfrag[t]->nodes) - 1;
    // start at first node in transfrag
    int32_t node_id = graph.transfrag[t]->nodes[i];
    double total_coverage = 0;
    while (i < last_idx) {
      // If next node in transfrag is adjacent -> there is only one way to
      // continue on transfrag path
      if (node_id + 1 == graph.transfrag[t]->nodes[i + 1]) {
        node_id = graph.transfrag[t]->nodes[i + 1];
        i++;
      } else { // there might be multiple ways to reach next node from node
        // If there is an edge between node and next_node in transfrag -> we
        // don't have discontinuity there
        if (graph.path_bit_layout.has_edge_bit(graph.transfrag[t]->path_bits,
            node_id,
                                graph.transfrag[t]->nodes[i + 1])) {
          node_id = graph.transfrag[t]->nodes[i + 1];
          i++;
        // Found a potential discontinuity in transfrag path_bits --> try to see
        // if there are different branches going on from here
        } else {
          CGraphnode *node = graph.nodes[node_id].get();
          int32_t next_node_id = graph.transfrag[t]->nodes[i + 1];
          // Total number of transfrags in node
          int32_t node_transfrag_count = count_of(node->transfrag_ids);
          // Children of node through which next_node can be reached
          std::vector<int32_t> branch_children;
          // Total abundance of transfrags going from node to children of node
          total_coverage = 0;

          // We check if there are different paths to reach next_node from node
          for (int32_t j = 0; j < count_of(node->child); j++) {
            // If child of node is after next_node -> we don't care
            if (node->child[j] > next_node_id) break;
            // We found a node on the path to next_node
            if (node->child[j] == next_node_id
                || graph.nodes[node->child[j]]
                       ->child_reachable_bits[next_node_id]) {
              float branch_abundance = 0;
              int32_t *edge_pos = graph.path_bit_layout.find_edge_bit(node_id,
                                                       node->child[j]);
              if (edge_pos)
                for (int32_t f = 0; f < node_transfrag_count; f++) {
                  // This is a transfrag that goes through node and child;
                  // ATTENTION: it might not reach next_node!! or be compatible
                  // with next_node
                  if (graph.transfrag[node->transfrag_ids[f]]
                          ->path_bits[*edge_pos]) {
                    branch_abundance
                        += graph.transfrag[node->transfrag_ids[f]]->abundance;
                    total_coverage
                        += graph.transfrag[node->transfrag_ids[f]]->abundance;
                  }
                }
              if (branch_abundance) branch_children.push_back(node->child[j]);
            }
          }
          // There is no way to continue
          if (!total_coverage || branch_children.empty()) return true;
          if (branch_children.size() == 1) { // only one way to go forward
            node_id = branch_children[0];
            if (node_id == next_node_id) i++;
          // Break at first discontinue point in transfrag -> might be many
          // other ways to continue but this functions doesn't care
          } else break;
        }
      }
    }
    if (i == last_idx) return true;

    // This implementation only assumes one break in transfrag
    return false;
  }

  // Extends the path of a guide's transfrag with the source before its first
  // node and the sink after its last node, setting the bits of both nodes and
  // both edges.
  void extend_path_to_source_and_sink(CTransfrag& transfrag,
                                      PathBitLayout& layout,
                                      int32_t node_count) {
    transfrag.nodes.insert(transfrag.nodes.begin(), 0);
    transfrag.path_bits[0] = 1;
    layout.set_edge_bit(transfrag.path_bits, 0, transfrag.nodes[1]);

    int32_t sink = node_count - 1;
    transfrag.nodes.push_back(sink);
    transfrag.path_bits[sink] = 1;
    layout.set_edge_bit(transfrag.path_bits,
                        transfrag.nodes[count_of(transfrag.nodes) - 2],
                        transfrag.nodes.back());
  }

  // Prepares the transfrags of one graph for find_transcripts: adds the guide
  // paths and the transfrags of uncovered edges, drops weak ones, sorts them by
  // reach, registers them with their nodes and resolves the incomplete ones.
  void process_transfrags(SpliceGraph& graph,
                          std::vector<CGuide>& matched_guides,
                          BundleData *bundle) {
    int32_t cell_count = count_of(bundle->cell_name);

    // Add each guide's path as a backbone transfrag (an incompatible transfrag
    // may connect its nodes) and extend the guide's path to source and sink.
    for (int32_t i = 0; i < count_of(matched_guides); i++) {
      CTransfrag& guide_transfrag = *matched_guides[i].transfrag;
      std::unique_ptr<CTransfrag> backbone;
      if (guide_transfrag.guide) {
        backbone = std::make_unique<CTransfrag>(
            cell_count, guide_transfrag.nodes, guide_transfrag.path_bits,
            min_transfrag_abundance * ERROR_PERC);
      }

      extend_path_to_source_and_sink(guide_transfrag, graph.path_bit_layout,
          graph.node_count);

      if (backbone) {
        backbone->guide = 1 + matched_guides[i].g;

        backbone->long_start = graph.nodes[backbone->nodes[0]]->start;
        backbone->long_end = graph.nodes[backbone->nodes.back()]->end;
        graph.transfrag.push_back(std::move(backbone));
      }
    }

    // Drop the weak transfrags that do not start at the source (noise).
    remove_low_abundance_transfrags(graph, min_transfrag_abundance);

    // Collect the nodes and edges that the transfrags cover.
    BitVec all_path_bits(graph.node_count + graph.edge_count);

    for (int32_t t = 0; t < count_of(graph.transfrag); t++)
      all_path_bits = all_path_bits | graph.transfrag[t]->path_bits;

    // Add a minimal transfrag for every edge that no transfrag covers.
    for (int32_t i = 1; i < graph.node_count - 1; i++) {
      CGraphnode *n = graph.nodes[i].get();

      for (int32_t c = 0; c < count_of(n->child); c++) {
        int32_t *edge_pos = graph.path_bit_layout.find_edge_bit(i, n->child[c]);
        if (edge_pos && !all_path_bits[*edge_pos]) {
          std::vector<int32_t> nodes;
          nodes.push_back(i);
          nodes.push_back(n->child[c]);
          BitVec edge_path_bits(graph.node_count + graph.edge_count);
          edge_path_bits[i] = 1;
          edge_path_bits[n->child[c]] = 1;
          edge_path_bits[*edge_pos] = 1;
          graph.transfrag.push_back(std::make_unique<CTransfrag>(
              cell_count, nodes, edge_path_bits, min_transfrag_abundance));
        }
      }
    }

    // Sort the transfrags: longest reach first, then more nodes, then larger
    // abundance; the transfrags with one node come last.
    sort_by(graph.transfrag, compare_transfrag_by_reach);

    // Remembers incomplete transfrags (the ones that don't have edges between
    // two consecutive nodes
    std::vector<int32_t> incomplete_transfrags;

    // Register each transfrag with the nodes it passes through. We go in the
    // sorted order, which matters later; a transfrag with gaps is also
    // registered with the nodes on its possible missing paths.
    for (int32_t transfrag_idx = 0; transfrag_idx < count_of(graph.transfrag);
         transfrag_idx++) {
      // Update nodes
      int32_t transfrag_node_count
          = count_of(graph.transfrag[transfrag_idx]->nodes);

      // Add transfrag to nodes' in and out; if a transfrag only has one node
      // then it is not added to a node; we might want to change this for the
      // computation of fpkm's
      if (transfrag_node_count > 1) {
        bool incomplete = false;

        // For all nodes in transfrag
        for (int32_t n = 0; n < transfrag_node_count; n++) {
          // Not first or last node
          if (n && n < count_of(graph.transfrag[transfrag_idx]->nodes) - 1) {
            // Add t1 to in and out of node
            graph.nodes[graph.transfrag[transfrag_idx]->nodes[n]]
                ->transfrag_ids.push_back(transfrag_idx);

            if (graph.transfrag[transfrag_idx]->nodes[n - 1]
                && graph.transfrag[transfrag_idx]->nodes[n] < graph.node_count
                    - 1) {
              // Check if transfrag t1 is incomplete between node[n-1] and node
              // [n]
              // There is no edge between node[n-1] and node[n]
              if (!graph.path_bit_layout.has_edge_bit(
                      graph.transfrag[transfrag_idx]->path_bits,
                      graph.transfrag[transfrag_idx]->nodes[n - 1],
                      graph.transfrag[transfrag_idx]->nodes[n]))
                // This ensures that we still have compatibilities when going
                // through a certain path: !!! THIS IS NEW COMPARED TO PERL
                incomplete = assign_incomplete_transfrag_to_nodes(transfrag_idx,
                    graph.transfrag[transfrag_idx]->nodes[n - 1],
                    graph.transfrag[transfrag_idx]->nodes[n],
                    graph.nodes) or incomplete;
            }
          } else if (n) { // last but not first node
            // Add t1 to in of node
            graph.nodes[graph.transfrag[transfrag_idx]->nodes[n]]
                ->transfrag_ids.push_back(transfrag_idx);

            if (graph.transfrag[transfrag_idx]->nodes[n - 1]
                && graph.transfrag[transfrag_idx]->nodes[n] < graph.node_count
                    - 1) {
              // Check if transfrag t1 is incomplete between node[n-1] and node
              // [n]
              // There is no edge between node[n-1] and node[n]
              if (!graph.path_bit_layout.has_edge_bit(
                      graph.transfrag[transfrag_idx]->path_bits,
                      graph.transfrag[transfrag_idx]->nodes[n - 1],
                      graph.transfrag[transfrag_idx]->nodes[n]))
                // This ensures that we still have compatibilities when going
                // through a certain path: !!! THIS IS NEW COMPARED TO PERL
                incomplete = assign_incomplete_transfrag_to_nodes(transfrag_idx,
                    graph.transfrag[transfrag_idx]->nodes[n - 1],
                    graph.transfrag[transfrag_idx]->nodes[n],
                    graph.nodes) or incomplete;
            }
          } else { // first node -> only add transfrag to out of node
            graph.nodes[graph.transfrag[transfrag_idx]->nodes[n]]
                ->transfrag_ids.push_back(transfrag_idx);
          }
        }

        if (incomplete) incomplete_transfrags.push_back(transfrag_idx);
        else graph.transfrag[transfrag_idx]->real = true;
      }
    }

    // Resolve the incomplete transfrags: those with a unique path are real.
    for (int32_t t = 0; t < count_of(incomplete_transfrags); t++)
      graph.transfrag[incomplete_transfrags[t]]->real
          = resolve_transfrag_paths(incomplete_transfrags[t], graph);

    // Set the abundance of the transfrag from the source to a node whose only
    // parent is the source to the abundance of the node's other transfrags.
    CGraphnode *source = graph.nodes[0].get();
    for (int32_t i = 0; i < count_of(source->child); i++) {
      double abundance = 0;
      int32_t link_transfrag = -1;
      // Source is the only parent of node
      if (count_of(graph.nodes[source->child[i]]->parent) == 1
          && !graph.nodes[source->child[i]]->parent[0]) {
        for (int32_t j = 0;
             j < count_of(graph.nodes[source->child[i]]->transfrag_ids); j++) {
          int32_t t = graph.nodes[source->child[i]]->transfrag_ids[j];
          if (graph.transfrag[t]->nodes.back() == source->child[i])
            link_transfrag = t;
          else abundance += graph.transfrag[t]->abundance;
        }
        // Found transfrag from source to node and the transfrag wasn't deleted
        if (link_transfrag > -1 && graph.transfrag[link_transfrag]->abundance) {
          graph.transfrag[link_transfrag]->abundance = abundance;
        }
      }
    }

    // Set the abundance of a transfrag that links a node to the sink (if it is
    // at most 1) from the node coverage minus its other transfrags.
    CGraphnode *sink = graph.nodes[graph.node_count - 1].get();
    for (int32_t i = 0; i < count_of(sink->parent); i++) {
      double abundance = 0;
      int32_t link_transfrag = -1;
      int32_t transfrag_count
          = count_of(graph.nodes[sink->parent[i]]->transfrag_ids);
      for (int32_t j = 0; j < transfrag_count; j++) {
        int32_t t = graph.nodes[sink->parent[i]]->transfrag_ids[j];
        if (graph.transfrag[t]->nodes.back() == graph.node_count
            - 1) link_transfrag = t;
        else abundance += graph.transfrag[t]->abundance;
      }
      if (link_transfrag > -1 && graph.transfrag[link_transfrag]->abundance
          && graph.transfrag[link_transfrag]->abundance <= 1) {
        graph.transfrag[link_transfrag]->abundance
            = graph.nodes[sink->parent[i]]->cov
            / graph.nodes[sink->parent[i]]->len() - abundance;
        if (graph.transfrag[link_transfrag]->abundance < 1)
          graph.transfrag[link_transfrag]->abundance = 1;
      }
    }
  }

  // Walks the segments of a read over the graph nodes of its groups, adds its
  // coverage to each node it overlaps and lists the touched nodes, per graph:
  // read_graph_idxs holds the graphs and read_nodes[i] the nodes in graph i.
  void get_read_graph_nodes(
      char strand, float read_coverage, std::vector<int32_t>& read_graph_idxs,
      float strand_prop, std::vector<std::vector<int32_t>>& read_nodes,
      ReadList& reads, int32_t read_idx,
      std::vector<std::vector<int32_t>>& read_group,
      std::vector<int32_t>& merge_map,
      PerStrand<std::vector<int32_t>>& group_to_pre_graph_segment,
      PerStrand<SegmentNodesList>& segment_nodes,
      PerStrand<std::vector<SpliceGraph>>& splice_graphs) {
    // Where we are in the read: the last node added and the current segment.
    int32_t last_node_id = -1;
    int32_t last_graph_idx = -1;
    int32_t segment_count = count_of(reads[read_idx]->segs);
    // Need to keep track of coordinates already added to coverages of
    // graphnodes
    int32_t seg_idx = 0;

    std::unordered_set<int32_t> seen_pre_graph_segments;

    // How can a read be associated to multiple groups? ---> If it is spliced
    for (int32_t i = 0; i < count_of(read_group[read_idx]); i++) {
      int32_t group_id = read_group[read_idx][i];
      while (merge_map[group_id] != group_id) group_id = merge_map[group_id];
      read_group[read_idx][i] = group_id;
      // Group was associated to a pre-graph
      int32_t pre_graph_segment_id
          = group_to_pre_graph_segment[strand][group_id];
      // Group has a pre-graph segment associated with it and the pre-graph was
      // processed
      if (pre_graph_segment_id > -1
          && count_of(
                 segment_nodes[strand][pre_graph_segment_id])) {
        // First time this pre-graph segment is seen
        if (seen_pre_graph_segments.insert(pre_graph_segment_id).second) {
          // Number of graph nodes in the pre-graph segment
          int32_t graph_node_count = count_of(
              segment_nodes[strand][pre_graph_segment_id]);
          int32_t j = 0;

          // Go through the graph nodes of this pre-graph segment, adding
          // coverage.
          while (j < graph_node_count && seg_idx < segment_count) {
            int32_t graph_idx
                = segment_nodes[strand][pre_graph_segment_id][j]
                      .n_graph;
            int32_t node_id
                = segment_nodes[strand][pre_graph_segment_id][j]
                      .node_no;
            CGraphnode *node
                = splice_graphs[strand][graph_idx].nodes[node_id].get();
            // Compute graphnode coverage by read here (assuming they come in
            // order on the genomic line)
            bool overlaps_node = false;
            // Found how much of reads[n] intersects node represented by the
            // pre-graph segment;
            // here we assume that read is always included in
            // segment_nodes
            while (seg_idx < segment_count) {
              if (reads[read_idx]->segs[seg_idx].end < node->start) seg_idx++;
              else {
                int32_t overlap_len
                    = reads[read_idx]->segs[seg_idx].overlapLen(node);
                if (overlap_len) {
                  overlaps_node = true;

                  node->cov += strand_prop * overlap_len * read_coverage;
                  if (reads[read_idx]->segs[seg_idx].end <= node->end)
                    seg_idx++;
                  else break;
                } else break;
              }
            }

            // Record the node in the node list of its graph.
            if (overlaps_node) { // read intersects node
              int32_t graph_pos = -1; // position of graph in list
              if (count_of(read_graph_idxs)) {
                if (graph_idx != read_graph_idxs.back()) {
                  // Count should be small
                  for (int32_t r = 0; r < count_of(read_graph_idxs) - 1; r++) {
                    if (graph_idx == read_graph_idxs[r]) {
                      graph_pos = r;
                      break;
                    }
                  }
                } else graph_pos = count_of(read_graph_idxs) - 1;
              }
              // Add n_graph to list of graphs if not seen before
              if (graph_pos < 0) {
                read_graph_idxs.push_back(graph_idx);
                graph_pos = count_of(read_graph_idxs) - 1;
              }

              // This is a new graph, or a new node
              if (last_graph_idx != graph_idx || last_node_id != node_id) {
                read_nodes[graph_pos].push_back(node_id);
              }

              last_node_id = node_id;
              last_graph_idx = graph_idx;
            }

            j++;
          }
        }
      }
    }
  }

  // Adds the abundance of a read (or pair) with this node path to its
  // transfrag, creating it if needed, and to the cell coverage of the nodes
  // it covers (and abundance x raw_scale to their raw coverage, see
  // weight_reads_by_molecule). Updates long_start and long_end; returns NULL
  // for one node.
  CTransfrag *add_read_abundance_to_transfrag(
      SpliceGraph& graph, std::vector<int32_t>& nodes, BitVec& path_bits,
      float abundance, float raw_scale, uint32_t read_start, uint32_t read_end,
      int32_t cell_count, int32_t cell, uint32_t gap_start, uint32_t gap_end) {
    // Add the read's coverage in this cell to each node, without the bases it
    // does not cover (before its start, after its end or in the mate gap).
    for (int32_t i = 0; i < count_of(nodes); i++) {
      CGraphnode *node = graph.nodes[nodes[i]].get();
      int32_t len = node->len();
      // Read starts after node
      if (node->start < read_start) len -= read_start - node->start;
      // Read ends before node
      if (node->end > read_end) len -= node->end - read_end;

      // Gap intersects node -> subtract from covered length
      if (gap_end && node->start <= gap_end && gap_start <= node->end) {
        len -= std::min(node->end, gap_end) - std::max(node->start,
            gap_start) + 1;
      }

      if (len < 0) len = 0;

      node->cell_cov[cell] += abundance * len;
      node->cell_cov_raw[cell] += abundance * raw_scale * len;
    }

    // A path with one node is not worth keeping as a transfrag.
    if (count_of(nodes) == 1) return (NULL);

    // Find the transfrag with this path, or create it and add it to the index.
    std::vector<int32_t> path_key
        = make_path_key(graph.path_bit_layout, nodes, path_bits);
    auto found = graph.transfrag_index.find(path_key);
    CTransfrag *tf = found == graph.transfrag_index.end() ? NULL
        : found->second;
    if (!tf) { // t is NULL
      graph.transfrag.push_back(
          std::make_unique<CTransfrag>(cell_count, nodes, path_bits, 0));
      tf = graph.transfrag.back().get();
      graph.transfrag_index[path_key] = tf;
    }

    // Add the abundance to the transfrag, pooled and for this cell.
    tf->abundance += abundance;
    tf->cell_abundance[cell] += abundance;

    // Extend long_start and long_end when the read ends inside the first or last
    // node.
    if (graph.nodes[nodes[0]]->start <= read_start)
      // Value of long_start inside node; the other extensions could come from
      // spliced nodes
      if (!tf->long_start || read_start < tf->long_start)
        tf->long_start = read_start;
    if (graph.nodes[nodes.back()]->end >= read_end)
      // Value of long_end inside node; the other extensions could come from
      // spliced nodes
      if (!tf->long_end || read_end > tf->long_end) tf->long_end = read_end;

    return (tf);
  }

  // Where a fragment goes: its reads, the strand and graphs of the tables it is
  // placed in, and the span and abundance of the transfrags it makes there.
  struct FragmentPlacement {
    ReadList& reads;
    GraphTables& tables;
    char strand;
    int32_t read_idx;
    int32_t mate_idx;
    float abundance; // the share of the fragment on this strand
    // Raw / deduplicated count of the fragment (1 without --sq3)
    float raw_scale;
    uint32_t read_start;
    uint32_t read_end;
    // Start/end of the inner gap between the two paired reads in the fragment
    uint32_t gap_start;
    uint32_t gap_end;
    int32_t cell_count;
    int32_t cell;
  };

  // Sets the share of a fragment (a read and its mate, if any) for each strand
  // (strand_prop): from the read strand or, if unstranded, from the neg_prop of
  // its groups.
  void compute_strand_prop(ReadList& reads, AssemblyState& state,
                           int32_t read_idx, int32_t mate_idx,
                           PerStrand<float>& strand_prop) {
    // Adds the neg_prop of the groups of a read to strand_prop['-'] and counts
    // the groups.
    auto add_group_neg_prop = [&](int32_t idx, int32_t& group_count) {
      for (int32_t i = 0; i < count_of(state.read_group[idx]); i++) {
        group_count++;
        // Read is unstranded => it should belong to one group only -> not
        // necessarily see above
        int32_t group_id = find_root(state.merge_map, state.read_group[idx][i]);
        strand_prop['-'] += state.group[group_id]->neg_prop;
      }
    };

    // Takes the strand share from one read of the fragment: the average of the
    // neg_prop of its groups if it is unstranded, else all on its strand.
    auto set_strand_prop_from_read = [&](int32_t idx) {
      // The paired read is not present otherwise it would have the same
      // strand from assign_read_to_groups
      if (reads[idx]->strand == '.') {
        // Due to errors read might not belong to any groups
        int32_t group_count = 0;
        add_group_neg_prop(idx, group_count);
        if (group_count) {
          strand_prop['-'] /= group_count;
          strand_prop['+'] = 1 - strand_prop['-'];
        }
      } else {
        if (reads[idx]->strand == '-') strand_prop['-'] = 1;
        else strand_prop['+'] = 1;
      }
    };

    // Both reads are unstranded: average the neg_prop of the groups of both.
    auto is_unstranded_with_hits = [&](int32_t idx) {
      return reads[idx]->nh && reads[idx]->strand == '.';
    };
    if (is_unstranded_with_hits(read_idx) && mate_idx > -1
        && is_unstranded_with_hits(mate_idx)) {
      // due to errors read might not belong to any groups
      int32_t group_count = 0;
      add_group_neg_prop(read_idx, group_count);
      add_group_neg_prop(mate_idx, group_count);

      if (group_count) {
        strand_prop['-'] /= group_count;
        strand_prop['+'] = 1 - strand_prop['-'];
      }
    } else {
      if (reads[read_idx]->nh) {
        set_strand_prop_from_read(read_idx);
      } else if (mate_idx > -1 && reads[mate_idx]->nh) { // reads[n] is deleted
        set_strand_prop_from_read(mate_idx);
      }
    }
  }

  // Adds a path_bits with its nodes as a transfrag of the graph.
  void add_path_bits(const FragmentPlacement& placement, int32_t graph,
                   BitVec& path_bits, std::vector<int32_t>& nodes) {
    GraphTables& tables = placement.tables;
    char strand = placement.strand;
    add_read_abundance_to_transfrag(
        tables.splice_graphs[strand][graph], nodes, path_bits,
        placement.abundance, placement.raw_scale, placement.read_start,
        placement.read_end,
        placement.cell_count, placement.cell, placement.gap_start,
        placement.gap_end);
  }

  // Sets the path_bits of the mate's nodes in the graph: the nodes, and the
  // edges between them that are continuous or spliced by a junction of the
  // mate.
  void build_mate_path_bits(const FragmentPlacement& placement, int32_t graph,
                          std::vector<int32_t>& nodes, BitVec& path_bits) {
    GraphTables& tables = placement.tables;
    char strand = placement.strand;
    CReadAln& mate = *placement.reads[placement.mate_idx];
    GraphNodeList& graph_node_list = tables.splice_graphs[strand][graph].nodes;
    PathBitLayout& layout = tables.splice_graphs[strand][graph].path_bit_layout;
    int32_t node_count = tables.splice_graphs[strand][graph].node_count;
    path_bits.clear();
    path_bits.resize(node_count
        + tables.splice_graphs[strand][graph].edge_count);
    int32_t i = 0;
    for (int32_t j = 0; j < count_of(nodes); j++) {
      path_bits[nodes[j]] = 1;
      if (j) { // add edge bit
        // Continuous node
        if (graph_node_list[nodes[j]]->start - 1
            == graph_node_list[nodes[j - 1]]->end) {
          layout.set_edge_bit(path_bits, nodes[j - 1], nodes[j]);
        } else { // non-continuous node
          while (i < count_of(mate.juncs)
                 && mate.juncs[i]->start
                 < graph_node_list[nodes[j - 1]]->end) i++;
          bool mate_junction_links_nodes = i < count_of(mate.juncs)
              && mate.juncs[i]->strand != '.'
              && mate.juncs[i]->start == graph_node_list[nodes[j - 1]]->end
              && mate.juncs[i]->end == graph_node_list[nodes[j]]->start;
          if (mate_junction_links_nodes) {
            layout.set_edge_bit(path_bits, nodes[j - 1], nodes[j]);
          }
        }
      }
    }
  }

  // Sets the path_bits of the read's nodes in the graph: the nodes, and the
  // edges between them that are continuous or spliced by a junction of the
  // read. If the read goes between two nodes without a valid junction, the
  // nodes so far are added as a transfrag of their own and the path_bits
  // restarts at the second node. The nodes that were added are removed from
  // nodes.
  void build_read_path_bits(const FragmentPlacement& placement, int32_t graph,
                          std::vector<int32_t>& nodes, BitVec& path_bits) {
    GraphTables& tables = placement.tables;
    char strand = placement.strand;
    CReadAln& read = *placement.reads[placement.read_idx];
    GraphNodeList& graph_node_list = tables.splice_graphs[strand][graph].nodes;
    PathBitLayout& layout = tables.splice_graphs[strand][graph].path_bit_layout;
    int32_t node_count = tables.splice_graphs[strand][graph].node_count;
    path_bits.clear();
    path_bits.resize(node_count
        + tables.splice_graphs[strand][graph].edge_count);
    int32_t i = 0;
    int32_t flush_start = 0;
    for (int32_t j = 0; j < count_of(nodes); j++) {
      path_bits[nodes[j]] = 1;
      if (j) { // add edge bit
        // Continuous node
        if (graph_node_list[nodes[j]]->start - 1
            == graph_node_list[nodes[j - 1]]->end) {
          layout.set_edge_bit(path_bits, nodes[j - 1], nodes[j]);
        } else { // non-continuous node
          // End of junction should be after start of node
          while (i < count_of(read.juncs)
                 && read.juncs[i]->end
                 < graph_node_list[nodes[j]]->start) i++;
          if (i < count_of(read.juncs)
              && read.juncs[i]->start
              <= graph_node_list[nodes[j - 1]]->end) {
            // Valid junction
            if (read.juncs[i]->strand != '.'
                && read.juncs[i]->start == graph_node_list[nodes[j - 1]]->end
                && read.juncs[i]->end == graph_node_list[nodes[j]]->start) {
              layout.set_edge_bit(path_bits, nodes[j - 1], nodes[j]);
            } else { // update abundance for what we have so far
              path_bits[nodes[j]] = 0;
              std::vector<int32_t> node(nodes.begin() + flush_start,
                                        nodes.begin() + j);
              add_path_bits(placement, graph, path_bits, node);
              path_bits.reset();
              // Restart the path_bits
              path_bits[nodes[j]] = 1;
              flush_start = j;
            }
          }
        }
      }
    }

    // Drop the nodes that we already added as separate transfrags.
    nodes.erase(nodes.begin(), nodes.begin() + flush_start);
  }

  // Turns a read, or a pair seen as a fragment, into transfrag abundance: finds
  // its strand(s), the nodes and edges it touches in each graph (read and mate
  // path_bits) and adds them as one transfrag if compatible, else separately.
  void assign_fragment_to_transfrags(BundleData *bundle,
                                     AssemblyState& state, GraphTables& tables,
                                     int32_t read_idx, int32_t mate_idx,
                                     float read_coverage, float raw_scale) {
    // The bundle's reads and groups, and the graphs they go into
    ReadList& reads = bundle->reads;
    int32_t cell_count = count_of(bundle->cell_name);

    // Span of the fragment and gap between the mates (if there is a mate).
    uint32_t read_start = reads[read_idx]->start;
    uint32_t read_end = reads[read_idx]->end;
    // Start/end of the inner gap between the two paired reads in the fragment
    uint32_t gap_start = read_end + 1;
    uint32_t gap_end = 0;
    if (mate_idx > -1) {
      if (reads[mate_idx]->start > read_end + 1)
        gap_end = reads[mate_idx]->start - 1;
      if (reads[mate_idx]->end > read_end) read_end = reads[mate_idx]->end;
    }

    // Share of the fragment for each strand
    PerStrand<float> strand_prop;
    compute_strand_prop(reads, state, read_idx, mate_idx, strand_prop);

    // For each strand that gets a share of the fragment, get the nodes that the
    // read and its mate touch in each graph.
    for (char strand : {'-', '+'})
      if (strand_prop[strand]) {
        FragmentPlacement placement = {reads, tables, strand, read_idx,
            mate_idx, strand_prop[strand] * read_coverage, raw_scale,
            read_start, read_end, gap_start, gap_end, cell_count,
            reads[read_idx]->bc};

        std::vector<int32_t> read_graphs;
        std::vector<std::vector<int32_t>> read_nodes(
            count_of(state.read_group[read_idx]));
        if (reads[read_idx]->nh)
          get_read_graph_nodes(strand, read_coverage, read_graphs,
                               strand_prop[strand], read_nodes, reads, read_idx,
                               state.read_group, state.merge_map,
                                   state.group_to_pre_graph_segment,
                               tables.pre_graph_segment_to_graph_nodes,
                                   tables.splice_graphs);

        std::vector<int32_t> mate_graphs;
        std::vector<std::vector<int32_t>> mate_nodes;
        // Get pair path_bits if pair exists and it hasn't been deleted
        if (mate_idx > -1 && reads[mate_idx]->nh) {
          mate_nodes.resize(count_of(state.read_group[mate_idx]));
          get_read_graph_nodes(strand, read_coverage, mate_graphs,
                               strand_prop[strand], mate_nodes, reads, mate_idx,
                               state.read_group, state.merge_map,
                                   state.group_to_pre_graph_segment,
                               tables.pre_graph_segment_to_graph_nodes,
                                   tables.splice_graphs);
        }

        // Build the read's path_bits in each graph and add it as a transfrag,
        // merged with the mate's path_bits if they are compatible.
        BitVec read_path_bits;
        BitVec mate_path_bits;
        int32_t mates_used = 0;

        for (int32_t graph_pos = 0; graph_pos < count_of(read_graphs);
             graph_pos++) {
          int32_t graph = read_graphs[graph_pos];
          std::vector<int32_t>& nodes = read_nodes[graph_pos];

          // Set read path_bits
          build_read_path_bits(placement, graph, nodes, read_path_bits);

          // Look for the mate's nodes in the same graph.
          int32_t p = 0;
          while (p < count_of(mate_graphs) && graph != mate_graphs[p]) p++;
          // We found a node that is within the same graph -> can only happen if
          // pair read is present;
          if (p < count_of(mate_graphs)) {
            mates_used++;

            // Set pair path_bits
            build_mate_path_bits(placement, graph, mate_nodes[p],
                                 mate_path_bits);

            // Check if there is a conflict of path_bits
            CGraphnode *gnode = tables.splice_graphs[strand][graph]
                                    .nodes[mate_nodes[p][0]].get();
            BitVec conflict_path_bits = gnode->parent_reachable_bits;
            conflict_path_bits[mate_nodes[p][0]] = 1;

            // There isn't a conflict -> pair parents should contain read
            // path_bits
            if (conflict_path_bits.contains(read_path_bits)) {
              int32_t i = 0;
              // Read and pair share a node
              if (mate_nodes[p][0] == nodes.back())
                i++;
              while (i < count_of(mate_nodes[p])) {
                nodes.push_back(mate_nodes[p][i]);
                i++;
              }
              read_path_bits = read_path_bits | mate_path_bits;
              add_path_bits(placement, graph, read_path_bits, nodes);
            } else { // update both path_bits separately
              add_path_bits(placement, graph, read_path_bits, nodes);
              add_path_bits(placement, mate_graphs[p], mate_path_bits,
                          mate_nodes[p]);
            }
            mate_graphs[p] = -1;
          } else { // pair has no valid path_bits in this graph
            add_path_bits(placement, graph, read_path_bits, nodes);
          }
        }

        // Add the mate path_bits in graphs where the read has no nodes.
        if (mates_used < count_of(mate_graphs))
          for (int32_t p = 0; p < count_of(mate_graphs); p++)
            if (mate_graphs[p] > -1) {
              // Set pair path_bits
              build_mate_path_bits(placement, mate_graphs[p], mate_nodes[p],
                                 mate_path_bits);
              add_path_bits(placement, mate_graphs[p], mate_path_bits,
                          mate_nodes[p]);
            }
      }
  }

  // Indexes the transfrags (except those from the source) by their path, to
  // find one without a linear scan (see make_path_key).
  TransfragIndex build_transfrag_index(SpliceGraph& graph) {
    TransfragIndex index;
    index.reserve(graph.transfrag.size());

    for (int32_t t = 0; t < count_of(graph.transfrag); t++)
      // Don't include transfrags from source -> not needed
      if (graph.transfrag[t]->nodes[0]) {
        // The nodes that are set in the path_bits, in order. If links from
        // source to nodes are desired source==1 and all nodes are treated as
        // +1
        std::vector<int32_t> nodes;
        for (int32_t node_id = 1; node_id < graph.node_count; node_id++) {
          if (graph.transfrag[t]->path_bits[node_id]) nodes.push_back(node_id);
        }
        index[make_path_key(graph.path_bit_layout, nodes,
                            graph.transfrag[t]->path_bits)]
                                = graph.transfrag[t].get();
      }

    return index;
  }

} // namespace ribbon
