#include "ribbon.h"
#include "GThreads.h"

// Anchor length for junction to be considered well supported
extern uint32_t junction_support;
extern bool trim_by_coverage;
extern int32_t allowed_nodes;

// Builds the splice graph of one pre-graph: nodes between junctions, the
// source and sink links, the guide boundaries and the parent/child bit vectors.

namespace ribbon {

  // Allocates a new CGraphnode(start,end,node_no) for the graph being built
  // and registers it in both the graph's nodes and in the segment's entry of
  // pre_graph_segment_to_graph_nodes (context.segment_nodes).
  CGraphnode *add_graph_node(const GraphBuildContext& context, uint32_t start,
                             uint32_t end, int32_t node_id,
                             PreGraphSegment *pre_graph_segment) {
    CGraphinfo graph_info(context.graph_idx, node_id);
    (*context.segment_nodes)[pre_graph_segment->bid].push_back(graph_info);
    GraphNodeList& graph_nodes = context.graph->nodes;
    graph_nodes.push_back(std::make_unique<CGraphnode>(context.cell_count,
                                                       start, end, node_id));

    return (graph_nodes.back().get());
  }

  // Splits `graph_node` at a guide transcript's start (`new_start`) and links
  // the new downstream node directly from source, using the coverage jump
  // around new_start to estimate the new source->guide transfrag's abundance.
  CGraphnode *link_source_to_guide_start(const GraphBuildContext& context,
      uint32_t new_start, uint32_t new_end, CGraphnode *graph_node,
      PendingTransfragList& pending_transfrags, int32_t& node_count,
      PreGraphSegment *pre_graph_segment, int32_t& edge_count) {
    // The new_start is very close to graph_node start
    if (graph_node->start + long_intron_anchor > new_start) {
      for (int32_t p = 0; p < count_of(graph_node->parent); p++)
        if (!graph_node->parent[p]) return (graph_node);
    }

    // Estimate the abundance of the source link from the coverage jump at
    // new_start.
    double coverage_left = 0;
    double coverage_right = 0;

    if (new_start > graph_node->start) {
      uint32_t window_start = graph_node->start;
      if (new_start - window_start > CHI_WIN) {
        window_start = new_start - CHI_WIN;
      }

      coverage_left = get_stranded_coverage(context.strand,
                                            window_start - context.ref_start,
                                            new_start - 1 - context.ref_start,
                                            *context.bp_coverage);
      coverage_left /= new_start - window_start;
    }
    if (new_start < new_end) {
      uint32_t window_end = new_end;
      if (new_end - new_start >= CHI_WIN) {
        window_end = new_start + CHI_WIN - 1;
      }

      coverage_right = get_stranded_coverage(context.strand,
                                             new_start - context.ref_start,
                                             window_end - context.ref_start,
                                             *context.bp_coverage);
      coverage_right /= window_end - new_start + 1;
    }

    double max_abundance = coverage_right - coverage_left;
    if (max_abundance < min_transfrag_abundance)
      max_abundance = min_transfrag_abundance;

    // Split the node at new_start if there is a part before it.
    if (graph_node->start <= new_start - 1) {
      uint32_t prev_end = graph_node->end;
      graph_node->end = new_start - 1;
      CGraphnode *prev_node = graph_node;
      graph_node = add_graph_node(context, new_start, prev_end, node_count,
                                  pre_graph_segment);
      node_count++;

      pending_transfrags.push_back({float(prev_node->node_id),
                                    float(graph_node->node_id),
                                    min_transfrag_abundance});
      // This node is the child of previous node
      prev_node->child.push_back(graph_node->node_id);
      // This node has as parent the previous node
      graph_node->parent.push_back(prev_node->node_id);
    }
    // This node is the child of source
    context.source->child.push_back(graph_node->node_id);
    // This node has source as parent
    graph_node->parent.insert(graph_node->parent.begin(),
                              context.source->node_id);
    pending_transfrags.push_back({0.0f, float(node_count - 1),
                                  float(max_abundance), true});
    // COUNT 1 EDGE HERE because the source to guide edge was already included
    // in our count
    edge_count++;

    return (graph_node);
  }

  // Splits `graph_node` at a guide transcript's end (`new_start`) and links the
  // upstream part directly to sink, using the coverage drop around new_start to
  // estimate the new guide->sink transfrag's abundance.
  CGraphnode *link_guide_end_to_sink(const GraphBuildContext& context,
      uint32_t new_start, uint32_t new_end, CGraphnode *graph_node,
      PendingTransfragList& pending_transfrags, int32_t& node_count,
      PreGraphSegment *pre_graph_segment, int32_t& edge_count) {
    // Estimate the abundance of the sink link from the coverage drop at
    // new_start.
    double coverage_left = 0;
    double coverage_right = 0;

    if (new_start >= graph_node->start) {
      uint32_t window_start = graph_node->start;
      if (new_start - window_start >= CHI_WIN) {
        window_start = new_start - CHI_WIN + 1;
      }

      coverage_left = get_stranded_coverage(context.strand,
                                            window_start - context.ref_start,
                                            new_start - context.ref_start,
                                                *context.bp_coverage);
      coverage_left /= new_start - window_start + 1;
    }
    if (new_start + 1 < new_end) {
      uint32_t window_end = new_end;
      if (new_end - new_start > CHI_WIN) {
        window_end = new_start + CHI_WIN;
      }

      coverage_right = get_stranded_coverage(context.strand,
                                             new_start + 1 - context.ref_start,
                                             window_end - context.ref_start,
                                             *context.bp_coverage);
      coverage_right /= window_end - new_start;
    }

    double max_abundance = coverage_left - coverage_right;
    if (max_abundance < min_transfrag_abundance)
      max_abundance = min_transfrag_abundance;

    // Split the node after new_start and link the left part to the sink.
    float prev_node_id = node_count - 1;
    uint32_t prev_end = graph_node->end;
    graph_node->end = new_start;
    CGraphnode *prev_node = graph_node;
    graph_node = add_graph_node(context, new_start + 1, prev_end, node_count,
                                pre_graph_segment);
    node_count++;
    // This node is the child of previous node
    prev_node->child.push_back(graph_node->node_id);
    // This node has as parent the previous node
    graph_node->parent.push_back(prev_node->node_id);
    // The prev_node is the parent of sink
    context.sink->parent.push_back(prev_node->node_id);

    pending_transfrags.push_back({prev_node_id, -1.0f, float(max_abundance),
                                  true});
    pending_transfrags.push_back({float(prev_node->node_id),
                                  float(graph_node->node_id),
                                  min_transfrag_abundance});
    // COUNT 1 EDGE HERE because the source to guide edge was already included
    // in our count
    edge_count++;

    return (graph_node);
  }

  // Makes child a child of parent and counts the edge between them.
  void link_nodes(CGraphnode *parent, CGraphnode *child, int32_t& edge_count) {
    parent->child.push_back(child->node_id);
    child->parent.push_back(parent->node_id);
    edge_count++;
  }

  // Walks the graph depth-first from node, passing down the edges of its
  // ancestors (parents) and collecting child_reachable_bits and
  // parent_reachable_bits on every visited node; adds source/sink transfrags
  // where needed. Returns reachable edges.
  BitVec build_reachable_bits_depth_first(
      const GraphBuildContext& context,
      CGraphnode *node, BitVec parents,
      int32_t node_count,
      std::vector<bool>& visited,
      int32_t& edge_count) {
    PathBitLayout& layout = context.graph->path_bit_layout;
    TransfragList& transfrag = context.graph->transfrag;
    // This graph's nodes, indexed by node id
    GraphNodeList& nodes
        = context.graph->nodes;

    if (visited[node->node_id]) {
      node->parent_reachable_bits = node->parent_reachable_bits | parents;
      for (int32_t other_id = 0; other_id < node_count; other_id++) {
        // Add the node's children to all parents of the node
        if (parents[other_id])
          nodes[other_id]->child_reachable_bits
              = nodes[other_id]->child_reachable_bits
                | node->child_reachable_bits;
        else if (node->child_reachable_bits[other_id])
          nodes[other_id]->parent_reachable_bits
              = nodes[other_id]->parent_reachable_bits
                | node->parent_reachable_bits;
      }
    } else {
      // First visit: record the node's ancestors, add any missing source and
      // sink transfrags, then visit its children.
      node->child_reachable_bits.resize(node_count + edge_count);
      node->parent_reachable_bits.resize(node_count + edge_count);
      node->parent_reachable_bits = node->parent_reachable_bits | parents;
      visited[node->node_id] = true;
      parents[node->node_id] = 1; // add the node to the parents

      // Node has source only as parent -> add transfrag from source to node
      if (count_of(node->parent) == 1 && !node->parent[0]) {
        BitVec edge_path_bits(node_count + edge_count);
        edge_path_bits[0] = 1;
        edge_path_bits[node->node_id] = 1;

        edge_path_bits[layout.get_or_create_edge_bit(0, node->node_id)] = 1;

        std::vector<int32_t> nodes;
        nodes.push_back(0);
        nodes.push_back(node->node_id);

        transfrag.push_back(
            std::make_unique<CTransfrag>(context.cell_count, nodes,
                                         edge_path_bits,
                                         min_transfrag_abundance));
      }

      // A node without children gets the sink as its only child.
      int32_t child_count = count_of(node->child);
      if (node != context.sink && !child_count) {
        // Add the sink to the node's children and the node to the sink's
        // parents.
        node->child.push_back(context.sink->node_id);
        context.sink->parent.push_back(node->node_id);
        // Create the transfrag that ends the node
        BitVec edge_path_bits(node_count + edge_count);
        edge_path_bits[node->node_id] = 1;
        edge_path_bits[node_count - 1] = 1;

        edge_path_bits[layout.get_or_create_edge_bit(node->node_id,
                                                   node_count - 1)] = 1;

        std::vector<int32_t> nodes;
        nodes.push_back(node->node_id);
        nodes.push_back(context.sink->node_id);

        transfrag.push_back(
            std::make_unique<CTransfrag>(context.cell_count, nodes,
                                         edge_path_bits,
                                         min_transfrag_abundance));
        child_count++;
      }

      // For each child, mark the edge for the child's ancestors and for this
      // node's descendants, then recurse into the child.
      for (int32_t i = 0; i < child_count; i++) { // for all children
        BitVec child_parents = parents;
        // Nodeid is always smaller than child node ?
        int32_t min_node = node->node_id;
        int32_t max_node = node->child[i];
        if (min_node > max_node) {
          max_node = node->node_id; // node_id is always smaller than child node ?
          min_node = node->child[i];
        }

        int32_t edge_bit = layout.get_or_create_edge_bit(min_node, max_node);
        // Add edge from node to child to the set of parents from child
        child_parents[edge_bit] = 1;
        // Add edge from node to child to the set of node children
        node->child_reachable_bits[edge_bit] = 1;

        node->child_reachable_bits = node->child_reachable_bits
            | build_reachable_bits_depth_first(context,
                nodes[node->child[i]].get(), child_parents, node_count,
                visited, edge_count);
      }
    }

    // Return the edges reachable from this node, including itself.
    BitVec children = node->child_reachable_bits;
    children[node->node_id] = 1;

    return (children);
  }

  // Splits graph_node at the guide starts and ends of this strand that lie
  // before limit, advancing guide_boundary_idx past them. A guide start links
  // its part of the node to the source, a guide end links the part before it
  // to the sink. Returns the last node.
  CGraphnode *split_node_at_guide_boundaries(
      const GraphBuildContext& context,
      std::vector<GuideBoundary>& guide_boundaries,
      int32_t& guide_boundary_idx, uint32_t limit, uint32_t curr_start,
      uint32_t pre_graph_segment_end, CGraphnode *graph_node,
      PendingTransfragList& pending_transfrags, int32_t& node_count,
      PreGraphSegment *pre_graph_segment, int32_t& edge_count) {
    int32_t boundary_count = count_of(guide_boundaries);

    // Skip the boundaries that are not after the start of the node.
    while (guide_boundary_idx < boundary_count
           && guide_boundaries[guide_boundary_idx].pos
           <= graph_node->start) guide_boundary_idx++;

    while (true) {
      while (guide_boundary_idx < boundary_count
             && guide_boundaries[guide_boundary_idx].strand
             != context.strand) guide_boundary_idx++;

      if (guide_boundary_idx >= boundary_count
          || guide_boundaries[guide_boundary_idx].pos >= limit) break;

      uint32_t segment_start = guide_boundaries[guide_boundary_idx].pos;
      uint32_t segment_end = limit;
      bool is_guide_start = guide_boundaries[guide_boundary_idx].pos
          < guide_boundaries[guide_boundary_idx].exon_other_end;
      guide_boundary_idx++;
      // Ignore a guide that does not reach into this segment.
      if (is_guide_start) {
        if (guide_boundaries[guide_boundary_idx - 1].exon_other_end
            > pre_graph_segment_end)
          continue;
      } else if (guide_boundaries[guide_boundary_idx - 1].exon_other_end
                 < curr_start)
        continue;

      // The part ends at the next boundary of this strand, if it is before the
      // limit.
      while (guide_boundary_idx < boundary_count
             && guide_boundaries[guide_boundary_idx].strand
             != context.strand) guide_boundary_idx++;
      if (guide_boundary_idx < boundary_count
          && guide_boundaries[guide_boundary_idx].pos < limit)
        segment_end = guide_boundaries[guide_boundary_idx].pos;

      if (is_guide_start) {
        graph_node = link_source_to_guide_start(context, segment_start,
            segment_end, graph_node, pending_transfrags, node_count,
            pre_graph_segment, edge_count);
      } else {
        graph_node = link_guide_end_to_sink(context, segment_start,
            segment_end, graph_node, pending_transfrags, node_count,
            pre_graph_segment, edge_count);
      }
    }

    return graph_node;
  }

  // Links the nodes to the source or sink where the coverage jumps up or drops
  // sharply compared to the nodes around them, and counts the edges of the
  // nodes without children (build_reachable_bits_depth_first links those to
  // the sink).
  void link_weak_nodes_to_source_and_sink(
      const GraphBuildContext& context, int32_t node_count,
      int32_t& edge_count, PendingTransfragList& pending_transfrags) {
    char strand = context.strand;
    int32_t ref_start = context.ref_start;
    BpCoverage& bp_coverage = *context.bp_coverage;
    GraphNodeList& nodes
        = context.graph->nodes;
    CGraphnode *source = context.source;
    CGraphnode *sink = context.sink;

    // The stranded coverage per bp of node idx.
    auto coverage_per_bp = [&](int32_t idx) {
      return get_stranded_coverage(strand, nodes[idx]->start - ref_start,
                                   nodes[idx]->end - ref_start, bp_coverage)
          / nodes[idx]->len();
    };

    for (int32_t i = 1; i < node_count; i++) {
      double node_coverage = 0;
      // Node i has parents, and does not come from source => might need to be
      // linked to source
      if (i > 1 && nodes[i]->parent[0]) {
        node_coverage = coverage_per_bp(i);
        double parents_coverage = 0; // all parents' coverage
        bool link_to_source = true;
        for (int32_t j = 0; j < count_of(nodes[i]->parent); j++) {
          if (!nodes[i]->parent[j]) { // parent is source
            link_to_source = false;
            break;
          }
          parents_coverage += coverage_per_bp(nodes[i]->parent[j]);
        }
        if (link_to_source
            && parents_coverage < node_coverage * ERROR_PERC * DROP) {
          nodes[i]->parent.insert(nodes[i]->parent.begin(), source->node_id);
          source->child.push_back(i);
          pending_transfrags.push_back({0.0f, float(i),
              float((node_coverage - parents_coverage) / DROP)});
          edge_count++;
        }
      }

      if (i < node_count - 1) {
        if (count_of(nodes[i]->child)) { // might need to be linked to sink
          bool link_to_sink = true;
          for (int32_t j = 0; j < count_of(sink->parent); j++) {
            if (sink->parent[j] == i) {
              link_to_sink = false;
              break;
            }
          }
          if (link_to_sink) {
            if (!node_coverage) node_coverage = coverage_per_bp(i);
            double children_coverage = 0; // all children coverages
            for (int32_t j = 0; j < count_of(nodes[i]->child); j++) {
              children_coverage += coverage_per_bp(nodes[i]->child[j]);
            }
            if (children_coverage < node_coverage * ERROR_PERC * DROP) {
              sink->parent.push_back(i); // prev_node is the parent of sink
              pending_transfrags.push_back({float(i), -1.0f,
                  float((node_coverage - children_coverage) / DROP)});
              edge_count++;
            }
          }
        // Node has no children so it might get linked to sink in
        // build_reachable_bits_depth_first function
        } else edge_count++;
      }
    }
  }

  // Turns the pending links into real transfrags, skipping the removed ones.
  // node_count is final: the sink is the last node.
  void add_pending_transfrags(const GraphBuildContext& context,
                              int32_t node_count, int32_t edge_count,
                              PendingTransfragList& pending_transfrags) {
    PathBitLayout& layout = context.graph->path_bit_layout;
    int32_t cell_count = context.cell_count;
    GraphNodeList& nodes = context.graph->nodes;
    CGraphnode *sink = context.sink;

    for (const PendingTransfrag& pending : pending_transfrags) {
      // Add links between node and sink
      int32_t from_node = int32_t(pending.from);
      if (from_node >= 0) {
        int32_t to_node = int32_t(pending.to);
        BitVec edge_path_bits(node_count + edge_count);
        edge_path_bits[from_node] = 1;
        std::vector<int32_t> transfrag_nodes;
        if (to_node < 0) {
          CGraphnode *node = nodes[from_node].get();

          int32_t scan_idx = from_node + 1;
          if (scan_idx < sink->node_id) {
            CGraphnode *prev_node = node;
            // This should be okay
            CGraphnode *next_node = nodes[scan_idx].get();
            uint32_t dist = next_node->len() - 1;
            bool keep_transfrag = true;
            while (next_node->start == prev_node->end + 1
                   && dist < long_intron_anchor) {
              for (int32_t c = 0; c < count_of(next_node->child); c++)
                if (next_node->child[c] == sink->node_id) {
                  // Do not continue with the pending_transfrags
                  keep_transfrag = false;
                  break;
                }
              scan_idx++;
              if (scan_idx == sink->node_id) {
                // No link to sink found -> can continue with pending_transfrags
                break;
              }
              prev_node = next_node;
              next_node = nodes[scan_idx].get();
              dist += next_node->len();
            }
            // Found link to sink; no need to keep pending_transfrags
            if (!keep_transfrag) continue;
          }

          node->child.push_back(sink->node_id);
          // Add node to sink transfrag
          edge_path_bits[node_count - 1] = 1;

          edge_path_bits[layout.get_or_create_edge_bit(from_node,
                                                     node_count - 1)] = 1;

          transfrag_nodes.push_back(from_node);
          transfrag_nodes.push_back(sink->node_id);
        } else {
          edge_path_bits[to_node] = 1;
          edge_path_bits[layout.get_or_create_edge_bit(from_node, to_node)] = 1;

          transfrag_nodes.push_back(from_node);
          transfrag_nodes.push_back(to_node);
        }
        TransfragList& transfrag = context.graph->transfrag;
        transfrag.push_back(
            std::make_unique<CTransfrag>(cell_count, transfrag_nodes,
                                         edge_path_bits,
                                         pending.abundance));
        transfrag.back()->made_by_guide = pending.made_by_guide;
      }
    }
  }

  // Checks if any guide boundary (a guide start or end, sorted by position)
  // on this strand falls inside the pre-graph; if so remembers it and rewinds
  // guide_boundary_idx, which the walk advances past handled edges.
  void SpliceGraphBuilder::find_guide_boundary_in_pre_graph() {
    PreGraphSegment *segment = pre_graph_segments[pre_graph->start_node].get();
    while (guide_boundary_idx < count_of(guide_boundaries)
           && segment != NULL) {
      uint32_t start = guide_boundaries[guide_boundary_idx].pos;
      uint32_t end = guide_boundaries[guide_boundary_idx].exon_other_end;
      if (start > end) std::swap(start, end);
      if (segment->end < start) segment = segment->next_node;
      // The segment->end>=start
      else if (guide_boundaries[guide_boundary_idx].strand == strand
               && segment->start <= end) {
        guide_boundary_idx = 0;
        pre_graph_has_guide_boundary = true;
        break;
      } else guide_boundary_idx++;
    }
  }

  // True if the next junction end (start) is at or before pos.
  bool SpliceGraphBuilder::next_junction_end_is_by(uint32_t pos) {
    return junc_end_idx < n_junctions
           && state.end_sorted_junctions[junc_end_idx]->end <= pos;
  }

  bool SpliceGraphBuilder::next_junction_start_is_by(uint32_t pos) {
    return junc_start_idx < n_junctions
           && junction[junc_start_idx]->start <= pos;
  }

  // Reads all junction ends at or before position -> assuming there are any
  // (at this point, smaller junction ends should not be relevant to the
  // pre-graph segment). Returns whether a junction on this strand ends exactly
  // at position.
  bool SpliceGraphBuilder::read_junction_ends_up_to(uint32_t position) {
    bool ends_here = false;
    while (junc_end_idx < n_junctions
           && state.end_sorted_junctions[junc_end_idx]->end <= position) {
      if (state.end_sorted_junctions[junc_end_idx]->end == position
          && state.end_sorted_junctions[junc_end_idx]->strand
              == context.strand) {
        ends_here = true;
      }
      junc_end_idx++;
    }
    return ends_here;
  }

  // Skips the junction ends (starts) that are not on this strand.
  void SpliceGraphBuilder::skip_ends_on_other_strands() {
    while (junc_end_idx < n_junctions
           && state.end_sorted_junctions[junc_end_idx]->strand
              != context.strand) junc_end_idx++;
  }

  void SpliceGraphBuilder::skip_starts_on_other_strands() {
    while (junc_start_idx < n_junctions
           && junction[junc_start_idx]->strand != context.strand)
      junc_start_idx++;
  }

  // Links node to the nodes that have a junction ending at pos; returns
  // whether there are any.
  bool SpliceGraphBuilder::link_from_junction_origins(uint32_t pos,
                                                      CGraphnode *node) {
    auto origin_it = junction_origins.find(pos);
    if (origin_it == junction_origins.end()) return false;
    for (int32_t i = 0; i < count_of(origin_it->second); i++) {
      link_nodes(nodes[origin_it->second[i]].get(), node, edge_count);
    }
    return true;
  }

  // Makes a node that starts at start and ends at the segment's end, and
  // links the open node to it; the new node becomes the open node.
  void SpliceGraphBuilder::open_next_node(uint32_t start) {
    CGraphnode *next_node = add_graph_node(context, start,
                                           pre_graph_segment_end, node_count,
                                           pre_graph_segment);
    node_count++;
    link_nodes(graph_node, next_node, edge_count);
    graph_node = next_node;
  }

  // Splits the open node at the guide boundaries and at the coverage drops
  // before end_pos, and ends it there. The splits may already have added
  // smaller nodes, and the open node becomes the last of them.
  void SpliceGraphBuilder::close_open_node_at(uint32_t end_pos) {
    // Add the guide starts and ends before end_pos.
    if (pre_graph_has_guide_boundary)
      graph_node = split_node_at_guide_boundaries(context, guide_boundaries,
          guide_boundary_idx, end_pos, curr_start, pre_graph_segment_end,
          graph_node, pending_transfrags, node_count, pre_graph_segment,
          edge_count);

    // Find intermediate nodes at coverage drops; alternatively this could be
    // done only for end nodes
    if (trim_by_coverage)
      graph_node = trim_node_at_all_points(context, end_pos, graph_node,
                                           pending_transfrags, node_count,
                                           pre_graph_segment, edge_count);

    graph_node->end = end_pos;
  }

  // Moves curr_start past a little hanging piece before a junction end that
  // is close to it, and reads the junction ends up to curr_start. Returns
  // whether a junction on this strand ends exactly at curr_start.
  bool SpliceGraphBuilder::skip_hanging_piece_at_segment_start() {
    // True if a junction on this strand ends exactly at curr_start
    bool junction_ends_here = read_junction_ends_up_to(curr_start);

    // See if we need to adjust the start to ignore little hanging pieces that
    // make no sense
    if (!junction_ends_here) {
      skip_ends_on_other_strands();
      // There is a junction ending soon here
      if ((junc_end_idx < n_junctions
           && state.end_sorted_junctions[junc_end_idx]->end - curr_start
           < junction_support)) {
        double coverage_left = get_region_coverage('.', curr_start - ref_start,
            state.end_sorted_junctions[junc_end_idx]->end - 1 - ref_start,
            bp_coverage);
        double coverage_right = get_region_coverage('.',
            state.end_sorted_junctions[junc_end_idx]->end - ref_start,
            2 * state.end_sorted_junctions[junc_end_idx]->end - curr_start - 1
            - ref_start, bp_coverage);
        // Adjust start here if needed
        if (coverage_left < coverage_right * (1 - ERROR_PERC)) {
          curr_start = state.end_sorted_junctions[junc_end_idx]->end;
          // We have to check ending junctions here again
          if (read_junction_ends_up_to(curr_start)) junction_ends_here = true;
        }
      }
    }
    return junction_ends_here;
  }

  // Creates the first node of the segment, numbered node_count, and links it
  // to the nodes that have a junction ending here, or to the source if there
  // are none.
  void SpliceGraphBuilder::open_first_node(bool junction_ends_here) {
    graph_node = add_graph_node(context, curr_start, pre_graph_segment_end,
                                node_count, pre_graph_segment);
    node_count++;

    if (junction_ends_here) {
      // We haven't seen nodes before that finish here (maybe due to error
      // correction?) => link to source
      if (!link_from_junction_origins(curr_start, graph_node))
        link_nodes(source, graph_node, edge_count);
    } else { // this node comes from source directly
      link_nodes(source, graph_node, edge_count);
    }
  }

  // Decides what to process next in the segment: 0 = a junction start, 1 = a
  // junction end, -1 = nothing left. Whichever comes first on the genome; on
  // a tie the end is processed first.
  int32_t SpliceGraphBuilder::next_junction_event() {
    int32_t event = -1;
    if (next_junction_end_is_by(pre_graph_segment_end)
        || next_junction_start_is_by(pre_graph_segment_end)) {
      if (junc_start_idx < n_junctions
          && (junction[junc_start_idx]->start <= pre_graph_segment_end)
          && junction[junc_start_idx]->end > pre_graph_end) junc_start_idx++;
      else {
        // There are still junctions endings
        if (junc_end_idx < n_junctions) {
          // There are still junctions starting
          if (junc_start_idx < n_junctions) {
            // One of them is clearly before the pre_graph_segment_end from the
            // initial if
            event = junction[junc_start_idx]->start
                >= state.end_sorted_junctions[junc_end_idx]->end ? 1 : 0;
          } else event = 1;
        } else event = 0;
      }
    }
    return event;
  }

  // True if the coverage drops after pos, in the short tail that is left of
  // the segment.
  bool SpliceGraphBuilder::tail_has_coverage_drop(uint32_t pos) {
    double coverage_left = get_region_coverage('.',
        2 * pos - pre_graph_segment_end + 1 - ref_start, pos - ref_start,
        bp_coverage);
    double coverage_right = get_region_coverage('.', pos + 1 - ref_start,
        pre_graph_segment_end - ref_start, bp_coverage);
    return coverage_right < coverage_left * (1 - ERROR_PERC);
  }

  // Closes the open node at a junction start and remembers where the
  // junction ends lead. If the segment goes on after it, opens the next node,
  // unless it is a short tail with a coverage drop.
  void SpliceGraphBuilder::process_junction_start() {
    uint32_t pos = junction[junc_start_idx]->start;
    close_open_node_at(pos);

    // Remember ends here
    while (junc_start_idx < n_junctions
           && junction[junc_start_idx]->start == pos) {
      if (junction[junc_start_idx]->strand == context.strand) {
        junction_origins[junction[junc_start_idx]->end].push_back(
            graph_node->node_id);
      }
      junc_start_idx++;
    }

    if (pos < pre_graph_segment_end) {
      // See if we should just skip node
      if (pre_graph_segment_end - pos < junction_support) {
        skip_starts_on_other_strands();
        // There are no more junctions starting within this pre-graph
        if (!next_junction_start_is_by(pre_graph_segment_end)
            && !next_junction_end_is_by(pre_graph_segment_end)) {
          // Adjust start here if needed
          if (tail_has_coverage_drop(pos)) pre_graph_segment_closed = true;
        }
      }
      if (!pre_graph_segment_closed) open_next_node(pos + 1);
    } else pre_graph_segment_closed = true;
  }

  // Reads the junction ends at the next end position and links the nodes that
  // have a junction ending there; the open node is split at that position if
  // it starts before it.
  void SpliceGraphBuilder::process_junction_end() {
    uint32_t pos = state.end_sorted_junctions[junc_end_idx]->end;
    // Read all junction ends at the current start
    while (junc_end_idx < n_junctions
           && state.end_sorted_junctions[junc_end_idx]->end == pos) {
      junc_end_idx++;
    }

    // Last created node starts before the position of the new node we want
    // to create
    if (graph_node->start < pos) {
      close_open_node_at(pos - 1);
      open_next_node(pos);
    }

    link_from_junction_origins(pos, graph_node);
  }

  // Walks the junction events inside the segment, closing and opening
  // nodes, and ends the last node at the segment's end with an edge to the
  // sink.
  void SpliceGraphBuilder::process_segment() {
    curr_start = pre_graph_segment->start;
    pre_graph_segment_end = pre_graph_segment->end;
    bool junction_ends_here = skip_hanging_piece_at_segment_start();
    open_first_node(junction_ends_here);

    pre_graph_segment_closed = false;
    do {
      skip_ends_on_other_strands();
      // Junctions that start before the current graph_node and we haven't
      // seen them before are part of a different pre-graph
      while (junc_start_idx < n_junctions
             && (junction[junc_start_idx]->strand != context.strand
                 || junction[junc_start_idx]->start < curr_start))
        junc_start_idx++;

      int32_t event = next_junction_event();
      if (event == 0) process_junction_start();
      else if (event == 1) process_junction_end();
    } while (next_junction_end_is_by(pre_graph_segment_end)
             || next_junction_start_is_by(pre_graph_segment_end));

    // The node is still open, so it ends here.
    if (!pre_graph_segment_closed) {
      close_open_node_at(pre_graph_segment_end);
      // COUNT EDGE HERE (this is an edge to sink)
      edge_count++;
    }
  }

  // Completes the graph once all segments are walked: links weak nodes,
  // prunes if there are too many nodes, adds the sink, and builds the
  // transfrags and path_bits. Returns the final node count.
  int32_t SpliceGraphBuilder::finish_graph() {
    link_weak_nodes_to_source_and_sink(context, node_count, edge_count,
                                       pending_transfrags);

    // The graph is complete, so node_count is known: prune it if it has too
    // many nodes (allowed_nodes).
    if (node_count > allowed_nodes) {
      node_count = prune_graph_to_node_limit(context, node_count,
                                             count_of(pre_graph_segments),
                                             junction, edge_count,
                                             pending_transfrags);
    }

    // Add the sink as the last node.
    sink->node_id = node_count;
    nodes.push_back(std::move(sink_owner));
    node_count++;

    // Nodes are from 0 to node_count-1, so the first "available" position in
    // BitVec is node_count
    graph.path_bit_layout.init(node_count);

    // Now that node_count is final, turn the pending links into real
    // transfrags (skipping the removed ones).
    add_pending_transfrags(context, node_count, edge_count,
                           pending_transfrags);

    // Finished reading the pre-graph: now build the parent and child path_bits.
    std::vector<bool> visited(node_count, false);
    BitVec parents(node_count + edge_count);

    build_reachable_bits_depth_first(context, source, parents, node_count,
                                     visited, edge_count);

    junction_origins.clear();

    return node_count;
  }

  // Sets up the source and sink nodes and checks whether any guide boundary
  // falls inside the pre-graph.
  SpliceGraphBuilder::SpliceGraphBuilder(
      BundleData *bundle_in, AssemblyState& state_in, GraphTables& tables_in,
      char strand_in, int32_t graph_idx_in,
      std::vector<GuideBoundary>& guide_boundaries_in)
      : bundle(bundle_in), state(state_in),
        graph(tables_in.splice_graphs[strand_in][graph_idx_in]),
        strand(strand_in), graph_idx(graph_idx_in),
        edge_count(graph.edge_count),
        guide_boundaries(guide_boundaries_in),
        ref_start(bundle_in->start), junction(bundle_in->junction),
        n_junctions(count_of(bundle_in->junction)),
        pre_graph(state_in.pre_graph_list[strand_in][graph_idx_in].get()),
        pre_graph_segments(state_in.pre_graph_segment_list[strand_in]),
        pre_graph_end(pre_graph_segments[pre_graph->last_node_id]->end),
        nodes(graph.nodes),
        // We might want to use a different type of data for bp_coverage to
        // save memory in the case of very long bundles
        bp_coverage(bundle_in->bp_coverage),
        cell_count(count_of(bundle_in->cell_name)) {
    nodes.push_back(std::make_unique<CGraphnode>(cell_count, 0, 0, 0));
    source = nodes.back().get();
    sink_owner = std::make_unique<CGraphnode>(cell_count, 0, 0, 0);
    sink = sink_owner.get();
    context = {strand, graph_idx, ref_start, cell_count, &bp_coverage, &graph,
               &tables_in.pre_graph_segment_to_graph_nodes[strand_in], source,
               sink};
    find_guide_boundary_in_pre_graph();
  }

  // Goes through the pre-graph segments left to right, creates the graph
  // nodes, and completes the graph.
  int32_t SpliceGraphBuilder::build() {
    pre_graph_segment = pre_graph_segments[pre_graph->start_node].get();
    while (pre_graph_segment != NULL) {
      process_segment();
      // advance to the next segment
      pre_graph_segment = pre_graph_segment->next_node;
    }
    return finish_graph();
  }

  // Builds the splice graph of one pre-graph on a strand: walks its segments
  // left to right, makes a node per stretch between junctions, splits them at
  // guide ends and trim points, prunes the graph and builds its path_bits.
  int32_t build_splice_graph(BundleData *bundle, AssemblyState& state,
                         GraphTables& tables, char strand, int32_t graph_idx,
                         std::vector<GuideBoundary>& guide_boundaries) {
    SpliceGraphBuilder builder(bundle, state, tables, strand, graph_idx,
                               guide_boundaries);
    return builder.build();
  }

} // namespace ribbon
