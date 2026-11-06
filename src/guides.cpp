#include "ribbon.h"
#include "GThreads.h"

// Matches the reference (guide) transcripts against a bundle's splice graph.
// match_guide_to_graph finds a guide's exact exon/intron path among the nodes;
// process_reference_guides matches and links every usable guide;
// is_within_guide_exon tests an interval against the guide exons. The guide's
// bundle status and its boundaries inside a pre-graph are also kept here.

namespace ribbon {

  // 2-bit mask of the bundle status in a guide's user flags (the values are
  // GuideBundleStatus, ribbon.h)
  static const byte tx_status_mask = 0x03;

  // Returns the 2-bit bundle-processing status stored in guide's user flags:
  // unused, in-bundle-but-not-fully-covered, all-introns-covered, or
  // stored-for-printing.
  GuideBundleStatus get_guide_status(GffObj *guide) {
    return (GuideBundleStatus)guide->getUserFlags(tx_status_mask);
  }

  // Sets guide's 2-bit bundle-processing status (see get_guide_status).
  void set_guide_status(GffObj *guide, GuideBundleStatus status) {
    guide->clearUserFlags(tx_status_mask);
    guide->setUserFlags(tx_status_mask & status); // can be 1, 2 or 0 (clear)
  }

  // True if reads cover every intron of the guide (compute_guide_boundaries),
  // so the guide is used to help build the graph.
  bool is_guide_supported_by_reads(GffObj *guide) {
    return get_guide_status(guide) >= GBST_ALL_INTR_COV;
  }

  // Records the start and end boundaries (guide_boundaries) of every guide
  // whose introns are all covered by reads (marking it GBST_ALL_INTR_COV) and
  // gives guide_prediction (the prediction index of each guide) one slot per
  // guide (-1 until find_transcripts sets it).
  void compute_guide_boundaries(
      std::vector<GffObj *>& guides,
          std::vector<GuideBoundary>& guide_boundaries,
      std::vector<int32_t>& guide_prediction) {
    // Adds a guide boundary to guide_boundaries. If one with the same position
    // and strand exists, keeps the narrower exon (smaller exon_other_end for a
    // start, larger for an end).
    auto add_boundary = [&](uint32_t pos, uint32_t exon_other_end,
                            char strand, bool is_start) {
      GuideBoundary boundary(pos, exon_other_end, strand);
      int32_t idx = find_sorted(guide_boundaries, boundary);
      if (idx < 0) insert_sorted(guide_boundaries, boundary);
      else if (is_start ? guide_boundaries[idx].exon_other_end > exon_other_end
                        : guide_boundaries[idx].exon_other_end < exon_other_end)
        guide_boundaries[idx].exon_other_end = exon_other_end;
    };

    for (int32_t guide = 0; guide < count_of(guides); guide++) {
      guide_prediction.push_back(-1);

      // Check whether every intron of the guide is covered by reads: use
      // TranscriptCounts::transcript_introns[]->read_count.
      bool introns_are_covered = true;
      TranscriptCounts *tdata = (TranscriptCounts *)(guides[guide]->uptr);
      for (int32_t i = 0; i < count_of(tdata->transcript_introns); i++) {
        if (!tdata->transcript_introns[i]->read_count) {
          introns_are_covered = false;
          break;
        }
      }
      if (!introns_are_covered) continue;

      // All introns are covered by long or short reads: add the guide start
      // and end boundaries, on both strands if the guide strand is unknown.
      set_guide_status(guides[guide], GBST_ALL_INTR_COV);
      bool unknown_strand = guides[guide]->strand != '+'
                            && guides[guide]->strand != '-';
      for (char strand : {'-', '+'}) {
        if (!unknown_strand && strand != guides[guide]->strand)
          continue;
        add_boundary(guides[guide]->start, guides[guide]->exons[0]->end,
                     strand, true);
        add_boundary(guides[guide]->end, guides[guide]->exons.Last()->start,
                     strand, false);
      }
    }
  }

  // Moves guide_cursor past the guides that end before the pre-graph starts,
  // then counts the guides from there on that start before the pre-graph ends,
  // are on the strand (or unstranded) and are supported by reads. Each of them
  // may add a source and a sink edge to the graph. The guides must be sorted by
  // start.
  int32_t count_guides_overlapping_pre_graph(
      std::vector<GffObj *>& guides, char strand_char, uint32_t pre_graph_start,
      uint32_t pre_graph_end, int32_t& guide_cursor) {
    int32_t guide_total = count_of(guides);
    while (guide_cursor < guide_total
           && guides[guide_cursor]->end < pre_graph_start) guide_cursor++;

    int32_t overlapping_guides = 0;
    for (int32_t guide_idx = guide_cursor;
         guide_idx < guide_total && guides[guide_idx]->start <= pre_graph_end;
         guide_idx++) {
      bool strand_matches = guides[guide_idx]->strand == strand_char
                            || guides[guide_idx]->strand == '.';
      if (strand_matches && is_guide_supported_by_reads(guides[guide_idx])) {
        overlapping_guides++;
      }
    }
    return overlapping_guides;
  }

  // Follows the guide's exons through the splice graph. Returns a transfrag
  // with the guide's exact path, or an empty pointer if the graph lacks a node
  // or an edge for part of it.
  std::shared_ptr<CTransfrag> match_guide_to_graph(
      GffObj *guide, SpliceGraph& graph, int32_t cell_count) {
    std::shared_ptr<CTransfrag> guide_transfrag;

    // Find the node where the guide starts: the last node of the contiguous,
    // unspliced run of nodes that overlap the guide's first exon.
    int32_t node_id = 1;
    while (node_id < graph.node_count - 1) {
      // First exon of guide overlaps this node
      if (graph.nodes[node_id]->overlap(guide->exons[0])) {
        int32_t run_end = node_id + 1;
        // True if the node overlaps the first exon, touches the previous node
        // and is its first child.
        auto extends_unspliced_run = [&](int32_t node) {
          auto& previous = graph.nodes[node - 1];
          return graph.nodes[node]->overlap(guide->exons[0])
                 && graph.nodes[node]->start == previous->end + 1
                 && count_of(previous->child) && previous->child[0] == node;
        };
        while (run_end < graph.node_count
            - 1 && extends_unspliced_run(run_end)) {
            run_end++;
          }

        // We are skipping nodes that are linked by introns to the last node
        // that is overlapped by the first exon in guide
        if (run_end < graph.node_count - 1
            && graph.nodes[run_end]->overlap(guide->exons[0])) node_id = run_end
            - 1;
        else break;
      }
      node_id++;
    }

    if (node_id < graph.node_count - 1) { // found start node
      // Start the guide's path at that node and follow its exons and introns.
      BitVec guide_path_bits(graph.node_count + graph.edge_count);
      std::vector<int32_t> path_nodes;
      guide_path_bits[node_id] = 1;
      path_nodes.push_back(node_id);
      CGraphnode *current_node = graph.nodes[node_id].get();
      int32_t exon_idx = 0; // the guide exon being matched
      int32_t exon_count = guide->exons.Count();
      while (exon_idx < exon_count) {
        if (guide->exons[exon_idx]->end < current_node->end) {
          // Guide ends inside the node: a match only if it is the last exon.
          // The abundance is estimated later, only if needed.
          if (exon_idx == exon_count - 1)
            guide_transfrag = std::make_shared<CTransfrag>(
                cell_count, path_nodes, guide_path_bits);
          return (guide_transfrag);
        } else if (guide->exons[exon_idx]->end == current_node->end) {
          // Exon ends exactly at the node's end: the guide splices to the next
          // exon here (or ends)
          if (exon_idx == exon_count - 1) {
            guide_transfrag = std::make_shared<CTransfrag>(
                cell_count, path_nodes, guide_path_bits);
            return (guide_transfrag);
          }
          exon_idx++;
          // The next exon must start exactly at a child of this node (a real
          // graph edge)
          int32_t child_pos = 0;
          int32_t child_count = count_of(current_node->child);
          while (child_pos < child_count
                 && graph.nodes[current_node->child[child_pos]]->start
                 != guide->exons[exon_idx]->start) child_pos++;
          // No such edge in the graph: no match
          if (child_pos == child_count) return (guide_transfrag);
          guide_path_bits[current_node->child[child_pos]] = 1;
          graph.path_bit_layout.set_edge_bit(guide_path_bits, node_id,
                              current_node->child[child_pos]);

          node_id = current_node->child[child_pos];
          path_nodes.push_back(node_id);
          current_node = graph.nodes[node_id].get();
        } else {
          // Exon continues past the node's end: it must continue into the
          // adjacent, unspliced next node
          if (node_id < graph.node_count - 1
              && (graph.nodes[node_id + 1]->start == current_node->end + 1)) {
            if (count_of(current_node->child)
                && current_node->child[0] == node_id + 1) {
              guide_path_bits[node_id + 1] = 1;
              graph.path_bit_layout.set_edge_bit(guide_path_bits, node_id,
                  node_id + 1);

              node_id++;
              path_nodes.push_back(node_id);
              current_node = graph.nodes[node_id].get();
            } else return (guide_transfrag);
          } else { // node doesn't continue
            // Abundance is estimated later, only if needed
            if (exon_idx == exon_count - 1)
              guide_transfrag = std::make_shared<CTransfrag>(
                  cell_count, path_nodes, guide_path_bits);
            return (guide_transfrag);
          }
        }
      }
    }
    return (guide_transfrag);
  }

  // Comparator: sorts CGuide by path_bits size (number of graph nodes/edges
  // set), largest first.
  int32_t compare_guide_by_path_bits_size(const void *p1, const void *p2) {
    CGuide *a = (CGuide *)p1;
    CGuide *b = (CGuide *)p2;
    // Longest guide takes precedence
    if (a->transfrag->path_bits.count()
        < b->transfrag->path_bits.count()) return 1;
    if (a->transfrag->path_bits.count()
        > b->transfrag->path_bits.count()) return -1;
    return 0;
  }

  // Links the end of matched guide g to the source (to_source) or the sink if
  // its transfrag does not reach it yet, adding a transfrag whose abundance
  // is the bp coverage step at the guide boundary.
  static void extend_guide_to_terminal(bool to_source, int32_t g,
                                       SpliceGraph& graph,
                                       std::vector<CGuide>& matched_guides,
                                       std::vector<GffObj *>& guides,
                                       BundleData *bundle, int32_t cell_count,
                                       uint32_t ref_start) {
    GffObj *guide = guides[matched_guides[g].g];
    // Find the guide's end node on the side being extended and the terminal
    // (source = 0, sink = last node) to link it to.
    int32_t end_node_id = to_source ? matched_guides[g].transfrag->nodes[0]
        : matched_guides[g].transfrag->nodes.back();
    int32_t terminal_id = to_source ? 0 : graph.node_count - 1;
    // Terminal is already present in guide's path_bits
    if (to_source ? !end_node_id : (end_node_id == terminal_id)) return;

    // Bit position of the edge (end node <-> terminal) in the path_bits;
    // allocate one if the edge is new
    int32_t edge_bit = to_source
        ? graph.path_bit_layout.get_or_create_edge_bit(0, end_node_id)
        : graph.path_bit_layout.get_or_create_edge_bit(end_node_id,
            terminal_id);

    // Check the end node's neighbors: stop if one is the terminal, else
    // remember the coverage on both sides of an unspliced neighbor.
    bool already_linked = false; // assume there is no extension to source/sink
    CGraphnode *end_node = graph.nodes[end_node_id].get();
    // Per-bp coverage just before the guide's boundary
    double left_coverage = 0;
    double right_coverage = 0; // per-bp coverage just after it
    int32_t neighbor = -1;
    std::vector<int32_t>& neighbors = to_source ? end_node->parent
        : end_node->child;
    for (int32_t k = 0; k < count_of(neighbors); k++) {
      if (neighbors[k] == terminal_id) {
        already_linked = true;
        break;
      } // found link to source/sink -> break from path

      bool adjacent = to_source
          ? (graph.nodes[neighbors[k]]->end == end_node->start - 1)
        // No splice site between neighbor and end_node
        : (graph.nodes[neighbors[k]]->start == end_node->end + 1);

      if (adjacent) {
        neighbor = neighbors[k];
        bool in_range = to_source ? (end_node->start >= guide->start)
          : (end_node->end <= guide->end);

        if (in_range) {
          if (to_source) {
            right_coverage = end_node->cov / end_node->len();
            left_coverage = graph.nodes[neighbor]->cov
                / graph.nodes[neighbor]->len();
          } else {
            left_coverage = end_node->cov / end_node->len();
            right_coverage = graph.nodes[neighbor]->cov
                / graph.nodes[neighbor]->len();
          }
        }
      }
    }
    // Source/sink is already a parent/child of the guide's terminal node
    if (already_linked) return;

    // If the guide starts/ends inside the node -> compute the coverage on
    // either side of the boundary from the bp coverage
    if (to_source) {
      if (end_node->start < guide->start) {
        // Cumulative bp coverage before the guide start
        left_coverage += get_region_coverage('.', end_node->start - ref_start,
          guide->start - 1 - ref_start, bundle->bp_coverage);
        left_coverage /= (guide->start - end_node->start);
        uint32_t overlap_end = end_node->end + 1;
        if (guide->end < end_node->end) overlap_end = guide->end + 1;

        // Cumulative bp coverage from the guide start to the end of the overlap
        if (overlap_end > guide->start) {
          right_coverage += get_region_coverage('.', guide->start - ref_start,
            overlap_end - 1 - ref_start, bundle->bp_coverage);
          right_coverage /= (overlap_end - guide->start);
        }
      }
    } else { // node is longer than guide
      if (end_node->end > guide->end) {
        uint32_t overlap_start = end_node->start;
        if (guide->start > end_node->start) overlap_start = guide->start;

        // Cumulative bp coverage from the overlap start to the guide end
        left_coverage += get_region_coverage('.', overlap_start - ref_start,
          guide->end - ref_start, bundle->bp_coverage);
        left_coverage /= (guide->end - overlap_start + 1);

        // Cumulative bp coverage after the guide end
        right_coverage += get_region_coverage('.', guide->end - ref_start,
          end_node->end - ref_start, bundle->bp_coverage);
        right_coverage /= (end_node->end - guide->end + 1);
      }
    }

    // Abundance of the new link: the coverage step at the boundary, but at
    // least the minimum transfrag abundance
    double max_abundance = min_transfrag_abundance;
    if (to_source ? (right_coverage > left_coverage)
        : (right_coverage < left_coverage)) {
      max_abundance = to_source ? (right_coverage - left_coverage)
          : (left_coverage - right_coverage);
    }

    // Add a source->end node (or end node->sink) transfrag
    std::vector<int32_t> path_nodes;
    if (to_source) {
      path_nodes.push_back(0);
      path_nodes.push_back(end_node_id);
    } else {
      path_nodes.push_back(end_node_id);
      path_nodes.push_back(terminal_id);
    }
    BitVec edge_path_bits(graph.node_count + graph.edge_count);
    edge_path_bits[to_source ? 0 : end_node_id] = 1;
    edge_path_bits[to_source ? end_node_id : terminal_id] = 1;
    edge_path_bits[edge_bit] = 1;
    graph.transfrag.push_back(std::make_unique<CTransfrag>(
        cell_count, path_nodes, edge_path_bits, max_abundance));
    graph.transfrag.back()->made_by_guide = true;

    // Record the source/sink link in the end node's parents/children and edge
    // path_bits
    if (to_source) {
      end_node->parent.insert(end_node->parent.begin(), 0);
      // Source should already be among the node's parents, but not the edge to
      // source
      end_node->parent_reachable_bits[edge_bit] = 1;
    } else {
      end_node->child.push_back(terminal_id);
      // The sink is not yet among the node's children; record the edge to sink
      end_node->child_reachable_bits[edge_bit] = 1;
    }
  }

  // Matches every usable guide of this strand to the graph, marks guides that
  // are contained in longer ones and extends each guide to the source and
  // sink. Returns the matched guides in matched_guides, longest first.
  void process_reference_guides(SpliceGraph& graph, char strand,
                                std::vector<CGuide>& matched_guides,
                                BundleData *bundle) {
    // Get the guides of this bundle.
    std::vector<GffObj *>& guides = bundle->guides;
    int32_t cell_count = count_of(bundle->cell_name);

    // Match each guide of this strand that is covered by reads to the graph.
    for (int32_t g = 0; g < count_of(guides); g++) {
      bool strand_matches
          = guides[g]->strand == strand || guides[g]->strand == '.';
      if (strand_matches && is_guide_supported_by_reads(guides[g])
          && guides[g]->overlap(graph.nodes[1]->start,
                                graph.nodes[graph.node_count - 2]->end)) {
        std::shared_ptr<CTransfrag> guide_transfrag = match_guide_to_graph(
            guides[g], graph, cell_count);
        if (guide_transfrag) { // the guide can be found among the graph nodes

          CGuide new_guide(guide_transfrag, g);
          guide_transfrag->real = true;
          guide_transfrag->guide = 1;
          matched_guides.push_back(new_guide);
        }
      }
    } // for g

    // Order the matched guides by path_bits size, longest first.
    sort_by(matched_guides, compare_guide_by_path_bits_size);

    // Go through the guides from the longest path_bits to the shortest.
    uint32_t ref_start = (uint32_t)bundle->start;
    int32_t g = 0;
    while (g < count_of(matched_guides)) {
      // Check if guide is included in previous guide paths
      int32_t longer_idx = 0;

      // Here we assume that the guides are sorted by path_bits
      while (longer_idx < g) {
        if (matched_guides[longer_idx].transfrag->path_bits.contains(
                matched_guides[g].transfrag->path_bits)) {
          // A guide contained in a longer one loses its real/guide status.
          matched_guides[g].transfrag->real = false;
          matched_guides[g].transfrag->guide = 0;
          break;
        }
        longer_idx++;
      }

      // Find if first node of matched_guides extends to source
      extend_guide_to_terminal(true, g, graph, matched_guides, guides, bundle,
                               cell_count, ref_start);

      // Find if we can extend matched_guides to sink
      extend_guide_to_terminal(false, g, graph, matched_guides, guides, bundle,
                               cell_count, ref_start);

      g++;
    }
  }

  // Returns whether [start,end] lies inside one exon of a covered guide on
  // the strand given by strand ('.' for either).
  bool is_within_guide_exon(std::vector<GffObj *>& guides, char strand,
                            uint32_t start, uint32_t end) {
    // Look for a covered guide on that strand with an exon that contains it.
    for (int32_t g = 0; g < count_of(guides); g++)
      if (is_guide_supported_by_reads(guides[g])) {
        // Guide overlaps, so look at its exons
        if ((strand == '.' || guides[g]->strand == strand)
            && guides[g]->overlap(start, end)) {
          for (int32_t i = 0; i < guides[g]->exons.Count(); i++) {
            // There won't be any further overlap
            if (end < guides[g]->exons[i]->start) break;
            // Overlap is before start of exon
            if (start < guides[g]->exons[i]->start) break;
            if (end <= guides[g]->exons[i]->end) {
              // Start is not before the exon's start and end is not past its
              // end: inside this exon
              return true;
            }
          }
        }
      }

    return false;
  }

} // namespace ribbon
