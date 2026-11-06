#include "ribbon.h"
#include "GThreads.h"

// Flow for reference (guide) transcripts: guides claim their evidence first,
// their ends can move to better supported nodes, and their cells are credited.

namespace ribbon {

  // Credits this cell with added_coverage of the guide prediction. A guide
  // has one prediction for all cells, so its cell is -1 and its cell_cov holds
  // the coverage of the cells credited with it (print_results merges it).
  // credit_guides_found_by_denovo later clears the cell's credit where de novo
  // does not find the transcript (the flow credit stays in flow_cov).
  void credit_guide_cell(const FlowContext& context, CPrediction *prediction,
                         double added_coverage) {
    prediction->cell_cov.resize(count_of(context.bundle->cell_name), 0.0);
    prediction->cell_cov[context.cell] += added_coverage;
  }

  // Makes the new prediction (just stored for a guide) the guide's
  // all-cells prediction, crediting the cell that stored it.
  void claim_guide_prediction(const FlowContext& context, CPrediction *prediction) {
    prediction->cell = -1;
    credit_guide_cell(context, prediction, prediction->cov);
  }

  // Adds the flow of path to the guide prediction predictions[p]: the
  // used fraction (node_flux) of each node's coverage goes to the prediction exons
  // that overlap it, and removes it from node_coverage. The cell is credited
  // with what it added.
  void update_guide_prediction(const FlowContext& context,
                         PredictionList& predictions,
                         int32_t p, std::vector<int32_t>& path,
                         std::vector<float>& node_flux,
                         std::vector<float>& node_coverage) {
    // Graph and cell data shared through the context.

    int32_t exon_idx = 0;
    int32_t exon_count = count_of(predictions[p]->exons);
    double cell_added_coverage = 0;

    for (int32_t i = 0; i < count_of(path); i++) {
      // Not first or last node in graph
      if (path[i] && path[i] < context.graph.node_count - 1) {
        CGraphnode *node = context.graph.nodes[path[i]].get();
        double added_coverage = node_flux[i] * node_coverage[path[i]];
        node_coverage[path[i]] -= added_coverage;
        // Addcov/=node->len(); this is not necessary because it was already
        // divided to node length
        // if exon end before start of node skip it
        while (exon_idx < exon_count
               && predictions[p]->exons[exon_idx].end < node->start) exon_idx++;
        // There is overlap between exon and node
        while (exon_idx < exon_count
               && predictions[p]->exons[exon_idx].start <= node->end) {
          int32_t overlap_len
              = node->overlapLen(predictions[p]->exons[exon_idx].start,
              predictions[p]->exons[exon_idx].end);
          double exon_coverage = added_coverage * overlap_len;
          predictions[p]->exon_cov[exon_idx] += exon_coverage
              / predictions[p]->exons[exon_idx].len();
          predictions[p]->cov += exon_coverage / abs(predictions[p]->tlen);
          cell_added_coverage += exon_coverage / abs(predictions[p]->tlen);
          exon_idx++;
        }
      }
    }

    credit_guide_cell(context, predictions[p].get(), cell_added_coverage);
  }

  // Prepares what the flow functions reuse of the guide for every cell: which
  // transfrags of the graph have their path_bits inside the path_bits of the
  // guide, and, for each node of its path, the transfrags with abundance that
  // are inside (the others only count with the abundance of the cell). This
  // depends on the guide and the transfrags only, not on the cell. Call it
  // again when the path of the guide changes.
  void prepare_guide_flow_cache(const FlowContext& context,
                                CTransfrag& guide) {
    auto cache = std::make_shared<GuideFlowCache>();
    int32_t path_len = count_of(guide.nodes);
    cache->contained_transfrags = BitVec(count_of(context.graph.transfrag));
    for (int32_t t = 0; t < count_of(context.graph.transfrag); t++)
      if (guide.path_bits.contains(context.graph.transfrag[t]->path_bits))
        cache->contained_transfrags[t] = 1;

    cache->path_transfrags.resize(path_len);
    cache->pooled_in.resize(path_len);
    cache->pooled_out.resize(path_len);
    for (int32_t i = 1; i < path_len - 1; i++) {
      int32_t path_node = guide.nodes[i];
      for (int32_t t : context.graph.nodes[path_node]->transfrag_ids) {
        const CTransfrag *transfrag = context.graph.transfrag[t].get();
        if (!transfrag->abundance) continue;
        if (cache->contained_transfrags[t]) {
          cache->path_transfrags[i].push_back(t);
        } else {
          if (path_node > transfrag->nodes[0])
            cache->pooled_in[i] += transfrag->abundance;
          if (path_node < transfrag->nodes.back())
            cache->pooled_out[i] += transfrag->abundance;
        }
      }
    }
    guide.flow_cache = cache;
  }

  // Like push_flow_and_deplete_transfrags, but along the nodes of guide g; a
  // transfrag shared by several overlapping guides is split between them in
  // proportion to their abundance, so that no guide claims another's evidence.
  float push_guide_flow_and_deplete_transfrags(const FlowContext& context,
      int32_t g, std::vector<CGuide>& matched_guides, BitVec& transfrag_on_path,
      std::vector<float>& node_flux) {
    // Graph and cell data shared through the context.

    // Positions of the guide's nodes and the entering/leaving abundances.
    int32_t path_len = count_of(matched_guides[g].transfrag->nodes);
    std::vector<int32_t> path_pos_of_node;
    path_pos_of_node.resize(context.graph.node_count, -1);

    for (int32_t i = 0; i < path_len; i++) {
      path_pos_of_node[matched_guides[g].transfrag->nodes[i]] = i;
      node_flux.push_back(0.0);
    }

    // How many transcripts compatible to path enter node
    std::vector<float> path_in_abundance;
    // How many transcripts compatible to path exit node
    std::vector<float> path_out_abundance;
    path_in_abundance.resize(path_len);
    path_out_abundance.resize(path_len);
    std::vector<float> total_in_abundance; // how many transcripts enter node
    std::vector<float> total_out_abundance; // how many transcripts exit node
    total_in_abundance.resize(path_len);
    total_out_abundance.resize(path_len);

    const GuideFlowCache& cache = *matched_guides[g].transfrag->flow_cache;

    // Compute capacities and sums for all nodes
    for (int32_t i = 1; i < path_len - 1; i++) {
      int32_t path_node = matched_guides[g].transfrag->nodes[i];

      // The transfrags inside the guide count with their pooled abundance,
      // shared with the earlier guides that contain them in proportion to
      // their abundance.
      for (int32_t t : cache.path_transfrags[i]) {
        transfrag_on_path[t] = 1;
        double total_abundance = matched_guides[g].transfrag->abundance;
        for (int32_t other_idx = g - 1; other_idx >= 0; other_idx--)
          if (matched_guides[other_idx].transfrag->flow_cache
                  ->contained_transfrags[t]) {
            total_abundance += matched_guides[other_idx].transfrag->abundance;
          }
        float share = 1;
        if (total_abundance > matched_guides[g].transfrag->abundance)
          share = matched_guides[g].transfrag->abundance / total_abundance;

        context.graph.transfrag[t]->use_path
            = share * context.graph.transfrag[t]->abundance;

        // Transfrag starts before this node
        if (context.graph.transfrag[t]->nodes[0] < path_node) {
          total_in_abundance[i] += context.graph.transfrag[t]->use_path;
          path_in_abundance[i] += context.graph.transfrag[t]->use_path;
        }
        // Transfrag ends after this node
        if (context.graph.transfrag[t]->nodes.back() > path_node) {
          total_out_abundance[i] += context.graph.transfrag[t]->use_path;
          path_out_abundance[i] += context.graph.transfrag[t]->use_path;
        }
      }

      // The other transfrags count with the abundance of this cell
      for (int32_t t : transfrags_of_cell(context, path_node)) {
        if (!context.graph.transfrag[t]->abundance
            || cache.contained_transfrags[t])
          continue;
        if (path_node > context.graph.transfrag[t]->nodes[0])
          total_in_abundance[i]
              += context.graph.transfrag[t]->cell_abundance[context.cell];
        if (path_node < context.graph.transfrag[t]->nodes.back())
          total_out_abundance[i]
              += context.graph.transfrag[t]->cell_abundance[context.cell];
      }

      apply_sparse_coverage_fallback(i, path_in_abundance, total_in_abundance,
                                     path_out_abundance, total_out_abundance);
    }

    // Compute flow -> deal with rounding errors here
    float running_flow = compute_path_flow(path_len, path_in_abundance,
                                           total_in_abundance,
                                           path_out_abundance,
                                           total_out_abundance);
    if (!running_flow) return (0);

    // Compute the fraction of each node's coverage that the flow uses.
    for (int32_t i = path_len - 2; i > 0; i--) {
      node_flux[i] = running_flow / total_out_abundance[i];
      // Because of rounding errors this could become more than 1
      if (node_flux[i] > 1) node_flux[i] = 1;
      path_out_abundance[i] = running_flow;
      running_flow = running_flow * total_in_abundance[i]
          / total_out_abundance[i];
    }

    // Subtract the flow from the transfrags that start at each path node.
    for (int32_t i = 1; i < path_len - 1; i++) {
      int32_t path_node = matched_guides[g].transfrag->nodes[i];
      const std::vector<int32_t>& node_transfrags
          = transfrags_of_cell(context, path_node);
      int32_t transfrag_count = count_of(node_transfrags);
      for (int32_t j = 0; j < transfrag_count; j++) {
        int32_t t = node_transfrags[j];
        if (transfrag_on_path[t]
            && context.graph.transfrag[t]->cell_abundance[context.cell]) {
          // use_path is set for every transfrag on the path in the capacity
          // loop above, and transfrag_on_path is reset after each guide, so it
          // is the abundance this guide takes from the transfrag. It shrinks
          // below as the transfrag gets depleted.
          double cell_abundance = context.graph.transfrag[t]->use_path;

          // Transfrag starts at this node
          if (context.graph.transfrag[t]->nodes[0] == path_node) {
            if (path_out_abundance[i] > cell_abundance) {
              path_out_abundance[i] -= cell_abundance;
              int32_t end_pos
                  = path_pos_of_node[context.graph.transfrag[t]->nodes.back()];
              for (int32_t k = i + 1; k < end_pos; k++) {
                path_out_abundance[k] -= cell_abundance;
              }
              if (!deplete_transfrag_abundance(context.graph.transfrag[t].get(),
                                               context.cell,
                                               cell_abundance, true))
                // We only need to delete this if we still have abundance left
                // in transfrag
                context.graph.transfrag[t]->use_path = 0;
            } else if (path_out_abundance[i]) {
              if (!deplete_transfrag_abundance(context.graph.transfrag[t].get(),
                                               context.cell,
                                               path_out_abundance[i], true)) {
                // The transfrag still has abundance, so this guide has taken
                // path_out_abundance[i] of its use_path. If the transfrag was
                // floored to 0 instead, its cell abundance is 0 and the check
                // above skips it from now on.
                context.graph.transfrag[t]->use_path -= path_out_abundance[i];
                if (context.graph.transfrag[t]->use_path < epsilon)
                  context.graph.transfrag[t]->use_path = 0;
              }
              int32_t end_pos
                  = path_pos_of_node[context.graph.transfrag[t]->nodes.back()];
              for (int32_t k = i + 1; k < end_pos; k++) {
                path_out_abundance[k] -= path_out_abundance[i];
              }
              path_out_abundance[i] = 0;
              // Break; we can not break because we might still have transfrag
              // that need the cleanup
            }
          }
        }
      }
    }

    // Subtract the flow from the first transfrag that starts at the source.
    const std::vector<int32_t>& source_transfrags
        = transfrags_of_cell(context, matched_guides[g].transfrag->nodes[0]);
    int32_t transfrag_count = count_of(source_transfrags);
    for (int32_t j = 0; j < transfrag_count; j++) {
      int32_t t = source_transfrags[j];
      if (transfrag_on_path[t]
          && context.graph.transfrag[t]->cell_abundance[context.cell]) {
        deplete_transfrag_abundance(context.graph.transfrag[t].get(),
            context.cell,
                                    running_flow,
                                    false);
        break;
      }
    } //*/

    return (node_flux[1]);
  }

  // Moves the last node of a guide's path (the one before the sink) to
  // new_last_node: drops the nodes after it, or adds the nodes up to it, and
  // fixes the node and edge bits of the path_bits.
  void move_guide_path_end(CTransfrag& transfrag, int32_t new_last_node,
                           PathBitLayout& layout, int32_t node_count) {
    int32_t sink = node_count - 1;
    auto set_bit = [&](int32_t from, int32_t to, bool value) {
      layout.set_edge_bit(transfrag.path_bits, from, to, value);
    };
    int32_t last_node = transfrag.nodes[count_of(transfrag.nodes) - 2];

    if (new_last_node < last_node) {
      // The end moves left: drop the nodes after the new last node.
      auto first_dropped = std::upper_bound(transfrag.nodes.begin() + 1,
                                            transfrag.nodes.end() - 1,
                                            new_last_node);
      set_bit(last_node, sink, false);
      for (auto node = first_dropped; node != transfrag.nodes.end() - 1;
           ++node) {
        transfrag.path_bits[*node] = 0;
        set_bit(*(node - 1), *node, false);
      }
      transfrag.nodes.erase(first_dropped, transfrag.nodes.end() - 1);
      set_bit(new_last_node, sink, true);
    } else if (new_last_node > last_node) {
      // The end moves right: add the nodes up to the new last node.
      set_bit(last_node, sink, false);
      transfrag.nodes.pop_back();
      for (int32_t node = last_node + 1; node <= new_last_node; node++) {
        transfrag.path_bits[node] = 1;
        set_bit(node - 1, node, true);
        transfrag.nodes.push_back(node);
      }
      transfrag.nodes.push_back(sink);
      set_bit(new_last_node, sink, true);
    }
  }

  // Moves the first node of a guide's path (the one after the source) to
  // new_first_node: drops the nodes before it, or adds the nodes down to it,
  // and fixes the node and edge bits of the path_bits.
  void move_guide_path_start(CTransfrag& transfrag, int32_t new_first_node,
                             PathBitLayout& layout) {
    auto set_bit = [&](int32_t from, int32_t to, bool value) {
      layout.set_edge_bit(transfrag.path_bits, from, to, value);
    };
    int32_t first_node = transfrag.nodes[1];

    if (new_first_node > first_node) {
      // The start moves right: drop the nodes before the new first node.
      set_bit(0, new_first_node, true);
      set_bit(0, first_node, false);
      auto first_kept = std::lower_bound(transfrag.nodes.begin() + 1,
                                         transfrag.nodes.end(),
                                         new_first_node);
      for (auto node = transfrag.nodes.begin() + 1; node != first_kept;
           ++node) {
        transfrag.path_bits[*node] = 0;
        set_bit(*node, *(node + 1), false);
      }
      transfrag.nodes.erase(transfrag.nodes.begin() + 1, first_kept);
    } else if (new_first_node < first_node) {
      // The start moves left: add the nodes back down to the new first node.
      set_bit(0, new_first_node, true);
      set_bit(0, first_node, false);
      std::vector<int32_t> added_nodes;
      for (int32_t node = new_first_node; node < first_node; node++) {
        transfrag.path_bits[node] = 1;
        set_bit(node, node + 1, true);
        added_nodes.push_back(node);
      }
      transfrag.nodes.insert(transfrag.nodes.begin() + 1, added_nodes.begin(),
                       added_nodes.end());
    }
  }

  // Tries to move the start and/or end of the multi-exon guide g to a better
  // supported first or last exon node, from the transfrags starting at the
  // source/sink. If it moves, pushes flow along the new path and returns it.
  float try_alternative_guide_ends(const FlowContext& context, int32_t g,
                                   std::vector<CGuide>& matched_guides,
                                   BitVec& transfrag_on_path,
                                   std::vector<float>& node_flux) {
    // Graph and cell data shared through the context.
    CTransfrag& transfrag = *matched_guides[g].transfrag;

    float flux = 0;

    // Find the first exon of the guide and, if it has several nodes, the start
    // node that the source transfrags of this cell support best.
    int32_t new_first_node = transfrag.nodes[1];
    int32_t scan_idx = 2;
    while (scan_idx < count_of(transfrag.nodes)) {
      // We know there is at least one node like this since this is multiexon
      if (context.graph.nodes[transfrag.nodes[scan_idx - 1]]->end + 1
          < context.graph.nodes[transfrag.nodes[scan_idx]]->start) {
        break;
      }
      scan_idx++;
    }
    // This shouldn't happen
    if (scan_idx == count_of(transfrag.nodes)) return (flux);
    if (scan_idx > 2) { // if there is another node in the first exon
      int32_t run_node = new_first_node;
      // Find continuous nodes
      while (run_node
             && context.graph.nodes[run_node - 1]->end + 1
             == context.graph.nodes[run_node]->start) run_node--;
      if (!run_node) run_node++;
      double best_abundance = 0;
      int32_t transfrag_count = count_of(context.graph.nodes[0]->transfrag_ids);
      for (int32_t j = 0; j < transfrag_count; j++) {
        int32_t t = context.graph.nodes[0]->transfrag_ids[j];
        if (context.graph.transfrag[t]->nodes.back() < transfrag.nodes[scan_idx]
            && context.graph.transfrag[t]->nodes.back() >= run_node
            && context.graph.transfrag[t]->cell_abundance[context.cell]
                > best_abundance) {
          new_first_node = context.graph.transfrag[t]->nodes.back();
          best_abundance
              = context.graph.transfrag[t]->cell_abundance[context.cell];
        }
      }
    }

    // Do the same for the last exon: the end node best supported by the
    // transfrags that go to the sink.
    int32_t new_last_node = transfrag.nodes[count_of(transfrag.nodes) - 2];
    scan_idx = count_of(transfrag.nodes) - 3;
    while (scan_idx > 0) {
      if (context.graph.nodes[transfrag.nodes[scan_idx]]->end + 1
          < context.graph.nodes[transfrag.nodes[scan_idx + 1]]->start) {
        break;
      }
      scan_idx--;
    }
    if (!scan_idx) return (flux); // this shouldn't happen
    if (scan_idx < count_of(transfrag.nodes) - 3) {
      int32_t run_node = new_last_node;
      // Find continuous nodes
      while (run_node < context.graph.node_count - 1
             && context.graph.nodes[run_node]->end + 1
                 == context.graph.nodes[run_node + 1]->start) run_node++;
      double best_abundance = 0;
      const std::vector<int32_t>& sink_transfrag_ids
          = context.graph.nodes[context.graph.node_count - 1]->transfrag_ids;
      int32_t transfrag_count = count_of(sink_transfrag_ids);
      for (int32_t j = 0; j < transfrag_count; j++) {
        int32_t t = sink_transfrag_ids[j];
        if (context.graph.transfrag[t]->nodes[0] > transfrag.nodes[scan_idx]
            && context.graph.transfrag[t]->nodes[0] <= run_node
            && context.graph.transfrag[t]->cell_abundance[context.cell]
                > best_abundance) {
          new_last_node = context.graph.transfrag[t]->nodes[0];
          best_abundance
              = context.graph.transfrag[t]->cell_abundance[context.cell];
        }
      }
    }

    // If the start or the end changed, rewrite the nodes and the path_bits of
    // the guide and push flow along the new path.
    if (new_first_node != transfrag.nodes[1]
        || new_last_node != transfrag.nodes[count_of(transfrag.nodes) - 2]) {
      move_guide_path_end(transfrag, new_last_node,
          context.graph.path_bit_layout,
                          context.graph.node_count);
      move_guide_path_start(transfrag, new_first_node,
          context.graph.path_bit_layout);
      prepare_guide_flow_cache(context, transfrag);

      // Push flow along the modified path.
      flux = push_guide_flow_and_deplete_transfrags(context, g, matched_guides,
          transfrag_on_path, node_flux);
      transfrag_on_path.reset();
    }

    return (flux);
  }

  // Estimates how much abundance the path of a guide can claim, without
  // depleting any transfrag, with the same in/out propagation as the other flow
  // functions. Returns the flow weighted by the coverage of each path node.
  float estimate_guide_flow(const FlowContext& context,
                            std::vector<int32_t>& path,
                            const GuideFlowCache& cache) {
    // Graph and cell data shared through the context.

    double guide_abundance = 0;

    // Positions of the path nodes and the entering/leaving abundances.
    int32_t path_len = count_of(path);
    std::vector<int32_t> path_pos_of_node;
    path_pos_of_node.resize(context.graph.node_count, -1);
    for (int32_t i = 0; i < path_len; i++) {
      path_pos_of_node[path[i]] = i;
    }
    // How many transcripts compatible to path enter node
    std::vector<float> path_in_abundance;
    // How many transcripts compatible to path exit node
    std::vector<float> path_out_abundance;
    path_in_abundance.resize(path_len);
    path_out_abundance.resize(path_len);
    std::vector<float> total_in_abundance; // how many transcripts enter node
    std::vector<float> total_out_abundance; // how many transcripts exit node
    total_in_abundance.resize(path_len);
    total_out_abundance.resize(path_len);

    // Compute capacities and sums for all nodes
    for (int32_t i = 1; i < path_len - 1; i++) {
      // The transfrags inside the guide count with the abundance of this cell
      for (int32_t t : transfrags_of_cell(context, path[i])) {
        if (!context.graph.transfrag[t]->abundance
            || !cache.contained_transfrags[t])
          continue;
        // Transfrag starts before this node
        if (context.graph.transfrag[t]->nodes[0] < path[i]) {
          total_in_abundance[i]
              += context.graph.transfrag[t]->cell_abundance[context.cell];
          path_in_abundance[i]
              += context.graph.transfrag[t]->cell_abundance[context.cell];
        }
        // Transfrag ends after this node
        if (context.graph.transfrag[t]->nodes.back() > path[i]) {
          total_out_abundance[i]
              += context.graph.transfrag[t]->cell_abundance[context.cell];
          path_out_abundance[i]
              += context.graph.transfrag[t]->cell_abundance[context.cell];
        }
      }
      // The others count with their pooled abundance
      total_in_abundance[i] += cache.pooled_in[i];
      total_out_abundance[i] += cache.pooled_out[i];

      apply_sparse_coverage_fallback(i, path_in_abundance, total_in_abundance,
                                     path_out_abundance, total_out_abundance);
    }

    // Compute flow
    float running_flow = compute_path_flow(path_len, path_in_abundance,
                                           total_in_abundance,
                                           path_out_abundance,
                                           total_out_abundance);
    if (!running_flow) return (0);

    // Weight the flow by the coverage of each node.
    for (int32_t i = path_len - 2; i > 0; i--) {
      guide_abundance += context.graph.nodes[path[i]]->cov * running_flow
          / total_out_abundance[i];
      running_flow = running_flow * total_in_abundance[i]
          / total_out_abundance[i];
    }

    return (guide_abundance);
  }

  // Comparator for CGuide: guides contained in a longer guide (transfrag->real
  // is false) come last; the rest by decreasing abundance, then longest
  // path_bits.
  int32_t compare_guide_by_abundance(const void *p1, const void *p2) {
    CGuide *a = (CGuide *)p1;
    CGuide *b = (CGuide *)p2;
    // This ensures included guides are treated last
    if (a->transfrag->real && !b->transfrag->real) return 1;
    if (!a->transfrag->real && b->transfrag->real) return -1;
    // Most abundant guide takes precedence
    if (a->transfrag->abundance < b->transfrag->abundance) return 1;
    if (a->transfrag->abundance > b->transfrag->abundance) return -1;
    // Longest guide takes precedence
    if (a->transfrag->path_bits.count()
        < b->transfrag->path_bits.count()) return 1;
    if (a->transfrag->path_bits.count()
        > b->transfrag->path_bits.count()) return -1;
    return 0;
  }

  // Stores the prediction of matched_guides[g] if its guide has none yet (and
  // notes the guide as printed if it has several exons), otherwise updates the
  // prediction of the guide with the flow that was just pushed.
  void store_or_update_guide_prediction(
      const FlowContext& context, int32_t g,
      std::vector<CGuide>& matched_guides, std::vector<float>& node_flux,
      std::vector<float>& node_coverage, int32_t& gene_count,
      PredictionList& predictions, BitVec& current_path_bits,
      bool& is_first_transcript, std::vector<GffObj *>& guides, std::vector<int32_t>& guide_prediction,
      std::vector<int32_t>& printed_guides) {
    CGuide& matched_guide = matched_guides[g];
    GffObj *guide = guides[matched_guide.g];
    if (guide_prediction[matched_guide.g] == -1) {
      bool include = true;
      // The prev_prediction_count is set to 0 here to address later
      store_transcript(context, 0, predictions, matched_guide.transfrag->nodes,
                       node_flux, node_coverage, gene_count,
                       is_first_transcript, include,
                       current_path_bits, context.bundle, guide);
      // NEED TO TEST: if this doesn't work for single genes we might want to
      // recombine with the previous prediction in store_transcript
      guide_prediction[matched_guide.g] = count_of(predictions) - 1;
      claim_guide_prediction(context, predictions.back().get());
      // The guide was printed
      if (guide->exons.Count() > 1) printed_guides.push_back(g);
    } else {
      update_guide_prediction(context, predictions,
                              guide_prediction[matched_guide.g],
                              matched_guide.transfrag->nodes, node_flux,
                              node_coverage);
    }
  }

  // Gives the guides that matched this graph first claim on the evidence: one
  // guide gets its flow directly; with several, we rank them by estimated flow
  // and push the weakest first. Multi-exon guides then try better ends.
  void assign_flow_to_guides(const FlowContext& context,
                             std::vector<CGuide>& matched_guides,
                             int32_t& gene_count, PredictionList& predictions,
                             std::vector<float>& node_coverage,
                             BitVec& transfrag_on_path,
                             BitVec& current_path_bits,
                             bool& is_first_transcript,
                             std::vector<GffObj *>& guides,
                             std::vector<int32_t>& guide_prediction,
                             std::vector<int32_t>& printed_guides) {
    // Graph and cell data shared through the context.

    int32_t guide_count = count_of(matched_guides);

    // The first cell of the graph finds out which transfrags each guide holds
    for (CGuide& guide : matched_guides)
      if (!guide.transfrag->flow_cache)
        prepare_guide_flow_cache(context, *guide.transfrag);

    // Stores or updates the prediction of matched_guides[g] after pushing flow.
    auto store_or_update = [&](int32_t g, std::vector<float>& node_flux) {
      store_or_update_guide_prediction(context, g, matched_guides, node_flux,
          node_coverage, gene_count, predictions, current_path_bits,
          is_first_transcript, guides, guide_prediction, printed_guides);
    };

    // Gives the multi-exon guide g a second chance with better supported ends.
    auto try_better_ends = [&](int32_t g) {
      std::vector<float> node_flux;
      float flux = try_alternative_guide_ends(context, g, matched_guides,
                                              transfrag_on_path, node_flux);
      if (flux > epsilon) store_or_update(g, node_flux);
    };

    if (guide_count == 1) { // if only one guide we do not need to do the 2 pass

      // Push the flow of the single guide and store (or update) its prediction.
      std::vector<float> node_flux;

      float flux = push_flow_and_deplete_transfrags(context,
          matched_guides[0].transfrag->nodes, transfrag_on_path, node_flux,
          matched_guides[0].transfrag->path_bits);
      transfrag_on_path.reset();

      if (flux > epsilon) store_or_update(0, node_flux);

      // Then give the guide a second chance with better supported ends.
      if (guides[matched_guides[0].g]->exons.Count() > 1) try_better_ends(0);

    } else if (guide_count) {
      // Estimate how much flow each guide can carry and rank the guides by it.
      for (int32_t g = 0; g < count_of(matched_guides); g++) {
        matched_guides[g].transfrag->abundance = estimate_guide_flow(context,
            matched_guides[g].transfrag->nodes,
            *matched_guides[g].transfrag->flow_cache);
        transfrag_on_path.reset();
      }
      sort_by(matched_guides, compare_guide_by_abundance);

      // Push the flow guide by guide, from the weakest to the strongest, and
      // store or update the prediction of each.
      for (int32_t g = guide_count - 1; g >= 0; g--)
        // Calculate maximum push flow for each guide starting from the less
        // covered one
        if (matched_guides[g].transfrag->abundance) {
          std::vector<float> node_flux;
          float flux = push_guide_flow_and_deplete_transfrags(context, g,
              matched_guides, transfrag_on_path, node_flux);

          transfrag_on_path.reset();

          if (flux > epsilon) {
            store_or_update(g, node_flux);

          // It's possible that this is a single exon gene included in a much
          // larger interval and this is why it didn't get predicted
          } else if (guides[matched_guides[g].g]->exons.Count() == 1) {
            // Single exon gene included in another prediction
            // Check if it overlaps other single exon genes

            bool overlap = false;
            for (int32_t r = guide_count - 1; r > g; r--)
              if (guide_prediction[matched_guides[r].g] > -1) {
                if (matched_guides[r].transfrag->path_bits.contains(
                        matched_guides[g].transfrag->path_bits)) {
                  if (guides[matched_guides[r].g]->exons.Count() > 1) {
                    overlap = true;
                    break;
                  // Overlaps single exon gene -> check if it's a true overlap
                  } else {
                    if (guides[matched_guides[r].g]->exons[0]->overlap(
                        guides[matched_guides[g].g]->exons[0])) {
                      overlap = true;
                      break;
                    }
                  }
                }
              }
            if (!overlap) { // no overlap detected
              for (int32_t z = 0;
                   z < count_of(matched_guides[g].transfrag->nodes); z++)
                node_flux.push_back(1.0);
              store_or_update(g, node_flux);
            }
          }
        }

      // Give each multi-exon guide a second chance with better supported ends.
      for (int32_t g = 0; g < guide_count; g++)
        if (matched_guides[g].transfrag->abundance
            && guides[matched_guides[g].g]->exons.Count() > 1) {
          try_better_ends(g);
        }
    }
  }

  // Runs assign_flow_to_guides, then returns the node with the highest
  // remaining coverage, which is where the greedy extraction continues.
  int32_t assign_flow_to_guides_and_pick_seed(
      const FlowContext& context, std::vector<CGuide>& matched_guides,
      int32_t& gene_count, PredictionList& predictions,
      std::vector<float>& node_coverage, BitVec& transfrag_on_path,
      BitVec& current_path_bits, bool& is_first_transcript,
      std::vector<GffObj *>& guides, std::vector<int32_t>& guide_prediction) {
    // Graph and cell data shared through the context.

    // Remembers the multi-exons matched_guides's that were printed here
    std::vector<int32_t> printed_guides;

    assign_flow_to_guides(context, matched_guides, gene_count, predictions,
                          node_coverage, transfrag_on_path, current_path_bits,
                          is_first_transcript, guides, guide_prediction,
                          printed_guides);

    // Pick the node with the most remaining coverage as the seed.
    int32_t seed_node = 1;
    // Node coverages:
    for (int32_t i = 1; i < context.graph.node_count - 1; i++)
      if (node_coverage[i] > node_coverage[seed_node]) seed_node = i;

    return (seed_node);
  }

} // namespace ribbon
