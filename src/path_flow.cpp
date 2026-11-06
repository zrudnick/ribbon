#include "ribbon.h"
#include "GThreads.h"

// minimum number for a transcript to be printed
extern int32_t min_transcript_len;

// Pushing flow along a path of the splice graph: how much a path carries,
// depleting the transfrags that support it, and storing the path as a
// transcript prediction.

namespace ribbon {

  // Node i has no compatible-transfrag capacity entering and/or leaving it.
  // Single-cell coverage is sparse, so force a minimal (1-unit) flow through
  // the node instead of giving up.
  void apply_sparse_coverage_fallback(int32_t i,
      std::vector<float>& path_in_abundance,
      std::vector<float>& total_in_abundance,
      std::vector<float>& path_out_abundance,
      std::vector<float>& total_out_abundance) {
    if (!path_in_abundance[i]) {
      path_in_abundance[i] = 1;
      if (total_in_abundance[i] < 1) total_in_abundance[i] = 1;
    }
    if (!path_out_abundance[i]) {
      path_out_abundance[i] = 1;
      if (total_out_abundance[i] < 1) total_out_abundance[i] = 1;
    }
  }

  // Propagates a single bottleneck flow along the path from node 1, given the
  // entering/leaving capacities and totals of each node, capping it at every
  // node by the smaller side. Returns 0 if the flow is blocked.
  float compute_path_flow(int32_t path_len,
                                 std::vector<float>& path_in_abundance,
                                 std::vector<float>& total_in_abundance,
                                 std::vector<float>& path_out_abundance,
                                 std::vector<float>& total_out_abundance) {
    float running_flow = path_in_abundance[1];
    for (int32_t i = 1; i < path_len - 1; i++) {
      float left_fraction = running_flow / total_in_abundance[i];
      float right_fraction = path_out_abundance[i] / total_out_abundance[i];
      if (right_fraction > left_fraction) { // more transfrags leave node
        right_fraction = left_fraction;
      }
      running_flow = right_fraction * total_out_abundance[i];
    }
    return running_flow;
  }

  // Subtracts amount from the cell abundance of transfrag t and floors it to 0
  // below epsilon, or when zero_if_pooled_depleted is set and the pooled
  // abundance is negligible too. Returns true when it was floored to 0.
  bool deplete_transfrag_abundance(CTransfrag *t, int32_t cell,
                                          double amount,
                                          bool zero_if_pooled_depleted) {
    t->cell_abundance[cell] -= amount;
    bool is_cell_depleted = t->cell_abundance[cell] < epsilon;
    bool is_pooled_depleted = zero_if_pooled_depleted
                              && t->abundance < epsilon;
    if (is_cell_depleted || is_pooled_depleted) {
      t->cell_abundance[cell] = 0;
      return true;
    }
    return false;
  }

  // The transfrags of node that can have abundance in this cell: the list of
  // the cell when there is one, else all the transfrags of the node. The loops
  // over it check the abundance of the cell themselves.
  const std::vector<int32_t>& transfrags_of_cell(
      const FlowContext& context, int32_t node) {
    if (context.cell_transfrags) return (*context.cell_transfrags)[node];
    return context.graph.nodes[node]->transfrag_ids;
  }

  // Pushes a bottleneck flow along path (node ids from source to sink) using
  // this cell's transfrag abundances and subtracts it from the transfrags that
  // support the path. Fills node_flux; returns node_flux[1] (0 if no flow).
  float push_flow_and_deplete_transfrags(const FlowContext& context,
                                         std::vector<int32_t>& path,
                                         BitVec& transfrag_on_path,
                                         std::vector<float>& node_flux,
                                         BitVec& current_path_bits) {
    // Graph and cell data shared through the context.

    // Positions of the path nodes and the entering/leaving abundances (a sparse
    // single cell may cover the path only partly).

    int32_t path_len = count_of(path);
    std::vector<int32_t> path_pos_of_node;
    path_pos_of_node.resize(context.graph.node_count, -1);
    for (int32_t i = 0; i < path_len; i++) {
      path_pos_of_node[path[i]] = i;
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

    // Compute capacities and sums for all nodes
    for (int32_t i = 1; i < path_len - 1; i++) {
      const std::vector<int32_t>& node_transfrags
          = transfrags_of_cell(context, path[i]);
      int32_t transfrag_count = count_of(node_transfrags);
      for (int32_t j = 0; j < transfrag_count; j++) {
        int32_t t = node_transfrags[j];
        // A transfrag without abundance in this cell adds nothing to the sums
        // below, and the loops that deplete need it to have some
        if (context.graph.transfrag[t]->abundance
            && context.graph.transfrag[t]->cell_abundance[context.cell]) {
          bool on_path = false;
          if (transfrag_on_path[t]) on_path = true;
          else if (!context.graph.transfrag[t]->nodes[0]) {
            if (context.graph.transfrag[t]->nodes.back() == path[1]) on_path
                = true;
          } else if (context.graph.transfrag[t]->nodes.back() == path.back()) {
            if (context.graph.transfrag[t]->nodes[0]
                == path[path_len - 2]) on_path
                = true;
          // Only need to check transfrag the first time we encounter it
          } else if (context.graph.transfrag[t]->nodes[0] == path[i]
                     && current_path_bits.contains(
                            context.graph.transfrag[t]->path_bits)) {
            on_path = true;
          }

          if (on_path) { // transcript on path
            transfrag_on_path[t] = 1;

            // Transfrag starts before this node
            if (context.graph.transfrag[t]->nodes[0] < path[i]) {
              // Consider the scell abundance instead of total abundance of
              // transfrag
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
          } else { // transfrag not on path
            if (path[i] > context.graph.transfrag[t]->nodes[0]) {
              total_in_abundance[i]
                  += context.graph.transfrag[t]->cell_abundance[context.cell];
            }

            if (path[i] < context.graph.transfrag[t]->nodes.back()) {
              total_out_abundance[i]
                  += context.graph.transfrag[t]->cell_abundance[context.cell];
            }
          }
        }
      }

      apply_sparse_coverage_fallback(i, path_in_abundance, total_in_abundance,
                                     path_out_abundance, total_out_abundance);
    }

    // Compute flow
    float running_flow = compute_path_flow(path_len, path_in_abundance,
                                           total_in_abundance,
                                           path_out_abundance,
                                           total_out_abundance);
    if (!running_flow) return (0);

    // Compute the fraction of each node's coverage that the flow uses.
    for (int32_t i = path_len - 2; i > 0; i--) {
      node_flux[i] = running_flow / total_out_abundance[i];
      path_out_abundance[i] = running_flow;
      running_flow = node_flux[i] * total_in_abundance[i];
    }

    // Subtract the flow from the transfrags that start at each path node.
    for (int32_t i = 1; i < path_len - 1; i++)
      if (path_out_abundance[i]) {
        const std::vector<int32_t>& node_transfrags
            = transfrags_of_cell(context, path[i]);
        int32_t transfrag_count = count_of(node_transfrags);
        for (int32_t j = 0; j < transfrag_count; j++) {
          int32_t t = node_transfrags[j];
          // Updates only scell abundances - if not empty
          if (transfrag_on_path[t]
              && context.graph.transfrag[t]->cell_abundance[context.cell]) {
            double cell_abundance
                = context.graph.transfrag[t]->cell_abundance[context.cell];

            // Transfrag starts at this node
            if (cell_abundance && context.graph.transfrag[t]->nodes[0]
                == path[i]) {
              if (path_out_abundance[i] > cell_abundance) {
                path_out_abundance[i] -= cell_abundance;
                int32_t end_pos
                    = path_pos_of_node[
                        context.graph.transfrag[t]->nodes.back()];
                for (int32_t k = i + 1; k < end_pos; k++) {
                  path_out_abundance[k] -= cell_abundance;
                }
                deplete_transfrag_abundance(context.graph.transfrag[t].get(),
                                            context.cell,
                                            cell_abundance, true);
              } else {
                deplete_transfrag_abundance(context.graph.transfrag[t].get(),
                                            context.cell,
                                            path_out_abundance[i], true);

                int32_t end_pos
                    = path_pos_of_node[
                        context.graph.transfrag[t]->nodes.back()];
                for (int32_t k = i + 1; k < end_pos; k++) {
                  path_out_abundance[k] -= path_out_abundance[i];
                }
                path_out_abundance[i] = 0;
                break;
              }
            }
          }
        }
      }

    // Subtract the flow from the first transfrag that starts at the source.
    const std::vector<int32_t>& source_transfrags
        = transfrags_of_cell(context, path[0]);
    int32_t transfrag_count = count_of(source_transfrags);
    for (int32_t j = 0; j < transfrag_count; j++) {
      int32_t t = source_transfrags[j];
      if (transfrag_on_path[t]
          && context.graph.transfrag[t]->cell_abundance[context.cell]) {
        deplete_transfrag_abundance(context.graph.transfrag[t].get(),
            context.cell,
                                    running_flow,
                                    false);
        // There is no point in updating more than one transfrag from source
        break;
      }
    }

    return (node_flux[1]);
  }

  // If the transcript (exons, exon_cov) equals prediction (same intron chain;
  // for one exon, an overlap of at least half), adds its coverage to prediction
  // and returns true; otherwise changes nothing and returns false.
  bool merge_if_equal_prediction(CPrediction *prediction, std::vector<GSeg>& exons,
                                 std::vector<float>& exon_cov) {
    // Different exon counts -> not equal
    if (count_of(prediction->exons) != count_of(exons)) return (false);

    // If prediction overlap at all
    if (exons[0].start < prediction->end && prediction->start < exons.back().end) {
      if (count_of(exons) == 1) { // if single exon
        int32_t len = (std::min(exons[0].end, prediction->end) - std::max(prediction->start,
            exons[0].start) + 1) / DROP;
        // If single exon make sure they overlap by a good amount (at least
        // half)
        if (len < (int32_t)prediction->len() && len
            < (int32_t)exons[0].len()) return (false);
        prediction->exon_cov[0] += exon_cov[0] * len / prediction->len();
        prediction->cov = prediction->exon_cov[0];
      } else { // multi exon
        // Make sure they have the same intron chain
        for (int32_t i = 0; i < count_of(exons) - 1; i++) {
          if (exons[i].end != prediction->exons[i].end
              || exons[i + 1].start != prediction->exons[i + 1].start) return false;
        }
        // For loop did not return false --> adjust coverages
        int32_t len = prediction->exons[0].len();
        if (exons[0].start > prediction->start) len -= exons[0].start - prediction->start;
        float cov = len * exon_cov[0];
        prediction->exon_cov[0] += cov / prediction->exons[0].len();
        int32_t last_exon_idx = count_of(exons) - 1;
        // Check internal exons, make sure they have the same intron chain
        for (int32_t i = 1; i < last_exon_idx; i++) {
          cov += exon_cov[i] * exons[i].len();
          prediction->exon_cov[i] += exon_cov[i];
        }
        // Adjust coverage of last exon
        len = prediction->exons[last_exon_idx].len();
        if (exons[last_exon_idx].end > prediction->end)
          len -= exons[last_exon_idx].end - prediction->end;
        cov += len * exon_cov[last_exon_idx];
        prediction->exon_cov[last_exon_idx] += cov / prediction->exons[last_exon_idx].len();
        prediction->cov += cov / prediction->tlen;
      }
      return (true);
    }
    return (false);
  }

  // Turns path (source to sink) into a prediction in predictions: exons are the
  // runs of adjacent nodes, with the flow's coverage, clipped to guide if given.
  // Merges copies from other cells. Returns the coverage stored, or 0.
  double store_transcript(const FlowContext& context, int32_t prev_prediction_count,
                          PredictionList& predictions, std::vector<int32_t>& path,
                          std::vector<float>& node_flux,
                          std::vector<float>& node_coverage,
                              int32_t& gene_count,
                          bool& is_first_transcript, bool& included,
                          BitVec& covered_path_bits, BundleData *bundle,
                          GffObj *guide) {
    // Graph and cell data shared through the context.

    // Collect the exons, their coverage and the transcript length.
    double cov = 0;
    // The same before PCR copies were folded (see weight_reads_by_molecule)
    double raw_cov = 0;
    int32_t len = 0;
    CGraphnode *prev_node = NULL;
    std::vector<GSeg> exons;
    std::vector<float> exon_cov;
    double exon_coverage = 0;

    uint32_t ref_start = 0;
    if (bundle) ref_start = (uint32_t)bundle->start;

    // Start after the source node and track the exon that we are building.
    int32_t first_node_idx = 0;
    if (!path[0]) first_node_idx = 1;

    bool at_first_exon = true;
    bool past_last_exon = false;

    // Go through the path nodes, adding each one to its exon.
    for (int32_t i = first_node_idx; i < count_of(path) - 1; i++) {
      // A node that no earlier path used makes this a new isoform (included =
      // false).
      if (!covered_path_bits[path[i]]) {
        included = false;
        covered_path_bits[path[i]] = 1;
      }
      // The same for an edge that no earlier path used.
      if (i) {
        int32_t *edge_pos
            = context.graph.path_bit_layout.find_edge_bit(path[i - 1], path[i]);
        if (edge_pos && !covered_path_bits[*edge_pos]) {
          included = false;
          covered_path_bits[*edge_pos] = 1;
        }
      }

      CGraphnode *node = context.graph.nodes[path[i]].get();

      // With a guide we skip the nodes outside it, so that they add no
      // coverage.
      if (guide && (node->end < guide->start || past_last_exon)) {
        prev_node = node;
        continue;
      }

      // Coverage that this path takes from the node.
      double used_coverage = node_coverage[path[i]] * node_flux[i]
          * (node->end - node->start + 1);
      // The flow takes the same share of the node's raw coverage: the node is
      // depleted by a share (node_coverage *= 1 - node_flux below), so its raw
      // / deduplicated ratio does not change. Without --sq3 it is 1.
      double node_raw_ratio = 1;
      if (node->cell_cov[context.cell] > 0)
        node_raw_ratio = node->cell_cov_raw[context.cell]
                         / node->cell_cov[context.cell];

      uint32_t exon_start = node->start;
      uint32_t exon_end = node->end;

      // Clip the first and last exon to the guide and remove the share of the
      // coverage that falls outside it.
      if (guide) {
        double first_cut_fraction = 0;
        double last_cut_fraction = 0;
        if (at_first_exon) {
          if (exon_start < guide->start && bundle) {
            double right_coverage = 0;

            // Cumulative bp coverage
            double left_coverage = get_region_coverage(context.strand,
                node->start - ref_start, guide->start - 1 - ref_start,
                bundle->bp_coverage);
            if (node->end > guide->start)
              right_coverage = get_region_coverage(context.strand,
                  guide->start - ref_start, node->end - ref_start,
                  bundle->bp_coverage);

            if (left_coverage)
              first_cut_fraction = left_coverage / (left_coverage
                  + right_coverage);
          }
          exon_start = guide->start;
        }
        if (guide->end <= exon_end || i == count_of(path) - 2) {
          past_last_exon = true;
          if (guide->end < node->end && bundle) {
            double left_coverage = 0;

            // Cumulative bp coverage
            if (node->start < guide->end) {
              left_coverage = get_region_coverage(context.strand,
                                                  node->start - ref_start,
                                                  guide->end - ref_start,
                                                  bundle->bp_coverage);
            }
            double right_coverage = get_region_coverage(context.strand,
                guide->end + 1 - ref_start, node->end - ref_start,
                bundle->bp_coverage);

            if (right_coverage)
              last_cut_fraction = right_coverage / (left_coverage
                  + right_coverage);
          }
          exon_end = guide->end;
        }
        if (first_cut_fraction || last_cut_fraction) {
          used_coverage -= (first_cut_fraction + last_cut_fraction)
              * used_coverage;
          node_flux[i] *= (1 - first_cut_fraction - last_cut_fraction);
        }
      }

      // Take the used coverage from the node (never below zero).
      node_coverage[path[i]] *= (1 - node_flux[i]);
      if (node_coverage[path[i]] < 0) node_coverage[path[i]] = 0;

      // This is a new exon
      if (!prev_node || at_first_exon || node->start > prev_node->end + 1) {
        if (prev_node && !at_first_exon) { // compute exon coverage
          exon_coverage /= exons.back().end - exons.back().start + 1;
          float exon_coverage_f = (float)exon_coverage;
          exon_cov.push_back(exon_coverage_f);
          exon_coverage = 0;
        }
        GSeg exon(exon_start, exon_end);
        exons.push_back(exon);
        at_first_exon = false;
      } else if (!at_first_exon) exons.back().end = exon_end;

      len += exon_end - exon_start + 1;

      cov += used_coverage;
      raw_cov += used_coverage * node_raw_ratio;
      exon_coverage += used_coverage;

      prev_node = node;
    }

    // Add last exon coverage
    if (prev_node) { // compute exon coverage
      exon_coverage /= exons.back().end - exons.back().start + 1;
      float exon_coverage_f = (float)exon_coverage;
      exon_cov.push_back(exon_coverage_f);
    }
    if (len) cov /= len;
    if (len) raw_cov /= len;

    // Store the transcript: always for a guide, otherwise only if it has
    // coverage and is at least min_transcript_len long.
    if (guide || (cov && len >= min_transcript_len)) {
      if (is_first_transcript) { gene_count++; }
      double transcript_coverage = cov;

      // We drop it if it equals a prediction that another cell already stored.
      bool keep_prediction = true;

      if (guide && guide->uptr) {
        set_guide_status(guide, GBST_STORED);
      // Check if prediction overlaps any previous predictions from same graph
      } else if (prev_prediction_count && count_of(predictions) > prev_prediction_count) {
        for (int32_t p = prev_prediction_count; p < count_of(predictions); p++) {
          // Check if this prediction is equal or included in a previous one
          double old_cov = predictions[p]->cov;
          if (merge_if_equal_prediction(predictions[p].get(), exons, exon_cov)) {
            // The merged raw ratio: the mean of both, weighted by coverage
            double added_cov = predictions[p]->cov - old_cov;
            if (added_cov > 0 && cov > 0)
              predictions[p]->raw_ratio = (old_cov * predictions[p]->raw_ratio
                  + added_cov * (raw_cov / cov)) / predictions[p]->cov;
            keep_prediction = false;
            break;
          }
        }
      }

      if (keep_prediction) {
        // Store the transcript prediction for this cell
        std::unique_ptr<CPrediction> p = std::make_unique<CPrediction>(
            gene_count - 1, guide, context.cell, exons[0].start,
            exons.back().end, transcript_coverage, context.strand, len);

        p->exons = exons;
        if (guide) {
          if (guide->exons.Count() == 1) exon_cov[0] = transcript_coverage;
        }
        p->exon_cov = exon_cov;
        p->path_nodes = path;
        if (cov > 0) p->raw_ratio = raw_cov / cov;

        predictions.push_back(std::move(p));
        is_first_transcript = false;
      }
    } else cov = 0;

    return (cov);
  }

} // namespace ribbon
