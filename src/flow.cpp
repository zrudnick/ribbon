#include "ribbon.h"
#include "GThreads.h"
#include <map>
#include <algorithm>

// Import globals from main program:
extern float isoform_fraction;
// See is_path_contradicted (off with --keep-contradicted-paths)
extern bool drop_contradicted_paths;
// See drop_shadowed_paths (off with --keep-shadowed-paths)
extern bool drop_shadowed;
// -c, applied later to the merged predictions (see drop_weak_predictions)
extern float min_read_coverage;

// Extracts transcripts from a bundle's splice graph (splice_graph.cpp) and its
// transfrags (transfrag.cpp). A transfrag is a read-supported path (nodes,
// path_bits, abundance pooled and per cell); a path is a source-to-sink node
// list that becomes a transcript; pushing flow along a path subtracts
// ("depletes") what it carries from the transfrags that support it. For each
// path node, the in/out abundance is that of the transfrags entering from the
// left or leaving to the right (path_* count only path-compatible ones,
// total_* all of them). low_cov_flags marks the nodes and edges to avoid
// because their pooled coverage is negligible
// (mark_low_coverage_nodes_and_edges). find_transcripts is the entry point,
// called once per cell and graph: guides claim their evidence first
// (assign_flow_to_guides), then extract_transcripts_greedily repeatedly builds
// the best supported path from the node with the most remaining coverage,
// stores it (store_transcript) and removes its evidence.

namespace ribbon {

  // Logistic function l_min + (l_max - l_min) / (1 + exp(-k * (x - x0))), used to
  // turn a log-coverage into a relative threshold (defaults: f(10) ~ 0.065,
  // f(50) ~ 0.159, f(100) ~ 0.26, approaching l_max = 0.3).
  float f_logistic(float x, float k = 0.03347, float x0 = 43.60,
                   float l_min = 0.05, float l_max = 0.3) {
    return l_min + (l_max - l_min) / (1.0 + (float)exp(-k * (x - x0)));
  }

  // The node next to node in transfrag tf (step 1: toward the sink, -1: toward
  // the source) if tf has the edge between them, else -1.
  int32_t next_node_in_transfrag(const FlowContext& context,
                                 const CTransfrag& tf, int32_t node,
                                 int32_t step) {
    for (int32_t j = 0; j < count_of(tf.nodes); j++) {
      if (tf.nodes[j] != node) continue;
      int32_t k = j + step;
      if (k < 0 || k >= count_of(tf.nodes)) return -1;
      int32_t *edge_pos = context.graph.path_bit_layout.find_edge_bit(
          std::min(node, tf.nodes[k]), std::max(node, tf.nodes[k]));
      return edge_pos && tf.path_bits[*edge_pos] ? tf.nodes[k] : -1;
    }
    return -1;
  }

  // True if the reads of all cells contradict two consecutive branch choices
  // of path: no transfrag has both chosen edges, but some transfrag has one of
  // them and takes another edge at the other branch node. Such paths are
  // mostly chimeras of real isoforms (pooled abundance, which flow does not
  // deplete, so the cell order does not matter).
  bool is_path_contradicted(const FlowContext& context,
                            const std::vector<int32_t>& path,
                            int32_t seed_node) {
    // The edges chosen where the walk had several neighbors, in path order
    // (from the seed the walk went toward the source, then toward the sink).
    struct Branch { int32_t from, to, edge_bit; };
    std::vector<Branch> branches;
    int32_t sink = context.graph.node_count - 1;
    bool toward_source = true;
    for (int32_t k = 1; k < count_of(path); k++) {
      int32_t a = path[k - 1], b = path[k];
      if (a == seed_node) toward_source = false;
      if (!a || b == sink) continue;
      const std::vector<int32_t>& neighbors = toward_source
          ? context.graph.nodes[b]->parent : context.graph.nodes[a]->child;
      if (count_of(neighbors) < 2) continue;
      int32_t *edge_pos = context.graph.path_bit_layout.find_edge_bit(a, b);
      if (edge_pos) branches.push_back({a, b, *edge_pos});
    }

    // For each pair of consecutive choices e1 = (u1, v1), e2 = (u2, v2), add
    // up the transfrags that link them and those that contradict them.
    for (int32_t k = 1; k < count_of(branches); k++) {
      const Branch& e1 = branches[k - 1];
      const Branch& e2 = branches[k];
      float link = 0, contradict = 0;

      // Through u2: with e1, then e2 (link) or another edge out of u2
      for (int32_t t : context.graph.nodes[e2.from]->transfrag_ids) {
        const CTransfrag& tf = *context.graph.transfrag[t];
        if (!tf.path_bits[e1.edge_bit]) continue;
        if (tf.path_bits[e2.edge_bit]) link += tf.abundance;
        else {
          int32_t next = next_node_in_transfrag(context, tf, e2.from, 1);
          if (next >= 0 && next != e2.to) contradict += tf.abundance;
        }
      }

      // Through v1: with e2 but entering v1 by another edge than e1
      for (int32_t t : context.graph.nodes[e1.to]->transfrag_ids) {
        const CTransfrag& tf = *context.graph.transfrag[t];
        if (tf.path_bits[e1.edge_bit] || !tf.path_bits[e2.edge_bit]) continue;
        int32_t prev = next_node_in_transfrag(context, tf, e1.to, -1);
        if (prev >= 0 && prev != e1.from) contradict += tf.abundance;
      }

      if (link <= epsilon && contradict > epsilon) return true;
    }
    return false;
  }

  // Extracts transcripts greedily in one cell (recursive): builds a path
  // through seed_node, pushes flow along it and stores it if it carries enough
  // flow, then repeats from the node with the most coverage left.
  void extract_transcripts_greedily(const FlowContext& context,
                                    int32_t seed_node,
                                    int32_t& gene_count,
                                    bool is_first_transcript,
                                    PredictionList& predictions,
                                    std::vector<float>& node_coverage,
                                    BitVec& transfrag_on_path,
                                    BitVec& used_node, float prev_best_coverage,
                                    BitVec& covered_path_bits,
                                    BitVec& low_cov_flags,
                                    int32_t prev_prediction_count) {
    // Graph and cell data shared through the context.

    // Start a path at the seed node and extend it to the source and the sink.
    std::vector<int32_t> path;
    path.push_back(seed_node);
    BitVec current_path_bits(context.graph.node_count
        + context.graph.edge_count);
    current_path_bits[seed_node] = 1;
    transfrag_on_path.reset();

    float flux = 0;
    std::vector<float> node_flux;

    // Only check if the path contains any non-empty transfrags, otherwise do
    // not compute flow
    bool path_empty = true;
    // Set if the path took a weakly supported junction from other cells' reads
    // (see check_weak_borrow); its flow is still pushed, but it is not stored
    bool weak_borrow = false;
    // The same if the reads contradict two of its choices (is_path_contradicted)
    bool contradicted = false;
    // V6
    if (extend_path(context, seed_node, true, path, current_path_bits,
                              node_coverage, path_empty, weak_borrow,
                              low_cov_flags)) {
      path.push_back(0); // include the source node
      // Back to source adds the nodes at the end to avoid pushing the list all
      // the time
      std::reverse(path.begin(), path.end());

      // V6
      if (extend_path(context, seed_node, false, path, current_path_bits,
                              node_coverage, path_empty, weak_borrow,
                              low_cov_flags)) {
        if (drop_contradicted_paths && !context.credit_pass)
          contradicted = is_path_contradicted(context, path, seed_node);
        // If we found at least one non-empty node on the path, we compute the
        // flow
        if (!path_empty)
          flux = push_flow_and_deplete_transfrags(context, path,
                                                  transfrag_on_path, node_flux,
                                                  current_path_bits);
      }
    }

    // Store the path as a transcript if it carried flow, borrowed no weak
    // junction and is not contradicted by the reads. The seed node is used up
    // if the path was already covered or its coverage is below
    // isoform_fraction times the best so far.
    if (flux > epsilon && !weak_borrow && !contradicted) {
      bool included = true;
      double cov = store_transcript(context, prev_prediction_count, predictions, path,
                                    node_flux, node_coverage, gene_count,
                                    is_first_transcript, included,
                                    covered_path_bits);

      if (included || cov < isoform_fraction * prev_best_coverage) {
        used_node[seed_node] = 1;
        seed_node = 0;
        prev_best_coverage = 0;
      } else if (cov > prev_best_coverage) prev_best_coverage = cov;

    } else {
      used_node[seed_node] = 1;
      seed_node = 0;
      prev_best_coverage = 0;
    }

    // Pick the next seed: the unused node with the most remaining coverage, or
    // the one with the most transfrag abundance passing through it.
    int32_t max_abundance_node = 0;
    float max_node_abundance = 0;
    // Node coverages:
    for (int32_t i = 1; i < context.graph.node_count - 1; i++)
      if (!used_node[i]) {
        if (node_coverage[i] > node_coverage[seed_node]) {
          if (!low_cov_flags[i]) // minimum cell coverage to consider A3
            seed_node = i;
        }
        // Only nodes that are not low coverage count their transfrag abundance
        // (low covered predictions could also be filtered out later).
        if (!low_cov_flags[i]) {
          float through_abundance = 0;
          const std::vector<int32_t>& node_transfrags
              = transfrags_of_cell(context, i);
          int32_t transfrag_count = count_of(node_transfrags);
          for (int32_t j = 0; j < transfrag_count; j++) {
            // Maybe we won't consider all transfrags in the future
            int32_t t = node_transfrags[j];
            if (!context.graph.transfrag[t]->path_bits[0]
                && !context.graph.transfrag[t]->path_bits[
                       context.graph.node_count - 1]) {
              through_abundance
                  += context.graph.transfrag[t]->cell_abundance[context.cell];
            }
          }
          if (through_abundance > max_node_abundance) {
            max_node_abundance = through_abundance;
            max_abundance_node = i;
          }
        }
      }

    // V24
    if (max_abundance_node && node_coverage[max_abundance_node] > ERROR_PERC) {
      seed_node = max_abundance_node;
    }

    // Repeat from the new seed while enough coverage is left.
    if (node_coverage[seed_node] >= 1
        || node_coverage[seed_node] > ERROR_PERC) {
      path.clear();
      node_flux.clear();
      extract_transcripts_greedily(context, seed_node, gene_count,
                                   is_first_transcript, predictions, node_coverage,
                                   transfrag_on_path, used_node,
                                   prev_best_coverage, covered_path_bits,
                                   low_cov_flags, prev_prediction_count);
    }
  }

  // Runs the de novo extraction of extract_transcripts_greedily on this cell
  // as if there were no guides, and returns the transcripts it finds. The
  // cell's transfrag abundances are restored afterwards, so the guides still
  // see all the evidence.
  PredictionList extract_cell_denovo_predictions(const FlowContext& context, 
                                                int32_t seed_node, 
                                                int32_t gene_count,
                                                std::vector<float> node_coverage, 
                                                BitVec current_path_bits, 
                                                BitVec& low_cov_flags) {
    // Save the cell's abundances. The links that guides made to the source and
    // sink are no evidence of this cell: we count them as pooled evidence
    // only, so that they do not end a path where this cell's reads go on.
    std::vector<float> saved_abundance;
    for (const auto& transfrag : context.graph.transfrag) {
      saved_abundance.push_back(transfrag->cell_abundance[context.cell]);
      if (transfrag->made_by_guide)
        transfrag->cell_abundance[context.cell] = 0;
    }

    // Mark the pass, so that contradicted paths are still found here
    FlowContext credit_context = context;
    credit_context.credit_pass = true;

    PredictionList denovo_predictions;
    if (node_coverage[seed_node] >= 1
        || node_coverage[seed_node] > ERROR_PERC) {
      BitVec transfrag_on_path(count_of(context.graph.transfrag));
      BitVec used_node(context.graph.node_count + context.graph.edge_count);
      extract_transcripts_greedily(credit_context, seed_node, gene_count, true,
                                   denovo_predictions, node_coverage,
                                   transfrag_on_path, used_node, 0,
                                   current_path_bits, low_cov_flags, 0);
    }

    for (int32_t t = 0; t < count_of(context.graph.transfrag); t++)
      context.graph.transfrag[t]->cell_abundance[context.cell]
          = saved_abundance[t];
    return denovo_predictions;
  }

  // True if the de novo transcript prediction has the intron chain of guide
  // (for a single-exon guide: overlaps it).
  bool is_guide_found_by_denovo(GffObj *guide, const CPrediction *prediction) {
    int32_t exon_count = guide->exons.Count();
    if (count_of(prediction->exons) != exon_count) return false;
    if (exon_count == 1)
      return prediction->exons[0].start <= guide->exons[0]->end
             && guide->exons[0]->start <= prediction->exons[0].end;
    for (int32_t e = 0; e < exon_count; e++) {
      if (e > 0 && prediction->exons[e].start != guide->exons[e]->start)
        return false;
      if (e + 1 < exon_count && prediction->exons[e].end != guide->exons[e]->end)
        return false;
    }
    return true;
  }

  // Sets this cell's credit for each guide to what de novo found for it: the
  // guide is a model, and the cell supports it exactly where
  // extract_cell_denovo_predictions finds its intron chain. Call it right after
  // the guides have claimed their flow. The coverage the flow gave the cell is
  // kept in flow_cov, for guides that no cell ends up with.
  void credit_guides_found_by_denovo(
      const FlowContext& context, const std::vector<CGuide>& matched_guides,
      std::vector<GffObj *>& guides, const PredictionList& predictions,
      const std::vector<int32_t>& guide_prediction,
      const PredictionList& denovo_predictions) {
    for (const CGuide& matched_guide : matched_guides) {
      int32_t p = guide_prediction[matched_guide.g];
      if (p < 0) continue;
      std::vector<float>& cell_cov = predictions[p]->cell_cov;
      cell_cov.resize(count_of(context.bundle->cell_name), 0.0);
      std::vector<float>& flow_cov = predictions[p]->flow_cov;
      flow_cov.resize(count_of(context.bundle->cell_name), 0.0);

      // The cell keeps the coverage the flow gave it only in flow_cov.
      flow_cov[context.cell] = cell_cov[context.cell];
      cell_cov[context.cell] = 0;
      for (const auto& denovo : denovo_predictions)
        if (is_guide_found_by_denovo(guides[matched_guide.g], denovo.get()))
          cell_cov[context.cell] = std::max(cell_cov[context.cell],
                                           (float)denovo->cov);
    }
  }

  // Clears (in cell_low_cov_flags) the low-coverage flags of the touching nodes
  // next to node i, on the left or on the right, when that neighbor has enough
  // coverage next to the coverage of i: its signal is then real for this cell.
  void rescue_low_coverage_neighbors(const FlowContext& context, 
                                     int32_t i,
                                     float node_coverage,
                                     float neighbor_coverage,
                                     float rescue_threshold, bool leftwards,
                                     const BitVec& low_cov_flags,
                                     BitVec& cell_low_cov_flags) {
    bool rescue_by_ratio = node_coverage > CHI_WIN
        && neighbor_coverage > rescue_threshold * node_coverage;
    bool rescue_by_minimum = node_coverage <= CHI_WIN
        && neighbor_coverage > SCELL_COV; // v3
    bool rescue = rescue_by_ratio || rescue_by_minimum;

    // True if node j is flagged and touches the node on its right (leftwards)
    // or on its left (rightwards).
    auto continues_flagged_run = [&](int32_t j) {
      if (leftwards)
        return j > 0 && low_cov_flags[j]
               && context.graph.nodes[j]->end
                   == context.graph.nodes[j + 1]->start;
      return j < context.graph.node_count - 1 && low_cov_flags[j]
             && context.graph.nodes[j]->start
                 == context.graph.nodes[j - 1]->end + 1;
    };

    int32_t j = leftwards ? i - 1 : i + 1;
    if (low_cov_flags[j] && rescue) {
      while (continues_flagged_run(j)) {
        cell_low_cov_flags[j] = 0;
        j += leftwards ? -1 : 1;
      }
    }
  }

  // The fraction of the pooled abundance of a transfrag from the source (or to
  // the sink) that this cell gets: the cell's share of the pooled coverage of
  // the node next to the source (or sink). Without reads there, the cell gets
  // none of it. The links of a node add up to the pooled abundance over cells.
  float cell_share_of_terminal_transfrag(const FlowContext& context,
                                         const CTransfrag *transfrag) {
    float share = 1;
    int32_t node_total = count_of(transfrag->nodes);
    if (node_total < 2) return share;

    std::vector<int32_t> next_to_terminal;
    if (transfrag->nodes[0] == 0)
      next_to_terminal.push_back(transfrag->nodes[1]);
    if (transfrag->nodes.back() == context.graph.node_count - 1)
      next_to_terminal.push_back(transfrag->nodes[node_total - 2]);

    for (int32_t node_id : next_to_terminal) {
      const CGraphnode *node = context.graph.nodes[node_id].get();
      // The source and sink have no coverage, so they keep the full abundance
      if (node->cov > 0)
        share = std::min(share, node->cell_cov[context.cell] / node->cov);
    }
    return share;
  }

  // Computes this cell's coverage of each node (node_coverage), rescues weak
  // nodes next to well covered ones (cell_low_cov_flags), and makes the
  // transfrags of the source and sink carry their full abundance in this cell.
  // Returns the node to start extracting from: the one with the most
  // transfrag abundance passing through it, else the one with the most
  // coverage.
  int32_t compute_cell_node_coverages(const FlowContext& context,
                                      const BitVec& low_cov_flags,
                                      BitVec& cell_low_cov_flags,
                                      std::vector<float>& node_coverage) {
    int32_t seed_node = 0; // node with maximum coverage
    // Node with maximum transcripts going through it
    int32_t max_abundance_node = 0;
    float max_node_abundance = 0;

    for (int32_t i = 0; i < context.graph.node_count; i++) {
      CGraphnode *current_node = context.graph.nodes[i].get();
      node_coverage.push_back(0.0);

      if (i) { // for all nodes but the source
        if (i < context.graph.node_count - 1 && current_node->len()) {
          // Node coverages are computed based on cell coverages
          // sink also has 0 coverage; node_coverage's have single cell
          // coverages
          node_coverage[i]
              = current_node->cell_cov[context.cell] / current_node->len();

          // Check if we can rescue some of the points
          // not a lowcov point -> check previous points ; What about the cases
          // where i is not covered but nearby nodes are quite high in
          // expression
          if (!low_cov_flags[i]) {
            // The threshold is 0.2 at coverage 100, 0.1 at 1000, 0.02 at
            // 10,000 and approaches 0.01 after that.
            float rescue_threshold = f_logistic(log10(node_coverage[i]),
                                                -2.7212954278522243,
                                                2.9063342233922893, 0.01,
                                                0.21612903225806462);

            rescue_low_coverage_neighbors(context, i, node_coverage[i],
                node_coverage[i - 1], rescue_threshold, true, low_cov_flags,
                cell_low_cov_flags);
            // The next node can be the sink, which has no length
            CGraphnode *next_node = context.graph.nodes[i + 1].get();
            float next_coverage = next_node->len()
                ? next_node->cell_cov[context.cell] / next_node->len() : 0;
            rescue_low_coverage_neighbors(context, i, node_coverage[i],
                next_coverage, rescue_threshold, false, low_cov_flags,
                cell_low_cov_flags);
          }

          if (!cell_low_cov_flags[i] || node_coverage[i] > SCELL_COV) {
            if (node_coverage[i] > node_coverage[seed_node]) seed_node = i;
            float through_abundance = 0;
            int32_t transfrag_count = count_of(current_node->transfrag_ids);
            for (int32_t j = 0; j < transfrag_count; j++) {
              // Maybe we won't consider all transfrags in the future
              int32_t t = current_node->transfrag_ids[j];
              if (!context.graph.transfrag[t]->path_bits[0]
                  && !context.graph.transfrag[t]->path_bits[
                       context.graph.node_count - 1]) {
                through_abundance
                    += context.graph.transfrag[t]->cell_abundance[context.cell];
              }
            }
            if (through_abundance > max_node_abundance) {
              max_node_abundance = through_abundance;
              max_abundance_node = i;
            }
          }
        }
      }

      // The transfrags that go from the source or to the sink carry this
      // cell's share of their pooled abundance
      if (i == 0 || i == context.graph.node_count - 1) {
        for (int32_t t : current_node->transfrag_ids) {
          context.graph.transfrag[t]->cell_abundance[context.cell]
              = context.graph.transfrag[t]->abundance
                * cell_share_of_terminal_transfrag(
                      context, context.graph.transfrag[t].get());
        }
      }
    }

    // Prefer the node with the most transfrag abundance passing through it.
    if (max_abundance_node && node_coverage[max_abundance_node] > ERROR_PERC) {
      seed_node = max_abundance_node;
    }
    return seed_node;
  }

  // A flagged junction is rescued for a cell with at least this many reads of
  // its own on it, and at least this share of the cell's strongest sibling
  // junction (the cell's version of the weak junction rule in flag_weak_edges)
  const float min_rescued_junction_reads = 2;

  const float min_rescued_junction_share = 0.1;

  // The abundance in this cell of the transfrags of node that cross the edge at
  // edge_pos.
  float cell_abundance_across_edge(const FlowContext& context,
                                   const CGraphnode *node, int32_t edge_pos) {
    float abundance = 0;
    for (int32_t t : node->transfrag_ids)
      if (context.graph.transfrag[t]->path_bits[edge_pos])
        abundance += context.graph.transfrag[t]->cell_abundance[context.cell];
    return abundance;
  }

  // Clears (in cell_low_cov_flags) the low-coverage flag of the junctions that
  // this cell has its own reads on. The flags come from the pooled coverage of
  // all cells, so a junction that is minor in the pool, but main in this cell,
  // would be closed to the cell otherwise.
  void rescue_junctions_with_cell_support(const FlowContext& context,
                                          BitVec& cell_low_cov_flags) {
    for (int32_t i = 1; i < context.graph.node_count - 1; i++) {
      const CGraphnode *from = context.graph.nodes[i].get();
      for (int32_t child : from->child) {
        const CGraphnode *to = context.graph.nodes[child].get();
        // Only the junctions: not the sink, not nodes that touch each other
        if (child == context.graph.node_count - 1 || to->start <= from->end + 1)
          continue;
        int32_t *edge_pos
            = context.graph.path_bit_layout.find_edge_bit(i, child);
        if (!edge_pos || !cell_low_cov_flags[*edge_pos]) continue;

        float support = cell_abundance_across_edge(context, from, *edge_pos);
        if (support < min_rescued_junction_reads) continue;

        // The strongest junction of this cell leaving from and entering to
        float strongest = 0;
        for (int32_t sibling : from->child) {
          int32_t *pos
              = context.graph.path_bit_layout.find_edge_bit(i, sibling);
          if (pos && sibling != context.graph.node_count - 1)
            strongest = std::max(strongest,
                cell_abundance_across_edge(context, from, *pos));
        }
        for (int32_t parent : to->parent) {
          int32_t *pos
              = context.graph.path_bit_layout.find_edge_bit(parent, child);
          if (pos && parent > 0)
            strongest = std::max(strongest, cell_abundance_across_edge(
                context, context.graph.nodes[parent].get(), *pos));
        }
        if (support >= min_rescued_junction_share * strongest)
          cell_low_cov_flags[*edge_pos] = 0;
      }
    }
  }

  // Lists, for each node, its transfrags that have abundance in this cell, in
  // the order of the node's transfrags. A cell has reads on only some of the
  // transfrags of the pooled graph, and its abundances only go down from here
  // (the copy run of the de novo extraction restores them to what they are
  // now), so the list holds every transfrag that can have abundance later.
  std::vector<std::vector<int32_t>> list_cell_transfrags(
      const FlowContext& context) {
    std::vector<std::vector<int32_t>> cell_transfrags(context.graph.node_count);
    for (int32_t i = 0; i < context.graph.node_count; i++)
      for (int32_t t : context.graph.nodes[i]->transfrag_ids)
        if (context.graph.transfrag[t]->cell_abundance[context.cell])
          cell_transfrags[i].push_back(t);
    return cell_transfrags;
  }

  // Extracts the transcripts of one cell from one graph into
  // bundle->predictions: computes the node coverage of this cell (rescuing weak
  // nodes with signal), lets guides claim evidence, extracts the rest greedily.
  // Returns gene_count.
  // True if transfrag tf fits a path (position_on_path: each node's index on
  // the path, -1 if not on it; first_node, last_node: the path's first and
  // last real nodes): its nodes are on the path, and each of its edges joins
  // nodes that are next to each other on the path (a gap between mates may
  // span several path nodes). It may run past the path's ends into touching
  // nodes (the path's end exon is shorter than the reads), not through a
  // splice.
  bool is_transfrag_on_whole_path(const FlowContext& context,
                                  const CTransfrag& tf,
                                  const std::vector<int32_t>& position_on_path,
                                  int32_t first_node, int32_t last_node) {
    auto touches = [&](int32_t a, int32_t b) {
      return context.graph.nodes[a]->end + 1 == context.graph.nodes[b]->start;
    };

    // Skip the overhangs before the first and after the last path node.
    int32_t first = 0, last = count_of(tf.nodes) - 1;
    while (first < last && position_on_path[tf.nodes[first]] < 0
           && touches(tf.nodes[first], tf.nodes[first + 1]))
      first++;
    if (first > 0 && tf.nodes[first] != first_node) return false;
    while (last > first && position_on_path[tf.nodes[last]] < 0
           && touches(tf.nodes[last - 1], tf.nodes[last]))
      last--;
    if (last < count_of(tf.nodes) - 1 && tf.nodes[last] != last_node)
      return false;

    int32_t prev_position = -1;
    for (int32_t j = first; j <= last; j++) {
      int32_t position = position_on_path[tf.nodes[j]];
      if (position < 0) return false;
      if (j > first) {
        int32_t *edge_pos = context.graph.path_bit_layout.find_edge_bit(
            tf.nodes[j - 1], tf.nodes[j]);
        if (edge_pos && tf.path_bits[*edge_pos]
            && position != prev_position + 1)
          return false;
      }
      prev_position = position;
    }
    return true;
  }

  // Called once all cells of the graph have their paths: takes a cell off a
  // novel path (its own copy of the prediction is dropped) when the cell's
  // reads do not single it out: each transfrag of the cell that fits the path
  // also fits another path of the graph with more of the cell's abundance. Uses the cell abundances from before the flow
  // (original_cell_abundance) and the paths of all cells and guides.
  void drop_shadowed_paths(
      const FlowContext& context,
      const std::vector<std::vector<float>>& original_cell_abundance,
      int32_t first_prediction) {
    PredictionList& predictions = context.bundle->predictions;
    int32_t prediction_count = count_of(predictions);
    int32_t sink = context.graph.node_count - 1;

    // The distinct paths of the graph's predictions, and the intron chain of
    // each: paths with the same chain (they differ only at their ends) are
    // one unit, so that a path is never shadowed by a copy of itself. A
    // single-exon path is a unit of its own.
    std::map<std::vector<int32_t>, int32_t> index_of_path;
    std::vector<int32_t> path_of_prediction(prediction_count, -1);
    std::vector<const std::vector<int32_t> *> paths;
    for (int32_t p = first_prediction; p < prediction_count; p++) {
      if (!predictions[p] || predictions[p]->path_nodes.empty()) continue;
      auto inserted = index_of_path.emplace(predictions[p]->path_nodes,
                                            count_of(paths));
      if (inserted.second) paths.push_back(&inserted.first->first);
      path_of_prediction[p] = inserted.first->second;
    }
    if (count_of(paths) < 1) return;
    std::map<std::vector<int32_t>, int32_t> index_of_chain;
    std::vector<int32_t> chain_of_path(count_of(paths));
    for (int32_t k = 0; k < count_of(paths); k++) {
      std::vector<int32_t> chain; // intron start and end coordinates
      const std::vector<int32_t>& path = *paths[k];
      for (int32_t j = 1; j < count_of(path); j++) {
        int32_t a = path[j - 1], b = path[j];
        if (!a || b == sink) continue;
        const CGraphnode& left = *context.graph.nodes[a];
        const CGraphnode& right = *context.graph.nodes[b];
        if (left.end + 1 < right.start) {
          chain.push_back(left.end);
          chain.push_back(right.start);
        }
      }
      if (chain.empty()) chain = path; // single exon: the path itself
      chain_of_path[k] = index_of_chain.emplace(chain,
          count_of(index_of_chain)).first->second;
    }
    int32_t chain_count = count_of(index_of_chain);

    // Only chains that will survive the coverage cutoff can shadow another:
    // their copies' coverage summed over the cells (as merge_cell_copies
    // does) reaches min_read_coverage, or a guide has the chain.
    std::vector<double> chain_coverage(chain_count, 0);
    std::vector<bool> chain_is_guide(chain_count, false);
    for (int32_t p = first_prediction; p < prediction_count; p++) {
      if (path_of_prediction[p] < 0 || !predictions[p]->flag) continue;
      int32_t u = chain_of_path[path_of_prediction[p]];
      chain_coverage[u] += predictions[p]->cov;
      if (predictions[p]->t_eq) chain_is_guide[u] = true;
    }
    auto can_shadow = [&](int32_t u) {
      return chain_is_guide[u] || chain_coverage[u] >= min_read_coverage;
    };

    // Which read transfrags (no source or sink links, none made by guides)
    // fit which chain (any of its paths)
    int32_t transfrag_count = count_of(context.graph.transfrag);
    std::vector<std::vector<int32_t>> transfrags_of_path(chain_count);
    std::vector<std::vector<int32_t>> paths_of_transfrag(transfrag_count);
    std::vector<int32_t> position_on_path(context.graph.node_count, -1);
    for (int32_t k = 0; k < count_of(paths); k++) {
      const std::vector<int32_t>& path = *paths[k];
      for (int32_t j = 0; j < count_of(path); j++) position_on_path[path[j]] = j;
      for (int32_t t = 0; t < transfrag_count; t++) {
        const CTransfrag& tf = *context.graph.transfrag[t];
        if (tf.made_by_guide || tf.nodes.empty() || !tf.nodes[0]
            || tf.nodes.back() == sink)
          continue;
        if (is_transfrag_on_whole_path(context, tf, position_on_path, path[1],
                                       path[count_of(path) - 2])) {
          int32_t u = chain_of_path[k];
          if (paths_of_transfrag[t].empty() || paths_of_transfrag[t].back() != u)
            if (std::find(paths_of_transfrag[t].begin(),
                          paths_of_transfrag[t].end(), u)
                == paths_of_transfrag[t].end()) {
              transfrags_of_path[u].push_back(t);
              paths_of_transfrag[t].push_back(u);
            }
        }
      }
      for (int32_t node : path) position_on_path[node] = -1;
    }

    // The abundance of a cell's transfrags that fit a chain (cached)
    std::map<std::pair<int32_t, int32_t>, float> support_cache;
    auto support = [&](int32_t k, int32_t cell) {
      auto it = support_cache.find({k, cell});
      if (it != support_cache.end()) return it->second;
      float total = 0;
      for (int32_t t : transfrags_of_path[k])
        total += original_cell_abundance[t][cell];
      support_cache[{k, cell}] = total;
      return total;
    };

    // Take each cell off the novel paths its reads do not single out.
    for (int32_t p = first_prediction; p < prediction_count; p++) {
      CPrediction *prediction = predictions[p].get();
      if (!prediction || !prediction->flag || prediction->t_eq
          || path_of_prediction[p] < 0 || prediction->cell < 0)
        continue;
      int32_t k = chain_of_path[path_of_prediction[p]];
      int32_t cell = prediction->cell;
      float own = support(k, cell);
      // Without any read transfrag of the cell on the whole path we cannot
      // tell: the flow may have used the cell's share of the source and sink
      // links (set in compute_cell_node_coverages), so the path is kept.
      if (own <= epsilon) continue;
      bool shadowed = true;
      {
        // Shadowed only if each of the cell's transfrags on the path is also
        // on a path with more support in the cell
        for (int32_t t : transfrags_of_path[k]) {
          if (original_cell_abundance[t][cell] <= epsilon) continue;
          bool on_better_path = false;
          for (int32_t q : paths_of_transfrag[t])
            if (q != k && can_shadow(q) && support(q, cell) > own) {
              on_better_path = true;
              break;
            }
          if (!on_better_path) {
            shadowed = false;
            break;
          }
        }
      }
      if (shadowed) prediction->flag = false;
    }
  }

  int32_t find_transcripts(const FlowContext& shared_context,
                       int32_t gene_count,
                       std::vector<CGuide>& matched_guides,
                       std::vector<GffObj *>& guides,
                       std::vector<int32_t>& guide_prediction,
                       BitVec& low_cov_flags) {
    // Graph and cell data shared through the context, with the transfrags of
    // this cell once we know them.
    FlowContext context = shared_context;

    PredictionList& predictions = context.bundle->predictions;

    // Remember how many predictions exist, to merge this cell's copies later.
    int32_t prev_prediction_count = count_of(predictions);

    // Compute this cell's coverage of each node, rescue weak nodes next to
    // well covered ones, and find the seed node and transfrag abundances.
    std::vector<float> node_coverage; // node coverages

    // Low-coverage flags of this cell: a copy that we clear for rescued nodes.
    BitVec cell_low_cov_flags(low_cov_flags);

    int32_t seed_node = compute_cell_node_coverages(context, low_cov_flags,
        cell_low_cov_flags, node_coverage);
    rescue_junctions_with_cell_support(context, cell_low_cov_flags);

    std::vector<std::vector<int32_t>> cell_transfrags
        = list_cell_transfrags(context);
    context.cell_transfrags = &cell_transfrags;

    // Set up the path bookkeeping; the guides claim their evidence first.
    BitVec transfrag_on_path(count_of(context.graph.transfrag));
    BitVec current_path_bits(context.graph.node_count
        + context.graph.edge_count);

    bool is_first_transcript = true;

    if (count_of(matched_guides)) {
      // What de novo would find in this cell, to credit the guides with.
      PredictionList denovo_predictions = extract_cell_denovo_predictions(
          context, seed_node, gene_count, node_coverage, current_path_bits,
          cell_low_cov_flags);
      seed_node = assign_flow_to_guides_and_pick_seed(context, matched_guides,
          gene_count, predictions, node_coverage, transfrag_on_path, current_path_bits,
          is_first_transcript, guides, guide_prediction);
      credit_guides_found_by_denovo(context, matched_guides, guides, predictions,
                                    guide_prediction, denovo_predictions);
    }

    // Extract the remaining transcripts greedily if enough coverage is left.
    if (node_coverage[seed_node] >= 1
        || node_coverage[seed_node] > ERROR_PERC) {
      BitVec used_node(context.graph.node_count + context.graph.edge_count);
      // A3
      extract_transcripts_greedily(context, seed_node, gene_count,
                                   is_first_transcript, predictions, node_coverage,
                                   transfrag_on_path, used_node, 0,
                                   current_path_bits, cell_low_cov_flags,
                                   prev_prediction_count);
    }

    return (gene_count);
  }

} // namespace ribbon
