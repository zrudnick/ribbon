#include "ribbon.h"

extern float multi_read_frac;
// minimum number for a transcript to be printed
extern int32_t min_transcript_len;
// See drop_shadowed_paths in flow.cpp (off with --keep-shadowed-paths)
extern bool drop_shadowed;

// The per-bundle assembly driver. infer_transcripts runs once for each bundle:
// it turns the reads into read groups and pre-graphs of groups, builds one
// splice graph per stranded pre-graph (splice_graph.cpp), fills the graphs with
// transfrags (transfrag.cpp) and extracts the transcripts (flow.cpp). A strand
// is '-' (negative), '+' (positive) or '.' (unknown); bp_coverage[strand] is
// the per-base coverage of the reads of that strand, with bp_coverage['.']
// holding the total (see get_region_coverage in coverage.cpp). state
// (AssemblyState, graph.h) carries the groups and pre-graphs between the steps.

namespace ribbon {

  // Adds a single-exon prediction for [start, end] that is not tied to a cell
  // or a strand, with the coverage per bp of the bundle, if it is longer than
  // min_transcript_len. Returns whether it was added.
  bool add_single_exon_prediction(BundleData *bundle, int32_t gene_count,
                                  uint32_t start, uint32_t end) {
    int32_t len = end - start + 1;
    if (len <= min_transcript_len) return false;

    int32_t ref_start = bundle->start;
    double cov = get_region_coverage('.', start - ref_start, end - ref_start,
                                     bundle->bp_coverage) / len;
    std::unique_ptr<CPrediction> p = std::make_unique<CPrediction>(
        gene_count - 1, nullptr, -1, start, end, cov, '.', len);
    p->exons.push_back(GSeg(start, end));
    p->exon_cov.push_back((float)cov);
    bundle->predictions.push_back(std::move(p));
    return true;
  }

  // Predicts single-exon transcripts from the coverage alone in unknown-strand
  // pre-graphs (state.pre_graph_list['.']): each segment is split at its trim
  // points and every piece above min_transcript_len is a prediction. Returns
  // the gene count.
  int32_t predict_unstranded_pre_graphs(BundleData *bundle,
                                 AssemblyState& state) {
    int32_t ref_start = bundle->start;
    std::vector<GffObj *>& guides = bundle->guides;
    BpCoverage& bp_coverage = bundle->bp_coverage;

    int32_t gene_count = 0;
    int32_t g = 0;

    // These are neutral pre-graphs that do not overlap any signed reads
    for (int32_t b = 0; b < count_of(state.pre_graph_list['.']); b++) {
      PreGraphSegment *first_segment
          = state.pre_graph_segment_list['.']
                [state.pre_graph_list['.'][b]->start_node].get();
      bool guide_ovlp = false;
      // Multi-exon guides and guides that end before the pre-graph are skipped.
      auto is_multi_exon_or_before_pre_graph = [&](int32_t idx) {
        return guides[idx]->exons.Count() > 1
               || guides[idx]->end < first_segment->start;
      };
      while (g < count_of(guides) && is_multi_exon_or_before_pre_graph(g)) {
        g++;
      }
      // Now guides[g] ends at or after the start of pre-graph b
      if (g < count_of(guides)
          && guides[g]->start <= first_segment->end) guide_ovlp = true;

      // && (count_of(guides) || adaptive || pre_graph_list[1][b]->len >=
      // min_transcript_len)) { // there might be small transfrags that are
      // worth showing, but here we are ignoring them
      auto& neutral_pre_graph = state.pre_graph_list['.'][b];
      bool has_coverage = neutral_pre_graph->cov;
      bool has_few_multi_reads = has_coverage
          && (neutral_pre_graph->multi / neutral_pre_graph->cov)
                 <= multi_read_frac * (1 - ERROR_PERC);
      if (has_coverage && (has_few_multi_reads || guide_ovlp)) {
        // The pre-graph might contain multiple fragments of a transcript but
        // since we
        // don't know the complete structure -> print only the pieces that are
        // well represented
        PreGraphSegment *curr_pre_graph_segment = first_segment;
        bool made_prediction = false;
        while (curr_pre_graph_segment != NULL) {
          // The first segment that makes no prediction yet counts as a gene
          if (!made_prediction) { gene_count++; }

          // Split the node at its trim points and predict a transcript
          // for each piece that is long enough.
          std::vector<CTrimPoint> trim_point;
          // Sign should not matter as we are in a totally neutral zone
          find_all_trim_points(ref_start, '-', curr_pre_graph_segment->start,
                               curr_pre_graph_segment->end, bp_coverage,
                                   trim_point);

          uint32_t prediction_start = curr_pre_graph_segment->start;
          uint32_t prediction_end = curr_pre_graph_segment->end;

          for (int32_t i = 0; i < count_of(trim_point); i++) {
            if (trim_point[i].pos) {
              if (trim_point[i].start) {
                if (add_single_exon_prediction(bundle, gene_count, prediction_start,
                                               trim_point[i].pos - CHI_WIN))
                  made_prediction = true;
                prediction_start = trim_point[i].pos;
              } else {
                if (add_single_exon_prediction(bundle, gene_count, prediction_start,
                                               trim_point[i].pos))
                  made_prediction = true;
                prediction_start = trim_point[i].pos + CHI_WIN;
              }
            }
          }

          // The last piece, after the final trim point.
          if (add_single_exon_prediction(bundle, gene_count, prediction_start,
                                         prediction_end))
            made_prediction = true;
          curr_pre_graph_segment = curr_pre_graph_segment->next_node;
        }
      }
    }

    return gene_count;
  }

  // A pre-graph gets a splice graph if it has coverage and is long enough with
  // few enough multi-hit reads, or if it overlaps a supported guide (there
  // might be small transfrags from source to sink worth keeping).
  static bool is_pre_graph_worth_processing(const PreGraph& pre_graph,
                                          int32_t overlapping_guides) {
    if (!pre_graph.cov) return false;
    bool has_few_multi_reads
        = (pre_graph.multi / pre_graph.cov) <= multi_read_frac;
    bool is_long_enough = pre_graph.len >= min_transcript_len;
    return (has_few_multi_reads && is_long_enough) || overlapping_guides;
  }

  // Builds the splice graph and the transfrag index of each pre-graph of one
  // strand ('-' or '+') that is worth processing.
  static void get_stranded_splice_graphs(
      BundleData *bundle, AssemblyState& state, GraphTables& tables,
      std::vector<GuideBoundary>& guide_boundaries, char strand) {
    std::vector<GffObj *>& guides = bundle->guides;
    PreGraphList& pre_graphs = state.pre_graph_list[strand];
    PreGraphSegmentList& segments = state.pre_graph_segment_list[strand];

    // The data of each graph g of this strand (see GraphTables)
    SegmentNodesList& pre_graph_segment_to_graph_nodes
        = tables.pre_graph_segment_to_graph_nodes[strand];
    std::vector<SpliceGraph>& splice_graphs = tables.splice_graphs[strand];

    pre_graph_segment_to_graph_nodes.resize(count_of(segments));
    if (!count_of(pre_graphs)) return;

    // One graph slot per pre-graph; only the pre-graphs worth processing get a
    // graph (the others keep a node count of 0).
    int32_t pre_graph_count = count_of(pre_graphs);
    splice_graphs.resize(pre_graph_count);

    int32_t guide_cursor = 0;
    for (int32_t i = 0; i < pre_graph_count; i++) {
      SpliceGraph& graph = splice_graphs[i];

      int32_t n_overlapping_guides = count_guides_overlapping_pre_graph(
          guides, strand, segments[pre_graphs[i]->start_node]->start,
          segments[pre_graphs[i]->last_node_id]->end, guide_cursor);
      graph.edge_count += 2 * n_overlapping_guides;

      // Guides are kept at any length, but novel transcripts still need the
      // minimum transcript length even next to a guide, to stay conservative.
      if (is_pre_graph_worth_processing(*pre_graphs[i], n_overlapping_guides)) {
        graph.node_count = build_splice_graph(bundle, state, tables, strand, i,
            guide_boundaries);
        if (!graph.node_count) continue;

        graph.transfrag_index = build_transfrag_index(graph);
      }
    }
  }

  // Builds the splice graph and the transfrag index of each stranded pre-graph
  // that is worth processing. The graphs are left in tables.
  void get_splice_graphs(BundleData *bundle, AssemblyState& state,
                    GraphTables& tables,
                    std::vector<GuideBoundary>& guide_boundaries) {
    for (char strand : {'-', '+'}) {
      get_stranded_splice_graphs(bundle, state, tables, guide_boundaries,
          strand);
    }
  }

  void free_pre_graphs(AssemblyState& state) {
    for (char strand : {'-', '.', '+'}) {
      state.pre_graph_segment_list[strand].clear();
      state.pre_graph_list[strand].clear();
    }
  }

  // Adds every read, or pair of reads, to the transfrags of its graph.
  void fill_transfrags(BundleData *bundle, AssemblyState& state,
                       GraphTables& tables) {
    ReadList& reads = bundle->reads;

    for (int32_t n = 0; n < count_of(reads); n++) {
      float single_count = reads[n]->read_count;
      for (int32_t j = 0; j < count_of(reads[n]->pair_idx); j++) {
        int32_t np = reads[n]->pair_idx[j];
        if (np > -1) {
          single_count -= reads[n]->pair_count[j];
          if (n < np) {
            assign_fragment_to_transfrags(bundle, state, tables, n, np,
                                          reads[n]->pair_count[j]);
          }
        }
      }
      if (single_count > epsilon) {
        assign_fragment_to_transfrags(bundle, state, tables, n, -1,
                                      single_count);
      }
    }
  }

  void free_read_groups(AssemblyState& state) {
    state.group.clear();
  }

  // Per graph: matches the guides, cleans the transfrags and extracts the
  // transcripts of each cell, then frees the graph. Returns the updated
  // gene_count.
  int32_t extract_transcripts_from_graphs(BundleData *bundle,
                                          GraphTables& tables,
                                          std::vector<int32_t>& guide_prediction, 
                                          int32_t gene_count) {
    std::vector<GffObj *>& guides = bundle->guides;
    int32_t cell_count = count_of(bundle->cell_name);
    // The data of each graph g on a strand s (see GraphTables)

    // Process each graph: guides, transfrags, low-coverage flags and the
    // transcripts of each cell.
    for (char strand : {'-', '+'}) {
      for (int32_t i = 0; i < count_of(tables.splice_graphs[strand]); i++) {
        SpliceGraph& graph = tables.splice_graphs[strand][i];
        if (graph.node_count) {
          // Include source to guide starts links
          std::vector<CGuide> matched_guides;
          if (count_of(guides))
            process_reference_guides(graph, strand, matched_guides, bundle);

          // Process transfrags to eliminate noise, and set compatibilities,
          // and node memberships
          process_transfrags(graph, matched_guides, bundle);

          // Find transcripts now
          // mark nodes/junctions in the graph that are below the (10%)
          // coverage threshold
          BitVec low_cov_flags(graph.node_count + graph.edge_count);
          FlowContext context = {graph, strand, bundle, 0};
          mark_low_coverage_nodes_and_edges(context, low_cov_flags);

          // Keep the cells' abundances from before the flow, and where this
          // graph's predictions start, for drop_shadowed_paths.
          std::vector<std::vector<float>> original_cell_abundance;
          if (drop_shadowed)
            for (const auto& transfrag : graph.transfrag)
              original_cell_abundance.push_back(transfrag->cell_abundance);
          int32_t first_prediction = count_of(bundle->predictions);

          // For all single cells in the bundle compute their transcripts
          // independently
          for (int32_t c = 0; c < cell_count; c++) {
            context.cell = c;
            gene_count = find_transcripts(context, gene_count,
                                          matched_guides, guides,
                                          guide_prediction, low_cov_flags);
          }

          // Take cells off the paths that other paths explain better, then
          // forget the graph paths.
          if (drop_shadowed)
            drop_shadowed_paths(context, original_cell_abundance,
                                first_prediction);
          for (int32_t p = first_prediction; p < count_of(bundle->predictions);
               p++)
            if (bundle->predictions[p]) bundle->predictions[p]->path_nodes.clear();
        }
        // Clean up what can be cleaned
        graph.transfrag_index = TransfragIndex();
        graph.transfrag.clear();
        graph.nodes.clear();
      }

      // Final clean up of the strand: the map from segments to graph nodes
      tables.pre_graph_segment_to_graph_nodes[strand].clear();
    }

    return gene_count;
  }

  // Assembles the transcripts of one bundle (called by process_bundle in
  // threads.cpp): finishes the junction support and coverage, then builds
  // the splice graphs and extracts the transcripts.
  void infer_transcripts(BundleData *bundle) {
    compute_junction_support(bundle);

    ReadList& reads = bundle->reads;
    JunctionList& junction = bundle->junction;
    std::vector<GffObj *>& guides = bundle->guides;

    std::vector<int32_t> guide_prediction;
    std::vector<GuideBoundary> guide_boundaries;
    compute_guide_boundaries(guides, guide_boundaries, guide_prediction);

    // Form groups on strands: all groups below are kept per strand ('-', '.'
    // unknown, '+')
    AssemblyState state(junction, count_of(reads));
    propagate_junction_support(bundle, state);
    assign_reads_to_groups(bundle, state);
    merge_close_groups(bundle, state);
    form_pre_graphs_from_groups(state);

    // Predict the unknown-strand pre-graphs from their coverage.
    int32_t gene_count = predict_unstranded_pre_graphs(bundle, state);

    // For the stranded pre-graphs: build a splice graph for each, put the reads
    // into the transfrags of the graphs, and extract the transcripts.
    bool has_stranded_groups = state.start_group['-'] != NULL
                               || state.start_group['+'] != NULL;
    if (has_stranded_groups) {
      GraphTables tables;
      get_splice_graphs(bundle, state, tables, guide_boundaries);
      free_pre_graphs(state);
      fill_transfrags(bundle, state, tables);
      free_read_groups(state);
      gene_count = extract_transcripts_from_graphs(bundle, tables,
                                                   guide_prediction, gene_count);
    } else {
      // No stranded groups: no graphs, only free what is left.
      free_pre_graphs(state);
      free_read_groups(state);
    }
  }

} // namespace ribbon
