#include "ribbon.h"
#include "GThreads.h"

// Reads at what distance should be considered part of separate bundles
extern uint32_t min_bundle_distance;

// Joins read groups into pre-graphs (chains of PreGraphSegments), one per
// stranded locus, that each become a splice graph.

namespace ribbon {

  // Starts a fresh pre-graph (and its first PreGraphSegment) from `group`, appends
  // it to `pre_graph_list`/`pre_graph_segments`, and returns the new pre-graph's
  // index.
  int32_t start_new_pre_graph(PreGraphList& pre_graph_list, CGroup *group,
                            PreGraphSegmentList& pre_graph_segments) {
    int32_t node_id = count_of(pre_graph_segments);
    int32_t pre_graph_idx = count_of(pre_graph_list);
    pre_graph_segments.push_back(std::make_unique<PreGraphSegment>(
        group->start, group->end, group->cov_sum, node_id));
    pre_graph_list.push_back(std::make_unique<PreGraph>(
        group->end - group->start + 1, group->cov_sum, group->multi, node_id,
        node_id));

    return (pre_graph_idx);
  }

  // Adds group to a pre-graph: far enough past the last segment it becomes a new
  // PreGraphSegment, otherwise it extends that segment. Updates the length,
  // coverage and multi-hit totals of the pre-graph.
  void add_group_to_pre_graph(CGroup *group, PreGraph *pre_graph,
                            PreGraphSegmentList& pre_graph_segments,
                            uint32_t merge_distance) {
    PreGraphSegment *last_node
        = pre_graph_segments[pre_graph->last_node_id].get();
    int32_t new_node_id = count_of(pre_graph_segments);

    // Group after last pre_graph_segments
    if (group->start > last_node->end + merge_distance) {
      pre_graph_segments.push_back(std::make_unique<PreGraphSegment>(
        group->start, group->end, group->cov_sum, new_node_id));
      last_node->next_node = pre_graph_segments.back().get();
      pre_graph->last_node_id = new_node_id;
      pre_graph->len += group->end - group->start + 1;
      pre_graph->cov += group->cov_sum;
      pre_graph->multi += group->multi;

    // Group overlaps the last segment within min_bundle_distance
    } else { 
      if (last_node->end < group->end) {
        pre_graph->len += group->end - last_node->end;
        last_node->end = group->end;
      }
      pre_graph->cov += group->cov_sum;
      pre_graph->multi += group->multi;
      last_node->cov += group->cov_sum;
    }
  }

  // Union-find merge for an unknown-strand group that overlaps the stranded
  // prev_group: their link ids, and those already equated to prev_group's,
  // end up with the same smallest link id in the root arrays.
  void propagate_strand_link_id(CGroup *prev_group, CGroup *group,
                                int32_t base_link_id,
                                std::vector<int32_t>& strand_link_id_root,
                                std::vector<int32_t>& link_id_root) {
    // Find the link id that prev_group's pre-graph is mapped to.
    int32_t mapped_link_id = strand_link_id_root[prev_group->link_id];
    if (mapped_link_id > -1) {
      while (strand_link_id_root[mapped_link_id] != -1
        && strand_link_id_root[mapped_link_id] != mapped_link_id) {
        mapped_link_id = strand_link_id_root[mapped_link_id];
      }
      int32_t saved_mapped_link_id = mapped_link_id;
      while (link_id_root[mapped_link_id] != mapped_link_id) {
        mapped_link_id = link_id_root[mapped_link_id];
      }
      strand_link_id_root[prev_group->link_id] = mapped_link_id;
      strand_link_id_root[saved_mapped_link_id] = mapped_link_id;

      // Make the larger of the two link ids point to the smaller one.
      if (mapped_link_id < base_link_id) {
        link_id_root[base_link_id] = mapped_link_id;
        group->link_id = mapped_link_id;
      } else if (base_link_id < mapped_link_id) {
        link_id_root[mapped_link_id] = base_link_id;
        strand_link_id_root[prev_group->link_id] = base_link_id;
      }
    } else {
      strand_link_id_root[prev_group->link_id] = base_link_id;
    }
  }

  // Returns the strand ('-', '.' or '+') whose current group starts first, so
  // the groups of the three strand tracks can be walked in genomic order as one
  // stream.
  char next_group_by_start(PerStrand<CGroup *>& curr_group) {
    char next_group = '-';
    if (curr_group['-'] != NULL) {
      if (curr_group['.'] != NULL) {
        if (curr_group['+'] != NULL) {
          char smaller_idx
            = curr_group['-']->start < curr_group['.']->start ? '-' : '.';

          next_group = curr_group[smaller_idx]->start < curr_group['+']->start
            ? smaller_idx : '+';
          return (next_group);
        } else {
          next_group = curr_group['-']->start < curr_group['.']->start
            ? '-' : '.';
          return (next_group);
        }
      } else {
        if (curr_group['+'] != NULL) {
          next_group = curr_group['-']->start < curr_group['+']->start
            ? '-' : '+';
          return (next_group);
        } else {
          return ('-');
        }
      }
    } else {
      if (curr_group['.'] != NULL) {
        if (curr_group['+'] != NULL) {
          next_group = curr_group['.']->start < curr_group['+']->start
              ? '.' : '+';
          return (next_group);
        } else {
          return ('.');
        }
      } else {
        return ('+');
      }
    }

    return (next_group);
  }

  // The part of other's coverage that falls in the stretch that group and other
  // share (at least one bp).
  float overlap_share(CGroup *group, CGroup *other) {
    uint32_t max_start = group->start > other->start ? group->start
                                                     : other->start;
    uint32_t min_end = group->end < other->end ? group->end 
                                               : other->end;

    // This can only happen if min_bundle_distance >0
    if (min_end < max_start) min_end = max_start;
    return other->cov_sum * (min_end - max_start + 1) / other->len();
  }

  // Resolves each group's link id (union-find root) and passes it on to the
  // stranded groups that an unknown-strand group overlaps, setting its neg_prop
  // (negative strand share). Then collects the groups of each pre-graph into a
  // PreGraph, chaining its covered stretches as PreGraphSegments.
  void form_pre_graphs_from_groups(AssemblyState& state) {

    // First do link_id assignment
    for (char strand : {'-', '.', '+'})
      state.curr_group[strand] = state.start_group[strand];
    PerStrand<CGroup *> prev_group;

    std::vector<int32_t> pos_strand_link_id_root;
    std::vector<int32_t> neg_strand_link_id_root;

    pos_strand_link_id_root.resize(count_of(state.link_id_root), -1);
    neg_strand_link_id_root.resize(count_of(state.link_id_root), -1);

    // Links an unknown-strand group to the groups of one stranded strand ('-'
    // or '+') that it overlaps: the previous one and the current
    // ones, which are consumed. The share of their coverage in the overlap is
    // added to share (for the positive strand only if the group already has a
    // negative share).
    auto link_overlapping_groups = [&](CGroup *unknown, char strand,
                                       std::vector<int32_t>& strand_root,
                                       float& share) {
      bool count_shares = strand == '-' || unknown->neg_prop;
      CGroup *&prev = prev_group[strand];
      CGroup *&curr = state.curr_group[strand];

      // Overlaps the previous group; this needs min_bundle_distance
      if (prev != NULL && unknown->start <= prev->end + min_bundle_distance) {
        propagate_strand_link_id(unknown, prev, prev->link_id, strand_root,
          state.link_id_root);
        if (count_shares) share += overlap_share(unknown, prev);
      }

      // Overlaps the current groups
      while (curr != NULL
        && unknown->start <= curr->end + min_bundle_distance
        && curr->start <= unknown->end + min_bundle_distance) {

        // Set smallest link_id for the group
        curr->link_id = find_root(state.link_id_root, curr->link_id);
        if (strand == '-') curr->neg_prop = 1;

        propagate_strand_link_id(unknown, curr, curr->link_id, strand_root,
                                 state.link_id_root);
        if (count_shares) share += overlap_share(unknown, curr);

        prev = curr;
        curr = curr->next_gr;
      }
    };

    // Each unknown-strand group remembers the share of the stranded groups it
    // overlaps (neg_prop), so that its reads can be split between the strands.

    // There are still groups to process
    while (state.curr_group['-'] != NULL || state.curr_group['.'] != NULL
      || state.curr_group['+'] != NULL) {

      // Gets the strand of curr_group with the left most begining
      char next_group = next_group_by_start(state.curr_group);
      CGroup *group = state.curr_group[next_group];

      // Set smallest link_id for curr_group[$next_group]
      group->link_id = find_root(state.link_id_root, group->link_id);

      if (next_group == '.') { // unknown strand group
        link_overlapping_groups(group, '-', neg_strand_link_id_root,
          group->neg_prop);

        float pos_prop = 0;
        link_overlapping_groups(group, '+', pos_strand_link_id_root, 
          pos_prop);

        if (pos_prop) {
          group->neg_prop /= (group->neg_prop + pos_prop);
        } else if (group->neg_prop)
          group->neg_prop = 1;
      } else if (next_group == '-') { // negative strand group
        group->neg_prop = 1;
      }

      prev_group[next_group] = state.curr_group[next_group];
      state.curr_group[next_group] = state.curr_group[next_group]->next_gr;
    }

    // Create the pre-graphs: a PreGraph collects the connected groups (same
    // link id).
    for (char strand : {'-', '.', '+'}) {
      state.curr_group[strand] = state.start_group[strand];
      prev_group[strand] = NULL;
    }

    for (char strand : {'-', '.', '+'}) {
      // For a given group id we get a certain pre_graph_segment id
      state.group_to_pre_graph_segment[strand].resize(count_of(state.group),
          -1);
    }

    // created_pre_graph_idx[link_id] is the index in pre_graph_list of the
    // PreGraph already created for that link id, or -1 if none yet. (Ids
    // are labels on the groups and have gaps; indices count the PreGraphs in
    // the order they are created.)
    std::vector<int32_t> created_pre_graph_idx;
    created_pre_graph_idx.resize(count_of(state.link_id_root), -1);

    // Adds a group to the pre-graph of a strand with this link id, creating
    // the pre-graph the first time, and remembers the pre-graph segment of the
    // group.
    auto add_to_pre_graph = [&](char strand, CGroup *group,
                              int32_t link_id) {
      int32_t pre_graph_idx = created_pre_graph_idx[link_id];
      if (pre_graph_idx > -1) { // pre-graph for group has been created before
        add_group_to_pre_graph(group,
                             state.pre_graph_list[strand][pre_graph_idx].get(),
                             state.pre_graph_segment_list[strand],
                             min_bundle_distance);
      } else { // create new pre-graph
        pre_graph_idx = start_new_pre_graph(
          state.pre_graph_list[strand], group,
          state.pre_graph_segment_list[strand]);
        created_pre_graph_idx[link_id] = pre_graph_idx;
      }
      state.group_to_pre_graph_segment[strand][group->grid]
        = state.pre_graph_list[strand][pre_graph_idx]->last_node_id;
    };

    // There are still groups to process
    while (state.curr_group['-'] != NULL || state.curr_group['.'] != NULL
      || state.curr_group['+'] != NULL) {

      // Next group based on starting position
      char next_group = next_group_by_start(state.curr_group);
      CGroup *group = state.curr_group[next_group];

      // Get group link_id; we need to redo this to ensure we equalize all
      // link_ids -> they could still be handled by
      // propagate_strand_link_id
      int32_t resolved_link_id = find_root(state.link_id_root, group->link_id);
      group->link_id = resolved_link_id;

      // Negative or positive strand pre-graph or unstranded pre-graph
      bool is_stranded_group = next_group == '-' || next_group == '+';
      bool is_unstranded_without_strand_link = next_group == '.'
          && neg_strand_link_id_root[resolved_link_id] == -1
          && pos_strand_link_id_root[resolved_link_id] == -1;
      if (is_stranded_group || is_unstranded_without_strand_link) {
        add_to_pre_graph(next_group, group, resolved_link_id);
        
      // Unknown strand : here is where we should compute positive and negative
      // proportions
      } else {
        if (neg_strand_link_id_root[resolved_link_id] != -1) {
          add_to_pre_graph('-', group,
              find_root(state.link_id_root,
                        neg_strand_link_id_root[resolved_link_id]));
        }
        if (pos_strand_link_id_root[resolved_link_id] != -1) {
          add_to_pre_graph('+', group,
              find_root(state.link_id_root,
                        pos_strand_link_id_root[resolved_link_id]));
        }
      }

      state.curr_group[next_group] = state.curr_group[next_group]->next_gr;
    }

    // Free the arrays that we no longer need.
    state.link_id_root.clear();
    pos_strand_link_id_root.clear();
    neg_strand_link_id_root.clear();
    created_pre_graph_idx.clear();
  }

} // namespace ribbon
