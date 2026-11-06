#include "ribbon.h"
#include "GThreads.h"

// Anchor length for junction to be considered well supported
extern uint32_t junction_support;
// Reads at what distance should be considered part of separate bundles
extern uint32_t min_bundle_distance;

// Places reads into groups (same-strand runs of overlapping exon coverage)
// and links the groups of read pairs.

namespace ribbon {

  // True if the first exon of the read is long enough or its first junction
  // has a strand.
  bool has_valid_first_exon(CReadAln& read) {
    return read.segs[0].len() >= junction_support
           || (count_of(read.juncs) && read.juncs[0]->strand != '.');
  }

  // True if the last exon of the read is long enough or its last junction has a
  // strand.
  bool has_valid_last_exon(CReadAln& read) {
    return read.segs.back().len() >= junction_support
           || (count_of(read.juncs) && read.juncs.back()->strand != '.');
  }

  // True if the mate was placed before the read as part of the same fragment: a
  // valid last exon of the mate (long enough or with a good junction) that is
  // close to a valid first exon of the read.
  bool is_placed_with_mate(ReadList& reads, int32_t read_idx, int32_t mate_idx,
                           uint32_t merge_distance) {
    if (mate_idx <= -1 || mate_idx >= read_idx) return false;
    CReadAln& mate = *reads[mate_idx];
    CReadAln& read = *reads[read_idx];
    bool mate_ends_before_read = mate.end < read.start;
    bool gap_is_small = mate.end + merge_distance > read.start;
    return mate_ends_before_read && gap_is_small && has_valid_last_exon(mate)
           && has_valid_first_exon(read);
  }

  // True if the mate comes after the read and the two form one fragment: a
  // valid last exon of the read that is close to a valid first exon of the
  // mate.
  bool continues_with_mate(ReadList& reads, int32_t read_idx, int32_t mate_idx,
                           uint32_t merge_distance) {
    if (mate_idx == -1 || mate_idx <= read_idx) return false;
    CReadAln& read = *reads[read_idx];
    CReadAln& mate = *reads[mate_idx];
    bool read_ends_before_mate = read.end < mate.start;
    bool gap_is_small = read.end + merge_distance > mate.start;
    return read_ends_before_mate && gap_is_small && has_valid_first_exon(mate);
  }

  // True if the mate was already placed in the groups (it comes before the read
  // and was not deleted).
  bool has_placed_mate(ReadList& reads, int32_t read_idx, int32_t mate_idx) {
    return mate_idx > -1 && reads[mate_idx]->nh && mate_idx < read_idx;
  }

  // False for an exon of the read that is too short, unless it has a good
  // junction.
  bool is_exon_kept(CReadAln& read, int32_t exon_idx) {
    if (read.segs[exon_idx].len() < junction_support) {
      if (exon_idx < count_of(read.juncs)) {
        if (read.juncs[exon_idx]->strand == '.') { // check first exon
          // See if first exon is big enough
          if (!exon_idx || read.juncs[exon_idx - 1]->strand == '.')
            return false;
        }
      } else if (count_of(read.juncs)) { // this is last exon
        if (read.juncs[exon_idx - 1]->strand == '.') return false;
      }
    }
    return true;
  }

  // Adds a new group with this extent, link id and coverage to the groups.
  CGroup *add_group(AssemblyState& state, uint32_t start, uint32_t end,
                    int32_t link_id, float cov_sum, float multi) {
    int32_t group_id = count_of(state.group);
    state.group.push_back(std::make_unique<CGroup>(start, end, link_id,
                                                   group_id, cov_sum, multi));
    state.merge_map.push_back(group_id);
    return state.group.back().get();
  }

  // Starts a new link id for the read, forgetting what was assigned to
  // last_link_id before. Returns it.
  int32_t start_new_link_id(AssemblyState& state,
                            PerStrand<int32_t>& last_link_id,
                            char strand) {
    last_link_id[strand] = state.link_id;
    state.link_id_root.push_back(state.link_id);
    return state.link_id++;
  }

  // Breaks the pairing link of a read and its mate (their ends need different
  // link ids) and gives read_link_id a fresh link id for this
  // strand (last_link_id if set, else the next free id).
  void split_read_pair_and_relabel_link_id(int32_t read_idx,
                                           int32_t mate_idx,
                                           int32_t pair_slot,
                                           char strand, ReadList& reads,
                                           int32_t& read_link_id,
                                           PerStrand<int32_t>& last_link_id,
                                           AssemblyState& state) {
    // Unlink the read and its mate.
    reads[read_idx]->pair_idx[pair_slot] = -1;
    for (int32_t j = 0; j < count_of(reads[mate_idx]->pair_idx); j++)
      if (reads[mate_idx]->pair_idx[j] == read_idx) {
        reads[mate_idx]->pair_idx[j] = -1;
        break;
      }

    // The read_link_id gets back the new link_id
    read_link_id = last_link_id[strand];
    if (read_link_id < 0) {
      last_link_id[strand] = state.link_id;
      read_link_id = state.link_id;
      state.link_id_root.push_back(state.link_id);
      state.link_id++;
    }
  }

  // Merges group2 into group1 (which comes first): group1 takes over group2's
  // end, coverage and multi and the smaller link id, skips group2 in its
  // chain and records the merge in merge_map so reads can be redirected.
  void merge_adjacent_groups(AssemblyState& state, CGroup *group1,
                             CGroup *group2) {
    // Get end of group (group1 is assumed to come before group2)
    group1->end = group2->end;

    // Get smallest link_id of group
    while (state.link_id_root[group1->link_id] != group1->link_id) {
      group1->link_id = state.link_id_root[group1->link_id];
    }
    while (state.link_id_root[group2->link_id] != group2->link_id) {
      group2->link_id = state.link_id_root[group2->link_id];
    }

    if (group1->link_id < group2->link_id) {
      state.link_id_root[group2->link_id] = group1->link_id;
    } else if (group1->link_id > group2->link_id) {
      state.link_id_root[group1->link_id] = group2->link_id;
      group1->link_id = group2->link_id;
    }

    // Add the coverage of group2 to group1 and remove group2 from the chain.
    group1->cov_sum += group2->cov_sum;
    // This is possible because group1->next_gr=group2
    group1->next_gr = group2->next_gr;

    state.merge_map[group2->grid] = group1->grid;

    group1->multi += group2->multi;

    // Delete group2
    state.group[group2->grid].reset();
  }

  // Merges the groups that follow the group and start at or before position
  // into it, then extends it to position.
  void extend_group_to(AssemblyState& state, CGroup *group,
                       uint32_t position) {
    CGroup *following_group = group->next_gr;
    while (following_group != NULL && position >= following_group->start) {
      merge_adjacent_groups(state, group, following_group);
      following_group = group->next_gr;
    }
    if (position > group->end) group->end = position;
  }

  // Moves curr_group to the group where the read starts (or the one before
  // it) and starts the walk there.
  void ReadSegmentPlacer::find_group_where_read_starts() {
    // While read start after the current group's end advance group
    auto advance_past_groups_before_read = [&]() {
      while (curr_group != NULL && reads[read_idx]->start > curr_group->end) {
        //-> might have more groups leaving from current group due to splicing
        prev_group = curr_group;
        curr_group = curr_group->next_gr;
      }
    };

    prev_group = NULL;
    advance_past_groups_before_read();

    // The read lies before curr_group, which is where the previous read was
    // placed. The chain has no back links, so we walk again from the start
    // of the chain to find the group before curr_group.
    if (prev_group == NULL && curr_group != NULL
        && reads[read_idx]->segs[0].end < curr_group->start
        && curr_group != state.start_group[strand]) {
      curr_group = state.start_group[strand];
      advance_past_groups_before_read();
    }

    // The curr_group is null only if we reached end of curr_group list.
    // Because curr_group is not NULL initially and it starts BEFORE read we
    // make sure read comes after curr_group. If there is no group before it
    // (prev_group is NULL) the read comes before the first group of the
    // chain, so we keep curr_group and the read gets a new group before it.
    if (prev_group != NULL
        && (curr_group == NULL
            || reads[read_idx]->segs[0].end < curr_group->start))
      curr_group = prev_group;

    matching_group = curr_group;
  }

  // For the first exon in a group: finds out if we saw this group of the
  // read before, and splits the pair if its link id didn't reach the group.
  // We need to split pairs here if link_id didn't reach this group: it means
  // there is a gap between these groups and no reads joining them.
  void ReadSegmentPlacer::check_pair_at_first_exon() {
    matching_group->grid = find_root(state.merge_map, matching_group->grid);

    for (int32_t g = 0; g < count_of(state.read_group[read_idx]); g++) {
      if (find_root(state.merge_map, state.read_group[read_idx][g])
          == matching_group->grid) {
        first_visit = false; // we did see this groups before
        break;
      }
    }

    // We only consider first exon here because the rest of the groups
    // need to get the same link_id if the read is still paired
    if (has_placed_mate(reads, read_idx, mate_idx)) {
      // Get smallest link_id
      int32_t resolved_link_id
          = find_root(state.link_id_root, matching_group->link_id);
      matching_group->link_id = resolved_link_id;

      // Pair link_id didn't reach this group
      if (resolved_link_id != read_link_id) {
        split_read_pair_and_relabel_link_id(read_idx, mate_idx, pair_slot,
                                            strand, reads, read_link_id,
                                            last_link_id, state);
      }
    }
  }

  // Makes the link id of the read and of the group the same: the smaller of
  // the two wins, unless the read is interrupted before this exon.
  void ReadSegmentPlacer::unify_link_ids() {
    // Get smallest link_id of group
    matching_group->link_id
        = find_root(state.link_id_root, matching_group->link_id);

    // Read link_id is different from group link_id
    if (read_link_id != matching_group->link_id) {
      // Bad junction before this exon
      if (exon_idx && reads[read_idx]->juncs[exon_idx - 1]->strand == '.') {
        // Before we equal read_link_id to matching_group->link_id, we need
        // to make sure that the read is not interrupted read_link_id for
        // next exons should still be the link_id of matching_group
        read_link_id = matching_group->link_id;
      } else {
        // Set group link_id to current read link_id
        if (read_link_id < matching_group->link_id) {
          state.link_id_root[matching_group->link_id] = read_link_id;
          matching_group->link_id = read_link_id;
        } else { // read link_id is bigger than group

          state.link_id_root[read_link_id] = matching_group->link_id;
          read_link_id = matching_group->link_id;
        }
      }
    }
  }

  // Adds the exon to the group it overlaps (or is close to): the group grows
  // to the exon, and takes the coverage and the link id of the read.
  void ReadSegmentPlacer::add_exon_to_matching_group() {
    // We need to split pairs here if link_id didn't reach this group
    if (!exon_idx) check_pair_at_first_exon();

    if (reads[read_idx]->nh > 1) {
      // Multi-hit reads add their coverage to the group's multi total.
      matching_group->multi += read_coverage
          * reads[read_idx]->segs[exon_idx].len();
    }

    // Extend the group to the exon, merging the groups it reaches.
    if (reads[read_idx]->segs[exon_idx].start < matching_group->start) {
      matching_group->start = reads[read_idx]->segs[exon_idx].start;
    }

    // Find end of new group
    extend_group_to(state, matching_group,
                    reads[read_idx]->segs[exon_idx].end);

    unify_link_ids();

    if (matching_group->grid != last_pushed_group_id) {
      // The read_group for read read_idx gets the id of group
      if (first_visit)
        state.read_group[read_idx].push_back(matching_group->grid);
      last_pushed_group_id = matching_group->grid;
    }

    // Coverage is different than number of reads
    matching_group->cov_sum += (reads[read_idx]->segs[exon_idx].end
                                - reads[read_idx]->segs[exon_idx].start
                                + 1) * read_coverage;
  }

  // The exon is at the end of the groups, or it doesn't overlap any group:
  // makes a new group for it, chained before matching_group.
  void ReadSegmentPlacer::add_exon_as_new_group() {
    // On the first exon we split the pair: no group overlaps the read,
    // so the link id cannot have reached this place.
    if (!exon_idx && has_placed_mate(reads, read_idx, mate_idx)) {
      split_read_pair_and_relabel_link_id(read_idx, mate_idx, pair_slot,
                                          strand, reads, read_link_id,
                                          last_link_id, state);
    // If exons are interrupted use new link_id
    } else if (exon_idx
               && reads[read_idx]->juncs[exon_idx - 1]->strand == '.') {
      // We have to use new link_id for read and we have no idea if
      // last_link_id was already used or not => we have to ignore
      // what was already assigned to last_link_id
      read_link_id = start_new_link_id(state, last_link_id, strand);
    }

    // Create a new group for the exon.
    float multi_cov = 0;
    if (reads[read_idx]->nh > 1)
      multi_cov = read_coverage * reads[read_idx]->segs[exon_idx].len();

    CGroup *new_group = add_group(state,
        reads[read_idx]->segs[exon_idx].start,
        reads[read_idx]->segs[exon_idx].end, read_link_id,
        reads[read_idx]->segs[exon_idx].len() * read_coverage, multi_cov);
    int32_t new_group_id = new_group->grid;
    // Link the new group into the chain; without a prev_group it becomes the
    // start of the chain.
    if (prev_group) prev_group->next_gr = new_group;
    else state.start_group[strand] = new_group;
    new_group->next_gr = matching_group;

    // Since group was created it clearly is new
    state.read_group[read_idx].push_back(new_group_id);
    last_pushed_group_id = new_group_id;
    matching_group = new_group;
  }

  // At the last exon of the read: if the mate continues the fragment, goes on
  // with the exons of the mate as if they were part of the read.
  void ReadSegmentPlacer::continue_with_mate_exons() {
    // Valid last exon and close
    if (continues_with_mate(reads, read_idx, mate_idx, merge_distance)) {
      // Next we will process read pair as it was part of the same fragment
      read_idx = mate_idx;
      mate_idx = -1;
      exon_count = count_of(reads[read_idx]->segs);
      first_visit = true;
      // We need to extend group
      if (reads[read_idx]->start > matching_group->end) {
        extend_group_to(state, matching_group, reads[read_idx]->start);
      }
      last_pushed_group_id = -1;
      exon_idx = -1;
    }
  }

  // The strand already has groups: each kept exon joins the group where
  // it fits or gets a new group.
  void ReadSegmentPlacer::place_in_existing_chain() {
    find_group_where_read_starts();

    // Now process each group of coordinates individually
    exon_count = count_of(reads[read_idx]->segs);
    last_pushed_group_id = -1;

    exon_idx = 0;
    while (exon_idx < exon_count) {
      // Skip an exon that is too short, unless it has a good junction.
      bool keep = is_exon_kept(*reads[read_idx], exon_idx);
      if (!keep && !exon_idx && prev_group) curr_group = prev_group;

      if (keep) {
        // Skip groups that are left behind
        // find the group where "exon" fits
        while (matching_group != NULL
               && reads[read_idx]->segs[exon_idx].start
               > matching_group->end) {
          prev_group = matching_group;
          matching_group = matching_group->next_gr;
        }

        // Read overlaps group or it's close enough -> it should get group
        // link_id
        if (matching_group
            && reads[read_idx]->segs[exon_idx].end + merge_distance
            >= matching_group->start) {
          add_exon_to_matching_group();
        // Read is at the end of groups, or read is not overlapping other
        // groups
        } else {
          add_exon_as_new_group();
        }

        // Process the mate's exons too if the pair forms one fragment.
        if (exon_idx == exon_count - 1) { // last segment in read
          continue_with_mate_exons();
        }
      }
      exon_idx++;
    }
  }

  // The strand has no groups yet: the read makes the first groups, one
  // for each kept exon, chained in order.
  void ReadSegmentPlacer::place_in_new_chain() {
    // We need to split pairs here because we have no overlap => link_id
    // didn't reach here for sure
    if (has_placed_mate(reads, read_idx, mate_idx)) {
      split_read_pair_and_relabel_link_id(read_idx, mate_idx, pair_slot,
                                          strand, reads, read_link_id,
                                          last_link_id, state);
      // Even if we split pairs we do not need to update single_count because
      // this read_coverage is taken care of
    }

    // Create a new group for each exon that is kept and chain them.
    exon_count = count_of(reads[read_idx]->segs);
    prev_group = NULL;
    float multi_cov = 0;
    // This will be adjusted by the length of each segment
    if (reads[read_idx]->nh > 1) multi_cov = read_coverage;
    exon_idx = 0;
    while (exon_idx < exon_count) {
      // Skip an exon that is too short, unless it has a good junction.
      bool keep = is_exon_kept(*reads[read_idx], exon_idx);

      if (!exon_idx || keep) {
        // Use a new link id if the read is interrupted before this exon.
        if (exon_idx && reads[read_idx]->juncs[exon_idx - 1]->strand == '.') {
          read_link_id = start_new_link_id(state, last_link_id, strand);
        }

        // Create the group and chain it to the previous one.
        CGroup *new_group = add_group(state,
            reads[read_idx]->segs[exon_idx].start,
            reads[read_idx]->segs[exon_idx].end, read_link_id,
            reads[read_idx]->segs[exon_idx].len() * read_coverage,
            reads[read_idx]->segs[exon_idx].len() * multi_cov);
        int32_t new_group_id = new_group->grid;
        if (prev_group != NULL) {
          prev_group->next_gr = new_group;
        } else {
          curr_group = new_group;
        }
        prev_group = new_group;
        if (keep) {
          // Add group to read so that we can remember what groups belong to a
          // read
          state.read_group[read_idx].push_back(new_group_id);

          // Process the mate's exons too if the pair forms one fragment.
          if (exon_idx == exon_count - 1) { // last segment in read
            // Valid last exon and close
            if (continues_with_mate(reads, read_idx, mate_idx,
                                    merge_distance)) {
              // Next we will process read pair as it was part of the same
              // fragment
              read_idx = mate_idx;
              mate_idx = -1;
              exon_count = count_of(reads[read_idx]->segs);
              // We need to extend group
              if (reads[read_idx]->segs[0].end > prev_group->end) {
                prev_group->end = reads[read_idx]->segs[0].end;
                // Add group to read so that we can remember what groups
                // belong to a read
                if (first_visit)
                  state.read_group[read_idx].push_back(new_group_id);
                // Coverage is different than number of reads
                prev_group->cov_sum += reads[read_idx]->segs[0].len()
                    * read_coverage;
              }
              exon_idx = 0;
            }
          }
        }
      }
      exon_idx++;
    }
  }

  ReadSegmentPlacer::ReadSegmentPlacer(ReadList& reads_in,
                                       AssemblyState& state_in,
                                       ReadPlacement placement,
                                       PerStrand<int32_t>& last_link_id_in)
      : reads(reads_in), state(state_in), pair_slot(placement.pair_slot),
        read_coverage(placement.read_coverage),
        strand(placement.strand), last_link_id(last_link_id_in),
        merge_distance(min_bundle_distance + long_intron_anchor),
        read_idx(placement.read_idx), mate_idx(placement.mate_idx),
        read_link_id(placement.read_link_id),
        // Start from the group where the previous read on this strand was
        // added.
        curr_group(state_in.curr_group[placement.strand]) {}

  void ReadSegmentPlacer::place() {
    // Skip a read that its mate already placed as part of the same fragment.
    // Valid last/first exon and close
    if (is_placed_with_mate(reads, read_idx, mate_idx, merge_distance)) {
      return;
    }

    if (curr_group != NULL) { // this type of group was created before
      place_in_existing_chain();
    } else { // create new group of this type
      place_in_new_chain();
    }

    // Remember where we are in this strand's group chain for the next read.
    state.curr_group[strand] = curr_group;

    if (state.start_group[strand] == NULL) {
      state.start_group[strand] = curr_group;
    }
  }

  // Places a read (and its mate) into the group chain of one strand:
  // each exon segment joins an overlapping group (merging the groups it
  // bridges) or starts a new CGroup.
  void assign_read_segments_to_group(ReadList& reads, AssemblyState& state,
                                     ReadPlacement placement,
                                     PerStrand<int32_t>& last_link_id) {
    ReadSegmentPlacer placer(reads, state, placement, last_link_id);
    placer.place();
  }

  // Places a read in the groups of its strand: for each pair we find the
  // shared link id (the mate's, if we saw it) and call
  // assign_read_segments_to_group, then we do it once for the unpaired part.
  void assign_read_to_groups(int32_t read_idx, ReadList& reads,
                             AssemblyState& state) {
    // Remembers the last link_id used for this read per strand, so that
    // we don't use different link_ids for the same read
    PerStrand<int32_t> last_link_id(-1);

    // Count of the read that is not part of any pair (we subtract the pairs).
    float single_count = reads[read_idx]->read_count;

    // Place the read once for each pair that it is part of.
    for (int32_t pair_slot = 0; pair_slot < count_of(reads[read_idx]->pair_idx);
         pair_slot++) {
      // Each pair is handled independently of the previous one, so we reset
      // strand every time: the mate may change it.
      char strand = reads[read_idx]->strand;

      // Check if we've seen read's pair and if yes get its read_link_id; at
      // the least get read's pair strand if available
      // pair read number
      int32_t mate_idx = reads[read_idx]->pair_idx[pair_slot];

      // Read pair still exists and it wasn't deleted -> false keep elminates
      // all these cases
      if (mate_idx > -1 && reads[mate_idx]->nh) {
        // See if we have the correct read strand
        // mate_strand is the strand of pair read
        char mate_strand = reads[mate_idx]->strand;
        // Different strands for read and pair
        if (strand != mate_strand) {
          // Read read_idx is on zero (neutral) strand, but pair has strand
          if (strand == '.') {
            strand = mate_strand; // assign strand of pair to read
          // Conflicting strands -> un-pair reads in the hope that one is right
          } else if (mate_strand != '.') {
            reads[read_idx]->pair_idx[pair_slot] = -1;
            // Also unpair read mate_idx to read_idx
            for (int32_t j = 0; j < count_of(reads[mate_idx]->pair_idx); j++)
              if (reads[mate_idx]->pair_idx[j] == read_idx) {
                reads[mate_idx]->pair_idx[j] = -1;
                break;
              }
            mate_idx = -1;
          }
        }

        if (mate_idx > -1) { // reads didn't get split

          // Update single count of read here
          single_count -= reads[read_idx]->pair_count[pair_slot];

          int32_t read_link_id = last_link_id[strand];

          // There is a pair and it came before the current read in sorted order
          // of position
          if (mate_idx < read_idx) {
            // First group of pair read is: read_group[mate_idx][0]

            int32_t mate_first_group = state.read_group[mate_idx][0];
            while (state.merge_map[mate_first_group] != mate_first_group) {
              mate_first_group = state.merge_map[mate_first_group];
            }
            state.read_group[mate_idx][0] = mate_first_group;

            // The read_link_id gets assigned the link_id of the pair's
            // first group --> this might not
            read_link_id
                = state.group[state.read_group[mate_idx][0]]->link_id;
            // Work if we split  a read into multiple groups
            // get smallest link_id
            while (state.link_id_root[read_link_id] != read_link_id) {
              read_link_id = state.link_id_root[read_link_id];
            }
            state.group[state.read_group[mate_idx][0]]->link_id
                = read_link_id;
          } else { // it's the first time we see the read in the fragment
            // We didn't use the link_id yet
            if (last_link_id[strand] < 0) {
              last_link_id[strand] = state.link_id;
              read_link_id = state.link_id;
              // Read might keep the new link_id so we need to add it to the
              // link_id_root
              state.link_id_root.push_back(state.link_id);
              state.link_id++;
            }
          }

          // Place the read, with the coverage of this pair, in the groups.
          double read_coverage = reads[read_idx]->pair_count[pair_slot];
          ReadPlacement placement = {read_idx, mate_idx, pair_slot,
                                     (float)read_coverage, strand,
                                     read_link_id};
          assign_read_segments_to_group(reads, state, placement,
                                        last_link_id);
        }
      }
    }

    // Place the part of the read count that is not in any pair.
    if (single_count > epsilon) { // our way of controlling for rounding errors
      int32_t read_link_id = last_link_id[reads[read_idx]->strand];
      if (read_link_id < 0) {
        read_link_id = state.link_id;
        last_link_id[reads[read_idx]->strand] = state.link_id;
        // Read might keep the new link_id so we need to add it to the
        // link_id_root
        state.link_id_root.push_back(state.link_id);
        state.link_id++;
      }
      double read_coverage = single_count;
      ReadPlacement placement = {read_idx, -1, -1, (float)read_coverage,
                                 reads[read_idx]->strand,
                                 read_link_id};
      assign_read_segments_to_group(reads, state, placement, last_link_id);
    }
  }

  // Unpairs a read from its mates, and unstrands a mate without junctions.
  void unpair_read(ReadList& reads, int32_t read_idx) {
    CReadAln& read = *reads[read_idx];
    for (int32_t p = 0; p < count_of(read.pair_idx); p++) {
      int32_t mate_idx = read.pair_idx[p];
      if (mate_idx > -1) {
        read.pair_idx[p] = -1;
        // Also unpair read np to n
        for (int32_t j = 0; j < count_of(reads[mate_idx]->pair_idx); j++)
          if (reads[mate_idx]->pair_idx[j] == read_idx) {
            reads[mate_idx]->pair_idx[j] = -1;
            break;
          }
        // Remove strand for pair read if single exon ? is strand
        // assigned before here ?
        if (!count_of(reads[mate_idx]->juncs)) reads[mate_idx]->strand = '.';
      }
    }
  }

  // True if any junction of the read has a strand.
  bool has_stranded_junction(const CReadAln& read) {
    for (int32_t j = 0; j < count_of(read.juncs); j++)
      if (read.juncs[j]->strand != '.') return true;
    return false;
  }

  // Decides if an unspliced, stranded read keeps its strand because a later
  // mate supports it: a mate without junctions on the same strand, or one with
  // a stranded, well supported junction. A later mate that does not support it
  // is unstranded.
  bool mates_keep_strand(ReadList& reads, int32_t read_idx, int32_t ref_start,
                         BpCoverage& bp_coverage) {
    CReadAln& read = *reads[read_idx];
    bool keep_strand = false;
    for (int32_t p = 0; p < count_of(read.pair_idx); p++) {
      int32_t mate_idx = read.pair_idx[p];
      if (mate_idx > -1 && read_idx < mate_idx) {
        if (!count_of(reads[mate_idx]->juncs)) {
          if (reads[mate_idx]->strand == read.strand) {
            keep_strand = true;
            break;
          }
        } else {
          for (int32_t j = 0; j < count_of(reads[mate_idx]->juncs); j++) {
            CJunction& junc = *(reads[mate_idx]->juncs[j]);
            if (junc.strand != '.'
                && is_junction_well_supported(junc, ref_start, bp_coverage)) {
              keep_strand = true;
              break;
            }
          }
          if (!keep_strand) reads[mate_idx]->strand = '.';
        }
      }
      if (keep_strand) break;
    }
    return keep_strand;
  }

  // The number of new fragments that a read adds: its count, minus the pairs
  // that were already counted with their first read (the mate was seen before
  // and was not deleted).
  double new_fragment_count(ReadList& reads, int32_t read_idx) {
    CReadAln& read = *reads[read_idx];
    double single_count = read.read_count;
    for (int32_t i = 0; i < count_of(read.pair_idx); i++) {
      if (read.pair_idx[i] != -1 && read_idx > read.pair_idx[i]
          && reads[read.pair_idx[i]]->nh) {
        single_count -= read.pair_count[i];
      }
    }
    return single_count;
  }

  // Reclassifies the junctions of each read (replacing those redirected by
  // propagate_junction_support, creating new ones if needed), unstrands reads
  // without reliable junctions, and groups them (assign_read_to_groups).
  void assign_reads_to_groups(BundleData *bundle, AssemblyState& state) {
    int32_t ref_start = bundle->start;
    ReadList& reads = bundle->reads;
    JunctionList& junction = bundle->junction;
    std::vector<GffObj *>& guides = bundle->guides;
    BpCoverage& bp_coverage = bundle->bp_coverage;

    // Go through the reads: fix their junctions, then assign them to groups.
    bool resort = false;
    int32_t original_junction_count = count_of(junction);
    for (int32_t read_idx = 0; read_idx < count_of(reads); read_idx++) {
      CReadAln& read = *(reads[read_idx]);

      // Check the junctions of the read; junctions_ok turns false at the
      // first bad one.
      bool junctions_ok = true;
      int32_t i = 0;

      while (i < count_of(read.juncs)) {
        CJunction& junc = *(read.juncs[i]);

        bool change_right = junc.end_redirect != 0;
        bool change_left = junc.start_redirect != 0;

        if (junc.strand != '.') {
          // Re-check a stranded junction: weak or redirected junctions
          // lose their strand.
          if (!is_junction_well_supported(junc, ref_start, bp_coverage))
            junc.rejected = true;
          if (change_left || change_right) junc.strand = '.';
          // We can't believe the junction because it's too low
          if (junc.rejected) junc.strand = '.';
        }

        if (junc.strand == '.') { // found a bad junction
          // Replace a redirected junction by the better one it points to,
          // or by a new one with the redirected coordinates.
          if (change_left || change_right) {
            uint32_t new_start = read.segs[i].end;
            uint32_t new_end = read.segs[i + 1].start;
            int32_t start_side_junc_idx = -1;
            int32_t end_side_junc_idx = -1;
            if (change_left) {
              start_side_junc_idx = junc.start_redirect;
              // Not a valid left support -> delete junction
              if (junction[start_side_junc_idx]->start_redirect) {
                new_start = read.segs[i].start - 1;
              } else {
                // Junction jk is good
                new_start = junction[start_side_junc_idx]->start;
                if (junction[start_side_junc_idx]->strand != '.')
                  start_side_junc_idx = -1;
              }
            }
            if (change_right) {
              end_side_junc_idx = junc.end_redirect;

              if (state.end_sorted_junctions[end_side_junc_idx]->end_redirect) {
                new_end = read.segs[i + 1].end + 1;
              } else {
                new_end = state.end_sorted_junctions[end_side_junc_idx]->end;
                if (state.end_sorted_junctions[end_side_junc_idx]->strand
                    != '.')
                  end_side_junc_idx = -1;
              }
            }
            // Junction inside read boundaries
            if (new_start >= read.segs[i].start
                && new_end <= read.segs[i + 1].end && new_start <= new_end) {
              read.segs[i].end = new_start; // adjust start
              read.segs[i + 1].start = new_end; // adjust end
              if (!junc.rejected && (start_side_junc_idx > 0
                                     || end_side_junc_idx > 0)) {
                // Use the junction that already has these coordinates, or add
                // one.
                CJunction *replacement = start_side_junc_idx > 0
                    ? find_junction_around(junction, start_side_junc_idx,
                                           false, new_start, new_end,
                                           read.strand)
                    : find_junction_around(state.end_sorted_junctions,
                                           end_side_junc_idx, true, new_start,
                                           new_end, read.strand);
                if (!replacement) {
                  replacement = find_added_junction(junction,
                      original_junction_count, new_start, new_end,
                      read.strand);
                }
                if (!replacement) {
                  resort = true;
                  junction.push_back(std::make_unique<CJunction>(
                      new_start, new_end, read.strand));
                  replacement = junction.back().get();
                  state.end_sorted_junctions.push_back(replacement);
                }
                read.juncs[i] = replacement;
              }
            }
          }

          // A read with a bad junction is unpaired the first time we see one;
          // a mate without junctions then loses its strand.
          if (junctions_ok) {
            unpair_read(reads, read_idx);
            junctions_ok = false;
          }
        } else if (count_of(guides)) { // need to remember boundary
          state.junction_start_positions.insert(junc.start);
          state.junction_end_positions.insert(junc.end);
        }
        i++;
      }

      // Read has bad junctions -> check if we need to unstrand it
      if (!junctions_ok) {
        if (read.strand != '.' && !has_stranded_junction(read))
          read.strand = '.';
      // Check if we need to unstrand current read due to future mapping
      } else if (!count_of(read.juncs) && read.strand != '.') {
        if (count_of(read.pair_idx)
            && !mates_keep_strand(reads, read_idx, ref_start, bp_coverage))
          read.strand = '.';
      }

      // Put the read in the groups of its strand.
      assign_read_to_groups(read_idx, reads, state);

      // Add the read to the fragment length and count the fragment once
      // (not for the second read of a pair).
      bundle->frag_len += read.len * read.read_count;
      // The same for the cell of the read.
      if (read.bc >= 0) {
        if (read.bc >= count_of(bundle->cell_frag_len))
          bundle->cell_frag_len.resize(read.bc + 1, 0.0);
        bundle->cell_frag_len[read.bc] += read.len * read.read_count;
      }
      double single_count = junctions_ok ? new_fragment_count(reads, read_idx)
                                         : read.read_count;
      if (single_count > epsilon) {
        bundle->num_fragments += single_count;
      }
    }

    if (resort) {
      sort_by(junction, compare_junction_by_position);
      sort_by(state.end_sorted_junctions, compare_junction_by_end);
    }
  }

  // Merges adjacent groups of each strand unless one side is a guide-derived
  // boundary (state.junction_start_positions/junction_end_positions); the gap
  // must be within min_bundle_distance or inside a guide exon.
  void merge_close_groups(BundleData *bundle, AssemblyState& state) {
    std::vector<GffObj *>& guides = bundle->guides;

    // Walk the groups of each strand in order and merge the close ones.
    if (min_bundle_distance || (count_of(guides))) {
      for (char strand : {'-', '.', '+'}) {
        CGroup *prev_group = NULL;
        CGroup *curr_group = state.start_group[strand];
        while (curr_group != NULL) {
          if (prev_group) {
            // We might consider separating the guides here and increasing the
            // edges for that case
            bool no_guide_boundary
                = !state.junction_start_positions.count(prev_group->end)
                  && !state.junction_end_positions.count(curr_group->start);
            bool gap_is_small
                = curr_group->start - prev_group->end <= min_bundle_distance;
            if (no_guide_boundary
                && (gap_is_small
                    || (count_of(guides)
                        && is_within_guide_exon(guides, strand,
                                                prev_group->end,
                                                curr_group->start)))) {
              merge_adjacent_groups(state, prev_group, curr_group);
              curr_group = prev_group->next_gr;
              continue;
            }
          }
          prev_group = curr_group;
          curr_group = curr_group->next_gr;
        }
      }
    }

    // Forget the guide-derived boundaries, which nothing needs from here on.
    if (count_of(guides)) {
      state.junction_start_positions.clear();
      state.junction_end_positions.clear();
    }
  }

} // namespace ribbon
