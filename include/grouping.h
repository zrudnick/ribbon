#pragma once

#include "ribbon.h"

// Places reads into groups (same-strand runs of overlapping exon coverage)
// and links the groups of read pairs.

namespace ribbon {
  // Places a read (and its mate, if they form one fragment) into the group
  // chain of one strand. Each exon segment joins an overlapping group
  // (merging the groups it bridges) or starts a new CGroup. If the strand
  // has no groups yet, the read starts the chain.
  class ReadSegmentPlacer {
    ReadList& reads;
    AssemblyState& state;
    int32_t pair_slot; // -1 if the read is placed alone
    float read_coverage;
    char strand;
    PerStrand<int32_t>& last_link_id;
    uint32_t merge_distance;

    int32_t read_idx; // the read being placed (the mate once we move on to it)
    int32_t mate_idx; // -1 if the read is placed alone
    int32_t read_link_id;
    bool first_visit = true; // false once the read is seen in a group twice
    // The group where the chain of this strand is left for the next read
    CGroup *curr_group;

    // The walk over the exons of the read
    CGroup *prev_group = NULL;
    CGroup *matching_group = NULL; // the group where the exon fits
    int32_t exon_count = 0;
    int32_t exon_idx = 0;
    int32_t last_pushed_group_id = -1;

    void find_group_where_read_starts();
    void check_pair_at_first_exon();

    void unify_link_ids();
    void add_exon_to_matching_group();

    void add_exon_as_new_group();
    void continue_with_mate_exons();

    void place_in_existing_chain();
    void place_in_new_chain();

  public:
    ReadSegmentPlacer(ReadList& reads_in, AssemblyState& state_in,
                      ReadPlacement placement,
                      PerStrand<int32_t>& last_link_id_in);

    void place();
  };

  void assign_read_segments_to_group(ReadList& reads, AssemblyState& state,
                                     ReadPlacement placement,
                                     PerStrand<int32_t>& last_link_id);

  void assign_read_to_groups(int32_t read_idx, ReadList& reads,
                             AssemblyState& state);

  void assign_reads_to_groups(BundleData *bundle, AssemblyState& state);

  void merge_close_groups(BundleData *bundle, AssemblyState& state);
} // namespace ribbon
