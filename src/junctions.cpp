#include "ribbon.h"
#include "GThreads.h"

// Anchor length for junction to be considered well supported
extern uint32_t junction_support;
extern uint32_t splice_site_error;
// Number of reads needed to support a particular junction
extern int32_t min_junction_coverage;

// Junction support: counts the reads that support each junction, rescues or
// redirects weak junctions and filters the ones that stay unsupported.

namespace ribbon {

  // Finishes the read statistics of one bundle before the graph is built: fills
  // bp_coverage from the reads (accumulate_read_coverage), turns it into the
  // per-block cumulative form, and sets the anchor support of each junction.
  void compute_junction_support(BundleData *bundle) {

    // Size the coverage arrays of the bundle.
    ReadList& reads = bundle->reads;
    BpCoverage& bp_coverage = bundle->bp_coverage;
    
    for (char strand : {'-', '.', '+'}) {
      bp_coverage[strand].resize(bundle->end - bundle->start + 3);
    }

    // For each read: add its coverage and the support of its junctions, from
    // the length of its exons on each side.
    for (int32_t read_idx = 0; read_idx < count_of(reads); read_idx++) {
      CReadAln& read = *(reads[read_idx]);
      float read_count = read.read_count;

      int32_t exon_count = count_of(read.segs);

      accumulate_read_coverage(reads, read_idx, bp_coverage, bundle);

      // Longest exon so far from each end of the read, for each junction.
      std::vector<uint32_t> left_support;
      std::vector<uint32_t> right_support;
      uint32_t max_left_support = 0;
      uint32_t max_right_support = 0;

      for (int32_t i = 0; i < exon_count; i++) {
        if (i) {
          if (read.segs[i - 1].len() > max_left_support)
            max_left_support = read.segs[i - 1].len();
          if (read.segs[exon_count - i].len() > max_right_support)
            max_right_support = read.segs[exon_count - i].len();
          left_support.push_back(max_left_support);
          right_support.push_back(max_right_support);
        }
      }

      // Add the read to the support of each junction on the sides where its
      // exon is at least as long as the anchor.
      for (int32_t i = 1; i < exon_count; i++) {
        uint32_t anchor = junction_support;
        CJunction& junc = *read.juncs[i - 1];
        bool is_long_intron = junc.len() > long_intron;
        bool is_weakly_supported
            = junc.n_reads - junc.nm < min_junction_coverage && !junc.mm;
        if (is_long_intron || is_weakly_supported) {
          if (anchor < long_intron_anchor) anchor = long_intron_anchor;
        }

        // Support only comes from spliced reads that are bigger than the anchor
        if (left_support[i - 1] >= anchor) {
          read.juncs[i - 1]->left_support += read_count;
          if (right_support[exon_count - i - 1] >= anchor) {
            read.juncs[i - 1]->right_support += read_count;
            read.juncs[i - 1]->n_reads_good += read_count;
          }
        } else if (right_support[exon_count - i - 1] >= anchor) {
          read.juncs[i - 1]->right_support += read_count;
        }
      }
    }

    // Turn the difference array into per-block running sums, so that interval
    // coverages can be found fast (get_region_coverage).
    int32_t base_count = count_of(bp_coverage['.']);
    int32_t last_full_block = int32_t((base_count - 1) / BSIZE);

    // This loop runs over every base, so it uses plain numbers instead of
    // strand characters: 0 is the positive strand, 1 the negative strand and 2
    // all reads. The three running sums are independent and advance together,
    // which keeps the loop fast.
    float *coverage[3] = {bp_coverage['+'].data(), bp_coverage['-'].data(),
                          bp_coverage['.'].data()};
    float prev_val[3] = {0, 0, 0};
    float prev_sum[3] = {0, 0, 0};
    int32_t block_idx = 0;
    while (block_idx < last_full_block) {
      int32_t end = (block_idx + 1) * BSIZE;
      for (int32_t i = block_idx * BSIZE; i < end; i++) {
        for (int32_t k = 0; k < 3; k++) {
          coverage[k][i] += prev_val[k];
          if (coverage[k][i] < 0) coverage[k][i] = 0;
          prev_val[k] = coverage[k][i];
          coverage[k][i] += prev_sum[k];
          prev_sum[k] = coverage[k][i];
        }
      }

      for (int32_t k = 0; k < 3; k++) prev_sum[k] = 0;

      block_idx++;
    }

    // Do the same for the last, partial block.
    for (int32_t i = block_idx * BSIZE; i < base_count; i++) {
      for (int32_t k = 0; k < 3; k++) {
        coverage[k][i] += prev_val[k];
        if (coverage[k][i] < 0) coverage[k][i] = 0;
        prev_val[k] = coverage[k][i];
        coverage[k][i] += prev_sum[k];
        prev_sum[k] = coverage[k][i];
      }
    }
  }

  // Tiered threshold for demoting a junction to a better supported neighbor:
  // the weaker the junction (fewer reads), the larger the fraction of the
  // neighbor's reads it may have and still be demoted.
  static inline double junc_ratio_threshold(double support) {
    if (support <= 25.0) return 1.0; // under 25
    else if (support <= 100.0) return 0.50; // 25-100
    else if (support <= 200.0) return 0.25; // 100-200 included here
    else return 0.10; // >200
  }

  // Should junction `cur` be redirected to the nearby junction `better`? Yes if
  // `better` has reads and cur's read count is at most
  // junc_ratio_threshold(cur's reads) times better's.
  static inline bool ok_to_demote(const CJunction *cur,
                                  const CJunction *better) {
    const double better_reads = (better->effective_n_reads() > 0.0
                                 ? better->effective_n_reads() : 0.0);
    // Can't demote to a neighbor with no support
    if (better_reads <= 0.0) return false;

    const double cur_reads = (cur->effective_n_reads() > 0.0
                                  ? cur->effective_n_reads() : 0.0);
    const double threshold = junc_ratio_threshold(cur_reads);
    return (cur_reads / better_reads) <= threshold;
  }

  // Looks at one side of junctions[i] (its start, or its end if end_side): if
  // it is a bad stranded junction, it may be rejected, and it is redirected to
  // a better supported junction nearby (start_redirect / end_redirect).
  static void redirect_to_better_junction(std::vector<CJunction *>& junctions,
                                          int32_t i, bool end_side,
                                          uint32_t search_window,
                                              float tolerance,
                                          int32_t ref_start,
                                          BpCoverage& bp_coverage) {
    CJunction *junc = junctions[i];

    // Accessors for the side we are looking at.
    auto pos = [&](const CJunction *x) -> uint32_t {
      return end_side ? x->end : x->start;
    };
    auto side_support = [&](const CJunction *x) {
      return end_side ? x->right_support : x->left_support;
    };
    auto redirect = [&](CJunction *x) -> int32_t& {
      return end_side ? x->end_redirect : x->start_redirect;
    };

    // Only stranded junctions that are bad (some reads with many mismatches,
    // no guide match, and mostly bad reads) are considered.
    if (junc->strand == '.') return;
    if (!junc->nm || junc->guide_match || junc->nm < junc->effective_n_reads())
      return;

    // Reject the junction if its good read support is too low, or if the
    // coverage inside the intron stays high next to the junction.
    if (junc->effective_n_reads_good() >= 0 && !junc->rejected) {
      // Threshold for bad junctions is higher. We do not reject short
      // junctions here: a 100 bp minimum removed real 85-90 bp introns on the
      // simulated cells.
      if (junc->effective_n_reads_good() < 1.25 * min_junction_coverage) {
        junc->rejected = true;
      } else { // check for coverage support of junction
        int32_t exon_point = pos(junc) - ref_start;
        int32_t intron_point = end_side ? exon_point - 1 : exon_point + 1;
        double exon_coverage = get_region_coverage('.', exon_point, exon_point,
                                                   bp_coverage);
        double intron_coverage = get_region_coverage('.', intron_point,
                                                     intron_point, bp_coverage);
        if (intron_coverage > tolerance * exon_coverage) junc->rejected = true;
      }
    }

    float support = 0; // support of the junction chosen so far to redirect to

    // Redirects junc to junctions[j] if it may be demoted to it.
    auto redirect_to = [&](int32_t j) {
      if (!ok_to_demote(junc, junctions[j])) return false;
      redirect(junc) = j;
      support = side_support(junctions[j]);
      return true;
    };
    // A guide junction, or one with mostly good reads.
    auto is_reliable = [&](const CJunction *x) {
      return x->guide_match || x->nm < x->effective_n_reads();
    };
    // An unreliable junction can still take over if it is close, better
    // supported than the junction chosen so far and than junc.
    auto is_better_alternative = [&](const CJunction *x, uint32_t distance) {
      return side_support(x) > support && distance < splice_site_error
          && side_support(x) * tolerance > side_support(junc);
    };

    // Look before the junction for a better supported one to redirect it to.
    int32_t j = i - 1;
    bool keep_searching = true;
    bool reliable = false;

    while (j > 0 && pos(junc) - pos(junctions[j]) < search_window) {
      CJunction *other = junctions[j];
      if (other->strand == junc->strand) {
        // Found a junction with the same coordinate -> we have already
        // searched it if it's bad
        if (pos(other) == pos(junc)) {
          if (redirect(other)
              && ok_to_demote(junc, junctions[redirect(other)])) {
            redirect(junc) = redirect(other);
            keep_searching = false;
          }
          break;
        // Nearby junction is much more reliable
        } else if (is_reliable(other)) {
          // The good junction is close enough
          if (side_support(other) > side_support(junc) * tolerance
              && redirect_to(j)) {
            reliable = true;
            break;
          }
        } else if (is_better_alternative(other, pos(junc) - pos(other))) {
          redirect_to(j);
        }
      }
      j--;
    }

    // Then look after it, unless it already follows a redirected junction.
    if (keep_searching) {
      j = i + 1;
      int32_t dist = search_window;
      if (redirect(junc)) {
        dist = pos(junc) - pos(junctions[redirect(junc)]);
      }

      while (j < count_of(junctions)
             && pos(junctions[j]) - pos(junc) < search_window) {
        CJunction *other = junctions[j];
        if (other->strand == junc->strand && pos(other) != pos(junc)) {
          int32_t d = (int32_t)(pos(other) - pos(junc));
          if (is_reliable(other)) {
            bool is_closer_or_stronger
                = d < dist || (d == dist && side_support(other) > support);
            bool is_much_stronger
                = side_support(other) > side_support(junc) * tolerance;
            if (is_closer_or_stronger && is_much_stronger && redirect_to(j)) {
              break;
            }
          // Junction is not best within window
          } else if (!reliable && is_better_alternative(other, (uint32_t)d)) {
            redirect_to(j);
          }
        }
        j++;
      }
    }
  }

  // True if the splice donor of an intron that starts after the exon base at
  // site (an index into the genomic sequence) is on this strand: GT or GC on
  // the positive strand, CT (the reverse complement of AG) on the negative one.
  bool has_donor_motif(const char *gseq, uint32_t site, char strand) {
    char first = tolower(gseq[site + 1]);
    char second = tolower(gseq[site + 2]);
    if (strand == '+') return first == 'g' && (second == 't' || second == 'c');
    return first == 'c' && second == 't';
  }

  // True if the splice acceptor of an intron that ends at the exon base at
  // site (an index into the genomic sequence) is on this strand: AG on the
  // positive strand, AC or GC (the reverse complements of GT and GC) on the
  // negative one.
  bool has_acceptor_motif(const char *gseq, uint32_t site, char strand) {
    char before = tolower(gseq[site - 2]);
    char last = tolower(gseq[site - 1]);
    if (strand == '+') return before == 'a' && last == 'g';
    return last == 'c' && (before == 'a' || before == 'g');
  }

  // Ends the run of junctions of list that share a start (or an end, if
  // by_end) and precede index i: weak junctions without the splice consensus
  // lose their strand, the others get the support summed over the run.
  template <class JunctionPtrList>
  void close_junction_run(JunctionPtrList& list, int32_t i, bool by_end,
                          uint32_t run_position, char consensus,
                          const PerStrand<double>& support) {
    for (int32_t j = i - 1; j >= 0; j--) {
      CJunction *junction = &*list[j];
      if ((by_end ? junction->end : junction->start) != run_position) break;
      if (!junction->guide_match && !consensus
          && junction->n_reads_good < DROP / ERROR_PERC) junction->strand = '.';
      if (junction->strand != '.') {
        (by_end ? junction->right_support : junction->left_support)
            = support[junction->strand];
      }
    }
  }

  // Adds up the support of junctions that share a start or an end, sets the
  // splice-site consensus, and redirects each unreliable junction to a better
  // supported nearby one (recorded in its start_redirect or end_redirect).
  void propagate_junction_support(BundleData *bundle,
                                  AssemblyState& state) {
    int32_t ref_start = bundle->start;
    JunctionList& junction = bundle->junction;
    BpCoverage& bp_coverage = bundle->bp_coverage;

    sort_by(state.end_sorted_junctions, compare_junction_by_end);

    // Walk the junctions in start order and, in parallel, in end order: when a
    // new start (end) begins, the previous run gets its summed left (right)
    // support and we set the splice-site consensus.
    uint32_t run_start = 0;
    uint32_t run_end = 0;
    PerStrand<double> left_support; // the sums of the junctions of each strand
    PerStrand<double> right_support;
    bool has_bad_junctions = false;
    char left_consensus = -1;
    char right_consensus = -1;
    for (int32_t i = 0; i < count_of(junction); i++) {
      // Remember whether some stranded junction is only supported by bad reads.
      if ((!has_bad_junctions) && junction[i]->strand != '.'
          && junction[i]->nm == junction[i]->n_reads
          && !junction[i]->guide_match) {
        has_bad_junctions = true;
      }

      if (junction[i]->start != run_start) { // new junction starting at i

        // Close the previous run of junctions with the same start: they get
        // the summed support; weak non-consensus ones lose their strand.
        close_junction_run(junction, i, false, run_start, left_consensus,
                           left_support);

        // Unstrand the weaker of two identical junctions on different strands.
        int32_t j = i + 1; // check if there is the same junction further ahead
        if (junction[i]->strand != '.')
          while (j < count_of(junction)
                 && junction[j]->start == junction[i]->start
                 && junction[j]->end == junction[i]->end) {
            if (junction[j]->strand != '.'
                && junction[i]->strand != junction[j]->strand) {
              // Possible missaligned junction --> delete junction strand
              if (junction[i]->n_reads > junction[j]->n_reads
                  && junction[i]->n_reads_good
                  > junction[j]->n_reads_good) junction[j]->strand = '.';
              if (junction[i]->n_reads < junction[j]->n_reads
                  && junction[i]->n_reads_good
                  < junction[j]->n_reads_good) junction[i]->strand = '.';
            }
            j++;
          }

        // Start the sums of the new run and set its splice-site consensus.
        left_support['-'] = 0;
        left_support['+'] = 0;
        // We might have deleted the i junction above
        if (junction[i]->strand != '.')
          left_support[junction[i]->strand]
              = junction[i]->left_support;
        run_start = junction[i]->start;

        if (bundle->gseq && junction[i]->strand != '.') {
          left_consensus = has_donor_motif(bundle->gseq,
              junction[i]->start - ref_start, junction[i]->strand);
        }
      } else if (junction[i]->strand != '.')
        left_support[junction[i]->strand]
            += junction[i]->left_support;

      // Do the same on the junction ends: close the previous run with the same
      // end and start a new one.
      CJunction *end_junction = state.end_sorted_junctions[i];
      if (end_junction->end != run_end) {
        close_junction_run(state.end_sorted_junctions, i, true, run_end,
                           right_consensus, right_support);
        // We do not check here for possible missalignments
        right_support['-'] = 0;
        right_support['+'] = 0;
        if (end_junction->strand != '.')
          right_support[end_junction->strand]
              = end_junction->right_support;
        run_end = end_junction->end;

        if (bundle->gseq && end_junction->strand != '.') {
          right_consensus = has_acceptor_motif(bundle->gseq,
              end_junction->end - ref_start, end_junction->strand);
        }
      } else if (end_junction->strand != '.') {
        right_support[end_junction->strand]
            += end_junction->right_support;
      }
    }

    // Some reads have very bad junctions: look for a better junction nearby
    // to redirect them to.
    if (has_bad_junctions) {
      uint32_t search_window = junction_support;

      float tolerance = 1 - ERROR_PERC;

      // Check each bad junction on the start side, then on the end side.
      std::vector<CJunction *> start_sorted_junctions;
      for (auto& j : junction) start_sorted_junctions.push_back(j.get());
      for (int32_t i = 1; i < count_of(junction); i++) {
        redirect_to_better_junction(start_sorted_junctions, i, false,
                                    search_window, tolerance, ref_start,
                                    bp_coverage);
        redirect_to_better_junction(state.end_sorted_junctions, i, true,
                                    search_window, tolerance, ref_start,
                                    bp_coverage);
      }
    } // if higherr
  }

  // Looks for the junction with this start, end and strand among the ones that
  // were added to the list after its first original_count junctions. Returns
  // NULL if there is none.
  CJunction *find_added_junction(JunctionList& junction,
                                 int32_t original_count, uint32_t start,
                                 uint32_t end, char strand) {
    for (int32_t k = original_count; k < count_of(junction); k++) {
      if (junction[k]->strand == strand && junction[k]->start == start
          && junction[k]->end == end) {
        return junction[k].get();
      }
    }
    return NULL;
  }

  // Returns whether a junction is well supported rather than an artifact: guide
  // junctions are always kept; others are rejected for low read support, a long
  // intron of bad reads, no coverage step or a small share of the coverage.
  bool is_junction_well_supported(CJunction& junction, int32_t ref_start,
                                  BpCoverage& bp_coverage) {
    // Guide junctions are always kept (a read that calls this function covers
    // them).
    if (junction.guide_match) return true;
    // ^ KEEP THIS IN MIND IF WE CHANGE HOW WE USE GOOD_JUNC

    // Reject junctions with less good read support than the minimum.
    if (junction.effective_n_reads_good() < min_junction_coverage) {
      junction.strand = '.';
      return false;
    }

    bool mismatch = false;

    // Flag a junction that only reads with many mismatches support.
    if (junction.nm
        && round(junction.nm) == round(junction.effective_n_reads())) {
      mismatch = true;
      // Don't believe long intron if all junctions are from bad reads unless
      // highly covered
      if (junction.end - junction.start > long_intron
          && junction.effective_n_reads() < CHI_WIN * ERROR_PERC) {
        junction.strand = '.';
        return false;
      }
    }

    // Reject a junction if the coverage next to it hardly changes and its
    // support is a small fraction of that coverage. We sum the coverage in a
    // window of window_len bp on each side of the junction boundaries.
    float support_boost = 1 / ERROR_PERC;
    int32_t window_len = 5;
    int32_t pos = junction.start - ref_start;

    double start_left_coverage = 0;
    double start_right_coverage = 0;
    if (pos + window_len + 1 < count_of(bp_coverage['.'])
        && pos >= window_len) {
      start_left_coverage = get_region_coverage('.', pos - window_len + 1, pos,
                                                bp_coverage);
      start_right_coverage = get_region_coverage('.', pos + 1, pos + window_len,
                                                 bp_coverage);
    }

    // Shouldn't coverage_left be normalized by window_len below when comparing
    // to left_support?
    // gave some boost (*10) to junctions here assumming not all of them are
    // captured
    bool start_support_is_tiny = start_left_coverage > 1 / ERROR_PERC
        && junction.left_support * support_boost
               < ERROR_PERC * start_left_coverage;
    bool start_coverage_hardly_changes = mismatch
        || start_right_coverage > start_left_coverage * (1 - ERROR_PERC);
    if (start_support_is_tiny && start_coverage_hardly_changes) {
      junction.strand = '.';
      return false;
    }

    pos = junction.end - ref_start - 1;

    // The same check at the end of the junction.
    double end_left_coverage = 0;
    double end_right_coverage = 0;
    if (pos + window_len + 1 < count_of(bp_coverage['.'])
        && pos >= window_len) {
      end_left_coverage = get_region_coverage('.', pos - window_len + 1, pos,
                                              bp_coverage);
      end_right_coverage = get_region_coverage('.', pos + 1, pos + window_len,
                                               bp_coverage);
    }

    // Gave some boost (*10) to junctions here assumming not all of them are
    // captured
    bool end_support_is_tiny = end_right_coverage > 1 / ERROR_PERC
        && junction.right_support * support_boost
               < ERROR_PERC * end_right_coverage;
    bool end_coverage_hardly_changes = mismatch
        || end_left_coverage > end_right_coverage * (1 - ERROR_PERC);
    if (end_support_is_tiny && end_coverage_hardly_changes) {
      junction.strand = '.';
      return false;
    }

    // Reject a junction whose support is a tiny share of all junctions at its
    // start or end (probably a spliceosome error).
    if (junction.effective_n_reads() * 10 < ERROR_PERC * junction.left_support
        || junction.effective_n_reads() * 10
        < ERROR_PERC * junction.right_support) {
      junction.strand = '.';
      return false;
    }

    return true;
  }

  // Comparator: orders junctions by end position, then by start.
  int32_t compare_junction_by_end(const void *p1, const void *p2) {
    CJunction *a = (CJunction *)p1;
    CJunction *b = (CJunction *)p2;
    if (a->end < b->end) return -1;
    if (a->end > b->end) return 1;
    if (a->start < b->start) return -1;
    if (a->start > b->start) return 1;
    return 0;
  }

} // namespace ribbon
