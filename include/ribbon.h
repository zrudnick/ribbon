#pragma once

#include <algorithm>
#include <filesystem>
#include <functional>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "GArgs.h"
#include "GList.hh"
#include "GSam.h"
#include "gff.h"
#include "time.h"
#include "utils.h" // helpers that do not depend on ribbon types

// A dynamically sized bit vector: the path_bits of nodes and edges that
// transfrags, paths and graph nodes carry. Missing bits count as 0 when
// comparing or combining vectors of different sizes.
class BitVec {
 public:
  class BitRef {
   public:
    BitRef(uint64_t& word, uint32_t bit) : word_(word), bit_(bit) {}
    BitRef& operator=(bool value) {
      if (value) word_ |= (uint64_t(1) << bit_);
      else word_ &= ~(uint64_t(1) << bit_);
      return *this;
    }
    BitRef& operator=(const BitRef& other) { return *this = bool(other); }
    operator bool() const { return (word_ >> bit_) & 1; }

   private:
    uint64_t& word_;
    uint32_t bit_;
  };

  BitVec() : size_(0) {}
  explicit BitVec(size_t size, bool value = false) : size_(0) {
    resize(size, value);
  }

  size_t size() const { return size_; }

  bool operator[](size_t i) const { return (words_[i / 64] >> (i % 64)) & 1; }
  BitRef operator[](size_t i) { return BitRef(words_[i / 64], i % 64); }

  // Removes all bits (size() becomes 0).
  void clear() {
    words_.clear();
    size_ = 0;
  }

  // Sets all bits to 0 (size() stays the same).
  void reset() { std::fill(words_.begin(), words_.end(), uint64_t(0)); }

  // Grows or shrinks to the given number of bits; new bits get the value.
  void resize(size_t size, bool value = false) {
    if (size > size_ && value && size_ % 64)
      words_[size_ / 64] |= ~uint64_t(0) << (size_ % 64);
    words_.resize((size + 63) / 64, value ? ~uint64_t(0) : uint64_t(0));
    size_ = size;
    if (size % 64) words_.back() &= (uint64_t(1) << (size % 64)) - 1;
  }

  // Number of bits that are set.
  size_t count() const {
    size_t n = 0;
    for (uint64_t word : words_) n += __builtin_popcountll(word);
    return n;
  }

  bool operator==(const BitVec& other) const {
    size_t common = std::min(words_.size(), other.words_.size());
    for (size_t i = 0; i < common; i++)
      if (words_[i] != other.words_[i]) return false;
    for (size_t i = common; i < words_.size(); i++)
      if (words_[i]) return false;
    for (size_t i = common; i < other.words_.size(); i++)
      if (other.words_[i]) return false;
    return true;
  }
  bool operator!=(const BitVec& other) const { return !(*this == other); }

  // True if every bit that is set in other is set here too. Same as
  // (*this & other) == other, without making a copy.
  bool contains(const BitVec& other) const {
    size_t common = std::min(words_.size(), other.words_.size());
    for (size_t i = 0; i < common; i++)
      if (other.words_[i] & ~words_[i]) return false;
    for (size_t i = common; i < other.words_.size(); i++)
      if (other.words_[i]) return false;
    return true;
  }

  // Bits only in this vector become 0; bits only in other are ignored.
  BitVec& operator&=(const BitVec& other) {
    size_t common = std::min(words_.size(), other.words_.size());
    for (size_t i = 0; i < common; i++) words_[i] &= other.words_[i];
    for (size_t i = common; i < words_.size(); i++) words_[i] = 0;
    return *this;
  }

  // Grows to the size of other first if it is longer.
  BitVec& operator|=(const BitVec& other) {
    if (size_ < other.size_) resize(other.size_);
    for (size_t i = 0; i < other.words_.size(); i++)
      words_[i] |= other.words_[i];
    return *this;
  }

 private:
  std::vector<uint64_t> words_;
  size_t size_;
};

inline BitVec operator&(BitVec a, const BitVec& b) { return a &= b; }
inline BitVec operator|(BitVec a, const BitVec& b) { return a |= b; }

#define MIN_EXON_OVERLAP 5 // minimum overlap for exons

namespace ribbon {
  // Read count data for exon or intron of a reference transcript
  struct FeatureCounts {
    // Feature id (>0), +1 to the index in
    // BundleCounts::local_exon_counts/local_intron_counts
    uint32_t id;
    std::vector<uint32_t> transcript_ids; // owners (guides index + 1)
    long start = 0;
    long end = 0; // genomic coordinates for the feature
    char strand = 0;
    // Number of read alignments that overlap this feature (more than 5 bp for
    // exons, exact coordinate match for introns).
    mutable uint32_t read_count = 0;

    FeatureCounts(long start0 = 0, long end0 = 0, char s = '.', uint32_t fid
                  = 0,
                  uint32_t tid = 0)
        : id(fid), transcript_ids(1), start(start0), end(end0), strand(s) {
      if (start > end) std::swap(start, end);
      if (tid > 0) transcript_ids.push_back(tid);
    }
    FeatureCounts(const FeatureCounts& seg)
        : id(seg.id), transcript_ids(seg.transcript_ids), start(seg.start),
          end(seg.end), strand(seg.strand) {
    }

    FeatureCounts(const FeatureCounts& seg, uint32_t tid)
        : id(seg.id), transcript_ids(1), start(seg.start), end(seg.end),
          strand(seg.strand) {
      if (start > end) std::swap(start, end);
      if (tid > 0) transcript_ids.push_back(tid);
    }

    // Orders by start, then end, then strand (a '.' strand matches any).
    bool operator<(const FeatureCounts& o) const {
      if (start != o.start) return (start < o.start);
      if (end != o.end) return (end < o.end);
      if (strand == '.' || o.strand == '.') return false;
      if (strand != o.strand) return (strand < o.strand);
      return false;
    }

    // Equal if start and end match and the strands match or one is '.'.
    bool operator==(const FeatureCounts& o) const {
      bool strands_match = strand == o.strand || strand == '.' || o.strand == '.';
      return start == o.start && end == o.end && strands_match;
    }
  };

  // The guide exons (or introns) of one bundle, in the order they were first
  // seen: feature id N is at index N-1. The table owns them.
  using FeatureTable = std::vector<std::unique_ptr<FeatureCounts>>;

  struct ExonOverlap {
    // Pointer to an item of BundleCounts::local_exon_counts
    FeatureCounts *feature = NULL;
    byte mate_overlap = 0; // =1 if the mate of this read overlaps the same exon
    long overlap_len = 0;

    // Orders by mate overlap, then longer overlap, then longer feature, then
    // strand and start.
    bool operator<(const ExonOverlap& overlap) const {
      if (mate_overlap != overlap.mate_overlap)
        return (mate_overlap > overlap.mate_overlap);
      if (overlap_len != overlap.overlap_len)
        return (overlap_len > overlap.overlap_len);
      if (feature->end - feature->start != overlap.feature->end
          - overlap.feature->start)
        return (feature->end - feature->start > overlap.feature->end
                - overlap.feature->start);
      if (feature->strand != overlap.feature->strand)
        return (feature->strand < overlap.feature->strand);
      return (feature->start < overlap.feature->start);
    } // operator <

    // Equal if the mate flag, overlap length and feature are all the same.
    bool operator==(const ExonOverlap& overlap) const {
      return (mate_overlap == overlap.mate_overlap
              && overlap_len == overlap.overlap_len
              && feature == overlap.feature);
    }
    ExonOverlap(FeatureCounts *f = NULL, long overlap_length = 0,
                byte m_overlap = 0)
        : feature(f), mate_overlap(m_overlap), overlap_len(overlap_length) {
    }
  };

  struct TranscriptCounts { // storing read-count data for a transcript (guide)
    // Read-count data of a reference transcript (guide). We create one for each
    // guide, even those not used in a bundle; it tracks the exon and intron
    // features used during assembly (see set_guide_status/get_guide_status).
    GffObj *reference_transcript;
    uint32_t transcript_id; // index in the global transcript-counts table
    long start;
    long end;
    std::vector<FeatureCounts *> transcript_exons;
    std::vector<FeatureCounts *> transcript_introns;
    // Scans the sorted feature_set from idx, to the left if search_left or to
    // the right otherwise. Returns false and sets found_id if an equal feature
    // exists, else true with idx at the place where new_feature belongs.
    static bool locate_feature(std::vector<FeatureCounts *>& feature_set,
                               FeatureCounts *new_feature, bool search_left,
                               int32_t& idx, int32_t& found_id) {
      bool eq = false;
      if (search_left) {
        while (*new_feature < *(feature_set[idx])
               || (eq = (*new_feature == *(feature_set[idx])))) {
          if (eq) {
            found_id = feature_set[idx]->id;
            return false;
          }
          // The new_feature < feature_set[idx]
          --idx;
          if (idx < 0) break; // new_feature should be inserted at 0
        } // while new_feature<feature_set[idx]
        ++idx;
        // The new_feature < feature_set[idx+1]
        // we'll insert new_feature at position idx+1
      } else { // new_feature >= *feature_set[idx]
        while (*(feature_set[idx]) < *new_feature
               || (eq = (*new_feature == *(feature_set[idx])))) {
          if (eq) {
            found_id = feature_set[idx]->id;
            return false;
          }
          ++idx;
          if (idx == count_of(feature_set)) {
            // The new_feature should be appended to the list
            break;
          }
        }
      }
      return true;
    }

    // Registers the exons and introns of the reference transcript as features.
    void collect_features(uint32_t& next_exon_id,
                          std::vector<FeatureCounts *>& exon_set,
                          FeatureTable& exon_table,
                          uint32_t& next_intron_id,
                          std::vector<FeatureCounts *>& intron_set,
                          FeatureTable& intron_table) {
      GffObj& transcript = *(reference_transcript);

      // Add each exon, and the intron before it, to the bundle's feature sets.
      int32_t exon_cache_idx = count_of(exon_set) - 1;
      int32_t intron_cache_idx = count_of(intron_set) - 1;
      for (int32_t i = 0; i < transcript.exons.Count(); ++i) {
        add_feature((long)transcript.exons[i]->start,
                    (long)transcript.exons[i]->end, transcript_exons,
                    next_exon_id, exon_set, exon_table, exon_cache_idx);
        if (i > 0) { // store intron
          add_feature(long(transcript.exons[i - 1]->end) + 1,
                      long(transcript.exons[i]->start) - 1, transcript_introns,
                      next_intron_id, intron_set, intron_table,
                      intron_cache_idx);
        } // for each intron
      } // for each exon
    }

    // Registers one exon or intron of this transcript in the bundle's feature
    // set and table. An identical feature gains this transcript as a parent;
    // otherwise we add a new feature. Either way it joins transcript_features.
    void add_feature(long feature_start, long feature_end,
                     std::vector<FeatureCounts *>& transcript_features,
                     uint32_t& next_id,
                     std::vector<FeatureCounts *>& feature_set,
                     FeatureTable& feature_table, int32_t& cache_idx) {
      // Create the feature and look for an equal one in the sorted set,
      // starting from the cached index (next_id is the largest id so far).
      bool add_new = true;
      // The new feature belongs to feature_table if it is not a duplicate.
      std::unique_ptr<FeatureCounts> new_feature_owner
          = std::make_unique<FeatureCounts>(feature_start, feature_end,
              reference_transcript->strand, 0, this->transcript_id);
      FeatureCounts *new_feature = new_feature_owner.get();
      int32_t idx = cache_idx < 0 ? count_of(feature_set) - 1 : cache_idx;
      int32_t found_id = -1;
      if (count_of(feature_set) > 0) {
        add_new = locate_feature(feature_set, new_feature,
                                 *new_feature < *(feature_set[idx]), idx,
                                 found_id);
      } // check existing set

      // Insert a new feature, or add this transcript to an equal feature.
      if (add_new) { // did not see this feature before
        new_feature->id = ++next_id;
        if (idx < 0) {
          feature_set.push_back(new_feature);
          idx = count_of(feature_set) - 1;
        } else feature_set.insert(feature_set.begin() + idx, new_feature);
        cache_idx = idx;
        feature_table.push_back(std::move(new_feature_owner));
        found_id = count_of(feature_table);

  #ifdef DEBUG
        if (count_of(feature_table) != (int32_t)next_id) {
          fprintf(stderr, "Error: feature_table.Count=%d, next_id=%d\n",
                  count_of(feature_table), next_id);
        }
  #endif
      } else { // feature seen before, only update its parent list
        feature_table[found_id
                      - 1]->transcript_ids.push_back(this->transcript_id);
      }

      // Add the resolved feature to the transcript's exons or introns.
      transcript_features.push_back(feature_table[found_id - 1].get());
    }
    TranscriptCounts(GffObj& transcript, uint32_t id = 0)
        : reference_transcript(&transcript), transcript_id(id),
          start(transcript.start), end(transcript.end), transcript_exons(),
          transcript_introns() {
    }
  };

  struct BundleCounts {
    long bundle_start;
    long cluster_end;
    // Count data of the reference transcripts of this bundle (element index =
    // TranscriptCounts::transcript_id - 1).
    std::vector<TranscriptCounts *> transcripts;
    // Guide exons of this bundle, sorted by start coordinate.
    std::vector<FeatureCounts *> guide_exons;
    // Set of guide introns in this bundle, sorted by start coordinate
    std::vector<FeatureCounts *> guide_introns;
    // Index of the first exon that overlaps exon_cache_pos.
    int32_t exon_cache_idx;
    // Left coordinate of last cached exon overlap query
    // (find_overlapping_exons())
    long exon_cache_pos;

    // Per-bundle tables of the guide features; the element index of a feature
    // is FeatureCounts::id - 1.
    FeatureTable local_exon_counts; // local exon read-count data
    FeatureTable local_intron_counts; // local intron read-count data
    // Local exon/intron ids within the bundle
    uint32_t exon_id_counter;
    uint32_t intron_id_counter;

    BundleCounts(long t_l = 0, long t_r = 0)
        : bundle_start(t_l), cluster_end(t_r), transcripts(), guide_exons(),
          guide_introns(), exon_cache_idx(0), exon_cache_pos(0),
          exon_id_counter(0), intron_id_counter(0) {}

    // Adds a guide to the bundle: extends the bundle range and collects the
    // exons and introns of the guide. Returns the guide's transcript_id.
    uint32_t add_transcript(GffObj& t) {
      // Extend the bundle range to cover the guide.
      if (bundle_start == 0 || bundle_start > (long)t.start) {
        bundle_start = t.start;
      }
      if (cluster_end == 0 || cluster_end < (long)t.end) {
        cluster_end = t.end;
      }

      // Collect the exons and introns of the guide as features.
      // We should always have a TranscriptCounts for each guide
      TranscriptCounts *tdata = (TranscriptCounts *)(t.uptr);
      transcripts.push_back(tdata);
      tdata->collect_features(exon_id_counter, guide_exons, local_exon_counts,
                              intron_id_counter, guide_introns,
                              local_intron_counts);
      return tdata->transcript_id;
    }

    // Finds the guide exons that overlap the query interval and adds them to
    // ex_overlaps (the strand and mate position refine the match). Returns
    // whether any exon overlaps.
    bool find_overlapping_exons(std::vector<ExonOverlap>& ex_overlaps,
                                long query_start, long query_end,
                                char strand = '.', int32_t mate_pos = 0) {
      // We start at the cached exon index if the last query was to the left of
      // this one.
      bool has_overlaps = false;
      if (count_of(guide_exons) == 0) return false;
      int32_t start_index = 0;
      bool no_cache = (exon_cache_pos == 0 || exon_cache_pos > query_start);
      if (no_cache) {
        exon_cache_idx = count_of(guide_exons) - 1;
        exon_cache_pos = 0;
      } else start_index = exon_cache_idx; // must have a valid value
      bool cache_update_pending = true;
      int32_t last_checked_exon = count_of(guide_exons) - 1;
      int32_t query_len = query_end - query_start + 1;
      int32_t min_overlap = GMIN(query_len, MIN_EXON_OVERLAP);

      // Check the exons that start before the end of the query.
      for (long i = start_index; i < count_of(guide_exons); ++i) {
        last_checked_exon = i;
        long exon_start = guide_exons[i]->start;
        long exon_end = guide_exons[i]->end;
        if (exon_start > query_end) break;
        if (query_start > exon_end) continue;
        // Exon overlap here
        long overlap_len = 0;
        if (query_start < exon_start) {
          overlap_len = (query_end < exon_end ? query_end - exon_start + 1
                         : exon_end - exon_start + 1);
        } else { // exon_start<=query_start
          overlap_len = (query_end < exon_end ? query_end - query_start + 1
                         : exon_end - query_start + 1);
        }
        if (cache_update_pending) {
          // Cache first overlap
          exon_cache_idx = i;
          cache_update_pending = false;
        }
        // Non-matching strand
        if (strand != '.' && strand != guide_exons[i]->strand) continue;

        // Keep the overlap if the mate read also overlaps this exon or the
        // overlap is long enough.
        byte mate_overlap = 0;
        if (mate_pos && mate_pos + 10 > exon_start && mate_pos + 5 < exon_end)
          mate_overlap = 1; // mate read likely overlaps this exon
        if (mate_overlap || overlap_len >= min_overlap) {
          has_overlaps = true;
          ExonOverlap overlap(guide_exons[i], overlap_len, mate_overlap);
          if (find_sorted(ex_overlaps, overlap) < 0)
            insert_sorted(ex_overlaps, overlap);
        }
      }

      // Remember where this query ended for the next one.
      // If no overlap was found, the cache points at the last exon checked.
      if (cache_update_pending) exon_cache_idx = last_checked_exon;
      exon_cache_pos = query_start;
      return has_overlaps;
    }

    // Returns the guide intron with exactly these coordinates and strand, or
    // NULL if there is none.
    FeatureCounts *find_intron(long query_start, long query_end, char strand) {
      FeatureCounts query(query_start, query_end, strand);
      auto it = std::lower_bound(guide_introns.begin(), guide_introns.end(),
                                 &query,
                                 [](const FeatureCounts *a,
                                    const FeatureCounts *b) {
                                   return *a < *b;
                                 });
      if (it != guide_introns.end() && !(query < **it)) return *it;
      return NULL;
    }
  }; // struct BundleCounts

  // Status of a guide within a bundle, stored in the guide's user flags.

  enum GuideBundleStatus {
    GBST_UNSET = 0,
    GBST_IN_BUNDLE, // 1: added to bundle
    GBST_ALL_INTR_COV, // 2: all introns covered by at least one read
    GBST_STORED, // 3: stored to be printed
  };
  GuideBundleStatus get_guide_status(GffObj *t);
  void set_guide_status(GffObj *t, GuideBundleStatus status);
  bool is_guide_supported_by_reads(GffObj *guide);

  // Collect all refguide transcripts for a single genomic sequence
  struct GRefData {
    // All guides on this genomic seq; they are added in order as sorted in
    // GffReader
    GList<GffObj> rnas;
    int32_t gseq_id; // chromosome index in GffObj::names->gseqs
    const char *gseq_name; // chromosome name
    // Rnas list: unsorted, dealloc on free, duplicates allowed
    GRefData(int32_t gid = -1)
        : rnas(false, true, false), gseq_id(gid), gseq_name(NULL) {
      gseq_id = gid;
      if (gseq_id >= 0)
        gseq_name = GffObj::names->gseqs.getName(gseq_id);
    }

    // Adds a guide to this chromosome's list (the first one sets the sequence).
    void add(GffReader *gffr, GffObj *t) {
      if (gseq_id >= 0) {
        if (gseq_id != t->gseq_id) {
          fprintf(stderr, "Error: invalid call to GRefData::add() - "
                          "different genomic sequence!\n");
          exit(1);
        }
      } else { // adding first transcript, initialize storage
        gseq_id = t->gseq_id;
        gseq_name = t->getGSeqName();
        if (gffr->gseqtable[gseq_id] == NULL) {
          fprintf(stderr, "Error: invalid genomic sequence data (%s)!\n",
                  gseq_name);
          exit(1);
        }
        rnas.setCapacity(gffr->gseqtable[gseq_id]->fcount);
      }
      rnas.Add(t);
      // Mark as used, to prevent deletion in GffReader destructor by
      // gflst::freeUnused()
      t->isUsed(true);
    }

    bool operator==(GRefData& d) {
      return gseq_id == d.gseq_id;
    }
    bool operator<(GRefData& d) {
      return (gseq_id < d.gseq_id);
    }
  };

  // The position of a strand in the per-strand arrays: 0 for '-', 1 for '.'
  // (unknown strand) and 2 for '+'.
  // The low two bits of '-', '.' and '+' are 1, 2 and 3, so no branch is needed.
  inline int32_t strand_index(char strand) {
    return (strand & 3) - 1;
  }

  // The other stranded strand: '+' for '-' and '-' for '+' ('.' stays '.').
  inline char opposite_strand(char strand) {
    return strand == '+' ? '-' : strand == '-' ? '+' : '.';
  }

  // One value for each strand: '-', '.' (unknown strand) and '+'. Only a char
  // can index it, so an old numeric strand index fails to compile.
  template <class T>
  class PerStrand {
    T values[3] = {};

  public:
    PerStrand() = default;
    explicit PerStrand(const T& initial) {
      for (T& value : values) value = initial;
    }
    T& operator[](char strand) { return values[strand_index(strand)]; }
    const T& operator[](char strand) const {
      return values[strand_index(strand)];
    }
    template <class I> T& operator[](I) = delete;
    template <class I> const T& operator[](I) const = delete;
  };

  // The coverage of each base of a bundle: bp_coverage['-'] and
  // bp_coverage['+'] for the reads of that strand, and bp_coverage['.'] for
  // all reads of any strand.
  using BpCoverage = PerStrand<std::vector<float>>;

  // A splice junction (intron) with its read support.
  struct CJunction : public GSeg {
    char strand; // '-', '+' or '.' (unknown)
    char guide_match; // exact match of a ref intron?
    double n_reads;
    double n_reads_good;
    // Set by propagate_junction_support when a better supported nearby
    // junction takes over this one: the index of that junction in the
    // start-sorted (start_redirect) or end-sorted (end_redirect) junction
    // list, or 0 if none. n_reads and n_reads_good keep their true counts.
    int32_t start_redirect;
    int32_t end_redirect;
    double left_support;
    double right_support;
    double nm; // number of reads with a high nm (high mismatch)
    // Number of reads that support a junction with both anchors bigger than
    // long_intron_anchor
    double mm;
    // Set when the junction is too weakly supported to be kept
    bool rejected;
    CJunction(int32_t s = 0, int32_t e = 0, char _strand = '.')
        : GSeg(s, e), strand(_strand), guide_match(0), n_reads(0),
          n_reads_good(0), start_redirect(0), end_redirect(0), left_support(0),
          right_support(0), nm(0), mm(0), rejected(false) {}

    // The read counts seen by the junction filter: a redirected junction
    // reports its redirect index as a negative count.
    double effective_n_reads() const {
      return start_redirect ? -start_redirect : n_reads;
    }
    double effective_n_reads_good() const {
      return end_redirect ? -end_redirect : n_reads_good;
    }

    // Orders by start, then end, then strand (reverse strand first).
    bool operator<(const CJunction& b) const {
      if (start < b.start) return true;
      if (start > b.start) return false;
      if (end < b.end) return true;
      if (end > b.end) return false;
      if (strand_index(strand) > strand_index(b.strand)) return true;
      return false;
    }

    // Equal if start, end and strand are the same.
    bool operator==(CJunction& b) {
      return (start == b.start && end == b.end && strand == b.strand);
    }
  };

  // The distinct junctions of one bundle, sorted by start, end and strand; the
  // list owns them (reads and end_sorted_junctions point to them).
  using JunctionList = std::vector<std::unique_ptr<CJunction>>;

  // Returns the junction in the sorted list that has the same start, end and
  // strand as new_junction. If there is none, moves new_junction into the list
  // at its place and returns it.
  inline CJunction *add_junction_if_new(
      JunctionList& junctions, std::unique_ptr<CJunction> new_junction) {
    auto it = std::lower_bound(
        junctions.begin(), junctions.end(), new_junction.get(),
        [](const std::unique_ptr<CJunction>& a, const CJunction *b) {
          return *a < *b;
        });
    if (it != junctions.end() && !(*new_junction < **it)) {
      return it->get();
    }
    CJunction *added = new_junction.get();
    junctions.insert(it, std::move(new_junction));
    return added;
  }

  // Compare function for sorting junctions by start, end and strand.
  inline int32_t compare_junction_by_position(const void *p1,
                                              const void *p2) {
    const CJunction& a = *(const CJunction *)p1;
    const CJunction& b = *(const CJunction *)p2;
    if (b < a) return 1;
    if (a < b) return -1;
    return 0;
  }
  struct CJunction;

  // A read alignment: its exon segments, junctions and read counts.
  struct CReadAln : public GSeg {
    char strand; // '-', '+' or '.' (unknown)
    int16_t nh;
    int32_t bc; // SCELL -1 if uninitiated
    std::string ub; // SCELL
    // Read 1 of a pair (BAM flag); with --sq3 the start of a read 1 with a UMI
    // is the 5' end of its molecule (see drop_unsupported_single_exons)
    bool is_read1 = false;
    uint32_t len;
    // Keeps count for all reads (including paired and unpaired)
    float read_count;
    std::vector<float> pair_count; // keeps count for all paired reads
    // With --sq3: read_count and pair_count before PCR copies were folded (see
    // weight_reads_by_molecule); read_count_raw < 0 when not set
    float read_count_raw = -1;
    std::vector<float> pair_count_raw;
    // Keeps indeces for all pairs in assembly mode, or all reads that were
    // collapsed in merge mode
    std::vector<int32_t> pair_idx;
    std::vector<GSeg> segs; //"exons"
    // Not owned: points to the shared junctions of BundleData::junction
    std::vector<CJunction *> juncs;

    // SCELL
    CReadAln(char _strand = '.', int16_t _nh = 0, int32_t _bc = -1,
             char *_ub = NULL, int32_t r_start = 0, int32_t r_end = 0)
        : GSeg(r_start, r_end), strand(_strand), nh(_nh), bc(_bc),
            ub(_ub ? _ub : ""), len(0), read_count(0), pair_count(), pair_idx(),
            segs(), juncs() {}
    CReadAln(CReadAln& rd) : GSeg(rd.start, rd.end) { // copy contructor
      strand = rd.strand;
      nh = rd.nh;
      bc = rd.bc; // SCELL
      ub = rd.ub; // SCELL
      is_read1 = rd.is_read1;
      len = rd.len;
      read_count = rd.read_count;
      pair_count = rd.pair_count;
      pair_idx = rd.pair_idx;
      juncs = rd.juncs;
    }

    ~CReadAln() {
    }
  };

  // The read alignments of one bundle, in the order they were read; the list
  // owns them.
  using ReadList = std::vector<std::unique_ptr<CReadAln>>;

  // An alignment from the BAM file with its strand, hit counts and introns.
  struct GReadAlnData {
    GSamRecord *brec;
    char strand; // '-', '+' or '.' (unknown)
    int32_t nh;
    int32_t hi;
    // The introns of this alignment (owned)
    std::vector<std::unique_ptr<CJunction>> juncs;

    GReadAlnData(GSamRecord *bam_rec = NULL, char strand_in = '.',
                 int32_t num_hits = 0,
                 int32_t hit_idx = 0)
        : brec(bam_rec), strand(strand_in), nh(num_hits), hi(hit_idx), juncs() {}
    ~GReadAlnData() {}
  };

  // The connected read groups of one strand, which become one splice graph:
  // pre-graph b of a strand is graph b of that strand. Only '-' and '+'
  // pre-graphs can become graphs, and only if they are worth processing
  // (is_pre_graph_worth_processing); the unknown-strand ('.') pre-graphs are
  // predicted directly and split between the strands.
  struct PreGraph {
    int32_t len;
    float cov;
    float multi;
    int32_t start_node; // id of start node in bundle of same strand
    int32_t last_node_id; // id of last node added to bundle
    PreGraph(int32_t _len = 0, float _cov = 0, float _multi
                = 0, int32_t _start = -1,
              int32_t _last = -1)
        : len(_len), cov(_cov), multi(_multi), start_node(_start),
          last_node_id(_last) {}
  };

  // A predicted transcript; one per cell, merged when the bundle is printed.
  struct CPrediction : public GSeg {
    int32_t gene_no;
    GffObj *t_eq; // equivalent reference transcript (guide)
    int32_t cell; // SCELL id; -1 if uninitiated
    double cov;
    char strand;
    int32_t tlen;
    bool flag;
    std::vector<GSeg> exons;
    std::vector<float> exon_cov;
    // SCELL this will stay empty until print_results
    std::vector<float> cell_cov;
    // Guides only: the coverage the flow gave each cell, kept when
    // credit_guides_found_by_denovo clears the cells where de novo does not find
    // the guide (see credit_best_cell_to_guides_without_cells)
    std::vector<float> flow_cov;
    // The graph path (node ids) the prediction came from, kept until the
    // graph's cells are done (see drop_shadowed_paths)
    std::vector<int32_t> path_nodes;
    // Coverage before PCR copies were folded, divided by cov (1 without --sq3)
    float raw_ratio = 1;
    // SCELL
    CPrediction(int32_t _gene_no = 0, GffObj *guide = NULL, int32_t _cell = -1,
                int32_t g_start = 0, int32_t g_end = 0, double _cov = 0,
                char _strand = '.', int32_t _len = 0, bool f = true)
        : GSeg(g_start, g_end), gene_no(_gene_no), t_eq(guide), cell(_cell),
          cov(_cov), strand(_strand), tlen(_len), flag(f), exons(),
          exon_cov(), cell_cov(), flow_cov(), path_nodes() {}
    // SCELL
    CPrediction(const CPrediction& c)
        : GSeg(c.start, c.end), gene_no(c.gene_no), t_eq(c.t_eq), cell(c.cell),
          cov(c.cov), strand(c.strand), tlen(c.tlen), flag(c.flag),
          exons(c.exons), exon_cov(c.exon_cov), cell_cov(c.cell_cov),
          flow_cov(c.flow_cov), path_nodes(c.path_nodes),
          raw_ratio(c.raw_ratio) {}
    ~CPrediction() { // GFREE(id);
    }
  };

  // The predicted transcripts of one bundle; the list owns them.
  using PredictionList = std::vector<std::unique_ptr<CPrediction>>;

  // A bundle: all the data parsed from the BAM file for one genomic region.
  struct BundleData {
    int32_t idx; // index in the main bundles array
    int32_t start;
    int32_t end;
    unsigned long num_reads; // number of reads in this bundle
    double num_fragments; // aligned read/pairs
    double frag_len;
    // Sum of all transcripts coverages --> needed to compute TPMs
    double sum_cov;
    // frag_len and the sum of the printed transcripts' coverage of each cell,
    // indexed by cell number (see cell_name); they give the per-cell FPKM and
    // TPM (see write_final_output)
    std::vector<double> cell_frag_len;
    std::vector<double> cell_sum_cov;

    std::string ref_seq; // reference sequence name
    char *gseq; // actual genomic sequence for the bundle
    ReadList reads;
    // Per-bp coverage of each strand and of all reads (see BpCoverage)
    BpCoverage bp_coverage;
    JunctionList junction;
    std::vector<GffObj *> guides; // guides in this bundle (not owned)
    PredictionList predictions;
    std::unique_ptr<BundleCounts> counts; // read count data for this bundle
    // Cell numbers by cell barcode (numbered in the order first seen), shared
    // by all the per-cell structures.
    std::unordered_map<std::string, int32_t> cell_name;
    BundleData()
        : idx(0), start(0), end(0), num_reads(0),
          num_fragments(0), frag_len(0), sum_cov(0), ref_seq(),
          gseq(NULL), reads(), junction(),
          predictions(), counts(), cell_name() {
      for (char strand : {'-', '.', '+'}) bp_coverage[strand].reserve(4096);
    }

    // Marks the bundle as ready to process and tags its guides.
    void get_ready(int32_t current_start, int32_t current_end) {
      // This is only called when the bundle is valid and ready to be processed
      start = current_start;
      end = current_end;
      // Tag all these guides
      for (int32_t i = 0; i < count_of(guides); ++i) {
        set_guide_status(guides[i], GBST_IN_BUNDLE);
      }
    }

    // Creates the count data of the bundle the first time it is needed.
    void init_counts(GffObj *t) {
      if (!counts) {
        counts = std::make_unique<BundleCounts>(t->start, t->end);
      }
    }
    // Registers guide as a reference transcript to keep for this bundle: lazily
    // initializes the bundle's per-transcript counters on first call and adds
    // guide to guides.
    void add_guide(GffObj *guide) {
      if (!counts) {
        init_counts(guide);
      }
      guides.push_back(guide);
      guide->udata = (int32_t)counts->add_transcript(*guide);
      // This also adds exon/intron info
    }

    // Checks the exons and introns of one read alignment against the guide
    // exons and introns of this bundle: counts the matched introns, collects
    // the alignment's junctions and refines splice_strand from the guides.
    void eval_read_aln(GReadAlnData& aln_data, char& splice_strand) {
      // Skip alignments that no guide can match.
      if (!counts) {
        return; // no ref transcripts available for this reads' region
      }

      GSamRecord& brec = *(aln_data.brec);
      int32_t mate_pos = brec.mate_start();

      if ((int32_t)brec.end < counts->bundle_start
          || (int32_t)brec.start > counts->cluster_end) {
        return; // hit outside coverage area
      }

      if (count_of(counts->guide_exons) == 0
          || counts->transcripts.empty()) {
        return; // nothing to do without transcripts
      }

      // Find the guide exons that each read exon overlaps and match the introns
      // between the exons to the guide introns.
      char strand_bits = 0;
      for (int32_t i = 0; i < brec.exons.Count(); i++) { // for each read exon

        std::vector<ExonOverlap> exon_overlaps;

        if (counts->find_overlapping_exons(exon_overlaps, brec.exons[i].start,
                                           brec.exons[i].end, splice_strand,
                                           mate_pos)) {
          // Largest guide exon overlap
          int32_t best_overlap = exon_overlaps[0].overlap_len;
          int32_t mate_overlap = exon_overlaps[0].mate_overlap;

          // Check all overlaps of similar length
          for (int32_t k = 0; k < count_of(exon_overlaps); ++k) {
            bool is_worse_overlap
                = exon_overlaps[k].mate_overlap < mate_overlap
                  || exon_overlaps[k].overlap_len + 5 < best_overlap;
            if (k && is_worse_overlap)
              break; // NOTE ignore further overlaps after a mate matched
            // Or if they are shorter than max_overlap-5
            if (exon_overlaps[k].feature->strand == '+') strand_bits |= 0x01;
            else if (exon_overlaps[k].feature->strand == '-') strand_bits |= 0x02;
          } // for each overlap
        } // has guide exon overlaps
        if (i > 0) { // adjacent intron check for match to guide introns
          int32_t j_l = brec.exons[i - 1].end + 1;
          int32_t j_r = brec.exons[i].start - 1;
          aln_data.juncs.push_back(std::make_unique<CJunction>(j_l, j_r));
          FeatureCounts *ri = counts->find_intron(j_l, j_r, splice_strand);
          if (ri) { // update guide intron counts
            ri->read_count++;
            aln_data.juncs.back()->guide_match = 1;
          }
        } // neighboring intron processing
      } // for each read exon

      // An unstranded alignment takes the strand of the guides it overlaps.
      if (splice_strand == '.' && strand_bits && strand_bits < 3) {
        splice_strand = (strand_bits == 1) ? '+' : '-';
      }
    }

    // Resets the bundle so that its slot can be reused.
    void clear() {
      // Clear the lists of guides, predictions, reads and junctions.
      guides.clear();
      predictions.clear();
      reads.clear();
      for (char strand : {'-', '.', '+'}) {
        bp_coverage[strand].clear();
        bp_coverage[strand].reserve(1024);
      }
      junction.clear();
      // Cell ids are local to one bundle; a recycled slot must not carry ids
      // over from earlier bundles
      cell_name.clear();

      // Reset the counters and free the genomic sequence and the count data.
      start = 0;
      end = 0;
      num_reads = 0;
      num_fragments = 0;
      frag_len = 0;
      sum_cov = 0;
      cell_frag_len.clear();
      cell_sum_cov.clear();
      counts.reset();
      GFREE(gseq);
    }

    ~BundleData() {
      clear();
    }
  };
} // namespace ribbon

#define MAX_NODE 1000000

#define DROP 0.5
#define ERROR_PERC 0.1

#define CHI_WIN 100
#define CHI_THR 50
#define SCELL_COV 20
#define SCELL_TOL 0.3 // scell tollerance for errors v17
// Number of single cells to accept as evidence for single exons
#define SCELL_NO 50

const double epsilon = 0.000001; //-E
const float min_transfrag_abundance = 1.0; // transfrag path_bits threshold

// Don't trust introns longer than this unless there is higher evidence; about
// 99% of all human annotated introns are shorter than this
const uint32_t long_intron = 100000;
// We need a higher anchor for long introns -> use a smaller value here? i.e.
// 20?
const uint32_t long_intron_anchor = 25;

const float mismatch_frac = 0.04; // v23 for single cell
const float low_iso_frac = 0.02;

// Maximum number of transfrag accepted so that the memory doesn't blow up
const int32_t max_transfrag_number = 40000;

namespace ribbon {
  // A segment of a bundle (a genomic interval) with its coverage, linked to
  // the next segment of the same bundle.
  struct PreGraphSegment : public GSeg {
    float cov;
    int32_t bid; // pre-graph segment id, to retrieve it easily
    PreGraphSegment *next_node; // next node in the same bundle
    PreGraphSegment(int32_t r_start = 0, int32_t r_end = 0, float _cov
                   = 0, int32_t _bid = -1,
                PreGraphSegment *_next_node = NULL)
        : GSeg(r_start, r_end), cov(_cov), bid(_bid), next_node(_next_node) {}
  };

  // What does not depend on the cell in the flow of a guide through the graph.
  struct GuideFlowCache {
    // Bit t is set if the path_bits of the guide contains the path_bits of
    // transfrag t of the graph
    BitVec contained_transfrags;
    // For each node of the guide: its transfrags with abundance that the guide
    // contains, in the order of the node's transfrags
    std::vector<std::vector<int32_t>> path_transfrags;
    // For each node of the guide: the sum of the abundance of the transfrags
    // with abundance that the guide does not contain, which enter the node
    // (pooled_in) or leave it (pooled_out)
    std::vector<float> pooled_in;
    std::vector<float> pooled_out;
  };

  // A path through the splice graph that reads support: its nodes, its
  // path_bits of nodes and edges, its abundance and the abundance in each cell.
  struct CTransfrag {
    int32_t cell_count;
    std::vector<int32_t> nodes;
    BitVec path_bits;
    float abundance;
    std::vector<float> cell_abundance; // SCELL -> abundance in each single cell
    // Share of the transfrag's abundance that the current guide path uses.
    float use_path;
    bool real : 1;
    // True for the link to the source or sink that a guide's start or end made,
    // as opposed to one that the reads support.
    bool made_by_guide = false;
    int32_t guide;
    // For the transfrag of a guide: what the flow functions reuse from cell to
    // cell (see prepare_guide_flow_cache in flow.cpp)
    std::shared_ptr<GuideFlowCache> flow_cache;
    // For long reads: min start of all longreads sharing transfrag
    uint32_t long_start;
    // for long reads: max end of all longreads sharing transfrag
    uint32_t long_end;
    // SCELL
    CTransfrag(int32_t _ncell, std::vector<int32_t>& _nodes, BitVec& bit,
               float abund = 0, bool treal = false, int32_t t_guide = 0)
        : cell_count(_ncell), nodes(_nodes), path_bits(bit), abundance(abund),
          use_path(-1), real(treal), guide(t_guide), long_start(0),
          long_end(0) {
      // SCELL : set all cell abundances to 0 initially
      cell_abundance.resize(_ncell, 0.0);
    }
    // SCELL
    CTransfrag(int32_t _ncell = 0, float abund = 0, bool treal = false,
               int32_t t_guide = 0)
        : cell_count(_ncell), nodes(), path_bits(), abundance(abund),
          cell_abundance(), use_path(-1), real(treal), guide(t_guide),
          long_start(0), long_end(0) {
      // SCELL : set all cell abundances to 0 initially
      cell_abundance.resize(cell_count, 0.0);
    }
  };

  // The transfrags of one splice graph; the list owns them.
  using TransfragList = std::vector<std::unique_ptr<CTransfrag>>;

  // A guide matched to the graph: its transfrag and its index in the guides.
  // CGuide is copied and sorted by value, so the copies share the transfrag.
  struct CGuide {
    std::shared_ptr<CTransfrag> transfrag;
    int32_t g; // stores guide index in guides instead of the actual pointer
    CGuide(std::shared_ptr<CTransfrag> _transfrag = nullptr, int32_t _g = -1)
        : transfrag(_transfrag), g(_g) {}
  };

  // A group of reads of the same strand that overlap or are close together.
  struct CGroup : public GSeg {
    // Label shared by the groups linked by a spliced read or a read pair; the
    // groups of one PreGraph end up with the same (root) link id
    int32_t link_id;
    int32_t grid;
    float cov_sum;
    float multi;
    // Proportion of negative reads assigned to group out of all positives and
    // negatives
    float neg_prop;
    CGroup *next_gr;
    CGroup(int32_t r_start = 0, int32_t r_end = 0, int32_t _link_id
           = -1, int32_t _grid = 0,
           float _cov_sum = 0, float _multi = 0, float _neg_prop = 0,
           CGroup *_next_gr = NULL)
        : GSeg(r_start, r_end), link_id(_link_id), grid(_grid),
          cov_sum(_cov_sum), multi(_multi), neg_prop(_neg_prop),
          next_gr(_next_gr) {}
  };

  // The read groups of one bundle, indexed by group id (CGroup::grid). The list
  // owns the groups; a merged-away group is reset to null instead of being
  // erased, so the other group ids stay valid.
  using GroupList = std::vector<std::unique_ptr<CGroup>>;

  // Position of a pre-graph segment in a graph: the graph and the node number.
  struct CGraphinfo {
    int32_t n_graph;
    int32_t node_no;
    CGraphinfo(int32_t ng = -1, int32_t nnode = -1) : n_graph(ng),
        node_no(nnode) {}
  };

  // Hash of a path key (see make_path_key).
  struct PathKeyHash {
    size_t operator()(const std::vector<int32_t>& key) const {
      size_t hash = key.size();
      for (int32_t value : key) {
        hash ^= std::hash<int32_t>()(value) + 0x9e3779b97f4a7c15ULL
            + (hash << 6) + (hash >> 2);
      }
      return hash;
    }
  };

  // Finds the transfrag of a graph with a given path (see make_path_key). It
  // does not own the transfrags: the TransfragList does.
  using TransfragIndex
      = std::unordered_map<std::vector<int32_t>, CTransfrag *, PathKeyHash>;

  struct GuideBoundary { // a start or end of a guide transcript
    // If pos < exon_other_end then this is a start; otherwise it is an end
    uint32_t pos; // genomic position of the guide's start or end
    // Other end of the guide's first exon (start) or last exon (end)
    uint32_t exon_other_end;
    char strand;
    bool operator<(const GuideBoundary& o) const {
      return (pos < o.pos
              || (pos == o.pos
                  && strand_index(strand) < strand_index(o.strand)));
    }
    bool operator==(const GuideBoundary& o) const {
      return (pos == o.pos && strand == o.strand);
    }
    GuideBoundary(uint32_t _pos = 0, uint32_t _exon_other_end = 0,
                  char _strand = '.')
        : pos(_pos), exon_other_end(_exon_other_end), strand(_strand) {}
  };

  // A node of the splice graph: an interval with its coverage, children,
  // parents and the transfrags that pass through it.
  struct CGraphnode : public GSeg {
    int32_t node_id;
    float cov;
    std::vector<float> cell_cov; // SCELL defines single cell coverage of node
    // The same before PCR copies were folded (equal to cell_cov without --sq3)
    std::vector<float> cell_cov_raw;
    std::vector<int32_t> child;
    std::vector<int32_t> parent;
    // The nodes and edges that can be reached from this node going toward
    // the sink (child) or toward the source (parent)
    BitVec child_reachable_bits;
    BitVec parent_reachable_bits;
    // The indices of the transfrags that pass the node
    std::vector<int32_t> transfrag_ids;
    // SCELL
    CGraphnode(int32_t cell_count = 0, int32_t s = 0, int32_t e = 0,
               uint32_t id = MAX_NODE, float node_cov = 0)
        : GSeg(s, e), node_id(id), cov(node_cov), child(), parent(),
          child_reachable_bits(), parent_reachable_bits(), transfrag_ids() {
      // SCELL : set all cell coverages to 0 initially
      cell_cov.resize(cell_count, 0.0);
      cell_cov_raw.resize(cell_count, 0.0);
    }
  };

  // The nodes of one splice graph, indexed by node id; the list owns them.
  using GraphNodeList = std::vector<std::unique_ptr<CGraphnode>>;

  // Assembles the transcripts of one bundle (defined in core.cpp).
  void infer_transcripts(BundleData *bundle);
} // namespace ribbon

extern std::string cram_ref; // reference genome FASTA for CRAM input

// Holds the next unread BAM record for one input file, so FileIO
// can pick, across all input files, whichever record comes first in
// genomic order.
struct InputRecord {
  std::unique_ptr<GSamRecord> brec; // owned here
  int32_t file_idx; // index of the source file in FileIO's files/readers
  bool operator<(const InputRecord& o) const {
    // Decreasing location sort
    GSamRecord& r1 = *brec;
    GSamRecord& r2 = *(o.brec);
    int32_t ref_cmp = r1.refId() - r2.refId();
    if (ref_cmp == 0) {
      // Higher coords first
      if (r1.start != r2.start)
        return (r1.start > r2.start);
      else {
        if (r1.end != r2.end)
          return (r1.end > r2.end);
        else if (file_idx == o.file_idx)
          return (strcmp(r1.name(), r2.name()) > 0);
        else
          return file_idx > o.file_idx;
      }
    } else { // use header order
      return (ref_cmp > 0);
    }
  }

  InputRecord(GSamRecord *b = NULL, int32_t i = 0) : brec(b), file_idx(i) {}
};

// Merges one or more input files into a single stream of BAM records
// ordered by genomic location (start() then repeated next() calls).
struct FileIO {
protected:
  // Record returned by the last next() call
  std::unique_ptr<InputRecord> current_rec;

public:
  // One reader per input file; owned here
  std::vector<std::unique_ptr<GSamReader>> readers;
  std::vector<std::string> files; // input file paths, same order as readers
  // Next unread record for each file, kept as a heap whose top is the record
  // that comes first in genomic order
  std::vector<std::unique_ptr<InputRecord>> pending_recs;
  FileIO() : current_rec(), readers(), files(), pending_recs() {}

  // Heap order of the pending records: the first in genomic order is on top.
  static bool record_less(const std::unique_ptr<InputRecord>& a,
                          const std::unique_ptr<InputRecord>& b) {
    return *a < *b;
  }

  // Adds a record to the heap of pending records.
  void push_pending(std::unique_ptr<InputRecord> rec) {
    pending_recs.push_back(std::move(rec));
    std::push_heap(pending_recs.begin(), pending_recs.end(), record_less);
  }

  // Registers an input file path (or "-" for stdin), after checking it exists.
  void add(const char *fn) {
    std::string sfn(fn);
    // Must exist and not be a directory
    bool is_stdin = sfn == "-";
    if (!is_stdin
        && (!std::filesystem::exists(fn) || std::filesystem::is_directory(fn))) {
      fprintf(stderr, "Error: input file %s cannot be found!\n", fn);
      exit(1);
    }
    files.push_back(sfn);
  }

  // Opens all registered input files and loads the first record from each.
  void start() {
    std::vector<std::string> bam_files = files;

    // Open every input file and load the first record of each
    for (int32_t i = 0; i < count_of(bam_files); ++i) {
      readers.push_back(std::make_unique<GSamReader>(
          bam_files[i].c_str(), cram_ref.empty() ? NULL : cram_ref.c_str()));
      GSamRecord *brec = readers.back()->next();
      if (brec)
        push_pending(std::make_unique<InputRecord>(brec, i));
    }
  }

  // Returns the next record across all input files in genomic order (a k-way
  // merge over each file's next unread record), refilling from whichever file
  // the returned record came from. Returns NULL once every file is exhausted.
  GSamRecord *next() {
    // Must free old current record first
    current_rec.reset();
    if (!pending_recs.empty()) {
      // Lowest coordinate
      std::pop_heap(pending_recs.begin(), pending_recs.end(), record_less);
      current_rec = std::move(pending_recs.back());
      pending_recs.pop_back();
      GSamRecord *next_bam_rec = readers[current_rec->file_idx]->next();
      if (next_bam_rec)
        push_pending(std::make_unique<InputRecord>(next_bam_rec,
                                                   current_rec->file_idx));
      return current_rec->brec.get();
    } else return NULL;
  }

  // Closes all input readers.
  void stop() {
    for (int32_t i = 0; i < count_of(readers); ++i) {
      readers[i]->bclose();
    }
  }
};

// The algorithmic components, split into their own files:
#include "path_bits.h" // which bit of a path stands for which node or edge
#include "graph.h" // types shared by the graph-building stages
#include "coverage.h" // bp coverage and coverage-drop trim points
#include "junctions.h" // junction support and filtering
#include "grouping.h" // reads into groups
#include "pre_graphs.h" // groups into pre-graphs
#include "graph_pruning.h" // pruning of graphs with too many nodes
#include "splice_graph.h" // splice graph construction
#include "flow.h" // types shared by the transcript extraction stages
#include "low_coverage.h" // low coverage nodes and edges
#include "path_flow.h" // flow along a path, storing transcripts
#include "guide_flow.h" // flow for reference (guide) transcripts
#include "path_extender.h" // greedy extension of a seed node into a path
#include "transfrag.h" // transfrag path_bits, simplification/stitching
#include "guides.h" // guide-transcript matching
#include "options.h" // command-line parsing
#include "read_input.h" // BAM records into reads
