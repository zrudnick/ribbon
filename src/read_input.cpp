#include "ribbon.h"

// Turns BAM records into reads (CReadAln): strand from the tags, exons and
// junctions from the CIGAR, junction checks against the MD tag and pairing of
// mates.

// --sq3: SmartSeq3 data, where read 1 with a UMI has the orientation of the
// transcript
extern bool smartseq3;
extern bool fr_strand;
extern bool rf_strand;
// Anchor length for junction to be considered well supported
extern uint32_t junction_support;

namespace ribbon {

  // To prevent repeated reallocation for each parsed read
  static std::string pair_key;

  // Result of checking a mismatch position against the exon boundaries.
  enum JunctionScanResult { JS_OK, JS_VIOLATION, JS_OUT_OF_SEGS };

  // Returns the strand ('+', '-' or '.') for a read: the tagged splice strand
  // if present, otherwise inferred from the library type (--fr/--rf), read
  // pairing and orientation.
  char get_strand(GSamRecord *brec) {
    char strand = brec->spliceStrand(); // tagged strand gets priority
    // Set strand if stranded library
    if (strand == '.' && (fr_strand || rf_strand)) {
      bool is_rev = brec->revStrand();
      // The orientation that means the + strand for a single read or for the
      // first read of a pair, given the library type.
      bool is_plus_orientation = (rf_strand && is_rev) || (fr_strand && !is_rev);
      if (brec->isPaired()) { // read is paired
        if (brec->pairOrder() == 1) { // first read in pair
          strand = is_plus_orientation ? '+' : '-';
        } else {
          strand = is_plus_orientation ? '-' : '+';
        }
      } else {
        strand = is_plus_orientation ? '+' : '-';
      }
    }
    return strand;
  }

  // Returns whether two exon-coordinate lists are identical (same count, same
  // start/end per exon).
  bool exon_match(std::vector<GSeg>& prev_exons, GVec<GSeg>& exons) {
    if (count_of(prev_exons) != exons.Count()) return false;

    for (int32_t i = 0; i < exons.Count(); i++) {
      if (prev_exons[i].end != exons[i].end
        || prev_exons[i].start != exons[i].start) {
          return false;
      }
    }

    return true;
  }

  // Checks whether a stored read has the same deletions next to its splice
  // sites as a new read (new_read_deletions: per intron, the deletion before
  // it in .start and after it in .end), so that such reads are not merged.
  bool junction_deletions_match(CReadAln *read,
                                GVec<GSeg>& new_read_deletions) {
    // Reads with a different number of junctions cannot share a deletion
    // pattern
    if (count_of(read->juncs) != new_read_deletions.Count()) return (false);

    // Compare the deletion lengths around each junction.
    for (int32_t i = 0; i < count_of(read->juncs); i++) {
      // Deletion length before junction i: gap between the end of exon i and
      // the junction start
      uint32_t start_j = read->segs[i].end - read->juncs[i]->start;
      // Deletion length after junction i: gap between the junction end and the
      // start of exon i+1
      uint32_t end_j = read->juncs[i]->end - read->segs[i + 1].start;
      if (start_j != new_read_deletions[i].start
          || end_j != new_read_deletions[i].end) return (false);
    }

    return (true);
  }

  // Records the exon segments of a new read as read's segments and registers
  // each junction between consecutive exons in the bundle's junction list
  // (adding the null junction first).
  static void process_exons(CReadAln *read, GSamRecord& brec,
                            JunctionList& junction, GReadAlnData& aln_data,
                            char strand) {
    for (int32_t i = 0; i < brec.exons.Count(); i++) {
      read->len += brec.exons[i].len();
      if (i) {
        // Register the intron before this exon in the bundle's junction list.
        if (!count_of(junction)) { // always add null junction first
          junction.push_back(std::make_unique<CJunction>(0, 0, '.'));
        }
        char j_strand = strand;
        uint32_t j_start = brec.exons[i - 1].end;
        uint32_t j_end = brec.exons[i].start;

        CJunction *shared_junction = add_junction_if_new(
            junction, std::make_unique<CJunction>(j_start, j_end, j_strand));
        if (count_of(aln_data.juncs))
          shared_junction->guide_match = aln_data.juncs[i - 1]->guide_match;
        if (shared_junction) {
          read->juncs.push_back(shared_junction);
        }
      }
      read->segs.push_back(brec.exons[i]);
    }
  }

  // Moves i (and rd_len) past the segments before parsed_len, then checks if the
  // segment at parsed_len lies within junction_support of its start or end;
  // extra_right shifts the end check for deletions.
  static JunctionScanResult advance_and_check_junction_support(CReadAln *rd,
      int32_t parsed_len, int32_t extra_right, int32_t& i, int32_t& rd_len) {
    while (i < count_of(rd->segs)
           && rd_len + (int32_t)rd->segs[i].len() < parsed_len) {
      rd_len += rd->segs[i].len();
      i++;
    }

    // Report a violation if the position is too close to an exon boundary.
    if (i == count_of(rd->segs)) return JS_OUT_OF_SEGS;
    bool too_close_to_start
        = i && parsed_len - rd->segs[i].start < junction_support;
    bool too_close_to_end = i < count_of(rd->segs) - 1
        && rd->segs[i].end + extra_right - parsed_len < junction_support;
    if (too_close_to_start || too_close_to_end) return JS_VIOLATION;
    return JS_OK;
  }

  // Returns true if the MD string of the read (mismatches and deletions)
  // has one within junction_support of an exon boundary.
  bool md_has_junction_violation(CReadAln *read, char *md_str,
                                 int32_t ref_start) {
    char *p = md_str;

    int32_t i = 0;
    int32_t parsed_len = ref_start;
    int32_t rd_len = 0;

    while (*p != '\0') {
      if (*p >= '0' && *p <= '9') {
        char *num_end = NULL;
        // Number of matching bases
        parsed_len += (uint32_t)strtoul(p, &num_end, 10);
        p = num_end;
        continue;
      }
      if (*p == '^') { // deletion --> found a problem with deletion

        int32_t del_length = 0; // tracking deletion length
        char del_base = *(++p);
        while (del_base >= 'A' && del_base <= 'Z') {
          del_length++;
          del_base = *(++p);
        }
        JunctionScanResult js = advance_and_check_junction_support(read,
            parsed_len, 1 - del_length, i, rd_len);
        if (js == JS_OUT_OF_SEGS) break;
        if (js == JS_VIOLATION) return true;
        parsed_len += del_length;
        continue;
      }
      if (*p >= 'A' && *p <= 'Z') { // mismatch
        JunctionScanResult js = advance_and_check_junction_support(read,
            parsed_len, 0, i, rd_len);
        if (js == JS_OUT_OF_SEGS) break;
        if (js == JS_VIOLATION) return true;
        parsed_len++;
      }
      p++;
    }

    return false;
  }

  // Returns true if the CIGAR of the read has an insertion within
  // junction_support of an exon boundary.
  bool cigar_has_junction_violation(CReadAln *read, bam1_t *b) {
    uint32_t *cigar = bam_get_cigar(b);
    int32_t i = 0;
    int32_t parsed_len = 0;
    int32_t rd_len = 0;

    for (uint32_t j = 0; j < b->core.n_cigar; ++j) {
      int32_t op = bam_cigar_op(cigar[j]);
      bool advances_reference = op == BAM_CMATCH || op == BAM_CEQUAL
                                || op == BAM_CDIFF || op == BAM_CDEL;
      if (advances_reference) {
        parsed_len += bam_cigar_oplen(cigar[j]);
      } else if (op == BAM_CINS) {
        JunctionScanResult js = advance_and_check_junction_support(read,
            parsed_len, 0, i, rd_len);
        if (js == JS_OUT_OF_SEGS) break;
        if (js == JS_VIOLATION) return true;
      }
    }

    return false;
  }

  // Returns true if the read has a mismatch, deletion or insertion (from MD
  // tag and CIGAR) within junction_support of an exon boundary, so that its
  // junctions are not trusted. Reads without an MD tag are not checked.
  bool mismatch_anchor(CReadAln *read, char *md_str, int32_t ref_start,
                       bam1_t *b) {
    if (md_str == NULL) return false;
    return md_has_junction_violation(read, md_str, ref_start)
        || cigar_has_junction_violation(read, b);
  }

  // Records that read `r` is paired with the read at `partner_idx`, adding
  // read_count to the existing pair_count entry if that partner is already
  // known, otherwise adding a new pair_idx/pair_count entry.
  static void add_or_update_pair(CReadAln *r, int32_t partner_idx,
                                 float read_count) {
    bool not_found = true;
    for (int32_t i = 0; i < count_of(r->pair_idx); i++)
      if (r->pair_idx[i] == partner_idx) {
        r->pair_count[i] += read_count;
        not_found = false;
        break;
      }
    if (not_found) { // we didn't see the pairing before
      r->pair_idx.push_back(partner_idx);
      r->pair_count.push_back(read_count);
    }
  }

  // Resolves pairing between read n and its mate: if the mate was already seen
  // (stored in hash_read), records the pair on both reads; otherwise stores read
  // n so its mate can find it later.
  static void process_pairs(ReadList& reads, int32_t read_idx,
                            GSamRecord& brec,
                            std::unordered_map<std::string, int32_t>& hash_read,
                            int32_t curr_start, int32_t read_start, int32_t nh,
                                int32_t hi,
                            float read_count) {
    // Only consider mate pairing data if mates are on the same
    // chromosome/contig and are properly paired
    if (brec.refId() == brec.mate_refId()) {
      int32_t mate_start = brec.mate_start();
      // If pairstart is in a previous bundle we don't care about it
      if (curr_start <= mate_start) {
        // Build the key that identifies this read and its mate.
        pair_key = brec.name();
        // If we've seen the pair already <- we might not have seen it yet
        // because the pair starts at the same place
        if (mate_start <= read_start) {
          // Look up the mate by its key (name, start and hit index).
          pair_key += '-';
          pair_key += std::to_string(mate_start);
          pair_key += ".=";
          // (!) this suffix actually speeds up the hash by improving
          // distribution!
          pair_key += std::to_string(hi);
          auto mate_it = hash_read.find(pair_key.c_str());
          const int32_t *mate_idx_ptr =
              mate_it == hash_read.end() ? NULL : &mate_it->second;
          // The pair was stored --> why wouldn't it be? : only in the case that
          // the pair starts at the same position
          if (mate_idx_ptr) {
            if (reads[*mate_idx_ptr]->nh > nh)
              read_count = float(1) / reads[*mate_idx_ptr]->nh;
            add_or_update_pair(reads[*mate_idx_ptr].get(), read_idx,
                               read_count);
            add_or_update_pair(reads[read_idx].get(), *mate_idx_ptr,
                               read_count);
            hash_read.erase(mate_it);
          }
        } else { // we might still see the pair in the future

          // Store this read under its own key so its mate can find it.
          pair_key += '-';
          pair_key += std::to_string(read_start); // this is the correct way
          pair_key += ".=";
          pair_key += std::to_string(hi);
          hash_read.emplace(pair_key.c_str(), read_idx);
        }
      }
    } //<-- if mate is mapped on the same chromosome
  }

  // Adds one BAM alignment to the bundle: builds its CReadAln (or updates the
  // identical read that it duplicates), registers its junctions and their
  // support and pairs it with its mate (process_pairs).
  void process_read_in(int32_t curr_start, int32_t curr_end, BundleData& bundle,
                       std::unordered_map<std::string, int32_t>& hash_read,
                       GReadAlnData& aln_data) {
    // Get the BAM record, the bundle lists and the properties of the alignment.
    GSamRecord& brec = *(aln_data.brec); // bam record
    ReadList& reads = bundle.reads; // list of reads gathered so far
    JunctionList& junction = bundle.junction; // junctions added so far
    char strand = aln_data.strand;
    int32_t nh = aln_data.nh;
    int32_t hi = aln_data.hi;
    int32_t read_start = brec.start;
    CReadAln *read = NULL; // read is initialized with NULL

    // Read the mismatch count and the cell barcode and UMI tags.
    double mismatch_count = (double)brec.tag_int("NM"); // read mismatch

    char *barcode = NULL;
    char *umi = NULL;
    // Read specific single cell data in process_read_in()
    barcode = brec.tag_str("BC"); // cell barcode
    umi = brec.tag_str("UB"); // UMI barcode

    if (!barcode || !strlen(barcode)) {
      fprintf(stderr, "No barcode present in read %s - ignoring\n",
              brec.name());
      return;
    }

    // Get the cell number from the cell barcode; a barcode not seen before in
    // this bundle gets the next free number
    auto cell_it = bundle.cell_name.find(barcode);
    if (cell_it == bundle.cell_name.end())
      cell_it = bundle.cell_name
                    .emplace(barcode, count_of(bundle.cell_name)).first;
    int32_t cell_id = cell_it->second;

    // In smartseq3 mode read 1 with umi has the same orientation as transcript
    if (smartseq3) {
      if (umi && strlen(umi)) { // read has UMI bar code
        char ustrand = '.';
        if (brec.flags() & BAM_FREAD1) { // read 1
          if (brec.flags() & BAM_FREVERSE) ustrand = '-';
          else ustrand = '+';
        }
        if (brec.flags() & BAM_FREAD2) { // read 2
          if (brec.flags() & BAM_FREVERSE) ustrand = '+';
          else ustrand = '-';
        }
        // Conflicting strands -> read is probably incorrect; do not use it
        if (strand != '.' && ustrand != strand) return;
        strand = ustrand;
      }
    }

    // Look for an identical read among the reads that start at the same place.
    bool match = false; // true if current read matches a previous read
    int32_t read_idx = count_of(reads) - 1;

    while (read_idx > -1 && reads[read_idx]->start == brec.start) {
      // Reads with the same barcode and UMI are considered the same
      if (reads[read_idx]->strand == strand
          && reads[read_idx]->bc == cell_id
          && reads[read_idx]->ub == (umi ? umi : "")) {
        match = exon_match(reads[read_idx]->segs, brec.exons);
        // DEL AWARE
        if (match)
          match = junction_deletions_match(reads[read_idx].get(),
                                           brec.juncsdel);
      }

      // This way we make sure that we keep the n of the matching reads
      if (match) break;
      read_idx--;
    }

    // Extend the bundle to this alignment and count it.
    if (bundle.end < curr_end) {
      bundle.start = curr_start;
      bundle.end = curr_end;
    }
    bundle.num_reads++; // number of reads gets increased no matter what

    // Set up a new read, or update the read that this alignment duplicates.
    if (!match) { // if this is a new read we are seeing we need to set it up

      std::unique_ptr<CReadAln> new_read = std::make_unique<CReadAln>(
          strand, nh, cell_id, umi, brec.start, brec.end);
      new_read->is_read1 = (brec.flags() & BAM_FREAD1) != 0;
      read = new_read.get();

      process_exons(read, brec, junction, aln_data, strand);
      reads.push_back(std::move(new_read));
      read_idx = count_of(reads) - 1; // reset n for the case there is no match
    } else { // redundant read alignment matching a previous alignment
      // Keep shortest nh so that we can see for each particular read the
      // multi-hit proportion
      if (nh < reads[read_idx]->nh) reads[read_idx]->nh = nh;
    }

    // Extend the end of the bundle if the alignment reaches past it.
    if ((int32_t)brec.end > curr_end) {
      curr_end = brec.end;
      bundle.end = curr_end;
    }

    // Add the alignment's weight (shared by its hits) to the read count.
    float read_count = (float)brec.tag_int("YC"); // alignment count
    if (!read_count) read_count = 1;

    read_count /= nh;
    // Increase single count just in case we don't see the pair
    reads[read_idx]->read_count += read_count;

    // Store the mismatch count per junction so that we can eliminate it later
    if (!mismatch_count) {
      // Paired mismatch : big problem with STAR alignments
      mismatch_count = (double)brec.tag_int("nM");
      if (brec.isPaired()) mismatch_count /= 2;
    }
    if (brec.clipL) mismatch_count++;
    if (brec.clipR) mismatch_count++;

    if (count_of(reads[read_idx]->juncs)) {
      bool mismatch = false;
      if (mismatch_count / reads[read_idx]->len > mismatch_frac) mismatch = true;
      else if (mismatch_count && count_of(reads[read_idx]->juncs)) {
        // Penalize mismatch that's too close to ss
        if (brec.clipL
            && reads[read_idx]->segs[0].len() < junction_support
            + brec.clipL) mismatch = true;
        else if (brec.clipR
                 && reads[read_idx]->segs.back().len() < junction_support
                 + brec.clipR) mismatch = true;
        // This line was not initially present in vs1 or vs3 but we noticed it
        // doesn't do any difference in real data, so far it only helped with
        // the SR in simulation -> we might want to take it out
        else if (mismatch_anchor(reads[read_idx].get(), brec.tag_str("MD"),
                                 curr_start, brec.get_b())) mismatch = true;
      }

      // If read is PacBio we might want to increase the mismatch fraction,
      // although the nm only gets used for longintrons
      for (int32_t i = 0; i < count_of(reads[read_idx]->juncs); i++) {
        if (mismatch || nh > 2) reads[read_idx]->juncs[i]->nm += read_count;
        if (reads[read_idx]->segs[i].len() > long_intron_anchor
            && reads[read_idx]->segs[i + 1].len() > long_intron_anchor)
          reads[read_idx]->juncs[i]->mm += read_count;
        reads[read_idx]->juncs[i]->n_reads += read_count;
      }
    }

    // Pair the read with its mate.
    process_pairs(reads, read_idx, brec, hash_read, curr_start, read_start, nh,
                  hi, read_count);
  }

} // namespace ribbon
