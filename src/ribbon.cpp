// Program entry point and front end of the pipeline: process_reads groups the
// sorted BAM records into bundles with their guides and hands them to the
// worker threads, and run drives everything. The command line is parsed in
// options.cpp and the BAM records are turned into reads in read_input.cpp.
#include "ribbon.h"
#include "out.h"

#include "CLI/CLI11.hpp"
#include "quill/Backend.h"
#include "quill/Frontend.h"
#include "quill/LogMacros.h"
#include "quill/Logger.h"
#include "quill/sinks/ConsoleSink.h"
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>
#ifndef NOTHREADS
#include "GThreads.h"
#include <atomic>
#endif
#include "threads.h"

//---- globals

// The output GTF (first the temporary results file, later the final output)
FILE *f_out = NULL;

// Command-line parameters (see process_options)

bool guided = false;
bool trim_by_coverage = true;
// --sq3: SmartSeq3 data, where read 1 with a UMI has the orientation of the
// transcript
bool smartseq3 = false;
// Weight of the per-cell abundance in the path-selection score (see
// extend_path)
double alpha = 1;
// Weight of the pooled (all-cells) abundance in the path-selection score
double beta = 0.75;
// --min-borrowed-share: a path is not stored in a cell if one of its junctions
// was chosen from other cells' reads with less than this share of the pooled
// support of the competing junctions (0: off; see check_weak_borrow)
float min_borrowed_share = 0.6;
// A path is not stored if the reads of all cells contradict two consecutive
// branch choices of it (is_path_contradicted); --keep-contradicted-paths: off
bool drop_contradicted_paths = true;
// A cell is taken off a novel path when its reads do not single the path out
// among the graph's paths (drop_shadowed_paths); --keep-shadowed-paths: off
bool drop_shadowed = true;

bool fr_strand = false;
bool rf_strand = false;
// minimum length for a transcript to be printed
int32_t min_transcript_len = 200;
// Anchor length for junction to be considered well supported (a shorter
// anchor might be worth trying)
uint32_t junction_support = 10;
// -E: how far (bp) from a poorly supported junction to look for a better
// supported one to replace it
uint32_t splice_site_error = 25;

// Minimum read coverage per bundle bp for a novel multi-exon transcript (-c;
// reference ones need min_reference_coverage, see drop_weak_predictions)
float min_read_coverage = 3;
// Minimum per-base read coverage a single-exon transcript needs to be trusted
// as real rather than noise
float min_single_exon_coverage = 4.75;
// Reads at what distance should be considered part of separate bundles
uint32_t min_bundle_distance = 50;
// Distance threshold to detect polymerase run-off artifacts
uint32_t runoff_distance = 200;
// Fraction of bundle allowed to be covered by multi-hit reads paper uses 1
float multi_read_frac = 1;
int32_t allowed_nodes = 1000;

// Number of reads needed to support a particular junction
int32_t min_junction_coverage = 1;

float isoform_fraction = 0.01;
std::string label("RBN");

bool QUIET = false;        // turn off progress logging, --quiet
bool RIBBON_DEBUG = false; // log processing details (bundles, warnings)
quill::Logger *logger = NULL;

int32_t gene_no = 0; //-- global "gene" counter
double total_fragments = 0; // global fragment counter (aligned pairs)
double total_frag_len = 0;
double total_cov_sum = 0;
// total_frag_len and the sum of the transcripts' coverage for each cell barcode, for
// the per-cell FPKM and TPM
std::map<std::string, double> frag_len_by_cell;
std::map<std::string, double> cov_sum_by_cell;
// With --demultiplex: the directory where we write one GTF per cell (empty:
// none), next to the output GTF
std::string demultiplexed_dir;

// Used as a dictionary for reference sequence names and ids
GffNames *gseq_names = NULL;

// number of reference sequences found in the guides file
int32_t n_ref_guides = 0;

// Chromosomes/contigs to exclude (e.g. chrM), from -x
std::unordered_set<std::string> excluded_seqs;

int32_t n_threads = 1;
std::unique_ptr<GFastaDb> g_fasta;
std::string guide_gff; // -G option
std::string out_fname;
std::string out_dir;
std::string tmp_path;
std::string cram_ref; // reference genome FASTA for CRAM input
std::string tmp_fname;
// Full command line, echoed in verbose log and output header
std::string cmd_line;

#ifndef NOTHREADS
// Threading: single producer, multiple consumers
// Main thread/program is always loading the producer

GMutex data_mutex; // guards the bundle slots that are ready to load
std::vector<int32_t> clear_data_pool; // indices of bundle slots free to load
GConditionVar have_bundles; // wakes a thread when a bundle is queued
char bundle_work = 1; // bit 0: bundles still being prepared; bit 1: queued

GMutex wait_mutex; // controls threads_waiting (idle threads counter)
std::atomic<uint8_t> threads_waiting; // idle worker threads
GConditionVar have_threads; // wakes the loader when a thread is idle

GConditionVar have_clear; // wakes the loader when a slot is free
GMutex queue_mutex; // controls bundle_queue and bundles access
GFastMutex print_mutex; // for writing the output to file
GFastMutex reading_mutex;
GFastMutex count_mutex;
#endif

bool no_more_bundles = false;

namespace ribbon {

  // Starts a bundle at the read at pos: skips the guides that end before it
  // and adds the guides that overlap the read, and transitively those that
  // overlap them, widening curr_start and curr_end. Updates the guide indices.
  static void add_guides_to_new_bundle(GList<GffObj> *guides,
                                       int32_t num_guides,
                                       int32_t pos, int32_t& curr_start,
                                       int32_t& curr_end,
                                       int32_t& first_guide_idx,
                                       int32_t& last_guide_idx,
                                       BundleData *bundle) {
    // Skip the guides that end before the read.
    first_guide_idx = last_guide_idx + 1;
    while (first_guide_idx < num_guides
           && (int32_t)(*guides)[first_guide_idx]->end < pos) {
      // For now, skip guides which have no overlap with current read
      first_guide_idx++;
    }
    int32_t overlap_guide_idx = first_guide_idx;
    // Add the guides that overlap the read, then the ones that overlap those,
    // and so on.
    while (overlap_guide_idx < num_guides
           && (int32_t)(*guides)[overlap_guide_idx]->start <= curr_end) {
      if (curr_start > (int32_t)(*guides)[overlap_guide_idx]->start)
        curr_start = (*guides)[overlap_guide_idx]->start;
      if (curr_end < (int32_t)(*guides)[overlap_guide_idx]->end)
        curr_end = (*guides)[overlap_guide_idx]->end;

      // For 1st guide overlap: needed to check previous guides for overlaps
      if (overlap_guide_idx == first_guide_idx && overlap_guide_idx > 0) {
        // First time only: we have to check back all possible transitive guide
        // overlaps
        // start from the overlapping guide, going backwards
        int32_t back_guide_idx = overlap_guide_idx;
        int32_t first_overlap_guide_idx = overlap_guide_idx;
        while (back_guide_idx > last_guide_idx + 1) {
          --back_guide_idx;

          // If overlap, set first_overlap_guide_idx=back_guide_idx and update
          // curr_start
          if (curr_start <= (int32_t)(*guides)[back_guide_idx]->end) {
            first_overlap_guide_idx = back_guide_idx;
            if (curr_start > (int32_t)(*guides)[back_guide_idx]->start)
              curr_start = (int32_t)(*guides)[back_guide_idx]->start;
          }
        } // while checking previous guides that could be pulled in this bundle

        // Add all guides that overlap the first added guide, in order
        for (int32_t bundle_guide_idx = first_overlap_guide_idx;
             bundle_guide_idx <= overlap_guide_idx; ++bundle_guide_idx) {
          bundle->add_guide((*guides)[bundle_guide_idx]);
        }
      } else { // for the rest, just add the guide to the bundle
        bundle->add_guide((*guides)[overlap_guide_idx]);
      }
      overlap_guide_idx++;
    } // while guide overlap

    // We update last_guide_idx even if no guide overlapped.
    last_guide_idx = overlap_guide_idx - 1;
  }

  // The current read extends the bundle to curr_end: adds any not-yet-assigned
  // guides that now overlap the bundle, extending curr_end further when a guide
  // reaches past it, until no more guides overlap. Updates last_guide_idx.
  static void add_guides_to_extended_bundle(GList<GffObj> *guides,
                                            int32_t num_guides,
                                                int32_t curr_start,
                                            int32_t& curr_end,
                                                int32_t& last_guide_idx,
                                            BundleData *bundle) {
    bool bundle_end_extended = true;
    while (bundle_end_extended) {
      bundle_end_extended = false;
      while (last_guide_idx + 1 < num_guides
             && (int32_t)(*guides)[last_guide_idx + 1]->start <= curr_end) {
        ++last_guide_idx;
        // More transcripts overlapping this bundle?
        if ((int32_t)(*guides)[last_guide_idx]->end >= curr_start) {
          // No transitive overlap check here, guide must overlap the bundle
          bundle->add_guide((*guides)[last_guide_idx]);
          if (curr_end < (int32_t)(*guides)[last_guide_idx]->end) {
            curr_end = (*guides)[last_guide_idx]->end;
            bundle_end_extended = true;
          }
        }
      }
    }
  }

  // Groups the alignments of all input files (read in genomic order) into
  // bundles, with their guides if guided, and queues each finished bundle for
  // the workers. Exits on unsorted input; no_ref_used turns false with guides.
  void process_reads(FileIO *io, std::vector<GRefData> &ref_guides,
                     BundleData *bundles, BundleData *bundle,
                     std::vector<BundleData *> *bundle_queue,
                     bool &no_ref_used) {
    // State of the loop: alignment counts per chromosome, the guides of the
    // current chromosome and the position of the current bundle.
    const char *ERR_BAM_SORT
        = "\nError: the input alignment file is not sorted!\n";

    // Keep track of the number of read alignments per chromosome [ref_id]
    std::vector<int32_t> aln_counts;
    // The read_name:pos:hit_index => reads index
    std::unordered_map<std::string, int32_t> hash_read;
    GList<GffObj> *guides = NULL; // list of transcripts on a specific reference
    int32_t curr_start = 0, curr_end = 0;
    int32_t first_guide_idx = 0;
    int32_t last_guide_idx = -1;
    int32_t num_guides = 0;
    std::string last_ref;
    int32_t last_ref_id = -1; // last seen ref_id

    GSamRecord *brec = NULL;
    bool more_alignments = true;
    int32_t prev_pos = 0;
    bool skip_gseq = false;
    while (more_alignments) {
      bool chr_changed = false;
      int32_t pos = 0;
      const char *ref_name = NULL;
      char strand = '.';
      int32_t nh = 1;
      int32_t hi = 0;
      int32_t ref_id = last_ref_id; // current chr id
      bool new_cluster = false;

      // Read the next alignment and check that it is mapped and sorted.
      if ((brec = io->next()) != NULL) {
        if (brec->isUnmapped()) continue;
        if (brec->start < 1 || brec->mapped_len < 10) {
          if (RIBBON_DEBUG)
            LOG_WARNING(logger,
                        "invalid mapping found for read {} (position={}, "
                        "mapped length={})",
                        brec->name(), brec->start, brec->mapped_len);
          continue;
        }

        ref_name = brec->refName();
        strand = get_strand(brec);

        pos = brec->start; // BAM is 0 based, but GBamRecord makes it 1-based
        chr_changed = (last_ref.empty() || last_ref != ref_name);

        // On a new chromosome, check that we have not seen it before.
        if (chr_changed) {
          skip_gseq = excluded_seqs.count(ref_name) > 0;
          ref_id = gseq_names->gseqs.addName(ref_name);
          if (guided) {
            if (ref_id >= n_ref_guides) {
              if (RIBBON_DEBUG)
                LOG_WARNING(logger,
                            "no reference transcripts found for genomic "
                            "sequence \"{}\"! (mismatched reference names?)",
                            ref_name);
            } else no_ref_used = false;
          }

          if (count_of(aln_counts) <= ref_id) {
            aln_counts.resize(ref_id + 1);
          } else if (aln_counts[ref_id] > 0) {
            fprintf(stderr, "%s\nAlignments (%d) already found for %s !\n",
                    ERR_BAM_SORT, aln_counts[ref_id], ref_name);
            exit(1);
          }
          prev_pos = 0;
        }
        if (pos < prev_pos) {
          fprintf(stderr,
                  "%s\nread %s (start %d) found at position %d on %s when "
                  "prev_pos=%d\n",
                  ERR_BAM_SORT, brec->name(), brec->start, pos, ref_name,
                  prev_pos);
          exit(1);
        }

        // Count the alignment and decide whether it starts a new bundle.
        prev_pos = pos;
        if (skip_gseq) continue;
        aln_counts[ref_id]++;
        nh = brec->tag_int("NH");
        if (nh == 0) nh = 1;
        hi = brec->tag_int("HI");

        if (!chr_changed && curr_end > 0
            && pos > curr_end + (int32_t)runoff_distance) {
          new_cluster = true;
        }
      } else { // no more alignments
        more_alignments = false;
        new_cluster = true; // fake a new start (end of last bundle)
      }

      if (new_cluster || chr_changed) { // read alignment starts a new bundle
        // Finish the previous bundle: get its sequence and queue it.
        hash_read.clear();
        if (count_of(bundle->reads) > 0) { // process reads in previous bundle

          bundle->get_ready(curr_start, curr_end);
          if (g_fasta) { // genomic sequence data requested
            GFaSeqGet *fa_seq = g_fasta->fetch(bundle->ref_seq.c_str());
            if (fa_seq == NULL) {
              fprintf(stderr,
                      "Error: could not retrieve sequence data for %s!\n",
                      bundle->ref_seq.c_str());
              exit(1);
            }
            bundle->gseq = fa_seq->copyRange(bundle->start, bundle->end, false,
                                            true);
          }
        } // have alignments to process

        // Enqueue for a worker thread, or clear it if empty
        push_bundle(bundle, bundle_queue);

        // Switch to the guides of the new chromosome.
        if (chr_changed) {
          if (guided) {
            num_guides = 0;
            guides = NULL;
            first_guide_idx = 0;
            last_guide_idx = -1;
            if (count_of(ref_guides) > ref_id
                && ref_guides[ref_id].rnas.Count() > 0) {
              guides = &(ref_guides[ref_id].rnas);
              num_guides = guides->Count();
            }
          }

          last_ref = ref_name;
          last_ref_id = ref_id;
          curr_end = 0;
        }

        // Stop after the last bundle: wait for the workers to finish.
        if (!more_alignments) {
          bool is_verbose = !QUIET || RIBBON_DEBUG;
          if (is_verbose && total_fragments) {
            LOG_INFO(logger, "{} aligned fragments found", total_fragments);
          }
          wait_for_bundles();
          break;
        }
#ifndef NOTHREADS

        // Get a free bundle slot for the new bundle.
        int32_t new_bidx = wait_for_data();
        if (new_bidx < 0) {
          // Should never happen!
          fprintf(stderr,
                  "Error: wait_for_data() returned invalid bundle index(%d)!\n",
                  new_bidx);
          exit(1);
          break;
        }
        bundle = &(bundles[new_bidx]);
#endif

        // Start the new bundle at this alignment with the guides that overlap.
        curr_start = pos;
        curr_end = brec->end;
        if (guides) { // guided and guides!=NULL
          add_guides_to_new_bundle(guides, num_guides, pos, curr_start,
                                   curr_end, first_guide_idx, last_guide_idx,
                                   bundle);
        } // guides present on the current chromosome
        bundle->ref_seq = last_ref;
        bundle->start = curr_start;
        bundle->end = curr_end;
      } //<---- read alignment started a new bundle

      // Read alignments extends the current bundle
      if (curr_end < (int32_t)brec->end) {
        // This might not happen if a longer guide had already been added to the
        // bundle
        curr_end = brec->end;
        if (guides) { // add any newly overlapping guides to bundle
          add_guides_to_extended_bundle(guides, num_guides, curr_start,
                                        curr_end, last_guide_idx, bundle);
        }
      } // adjusted curr_end and checked for overlapping reference transcripts

      // Check the alignment against the guides (they may set its strand), then
      // add it to the bundle.
      GReadAlnData aln_data(brec, '.', nh, hi);
      bundle->eval_read_aln(aln_data, strand);

      aln_data.strand = strand;
      process_read_in(curr_start, curr_end, *bundle, hash_read, aln_data);
    } // for each read alignment
  }

  // Closes and deletes the temporary results file and its directory, for the
  // runs that stop before write_final_output.
  void remove_temp_output() {
    fclose(f_out);
    remove(tmp_fname.c_str());
    while (!tmp_path.empty() && tmp_path.back() == '/') tmp_path.pop_back();
    std::error_code remove_error;
    std::filesystem::remove(tmp_path.c_str(), remove_error);
  }

  // Program driver called from main: parses the options, loads the guides if
  // given, runs process_reads so the workers assemble the bundles, writes the
  // final output (write_final_output) and cleans up. Returns the exit code.
  int32_t run(int32_t argc, char *argv[]) {
    // Start the logger and parse the command line.
    quill::Backend::start();
    quill::ConsoleSinkConfig sink_config;
    // Stdout carries the GTF when no -o is given
    sink_config.set_stream("stderr");
    auto console_sink
        = quill::Frontend::create_or_get_sink<quill::ConsoleSink>("sink_id_1",
        sink_config);
    logger = quill::Frontend::create_or_get_logger("root",
        std::move(console_sink));

    std::vector<std::string> in_files;
    process_options(argc, argv, in_files);

    // Register the input files.
    std::unique_ptr<FileIO> io = std::make_unique<FileIO>();
    for (const std::string& guide_file : in_files) io->add(guide_file.c_str());

    // Announce the version.
    if (!QUIET || RIBBON_DEBUG) {
      printf("\n[ribbon] starting version: %s\n", VERSION);
    }

    // Owns the count data of the guides (each guide points to its own).
    std::vector<std::unique_ptr<TranscriptCounts>> guide_transcript_counts;

    // The guides of each chromosome (only if guided).
    std::vector<GRefData> ref_guides;
    if (guided) { // read guiding transcripts from input gff file
      if (!QUIET || RIBBON_DEBUG) {
        LOG_INFO(logger, "loading reference annotation...");
      }

      // Open GFF/GTF file
      FILE *guide_file = fopen(guide_gff.c_str(), "r");
      if (guide_file == NULL) {
        LOG_ERROR(logger, "error: could not open reference annotation file {}",
          guide_gff.c_str());
        remove_temp_output();
        return 1;
      }

      // GffReader: transcripts only, sort by location
      // loading only recognizable transcript features
      GffReader gff_reader(guide_file, true, true);
      gff_reader.setRefAlphaSorted(); // alphabetical sorting of RefSeq IDs
      gff_reader.showWarnings(RIBBON_DEBUG);

      // Keep attributes, merge close exons, no exon attributes
      gff_reader.readAll(false, true, true);

      // The guides are in gff_reader.gflst, sorted by chromosome and start-end
      // coordinates. Check that the reader found some.
      n_ref_guides = gff_reader.gseqtable.Count();
      if (n_ref_guides == 0 || gff_reader.gflst.Count() == 0) {
        LOG_ERROR(logger,
                  "error: could not find valid reference transcripts in {}\n",
                  guide_gff.c_str());
        remove_temp_output();
        return 1;
      }

      // Make room for the guides of every chromosome.
      ref_guides.resize(n_ref_guides); // maximum reference guide ID

      uint32_t transcript_counter = 0; // counter for all transcripts
      // Assign unique transcript IDs based on the sorted order
      int32_t last_guide_gseq_id = -1;
      bool skip_gseq = false;
      for (int32_t i = 0; i < gff_reader.gflst.Count(); i++) {
        GffObj *guide = gff_reader.gflst[i];
        if (last_guide_gseq_id != guide->gseq_id) {
          // Chromosome switch
          last_guide_gseq_id = guide->gseq_id;
          skip_gseq = excluded_seqs.count(guide->getGSeqName()) > 0;
        }
        // Sanity check: make sure there are no exonless "genes" or other
        if (skip_gseq) continue;
        if (guide->exons.Count() == 0) {
          if (RIBBON_DEBUG)
            LOG_WARNING(logger,
                        "exonless GFF {} feature with ID {} found, added "
                        "implicit exon {}-{}",
                        guide->getFeatureName(), guide->getID(), guide->start,
                        guide->end);
          guide->addExon(guide->start, guide->end); // should never happen!
        }
        // Give each guide its own count data, reachable through uptr.
        guide_transcript_counts.push_back(
            std::make_unique<TranscriptCounts>(*guide, ++transcript_counter));
        guide->uptr = guide_transcript_counts.back().get();
        GRefData& g_ref_data = ref_guides[guide->gseq_id];
        // Add transcripts already sorted by location
        g_ref_data.add(&gff_reader, guide);
      }
      if (!QUIET || RIBBON_DEBUG) {
        LOG_INFO(logger,
            "reference annotation loaded! {} reference transcripts were found",
            gff_reader.gflst.Count());
      }
    }

    // Get the table of reference sequence names (shared with the guides).
    gseq_names = GffObj::names; // might have been populated already by gff data
    gffnames_ref(gseq_names); // initialize the names collection if not guided

    // Open the input files.

    io->start(); // setup and open input files

    std::vector<BundleData *> bundle_queue;

#ifndef NOTHREADS

#define DEF_TSTACK_SIZE 8388608
    size_t def_stack_size = DEF_TSTACK_SIZE;

#ifdef _GTHREADS_POSIX_
    size_t thread_stack_size = GThread::defaultStackSize();
    if (thread_stack_size < DEF_TSTACK_SIZE) def_stack_size = DEF_TSTACK_SIZE;
    if (RIBBON_DEBUG) {
      if (thread_stack_size < def_stack_size) {
        LOG_INFO(logger, "default stack size for threads: {}, increased to {}",
                 thread_stack_size, def_stack_size);
      } else {
        LOG_INFO(logger, "default stack size for threads: {}", thread_stack_size);
      }
    }
#endif

    // Allocate the worker threads and the bundle slots (one more than threads).
    // Threads for processing bundles
    std::unique_ptr<GThread[]> threads = std::make_unique<GThread[]>(n_threads);
    // Redef with more bundles
    std::unique_ptr<BundleData[]> bundle_slots
        = std::make_unique<BundleData[]>(n_threads + 1);
    BundleData *bundles = bundle_slots.get();

    clear_data_pool.reserve(n_threads + 1);

    // Start worker threads
    for (int32_t b = 0; b < n_threads; b++) {
      threads[b].kickStart(worker_thread, (void *)&bundle_queue,
                           def_stack_size);
      bundles[b + 1].idx = b + 1;
      clear_data_pool.push_back(b);
    }
    BundleData *bundle = &(bundles[n_threads]);

    if (RIBBON_DEBUG) {
      LOG_INFO(logger, "threads started");
    }

#else

    // Just put everything into the same bundle
    BundleData bundles[1];
    BundleData *bundle = &(bundles[0]);
#endif

    // Process all the alignments.
    if (!QUIET || RIBBON_DEBUG) {
      LOG_INFO(logger, "processing alignments :-)");
    }

    bool no_ref_used = true;
    process_reads(io.get(), ref_guides, bundles, bundle, &bundle_queue,
                  no_ref_used);

    // Wait for the workers and free the threads and the bundles.
#ifndef NOTHREADS
    for (int32_t t = 0; t < n_threads; t++) {
      threads[t].join();
    }
    if (RIBBON_DEBUG) {
      LOG_INFO(logger, "all threads finished");
    }
    threads.reset();
    bundle_slots.reset();
#endif

    // Close the input files.
    io->stop(); // close all BAM files
    io.reset();
    g_fasta.reset();

    // Warn if none of the guides matched the sequences of the alignments.
    if (guided && no_ref_used) {
      fprintf(stderr,
              "WARNING: no reference transcripts were found for the genomic "
              "sequences where reads were mapped!\n"
              "Please make sure the -G annotation file uses the same naming "
              "convention for the genome sequences.\n");
    }

    // Close the temporary results file and log the totals.
    fclose(f_out);

    if (!QUIET || RIBBON_DEBUG) {
      LOG_INFO(logger, "total count of aligned fragments: {}", total_fragments);
      if (total_fragments)
        LOG_INFO(logger, "fragment coverage length: {}",
                 total_frag_len / total_fragments);
    }

    // Open the final output file and write its header.
    f_out = stdout;
    if (out_fname != "stdout") {
      f_out = fopen(out_fname.c_str(), "w");
      if (f_out == NULL) {
        fprintf(stderr, "Error creating output file %s\n", out_fname.c_str());
        exit(1);
      }
    }

    fprintf(f_out, "# ");
    fprintf(f_out, "%s\n", cmd_line.c_str());
    fprintf(f_out, "# ribbon version %s\n", VERSION);

    // Add FPKM and TPM to the transcripts and write them to the output.
    write_final_output();

    // Remove the temporary directory.
    while (!tmp_path.empty() && tmp_path.back() == '/') tmp_path.pop_back();
    // Check the results of remove() and show error if it fails
    std::error_code remove_error;
    if (!std::filesystem::remove(tmp_path.c_str(), remove_error)
        || remove_error) {
      fprintf(stderr, "Error removing temporary directory: %s\n",
              tmp_path.c_str());
    }

    // Free the table of sequence names.
    gffnames_unref(gseq_names); // deallocate names collection

    return 0;
  }

} // namespace ribbon

int main(int32_t argc, char *argv[]) {
  return ribbon::run(argc, argv);
}