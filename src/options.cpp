#include "ribbon.h"

#include "CLI/CLI11.hpp"
#include "quill/LogMacros.h"
#include "quill/Logger.h"
#include <filesystem>

// Command-line parsing: process_options reads the options into the global
// parameters of the run.

// The output GTF (first the temporary results file, later the final output)
extern FILE *f_out;
extern bool guided;
extern bool trim_by_coverage;
// --sq3: SmartSeq3 data, where read 1 with a UMI has the orientation of the
// transcript
extern bool smartseq3;
extern bool fr_strand;
extern bool rf_strand;
// minimum length for a transcript to be printed
extern int32_t min_transcript_len;
// Anchor length for junction to be considered well supported
extern uint32_t junction_support;
// -E: how far (bp) from a poorly supported junction to look for a better
// supported one to replace it
extern uint32_t splice_site_error;
// Minimum read coverage per bundle bp for a novel multi-exon transcript (-c;
// reference ones need min_reference_coverage, see drop_weak_predictions)
extern float min_read_coverage;
// Minimum per-base read coverage a single-exon transcript needs to be trusted
// as real rather than noise
extern float min_single_exon_coverage;
// Reads at what distance should be considered part of separate bundles
extern uint32_t min_bundle_distance;
// Distance threshold to detect polymerase run-off artifacts
extern uint32_t runoff_distance;
// Fraction of bundle allowed to be covered by multi-hit reads paper uses 1
extern float multi_read_frac;
// Number of reads needed to support a particular junction
extern int32_t min_junction_coverage;
extern float isoform_fraction;
// --min-borrowed-share (see check_weak_borrow in path_extender.cpp)
extern float min_borrowed_share;
// --keep-contradicted-paths turns it off (see is_path_contradicted in flow.cpp)
extern bool drop_contradicted_paths;
// --keep-shadowed-paths turns it off (see drop_shadowed_paths in flow.cpp)
extern bool drop_shadowed;
extern std::string label;
extern bool QUIET;
extern bool RIBBON_DEBUG;
extern quill::Logger *logger;
// With --demultiplex: the directory where we write one GTF per cell (empty:
// none), next to the output GTF
extern std::string demultiplexed_dir;
// Chromosomes/contigs to exclude (e.g. chrM), from -x
extern std::unordered_set<std::string> excluded_seqs;
extern int32_t n_threads;
extern std::unique_ptr<GFastaDb> g_fasta;
extern std::string guide_gff;
extern std::string out_fname;
extern std::string out_dir;
extern std::string tmp_path;
extern std::string cram_ref;
extern std::string tmp_fname;
// Full command line, echoed in verbose log and output header
extern std::string cmd_line;

namespace ribbon {

  // Formats the current time as MM_DD_HH:MM:SS, used to build a unique
  // temp output filename when writing to stdout (see process_options()).
  char *sprint_time() {
    static char sbuf[32];
    time_t ltime = time(NULL);
    struct tm *t = localtime(&ltime);
    sprintf(sbuf, "%02d_%02d_%02d:%02d:%02d", t->tm_mon + 1, t->tm_mday,
            t->tm_hour, t->tm_min, t->tm_sec);
    return (sbuf);
  }

  // Defines and parses the command line (CLI11) into the global settings,
  // checks that they make sense and returns the input alignment paths in
  // in_files. Exits after printing help, the version or an error.
  void process_options(int32_t argc, char *argv[],
                       std::vector<std::string>& in_files) {
    // Create the command-line parser.
    CLI::App app{"ribbon v" VERSION
                 ": assemble RNA-Seq alignments into potential transcripts",
                 "ribbon"};

    // Values that we check or combine before storing them in the settings.
    std::string guide_path, out_path, label_str, fasta_path, exclude_str,
        cram_ref_path;
    bool no_trim = false;
    bool demultiplex = false;
    bool keep_contradicted_paths = false;
    bool keep_shadowed_paths = false;

    // Declare the command-line options.
    app.add_option("in.bam", in_files,
        "input alignment file(s) (BAM/CRAM, or - for stdin)")->required();
    app.add_flag("--rf", rf_strand, "assume stranded library fr-firststrand");
    app.add_flag("--fr", fr_strand, "assume stranded library fr-secondstrand");
    app.add_flag("--sq3", smartseq3, "SmartSeq3 single-cell mode");
    app.add_flag("-t", no_trim,
                 "disable trimming of predicted transcripts based on coverage");
    app.add_flag("--quiet", QUIET, "turn off verbose (log processing details)");
    app.add_option("-G", guide_path,
                   "reference annotation to use for guiding the assembly "
                   "process (GTF/GFF)")
        ->check(CLI::ExistingFile);
    app.add_option("-o", out_path,
                   "output path/file name for the assembled transcripts GTF "
                   "(default: stdout)");
    app.add_flag("--demultiplex", demultiplex,
                 "also write one GTF per cell in the directory 'demultiplexed' "
                 "next to the output GTF, with the cell's own coverage, FPKM "
                 "and TPM for each transcript");
    app.add_option("-l", label_str,
                   "name prefix for output transcripts (default: RBN)");
    app.add_option("-f", isoform_fraction,
                   "minimum isoform abundance as a fraction of the most "
                   "abundant transcript of its gene (coverage pooled over "
                   "cells); weaker transcripts, including reference ones, "
                   "are dropped (default: 0.01)");
    app.add_flag("--keep-contradicted-paths", keep_contradicted_paths,
                 "keep a path although no read has two of its consecutive "
                 "branch choices together while some reads take another way "
                 "at one of them (default: such paths are dropped)");
    app.add_flag("--keep-shadowed-paths", keep_shadowed_paths,
                 "keep a novel path in a cell although each of the cell's "
                 "reads on it also fits a path with more of the cell's reads "
                 "(default: the cell is taken off such paths)");
    app.add_option("--min-borrowed-share", min_borrowed_share,
                   "in each cell, drop a novel path that has a junction the "
                   "cell has no reads for, chosen from other cells' reads "
                   "among competing junctions with less than this share of "
                   "their reads (default: 0.6; 0 turns it off)");
    app.add_option("-m", min_transcript_len,
                   "minimum assembled transcript length (default: 200)");
    app.add_option("-a", junction_support,
                   "minimum anchor length for junctions (default: 10)");
    app.add_option("-j", min_junction_coverage,
                   "minimum junction coverage (default: 1)");
    app.add_option("-E", splice_site_error,
                   "window around splice sites used to build consensus for "
                   "long reads (default: 25)");
    app.add_option("-c", min_read_coverage,
                   "minimum reads per bp coverage for a novel multi-exon "
                   "transcript; reference ones need 1 (default: 3)");
    app.add_option("-s", min_single_exon_coverage,
                   "minimum reads per bp coverage to consider for "
                   "single-exon transcript (default: 4.75)");
    app.add_option("-g", min_bundle_distance,
                   "maximum gap allowed between read mappings (default: 50)");
    app.add_option("-M", multi_read_frac,
                   "fraction of bundle allowed to be covered by multi-hit "
                   "reads (default: 1)");
    app.add_option("-p", n_threads,
                   "number of threads (CPUs) to use (default: 1)");
    app.add_option("-x", exclude_str,
                   "do not assemble any transcripts on the given reference "
                   "sequence(s), comma-separated");
    app.add_option("-S,--rseq", fasta_path,
                   "reference genome sequence (FASTA)");
    app.add_option("--ref,--cram-ref", cram_ref_path,
                   "reference genome FASTA file for CRAM input");
    auto version_callback = [](int32_t) {
      std::cout << VERSION << "\n";
      std::exit(EXIT_SUCCESS);
    };
    app.add_flag_function("--version", version_callback,
                          "print just the version at stdout and exit");

    // Keep the command line to echo it in the log and in the output.
    for (int32_t i = 0; i < argc; i++) {
      cmd_line += argv[i];
      cmd_line += (i == argc - 1) ? "" : " ";
    }

    // Parse the command line; exit after help, the version or a parse error.
    try {
      app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
      std::exit(app.exit(e));
    }

    if (no_trim) trim_by_coverage = false;
    if (keep_contradicted_paths) drop_contradicted_paths = false;
    if (keep_shadowed_paths) drop_shadowed = false;

    // Check the combination of options and log the command line.
    if (fr_strand && rf_strand) {
      fprintf(stderr, "Error: --fr and --rf options are incompatible.\n");
      exit(1);
    }

    if (!QUIET || RIBBON_DEBUG) {
      LOG_INFO(logger, "starting ribbon version {}, command line: {}", VERSION,
               cmd_line.c_str());
    }

    // Check the option values and store them in the global settings.
    if (min_transcript_len < 30) {
      fprintf(stderr, "Error: invalid -m value, must be >=30)\n");
      exit(1);
    }
    if (!fasta_path.empty())
      g_fasta = std::make_unique<GFastaDb>(fasta_path.c_str());
    if (!cram_ref_path.empty()) cram_ref = cram_ref_path.c_str();
    if (!exclude_str.empty()) {
      // Split by comma and populate excludeGSeqs
      size_t token_start = exclude_str.find_first_not_of(" ,\t");
      while (token_start != std::string::npos) {
        size_t token_end = exclude_str.find_first_of(" ,\t", token_start);
        excluded_seqs.insert(exclude_str.substr(token_start,
                                                token_end - token_start));
        token_start = exclude_str.find_first_not_of(" ,\t", token_end);
      }
    }
    if (n_threads <= 0) n_threads = 1;
    if (min_read_coverage < 0.001) {
      fprintf(stderr, "Error: invalid -c value, must be >=0.001)\n");
      exit(1);
    }
    if (min_single_exon_coverage < 0.001) {
      fprintf(stderr, "Error: invalid -s value, must be >=0.001)\n");
      exit(1);
    }
    if (min_bundle_distance > runoff_distance)
      runoff_distance = min_bundle_distance;
    // Increase gap length
    if (min_bundle_distance < 100) min_bundle_distance = 100;
    if (!label_str.empty()) label = label_str.c_str();
    if (min_borrowed_share < 0 || min_borrowed_share > 1) {
      fprintf(stderr, "Error: invalid --min-borrowed-share value, must be "
              "between 0 and 1\n");
      exit(1);
    }
    if (isoform_fraction >= 1) {
      fprintf(stderr,
              "Miminum isoform fraction (-f coefficient: %f) needs to be "
              "less than 1\n",
              isoform_fraction);
      exit(1);
    }
    if (!guide_path.empty()) {
      guide_gff = guide_path.c_str();
      guided = true;
    }
    tmp_fname = out_path.c_str();

    // Deferred creation of output path
    out_fname = "stdout";
    out_dir = "./";
    if (!tmp_fname.empty() && tmp_fname != "-") {
      if (tmp_fname[0] == '.' && tmp_fname[1] == '/')
        tmp_fname.erase(0, 2);
      out_fname = tmp_fname;
      size_t pidx = out_fname.rfind('/');
      if (pidx != std::string::npos) { // path given
        out_dir = out_fname.substr(0, pidx + 1);
        tmp_fname = out_fname.substr(pidx + 1);
      }
    } else { // stdout
      tmp_fname = out_fname;
      char *stime = sprint_time();
      std::replace(tmp_fname.begin(), tmp_fname.end(), ':', '-');
      tmp_fname += '.';
      tmp_fname += stime;
    }

    // Create the output directory if it does not exist.
    if (out_dir != "./") {
      if (!std::filesystem::exists(out_dir.c_str())) {
        // Directory does not exist, create it
        if (Gmkdir(out_dir.c_str())
            && !std::filesystem::exists(out_dir.c_str())) {
          fprintf(stderr, "Error: cannot create directory %s!\n",
                  out_dir.c_str());
          exit(1);
        }
      }
    }

    // Create the directory of the per-cell GTFs now, so that a bad path fails
    // before the assembly runs. Files of an earlier run stay in it.
    if (demultiplex) {
      demultiplexed_dir = out_dir + "demultiplexed";
      std::error_code create_error;
      std::filesystem::create_directories(demultiplexed_dir, create_error);
      if (create_error) {
        fprintf(stderr, "Error creating directory %s: %s\n",
                demultiplexed_dir.c_str(), create_error.message().c_str());
        exit(1);
      }
    }

    // Create the temporary directory for the results of the worker threads.
    { // prepare temp path
      std::string templ = out_dir;
      while (!templ.empty() && templ.back() == '/') templ.pop_back();
      templ += "/tmp_XXXXXX";
      if (mkdtemp(&templ[0]) == NULL) {
        fprintf(stderr, "Error creating temp dir %s!(%s)\n", templ.c_str(),
                strerror(errno));
        exit(1);
      }
      tmp_path = templ;
      tmp_path += '/';
    }

    // Open the temporary results file.
    tmp_fname = tmp_path + tmp_fname;
    tmp_fname += ".tmp";
    f_out = fopen(tmp_fname.c_str(), "w");
    if (f_out == NULL) {
      fprintf(stderr, "Error creating output file %s\n", tmp_fname.c_str());
      exit(1);
    }
  }

} // namespace ribbon
