#include "out.h"
#include "predictions.h"
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include "ribbon.h"
#include "GThreads.h"

// Writes the predictions of a finished bundle (print_results) and completes
// the output file (write_final_output). A transcript appears once, with its
// coverage in every cell (`cell "name:coverage"`) and the number of cells.
// With --demultiplex it also appears in the GTF of every cell that has it,
// with the cell's own coverage, FPKM and TPM.
// Which predictions are kept was already decided (see predictions.cpp).

// Import globals from main program:

extern FILE *f_out;
extern std::string label;
extern std::string tmp_fname; // temporary results file written by worker threads
extern double total_frag_len;
extern double total_cov_sum;
extern std::map<std::string, double> frag_len_by_cell;
extern std::map<std::string, double> cov_sum_by_cell;
extern std::string demultiplexed_dir;

namespace ribbon {

  // A predicted gene: the union of the exons of its transcripts.
  struct CGene : public GSeg {
    char strand;
    double cov; // this is the actual gene coverage
    // This is a sum of transcripts coverages -> this is what we need for FPKM
    // and TPM estimations
    double cov_sum;
    // All possible exons in gene (those are bnodes in bundle)
    std::vector<GSeg> exons;
    CGene(int32_t g_start = 0, int32_t g_end = 0, char _strand = '.')
        : GSeg(g_start, g_end), strand(_strand), exons() {
      cov = 0;
      cov_sum = 0;
    }
  };

  // Reads the cells of a transcript from its first line: the `cell "name:cov"`
  // attributes, in the order they are printed.
  std::vector<std::pair<std::string, double>> parse_cell_coverages(
      const std::string& transcript_line) {
    std::vector<std::pair<std::string, double>> cells;
    const std::string attribute = " cell \"";
    size_t pos = 0;
    while ((pos = transcript_line.find(attribute, pos)) != std::string::npos) {
      pos += attribute.size();
      size_t end_quote = transcript_line.find('"', pos);
      if (end_quote == std::string::npos) break;
      size_t colon = transcript_line.rfind(':', end_quote);
      if (colon == std::string::npos || colon < pos) break;
      cells.emplace_back(transcript_line.substr(pos, colon - pos),
                         atof(transcript_line.substr(colon + 1,
                                                     end_quote - colon - 1).c_str()));
      pos = end_quote;
    }
    return cells;
  }

  // The GTF files of the cells (--demultiplex). A cell's file is created the
  // first time we write to it. To stay under the limit of open files, we close
  // them all when too many are open and append when we reopen one.
  struct CellGtfWriter {
    std::string dir;
    std::map<std::string, FILE *> open_files;
    std::set<std::string> created;
    size_t max_open_files = 500;

    explicit CellGtfWriter(const std::string& _dir) : dir(_dir) {}
    ~CellGtfWriter() { close_all(); }

    std::string path_of(const std::string& cell) const {
      std::string name = cell;
      for (char& ch : name)
        if (ch == '/') ch = '_';
      return dir + "/" + name + ".gtf";
    }

    void close_all() {
      for (auto& open_file : open_files) fclose(open_file.second);
      open_files.clear();
    }

    // The open file of the cell, creating or reopening it if needed.
    FILE *file_of(const std::string& cell) {
      auto it = open_files.find(cell);
      if (it != open_files.end()) return it->second;
      if (open_files.size() >= max_open_files) close_all();
      std::string path = path_of(cell);
      FILE *file = fopen(path.c_str(), created.count(cell) ? "a" : "w");
      if (!file) {
        perror(path.c_str());
        exit(1);
      }
      created.insert(cell);
      open_files[cell] = file;
      return file;
    }

    // Creates the empty file of every cell that has no transcript.
    void create_missing_files(const std::map<std::string, double>& all_cells) {
      for (const auto& cell : all_cells)
        if (!created.count(cell.first)) file_of(cell.first);
    }
  };

  // Writes a transcript (its first line and its exon lines) to the GTF of
  // each of its cells. The first line gets the cell's own coverage, FPKM and
  // TPM. The TPM of a cell is normalized over the transcripts that we report
  // for that cell, so it differs from the pooled TPM, whose total also counts
  // the coverage that no transcript took.
  void write_transcript_to_cell_gtfs(CellGtfWriter& writer,
                                     const std::vector<std::string>& lines) {
    std::vector<std::pair<std::string, double>> cells
        = parse_cell_coverages(lines[0]);
    // Everything of the first line before the pooled coverage
    size_t cov_start = lines[0].find(" cov \"");
    if (cells.empty() || cov_start == std::string::npos) return;
    std::string line_start = lines[0].substr(0, cov_start);

    for (const auto& cell : cells) {
      double cell_cov = cell.second;
      auto cell_frag_len = frag_len_by_cell.find(cell.first);
      auto cell_cov_sum = cov_sum_by_cell.find(cell.first);
      double fpkm = 0;
      double tpm = 0;
      if (cell_frag_len != frag_len_by_cell.end() && cell_frag_len->second > 0.001)
        fpkm = cell_cov * 1000000000 / cell_frag_len->second;
      if (cell_cov_sum != cov_sum_by_cell.end() && cell_cov_sum->second > 0.00001)
        tpm = cell_cov * 1000000 / cell_cov_sum->second;

      FILE *file = writer.file_of(cell.first);
      fprintf(file, "%s cov \"%.6f\"; FPKM \"%.6f\"; TPM \"%.6f\";\n",
              line_start.c_str(), cell_cov, fpkm, tpm);
      for (size_t i = 1; i < lines.size(); i++)
        fprintf(file, "%s\n", lines[i].c_str());
    }
  }

  // Final pass over the temporary results file: copies each transcript's lines
  // to f_out, adding FPKM and TPM to the first one (and writes it to the GTF of
  // each of its cells with --demultiplex), then closes f_out and deletes the
  // file.
  void write_final_output() {
    std::unique_ptr<CellGtfWriter> cell_writer;
    if (!demultiplexed_dir.empty())
      cell_writer = std::make_unique<CellGtfWriter>(demultiplexed_dir);
    // Read the temporary file that the worker threads wrote.
    std::ifstream tmp_in(tmp_fname.c_str());
    if (tmp_in) {
      std::string line;
      // Reads the next line (lines end at \n, \r or \r\n); returns false at end
      // of file
      auto next_line = [&tmp_in, &line]() {
        if (!std::getline(tmp_in, line)) return false;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        return true;
      };
      int32_t line_count;
      int32_t is_transcript;
      int32_t transcript_len;
      double cov; // increased precision for downstream calculations
      double fpkm;
      double tpm;
      int32_t transcript_no;

      // Compute FPKM and TPM of each transcript from its coverage and copy its
      // lines, adding FPKM and TPM to the first line.
      while (next_line()) {
        if (sscanf(line.c_str(), "%d %d %d %d %lg", &is_transcript,
                   &line_count, &transcript_len, &transcript_no, &cov) != 5) {
          fprintf(stderr, "Malformed line in temporary file %s: %s\n",
                  tmp_fname.c_str(), line.c_str());
          exit(1);
        }
        if (cov < 0) cov = 0;
        if (total_frag_len > 0.001) fpkm = cov * 1000000000 / total_frag_len;
        else fpkm = 0.0;
        if (total_cov_sum > 0.00001) tpm = cov * 1000000 / total_cov_sum;
        else tpm = 0.0;
        if (is_transcript) { // this is a transcript
          std::vector<std::string> transcript_lines;
          for (int32_t i = 0; i < line_count; i++) {
            if (!next_line()) {
              fprintf(stderr,
                      "Temporary file %s ends in the middle of a transcript\n",
                      tmp_fname.c_str());
              exit(1);
            }
            if (cell_writer) transcript_lines.push_back(line);
            if (!i) {
              fprintf(f_out, "%s", line.c_str());
              fprintf(f_out, " FPKM \"%.6f\";", fpkm);
              fprintf(f_out, " TPM \"%.6f\";", tpm);
              fprintf(f_out, "\n");
            } else fprintf(f_out, "%s\n", line.c_str());
          }
          if (cell_writer)
            write_transcript_to_cell_gtfs(*cell_writer, transcript_lines);
        }
      }

      // Every cell gets a file, even one without transcripts.
      if (cell_writer) cell_writer->create_missing_files(frag_len_by_cell);

      // Close and delete the temporary file.
      fclose(f_out);
      tmp_in.close();
      int32_t rv = remove(tmp_fname.c_str());
      if (rv) {
        perror(tmp_fname.c_str());
      }
    } else {
      fclose(f_out);
      fprintf(stderr, "No temporary file %s present!\n", tmp_fname.c_str());
      exit(1);
    }
  }

  // Merges the exons of a transcript into the exon list of its gene: exons that
  // overlap nothing are inserted, the ones that overlap are extended or
  // merged.
  void merge_exons(CGene& gene, const std::vector<GSeg>& exons) {
    int32_t ig = 0;
    int32_t ie = 0;
    while (ie < count_of(exons)) {
      // Insert an exon that lies before the next gene exon.
      if (ig == count_of(gene.exons) || exons[ie].end < gene.exons[ig].start) {
        gene.exons.insert(gene.exons.begin() + ig, exons[ie]);
        ie++;
        ig++;
        continue;
      }

      // Skip the gene exons that end before this exon, then extend or merge the
      // gene exons that it overlaps.
      while (ig < count_of(gene.exons) && exons[ie].start > gene.exons[ig].end)
        ig++;
      // Here exons[ie]->start<=gene.exons[ig].end and
      // exons[ie]->end>=gene.exons[ig].start
      if (ig < count_of(gene.exons)) {
        if (exons[ie].start <= gene.exons[ig].start)
          gene.exons[ig].start = exons[ie].start;
        if (exons[ie].end >= gene.exons[ig].end) {
          gene.exons[ig].end = exons[ie].end;
          ig++;
          while (ig < count_of(gene.exons)
                 && exons[ie].end >= gene.exons[ig].start) {
            if (gene.exons[ig].end > exons[ie].end)
              gene.exons[ig - 1].end = gene.exons[ig].end;
            gene.exons.erase(gene.exons.begin() + ig);
          }
        }
        ie++;
      }
    }
  }

  // Numbers the genes: all predictions of a gene get the same gene number.
  // Returns the gene number of each prediction (-1 for the unflagged ones) and
  // adds an entry to transcripts_per_gene for every gene.
  std::vector<int32_t> assign_gene_numbers(
      PredictionList& predictions, std::vector<int32_t>& gene_parent,
      std::vector<int32_t>& transcripts_per_gene) {
    int32_t prediction_count = count_of(predictions);
    std::vector<int32_t> gene_of_prediction(prediction_count, -1);
    for (int32_t p = 0; p < prediction_count; p++)
      // No gene was assigned to this prediction
      if (predictions[p]->flag && gene_of_prediction[p] == -1) {
        int32_t root_n = find_gene_root(gene_parent, p);
        gene_of_prediction[p] = count_of(transcripts_per_gene);
        transcripts_per_gene.push_back(0);
        uint32_t current_end = predictions[p]->end;
        int32_t m = p + 1;
        while (m < prediction_count && predictions[m]->start <= current_end) {
          if (predictions[m]->flag) {
            int32_t root_m = find_gene_root(gene_parent, m);
            if (root_n == root_m) {
              if (gene_of_prediction[m] > -1) {
                fprintf(stderr,
                        "Overlapping predictions with different gene "
                        "number!\n");
                exit(1);
              }
              gene_of_prediction[m] = gene_of_prediction[p];
            }
            if (predictions[m]->end > current_end) current_end = predictions[m]->end;
          }
          m++;
        }
      }

    return gene_of_prediction;
  }

  // Writes a transcript and its exons in GTF format, with its coverage and the
  // coverage of every cell where we found it. The line before them is read back
  // by write_final_output.
  void print_transcript(const CPrediction& p, const std::string& gene_id,
                        const std::string& transcript_id,
                        const std::string& ref_name,
                        const std::vector<const std::string *>& cell_names) {
    fprintf(f_out, "1 %d %d 0 %.6f\n", count_of(p.exons) + 1, p.tlen, p.cov);

    fprintf(f_out,
            "%s\tribbon\ttranscript\t%d\t%d\t1000\t%c\t.\t"
            "gene_id \"%s\"; transcript_id \"%s\"; ",
            ref_name.c_str(), p.start, p.end, p.strand, gene_id.c_str(),
            transcript_id.c_str());

    fprintf(f_out, "cov \"%.6f\";", p.cov);

    // Print cell coverages
    int32_t cells_with_coverage = 0;
    for (int32_t cell = 0; cell < count_of(cell_names); ++cell) {
      // If there is coverage in this cell (a cell with zero or negative
      // coverage does not have the transcript)
      if (p.cell_cov[cell] > 0) {
        fprintf(f_out, " cell \"%s:%f\";", cell_names[cell]->c_str(),
                p.cell_cov[cell]);
        cells_with_coverage++;
      }
    }
    if (cells_with_coverage) {
      fprintf(f_out, " no_cells \"%d\";", cells_with_coverage);
    }

    fprintf(f_out, "\n");
    for (int32_t j = 0; j < count_of(p.exons); j++) {
      fprintf(f_out,
              "%s\tribbon\texon\t%d\t%d\t1000\t%c\t.\t"
              "gene_id \"%s\"; transcript_id \"%s\"; exon_number \"%d\"; ",
              ref_name.c_str(), p.exons[j].start, p.exons[j].end, p.strand,
              gene_id.c_str(), transcript_id.c_str(), j + 1);
      fprintf(f_out, "\n");
    }
  }

  // Adds a printed transcript to its predicted gene (used for the
  // leftover-coverage computation): creates the gene the first time and
  // extends it afterwards.
  void add_to_predicted_gene(std::vector<CGene>& predicted_genes,
                             int32_t gene_idx, const CPrediction& p) {
    // We did not see this gene before
    if (gene_idx >= count_of(predicted_genes)) {
      CGene g(p.start, p.end, p.strand);
      for (int32_t i = 0; i < count_of(p.exons); i++) {
        g.exons.push_back(p.exons[i]);
      }
      g.cov += p.cov * p.tlen;
      g.cov_sum += p.cov;
      predicted_genes.push_back(g);
    } else { // we've seen this gene before
      if (p.start < predicted_genes[gene_idx].start)
        predicted_genes[gene_idx].start = p.start;
      if (p.end > predicted_genes[gene_idx].end)
        predicted_genes[gene_idx].end = p.end;
      merge_exons(predicted_genes[gene_idx], p.exons);
      predicted_genes[gene_idx].cov += p.cov * p.tlen;
      predicted_genes[gene_idx].cov_sum += p.cov;
    }
  }

  // Computes the coverage that is left for each gene after its transcripts,
  // and adds it to sum_cov for the TPM normalization.
  void add_leftover_gene_coverage(std::vector<CGene>& predicted_genes,
                                  BundleData *bundle) {
    BpCoverage& bp_coverage = bundle->bp_coverage;
    for (int32_t i = 0; i < count_of(predicted_genes); i++) {
      float cov = 0;
      char s = '.'; // strand of gene
      if (predicted_genes[i].strand == '+') s = '+';
      else if (predicted_genes[i].strand == '-') s = '-';
      int32_t gene_len = 0;
      // Evaluate unused coverage
      for (int32_t j = 0; j < count_of(predicted_genes[i].exons); j++) {
        gene_len += predicted_genes[i].exons[j].len();
        int32_t start = (int32_t)predicted_genes[i].exons[j].start - bundle->start;
        int32_t end = (int32_t)predicted_genes[i].exons[j].end - bundle->start + 1;
        // Predgene start and end might have been also adjusted to reflect the
        // annotation
        if (start < 0) start = 0;
        if (end >= count_of(bp_coverage['.']))
          end = count_of(bp_coverage['.']) - 1;

        // Cummulative bp_coverage
        PerStrand<float> local_coverage;
        for (char strand : {'-', '.', '+'}) {
          if (end)
            local_coverage[strand] = get_region_coverage(strand, start,
                end - 1, bp_coverage);
        }

        float minus_cov = local_coverage['-'];
        float all_cov = local_coverage['.'];
        float plus_cov = local_coverage['+'];
        switch (s) {
          case '-':
            if (plus_cov)
              cov += minus_cov + (all_cov - minus_cov - plus_cov)
                                 * minus_cov / (minus_cov + plus_cov);
            else cov += all_cov;
            break;
          case '.':
            cov += all_cov - plus_cov - minus_cov;
            break;
          case '+':
            if (minus_cov)
              cov += plus_cov + (all_cov - minus_cov - plus_cov)
                                * plus_cov / (minus_cov + plus_cov);
            else cov += all_cov;
            break;
        }
      }
      // THIS is the read coverage that is left for the genes after all the
      // predictions were taken into account
      predicted_genes[i].cov = cov - predicted_genes[i].cov;
      if (predicted_genes[i].cov > epsilon)
        predicted_genes[i].cov_sum += predicted_genes[i].cov / gene_len;
      bundle->sum_cov += predicted_genes[i].cov_sum;
    }
  }

  // Prints the predictions of one bundle that finalize_predictions kept:
  // groups them into genes and writes them. Returns the updated gene count.
  int32_t print_results(BundleData *bundle, int32_t gene_count,
                    const std::string& ref_name) {
    // Get the predictions of this bundle (sorted by position) and the cells'
    // names, indexed by cell number (cells are numbered in the order they were
    // first seen).
    PredictionList& predictions = bundle->predictions;
    int32_t prediction_count = count_of(predictions);
    int32_t cell_count = count_of(bundle->cell_name);
    std::vector<const std::string *> cell_names(cell_count);
    for (const auto& cell : bundle->cell_name)
      cell_names[cell.second] = &cell.first;

    // Group the remaining predictions into genes.
    std::vector<int32_t> gene_parent
        = group_predictions_into_genes(predictions);
    // For each gene remembers how many transcripts were printed
    std::vector<int32_t> transcripts_per_gene;
    std::vector<int32_t> gene_of_prediction
        = assign_gene_numbers(predictions, gene_parent, transcripts_per_gene);

    // Print each remaining transcript and add it to its predicted gene.
    std::vector<CGene> predicted_genes;
    for (int32_t p = 0; p < prediction_count; p++)
      if (predictions[p]->flag) {
        transcripts_per_gene[gene_of_prediction[p]]++;
        predictions[p]->gene_no = gene_of_prediction[p] + gene_count + 1;

        std::string gene_id = label + '.'
            + std::to_string(predictions[p]->gene_no);
        std::string transcript_id =
            gene_id + '.'
            + std::to_string(transcripts_per_gene[gene_of_prediction[p]]);

        print_transcript(*predictions[p], gene_id, transcript_id, ref_name,
                         cell_names);
        // Add the coverage of each cell to the sum for its TPM normalization.
        if (count_of(bundle->cell_sum_cov) < cell_count)
          bundle->cell_sum_cov.resize(cell_count, 0.0);
        for (int32_t c = 0; c < count_of(predictions[p]->cell_cov); c++)
          if (predictions[p]->cell_cov[c] > 0)
            bundle->cell_sum_cov[c] += predictions[p]->cell_cov[c];
        add_to_predicted_gene(predicted_genes, gene_of_prediction[p],
                              *predictions[p]);
      }

    // Update the gene count.
    gene_count += count_of(transcripts_per_gene);

    add_leftover_gene_coverage(predicted_genes, bundle);

    return (gene_count);
  }

} // namespace ribbon

