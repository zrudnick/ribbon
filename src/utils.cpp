#include "utils.h"
#include <cctype>
#include <string>

// General helpers used by many files. The templates are in utils.h.

// --mito: the name of the mitochondrial chromosome, if not a usual one
extern std::string mito_name;

namespace ribbon {

  // Returns the root of x in a union-find array (parent[x] == x for a root).
  int32_t find_root(const std::vector<int32_t>& parent, int32_t x) {
    while (parent[x] != x) x = parent[x];
    return x;
  }

  // Erases the first occurrence of value from v; false if there is none.
  bool erase_value(std::vector<int32_t>& v, int32_t value) {
    auto it = std::find(v.begin(), v.end(), value);
    if (it == v.end()) return false;
    v.erase(it);
    return true;
  }

  // Inserts value into the sorted vector v, unless it is already there.
  void insert_sorted_unique(std::vector<int32_t>& v, int32_t value) {
    auto it = std::lower_bound(v.begin(), v.end(), value);
    if (it == v.end() || *it != value) v.insert(it, value);
  }


  // True if sequence_name is the mitochondrial chromosome: --mito, or one of
  // the usual names (Ensembl MT, UCSC chrM, Ensembl yeast Mito, SGD chrmt,
  // FlyBase mitochondrion_genome, WormBase MtDNA, Ensembl Plants Mt, ...), case
  // ignored, or the RefSeq accession (any version) of the mitochondrial genome
  // of a common model organism.
  bool is_mitochondrial(const std::string& sequence_name) {
    if (!mito_name.empty() && sequence_name == mito_name) return true;
    std::string name;
    for (char c : sequence_name) name += (char)std::tolower((unsigned char)c);
    for (const char *usual : {"mt", "m", "chrm", "chrmt", "mito", "chrmito",
                              "mitochondrion", "mitochondrion_genome",
                              "mtdna", "chrmtdna"})
      if (name == usual) return true;
    // Human, mouse, rat, zebrafish, fruit fly, C. elegans, yeast
    std::string accession = name.substr(0, name.find('.'));
    for (const char *refseq : {"nc_012920", "nc_005089", "nc_001665",
                               "nc_002333", "nc_024511", "nc_001328",
                               "nc_001224"})
      if (accession == refseq) return true;
    return false;
  }
} // namespace ribbon
