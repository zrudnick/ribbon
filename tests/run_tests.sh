#!/usr/bin/env bash
# Rebuild ribbon, then run the single-cell regression tests against the baselines in
# ribbon_data/test/scell_baseline (3 smallest cells from ribbon_data/benchmark/simulation/bam):
#   - de novo and guided (-G Ensembl 104): the output body (everything but the '#' header
#     line) must be byte-identical to the baseline
# Runs use -p 1 and -c 1: the default -c 3 is a cutoff on coverage pooled over the cells and
# would remove most transcripts of 3 cells (4891 -> 854 de novo), leaving little to test.
# With -p > 1 the transcripts are the same, but gene/transcript id numbering,
# line order and the last digits of FPKM/TPM vary with thread scheduling.
# usage: tests/run_tests.sh [quick]     quick = skip the guided test (~35s)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RIBBON_DIR="$(dirname "$HERE")"
DATA_DIR="${RIBBON_DATA:-$RIBBON_DIR/../ribbon_data}"
OUT="${TMPDIR:-/tmp}/ribbon_tests"
ninja -C "$RIBBON_DIR/build" >"$OUT.build.log" 2>&1 || { grep -E 'error' "$OUT.build.log" | head; echo BUILD FAILED; exit 1; }
mkdir -p "$OUT" && rm -f "$OUT"/*.gtf
cd "$DATA_DIR" || exit 1
RIBBON="$RIBBON_DIR/build/ribbon"
BAMS="benchmark/simulation/bam/155.bam benchmark/simulation/bam/63.bam benchmark/simulation/bam/96.bam"
ok=1
$RIBBON $BAMS -c 1 -p 1 -o "$OUT/denovo.gtf" >/dev/null 2>&1
if diff <(grep -v "^#" "$OUT/denovo.gtf") <(grep -v "^#" test/scell_baseline/scell_denovo.gtf) >/dev/null; then echo "scell denovo OK"; else echo "scell denovo DIFFERS"; ok=0; fi
if [ "$1" != quick ]; then
  $RIBBON $BAMS -G data/ref/Homo_sapiens.GRCh38.104.gtf -c 1 -p 1 -o "$OUT/guided.gtf" >/dev/null 2>&1
  if diff <(grep -v "^#" "$OUT/guided.gtf") <(grep -v "^#" test/scell_baseline/scell_guided.gtf) >/dev/null; then echo "scell guided OK"; else echo "scell guided DIFFERS"; ok=0; fi
fi
[ $ok = 1 ] && echo ALL PASS || { echo FAIL; exit 1; }
