#!/usr/bin/env bash
# THE CUT-BATCH PER-TERRITORY FAIL-CLOSED ISOLATION - the regression gate:
# the 1-territory L1 batch re-run on the patched binary must reproduce the
# committed cutbatch-gate-L1 receipt BYTE-IDENTICALLY (census + cache +
# reads fasta + the cut-stats content fields) - the emissions-loop change
# must be inert on the success path.
set -uo pipefail
RV=/home/erikg/hprcv2/wgs35x-fleet/runtime-verification
IMPG=/home/erikg/impg-genome-inference/target/release/impg
BIN=/home/erikg/impg-genome-inference/target/release/examples/sample_projection
G=/home/erikg/hprcv2/generic-sampling/L1
SYNG_GLOBAL=/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng
AGC=/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc
SHARDS=/home/erikg/hprcv2/hg002-projection/full
COMMITTED=/home/erikg/hprcv2/hg002-projection/cutbatch-gate-L1
S=$RV/regression-L1
AXIS_PATH=HG002#1#chr1:4632437-4883593
mkdir -p "$S/results"
[ -f "$S/syng.done" ] || {
  "$IMPG" syng -f "$G/locus.fa" --syncmer-length 63 --smer-length 8 \
    --syncmer-seed 7 --position-sample-rate 64 -o "$S/locus-pos64.syng" \
    > "$S/syng.log" 2>&1 && touch "$S/syng.done" || { echo "syng FAILED"; exit 1; }
}
[ -f "$S/beds.done" ] || {
  printf '%s\n' "$AXIS_PATH" > "$S/starting-sequences.txt"
  "$IMPG" partition -a "$S/locus-pos64.syng" -w 10000 -d 1000 \
    --starting-sequences-file "$S/starting-sequences.txt" \
    --selection-mode haplotype --syng-min-chain-anchors 5 --syng-min-chain-fraction 0.5 \
    --no-rehome-singletons --sequence-files "$G/locus.fa" \
    -o bed --separate-files --output-folder "$S/results" -t 16 \
    > "$S/beds.log" 2>&1 && touch "$S/beds.done" || { echo "beds FAILED"; exit 1; }
}
cat > "$S/manifest.jsonl" <<EOF
{"name": "L1", "panel": "$S/locus-pos64.syng", "axis": "$G/reference-axis.json", "bed_directory": "$S/results", "component": "$AXIS_PATH", "sources": "$G/locus.fa", "out_dir": "$S/cut"}
EOF
rm -rf "$S/cut"; mkdir -p "$S/cut"
start=$(date +%s)
RAYON_NUM_THREADS=96 IMPG_CUT_LOCAL_DERIVE=1 nice -n 10 "$BIN" cut-batch \
  --global-panel "$SYNG_GLOBAL" --agc "$AGC" --shards-dir "$SHARDS" \
  --manifest "$S/manifest.jsonl" > "$S/batch.log" 2> "$S/batch.err"
st=$?
echo "$st" > "$S/batch.exit"
echo "$(( $(date +%s) - start ))" > "$S/batch.wall"
[ "$st" -ne 0 ] && { echo "BATCH FAILED"; tail -5 "$S/batch.err"; exit 1; }
fail=0
for f in cut-census.jsonl cut-cache.bin cut-reads.fasta; do
  if cmp -s "$S/cut/$f" "$COMMITTED/$f"; then
    echo "IDENTICAL: $f"
  else
    echo "DIFF: $f"; fail=1
  fi
done
python3 - "$S/cut/cut-stats.json" "$COMMITTED/cut-stats.json" >> "$S/batch.log" 2>&1 <<'EOF'
import json, sys
a, b = (json.load(open(p)) for p in sys.argv[1:3])
diff = {k: (a.get(k), b.get(k)) for k in set(a) | set(b) if a.get(k) != b.get(k)}
print("cut-stats content diff:", diff if diff else "IDENTICAL")
sys.exit(1 if diff else 0)
EOF
[ $? -ne 0 ] && fail=1
[ "$fail" -eq 0 ] && { touch "$S/gate.done"; echo "=== REGRESSION GATE PASS: the patched cut-batch byte-identical on the success path"; } || echo "=== REGRESSION GATE FAIL"
exit "$fail"
