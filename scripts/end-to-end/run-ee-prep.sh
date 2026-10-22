#!/usr/bin/env bash
# THE END-TO-END RECEIPT, STEP 1 (per territory) - THE SUBSTRATE CHAIN.
# The export-profile prep recipe VERBATIM (the committed runtime-
# verification chain: extract -> syng -> partition beds -> export
# graphs -> panel walks -> locality maps, the c-rider K=63 export
# profile - THE RUNTIME PROFILE, honestly stated), re-derived into the
# receipt's own workspace from the campaign's surviving smalls (the
# definitional per-territory inputs live in the RV prep dirs; the
# heavy substrate was reaped at the campaign close per the stream
# law). The SUBSTRATE INSTRUMENTS stay pinned to the committed
# target/release binaries (unchanged since Oct 19) - the era split,
# stated. Marker-gated per step; wall + core-s receipts at write
# time; resumable. usage: run-ee-prep.sh <TID>  (env: THREADS=4)
source /home/erikg/impg-genome-inference/scripts/end-to-end/ee-env.sh
T=$1
P=$EE/prep/$T
THREADS=${THREADS:-4}
mkdir -p "$P"
[ -f "$P/substrate.done" ] && exit 0
ee_budget_check || exit 1

# seed the definitional smalls from the campaign's surviving prep dir
for f in axis-path.txt locus.bed locus.members.bed locus.refined.bed extract.specs.txt split.info.json; do
  [ -f "$P/$f" ] || cp "$RV/prep/$T/$f" "$P/$f" 2>/dev/null || true
done
[ -f "$P/axis-path.txt" ] || { echo "$T no-axis (not accepted)" > "$P/substrate.skipped"; exit 0; }
AXIS_PATH=$(cat "$P/axis-path.txt")

# 1. extract (the AGC getctg chain, one xargs per territory)
if [ ! -f "$P/extract.done" ]; then
  ee_step "$P" extract python3 "$RV/rv-extract.py" "$EE" "$T" || exit 1
fi

# 2. the local syng (k63 s8 seed7 pos64, the committed recipe)
ee_step "$P" syng "$IMPG_SUB" syng -f "$P/locus.fa" --syncmer-length 63 \
  --smer-length 8 --syncmer-seed 7 --position-sample-rate 64 \
  -o "$P/locus-pos64.syng" || exit 1

# 3. the partition beds (w10k d1k haplotype, the committed recipe)
BEDS="$P/partition-w10k-d1k/results"
if [ ! -f "$P/beds.done" ]; then
  mkdir -p "$P/partition-w10k-d1k" "$BEDS"
  printf '%s\n' "$AXIS_PATH" > "$P/partition-w10k-d1k/starting-sequences.txt"
  ee_step "$P" beds nice -n 10 "$IMPG_SUB" partition \
    -a "$P/locus-pos64.syng" -w 10000 -d 1000 \
    --starting-sequences-file "$P/partition-w10k-d1k/starting-sequences.txt" \
    --selection-mode haplotype \
    --syng-min-chain-anchors 5 --syng-min-chain-fraction 0.5 \
    --no-rehome-singletons \
    --sequence-files "$P/locus.fa" \
    -o bed --separate-files --output-folder "$BEDS" \
    -t "$THREADS" || exit 1
fi
NP=$(ls "$BEDS" 2>/dev/null | grep -c '\.bed$')
[ "${NP:-0}" -gt 0 ] || { echo "empty partition at $T" > "$P/beds.empty"; exit 1; }

# 4. the export graphs (the syng-native range GFAs)
if [ ! -f "$P/export.done" ]; then
  mkdir -p "$P/export-graphs"
  IDS=$(seq 0 $((NP-1)) | tr '\n' ' ')
  ee_step "$P" export "$EX_SUB/partition_graph_export" export "$P/locus-pos64.syng" \
    "$P/locus.fa" "$BEDS" "$P/export-graphs" $IDS || exit 1
fi

# 5. the panel walks
if [ ! -f "$P/walks.done" ]; then
  ee_step "$P" walkreq python3 "$RV/rv-mk-walk-requests.py" "$BEDS" \
    "$P/locus-pos64.syng.names" "$P/walk-requests.tsv" || exit 1
  ee_step_out "$P" walks "$P/export-walks.jsonl" nice -n 10 \
    "$EX_SUB/syng_homology_dump" walks "$P/locus-pos64.syng" \
    "$P/walk-requests.tsv" || exit 1
fi

# 6. the export locality maps (the c-rider locality-maps-export at K=63)
if [ ! -f "$P/maps.done" ]; then
  rm -rf "$P/locality-graphs"
  ee_step "$P" maps python3 "$CRIDER/locality-maps-export.py" \
    "$AXIS_PATH" "$BEDS" "$P/export-graphs" "$P/export-walks.jsonl" \
    "$P/locality-graphs" "$P/reference-axis.json" "$P" \
    "$P/locality-domain-structure.txt" 63 || exit 1
fi
touch "$P/substrate.done"
echo "$(ee_now) $T substrate done"
