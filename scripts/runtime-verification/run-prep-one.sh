#!/usr/bin/env bash
# THE WHOLE-GENOME RUNTIME VERIFICATION - STAGE 2: the per-territory
# EXPORT-PROFILE substrate chain (the c-rider recipe at K=63, the syng-native
# export maps - THE RUNTIME PROFILE, honestly degraded vs pggb per the
# k-sweep receipt): extract -> syng -> partition beds (the committed
# -w10k -d1k haplotype recipe) -> export graphs -> panel walks ->
# locality-maps-export. Marker-gated per step; resumable; per-step walls.
# usage: run-prep-one.sh <TID>   (env: RV, IMPG, EXAMPLES, THREADS)
set -uo pipefail
RV=${RV:-/home/erikg/hprcv2/wgs35x-fleet/runtime-verification}
IMPG=${IMPG:-/home/erikg/impg-genome-inference/target/release/impg}
EXAMPLES=${EXAMPLES:-/home/erikg/impg-genome-inference/target/release/examples}
THREADS=${THREADS:-4}
T=$1
P=$RV/prep/$T
[ -f "$P/substrate.done" ] && exit 0
[ -f "$P/axis-path.txt" ] || { echo "$T no-axis (split not accepted)" > "$P/substrate.skipped"; exit 0; }
AXIS_PATH=$(cat "$P/axis-path.txt")

step() { # step <name> <cmd...>
  local name=$1; shift
  [ -f "$P/$name.done" ] && return 0
  local start
  start=$(date +%s)
  "$@" > "$P/$name.log" 2> "$P/$name.err" \
    && { touch "$P/$name.done"; echo "$(( $(date +%s) - start ))" > "$P/$name.wall"; return 0; } \
    || { echo 1 > "$P/$name.exit"; echo "$(( $(date +%s) - start ))" > "$P/$name.wall"; return 1; }
}

# 1. extract
if [ ! -f "$P/extract.done" ]; then
  start=$(date +%s)
  python3 "$RV/rv-extract.py" "$RV" "$T" > "$P/extract.log" 2> "$P/extract.err" \
    && { touch "$P/extract.done"; echo "$(( $(date +%s) - start ))" > "$P/extract.wall"; } \
    || { echo 1 > "$P/extract.exit"; echo "$(( $(date +%s) - start ))" > "$P/extract.wall"; exit 1; }
fi

# 2. the local syng (the committed recipe: k63 s8 seed7 pos64)
step syng "$IMPG" syng -f "$P/locus.fa" --syncmer-length 63 --smer-length 8 \
  --syncmer-seed 7 --position-sample-rate 64 -o "$P/locus-pos64.syng" || exit 1

# 3. the partition beds (the committed recipe: w10k d1k haplotype)
BEDS="$P/partition-w10k-d1k/results"
if [ ! -f "$P/beds.done" ]; then
  mkdir -p "$P/partition-w10k-d1k" "$BEDS"
  printf '%s\n' "$AXIS_PATH" > "$P/partition-w10k-d1k/starting-sequences.txt"
  step beds nice -n 10 "$IMPG" partition \
    -a "$P/locus-pos64.syng" -w 10000 -d 1000 \
    --starting-sequences-file "$P/partition-w10k-d1k/starting-sequences.txt" \
    --selection-mode haplotype \
    --syng-min-chain-anchors 5 --syng-min-chain-fraction 0.5 \
    --no-rehome-singletons \
    --sequence-files "$P/locus.fa" \
    -o bed --separate-files --output-folder "$BEDS" \
    -t "$THREADS" || exit 1
fi
NP=$(ls "$BEDS" | grep -c '\.bed$')
[ "$NP" -gt 0 ] || { echo "empty partition at $T" > "$P/beds.empty"; exit 1; }

# 4. the export graphs (the syng-native range GFAs)
if [ ! -f "$P/export.done" ]; then
  mkdir -p "$P/export-graphs"
  IDS=$(seq 0 $((NP-1)) | tr '\n' ' ')
  step export "$EXAMPLES/partition_graph_export" export "$P/locus-pos64.syng" \
    "$P/locus.fa" "$BEDS" "$P/export-graphs" $IDS || exit 1
fi

# 5. the panel walks (ONE syng load; the panel's own stored walks)
if [ ! -f "$P/walks.done" ]; then
  step walkreq python3 "$RV/rv-mk-walk-requests.py" "$BEDS" "$P/locus-pos64.syng.names" \
    "$P/walk-requests.tsv" || exit 1
  start=$(date +%s)
  nice -n 10 "$EXAMPLES/syng_homology_dump" walks "$P/locus-pos64.syng" "$P/walk-requests.tsv" \
    > "$P/export-walks.jsonl" 2> "$P/walks.err" \
    && { touch "$P/walks.done"; echo "$(( $(date +%s) - start ))" > "$P/walks.wall"; } \
    || { echo 1 > "$P/walks.exit"; echo "$(( $(date +%s) - start ))" > "$P/walks.wall"; exit 1; }
fi

# 6. the export locality maps (the c-rider locality-maps-export at K=63)
if [ ! -f "$P/maps.done" ]; then
  rm -rf "$P/locality-graphs"
  step maps python3 /home/erikg/hprcv2/wgs35x-fleet/c-rider/locality-maps-export.py \
    "$AXIS_PATH" "$BEDS" "$P/export-graphs" "$P/export-walks.jsonl" \
    "$P/locality-graphs" "$P/reference-axis.json" "$P" \
    "$P/locality-domain-structure.txt" 63 || exit 1
fi
touch "$P/substrate.done"
echo "$(date -u +%FT%TZ) $T substrate done"
