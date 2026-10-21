#!/usr/bin/env bash
# THE COMPLEXITY LEDGER'S GRAIN RE-MEASURE (the owner's bigger-loci lever):
# the window-grain cost curve re-measured under the CURRENT binary (the
# target-frontier fix build: bounded class scan + anchor hoist + the
# committed defaults), at the five-locality chain's own committed recipe
# (run-grain-chain.sh verbatim in stage shape), for ONE locality (L1)
# at 10kb/20kb/40kb/50kb. EVERY grain's substrate is REBUILT whole into
# this work dir (the committed chain's and the 2026-10-09 sweep's own
# beds/routes/census were deleted by their disk disciplines - only the
# receipts survive), so the four points share ONE binary and ONE
# substrate era (the T00041-era lesson: comparisons need one era; the
# committed 10kb receipts stay the cited accuracy baseline beside it).
# The CLI calls run SERIAL (RAYON_NUM_THREADS=1: wall = core-seconds,
# the autopsy's own method) with --truth-qv (the accuracy sidecar). The
# far-bin prune requires --global-frame and is therefore INERT on this
# local-frame chain - stated in the ledger doc.
# usage: run-grain-L1.sh <10000|20000|40000|50000>
set -uo pipefail
W=$1
G=/home/erikg/hprcv2/generic-sampling
IMPG=/home/erikg/impg-genome-inference/target-frontier/release/impg
# substrate tools (partition/export/walks/census/cut/projection) UNCHANGED
# since the Oct-19 shared build (git diff 82583eb..HEAD touches only
# commands/genome_infer.rs, realign_score.rs, examples/class_bound_probe.rs);
# the NFS wedge (pkg-config on -sys crates against the dead /gnu/store)
# makes a target-frontier examples rebuild a one-shot gamble - the
# committed target/release/examples are used instead. ALL call steps run
# the CURRENT binary (target-frontier/release/impg, the fix build at HEAD).
EXAMPLES=/home/erikg/impg-genome-inference/target/release/examples
QVBIWFA=/home/erikg/impg-genome-inference/genome/instrumented/qv-biwfa/target/release/qv-biwfa
BASE=/home/erikg/hprcv2/global-panel-call-mode/complexity-ledger
P0=$G/L1
PANEL_ROOT=$G/L1
READS="$P0/hg002.L1.300x.all.fq.gz"
OUT=$BASE/L1-w$((W/1000))k
P=$OUT
mkdir -p "$OUT"
AXIS_PATH=$(awk -F'\t' '$1 ~ /^HG002#1#/ {print $1":"$2"-"$3; exit}' "$P0/locus.refined.bed")
TRUTH2_PATH=$(awk -F'\t' '$1 ~ /^HG002#2#/ {print $1":"$2"-"$3; exit}' "$P0/locus.refined.bed")
PANEL="$P/locus-pos64.syng"
echo "=== L1 w$W: axis $AXIS_PATH truth2 $TRUTH2_PATH"

guard_poll() { local pid=$1 rssf=$2 cap=${3:-$((64*1024*1024))}
  ( while kill -0 "$pid" 2>/dev/null; do
      rss_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
      echo "$(date -u +%FT%TZ) $pid ${rss_kb:-0}" >> "$rssf"
      if [ "${rss_kb:-0}" -gt "$cap" ]; then echo "$(date -u +%FT%TZ) GUARD kill ${rss_kb}kB" >> "$rssf"; kill -9 "$pid" 2>/dev/null; break; fi
      sleep 10
    done ) &
}
runjob() { local name=$1 cap=$2; shift 2
  if [ -f "$OUT/$name.done" ]; then echo "=== $name done (marker)"; return 0; fi
  local start pid st; start=$(date +%s)
  nice -n 10 stdbuf -oL -eL "$@" > "$OUT/$name.log" 2> "$OUT/$name.err" &
  pid=$!; guard_poll "$pid" "$OUT/$name.rss" "$cap"; local poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$OUT/$name.exit"; echo "$(( $(date +%s) - start ))" > "$OUT/$name.wall"
  [ "$st" -eq 0 ] && touch "$OUT/$name.done" || { echo "=== $name FAILED ($st)"; return 1; }
  return 0
}

BEDS="$P/partition-w$((W/1000))k-d$((W/10000))k/results"; ROUTES="$P/routes"; LGRAPHS="$P/locality-graphs"
RAX="$P/reference-axis.json"
CENSUS="$P/cosine-multi-census.jsonl"; DCACHE="$P/realign-derive-cache.bin"

# stage 0: the locality syng (the committed recipe, rebuilt at this grain wave)
runjob s0-syng $((32*1024*1024)) "$IMPG" syng -f "$P0/locus.fa" --syncmer-length 63 --smer-length 8 \
  --syncmer-seed 7 --position-sample-rate 64 -o "$PANEL" || exit 1

# stage 1: the partition at this grain (the committed options verbatim)
mkdir -p "$BEDS"
printf '%s\n' "$AXIS_PATH" > "$P/partition-starting-sequences.txt"
runjob s1-partition $((32*1024*1024)) "$IMPG" partition \
  -a "$PANEL" -w "$W" -d $((W/10)) \
  --starting-sequences-file "$P/partition-starting-sequences.txt" \
  --selection-mode haplotype \
  --syng-min-chain-anchors 5 --syng-min-chain-fraction 0.5 \
  --no-rehome-singletons \
  --sequence-files "$P0/locus.fa" \
  -o bed --separate-files --output-folder "$BEDS" \
  -t 16 || exit 1
NP=$(ls "$BEDS" | grep -c '\.bed$')
echo "=== partition done: $NP partitions"

# stage 2: the catalog + the routes from this grain's beds
runjob s2-catalog $((32*1024*1024)) "$IMPG" genome-infer build-catalog \
  --bed-dir "$BEDS" --panel "$PANEL" --out-dir "$P/catalog" || exit 1
runjob s2b-routes $((32*1024*1024)) "$IMPG" genome-infer build-panel-routes \
  --panel "$PANEL" --catalog "$P/catalog/catalog.json" \
  --sources "$P0/locus.fa" --read-lengths 148 --out-dir "$ROUTES" || exit 1
rm -f "$P/catalog/catalog.json"

# stages 3-5: export graphs, the axis-only pggb lanes, the projection
IDS=$(seq 0 $((NP-1)) | tr '\n' ' ')
mkdir -p "$P/export-graphs"
runjob s3-export $((32*1024*1024)) "$EXAMPLES/partition_graph_export" export "$PANEL" "$P0/locus.fa" "$BEDS" "$P/export-graphs" $IDS || exit 1
AXIS_PIDS=$(for i in $(seq 0 $((NP-1))); do grep -qF "$AXIS_PATH" "$BEDS/partition$i.bed" 2>/dev/null && echo $i; done | tr '\n' ' ')
NAX=$(echo $AXIS_PIDS | wc -w); echo "=== axis partitions: $NAX of $NP ($AXIS_PIDS)"
if [ ! -f "$P/s4-pggb.done" ]; then
  start=$(date +%s)
  LOCUS_DIR="$P0" BED_ROOT="$BEDS" PGGB_OUT="$P/pggb-graphs" \
    bash "$G/run-pggb-builds-generic.sh" $AXIS_PIDS > "$P/s4-pggb.log" 2>&1
  bad=0; for i in $AXIS_PIDS; do
    [ "$(cat "$P/pggb-graphs/partition$i/graph.gfa.exit" 2>/dev/null)" = "0" ] || bad=1; done
  [ "$bad" = 0 ] && { touch "$P/s4-pggb.done"; echo "$(( $(date +%s) - start ))" > "$P/s4-pggb.wall"; } \
    || { echo "=== pggb FAILED (the sliver-class crash?) - see $P/s4-pggb.log"; exit 1; }
fi
for i in $AXIS_PIDS; do cp "$BEDS/partition$i.bed" "$P/pggb-graphs/partition$i/members.bed" 2>/dev/null; done
mkdir -p "$P/pggb-projection"
runjob s5-projection $((32*1024*1024)) "$EXAMPLES/pggb_projection" --panel "$PANEL" --agc "$P0/locus.fa" \
  --export-graphs "$P/export-graphs" --pggb-dir "$P/pggb-graphs" \
  --out-dir "$P/pggb-projection" --partitions "$(echo $AXIS_PIDS | tr ' ' ',')" || exit 1

# stage 6: locality maps (the committed pggb-era domain machinery)
if [ ! -f "$P/s6-maps.done" ]; then
  start=$(date +%s)
  python3 "$G/locality-maps.py" "$AXIS_PATH" "$BEDS" "$P/export-graphs" \
    "$P/pggb-graphs" "$P/pggb-projection" "$LGRAPHS" \
    "$RAX" "$P" "$P/locality-domain-structure.txt" \
    > "$P/s6-maps.log" 2> "$P/s6-maps.err" \
    && { touch "$P/s6-maps.done"; echo "$(( $(date +%s) - start ))" > "$P/s6-maps.wall"; } || { echo "=== maps FAILED"; exit 1; }
fi
NW=$(ls "$LGRAPHS" | grep -c 'window[0-9]*.gfa.map.json')
echo "=== locality maps done: $NW windows"

# stage 7: the deep census (the committed env; the committed depth receipt 137.3x)
if [ ! -f "$P/s7-census.done" ]; then
  start=$(date +%s)
  IMPG_FORCE_DIPLOID_DIAGNOSTIC=1 IMPG_MULTIPLICITY_VARIANT=gentle-beta \
  IMPG_COSINE_READ_MATCHED=1 \
  IMPG_COSINE_MULTI_CENSUS="$CENSUS" \
  IMPG_COSINE_DOMAIN_PILOT=1 \
    taskset -c 0-63 nice -n 10 \
    stdbuf -oL -eL "$EXAMPLES/panel_route_mem_routed" \
    --panel "$PANEL" --routes "$ROUTES" \
    --sample "$P0/sample/sample.membwt" \
    --reads "$READS" \
    --axis "$RAX" \
    --bed-directory "$BEDS" \
    --component "$AXIS_PATH" \
    --depth 137.3 --read-length 148 \
    --routing-universe component --spine --spine-phasing \
    --reference-pair "truth:$P0/truth-routes/HG002-1.json:$P0/truth-routes/HG002-2.json" \
    > "$P/s7-census.log" 2> "$P/s7-census.err" &
  cpid=$!; guard_poll "$cpid" "$P/s7-census.rss" $((64*1024*1024)); cpoll=$!
  cst=0; wait "$cpid" || cst=$?; kill "$cpoll" 2>/dev/null; wait "$cpoll" 2>/dev/null
  echo "$(( $(date +%s) - start ))" > "$P/s7-census.wall"
  [ "$cst" -eq 0 ] && touch "$P/s7-census.done" || { echo "=== census FAILED ($cst)"; exit 1; }
fi

# stage 8: the deep derive cache
if [ ! -f "$P/s8-derive-cache.done" ]; then
  python3 -c "
import sys
for i in range($NW): print('{\"locus\": %d, \"records\": []}' % i)" > "$P/ingredients-stub.jsonl"
  start=$(date +%s)
  nice -n 10 stdbuf -oL -eL "$EXAMPLES/partition_realign_likelihood" \
    --panel "$PANEL" --routes "$ROUTES" \
    --partition-graphs "$LGRAPHS" \
    --component "$AXIS_PATH" \
    --census "$CENSUS" \
    --ingredients "$P/ingredients-stub.jsonl" \
    --reads "$READS" \
    --derive-cache "$DCACHE" \
    --read-length 148 \
    --out "$P/realign-stub-out.jsonl" \
    > "$P/s8-derive-cache.log" 2> "$P/s8-derive-cache.err" &
  pid=$!; guard_poll "$pid" "$P/s8-derive-cache.rss" $((64*1024*1024)); poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$(( $(date +%s) - start ))" > "$P/s8-derive-cache.wall"
  if [ "$st" -eq 0 ] || { [ -s "$DCACHE" ] && grep -q 'the read hull differs from the placed path segment' "$P/s8-derive-cache.err"; }; then
    touch "$P/s8-derive-cache.done"
  else echo "=== derive cache FAILED ($st)"; exit 1; fi
fi

# stage 9: THE CLI PRODUCT RUN at the CURRENT binary, SERIAL, with --truth-qv
LOCI=$(python3 -c "
import json
a=json.load(open('$RAX'))
n=sum(1 for iv in a['intervals'] if iv['component'] == '$AXIS_PATH')
print(','.join(str(i) for i in range(n)))")
echo "=== CLI loci: $LOCI"
if [ ! -f "$P/s9-cli-deep.done" ]; then
  rm -rf "$P/cli-deep"   # the call refuses an existing out-dir (retry-clean)
  start=$(date +%s)
  RAYON_NUM_THREADS=1 IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS=1 IMPG_PANEL_LEGACY_V1_OK=1 \
    nice -n 10 stdbuf -oL -eL "$IMPG" genome-infer call-diplotypes \
    --panel "$PANEL" \
    --routes "$ROUTES" \
    --partition-graphs "$LGRAPHS" \
    --component "$AXIS_PATH" \
    --census "$CENSUS" \
    --reads "$READS" \
    --derive-cache "$DCACHE" \
    --loci "$LOCI" \
    --truth-qv "$AXIS_PATH" "$TRUTH2_PATH" \
    --far-bin-prune=false \
    --qv-biwfa "$QVBIWFA" \
    --out-dir "$P/cli-deep" \
    > "$P/s9-cli-deep.log" 2> "$P/s9-cli-deep.err" &
  pid=$!; guard_poll "$pid" "$P/s9-cli-deep.rss" $((64*1024*1024)); poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$P/s9-cli-deep.exit"; echo "$(( $(date +%s) - start ))" > "$P/s9-cli-deep.wall"
  [ "$st" -eq 0 ] && touch "$P/s9-cli-deep.done" || { echo "=== CLI FAILED ($st)"; exit 1; }
fi
echo "=== L1 w$W COMPLETE: $(grep -m1 'phase summary' "$P/s9-cli-deep.err")"
