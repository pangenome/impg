#!/usr/bin/env bash
# THE COMPLEXITY LEDGER'S GRAIN RE-MEASURE, the deep-stratum territory leg:
# T03226 (the S_EXP monster: 43.4% far-bin share, the receipted 1,486.8s
# serial loci phase at 10kb) re-grained to 20kb/40kb under the CURRENT
# architecture (the global-frame production stack: the target-frontier
# fix build, the far-bin prune DEFAULT, the bounded class scan, the
# anchor hoist) by replicating the tranche's own prep chain
# (run-prep-one.sh VERBATIM in stage shape, -w changed) + the single
# restricted cut (run-cutA.sh's own command) + the tranche callC recipe
# (assessment form, --truth-qv, serial RAYON=1 = core-seconds) + the
# production form (no --truth-qv, T16, the bounded scan engaged).
# The restriction basis (rvcut's cut-reads.fasta) and the reads are
# GRAIN-INDEPENDENT and REUSED read-only from the tranche; the tranche's
# own dirs are never written. The 10kb point re-runs the SAME recipe over
# the tranche's committed 10kb substrate (read-only) at the current
# binary for wave comparability with the new grains.
# usage: run-grain-T03226.sh <10000|20000|40000>
set -uo pipefail
W=$1
T=T03226
WT1=/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1   # READ-ONLY (the masters own it)
RV=/home/erikg/hprcv2/wgs35x-fleet/runtime-verification        # READ-ONLY (the campaign's ledgers)
HPRC=/home/erikg/hprcv2
BASE=/home/erikg/hprcv2/global-panel-call-mode/complexity-ledger
IMPG=/home/erikg/impg-genome-inference/target-frontier/release/impg
# substrate tools (partition/export/walks/census/cut/projection) UNCHANGED since the Oct-19 shared
# build (git diff 82583eb..HEAD touches only commands/genome_infer.rs, realign_score.rs,
# examples/class_bound_probe.rs); the NFS wedge (pkg-config on -sys crates against the dead
# /gnu/store) makes a target-frontier examples rebuild a one-shot gamble - the committed
# target/release/examples are used instead. ALL call steps run the CURRENT binary
# (target-frontier/release/impg, the fix build at HEAD).
EXAMPLES=/home/erikg/impg-genome-inference/target/release/examples
QVBIWFA=/home/erikg/impg-genome-inference/genome/instrumented/qv-biwfa/target/release/qv-biwfa
ROUTES=$HPRC/global-panel-call-mode/global-routes-par
AXIS_PATH=$(cat "$WT1/prep/$T/axis-path.txt")
TRUTH2_PATH=HG002#2#chr5:9034962-9284684
OUT=$BASE/T03226-w$((W/1000))k
P=$OUT/prep
mkdir -p "$P" "$OUT"

guard_poll() { local pid=$1 rssf=$2 cap=${3:-$((64*1024*1024))}
  ( while kill -0 "$pid" 2>/dev/null; do
      rss_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
      echo "$(date -u +%FT%TZ) $pid ${rss_kb:-0}" >> "$rssf"
      if [ "${rss_kb:-0}" -gt "$cap" ]; then echo "$(date -u +%FT%TZ) GUARD kill ${rss_kb}kB" >> "$rssf"; kill -9 "$pid" 2>/dev/null; break; fi
      sleep 10
    done ) &
}
step() { local name=$1 cap=$2; shift 2
  if [ -f "$P/$name.done" ]; then echo "=== $name done (marker)"; return 0; fi
  local start pid st; start=$(date +%s)
  nice -n 10 stdbuf -oL -eL "$@" > "$P/$name.log" 2> "$P/$name.err" &
  pid=$!; guard_poll "$pid" "$P/$name.rss" "$cap"; local poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$P/$name.exit"; echo "$(( $(date +%s) - start ))" > "$P/$name.wall"
  [ "$st" -eq 0 ] && touch "$P/$name.done" || { echo "=== $name FAILED ($st)"; return 1; }
  return 0
}

if [ "$W" = 10000 ]; then
  # the committed 10kb substrate, reused READ-ONLY
  BEDS="$WT1/prep/$T/partition-w10k-d1k/results"
  LGRAPHS="$WT1/prep/$T/locality-graphs"
  RAX="$WT1/prep/$T/reference-axis.json"
  CENSUS="$WT1/cutA/$T/cut-census.jsonl"; DCACHE="$WT1/cutA/$T/cut-cache.bin"
else
  BEDS="$P/partition-w$((W/1000))k-d$((W/10000))k/results"; LGRAPHS="$P/locality-graphs"
  RAX="$P/reference-axis.json"
  CENSUS="$OUT/cutA/cut-census.jsonl"; DCACHE="$OUT/cutA/cut-cache.bin"
  # the grain-independent inputs, symlinked read-only from the tranche
  for f in locus.fa locus.fa.fai locus-pos64.syng locus-pos64.syng.1gbwt \
           locus-pos64.syng.1khash locus-pos64.syng.meta locus-pos64.syng.names \
           locus-pos64.syng.pstep locus-pos64.syng.spos; do
    [ -e "$P/$f" ] || ln -s "$WT1/prep/$T/$f" "$P/$f"
  done
  for f in locus.bed locus.members.bed locus.refined.bed split.info.json; do
    [ -e "$P/$f" ] || cp "$WT1/prep/$T/$f" "$P/$f"
  done
  # ---- the tranche prep chain VERBATIM, -w changed ----
  mkdir -p "$BEDS"
  printf '%s\n' "$AXIS_PATH" > "$P/starting-sequences.txt"
  step beds $((32*1024*1024)) "$IMPG" partition \
    -a "$P/locus-pos64.syng" -w "$W" -d $((W/10)) \
    --starting-sequences-file "$P/starting-sequences.txt" \
    --selection-mode haplotype \
    --syng-min-chain-anchors 5 --syng-min-chain-fraction 0.5 \
    --no-rehome-singletons \
    --sequence-files "$P/locus.fa" \
    -o bed --separate-files --output-folder "$BEDS" \
    -t 8 || exit 1
  NP=$(ls "$BEDS" | grep -c '\.bed$'); echo "=== $T w$W: $NP partitions"
  IDS=$(seq 0 $((NP-1)) | tr '\n' ' ')
  mkdir -p "$P/export-graphs"
  step export $((32*1024*1024)) "$EXAMPLES/partition_graph_export" export "$P/locus-pos64.syng" \
    "$P/locus.fa" "$BEDS" "$P/export-graphs" $IDS || exit 1
  if [ ! -f "$P/walks.done" ]; then
    python3 "$RV/rv-mk-walk-requests.py" "$BEDS" "$P/locus-pos64.syng.names" \
      "$P/walk-requests.tsv" || exit 1
    start=$(date +%s)
    nice -n 10 "$EXAMPLES/syng_homology_dump" walks "$P/locus-pos64.syng" "$P/walk-requests.tsv" \
      > "$P/export-walks.jsonl" 2> "$P/walks.err" \
      && { touch "$P/walks.done"; echo "$(( $(date +%s) - start ))" > "$P/walks.wall"; } \
      || { echo 1 > "$P/walks.exit"; exit 1; }
  fi
  if [ ! -f "$P/maps.done" ]; then
    rm -rf "$LGRAPHS"
    start=$(date +%s)
    python3 /home/erikg/hprcv2/wgs35x-fleet/c-rider/locality-maps-export.py \
      "$AXIS_PATH" "$BEDS" "$P/export-graphs" "$P/export-walks.jsonl" \
      "$LGRAPHS" "$RAX" "$P" "$P/locality-domain-structure.txt" 63 \
      > "$P/maps.log" 2> "$P/maps.err" \
      && { touch "$P/maps.done"; echo "$(( $(date +%s) - start ))" > "$P/maps.wall"; } \
      || { echo "=== maps FAILED"; exit 1; }
  fi
  # ---- the single restricted cut (run-cutA.sh's own command, RAYON=32 co-tenant courtesy) ----
  if [ ! -f "$OUT/cutA.done" ]; then
    start=$(date +%s)
    IMPG_CUT_READSET_FASTA="$WT1/rvcut/$T/$T/cut-reads.fasta" \
    RAYON_NUM_THREADS=32 nice -n 10 stdbuf -oL -eL "$EXAMPLES/sample_projection" cut --frame global \
      --record-origin 1 \
      --global-panel $HPRC/HPRC_r2_assemblies_0.6.1.syng \
      --axis "$RAX" --bed-directory "$BEDS" --component "$AXIS_PATH" \
      --agc $HPRC/HPRC_r2_assemblies_0.6.1.agc \
      --shards-dir $HPRC/hg002-projection/full \
      --out-dir "$OUT/cutA" \
      > "$OUT/cutA.log" 2> "$OUT/cutA.err" &
    pid=$!; guard_poll "$pid" "$OUT/cutA.rss" $((120*1024*1024)); poll=$!
    st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
    echo "$st" > "$OUT/cutA.exit"; echo "$(( $(date +%s) - start ))" > "$OUT/cutA.wall"
    { [ "$st" -eq 0 ] && [ -f "$OUT/cutA/cut-stats.json" ]; } && touch "$OUT/cutA.done" \
      || { echo "=== cutA FAILED ($st)"; exit 1; }
  fi
fi

NW=$(ls "$LGRAPHS" | grep -c 'window[0-9]*.gfa.map.json')
LOCI=$(seq -s, 0 $((NW-1)))
echo "=== $T w$W: $NW windows, loci $LOCI"

# ---- the C-form call (assessment, --truth-qv, far-bin prune, SERIAL = core-s) ----
if [ ! -f "$OUT/callC.done" ]; then
  start=$(date +%s)
  RAYON_NUM_THREADS=1 nice -n 10 stdbuf -oL -eL "$IMPG" genome-infer call-diplotypes \
    --panel $HPRC/HPRC_r2_assemblies_0.6.1.syng \
    --routes "$ROUTES" \
    --partition-graphs "$LGRAPHS" \
    --component "$AXIS_PATH" \
    --census "$CENSUS" \
    --reads "$WT1/reads/$T.fq.gz" \
    --derive-cache "$DCACHE" \
    --loci "$LOCI" \
    --truth-qv "$AXIS_PATH" "$TRUTH2_PATH" \
    --far-bin-prune \
    --qv-biwfa "$QVBIWFA" \
    --global-frame \
    --record-origin 1 \
    --out-dir "$OUT/callC" \
    > "$OUT/callC.log" 2> "$OUT/callC.err" &
  pid=$!; guard_poll "$pid" "$OUT/callC.rss" $((300*1024*1024)); poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$OUT/callC.exit"; echo "$(( $(date +%s) - start ))" > "$OUT/callC.wall"
  [ "$st" -eq 0 ] && touch "$OUT/callC.done" || echo "=== callC FAILED ($st) - see $OUT/callC.err"
fi

# ---- the production form (no --truth-qv, T16, the bounded scan engaged) ----
if [ ! -f "$OUT/callProd.done" ]; then
  start=$(date +%s)
  RAYON_NUM_THREADS=16 nice -n 10 stdbuf -oL -eL "$IMPG" genome-infer call-diplotypes \
    --panel $HPRC/HPRC_r2_assemblies_0.6.1.syng \
    --routes "$ROUTES" \
    --partition-graphs "$LGRAPHS" \
    --component "$AXIS_PATH" \
    --census "$CENSUS" \
    --reads "$WT1/reads/$T.fq.gz" \
    --derive-cache "$DCACHE" \
    --loci "$LOCI" \
    --far-bin-prune \
    --qv-biwfa "$QVBIWFA" \
    --global-frame \
    --record-origin 1 \
    --out-dir "$OUT/callProd" \
    > "$OUT/callProd.log" 2> "$OUT/callProd.err" &
  pid=$!; guard_poll "$pid" "$OUT/callProd.rss" $((300*1024*1024)); poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$OUT/callProd.exit"; echo "$(( $(date +%s) - start ))" > "$OUT/callProd.wall"
  [ "$st" -eq 0 ] && touch "$OUT/callProd.done" || echo "=== callProd FAILED ($st)"
fi
echo "=== $T w$W COMPLETE"
grep -m1 'phase summary' "$OUT/callC.err" || true
