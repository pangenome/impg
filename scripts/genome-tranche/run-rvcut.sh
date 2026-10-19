#!/usr/bin/env bash
# THE GENOME TRANCHE - stage 3: the RV local-derive cut (the committed
# kept-read set = the restriction basis). The campaign's own cut form
# (run-cuts.sh VERBATIM: sample_projection cut-batch, IMPG_CUT_LOCAL_DERIVE=1,
# one global syng+AGC load + one shard scan per W=8 batch), re-run over
# the tranche's regenerated substrates to produce each territory's
# cut-reads.fasta (the committed kept-read set the restrict census
# profile states via IMPG_CUT_READSET_FASTA) - and a free determinism
# gate: where the campaign banked a cut-stats receipt, the regenerated
# reads_kept/census_mult/routed_patterns must reproduce it exactly.
# usage: run-rvcut.sh   (env: WT, BATCH_W=8, N_CONCURRENT=2, RAYON=24)
set -uo pipefail
WT=${WT:-/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1}
RV=/home/erikg/hprcv2/wgs35x-fleet/runtime-verification
BIN=/home/erikg/impg-genome-inference/target/release/examples/sample_projection
SYNG=/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng
AGC=/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc
SHARDS=/home/erikg/hprcv2/hg002-projection/full
W=${BATCH_W:-8}
C=${N_CONCURRENT:-2}
RAYON=${RAYON:-24}
CUTS=$WT/rvcut
mkdir -p "$CUTS" "$WT/logs"

one_batch_w() {
  local B=$1; local w=$2; shift 2
  local TODO=("$@")
  mkdir -p "$B"
  : > "$B/manifest.jsonl"
  local T P
  for T in "${TODO[@]}"; do
    [ -n "$T" ] || continue
    P=$WT/prep/$T
    rm -rf "$CUTS/$T/$T"; mkdir -p "$CUTS/$T/$T"
    printf '{"name": "%s", "panel": "%s/prep/%s/locus-pos64.syng", "axis": "%s/prep/%s/reference-axis.json", "bed_directory": "%s/prep/%s/partition-w10k-d1k/results", "component": "%s", "sources": "%s/prep/%s/locus.fa", "out_dir": "%s/rvcut/%s/%s"}\n' \
      "$T" "$WT" "$T" "$WT" "$T" "$WT" "$T" "$(cat "$P/axis-path.txt")" "$WT" "$T" "$WT" "$T" "$T" >> "$B/manifest.jsonl"
  done
  echo "$(date -u +%FT%TZ) batch $B: ${#TODO[@]} territories" >> "$CUTS/campaign.log"
  local start
  start=$(date +%s)
  RAYON_NUM_THREADS=$RAYON IMPG_CUT_LOCAL_DERIVE=1 \
    nice -n 10 stdbuf -oL -eL "$BIN" cut-batch \
    --global-panel "$SYNG" --agc "$AGC" --shards-dir "$SHARDS" \
    --manifest "$B/manifest.jsonl" > "$B/batch.log" 2> "$B/batch.err" &
  local pid=$!
  ( while kill -0 "$pid" 2>/dev/null; do
      rss_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
      echo "$(date -u +%FT%TZ) $rss_kb" >> "$B/batch.rss"
      if [ "${rss_kb:-0}" -gt $((64*1024*1024)) ]; then
        echo "$(date -u +%FT%TZ) GUARD kill ${rss_kb}kB" >> "$B/batch.rss"
        kill -9 "$pid" 2>/dev/null; break
      fi
      sleep 10
    done ) &
  local poll=$!
  local st=0; wait "$pid" || st=$?
  kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$B/batch.exit"
  echo "$(( $(date +%s) - start ))" > "$B/batch.wall"
  for T in "${TODO[@]}"; do
    local D=$CUTS/$T/$T
    if [ -f "$D/cut-stats.json" ]; then
      touch "$WT/prep/$T/rvcut.done"
      # THE FREE DETERMINISM GATE vs the campaign's banked receipt
      if [ -f "$RV/receipts/$T/cut-stats.json" ]; then
        python3 - "$D/cut-stats.json" "$RV/receipts/$T/cut-stats.json" >> "$CUTS/determinism-gate.txt" <<'EOF' || echo "$T GATE-DIFFERS" >> "$CUTS/determinism-gate.txt"
import json, sys
a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2]))
keys = ("reads_kept", "census_mult", "routed_patterns")
assert all(a.get(k) == b.get(k) for k in keys), (sys.argv[1], {k: (a.get(k), b.get(k)) for k in keys})
print("IDENTICAL " + sys.argv[1].split("/")[-2], {k: a[k] for k in keys})
EOF
      fi
      # THE JIT REAP: only the kept-read set (the restriction basis)
      # and the stats ride forward; the census/cache stay on disk until
      # the batch verifies, then go (the disk floor is tight).
      rm -f "$D/cut-census.jsonl" "$D/cut-cache.bin"
    elif [ "$st" -eq 137 ]; then
      # THE GUARD KILL IS A BATCH-PACKAGING EVENT, never a territory
      # strike - the width HALVES for this batch's territories (the
      # guard is the only oracle: the census working set is not
      # predictable pre-batch)
      echo $(( w > 1 ? w / 2 : 1 )) > "$WT/prep/$T/cut.width"
      echo "$(date -u +%FT%TZ) $T batch GUARD-KILLED at width $w (packaging - no strike, width -> $(( w > 1 ? w / 2 : 1 )))" >> "$CUTS/campaign.log"
    else
      echo "$T no cut-stats after batch (exit $st)" >> "$CUTS/campaign.log"
    fi
  done
  echo "$(date -u +%FT%TZ) batch $B done exit=$st wall=$(cat "$B/batch.wall")s" >> "$CUTS/campaign.log"
}

BID=0
W_START=${W_START:-8}
# THE ADAPTIVE WIDTH (the campaign's own rule, run-cuts.sh VERBATIM in
# kind): the guard is the ONLY oracle; the territories of a GUARD-KILLED
# batch (exit 137) retry at HALF the width (prep/<T>/cut.width, min 1) -
# a guard kill never strikes a territory (the width's own failure).
while true; do
  mapfile -t PENDING < <(grep -v '^#' "$WT/tranche.tsv" | cut -f1 | while read -r T; do
      [ -f "$WT/prep/$T/rvcut.done" ] && continue
      [ -f "$WT/prep/$T/substrate.done" ] || continue
      grep -q "^$T\b" "$WT/tranche-failed.tsv" 2>/dev/null && continue
      echo "$T"
    done)
  [ "${#PENDING[@]}" -eq 0 ] && break
  # the campaign's own composition rule: the pending list is
  # DETERMINISTICALLY HASH-INTERLEAVED before batch composition, so
  # genome-adjacent density clumps spread across batches (no
  # cherry-picking - a pure function of the TID)
  mapfile -t PENDING < <(printf '%s\n' "${PENDING[@]}" | sort -k1,1)
  mapfile -t PENDING < <(for T in "${PENDING[@]}"; do echo "$(printf '%s' "$T" | md5sum | cut -c1-12) $T"; done | sort | awk '{print $2}')
  # pack per width class (the widest first), W_pending per batch, C
  # concurrent batches
  JOBS=()
  declare -A BYW
  for T in "${PENDING[@]}"; do
    w=$(cat "$WT/prep/$T/cut.width" 2>/dev/null || echo $W_START)
    BYW[$w]="${BYW[$w]:-} $T"
  done
  for w in $(printf '%s\n' "${!BYW[@]}" | sort -rn); do
    mapfile -t WL < <(printf '%s\n' ${BYW[$w]})
    while [ "${#WL[@]}" -gt 0 ] && [ "${#JOBS[@]}" -lt "$(( C * 3 ))" ]; do
      B=$CUTS/batch-$BID; BID=$((BID+1))
      mapfile -t SLICE < <(printf '%s\n' "${WL[@]}" | head -"$w")
      JOBS+=("$B $w ${SLICE[*]}")
      WL=("${WL[@]:$w}")
    done
    [ "${#JOBS[@]}" -ge "$C" ] && break
  done
  while [ "${#JOBS[@]}" -gt 0 ]; do
    for j in "${JOBS[@]:0:$C}"; do
      one_batch_w $j &
    done
    wait
    JOBS=("${JOBS[@]:$C}")
  done
  # no-progress guard: if nothing new completed this round, stop
  nd1=$(ls $CUTS/*/*/cut-stats.json 2>/dev/null | wc -l)
  [ "${nd1:-0}" -eq "${ND0:-0}" ] && echo "$(date -u +%FT%TZ) rvcut NO PROGRESS ($nd1 done)" >> "$CUTS/campaign.log" && break
  ND0=$nd1
done
echo "$(date -u +%FT%TZ) rvcut stage complete: $(ls $WT/prep/*/rvcut.done 2>/dev/null | wc -l) territories"
