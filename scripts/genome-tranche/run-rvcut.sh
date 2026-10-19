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

one_batch() {
  local B=$1; shift
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
    else
      echo "$T no cut-stats after batch (exit $st)" >> "$CUTS/campaign.log"
    fi
  done
  echo "$(date -u +%FT%TZ) batch $B done exit=$st wall=$(cat "$B/batch.wall")s" >> "$CUTS/campaign.log"
}

BID=0
while true; do
  mapfile -t PENDING < <(grep -v '^#' "$WT/tranche.tsv" | cut -f1 | while read -r T; do
      [ -f "$WT/prep/$T/substrate.done" ] || continue
      [ -f "$WT/prep/$T/rvcut.done" ] && continue
      echo "$T"
    done)
  [ "${#PENDING[@]}" -eq 0 ] && break
  # pack W pending per batch, run C batches concurrently
  JOBS=()
  while [ "${#PENDING[@]}" -gt 0 ] && [ "${#JOBS[@]}" -lt "$C" ]; do
    B=$CUTS/batch-$BID; BID=$((BID+1))
    mapfile -t SLICE < <(printf '%s\n' "${PENDING[@]}" | head -"$W")
    JOBS+=("$B ${SLICE[*]}")
    PENDING=("${PENDING[@]:$W}")
  done
  for j in "${JOBS[@]}"; do
    one_batch $j &
  done
  wait
done
echo "$(date -u +%FT%TZ) rvcut stage complete: $(ls $WT/prep/*/rvcut.done 2>/dev/null | wc -l) territories"
