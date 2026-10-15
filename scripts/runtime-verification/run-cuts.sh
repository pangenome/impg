#!/usr/bin/env bash
# THE WHOLE-GENOME RUNTIME VERIFICATION - STAGE 3: THE CUT CAMPAIGN (the
# sample-side map step): sample_projection cut-batch from the COMMITTED 35x
# projection shards - batches of W territories, ONE global syng + AGC load +
# ONE shard scan per batch (the amortized-load shape; the per-territory
# fail-closed bind failures ISOLATED by the patched emissions loop - the
# poison ledger names them). C concurrent batches; the disk floor gates every
# launch (the streaming reaper must keep up). Marker-gated, resumable.
# usage: run-cuts.sh   (env: RV, BATCH_W=32, N_CONCURRENT=2, RAYON=64,
#                       DISK_FLOOR_GB=1500)
set -uo pipefail
RV=${RV:-/home/erikg/hprcv2/wgs35x-fleet/runtime-verification}
BIN=/home/erikg/impg-genome-inference/target/release/examples/sample_projection
SYNG=/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng
AGC=/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc
SHARDS=/home/erikg/hprcv2/hg002-projection/full
W=${BATCH_W:-32}
C=${N_CONCURRENT:-2}
RAYON=${RAYON:-64}
FLOOR=${DISK_FLOOR_GB:-1500}
CUTS=$RV/cuts
mkdir -p "$CUTS"

avail_gb() { df -BG --output=avail / | tail -1 | tr -dc '0-9'; }

guard_poll() {
  local pid=$1 rssf=$2
  ( while kill -0 "$pid" 2>/dev/null; do
      rss_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
      echo "$(date -u +%FT%TZ) ${rss_kb:-0}" >> "$rssf"
      if [ "${rss_kb:-0}" -gt $((64*1024*1024)); then
        echo "$(date -u +%FT%TZ) GUARD kill ${rss_kb}kB" >> "$rssf"
        kill -9 "$pid" 2>/dev/null; break
      fi
      sleep 10
    done ) &
}

one_batch() {
  local B=$1; shift
  local TODO=("$@")
  mkdir -p "$B"
  : > "$B/manifest.jsonl"
  local T P
  for T in "${TODO[@]}"; do
    P=$RV/prep/$T
    rm -rf "$RV/map/$T/$T"; mkdir -p "$RV/map/$T/$T"
    printf '{"name": "%s", "panel": "%s/prep/%s/locus-pos64.syng", "axis": "%s/prep/%s/reference-axis.json", "bed_directory": "%s/prep/%s/partition-w10k-d1k/results", "component": "%s", "sources": "%s/prep/%s/locus.fa", "out_dir": "%s/map/%s/%s"}\n' \
      "$T" "$RV" "$T" "$RV" "$T" "$RV" "$T" "$(cat "$P/axis-path.txt")" "$RV" "$T" "$RV" "$T" "$T" >> "$B/manifest.jsonl"
  done
  echo "$(date -u +%FT%TZ) batch $B: ${#TODO[@]} territories (df $(avail_gb)G)" >> "$CUTS/campaign.log"
  local start
  start=$(date +%s)
  RAYON_NUM_THREADS=$RAYON IMPG_CUT_LOCAL_DERIVE=1 \
    nice -n 10 stdbuf -oL -eL "$BIN" cut-batch \
    --global-panel "$SYNG" --agc "$AGC" --shards-dir "$SHARDS" \
    --manifest "$B/manifest.jsonl" > "$B/batch.log" 2> "$B/batch.err" &
  local pid=$!; guard_poll "$pid" "$B/batch.rss"; local poll=$!
  local st=0; wait "$pid" || st=$?
  kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$B/batch.exit"
  echo "$(( $(date +%s) - start ))" > "$B/batch.wall"
  for T in "${TODO[@]}"; do
    local D=$RV/map/$T/$T
    if [ -f "$D/cut-stats.json" ]; then
      touch "$RV/prep/$T/cut.done"
    elif [ -f "$D/CUT-FAILED.txt" ]; then
      echo -e "$T\t$(head -c 300 "$D/CUT-FAILED.txt" | tr '\n' ' ')" >> "$CUTS/poison.tsv"
      touch "$RV/prep/$T/cut.poison"
    else
      # batch-level failure or a missing outcome: count the attempt;
      # 3 batch-level failures -> terminal (named in the ledger)
      n=$(cat "$RV/prep/$T/cut.attempts" 2>/dev/null || echo 0)
      n=$((n+1)); echo "$n" > "$RV/prep/$T/cut.attempts"
      if [ "$n" -ge 3 ]; then
        touch "$RV/prep/$T/cut.exit"
        echo "$(date -u +%FT%TZ) $T TERMINAL after $n batch-level failures (last batch exit $st)" >> "$CUTS/batch-failures.log"
      else
        echo "$(date -u +%FT%TZ) $T batch-level failure (batch exit $st, attempt $n)" >> "$CUTS/batch-failures.log"
      fi
    fi
  done
  echo "$(date -u +%FT%TZ) batch $B done exit=$st wall=$(cat "$B/batch.wall")s" >> "$CUTS/campaign.log"
}

BATCH_ID=0
while true; do
  # the pending set: substrate done, not cut, not poisoned
  mapfile -t PENDING < <(while read -r T; do
      [ -f "$RV/prep/$T/substrate.done" ] || continue
      [ -f "$RV/prep/$T/cut.done" ] && continue
      [ -f "$RV/prep/$T/cut.poison" ] && continue
      [ -f "$RV/prep/$T/cut.exit" ] && continue
      echo "$T"
    done < "$RV/prep/accepted.txt")
  NPEND=${#PENDING[@]}
  NDONE=$(while read -r T; do [ -f "$RV/prep/$T/cut.done" ] && echo x; done < "$RV/prep/accepted.txt" | wc -l)
  NTOTAL=$(wc -l < "$RV/prep/accepted.txt")
  if [ "$NPEND" -eq 0 ]; then
    [ "$NDONE" -ge "$NTOTAL" ] && break
    echo "$(date -u +%FT%TZ) waiting for substrates ($NDONE/$NTOTAL cut)" >> "$CUTS/campaign.log"
    sleep 120; continue
  fi
  if [ "$(avail_gb)" -lt "$FLOOR" ]; then
    echo "$(date -u +%FT%TZ) disk floor: $(avail_gb)G < ${FLOOR}G - waiting for the reaper" >> "$CUTS/campaign.log"
    sleep 120; continue
  fi
  # compose up to C batches of W
  pids=()
  idx=0
  while [ "$idx" -lt "$NPEND" ] && [ "${#pids[@]}" -lt "$C" ]; do
    TODO=()
    for ((j=0; j<W && idx+j<NPEND; j++)); do TODO+=("${PENDING[$((idx+j))]}"); done
    one_batch "$CUTS/batch-$(date +%s)-$BATCH_ID" "${TODO[@]}" &
    pids+=("$!")
    BATCH_ID=$((BATCH_ID+1)); idx=$((idx+W))
  done
  for pid in "${pids[@]}"; do wait "$pid"; done
done
NCUT=$(while read -r T; do [ -f "$RV/prep/$T/cut.done" ] && echo x; done < "$RV/prep/accepted.txt" | wc -l)
NPOISON=$(ls "$RV"/prep/*/cut.poison 2>/dev/null | wc -l)
echo "=== cut campaign complete: $NCUT cut, $NPOISON poisoned $(date -u +%FT%TZ)" >> "$CUTS/campaign.log"
touch "$CUTS/campaign.done"
