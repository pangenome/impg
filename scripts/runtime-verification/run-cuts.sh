#!/usr/bin/env bash
# THE WHOLE-GENOME RUNTIME VERIFICATION - STAGE 3: THE CUT CAMPAIGN (the
# sample-side map step): sample_projection cut-batch from the COMMITTED 35x
# projection shards - batches of W territories, ONE global syng + AGC load +
# ONE shard scan per batch (the amortized-load shape; the per-territory
# fail-closed bind failures ISOLATED by the patched emissions loop - the
# poison ledger names them). C concurrent batches; the disk floor gates every
# launch (the streaming reaper must keep up). Marker-gated, resumable.
# THE WIDTH vs THE 64GiB GUARD (measured 2026-10-15): W=32 batches were
# SIGKILLed by the per-process 64GiB guard at 65-69GB RSS (4/4 batches, 0
# territories cut) - the batch memory grows with the batch's aggregate
# routing-index nodes (779k nodes -> ~63GB at scan end; the dense batch
# 1792039134-3 carried 1.76M nodes). THE GUARD IS LAW; the width yields:
# BATCH_W=16. Two attribution rules fixed with it: (a) a guard kill
# (exit 137) is a BATCH-PACKAGING event, never a per-territory strike - the
# 54 false cut.exit terminals from the W=32 era were reset (poison.tsv
# empty: zero intrinsic territory failures); (b) the pending list is
# deterministically hash-interleaved before batch composition so genome-
# adjacent density clumps spread across batches (per-batch memory ~ the
# genome average, not the local clump).
# usage: run-cuts.sh   (env: RV, BATCH_W=16, N_CONCURRENT=2, RAYON=48,
#                       DISK_FLOOR_GB=1500)
set -uo pipefail
RV=${RV:-/home/erikg/hprcv2/wgs35x-fleet/runtime-verification}
BIN=/home/erikg/impg-genome-inference/target/release/examples/sample_projection
SYNG=/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng
AGC=/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc
SHARDS=/home/erikg/hprcv2/hg002-projection/full
W=${BATCH_W:-16}
C=${N_CONCURRENT:-2}
RAYON=${RAYON:-48}
FLOOR=${DISK_FLOOR_GB:-1500}
CUTS=$RV/cuts
mkdir -p "$CUTS"

avail_gb() { df -BG --output=avail / | tail -1 | tr -dc '0-9'; }

guard_poll() {
  local pid=$1 rssf=$2
  ( while kill -0 "$pid" 2>/dev/null; do
      rss_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
      echo "$(date -u +%FT%TZ) ${rss_kb:-0}" >> "$rssf"
      if [ "${rss_kb:-0}" -gt $((64*1024*1024)) ]; then
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
      if [ "$st" -eq 137 ]; then
        # THE GUARD KILL IS A BATCH-PACKAGING EVENT, never a territory
        # strike (the width's own failure, not the territory's)
        echo "$(date -u +%FT%TZ) $T batch GUARD-KILLED (exit 137, packaging - no strike)" >> "$CUTS/batch-failures.log"
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
    fi
  done
  echo "$(date -u +%FT%TZ) batch $B done exit=$st wall=$(cat "$B/batch.wall")s" >> "$CUTS/campaign.log"
}

Q=$RV/query
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
  # THE DETERMINISTIC HASH INTERLEAVE: genome-adjacent territories carry
  # correlated density (the segdup clumps); shuffling the pending list by a
  # fixed hash makes every batch a uniform sample of the density profile -
  # per-batch memory ~ the genome average, never the local clump
  mapfile -t PENDING < <(printf '%s\n' "${PENDING[@]}" | while read -r T; do
      printf '%s %s\n' "$(printf '%s' "$T" | md5sum | cut -c1-12)" "$T"
    done | sort | awk '{print $2}')
  NPEND=${#PENDING[@]}
  NDONE=$(while read -r T; do [ -f "$RV/prep/$T/cut.done" ] && echo x; done < "$RV/prep/accepted.txt" | wc -l)
  NTOTAL=$(wc -l < "$RV/prep/accepted.txt")
  if [ "$NPEND" -eq 0 ]; then
    # terminal only when the upstream is terminal too: the query stage done,
    # the final reconciliation split done, and every accepted territory has a
    # terminal cut state (the accepted set grows while the query runs)
    if [ -f "$Q/queries.done" ] && [ -f "$RV/prep/final-split.done" ]; then
      NTERM=$((NDONE + $(ls "$RV"/prep/*/cut.poison 2>/dev/null | wc -l) + $(ls "$RV"/prep/*/cut.exit 2>/dev/null | wc -l)))
      [ "$NTERM" -ge "$NTOTAL" ] && break
    fi
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
