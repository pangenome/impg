#!/usr/bin/env bash
# THE WHOLE-GENOME RUNTIME VERIFICATION - STAGE 1: THE QUERY (the committed
# stage-B recipe verbatim: `impg query -d 10000` over the 9,891 bridged
# windows, sharded N ways - each process ONE global syng + AGC load
# amortized over its shard's windows). Marker-gated per shard; resumable.
# usage: run-query.sh [N_SHARDS]   (env: RV, IMPG, SYNG, AGC)
set -uo pipefail
RV=${RV:-/home/erikg/hprcv2/wgs35x-fleet/runtime-verification}
IMPG=${IMPG:-/home/erikg/impg-genome-inference/target/release/impg}
SYNG=${SYNG:-/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng}
AGC=${AGC:-/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc}
N=${1:-14}
Q=$RV/query
mkdir -p "$Q"

# 1. shard the windows round-robin (chromosome-balanced by construction)
if [ ! -f "$Q/shards.done" ]; then
  rm -f "$Q"/shard-*.bed
  python3 - "$Q" "$RV/query-input.bed" "$N" <<'EOF'
import sys
Q, inp, n = sys.argv[1:4]
n = int(n)
files = [open(f"{Q}/shard-{i:02d}.bed", "w") for i in range(n)]
for i, line in enumerate(open(inp)):
    if line.strip():
        files[i % n].write(line)
for f in files: f.close()
EOF
  touch "$Q/shards.done"
fi

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

one_shard() {
  local S=$1
  [ -f "$Q/shard-$S.done" ] && return 0
  rm -f "$Q/shard-$S.out.bed"
  local start
  start=$(date +%s)
  nice -n 10 stdbuf -oL -eL "$IMPG" query -a "$SYNG" \
    -b "$Q/shard-$S.bed" -d 10000 --sequence-files "$AGC" \
    -o bed -t 8 > "$Q/shard-$S.out.bed" 2> "$Q/shard-$S.err" &
  local pid=$!; guard_poll "$pid" "$Q/shard-$S.rss"; local poll=$!
  local st=0; wait "$pid" || st=$?
  kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$Q/shard-$S.exit"
  echo "$(( $(date +%s) - start ))" > "$Q/shard-$S.wall"
  # completion check: exit 0 AND every input window present in the output
  if [ "$st" -eq 0 ]; then
    local n_in n_out
    n_in=$(grep -c . "$Q/shard-$S.bed")
    n_out=$(cut -f4 "$Q/shard-$S.out.bed" | sort -u | wc -l)
    if [ "$n_in" -eq "$n_out" ]; then
      touch "$Q/shard-$S.done"
    else
      echo "shard-$S: $n_out/$n windows answered" > "$Q/shard-$S.short"
    fi
  fi
}

pids=()
for i in $(seq 0 $((N-1))); do
  one_shard "$(printf "%02d" $i)" &
  pids+=("$!")
done
for pid in "${pids[@]}"; do wait "$pid"; done
NDONE=$(ls "$Q"/shard-*.done 2>/dev/null | wc -l)
echo "$NDONE" > "$Q/queries.done.count"
[ "$NDONE" -eq "$N" ] && touch "$Q/queries.done"
echo "=== query stage: $NDONE/$N shards done $(date -u +%FT%TZ)"
