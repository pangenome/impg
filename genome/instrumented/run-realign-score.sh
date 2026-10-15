#!/usr/bin/env bash
# The realignment scoring layer (slice B), one component (assessment-side).
# Usage: run-realign-score.sh chrI "2,4,7" [TAG]
# Supplies the house runner's markers, external RSS poller and 64GiB guard.
set -uo pipefail
C="$1"; LOCI="$2"; TAG="${3:-realignscore}"
source /home/erikg/impg-genome-inference/build-env.sh
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
cd "$D" || exit 1
P="$D/run-$TAG-$C"
rm -f "$P.done" "$P.exit"
start=$(date +%s)
IMPG_CORES="${IMPG_CORES:-0-255}" taskset -c "${IMPG_CORES:-0-255}" nice -n 10 \
  stdbuf -oL -eL /home/erikg/impg-genome-inference/target/release/examples/partition_realign_score \
  --panel /home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng \
  --routes /home/erikg/yeast/genome-panel-routes-rebuild-v1/routes \
  --partition-graphs "$D/partition-graphs" \
  --component "S288C#0#$C" \
  --census "$D/cosine-diagnostic-scratch/$C/cosine-multi-census.jsonl" \
  --reads "$D/reads.fastq.gz" \
  --derive-cache "$D/realign-derive-cache-chrMT.bin" \
  --loci "$LOCI" \
  --out "$D/realign-score-$C.jsonl" \
  --exactness-out "$D/realign-score-$C.exactness.jsonl" \
  > "$P.log" 2> "$P.err" &
pid=$!
( while kill -0 "$pid" 2>/dev/null; do
    rss_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
    echo "$(date -u +%FT%TZ) $pid ${rss_kb:-0}" >> "$P.rss"
    sleep 5
  done ) & poller=$!
( while kill -0 "$pid" 2>/dev/null; do
    ticks=$(awk '{print $14+$15}' /proc/$pid/stat 2>/dev/null) || continue
    stage=$(grep -a -oE '^\[[^]]+\][^:]{0,60}' "$P.err" 2>/dev/null | tail -1)
    echo "$(date -u +%FT%TZ) cpu_s=$(( ${ticks:-0} / 100 )) stage=$stage" >> "$P.stages"
    sleep 30
  done ) & timer=$!
status=0; wait "$pid" || status=$?
kill "$poller" "$timer" 2>/dev/null; wait "$poller" "$timer" 2>/dev/null
echo "$status" > "$P.exit"
echo "$(( $(date +%s) - start ))" > "$P.wall"
touch "$P.done"
exit "$status"
