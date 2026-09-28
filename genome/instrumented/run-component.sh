#!/usr/bin/env bash
# genome-genotyping-v1-item1: one component's campaign run under the COLLAPSED
# pipeline (STEPS 1-4b: canonical anchors + streamed bounded routing + admissible
# pruning + the folded DP + posterior). The SAME invocation as the main lane.
# Usage: run-component.sh <S288C#0#chrX> [LANE_TAG]
set -uo pipefail
COMP="$1"; TAG="${2:-main}"
source /home/erikg/impg-genome-inference/build-env.sh
D=/home/erikg/impg/target/experiments/genome-mem-bwt-pipeline/genome-genotyping-v1-item1
cd "$D"
status=0
P="run-$TAG-$COMP"
rm -f "$P.exit" "$P.done" "$P.rss"
start=$(date +%s)
IMPG_MULTIPLICITY_VARIANT=gentle-beta taskset -c 240-243 nice -n 10 \
  stdbuf -oL -eL \
  /home/erikg/impg-genome-inference/target/release/examples/panel_route_mem_routed \
  --panel /home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng \
  --routes /home/erikg/yeast/genome-panel-routes-rebuild-v1/routes \
  --sample /home/erikg/yeast/genome-blind-mosaic-20260914T121834Z/sample/sample.membwt \
  --reads /home/erikg/yeast/genome-blind-mosaic-20260914T121834Z/reads.fastq.gz \
  --axis /home/erikg/yeast/genome-mem-bwt-bootstrap-20260911T055203Z/reference-axis-v1.json \
  --bed-directory /home/erikg/yeast/partition-pos64-w10k-d1k-to-completion/results \
  --component "$COMP" \
  --routing-universe component --spine --spine-phasing \
  --reference-pair "truth:$D/routes/truth-$COMP.json:$D/routes/native-$COMP.json" \
  --reference-pair "native2:$D/routes/native-$COMP.json:$D/routes/native-$COMP.json" \
  --reference-pair "haploid_truth:$D/routes/truth-$COMP.json:$D/routes/empty.json" \
  > "$P.log" 2> "$P.err" &
pid=$!
( while kill -0 "$pid" 2>/dev/null; do
    rss_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
    echo "$(date -u +%FT%TZ) $pid ${rss_kb:-0}" >> "$P.rss"
    sleep 5
  done ) &
poller=$!
# STAGE TIMER + CPU HEARTBEAT (the watchdog contract: every 30 s, the
# run's CPU seconds and its LATEST STAGE MARKER land in $P.stages — a
# watchdog sees progress even when a single stage runs for hours; the
# RSS poller above carries the memory side at 5 s resolution).
( while kill -0 "$pid" 2>/dev/null; do
    ticks=$(awk '{print $14+$15}' /proc/$pid/stat 2>/dev/null) || continue
    cpu_s=$(( ${ticks:-0} / 100 ))
    stage=$(grep -a -oE '^\[[^]]+\][^:]{0,40}' "$P.err" 2>/dev/null | tail -1)
    err_age=$(( $(date +%s) - $(stat -c %Y "$P.err" 2>/dev/null || date +%s) ))
    echo "$(date -u +%FT%TZ) cpu_s=$cpu_s stage=$stage (last-flushed ${err_age}s ago)" >> "$P.stages"
    sleep 30
  done ) &
timer=$!
wait "$pid" || status=$?
kill "$poller" 2>/dev/null; wait "$poller" 2>/dev/null
kill "$timer" 2>/dev/null; wait "$timer" 2>/dev/null
echo "$status" > "$P.exit"
echo "$(( $(date +%s) - start ))" > "$P.wall"
touch "$P.done"
