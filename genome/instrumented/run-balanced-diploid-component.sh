#!/usr/bin/env bash
# Balanced 15x/15x sample under a given diploid ploidy prior, one component.
# Usage: run-balanced-diploid-component.sh chrMT [TAG]
set -uo pipefail
C="$1"; TAG="${2:-balanced}"
source /home/erikg/impg-genome-inference/build-env.sh
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
cd "$D"
P="$D/run-$TAG-$C"
rm -f "$P.done" "$P.exit"
start=$(date +%s)
# The existing diagnostic switch enforces the experiment's GIVEN two-real-slot
# ploidy prior. --depth is per homolog, matching 15x on EACH pure route.
IMPG_FORCE_DIPLOID_DIAGNOSTIC=1 IMPG_MULTIPLICITY_VARIANT=gentle-beta \
  taskset -c "${IMPG_CORES:-0-255}" nice -n 10 \
  stdbuf -oL -eL /home/erikg/impg-genome-inference/target/release/examples/panel_route_mem_routed \
  --panel /home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng \
  --routes /home/erikg/yeast/genome-panel-routes-rebuild-v1/routes \
  --sample "$D/sample/sample.membwt" --reads "$D/reads.fastq.gz" \
  --axis /home/erikg/yeast/genome-mem-bwt-bootstrap-20260911T055203Z/reference-axis-v1.json \
  --bed-directory /home/erikg/yeast/partition-pos64-w10k-d1k-to-completion/results \
  --component "S288C#0#$C" --depth 15 --routing-universe component --spine --spine-phasing \
  --reference-pair "truth:$D/routes/S288C-$C.json:$D/routes/SK1-$C.json" \
  --reference-pair "native2:$D/routes/SK1-$C.json:$D/routes/SK1-$C.json" \
  --reference-pair "haploid_truth:$D/routes/S288C-$C.json:$D/routes/empty.json" \
  > "$P.log" 2> "$P.err" &
pid=$!
( while kill -0 "$pid" 2>/dev/null; do
    rss_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
    echo "$(date -u +%FT%TZ) $pid ${rss_kb:-0}" >> "$P.rss"
    sleep 5
  done ) & poller=$!
( while kill -0 "$pid" 2>/dev/null; do
    ticks=$(awk '{print $14+$15}' /proc/$pid/stat 2>/dev/null) || continue
    stage=$(grep -a -oE '^\[[^]]+\][^:]{0,40}' "$P.err" 2>/dev/null | tail -1)
    err_age=$(( $(date +%s) - $(stat -c %Y "$P.err" 2>/dev/null || date +%s) ))
    echo "$(date -u +%FT%TZ) cpu_s=$(( ${ticks:-0} / 100 )) stage=$stage (last-flushed ${err_age}s ago)" >> "$P.stages"
    sleep 30
  done ) & timer=$!
status=0; wait "$pid" || status=$?
kill "$poller" "$timer" 2>/dev/null; wait "$poller" "$timer" 2>/dev/null
echo "$status" > "$P.exit"
echo "$(( $(date +%s) - start ))" > "$P.wall"
touch "$P.done"
exit "$status"
