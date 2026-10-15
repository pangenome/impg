#!/usr/bin/env bash
# Context-domain p1 identity gate under the same 4-core baseline configuration.
set -uo pipefail
source /home/erikg/impg-genome-inference/build-env.sh
D=/home/erikg/impg/target/experiments/genome-mem-bwt-pipeline/context-domains-v1
cd "$D"
P=run-interface-p1
rm -f "$P.done" "$P.exit"
A=/home/erikg/impg/target/experiments/genome-mem-bwt-pipeline/chrIII-assessment-v1/routes
EMPTY=/home/erikg/impg/target/experiments/genome-mem-bwt-pipeline/chrIII-multiplicity-background-v1/empty-route.json
start=$(date +%s)
IMPG_MULTIPLICITY_VARIANT=gentle-beta taskset -c 240-243 nice -n 10 \
  stdbuf -oL -eL /home/erikg/impg-genome-inference/target/release/examples/panel_route_mem_routed \
  --panel /home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng \
  --routes /home/erikg/yeast/genome-panel-routes-rebuild-v1/routes \
  --sample /home/erikg/yeast/genome-blind-mosaic-20260914T121834Z/sample/sample.membwt \
  --reads /home/erikg/yeast/genome-blind-mosaic-20260914T121834Z/reads.fastq.gz \
  --axis /home/erikg/yeast/genome-mem-bwt-bootstrap-20260911T055203Z/reference-axis-v1.json \
  --bed-directory /home/erikg/yeast/partition-pos64-w10k-d1k-to-completion/results \
  --component S288C#0#chrIII --routing-universe component --spine --spine-phasing \
  --context-aware-domains --locus-range 12 27 \
  --reference-pair "truth:$A/truth-mosaic-chrIII.json:$A/native-chrIII.json" \
  --reference-pair "native2:$A/native-chrIII.json:$A/native-chrIII.json" \
  --reference-pair "haploid_truth:$A/truth-mosaic-chrIII.json:$EMPTY" \
  > "$P.log" 2> "$P.err" &
pid=$!
( while kill -0 "$pid" 2>/dev/null; do
    echo "$(date -u +%FT%TZ) $pid $(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null)" >> "$P.rss"
    sleep 5
  done ) & poller=$!
( while kill -0 "$pid" 2>/dev/null; do
    ticks=$(awk '{print $14+$15}' /proc/$pid/stat 2>/dev/null) || continue
    stage=$(grep -a -oE '^\[[^]]+\][^:]{0,40}' "$P.err" 2>/dev/null | tail -1)
    echo "$(date -u +%FT%TZ) cpu_s=$(( ${ticks:-0} / 100 )) stage=$stage" >> "$P.stages"
    sleep 30
  done ) & timer=$!
status=0; wait "$pid" || status=$?
kill "$poller" "$timer" 2>/dev/null; wait "$poller" "$timer" 2>/dev/null
echo "$status" > "$P.exit"
echo "$(( $(date +%s) - start ))" > "$P.wall"
touch "$P.done"
exit "$status"
