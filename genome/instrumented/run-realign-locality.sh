#!/usr/bin/env bash
# THE LOCALITY-DOMAIN REALIGNMENT RUN, one component (the gate-2
# runner of the alignment-induced locality domain). Identical to
# run-realign-run.sh except: --partition-graphs points at the
# locality map dir (locality-domain-maps.py's emission, the window-
# keyed locality folds), and IMPG_REALIGN_LOCALITY_DOMAIN=1 switches
# the record/occurrence filters to the territory-touch convention on
# the locality extents. The machinery (territory images, marginal
# realignment likelihood, DP, QUAL) stands; the receipts carry the
# locality model string.
# Usage: run-realign-locality.sh WIDTH|serial COMP "0,1,..." OUTPREFIX [TAG]
set -uo pipefail
W="$1"; C="$2"; LOCI="$3"; OUT="$4"; TAG="${5:-realignlocality-${W}-$C}"
source /home/erikg/impg-genome-inference/build-env.sh
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
cd "$D" || exit 1
P="$D/run-$TAG"
rm -f "$P.done" "$P.exit" "$P.rss" "$P.stages"
start=$(date +%s)
SEAMENV=()
if [ "$W" != "serial" ]; then
  SEAMENV=(IMPG_REALIGN_PARALLEL_SEAMS=1 IMPG_REALIGN_SEAM_WIDTH="$W")
fi
IMPG_CORES="${IMPG_CORES:-0-255}" taskset -c "${IMPG_CORES:-0-255}" nice -n 10 \
  env "${SEAMENV[@]}" IMPG_REALIGN_LOCALITY_DOMAIN=1 \
  stdbuf -oL -eL /home/erikg/impg-genome-inference/target/release/examples/partition_realign_score \
  --panel /home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng \
  --routes /home/erikg/yeast/genome-panel-routes-rebuild-v1/routes \
  --partition-graphs "$D/locality-graphs-$C" \
  --component "S288C#0#$C" \
  --census "$D/cosine-diagnostic-scratch/$C/cosine-multi-census.jsonl" \
  --reads "$D/reads.fastq.gz" \
  --derive-cache "$D/realign-derive-cache-chrMT.bin" \
  --loci "$LOCI" \
  --out "$D/$OUT.jsonl" \
  --exactness-out "$D/$OUT.exactness.jsonl" \
  --skeleton-out "$D/$OUT.skeleton.jsonl" \
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
