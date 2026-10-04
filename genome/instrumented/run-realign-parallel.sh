#!/usr/bin/env bash
# THE PARALLELISM PILOT rung (phase 3 of the runtime plan): the
# realignment instrument with the seam-parallel body (seam 1 = the
# per-record binding loop, seam 2 = the per-locus loop) at a given
# rayon pool width, one component. THE RACE DETECTOR IS THE IDENTITY
# GATE: the receipts and all four sidecars must be byte-identical to
# the serial receipts on every semantic field (check-realign-scoring.py
# --exhaustive <comp> --timered-base <base> --run-tag <tag>, phase 8t).
# Usage: run-realign-parallel.sh WIDTH|serial chrI "0,1,...,20" [TAG]
# WIDTH "serial" runs the seam code on the committed serial path (no
# env vars) — the refactor's own identity re-baseline; WIDTH N sets
# RAYON_NUM_THREADS=N with IMPG_REALIGN_PARALLEL_SEAMS=1 (the ladder
# rung). Supplies the house runner's markers, external RSS poller and
# 64GiB guard; the box's other work respected (nice 10).
set -uo pipefail
W="$1"; C="$2"; LOCI="$3"; TAG="${4:-realignpar-${W}-$C-$C}"
source /home/erikg/impg-genome-inference/build-env.sh
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
cd "$D" || exit 1
BASE="realign-par-${W}-$C"
P="$D/run-$TAG"
rm -f "$P.done" "$P.exit"
start=$(date +%s)
SEAMENV=()
if [ "$W" != "serial" ]; then
  SEAMENV=(IMPG_REALIGN_PARALLEL_SEAMS=1 IMPG_REALIGN_SEAM_WIDTH="$W")
fi
IMPG_CORES="${IMPG_CORES:-0-255}" taskset -c "${IMPG_CORES:-0-255}" nice -n 10 \
  env "${SEAMENV[@]}" \
  stdbuf -oL -eL /home/erikg/impg-genome-inference/target/release/examples/partition_realign_score \
  --panel /home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng \
  --routes /home/erikg/yeast/genome-panel-routes-rebuild-v1/routes \
  --partition-graphs "$D/partition-graphs" \
  --component "S288C#0#$C" \
  --census "$D/cosine-diagnostic-scratch/$C/cosine-multi-census.jsonl" \
  --reads "$D/reads.fastq.gz" \
  --derive-cache "$D/realign-derive-cache-chrMT.bin" \
  --loci "$LOCI" \
  --out "$D/$BASE.jsonl" \
  --exactness-out "$D/$BASE.exactness.jsonl" \
  --skeleton-out "$D/$BASE.skeleton.jsonl" \
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
