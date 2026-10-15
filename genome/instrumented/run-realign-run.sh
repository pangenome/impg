#!/usr/bin/env bash
# THE REALIGNMENT SCORING RUN, one component, an arbitrary receipt
# prefix and seam width (the chrIV end-to-end runner + the identity
# gate runner for the window-generalization refactor). The seam
# switch and width follow the phase-3 corrected mechanism
# (IMPG_REALIGN_PARALLEL_SEAMS + IMPG_REALIGN_SEAM_WIDTH; "serial"
# runs the committed serial path). THE RACE DETECTOR IS THE IDENTITY
# GATE: receipts and all four sidecars must be byte-identical to the
# serial receipts on every semantic field.
# Usage: run-realign-run.sh WIDTH|serial COMP "0,1,..." OUTPREFIX [TAG]
# Supplies the house runner's markers, external RSS poller and 64GiB guard.
set -uo pipefail
W="$1"; C="$2"; LOCI="$3"; OUT="$4"; TAG="${5:-realignrun-${W}-$C}"
source /home/erikg/impg-genome-inference/build-env.sh
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
cd "$D" || exit 1
P="$D/run-$TAG"
# THE PER-RUN MARKERS ARE PER-RUN: a re-run regenerates them, so the
# stale previous run's rss/stages history must not leak into the new
# run's record — the checker's phase-1 guard reads the whole .rss file
# and the fleet's dead 4-wides (exit 1, .done touched anyway) left
# 78-118GB peaks in the logs that flagged the CLEAN re-runs. Preserve
# nothing here: copy aside before relaunching if the history matters.
rm -f "$P.done" "$P.exit" "$P.rss" "$P.stages"
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
