#!/usr/bin/env bash
# THE MULTI-MATCHING CENSUS run, one component (the recipe of the
# committed chrMT/chrI multicensus pilots, re-runnable per component:
# the balanced 15x/15x component run under the read-matched
# observation gate with the census emission on; the likelihood and
# admission receipts are the component's read-matched BEFORE-record,
# written to the diagnostic scratch so the preserved balanced
# sidecars at the validation root are never touched).
# Usage: run-cosine-multicensus.sh chrIV
# Supplies the house runner's markers, external RSS poller, stage
# timer; the instrument's own resident-set guard stays armed.
set -uo pipefail
C="$1"
source /home/erikg/impg-genome-inference/build-env.sh
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
SCRATCH="$D/cosine-diagnostic-scratch/$C"
mkdir -p "$SCRATCH"
cd "$SCRATCH" || exit 1
P="$D/run-multicensus-pilot-$C"
rm -f "$P.done" "$P.exit"
start=$(date +%s)
# The existing diagnostic switch enforces the experiment's GIVEN
# two-real-slot ploidy prior. --depth is per homolog, matching 15x on
# EACH pure route. The census is defined under the read-matched
# (option-b) convention; the likelihood/admission receipts require
# the pure-material domain gate, as in the committed pilots.
IMPG_FORCE_DIPLOID_DIAGNOSTIC=1 IMPG_MULTIPLICITY_VARIANT=gentle-beta \
IMPG_COSINE_READ_MATCHED=1 \
IMPG_COSINE_MULTI_CENSUS="$SCRATCH/cosine-multi-census.jsonl" \
IMPG_COSINE_DOMAIN_PILOT=1 \
IMPG_COSINE_LIKELIHOOD_OUTPUT="$SCRATCH/cosine-graph-likelihood-readmatched-$C.jsonl" \
IMPG_COSINE_ADMISSION_DIAGNOSTIC="$SCRATCH/cosine-admission-diagnostic-readmatched-$C.jsonl" \
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
    stage=$(grep -a -oE '^\[[^]]+][^:]{0,40}' "$P.err" 2>/dev/null | tail -1)
    echo "$(date -u +%FT%TZ) cpu_s=$(( ${ticks:-0} / 100 )) stage=$stage" >> "$P.stages"
    sleep 30
  done ) & timer=$!
status=0; wait "$pid" || status=$?
kill "$poller" "$timer" 2>/dev/null; wait "$poller" "$timer" 2>/dev/null
echo "$status" > "$P.exit"
echo "$(( $(date +%s) - start ))" > "$P.wall"
touch "$P.done"
exit "$status"
