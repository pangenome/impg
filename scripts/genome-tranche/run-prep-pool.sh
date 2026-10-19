#!/usr/bin/env bash
# THE GENOME TRANCHE - stage 1 POOL driver: the prep mirror + regen for
# every tranche territory (POOL_W chains, marker-gated per territory via
# the chain's own .done files; the segdup-tail beds steps are the honest
# straggler risk - the pool drains ordinary members first by launching
# in tranche order, the walls receipted per step).
# usage: run-prep-pool.sh   (env: WT, POOL_W=8, THREADS=4)
set -uo pipefail
WT=${WT:-/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1}
POOL_W=${POOL_W:-8}
export WT
mkdir -p "$WT/logs"
grep -v '^#' "$WT/tranche.tsv" | cut -f1 | \
  xargs -P "$POOL_W" -I{} bash -c '
    T={}
    RV=$WT bash "$WT/run-prep-mirror.sh" "$T" > "$WT/logs/prep-$T.log" 2>&1 \
      && touch "$WT/prep/$T/mirror.done" \
      || echo "$T prep FAILED (see logs/prep-$T.log)" >> "$WT/logs/prep-failures.txt"
  '
n=$(grep -v '^#' "$WT/tranche.tsv" | cut -f1 | while read -r T; do
      [ -f "$WT/prep/$T/substrate.done" ] && echo x; done | wc -l)
echo "$(date -u +%FT%TZ) prep pool complete: $n substrate.done of $(grep -v '^#' "$WT/tranche.tsv" | wc -l)"
