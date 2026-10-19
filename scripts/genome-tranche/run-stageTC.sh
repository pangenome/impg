#!/usr/bin/env bash
# THE GENOME TRANCHE - stage TC (the steering's third form, run as a
# STANDALONE loop beside the master so the running master is never
# edited mid-flight): the NEW PRODUCT DEFAULT calls - the restricted
# census (form A's own cuts, reused; the prune is call-side) with
# --far-bin-prune, assessed with the official QV machinery. Trails
# stage A (form C needs cutA); exits when every runnable member has
# callC.done (or stage A is complete with nothing left).
# usage: run-stageTC.sh   (env: WT, LANES=6)
set -uo pipefail
WT=${WT:-/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1}
cd "$WT"
while true; do
  LANES=${LANES:-6} bash run-truthcalls.sh C >> "$WT/logs/stage-TC.log" 2>&1
  nd=$(ls prep/*/callC.done 2>/dev/null | wc -l)
  if [ "$nd" -ge "$(( $(grep -v '^#' tranche.tsv | wc -l) - $( [ -f tranche-failed.tsv ] && wc -l < tranche-failed.tsv || echo 0 ) ))" ]; then
    echo "$(date -u +%FT%TZ) stage TC complete ($nd)" >> "$WT/logs/stage-TC.log"
    touch "$WT/logs/stage-TC.done"; exit 0
  fi
  if [ -f "$WT/logs/stage-A.done" ]; then
    # stage A closed: whatever has no callC now is a cutA-failure class,
    # named by the stage-A logs; TC is done with what landed
    echo "$(date -u +%FT%TZ) stage TC complete at stage-A close ($nd)" >> "$WT/logs/stage-TC.log"
    touch "$WT/logs/stage-TC.done"; exit 0
  fi
  sleep 120
done
