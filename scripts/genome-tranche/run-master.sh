#!/usr/bin/env bash
# THE GENOME TRANCHE - the master driver (marker-gated stages, sequenced;
# every stage's own script is resumable, so this loop is killable and
# restartable at any point):
#   P  the prep mirror + regen pool (all 56 substrates)
#   R  the reads fetch pool
#   K  the RV local-derive cuts (the restriction basis; W=8 batches, C=2)
#   F  the fleet pipeline run (fork driver, profile-B cuts + the
#      4-concurrent T16 in-child production calls)
#   A  the profile-A restricted cuts (4-way, the committed single-cut
#      shape; the only committed restricted path)
#   TA the official truth assessment, profile A (serial lanes)
#   TB the official truth assessment, profile B
# The co-tenant state is receipted per stage in co-tenant.log (another
# worker's gate fleet has been on the box since before the tranche
# opened - every wall below is a CONTENDED wall).
set -uo pipefail
WT=${WT:-/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1}
cd "$WT"
N=$(grep -v '^#' tranche.tsv | wc -l)
NF=$( [ -f tranche-failed.tsv ] && wc -l < tranche-failed.tsv || echo 0 )
NOK=$(( N - NF ))

receipt_box() {
  echo "$(date -u +%FT%TZ) load=$(cut -d' ' -f1 /proc/loadavg) memavail=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)GiB df=$(df -BG --output=avail / | tail -1 | tr -dc '0-9')G $1" >> "$WT/logs/co-tenant.log"
}

stage_done() { [ -f "$WT/logs/stage-$1.done" ]; }
mark()       { touch "$WT/logs/stage-$1.done"; receipt_box "STAGE $1 COMPLETE"; }

receipt_box "OPEN tranche n=$N"

# P: prep
if ! stage_done P; then
  POOL_W=8 THREADS=4 bash run-prep-pool.sh >> "$WT/logs/stage-P.log" 2>&1
  nd=$(ls prep/*/substrate.done 2>/dev/null | wc -l)
  [ "$nd" -eq "$NOK" ] && mark P || echo "$(date -u +%FT%TZ) P incomplete ($nd/$N)" >> "$WT/logs/master.log"
fi
# R: reads
if ! stage_done R; then
  bash run-reads.sh >> "$WT/logs/stage-R.log" 2>&1
  nd=$(ls prep/*/reads.done 2>/dev/null | wc -l)
  [ "$nd" -eq "$NOK" ] && mark R || echo "$(date -u +%FT%TZ) R incomplete ($nd/$N)" >> "$WT/logs/master.log"
fi
# K: the restriction-basis cuts (loop until all done or no progress)
if ! stage_done K; then
  while true; do
    bash run-rvcut.sh >> "$WT/logs/stage-K.log" 2>&1
    nd=$(ls prep/*/rvcut.done 2>/dev/null | wc -l)
    [ "$nd" -eq "$NOK" ] && { mark K; break; }
    pend=$(grep -v '^#' tranche.tsv | cut -f1 | while read -r T; do
        [ -f "prep/$T/rvcut.done" ] && continue
        [ -f "prep/$T/substrate.done" ] && echo x; done | wc -l)
    [ "$pend" -eq 0 ] && { echo "$(date -u +%FT%TZ) K stopped ($nd done, 0 pending)" >> "$WT/logs/master.log"; break; }
  done
fi
# F: the fleet pipeline (profile-B cuts + production calls) - needs ALL prep
if ! stage_done F && [ "$(ls prep/*/substrate.done 2>/dev/null | wc -l)" -eq "$NOK" ]; then
  bash run-fleetB.sh >> "$WT/logs/stage-F.log" 2>&1
  nd=$(ls prep/*/cutB.done 2>/dev/null | wc -l)
  [ "$nd" -eq "$NOK" ] && mark F || echo "$(date -u +%FT%TZ) F incomplete ($nd/$N)" >> "$WT/logs/master.log"
fi
# A: the profile-A restricted cuts
if ! stage_done A && stage_done K; then
  POOL_W=4 RAYON=32 bash run-cutA.sh >> "$WT/logs/stage-A.log" 2>&1
  nd=$(ls prep/*/cutA.done 2>/dev/null | wc -l)
  [ "$nd" -eq "$NOK" ] && mark A || echo "$(date -u +%FT%TZ) A incomplete ($nd/$N)" >> "$WT/logs/master.log"
fi
# TA / TB: the official truth assessment
if stage_done A || [ -n "$(ls prep/*/cutA.done 2>/dev/null)" ]; then
  if ! stage_done TA; then
    LANES=6 bash run-truthcalls.sh A >> "$WT/logs/stage-TA.log" 2>&1
    nd=$(ls prep/*/callA.done 2>/dev/null | wc -l)
    [ "$nd" -eq "$NOK" ] && mark TA || echo "$(date -u +%FT%TZ) TA incomplete ($nd/$N)" >> "$WT/logs/master.log"
  fi
fi
if [ -n "$(ls prep/*/cutB.done 2>/dev/null)" ]; then
  if ! stage_done TB; then
    LANES=6 bash run-truthcalls.sh B >> "$WT/logs/stage-TB.log" 2>&1
    nd=$(ls prep/*/callB.done 2>/dev/null | wc -l)
    [ "$nd" -eq "$NOK" ] && mark TB || echo "$(date -u +%FT%TZ) TB incomplete ($nd/$N)" >> "$WT/logs/master.log"
  fi
fi
receipt_box "CLOSE"
echo "$(date -u +%FT%TZ) master driver exiting" >> "$WT/logs/master.log"
