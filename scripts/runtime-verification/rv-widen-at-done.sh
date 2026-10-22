#!/usr/bin/env bash
# THE WIDEN WATCHER (detached; rides the campaign through the queries.done
# transition): waits for the query stage's queries.done + the prep master's
# final reconciliation split, then widens the pipeline for the clean
# post-query box - THE WIDENING AXIS IS CONCURRENCY, NEVER THE BATCH WIDTH
# (the 64GiB guard caps the width; W=8 measured: peaks 64.4-65.3GB).
# THE KILL IS THE GROUP KILL (every master was setsid'd: pgid == pid; the
# lessons of this session: killing a master's children first hands it a 1s
# window to compose replacements; killing only the master orphans its pools
# mid-flight - the group kill takes master + pools + stage children at once,
# and every stage is marker-gated, so the relaunched master resumes exactly
# where the markers stand: cut batches die unaccounted (no false strikes),
# calls die before call.done (clean retry, no terminal call.exit), prep
# chains die mid-step (the orphan sweep then marks/cleans the partials)).
# usage: rv-widen-at-done.sh  (env: RV, WIDEN_MARK=widen-done)
set -uo pipefail
RV=${RV:-/home/erikg/hprcv2/wgs35x-fleet/runtime-verification}
MARK=$RV/${WIDEN_MARK:-widen-done}
[ -f "$MARK" ] && exit 0
log() { echo "$(date -u +%FT%TZ) $*" >> "$RV/widen.log"; }
log "widen watcher armed: waiting for queries.done + prep/final-split.done"

master_pid() { ps -eo pid,args | awk -v s="$1" '$2=="bash" && $3 ~ ("^\\./" s "$|^" s "$") {print $1; exit}'; }
group_kill() { # group_kill <script> : TERM the whole process group of its master
  local pid; pid=$(master_pid "$1")
  if [ -n "${pid:-}" ]; then
    kill -TERM -- "-$pid" 2>/dev/null && log "group-killed $1 master pid $pid" \
      || { kill -TERM "$pid" && log "killed $1 master pid $pid (no group)"; }
  else
    log "$1 master not found (already gone)"
  fi
}

# 1. the upstream terminal conditions
while [ ! -f "$RV/query/queries.done" ]; do sleep 60; done
log "queries.done landed"
while [ ! -f "$RV/prep/final-split.done" ]; do sleep 60; done
log "final-split.done landed (the bridge-rule reconciliation is the ledger)"

# 2. THE PREP POOL (the beds segdup tail is the post-query pole: measured
#    mean 385s/territory at -t4, p90 1,262s, max 5,672s): POOL_W 32 -> 64
group_kill run-prep-master.sh
sleep 20
python3 "$RV/rv-orphan-sweep.py" "$RV" >> "$RV/orphan-watcher.log" 2>&1
sleep 10
python3 "$RV/rv-orphan-sweep.py" "$RV" >> "$RV/orphan-watcher.log" 2>&1
setsid env POOL_W=64 THREADS=4 bash "$RV/run-prep-master.sh" >> "$RV/prep-master.nohup" 2>&1 < /dev/null &
log "prep master relaunched at POOL_W=64 THREADS=4 (orphan partials swept)"

# 3. THE CUT CAMPAIGN: C 2 -> 8 (8 x ~50GB ~ 400GB; memavail post-query ~900GB;
#    if memavail drops under 150G, the successor narrows C - see the handoff)
group_kill run-cuts.sh
sleep 5
setsid env BATCH_W=8 N_CONCURRENT=8 RAYON=24 bash "$RV/run-cuts.sh" >> "$RV/cuts-campaign.nohup" 2>&1 < /dev/null &
log "cut campaign relaunched at BATCH_W=8 N_CONCURRENT=8 RAYON=24"

# 4. THE WAVE: LANES 32 -> 128 (measured per-call peak RSS ~2.0GB)
group_kill run-wave.sh
sleep 20
# a lane killed mid-gpartition leaves partitions/<T> that the retry would
# REFUSE (the discarded-run rule) - clean the unmarked ones so the retry
# is clean (gpartition.done is the marker)
for d in "$RV"/partitions/T*; do
  [ -d "$d" ] || continue
  T=$(basename "$d")
  [ -f "$RV/prep/$T/gpartition.done" ] || rm -rf "$d"
done
log "cleaned unmarked partitions dirs (mid-gpartition kills)"
setsid env LANES=128 bash "$RV/run-wave.sh" >> "$RV/wave.nohup" 2>&1 < /dev/null &
log "wave relaunched at LANES=128"

touch "$MARK"
log "WIDEN COMPLETE: prep POOL_W=64, cuts C=8 W=8 RAYON=24, wave LANES=128"
