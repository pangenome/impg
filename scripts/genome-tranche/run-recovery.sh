#!/usr/bin/env bash
# THE GENOME TRANCHE - the T2-session recovery driver (the empty-basis
# class repair + the master relaunch loop). THE CLASS, named: the pre-fix
# rvcut batch-composition bug (stale members re-entering slices) composed
# batches containing members whose rvcut.done ALREADY existed; one_batch_w
# does `rm -rf $CUTS/$T/$T` at composition, and where such a batch was
# then guard-killed, the member's cut-reads.fasta (the restriction basis)
# was DESTROYED while its rvcut.done stayed stale-set - cutA can then
# never run (Os NotFound on the basis). The 8-member class:
# T02596 T04149 T09925 T00041 T00130 T07448 T04006 T02594
# (T09025 is the separate tranche-failed extract-length-mismatch member.)
# THE FIX: clear the stale rvcut.done markers so the as-is run-rvcut.sh
# PENDING loop re-cuts them (the determinism gate re-receipts each),
# then the as-is run-cutA.sh, then a bounded A pass, then the master
# relaunch loop until TA and TB are complete (T07036 callA is also
# cleared here: the one TRUE concurrent double-write of the session -
# master TA and a stale bounded lane ran it simultaneously, 01:11-01:19Z;
# T08776 callA's os-error-2 failure was the stale lane's rm -rf of the
# master's out-dir mid-run - both redo cleanly via the markers).
set -uo pipefail
WT=/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1
cd "$WT"
echo "$(date -u +%FT%TZ) recovery: clearing stale rvcut.done for the empty-basis 8" >> logs/stage.log
for T in T02596 T04149 T09925 T00041 T00130 T07448 T04006 T02594; do
  rm -f prep/$T/rvcut.done prep/$T/cut.width
done
# T07036 callA contamination redo (concurrent double-write, named above)
rm -rf callA/T07036 prep/T07036/callA.done
echo "$(date -u +%FT%TZ) recovery: T07036 callA cleared for clean redo" >> logs/stage.log

# 1: the rvcut re-cut (as-is script; only the 8 are pending now)
bash run-rvcut.sh >> logs/stage-K.log 2>&1
n8=0
for T in T02596 T04149 T09925 T00041 T00130 T07448 T04006 T02594; do
  [ -f rvcut/$T/$T/cut-reads.fasta ] && n8=$((n8+1))
done
echo "$(date -u +%FT%TZ) recovery: rvcut re-cut $n8/8 bases present" >> logs/stage.log
[ "$n8" -eq 8 ] || { echo "$(date -u +%FT%TZ) recovery: BASIS RE-CUT INCOMPLETE ($n8/8) - STOPPING" >> logs/stage.log; exit 1; }

# 2: the cutA pass (as-is script, marker-gated - only the 8 + T07036 redo lack it)
POOL_W=4 RAYON=32 bash run-cutA.sh >> logs/stage-A.log 2>&1
echo "$(date -u +%FT%TZ) recovery: cutA pass done, $(ls prep/*/cutA.done | wc -l) cutA.done" >> logs/stage.log

# 3: bounded A pass over the recovered members (the master's TA already
# passed them by; the master relaunch below will also re-check everyone)
: > lists-recovery.txt
for T in T02596 T04149 T09925 T00041 T00130 T07448 T04006 T02594 T07036 T08776; do
  echo "$T" >> lists-recovery.txt
done
LANES=3 bash run-truth-bounded.sh A lists-recovery.txt >> logs/bounded-A-rec.nohup 2>&1

# 4: the master relaunch loop until TA and TB are complete (capped: the
# deterministic exit-1 fail-closed members must not loop forever)
REL=0
while [ "$REL" -lt 6 ]; do
  [ -f logs/stage-TA.done ] && [ -f logs/stage-TB.done ] && break
  if ! pgrep -f "bash run-master.sh" >/dev/null; then
    # never relaunch the master while a bounded pass could still hold a
    # member mid-call (the same-form double-write class this session
    # already had to repair once - see the header)
    while pgrep -f "run-truth-bounded.sh [ABC] $WT/lists-recovery.txt" >/dev/null; do
      sleep 120
    done
    echo "$(date -u +%FT%TZ) recovery: relaunching master #$REL (TA=$(ls prep/*/callA.done|wc -l) TB=$(ls prep/*/callB.done|wc -l))" >> logs/stage.log
    nohup bash run-master.sh >> logs/master.nohup 2>&1 &
    disown
    REL=$((REL+1))
  fi
  sleep 300
done
echo "$(date -u +%FT%TZ) recovery loop exiting (TA=$(ls prep/*/callA.done|wc -l) TB=$(ls prep/*/callB.done|wc -l) relaunches=$REL)" >> logs/stage.log
