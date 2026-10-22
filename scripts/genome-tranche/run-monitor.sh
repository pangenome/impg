#!/usr/bin/env bash
# THE TRANCHE MONITOR - the master relaunch loop, WORK-DIR QUALIFIED (the
# 2026-10-20 fix: the first recovery/t2-monitor pair shared the bare
# pgrep pattern "bash run-master.sh", so T2's monitor was PACIFIED by
# T1's live master and never relaunched T2's (T2's master had
# self-exited at 06:11Z on a stale NOK - it read tranche-failed.tsv
# before T00762 was named - so stage K never marked and F never ran;
# 2h of fleet time lost, receipted). This monitor checks the master by
# CWD, not by name, and is parameterized by WT.
set -uo pipefail
WT=${1:?usage: run-monitor.sh <WT>}
MAX=${2:-8}
cd "$WT"
REL=0
while [ "$REL" -lt "$MAX" ]; do
  [ -f logs/stage-TA.done ] && [ -f logs/stage-TB.done ] && break
  alive=""
  for p in $(pgrep -f 'bash run-master.sh' 2>/dev/null); do
    [ "$(readlink /proc/$p/cwd 2>/dev/null)" = "$WT" ] && alive=1
  done
  if [ -z "$alive" ]; then
    while pgrep -f "run-truth-bounded.sh [ABC] $WT/" >/dev/null; do sleep 120; done
    echo "$(date -u +%FT%TZ) monitor: relaunching master #$REL (substrate=$(ls prep/*/substrate.done 2>/dev/null|wc -l) cutA=$(ls prep/*/cutA.done 2>/dev/null|wc -l) TA=$(ls prep/*/callA.done 2>/dev/null|wc -l) TB=$(ls prep/*/callB.done 2>/dev/null|wc -l))" >> logs/stage.log
    nohup bash run-master.sh >> logs/master.nohup 2>&1 &
    disown
    REL=$((REL+1))
  fi
  sleep 300
done
echo "$(date -u +%FT%TZ) monitor exiting (TA=$(ls prep/*/callA.done 2>/dev/null|wc -l) TB=$(ls prep/*/callB.done 2>/dev/null|wc -l) relaunches=$REL)" >> logs/stage.log
