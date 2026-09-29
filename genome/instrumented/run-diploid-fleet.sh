#!/usr/bin/env bash
# Resume-safe serial full-box diploid fleet. chrIV uses the proven 64-core form.
set -uo pipefail
D=/home/erikg/yeast/genome-diploid-validation-20260929
cd /home/erikg/impg-genome-inference
for C in chrVI chrIII chrIX chrVIII chrV chrXI chrX chrII chrXIV chrXIII chrXVI chrXII chrXV chrVII chrIV; do
  P="$D/run-fleet-$C"
  if [ -e "$P.done" ] && [ "$(cat "$P.exit" 2>/dev/null)" = 0 ]; then
    echo "skip $C (clean)" >> "$D/fleet-lane.log"
    continue
  fi
  echo "start $C $(date -u +%FT%TZ)" >> "$D/fleet-lane.log"
  if [ "$C" = chrIV ]; then
    IMPG_CORES=0-63 genome/instrumented/run-diploid-component.sh "$C" fleet
  else
    IMPG_CORES=0-255 genome/instrumented/run-diploid-component.sh "$C" fleet
  fi
  status=$?
  echo "end $C exit=$status wall=$(cat "$P.wall" 2>/dev/null) $(date -u +%FT%TZ)" >> "$D/fleet-lane.log"
  if [ "$status" != 0 ]; then
    echo "$C" > "$D/fleet-lane.failed"
    exit "$status"
  fi
  python3 genome/instrumented/score-diploid.py "$C" fleet > "$D/score-fleet-$C.json" 2> "$D/score-fleet-$C.err"
  status=$?
  echo "$status" > "$D/score-fleet-$C.exit"
  if [ "$status" != 0 ]; then
    echo "$C score" > "$D/fleet-lane.failed"
    exit "$status"
  fi
done
touch "$D/fleet-lane.done"
