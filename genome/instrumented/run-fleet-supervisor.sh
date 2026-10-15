#!/usr/bin/env bash
# THE WHOLE-GENOME FLEET supervisor: walks the 14 fleet components in
# size order (small first — the small chains close earliest), launching
# run-fleet-heavy.sh for each component as its deps land, keeping at
# most MAXHEAVY heavy chains live (the multicensus runs already ride
# the box; the heavy chains add anchor/scoring/checker load).
# Usage: run-fleet-supervisor.sh [MAXHEAVY]
set -uo pipefail
cd /home/erikg/impg-genome-inference
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
MAXHEAVY="${1:-4}"
ORDER=(chrVI chrIII chrIX chrV chrVIII chrXI chrX chrXIV chrII chrXIII chrXVI chrXII chrXV chrVII)
declare -A LIVE
while true; do
  n_live=0
  for C in "${!LIVE[@]}"; do
    if ! kill -0 "${LIVE[$C]}" 2>/dev/null; then
      unset "LIVE[$C]"
    else
      n_live=$((n_live+1))
    fi
  done
  all_done=1
  for C in "${ORDER[@]}"; do
    [ -f "$D/fleet-heavy-$C.done" ] && continue
    all_done=0
    [ -n "${LIVE[$C]:-}" ] && continue
    [ "$n_live" -ge "$MAXHEAVY" ] && break
    nohup nice -n 10 bash genome/instrumented/run-fleet-heavy.sh "$C" \
      > "$D/fleet-heavy-$C.log" 2>&1 &
    LIVE[$C]=$!
    echo "$(date -u +%FT%TZ) launched heavy chain $C pid ${LIVE[$C]}"
    n_live=$((n_live+1))
  done
  [ "$all_done" -eq 1 ] && break
  sleep 30
done
echo "$(date -u +%FT%TZ) ALL 14 FLEET HEAVY CHAINS COMPLETE"
