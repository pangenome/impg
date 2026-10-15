#!/usr/bin/env bash
# THE WHOLE-GENOME FLEET, the heavy chain for ONE component (the chrIV
# recipe exactly): wait for the component's multicensus + light
# phases, then anchor projection -> serial identity base -> 4-wide
# exhaustive run of record -> the anchor checker + the realign
# checker (sliced per the chrIV long-run discipline) -> the tables.
# Every phase is the committed generic runner with its own markers;
# this driver adds only the fleet-level done-marker.
# Usage: run-fleet-heavy.sh COMP
set -uo pipefail
cd /home/erikg/impg-genome-inference
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
C="$1"
PY=python3

# ---- wait for the deps (the multicensus run + the light phases)
while [ ! -f "$D/run-multicensus-pilot-$C.done" ] || [ ! -f "$D/fleet-light-$C.done" ]; do
  sleep 20
done
echo "=== $C deps ready $(date -u +%FT%TZ)"

# ---- the loci list: every window of the component (the axis file's
# S288C#0#COMP intervals; asserted 1:1 with the balanced receipt's rows)
read -r NLOCI BALN <<< "$($PY - "$C" << 'EOF'
import json, sys
c = sys.argv[1]
D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
axis = json.load(open(
    "/home/erikg/yeast/genome-mem-bwt-bootstrap-20260911T055203Z/reference-axis-v1.json"
))["intervals"]
n = sum(1 for iv in axis if iv["component"] == f"S288C#0#{c}")
bal = json.load(open(f"{D}/genotype-distance-balanced-{c}.json"))["rows"]
assert n == len(bal), (n, len(bal))
print(n, len(bal))
EOF
)"
LOCI="$($PY -c "print(','.join(str(i) for i in range($NLOCI)))")"
echo "=== $C: $NLOCI loci"

# ---- the anchor-projection receipts (the committed receipt names ARE
# the fleet's receipts of record; the committed marker convention is
# the bare anchorproj tag; the two earliest chains ran under the
# per-component tag anchorproj-<C> before the convention was adopted -
# both accepted on resume)
if [ ! -f "$D/run-anchorproj-$C.done" ] && [ ! -f "$D/run-anchorproj-$C-$C.done" ]; then
  bash genome/instrumented/run-anchor-projection.sh "$C" "anchorproj" \
    > "$D/fleet-heavy-$C-anchor.log" 2>&1 || { echo "$C anchor FAILED"; exit 1; }
fi
echo "=== $C anchor done $(date -u +%FT%TZ)"

# ---- the serial identity base
if [ ! -f "$D/run-realignpar-serial-$C-$C.done" ]; then
  bash genome/instrumented/run-realign-run.sh serial "$C" "$LOCI" \
    "realign-par-serial-$C" "realignpar-serial-$C-$C" \
    > "$D/fleet-heavy-$C-serial.log" 2>&1 || { echo "$C serial FAILED"; exit 1; }
fi
echo "=== $C serial done $(date -u +%FT%TZ)"

# ---- the 4-wide exhaustive run of record
if [ ! -f "$D/run-realignexhaustive-$C-$C.done" ]; then
  bash genome/instrumented/run-realign-run.sh 4 "$C" "$LOCI" \
    "realign-exhaustive-$C" "realignexhaustive-$C-$C" \
    > "$D/fleet-heavy-$C-4wide.log" 2>&1 || { echo "$C 4wide FAILED"; exit 1; }
fi
echo "=== $C 4-wide done $(date -u +%FT%TZ)"

# ---- the anchor checker
if [ ! -f "$D/check-anchor-projection-$C.done" ]; then
  rm -f "$D/check-anchor-projection-$C.log"
  if $PY genome/instrumented/check-anchor-projection.py "$C" \
      > "$D/check-anchor-projection-$C.log" 2>&1; then
    touch "$D/check-anchor-projection-$C.done"
  else
    echo "$C anchor checker FAILED"; exit 1
  fi
fi
echo "=== $C anchor checker done $(date -u +%FT%TZ)"

# ---- the realign checker (the corrected checker; sliced for the
# long components per the chrIV discipline)
if [ ! -f "$D/check-realign-$C.done" ]; then
  rm -f "$D/check-realign-$C-s0.log" "$D/check-realign-$C-s1.log" "$D/check-realign-$C-s2.log"
  if [ "$NLOCI" -le 60 ]; then
    SLICES=("0-$((NLOCI-1))")
    LOGS=(s0)
  else
    T1=$(( NLOCI / 3 )); T2=$(( 2 * NLOCI / 3 ))
    SLICES=("0-$((T1-1))" "$T1-$((T2-1))" "$T2-$((NLOCI-1))")
    LOGS=(s0 s1 s2)
  fi
  i=0
  for SL in "${SLICES[@]}"; do
    if ! $PY genome/instrumented/check-realign-scoring.py --exhaustive "$C" \
        --timered-base "realign-exhaustive-$C" --run-tag "realignexhaustive-$C-$C" \
        --loci-slice "$SL" > "$D/check-realign-$C-${LOGS[$i]}.log" 2>&1; then
      echo "$C realign checker slice $SL FAILED"; exit 1
    fi
    i=$((i+1))
  done
  touch "$D/check-realign-$C.done"
fi
echo "=== $C realign checker done $(date -u +%FT%TZ)"

# ---- the partition-graphs checker (the survey/build/census receipts)
if [ ! -f "$D/check-partition-graphs-$C.done" ]; then
  rm -f "$D/check-partition-graphs-$C.log"
  if $PY genome/instrumented/check-partition-graphs-chrIV.py "$C" \
      > "$D/check-partition-graphs-$C.log" 2>&1; then
    touch "$D/check-partition-graphs-$C.done"
  else
    echo "$C partition-graphs checker FAILED"; exit 1
  fi
fi
echo "=== $C partition-graphs checker done $(date -u +%FT%TZ)"

# ---- the tables
$PY genome/instrumented/realign-chrIV-tables.py "$C" > "$D/realign-$C-tables.stdout" 2>&1 \
  || { echo "$C tables FAILED"; exit 1; }
echo "=== $C tables done $(date -u +%FT%TZ)"

touch "$D/fleet-heavy-$C.done"
echo "=== $C HEAVY CHAIN COMPLETE $(date -u +%FT%TZ)"
