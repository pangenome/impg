#!/usr/bin/env bash
# THE TERRITORY-NORMALIZED FLEET, one component (stage 2 of the owner's
# approved four): the extent-normalization convention's scoring chain
# over the committed artifacts — the census, anchor receipts, partition
# graphs and expressibility all STAND (the convention changes only the
# comparison domain); the chain is serial identity base -> 4-wide
# exhaustive territory run of record -> the realign checker --territory
# (sliced per the chrIV long-run discipline). The anchor/graph/census
# checkers do not re-run (their artifacts are unchanged).
# Usage: run-fleet-territory.sh COMP
set -uo pipefail
cd /home/erikg/impg-genome-inference
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
C="$1"
PY=python3

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

if [ ! -f "$D/run-territory-serial-$C-$C.done" ]; then
  IMPG_REALIGN_TERRITORY_NORMALIZED=1 \
  bash genome/instrumented/run-realign-run.sh serial "$C" "$LOCI" \
    "realign-territory-serial-$C" "territory-serial-$C-$C" \
    > "$D/fleet-territory-$C-serial.log" 2>&1 || { echo "$C serial FAILED"; exit 1; }
fi
echo "=== $C serial done $(date -u +%FT%TZ)"

if [ ! -f "$D/run-territory-$C-$C.done" ]; then
  IMPG_REALIGN_TERRITORY_NORMALIZED=1 \
  bash genome/instrumented/run-realign-run.sh 4 "$C" "$LOCI" \
    "realign-exhaustive-territory-$C" "territory-$C-$C" \
    > "$D/fleet-territory-$C-4wide.log" 2>&1 || { echo "$C 4wide FAILED"; exit 1; }
fi
echo "=== $C 4-wide done $(date -u +%FT%TZ)"

if [ ! -f "$D/check-territory-$C.done" ]; then
  rm -f "$D/check-territory-$C-s0.log" "$D/check-territory-$C-s1.log" "$D/check-territory-$C-s2.log"
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
        --territory --timered-base "realign-exhaustive-territory-$C" \
        --run-tag "territory-$C-$C" --loci-slice "$SL" \
        > "$D/check-territory-$C-${LOGS[$i]}.log" 2>&1; then
      echo "$C territory checker slice $SL FAILED"; exit 1
    fi
    i=$((i+1))
  done
  touch "$D/check-territory-$C.done"
fi
echo "=== $C territory checker done $(date -u +%FT%TZ)"
touch "$D/fleet-territory-$C.done"
echo "=== $C TERRITORY CHAIN DONE $(date -u +%FT%TZ)"
