#!/usr/bin/env bash
# THE LOCALITY-DOMAIN FLEET, one component (gate 3 of the
# alignment-induced locality domain): the locality scoring chain
# over the fleet pggb substrate - the locality maps
# (locality-domain-maps.py over the gate+fleet pggb builds and
# projection sidecars) -> the serial identity base -> the 4-wide
# exhaustive locality run of record (run-realign-locality.sh, the
# committed gate-2 runner verbatim) -> the serial/4-wide identity
# gate (0 semantic diffs + all four sidecars md5-identical) -> the
# realign checker --territory --locality (sliced per the committed
# chrIV long-run discipline), phase 8N's BEFORE pointed at the
# committed territory receipts of record via --committed-receipt.
# Marker-idempotent per step; the fleet done marker lands last.
# Usage: run-fleet-locality.sh COMP
set -uo pipefail
cd /home/erikg/impg-genome-inference || exit 1
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
C="$1"
PY=python3

read -r NLOCI <<< "$($PY - "$C" << 'EOF'
import json, sys
c = sys.argv[1]
axis = json.load(open(
    "/home/erikg/yeast/genome-mem-bwt-bootstrap-20260911T055203Z/reference-axis-v1.json"
))["intervals"]
print(sum(1 for iv in axis if iv["component"] == f"S288C#0#{c}"))
EOF
)"
LOCI="$($PY -c "print(','.join(str(i) for i in range($NLOCI)))")"
echo "=== $C: $NLOCI loci $(date -u +%FT%TZ)"

if [ ! -f "$D/pggb-fleet-projection.done" ]; then
  echo "$C: the fleet pggb projection is not closed" >&2
  exit 1
fi

# step 1: the locality maps (the window-keyed fold maps + the
# structure receipt)
if [ ! -f "$D/locality-maps-$C.done" ]; then
  $PY genome/instrumented/locality-domain-maps.py "$C" \
    > "$D/locality-maps-$C.log" 2>&1 \
    || { echo "$C locality maps FAILED"; exit 1; }
  touch "$D/locality-maps-$C.done"
fi
echo "=== $C locality maps done $(date -u +%FT%TZ)"

# step 2: the serial identity base (the .exit guard: the runner
# touches .done even on failure - a failed run must re-run)
if [ ! -f "$D/run-realignlocality-serial-$C.done" ] \
   || [ "$(cat "$D/run-realignlocality-serial-$C.exit" 2>/dev/null)" != "0" ]; then
  bash genome/instrumented/run-realign-locality.sh serial "$C" "$LOCI" \
    "realign-locality-serial-$C" \
    > "$D/fleet-locality-$C-serial.log" 2>&1 \
    || { echo "$C serial FAILED"; exit 1; }
fi
echo "=== $C serial done $(date -u +%FT%TZ)"

# step 3: the 4-wide exhaustive locality run of record
if [ ! -f "$D/run-realignlocality-4-$C.done" ] \
   || [ "$(cat "$D/run-realignlocality-4-$C.exit" 2>/dev/null)" != "0" ]; then
  bash genome/instrumented/run-realign-locality.sh 4 "$C" "$LOCI" \
    "realign-locality-exhaustive-$C" \
    > "$D/fleet-locality-$C-4wide.log" 2>&1 \
    || { echo "$C 4wide FAILED"; exit 1; }
fi
echo "=== $C 4wide done $(date -u +%FT%TZ)"

# step 4: the identity gate (serial vs 4-wide: 0 semantic diffs,
# all four sidecars md5-identical)
if [ ! -f "$D/locality-identity-$C.done" ]; then
  $PY - "$C" << 'EOF' > "$D/locality-identity-$C.log" 2>&1 \
    || { echo "$C identity gate FAILED"; exit 1; }
import hashlib, json, sys
C = sys.argv[1]
D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
def strip(d):
    return {k: v for k, v in d.items() if k not in ("rss_kb", "walls")}
n = 0
with open(f"{D}/realign-locality-serial-{C}.jsonl") as fs, \
     open(f"{D}/realign-locality-exhaustive-{C}.jsonl") as fe:
    for ls, le in zip(fs, fe):
        rs, re_ = strip(json.loads(ls)), strip(json.loads(le))
        assert rs == re_, f"locus {rs['locus']} semantic diff"
        n += 1
print(f"identity: {n} loci, 0 semantic diffs")
# (the committed gate pattern: the FOUR SIDECARS are md5-identical;
# the main receipt carries per-locus rss_kb/walls timing fields and
# is proven by the semantic comparison above)
for ext in (".exactness.jsonl", ".skeleton.jsonl",
            ".jsonl.ingredients.jsonl", ".jsonl.records.jsonl"):
    a = hashlib.md5(open(f"{D}/realign-locality-serial-{C}{ext}", "rb")
                    .read()).hexdigest()
    b = hashlib.md5(open(f"{D}/realign-locality-exhaustive-{C}{ext}", "rb")
                    .read()).hexdigest()
    assert a == b, f"sidecar md5 differs: {ext}"
    print(f"md5 {ext}: {a} identical")
EOF
  touch "$D/locality-identity-$C.done"
fi
echo "=== $C identity gate done $(date -u +%FT%TZ)"

# step 5: the checker slices (the committed slice discipline)
if [ ! -f "$D/check-locality-$C.done" ]; then
  rm -f "$D/check-locality-$C-s0.log" "$D/check-locality-$C-s1.log" \
        "$D/check-locality-$C-s2.log"
  if [ "$NLOCI" -le 60 ]; then
    SLICES=("0-$((NLOCI-1))"); LOGS=(s0)
  else
    T1=$(( NLOCI / 3 )); T2=$(( 2 * NLOCI / 3 ))
    SLICES=("0-$((T1-1))" "$T1-$((T2-1))" "$T2-$((NLOCI-1))")
    LOGS=(s0 s1 s2)
  fi
  i=0
  for SL in "${SLICES[@]}"; do
    $PY genome/instrumented/check-realign-scoring.py --exhaustive "$C" \
        --territory --locality \
        --timered-base "realign-locality-exhaustive-$C" \
        --run-tag "realignlocality-4-$C" \
        --committed-receipt "$D/realign-exhaustive-territory-$C.jsonl" \
        --loci-slice "$SL" \
        > "$D/check-locality-$C-${LOGS[$i]}.log" 2>&1
    st=$?
    if [ $st -ne 0 ] || ! grep -q 'ALL PHASES PASS' "$D/check-locality-$C-${LOGS[$i]}.log"; then
      echo "$C locality checker slice $SL FAILED"; exit 1
    fi
    i=$((i+1))
  done
  touch "$D/check-locality-$C.done"
fi
echo "=== $C locality checker done $(date -u +%FT%TZ)"
touch "$D/fleet-locality-$C.done"
echo "=== $C LOCALITY CHAIN DONE $(date -u +%FT%TZ)"
