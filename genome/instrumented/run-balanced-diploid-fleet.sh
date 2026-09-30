#!/usr/bin/env bash
# Resume-safe serial fleet after BOTH chrMT/chrI balanced diploid-prior smoke gates.
set -uo pipefail
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
ROOT=/home/erikg/impg-genome-inference
for C in chrMT chrI; do
  P="$D/run-balanced-$C"
  if [ ! -e "$P.done" ] || [ "$(cat "$P.exit" 2>/dev/null)" != 0 ]; then
    echo "smoke $C not clean; fleet blocked" >&2
    exit 1
  fi
done
for C in chrVI chrIII chrIX chrVIII chrV chrXI chrX chrII chrXIV chrXIII chrXVI chrXII chrXV chrVII chrIV; do
  P="$D/run-balanced-$C"
  if [ -e "$P.done" ] && [ "$(cat "$P.exit" 2>/dev/null)" = 0 ]; then
    echo "skip $C (clean)" >> "$D/fleet-lane.log"
  else
    echo "start $C $(date -u +%FT%TZ)" >> "$D/fleet-lane.log"
    if [ "$C" = chrIV ]; then
      IMPG_CORES=0-63 "$ROOT/genome/instrumented/run-balanced-diploid-component.sh" "$C" balanced
    else
      IMPG_CORES=0-255 "$ROOT/genome/instrumented/run-balanced-diploid-component.sh" "$C" balanced
    fi
    status=$?
    echo "end $C exit=$status wall=$(cat "$P.wall" 2>/dev/null) $(date -u +%FT%TZ)" >> "$D/fleet-lane.log"
    if [ "$status" != 0 ]; then
      echo "$C" > "$D/fleet-lane.failed"
      exit "$status"
    fi
  fi
  A="$D/genotype-distance-balanced-$C"
  rm -f "$A.done" "$A.exit"
  start=$(date +%s)
  IMPG_DIPLOID_VALIDATION_DIR="$D" python3 -B "$ROOT/genome/instrumented/score-diploid-genotype.py" "$C" balanced \
    > "$A.json.tmp" 2> "$A.err"
  status=$?
  echo "$status" > "$A.exit"
  echo "$(( $(date +%s) - start ))" > "$A.wall"
  if [ "$status" != 0 ]; then
    rm -f "$A.json.tmp"
    touch "$A.done"
    echo "$C" > "$D/fleet-lane.failed"
    exit "$status"
  fi
  mv "$A.json.tmp" "$A.json"
  touch "$A.done"
done
IMPG_DIPLOID_VALIDATION_DIR="$D" IMPG_DIPLOID_GENOTYPE_TAG=balanced \
  bash "$ROOT/genome/instrumented/run-diploid-genotype.sh" \
  > "$D/genotype-distance-fleet-summary.out" 2> "$D/genotype-distance-fleet-summary.err"
status=$?
if [ "$status" != 0 ]; then
  echo "genotype-summary" > "$D/fleet-lane.failed"
  exit "$status"
fi
touch "$D/fleet-lane.done"
