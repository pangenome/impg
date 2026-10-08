#!/usr/bin/env bash
# THE WHOLE-GENOME FLEET, the light phases per component (survey
# requests -> the partition-graph build + dumps -> the expressibility
# census), serial within the lane, per-component done-markers at the
# validation dir so the state survives the window. The heavy phases
# (multicensus, anchor projection, scoring) are launched per component
# by run-fleet-heavy.sh once the light phases land.
# Usage: run-fleet-light.sh COMP...
set -uo pipefail
cd /home/erikg/impg-genome-inference
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
for C in "$@"; do
  if [ -f "$D/fleet-light-$C.done" ]; then echo "$C light already done"; continue; fi
  echo "=== $C: requests $(date -u +%FT%TZ)"
  if ! python3 genome/instrumented/partition-expressibility-census-chrIV.py \
      requests "$C" > "$D/fleet-light-$C-requests.log" 2>&1; then
    echo "$C requests FAILED"; exit 1
  fi
  echo "=== $C: graphs + dumps $(date -u +%FT%TZ)"
  if ! bash genome/instrumented/run-chrIV-partition-graphs.sh "$C" \
      > "$D/fleet-light-$C-build.log" 2>&1; then
    echo "$C build FAILED"; exit 1
  fi
  echo "=== $C: census $(date -u +%FT%TZ)"
  if ! python3 genome/instrumented/partition-expressibility-census-chrIV.py \
      census "$C" > "$D/fleet-light-$C-census.log" 2>&1; then
    echo "$C census FAILED"; exit 1
  fi
  touch "$D/fleet-light-$C.done"
  echo "=== $C light complete $(date -u +%FT%TZ)"
done
echo "lane complete: $*"
