#!/usr/bin/env bash
# GATE 3 of the alignment-induced locality domain: THE FLEET PGBB
# BUILDS - the pggb graphs for every remaining component's build set
# (the union of the 15 components' committed partition-graph build
# lists, minus the 56 already built for the chrMT/chrI gate), into the
# shared fleet dir. Marker-idempotent per partition (the committed
# driver skips partition<P>/graph.gfa.done); the union todo list is
# the committed receipt pggb-fleet-build-union-todo.txt.
# Usage: run-pggb-fleet-builds.sh [OUTDIR-NAME]
set -uo pipefail
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
OUTNAME="${1:-pggb-partition-graphs-fleet}"
cd /home/erikg/impg-genome-inference || exit 1
mapfile -t TODO < "$D/pggb-fleet-build-union-todo.txt"
if [ "${#TODO[@]}" -eq 0 ]; then
  echo "nothing to build" >&2
  exit 0
fi
echo "=== fleet pggb build: ${#TODO[@]} partitions -> $D/$OUTNAME $(date -u +%FT%TZ)"
bash genome/instrumented/run-pggb-partition-graphs.sh "$OUTNAME" "${TODO[@]}"
status=$?
echo "=== fleet pggb build exit $status $(date -u +%FT%TZ)"
exit "$status"
