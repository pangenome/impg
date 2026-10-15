#!/usr/bin/env bash
# GATE 3 of the alignment-induced locality domain: THE FLEET PGGB
# PROJECTION - the committed stage-3 tool (examples/pggb_projection,
# the spell-equality proof + the BED-completeness check) over every
# partition of the 1358-partition union todo that does not already
# carry a projection sidecar. The sidecars (per-row absolute-bp walks
# + per-segment windows + the scorer-interface maps) land at
# pggb-partition-graphs-fleet/projection/, beside the committed
# chrMT/chrI gate receipts at pggb-partition-graphs-chrMT-chrI/
# projection/ (the 56 gate partitions are NOT re-projected - the
# committed receipts of record stand). Marker-idempotent per
# partition (the sidecar is the marker); the fleet-level done marker
# lands only when every todo partition carries a sidecar.
# Usage: run-pggb-fleet-projection.sh
set -uo pipefail
cd /home/erikg/impg-genome-inference || exit 1
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
PROJ="$D/pggb-partition-graphs-fleet/projection"
PGGB="$D/pggb-partition-graphs-fleet"
TODO="$D/pggb-fleet-build-union-todo.txt"
P=pggb-fleet-projection
mkdir -p "$PROJ"
rm -f "$D/$P.done" "$D/$P.exit"
start=$(date +%s)
REMAIN=()
while read -r p || [ -n "$p" ]; do
  [ -n "$p" ] || continue
  [ -f "$PROJ/partition$p.pggb.projection.jsonl" ] || REMAIN+=("$p")
done < "$TODO"
echo "=== fleet pggb projection: ${#REMAIN[@]} of $(grep -c . "$TODO") partitions remaining $(date -u +%FT%TZ)"
if [ "${#REMAIN[@]}" -gt 0 ]; then
  IMPG_CORES=0-31 taskset -c 0-31 nice -n 10 \
    stdbuf -oL -eL target/release/examples/pggb_projection \
    --panel /home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng \
    --agc /home/erikg/yeast/yeast235.agc \
    --export-graphs "$D/partition-graphs" \
    --pggb-dir "$PGGB" \
    --out-dir "$PROJ" \
    --partitions "$(IFS=,; echo "${REMAIN[*]}")" \
    > "$D/$P.log" 2> "$D/$P.err" &
  pid=$!
  ( while kill -0 "$pid" 2>/dev/null; do
      rss_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
      echo "$(date -u +%FT%TZ) $pid ${rss_kb:-0}" >> "$D/$P.rss"
      sleep 10
    done ) & poller=$!
  status=0; wait "$pid" || status=$?
  kill "$poller" 2>/dev/null; wait "$poller" 2>/dev/null
  echo "$status" > "$D/$P.exit"
  if [ "$status" -ne 0 ]; then
    echo "=== fleet pggb projection FAILED exit $status $(date -u +%FT%TZ)"
    exit "$status"
  fi
fi
# the closure check: every todo partition must now carry a sidecar
MISSING=0
while read -r p || [ -n "$p" ]; do
  [ -n "$p" ] || continue
  [ -f "$PROJ/partition$p.pggb.projection.jsonl" ] || { echo "MISSING sidecar partition $p"; MISSING=1; }
done < "$TODO"
if [ "$MISSING" -ne 0 ]; then
  echo "=== fleet pggb projection INCOMPLETE $(date -u +%FT%TZ)"
  exit 1
fi
spell_total=$(grep -c 'spell_diff 0,' "$D/$P.log" || true)
echo "=== fleet pggb projection DONE: $spell_total partitions spell-equal, 0 diffs expected $(date -u +%FT%TZ) wall $(( $(date +%s) - start ))s"
touch "$D/$P.done"
