#!/usr/bin/env bash
# THE WHOLE-GENOME RUNTIME VERIFICATION - STAGE 2 MASTER: the prep pool
# behind the query stage. As query shards complete, the refined-bed split
# re-runs over the completed shards' outputs (idempotent - accepted.txt
# grows monotonically), and the per-territory EXPORT-PROFILE substrate chain
# (run-prep-one.sh) runs at POOL_W wide over accepted territories.
# Marker-gated everywhere; resumable; killable.
# usage: run-prep-master.sh   (env: RV, POOL_W=64, THREADS=4)
set -uo pipefail
RV=${RV:-/home/erikg/hprcv2/wgs35x-fleet/runtime-verification}
Q=$RV/query
POOL_W=${POOL_W:-64}
THREADS=${THREADS:-4}
export RV THREADS
PREP_ONE=$RV/run-prep-one.sh
export -f printf >/dev/null 2>&1 || true

split_done_count() { ls "$Q"/shard-*.done 2>/dev/null | wc -l; }

last_split=-1
while true; do
  n=$(split_done_count)
  all_in=$(ls "$Q"/shard-*.bed 2>/dev/null | wc -l)
  if [ "$n" -gt "$last_split" ]; then
    # the split over every COMPLETED shard's output (idempotent recompute)
    outs=()
    for d in "$Q"/shard-*.done; do
      [ -e "$d" ] || continue
      S=$(basename "$d" .done)
      outs+=("$Q/$S.out.bed")
    done
    if [ "${#outs[@]}" -gt 0 ]; then
      python3 "$RV/rv-split-refined-beds.py" "$RV" "$RV/../grid/bridge-accepted.bed" "${outs[@]}" \
        >> "$RV/prep-split.log" 2>&1
    fi
    last_split=$n
    echo "$(date -u +%FT%TZ) split at $n shards: $(tail -1 "$RV/prep-split.log")" >> "$RV/prep-master.log"
  fi
  # the prep pool over accepted, not-yet-done territories
  mapfile -t TODO < <(while read -r T; do
      [ -f "$RV/prep/$T/substrate.done" ] && continue
      [ -f "$RV/prep/$T/substrate.skipped" ] && continue
      [ -f "$RV/prep/$T/split.rejected" ] && continue
      echo "$T"
    done < "$RV/prep/accepted.txt")
  if [ "${#TODO[@]}" -gt 0 ]; then
    printf '%s\n' "${TODO[@]}" | xargs -P "$POOL_W" -I{} bash "$PREP_ONE" {}
  fi
  # terminal when the query is complete and every accepted territory is done
  if [ -f "$Q/queries.done" ]; then
    nleft=$(while read -r T; do
        [ -f "$RV/prep/$T/substrate.done" ] || [ -f "$RV/prep/$T/substrate.skipped" ] || echo x
      done < "$RV/prep/accepted.txt" | wc -l)
    if [ "$nleft" -eq 0 ]; then break; fi
  fi
  sleep 60
done
echo "=== prep master complete $(date -u +%FT%TZ)" >> "$RV/prep-master.log"
touch "$RV/prep-master.done"
