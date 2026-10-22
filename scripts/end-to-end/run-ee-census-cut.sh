#!/usr/bin/env bash
# THE END-TO-END RECEIPT, STEP 1 PIECE 0 - THE CENSUS-PRODUCING CUT.
# The whole-panel routes inventory (derive-routes) demands a
# global-frame census (the routing bind: every occurrence's path id
# resolves in the global panel's name map); the committed global-frame
# cut censuses were reaped with the gates-era workspace, so the receipt
# re-derives ONE here: the committed global-frame cut recipe at L1,
# d=1 (the uniform substrate convention, the delta=1 receipts), the
# GLOBAL-NATIVE census profile (the record form, the same class the
# committed routes artifact was bound to). Marker-gated, guarded,
# receipts at write time.
source /home/erikg/impg-genome-inference/scripts/end-to-end/ee-env.sh
OUT=$EE/cut-census-L1
RCP=$EE/receipts/census-cut-L1.txt
mkdir -p "$OUT"
if [ -f "$EE/ledger/census-cut.done" ]; then
  echo "census cut already done (marker)" ; exit 0
fi
: > "$RCP"
echo "$(ee_now) the L1 global-native census cut (the routing-bind input for the routes inventory)" >> "$RCP"
ee_budget_check || { echo "budget" >> "$RCP"; exit 1; }

# the cut is a whole-box-shaped one-time step: RAYON at 32 (the honest
# co-tenant form), guarded at the ONCE bound, RSS tracked.
export RAYON_NUM_THREADS=${RAYON_NUM_THREADS:-32}
rm -f "$OUT"/cut-census.jsonl "$OUT"/cut-cache.bin "$OUT"/cut-reads.fasta "$OUT"/cut-stats.json
"$SP_NEW" cut \
  --frame global --record-origin 1 --census-profile global-native \
  --global-panel "$PANEL" --agc "$AGC" \
  --axis "$L1_AXIS" --bed-directory "$L1_BEDS" \
  --component "$L1_COMPONENT" \
  --shards-dir "$SHARDS" --out-dir "$OUT" &
cutpid=$!
gp=$(ee_guard "$cutpid" "$ONCE_GUARD_GIB" "$RCP" census-cut)
tp=$(ee_rss_track "$cutpid" "$RCP" census-cut)
wait "$cutpid"; st=$?
kill "$gp" "$tp" 2>/dev/null; wait "$gp" "$tp" 2>/dev/null
wall=$SECONDS
if [ "$st" -ne 0 ]; then
  echo "$(ee_now) FAILED exit=$st (guard events above, if any)" >> "$RCP"
  echo "$st" > "$EE/ledger/census-cut.exit"
  exit "$st"
fi
# the receipt: the cut's own stats + the walls + the core-s
python3 - "$OUT/cut-stats.json" "$RCP" <<'PY'
import json, sys
s = json.load(open(sys.argv[1]))
keys = ("reads_scanned","reads_kept","census_mult","census_records","routed_patterns",
        "territory_rows","territory_entries","load_seconds","scan_seconds")
line = {k: s.get(k) for k in keys if k in s}
print(f"cut-stats: {json.dumps(line)}")
open(sys.argv[2], "a").write(f"cut-stats: {json.dumps(line)}\n")
PY
touch "$EE/ledger/census-cut.done"
echo "$(ee_now) DONE wall=${wall}s (see $RCP for the RSS peak; core-s below)" >> "$RCP"
echo "census cut done wall=${wall}s"
