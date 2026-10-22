#!/usr/bin/env bash
# THE END-TO-END RECEIPT, STEP 1 - PANEL-PREP'S PANEL-SIDE PIECE:
# THE WHOLE-PANEL ROUTES INVENTORY (derive-routes over the ONE
# committed global panel, the PARALLEL LANE-HASH form - the committed
# default; the serial form's ~7.5h wall is receipted in SS10.6.2 and
# the parallel form is the product's own). The reaped routes artifact
# re-derives here, once per panel (~25min). Receipt: wall, core-s,
# peak RSS, and the byte-comparison against the committed identity
# receipt IF AVAILABLE (the committed serial-vs-parallel receipt was
# reaped with the gates-era workspace - stated honestly when absent;
# the artifact's own manifest fields are checked against the committed
# receipted values: 38,790 lanes, 1.403 Tbp total source bp).
source /home/erikg/impg-genome-inference/scripts/end-to-end/ee-env.sh
CENSUS=$EE/cut-census-L1/cut-census.jsonl
OUT=$EE/global-routes
RCP=$EE/receipts/routes.txt
if [ -f "$EE/ledger/routes.done" ]; then echo done; exit 0; fi
: > "$RCP"
[ -f "$CENSUS" ] || { echo "$(ee_now) no census at $CENSUS (run run-ee-census-cut.sh first)" >> "$RCP"; exit 1; }
ee_budget_check || { echo budget >> "$RCP"; exit 1; }
mkdir -p "$OUT"

echo "$(ee_now) derive-routes over the ONE global panel (the parallel lane-hash form, the whole box's rayon)" >> "$RCP"
export RAYON_NUM_THREADS=${RAYON_NUM_THREADS:-192}
t0=$(date +%s)
cpu=$EE/scratch/routes.cpu
( TIMEFORMAT='%U %S'; time (
  "$IMPG_NEW" genome-infer derive-routes \
    --panel "$PANEL" --sources "$AGC" --read-lengths 148 \
    --census "$CENSUS" --out-dir "$OUT"
  ) ) > "$EE/scratch/routes.log" 2> "$cpu" &
dp=$!
gp=$(ee_guard "$dp" "$ONCE_GUARD_GIB" "$RCP" derive-routes)
tp=$(ee_rss_track "$dp" "$RCP" derive-routes)
wait "$dp"; st=$?
kill "$gp" "$tp" 2>/dev/null; wait "$gp" "$tp" 2>/dev/null
wall=$(( $(date +%s) - t0 ))
echo "$(ee_now) derive-routes exit=$st wall=${wall}s cpu='$(cat $cpu 2>/dev/null)' (committed parallel baseline: 1505s)" >> "$RCP"
[ "$st" -ne 0 ] && { echo "$st" > "$EE/ledger/routes.exit"; exit "$st"; }

# ---- THE RECEIPT: the manifest's own counters vs the committed
# receipted values (the committed serial-vs-parallel BYTE receipt is
# reaped with the gates-era workspace - the honest statement) ------
python3 - "$OUT/graph.json" "$RCP" <<'PY'
import json, sys
m = json.load(open(sys.argv[1]))
lanes = len(m["lanes"]) if isinstance(m.get("lanes"), list) else m.get("lanes")
bp = m.get("total_source_bp")
fam = len(m["families"]) if isinstance(m.get("families"), list) else m.get("families")
print(f"routes graph: lanes={lanes} families={fam} total_source_bp={bp}")
open(sys.argv[2], "a").write(
    f"routes graph: lanes={lanes} families={fam} total_source_bp={bp}\n")
if lanes is not None:
    ok = (int(lanes) == 38790)
    open(sys.argv[2], "a").write(
        ("routes lane count MATCHES the committed receipt (38,790 lanes)\n" if ok
         else f"routes lane count DIFFERS from the committed 38,790: {lanes}\n"))
if bp is not None:
    open(sys.argv[2], "a").write(
        f"routes total_source_bp {float(bp)/1e12:.3f} Tbp (the committed receipt: 1.403 Tbp)\n")
PY
echo "$(ee_now) the committed serial-vs-parallel byte-identity receipt (serial-vs-par.routes-cmp-receipt.txt) was REAPED with the gates-era workspace - NOT AVAILABLE for byte comparison; the manifest-field comparison above is the available identity receipt" >> "$RCP"
# the machinery's own manifest receipt (elapsed, the legacy-v1 opt-in provenance)
python3 - "$OUT/manifest.json" >> "$RCP" 2>&1 <<'PY'
import json, sys
m = json.load(open(sys.argv[1]))
s = m.get("stats", {})
print(f"routes manifest: elapsed={m.get('elapsed_seconds')}s "
      f"panel_metadata_version={s.get('panel_metadata_version')} "
      f"panel_legacy_opt_in={s.get('panel_legacy_opt_in')} "
      f"census_records={s.get('census_records')} "
      f"census_occurrences={s.get('census_occurrences')}")
PY
cp "$EE/scratch/routes.log" "$EE/receipts/routes.runlog.txt" 2>/dev/null || true
touch "$EE/ledger/routes.done"
echo "routes done wall=${wall}s"
