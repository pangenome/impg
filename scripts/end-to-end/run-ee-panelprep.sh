#!/usr/bin/env bash
# THE END-TO-END RECEIPT, STEP 1 (per territory) - PANEL-PREP, the
# three-step product's own STEP 1 invocation at the CURRENT STACK:
# the whole-panel routes inventory (REUSED loudly, built once per
# panel) + the PARTITION global form (the window beds regenerated
# from the locality maps, the registry recording the global pair, the
# content origin d=1 - the uniform substrate convention, the delta=1
# receipts at both substrate classes) + the panel-side per-locus
# precompute asset (slices 1+2, the current stack's own panel-prep
# third piece; its once-per-panel price is receipted per territory and
# the asset is REAPED with the territory's substrate after its call -
# the stream law's named exception: consumer = the call's serve, reap
# = the leg's completion, budget = the workdir bound).
# usage: run-ee-panelprep.sh <TID>
source /home/erikg/impg-genome-inference/scripts/end-to-end/ee-env.sh
T=$1
P=$EE/prep/$T
L=$EE/ledger/$T
mkdir -p "$L"
[ -f "$L/panelprep.done" ] && exit 0
[ -f "$P/substrate.done" ] || { echo "$T substrate not done" >&2; exit 1; }
[ -f "$EE/ledger/routes.done" ] || { echo "routes inventory not done" >&2; exit 1; }
ee_budget_check || exit 1
CENSUS=$EE/cut-census-L1/cut-census.jsonl

# the BED4 row (the committed bed4-genome.tsv is the universe ledger)
read -r CHROM START END < <(awk -v t="$T" '$4==t{print $1, $2, $3; exit}' \
  "$RV/bed4-genome.tsv")
[ -n "$CHROM" ] || { echo "$T not in bed4-genome.tsv" >&2; exit 1; }

# the panel-prep invocation pays the global panel load for the
# precompute asset build (the ~200s class, receipted per territory)
export RAYON_NUM_THREADS=${PANELPREP_T:-8}
ee_step_out "$L" panelprep "$L/panelprep.log" \
  env RAYON_NUM_THREADS=$RAYON_NUM_THREADS "$IMPG_NEW" genome-infer panel-prep \
  --name "$T" --chrom "$CHROM" --start "$START" --end "$END" \
  --sources "$P/locus.fa" \
  --component "$(cat "$P/axis-path.txt")" \
  --axis "$P/reference-axis.json" \
  --partition-graphs "$P/locality-graphs" \
  --census "$CENSUS" --read-lengths 148 \
  --global-panel "$PANEL" --global-agc "$AGC" \
  --global-routes "$EE/global-routes" \
  --record-origin 1 \
  --out-dir "$EE/partitions" || exit 1
touch "$L/panelprep.done"
echo "$(ee_now) $T panel-prep done wall=$(cat $L/panelprep.wall)s"
