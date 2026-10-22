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
[ -f "$L/panelprep.poison" ] && exit 0
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
# (THE 2-STRIKE POISON RULE at panel-prep shape: a panelprep that
# fails on a RE-RUN is a DETERMINISTIC pp-terminal - marked
# panelprep.poison, named in the fail-closed census, excluded from
# its chunk by the chunk-readiness gate, instead of retrying forever)
PP_PRIOR=0; [ -f "$L/panelprep.exit" ] && PP_PRIOR=1
# THE CITIZENSHIP GOVERNOR (the window-2 steering): the panel load
# is a ~34GB-per-lane appetite - hold while the floor is breached
ee_mem_governor "panelprep $T"
echo "$(ee_now) $T panelprep START weather=$(ee_weather)" >> "$L/panelprep.weather.txt"
# THE CASCADE-LEFTOVER ADJUDICATION (window-2 law): a leftover
# partition dir carrying a COMPLETE registry (partition.json +
# non-empty beds + the partitions.bed row) is the product's own
# completed output from a killed-mid-run marker gap (the EEXIST
# class: the re-run refuses on the name-dup check) - marked done,
# never re-run; a PARTIAL dir (no registry) is a corpse - removed
if [ -d "$EE/partitions/$T" ]; then
  if [ -f "$EE/partitions/$T/partition.json" ] && [ -n "$(ls "$EE/partitions/$T/beds" 2>/dev/null)" ] \
     && grep -q "	$T$" "$EE/partitions/partitions.bed"; then
    touch "$L/panelprep.done"
    rm -f "$P"/locus-pos64.syng* "$P/export-walks.jsonl" "$P/walk-requests.tsv" 2>/dev/null
    rm -rf "$P/partition-w10k-d1k" "$P/export-graphs" 2>/dev/null
    echo "$(ee_now) $T panelprep ADJUDICATED-COMPLETE (the cascade-leftover class; the artifact verified) weather=$(ee_weather)" >> "$L/panelprep.weather.txt"
    exit 0
  fi
  rm -rf "$EE/partitions/$T"
fi
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
  --out-dir "$EE/partitions"; pst=$?
if [ "$pst" -ne 0 ]; then
  w=$(ee_weather)
  if [ "$pst" -eq 137 ] || [ "$pst" -eq 9 ]; then
    # SIGKILL with no guard over this step: the kernel's
    # memory-pressure killer - KERNEL-WEATHER, never a code failure
    echo "$(ee_now) $T panelprep KERNEL-WEATHER kill (exit $pst; weather $w) - the box's memory-pressure killer; leg re-queues under the governor, NOT a code failure" >> "$EE/logs/failures.log"
    echo "$(ee_now) $T panelprep KERNEL-WEATHER exit=$pst weather=$w" >> "$L/panelprep.weather.txt"
  elif [ "$PP_PRIOR" -eq 1 ]; then
    echo 1 > "$L/panelprep.poison"
    echo "$(ee_now) $T PANELPREP-TERMINAL POISON: panelprep failed again (strike 2, exit $pst; weather $w) - named per the doctrine" >> "$EE/logs/failures.log"
  else
    echo "$(ee_now) $T panelprep FAILED strike 1 (exit $pst; weather $w)" >> "$EE/logs/failures.log"
  fi
  exit 1
fi
touch "$L/panelprep.done"
echo "$(ee_now) $T panelprep END weather=$(ee_weather)" >> "$L/panelprep.weather.txt"

# THE EARLY SUBSTRATE REAP (the budget law's own valve, a window-2
# runner fix): once the partition step has consumed the substrate,
# the dead-weight heavy pieces are the reaped class. THE TWO KEPT
# EXCEPTIONS, each named by its consumer: locus.fa (the MAP cut's
# own --sources input - reaped at the map leg's completion) and
# locality-graphs (the per-locus precompute asset, the call's serve
# - reaped by the chunk's own JIT reaper). The reaped pieces are
# re-derivable from the campaign's surviving smalls (the prep
# chain's own re-derivation receipt) if a product-side re-run ever
# demands it.
rm -f "$P"/locus-pos64.syng* "$P/export-walks.jsonl" "$P/walk-requests.tsv" 2>/dev/null
rm -rf "$P/partition-w10k-d1k" "$P/export-graphs" 2>/dev/null
echo "$(ee_now) $T panel-prep done wall=$(cat $L/panelprep.wall)s"
