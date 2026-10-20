#!/usr/bin/env bash
# THE GENOME TRANCHE - stage 6: THE OFFICIAL TRUTH ASSESSMENT (the same
# assessment as every banked accuracy receipt): call-diplotypes per
# territory per profile (A = the restricted cut, B = the global-native
# cut) with --truth-qv <axis> <truth2>, IMPG_TRUTH_QV_ALLWAVE=1, the
# likegt max-qv convention (cap 50), RAYON_NUM_THREADS=1 (the banked
# serial recipe; the truth-QV flag adds the test artifact, never
# changes the product output). The truth pair: hap1 = the bridge-fixed
# axis row (prep/<T>/axis-path.txt, the committed bridge output), hap2 =
# the truth2.tsv row (the HG002#2 leg of the same bridge machinery;
# fail-closed classes named). Lanes run LANES-wide over the tranche.
# usage: run-truthcalls.sh <A|B|C>   (env: WT, LANES=6)
# PROFILE C = THE NEW PRODUCT DEFAULT (the owner's far-bin ruling, commit
# 2a3dd46): the RESTRICTED census (form A's own cut, reused - the prune is
# a CALL-side census-view flag, no new cut) + --far-bin-prune - the pure
# window-geometry rule (touched-window diameter <= 1; boundary-straddlers
# keep both adjacent bins). Assessed with the same official QV machinery.
set -uo pipefail
WT=${WT:-/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1}
PROFILE=$1     # A (restricted) | B (global-native)
IMPG=/home/erikg/impg-genome-inference/target/release/impg
QVBIWFA=/home/erikg/impg-genome-inference/genome/instrumented/qv-biwfa/target/release/qv-biwfa
HPRC=/home/erikg/hprcv2
ROUTES=/home/erikg/hprcv2/global-panel-call-mode/global-routes-par
LANES=${LANES:-6}
CUTDIR=$WT/cutA; EXTRAS="--far-bin-prune"; CUTMARK=cutA.done
# (THE DEFAULTS-FLIP FIX, commit 431d55a: the 22:10 impg defaults the
# prune ON in call-diplotypes, so every form states its own census view
# EXPLICITLY: A/B = the pre-ruling committed UNPRUNED view
# (--far-bin-prune=false), C = the new product default, ON.)
if [ "$PROFILE" = "A" ]; then EXTRAS="--far-bin-prune=false"; fi
if [ "$PROFILE" = "B" ]; then EXTRAS="--far-bin-prune=false"; CUTDIR=$WT/cutB; CUTMARK=cutB.done; fi
mkdir -p "$WT/call$PROFILE" "$WT/logs"
export WT IMPG QVBIWFA HPRC ROUTES PROFILE CUTDIR CUTMARK EXTRAS LANES

# the work list: territories with both cuts done and an ok truth pair
grep -v '^#' "$WT/tranche.tsv" | cut -f1 | while read -r T; do
  [ -f "$WT/prep/$T/cut${PROFILE}.done" ] || continue
  [ -f "$WT/prep/$T/call${PROFILE}.done" ] && continue
  cls=$(awk -F'\t' -v t="$T" '$1==t {print $4}' "$WT/truth2.tsv")
  if [ "$cls" != "ok" ]; then
    [ -f "$WT/prep/$T/call${PROFILE}.skipped" ] || \
      echo "$(date -u +%FT%TZ) $T call$PROFILE SKIPPED truth2 class $cls" >> "$WT/logs/stage.log"
    touch "$WT/prep/$T/call${PROFILE}.skipped"
  fi
done | true

run_one() {
  local T=$1
  [ -f "$WT/prep/$T/$CUTMARK" ] || return 0
  [ -f "$WT/prep/$T/call${PROFILE}.done" ] && return 0
  [ -f "$WT/prep/$T/call${PROFILE}.skipped" ] && return 0
  local TRUTH2 AXIS OUT
  TRUTH2=$(awk -F'\t' -v t="$T" '$1==t {print $3}' "$WT/truth2.tsv")
  AXIS=$(cat "$WT/prep/$T/axis-path.txt")
  [ -n "$TRUTH2" ] || return 0
  OUT=$WT/call$PROFILE/$T
  rm -rf "$OUT"
  local LOCI
  LOCI=$(python3 -c "
import json
a=json.load(open('$WT/prep/$T/reference-axis.json'))
n=sum(1 for iv in a['intervals'] if iv['component'] == '$AXIS')
print(','.join(str(j) for j in range(n)))")
  local start
  start=$(date +%s)
  IMPG_PANEL_LEGACY_V1_OK=1 \
  IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS=1 \
  IMPG_TRUTH_QV_ALLWAVE=1 \
  RAYON_NUM_THREADS=1 nice -n 10 stdbuf -oL -eL "$IMPG" genome-infer call-diplotypes \
    --panel $HPRC/HPRC_r2_assemblies_0.6.1.syng \
    --routes "$ROUTES" \
    --partition-graphs "$WT/prep/$T/locality-graphs" \
    --component "$AXIS" \
    --census "$CUTDIR/$T/cut-census.jsonl" \
    --reads "$WT/reads/$T.fq.gz" \
    --derive-cache "$CUTDIR/$T/cut-cache.bin" \
    --loci "$LOCI" \
    --truth-qv "$AXIS" "$TRUTH2" $EXTRAS \
    --qv-biwfa "$QVBIWFA" \
    --global-frame --record-origin 1 \
    --out-dir "$OUT" \
    > "$WT/logs/call$PROFILE-$T.log" 2> "$WT/logs/call$PROFILE-$T.err" &
  local pid=$!
  ( while kill -0 "$pid" 2>/dev/null; do
      rss_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
      echo "$(date -u +%FT%TZ) $rss_kb" >> "$WT/logs/call$PROFILE-$T.rss"
      sleep 15
    done ) &
  local poll=$!
  local st=0; wait "$pid" || st=$?
  kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$WT/logs/call$PROFILE-$T.exit"
  echo "$(( $(date +%s) - start ))" > "$WT/logs/call$PROFILE-$T.wall"
  if [ "$st" -eq 0 ] && [ -s "$OUT/calls.jsonl" ]; then
    touch "$WT/prep/$T/call${PROFILE}.done"
    echo "$(date -u +%FT%TZ) call$PROFILE $T OK wall $(cat $WT/logs/call$PROFILE-$T.wall)s" >> "$WT/logs/stage.log"
  else
    echo "$st" > "$WT/prep/$T/call${PROFILE}.exit.code"
    echo "$(date -u +%FT%TZ) call$PROFILE $T FAILED exit $st: $(tail -1 $WT/logs/call$PROFILE-$T.err | head -c 300)" >> "$WT/logs/stage.log"
  fi
}
export -f run_one

grep -v '^#' "$WT/tranche.tsv" | cut -f1 | \
  xargs -P "$LANES" -I{} bash -c 'run_one {}'
echo "$(date -u +%FT%TZ) call$PROFILE stage complete: $(ls $WT/prep/*/call${PROFILE}.done 2>/dev/null | wc -l) ok"
