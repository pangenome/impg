#!/usr/bin/env bash
# THE GENOME TRANCHE - the BOUNDED truth-call pass (the partial-deliverable
# accelerator): run the official truth assessment for a FIXED member list
# and one profile, then EXIT (never loops, never waits - the master's own
# TA/TB stages own everything outside the list, and the per-territory
# callX.done markers make the master skip exactly what landed here).
# usage: run-truth-bounded.sh <A|B|C> <list-file>   (env: WT, LANES=4)
set -uo pipefail
WT=${WT:-/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1}
PROFILE=$1
LIST=$2
IMPG=/home/erikg/impg-genome-inference/target/release/impg
QVBIWFA=/home/erikg/impg-genome-inference/genome/instrumented/qv-biwfa/target/release/qv-biwfa
HPRC=/home/erikg/hprcv2
ROUTES=/home/erikg/hprcv2/global-panel-call-mode/global-routes-par
LANES=${LANES:-4}
CUTDIR=$WT/cutA; EXTRAS="--far-bin-prune"; CUTMARK=cutA.done
# (THE DEFAULTS-FLIP FIX, commit 431d55a: the 22:10 impg defaults the
# prune ON in call-diplotypes, so every form states its own census view
# EXPLICITLY: A/B = the pre-ruling committed UNPRUNED view
# (--far-bin-prune=false), C = the new product default, ON.)
if [ "$PROFILE" = "A" ]; then EXTRAS="--far-bin-prune=false"; fi
if [ "$PROFILE" = "B" ]; then EXTRAS="--far-bin-prune=false"; CUTDIR=$WT/cutB; CUTMARK=cutB.done; fi
export WT IMPG QVBIWFA HPRC ROUTES PROFILE CUTDIR CUTMARK EXTRAS

run_one() {
  local T=$1
  [ -f "$WT/prep/$T/$CUTMARK" ] || return 0
  [ -f "$WT/prep/$T/call${PROFILE}.done" ] && return 0
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

xargs -P "$LANES" -I{} bash -c 'run_one {}' < "$LIST"
echo "$(date -u +%FT%TZ) bounded $PROFILE pass complete over $(wc -l < "$LIST") listed members"
