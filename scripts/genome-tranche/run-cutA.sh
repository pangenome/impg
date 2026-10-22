#!/usr/bin/env bash
# THE GENOME TRANCHE - stage 4: the PROFILE-A (restricted) global-frame
# cuts, the stacked-rulings production form (Ruling 1 out-of-production
# = the default binary path; restrict-default census; exact scan). The
# single `cut --frame global` is the ONLY committed restricted path
# (the batch forms' scan is the full routed mass - the named follow-up
# carries the profile into the batch manifests), so each territory pays
# its own panel load; the runs go 4-wide at the committed 64-thread
# single-cut shape, IMPG_CUT_READSET_FASTA naming the territory's
# regenerated committed kept-read set (stage 3's rvcut cut-reads.fasta).
# usage: run-cutA.sh   (env: WT, POOL_W=4, RAYON=64)
set -uo pipefail
WT=${WT:-/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1}
BIN=/home/erikg/impg-genome-inference/target/release/examples/sample_projection
HPRC=/home/erikg/hprcv2
POOL_W=${POOL_W:-4}
RAYON=${RAYON:-64}
mkdir -p "$WT/cutA" "$WT/logs"
export WT BIN HPRC RAYON
grep -v '^#' "$WT/tranche.tsv" | cut -f1 | \
  xargs -P "$POOL_W" -I{} bash -c '
    T={}
    [ -f "$WT/prep/$T/cutA.done" ] && exit 0
    [ -f "$WT/prep/$T/rvcut.done" ] || exit 0
    OUT=$WT/cutA/$T
    rm -rf "$OUT"
    AXIS_JSON=$WT/prep/$T/reference-axis.json
    BEDS=$WT/prep/$T/partition-w10k-d1k/results
    AXIS_PATH=$(cat "$WT/prep/$T/axis-path.txt")
    BASIS=$WT/rvcut/$T/$T/cut-reads.fasta
    start=$(date +%s)
    IMPG_CUT_READSET_FASTA="$BASIS" \
    RAYON_NUM_THREADS=$RAYON nice -n 10 stdbuf -oL -eL "$BIN" cut --frame global \
      --record-origin 1 \
      --global-panel $HPRC/HPRC_r2_assemblies_0.6.1.syng \
      --axis "$AXIS_JSON" --bed-directory "$BEDS" --component "$AXIS_PATH" \
      --agc $HPRC/HPRC_r2_assemblies_0.6.1.agc \
      --shards-dir $HPRC/hg002-projection/full \
      --out-dir "$OUT" \
      > "$WT/logs/cutA-$T.log" 2> "$WT/logs/cutA-$T.err" &
    pid=$!
    st=0; wait "$pid" || st=$?
    echo "$st" > "$WT/logs/cutA-$T.exit"
    echo "$(( $(date +%s) - start ))" > "$WT/logs/cutA-$T.wall"
    if [ "$st" -eq 0 ] && [ -f "$OUT/cut-stats.json" ]; then
      touch "$WT/prep/$T/cutA.done"
      echo "$(date -u +%FT%TZ) cutA $T wall $(cat $WT/logs/cutA-$T.wall)s $(python3 -c "import json;s=json.load(open(\"$OUT/cut-stats.json\"));print({k:s.get(k) for k in (\"census_profile\",\"reads_kept\",\"census_mult\")})")" >> "$WT/logs/stage.log"
    else
      echo "$(date -u +%FT%TZ) cutA $T FAILED exit $st" >> "$WT/logs/stage.log"
      echo "$T cutA failed exit $st: $(tail -2 $WT/logs/cutA-$T.err | head -c 300)" >> "$WT/logs/cutA-failures.txt"
    fi
  '
echo "$(date -u +%FT%TZ) cutA stage complete: $(ls $WT/prep/*/cutA.done 2>/dev/null | wc -l)"
