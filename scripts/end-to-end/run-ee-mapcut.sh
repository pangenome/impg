#!/usr/bin/env bash
# THE END-TO-END RECEIPT, STEP 2 (per territory) - THE MAP LEG: the
# committed reads fetch rule (the 35x region fetch from the local BAM,
# +-50kb, -F 0x904) + the three-step product MAP step (the local syng
# built into the map's own partition dir + the committed LOCAL-frame
# cut against the ONE global panel - the territory's evidence basis
# cut-reads.fasta for the call step's restricted census profile).
# Marker-gated; wall + core-s receipts; the 64GiB worker guard over
# the whole invocation tree (the cut child included).
# usage: run-ee-mapcut.sh <TID>
source /home/erikg/impg-genome-inference/scripts/end-to-end/ee-env.sh
T=$1
P=$EE/prep/$T
L=$EE/ledger/$T
mkdir -p "$L"
[ -f "$L/map.done" ] && exit 0
[ -f "$L/panelprep.done" ] || { echo "$T panel-prep not done" >&2; exit 1; }
ee_budget_check || exit 1

read -r CHROM START END < <(awk -v t="$T" '$4==t{print $1, $2, $3; exit}' \
  "$RV/bed4-genome.tsv")
[ -n "$CHROM" ] || { echo "$T not in bed4-genome.tsv" >&2; exit 1; }

# 1. the reads fetch (the committed fetch rule, verbatim)
if [ ! -f "$L/reads.done" ]; then
  ee_step_out "$L" reads "$EE/reads/$T.fq.gz" \
    bash -c 'samtools view -u -F 0x904 "$1" "$2" | samtools fastq - 2>"$3" | gzip -1' \
    _ "$BAM" "$CHROM:$((START-50000))-$((END+50000))" "$L/reads.samtools.err" || exit 1
  touch "$L/reads.done"
fi

# 2. the three-step MAP step (one territory per invocation: the map
# root is the call step's own per-territory interface; the chunk roots
# for the fleet call are assembled from these outputs by the chunk
# runner - the batch-composition invariance receipt is their warrant)
# (core-s receipt: TIMEFORMAT written DIRECTLY to map.cpu - no $()
# pipe-hold; the product's raw stdout/stderr stream is the reaped
# class, kept only on failure for the fail-closed diagnosis)
export RAYON_NUM_THREADS=${MAPCUT_T:-32}
OUT=$EE/map/$T
( TIMEFORMAT='%U %S'; time ( "$IMPG_NEW" genome-infer map \
    --panel "$PANEL" --agc "$AGC" --shards-dir "$SHARDS" \
    --partitions "$EE/partitions" \
    --reads "$EE/reads/$T.fq.gz" \
    --region "$T" --out-dir "$OUT" > "$L/map.log" 2>&1 ) ) 2> "$L/map.cpu" &
mp=$!
gpid=$(ee_guard_tree "$mp" "$GUARD_GIB" "$L/map.guard.txt" mapcut)
tpid=$(ee_rss_track "$mp" "$L/map.guard.txt" mapcut-parent)
t0=$(date +%s)
wait "$mp"; st=$?
kill "$gpid" "$tpid" 2>/dev/null; wait "$gpid" "$tpid" 2>/dev/null
echo "$(( $(date +%s) - t0 ))" > "$L/map.wall"
if [ "$st" -ne 0 ]; then
  echo "$st" > "$L/map.exit"
  echo "$(ee_now) $T map FAILED exit=$st (guard: $(cat $L/map.guard.txt 2>/dev/null | tail -1); log kept: $L/map.log)" >> "$EE/logs/failures.log"
  exit "$st"
fi
touch "$L/map.done"
rm -f "$L/map.log"
# the map cut was locus.fa's last in-campaign consumer (the call's
# restricted cut is global-frame: no locus.fa) - reaped here with the
# territory's own reads fetch (the chunk call fetches its own
# multi-region stream), per the stream law, ahead of the chunk's JIT
# reaper
rm -f "$P/locus.fa" "$P/locus.fa.fai" "$EE/reads/$T.fq.gz" 2>/dev/null
echo "$(ee_now) $T map done wall=$(cat $L/map.wall)s"
