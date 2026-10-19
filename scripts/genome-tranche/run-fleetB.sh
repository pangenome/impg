#!/usr/bin/env bash
# THE GENOME TRANCHE - stage 5: the FLEET PIPELINE at genome scale (its
# first use beyond the four gate localities): fork-cut-call-batch-global,
# the parent loading the global panel ONCE (the shared-panel-load
# receipt), width-32 child batches (--rows-per-child 32), the children's
# in-child calls at the receipted T16 shape with the committed
# 4-concurrent call rung (the first 4 rows of each child carry the call
# object; the other rows' calls are the truth-assessment stage's own
# runs). The batch form's census is the GLOBAL-NATIVE full mass (the
# machinery's own limit, honestly receipted: the read-set instrument
# does not apply on the batch path) - these cuts are the tranche's
# profile-B (global-native) receipts.
# usage: run-fleetB.sh   (env: WT, RPC=32, RAYON=16, TREE_GUARD_GB=700)
set -uo pipefail
WT=${WT:-/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1}
BIN=/home/erikg/impg-genome-inference/target/release/examples/sample_projection
HPRC=/home/erikg/hprcv2
ROUTES=/home/erikg/hprcv2/global-panel-call-mode/global-routes-par
RPC=${RPC:-32}
RAYON=${RAYON:-16}
GUARD=${TREE_GUARD_GB:-700}
BATCH=$WT/fleetB
mkdir -p "$BATCH" "$WT/cutB" "$WT/prod" "$WT/logs"

# the manifest: every tranche row cut (out_dir cutB/<T>); the FIRST 4
# rows of every RPC-sized chunk carry the in-child call (the committed
# 4-concurrent rung). The call rows use the committed 35x region reads.
: > "$BATCH/manifest.jsonl"
i=0
while IFS=$'\t' read -r T STRAT CHROM S E REST; do
  [ -f "$WT/prep/$T/cutB.done" ] && { i=$(( i + 1 )); continue; }
  grep -q "^$T\b" "$WT/tranche-failed.tsv" 2>/dev/null && { i=$(( i + 1 )); continue; }
  LOCI=$(python3 -c "
import json
a=json.load(open('$WT/prep/$T/reference-axis.json'))
n=sum(1 for iv in a['intervals'] if iv['component'] == '$(cat "$WT/prep/$T/axis-path.txt")')
print(','.join(str(j) for j in range(n)))")
  CHUNK=$(( i / RPC ))
  if [ $(( i % RPC )) -lt 4 ]; then
    python3 - "$T" "$LOCI" >> "$BATCH/manifest.jsonl" <<'EOF'
import json, sys
T, loci = sys.argv[1], sys.argv[2]
WT = "/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1"
print(json.dumps({"name": T, "axis": f"{WT}/prep/{T}/reference-axis.json",
  "bed_directory": f"{WT}/prep/{T}/partition-w10k-d1k/results",
  "component": open(f"{WT}/prep/{T}/axis-path.txt").read().strip(),
  "out_dir": f"{WT}/cutB/{T}", "record_origin": 1,
  "call": {"reads": f"{WT}/reads/{T}.fq.gz",
           "partition_graphs": f"{WT}/prep/{T}/locality-graphs",
           "loci": loci, "out_dir": f"{WT}/prod/call-{T}"}}))
EOF
  else
    python3 - "$T" >> "$BATCH/manifest.jsonl" <<'EOF'
import json, sys
T = sys.argv[1]
WT = "/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1"
print(json.dumps({"name": T, "axis": f"{WT}/prep/{T}/reference-axis.json",
  "bed_directory": f"{WT}/prep/{T}/partition-w10k-d1k/results",
  "component": open(f"{WT}/prep/{T}/axis-path.txt").read().strip(),
  "out_dir": f"{WT}/cutB/{T}", "record_origin": 1}))
EOF
  fi
  i=$(( i + 1 ))
done < <(grep -v '^#' "$WT/tranche.tsv")
N=$(wc -l < "$BATCH/manifest.jsonl")
[ "$N" -gt 0 ] || { echo "nothing to do"; exit 0; }
echo "$(date -u +%FT%TZ) fleetB: $N rows, rows-per-child $RPC, T$RAYON" | tee "$WT/logs/fleetB.log"

start=$(date +%s)
IMPG_PANEL_LEGACY_V1_OK=1 \
  IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS=1 \
  IMPG_REALIGN_PARALLEL_SEAMS=1 \
  IMPG_REALIGN_SEAM_WIDTH=16 \
  RAYON_NUM_THREADS=$RAYON nice -n 10 stdbuf -oL -eL "$BIN" fork-cut-call-batch-global \
  --global-panel $HPRC/HPRC_r2_assemblies_0.6.1.syng \
  --agc $HPRC/HPRC_r2_assemblies_0.6.1.agc \
  --shards-dir $HPRC/hg002-projection/full \
  --manifest "$BATCH/manifest.jsonl" \
  --routes "$ROUTES" \
  --rows-per-child "$RPC" \
  > "$BATCH/fleet.log" 2> "$BATCH/fleet.err" &
pid=$!
( while kill -0 "$pid" 2>/dev/null; do
    total_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || total_kb=0
    for c in $(pgrep -P "$pid" 2>/dev/null); do
      ck=$(awk '/VmRSS/{print $2}' /proc/$c/status 2>/dev/null) || ck=0
      total_kb=$((total_kb + ck))
    done
    echo "$(date -u +%FT%TZ) tree-rss=$(( total_kb / 1024 ))MB children=$(pgrep -P "$pid" 2>/dev/null | wc -l)" >> "$BATCH/fleet.rss"
    if [ "${total_kb:-0}" -gt $((GUARD*1024*1024)) ]; then
      echo "$(date -u +%FT%TZ) TREE GUARD kill ${total_kb}kB" >> "$BATCH/fleet.rss"
      kill -9 "$pid" 2>/dev/null; break
    fi
    sleep 10
  done ) & poll=$!
st=0; wait "$pid" || st=$?
kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
echo "$st" > "$BATCH/fleet.exit"
echo "$(( $(date +%s) - start ))" > "$BATCH/fleet.wall"
echo "$(date -u +%FT%TZ) fleetB exit $st wall $(cat "$BATCH/fleet.wall")s peak tree RSS $(sed -n 's/.*tree-rss=\([0-9]*\)MB.*/\1/p' "$BATCH/fleet.rss" | sort -n | tail -1)MB" | tee -a "$WT/logs/fleetB.log"
# mark the cut-done per territory
while read -r row; do
  T=$(python3 -c "import json,sys;print(json.loads(sys.argv[1])['name'])" "$row")
  [ -f "$WT/cutB/$T/cut-stats.json" ] && touch "$WT/prep/$T/cutB.done"
done < "$BATCH/manifest.jsonl"
grep -E "\[fork-cut-call-batch-global\] (panel loaded|child .* cut phase done)" "$BATCH/fleet.err" | cut -c1-160 >> "$WT/logs/fleetB.log"
exit "$st"
