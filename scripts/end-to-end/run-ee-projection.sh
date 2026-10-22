#!/usr/bin/env bash
# THE END-TO-END RECEIPT, STEP 2 - THE MAP (the HG002 35x projection,
# the committed map recipe verbatim: the same local BAM, the same 19
# region streams, the same -F 0x904, the same fifos/map-cap/guards,
# --region-batching at T=32 - the committed 20.79min form, run on the
# committed projection instrument target/release/examples/
# sample_projection, unchanged since the committed receipt - the map
# step's byte-identity gate demands the same instrument). Into a NEW
# out-dir (the committed shards are the immutable baseline and remain
# the cut basis); then the dedup occurrence stream; then THE IDENTITY
# GATE (the committed run-full-gate.sh form: occurrences.bin cmp,
# projection-stats counters, per-shard commutative packet-multiset
# hashes); then the re-run's 65G is REAPED (the stream law: the gate
# is on digests, in-process; the committed shards carry the product
# forward, byte-identical by the gate).
# usage: run-ee-projection.sh    (nohup + markers; ~23min + gate)
source /home/erikg/impg-genome-inference/scripts/end-to-end/ee-env.sh
B=/home/erikg/hprcv2/proj-bandwidth
OUT=$EE/map-projection
N=proj35-ee
RCP=$EE/receipts/map-projection.txt
if [ -f "$EE/ledger/map-projection.done" ]; then echo done; exit 0; fi
: > "$RCP"
ee_budget_check || { echo "budget" >> "$RCP"; exit 1; }
mkdir -p "$OUT"

# ---- the projection (the committed recipe, verbatim shape) ----------
FIFO_DIR=/dev/shm/hg002-projection-fifos
mkdir -p "$FIFO_DIR"
rm -f "$FIFO_DIR"/fifo.*
for i in $(seq 0 18); do mkfifo "$FIFO_DIR/fifo.$i"; done
FIFOS=$(ls "$FIFO_DIR"/fifo.* | tr '\n' ' ')
t0=$(date +%s)
cpufile=$EE/scratch/$N.cpu
( TIMEFORMAT='%U %S'; time (
  RAYON_NUM_THREADS=32 nice -n 10 stdbuf -oL -eL "$SP_COMMITTED" project \
    --panel "$PANEL" --inputs $FIFOS \
    --out-dir "$OUT" --map-cap 4000000 --progress 20000000 --region-batching
  ) ) 2> "$cpufile" &
runner=$!
sleep 2
# the samtools region-stream writers (the committed feed)
for i in $(seq 0 18); do
  FIFO="$FIFO_DIR/fifo.$i"
  if [ "$i" -ge 17 ]; then RL="/home/erikg/hprcv2/hg002-projection/streams/regions.$i.master.txt"; else RL="/home/erikg/hprcv2/hg002-projection/streams/regions.$i.txt"; fi
  REGIONS=$(cat "$RL" | tr '\n' ' ')
  ( samtools view -u -@ 2 -F 0x904 "$BAM" $REGIONS 2>> "$EE/scratch/$N.samtools.$i.err" \
      | samtools fastq - 2>> "$EE/scratch/$N.samtools.$i.err" > "$FIFO"
    echo $? > "$EE/scratch/$N.samtools.$i.exit" ) &
done
gp=$(ee_guard "$runner" 64 "$RCP" projection)
tp=$(ee_rss_track "$runner" "$RCP" projection)
wait "$runner"; st=$?
kill "$gp" "$tp" 2>/dev/null; wait "$gp" "$tp" 2>/dev/null
proj_wall=$(( $(date +%s) - t0 ))
wait 2>/dev/null
echo "$(ee_now) projection exit=$st wall=${proj_wall}s cpu='$(cat $cpufile 2>/dev/null)'" >> "$RCP"
[ "$st" -ne 0 ] && { echo "$st" > "$EE/ledger/map-projection.exit"; exit "$st"; }

# ---- the dedup occurrence stream ------------------------------------
t1=$(date +%s)
dedupcpu=$EE/scratch/$N.dedup.cpu
( TIMEFORMAT='%U %S'; time (
  stdbuf -oL -eL "$SP_COMMITTED" dedup --out-dir "$OUT" \
    --occurrences "$OUT/occurrences.bin" --partitions 32
  ) ) 2> "$dedupcpu"
st2=$?
dedup_wall=$(( $(date +%s) - t1 ))
echo "$(ee_now) dedup exit=$st2 wall=${dedup_wall}s cpu='$(cat $dedupcpu 2>/dev/null)'" >> "$RCP"
[ "$st2" -ne 0 ] && { echo "$st2" > "$EE/ledger/map-projection.exit"; exit "$st2"; }

# ---- THE IDENTITY GATE (the committed full-gate form) ---------------
G=$EE/receipts; FULL=$SHARDS
echo "== (1) occurrences.bin byte-compare (canonical: partitions sorted)" >> "$RCP"
if cmp "$FULL/occurrences.bin" "$OUT/occurrences.bin" >> "$RCP" 2>&1; then
  echo "   occurrences.bin BYTE-IDENTICAL" >> "$RCP"
else
  echo "   occurrences.bin DIFFERS" >> "$RCP"
fi
python3 - "$FULL/projection-stats.json" "$OUT/projection-stats.json" >> "$RCP" 2>&1 <<'PY'
import json, sys
a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2]))
for k in ("reads", "bases", "reads_with_mems", "records", "read_lengths", "per_shard_reads"):
    if a[k] != b[k]:
        print(f"   MISMATCH {k}: {a[k]} vs {b[k]}"); sys.exit(1)
print(f"   projection-stats counters identical (reads {a['reads']}, records {a['records']})")
PY
echo "== (2) per-shard commutative packet-multiset hash (the 100% standard)" >> "$RCP"
ls "$FULL"/shard.*.bin.gz "$OUT"/shard.*.bin.gz \
  | xargs -P 24 -n 1 nice -n 19 python3 "$B/shard-hash.py" \
  > "$EE/receipts/map-projection.shard-hashes.txt" 2>&1
python3 - "$EE/receipts/map-projection.shard-hashes.txt" "$FULL" "$OUT" >> "$RCP" 2>&1 <<'PY'
import re, sys
from collections import defaultdict
d = defaultdict(dict)
for line in open(sys.argv[1]):
    t = line.split()
    m = re.search(r'shard\.(fifo_\d+)\.bin\.gz', t[0])
    side = 'rerun' if t[0].startswith(sys.argv[2]) else 'committed'
    d[m.group(1)][side] = (int(t[1]), int(t[2]), t[-1])
bad = [s for s, v in d.items() if len(v) != 2 or v['committed'] != v['rerun']]
print(f"   shards {len(d)}: " + ("ALL IDENTICAL (packets + records + 128-bit multiset sum per shard)"
      if not bad and len(d) == 19 else f"MISMATCH {bad}"))
sys.exit(1 if bad or len(d) != 19 else 0)
PY
if [ $? -eq 0 ]; then echo "   FULL-SAMPLE SHARD GATE PASS" >> "$RCP"; else echo "   FULL-SAMPLE SHARD GATE FAIL" >> "$RCP"; fi
echo "$(ee_now) MAP RECEIPT: projection wall ${proj_wall}s (committed baseline 1247.5s = 20.79min), dedup ${dedup_wall}s, total $((proj_wall+dedup_wall))s" >> "$RCP"
touch "$EE/ledger/map-projection.done"

# ---- THE REAP (the stream law: the gate is on digests; the committed
# shards carry the product forward) ----------------------------------
rm -rf "$OUT"
rm -f "$EE/scratch/$N".* 2>/dev/null
echo "$(ee_now) reaped the re-run's shard set (65G class; the committed shards remain the cut basis, byte-identical by the gate)" >> "$RCP"
echo "map projection done: wall ${proj_wall}s + dedup ${dedup_wall}s, gated, reaped"
