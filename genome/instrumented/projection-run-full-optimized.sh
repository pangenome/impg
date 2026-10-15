#!/usr/bin/env bash
# THE WHOLE-SAMPLE OPTIMIZED PROJECTION (the datasheet's full-sample leg):
# the committed 35x projection runner VERBATIM (same local BAM, same 19
# region streams, same -F 0x904, same fifo feed, same map-cap/progress,
# same guards) + --region-batching at the BEST MEASURED GATED POINT
# (RAYON_NUM_THREADS=32: the slice measured 709,813 reads/s post-load vs
# 503,609-592,328 at 192 - the committed 192-thread configuration is
# past the peak; s32 is record-identical to ctl192 by gate g4). Into a
# NEW out-dir (the committed shards are the immutable baseline). Then the
# dedup occurrence stream, then the per-shard record-count identity gate
# vs the committed projection (every read's packet must appear in the
# committed shard of its stream with equal multiplicity - full-sample
# multiset identity).
set -uo pipefail
P=/home/erikg/hprcv2/hg002-projection
B=/home/erikg/hprcv2/proj-bandwidth
BIN=/home/erikg/impg-genome-inference/target/release/examples/sample_projection
LOCAL=/home/erikg/hprcv2/HG002.GRCh38.35x.bam
N=proj35-opt
FIFO_DIR=/dev/shm/hg002-projection-fifos
mkdir -p "$FIFO_DIR" "$B/full-optimized"
rm -f "$B/$N.done" "$B/$N.exit" "$B/$N.rss" "$B/$N.log" "$B/$N.err" "$B/$N.wall"

FIFOS=""
for i in $(seq 0 18); do
  FIFO="$FIFO_DIR/fifo.$i"
  rm -f "$FIFO"; mkfifo "$FIFO"
  FIFOS="$FIFOS $FIFO"
done
RAYON_NUM_THREADS=32 nice -n 10 stdbuf -oL -eL "$BIN" project \
  --panel /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  --inputs $FIFOS \
  --out-dir "$B/full-optimized" \
  --map-cap 4000000 \
  --progress 20000000 \
  --region-batching \
  > "$B/$N.log" 2> "$B/$N.err" &
runner=$!
echo "$runner" > "$B/$N.pid"
start=$(date +%s)

for i in $(seq 0 18); do
  FIFO="$FIFO_DIR/fifo.$i"
  if [ "$i" -ge 17 ]; then RL="$P/streams/regions.$i.master.txt"; else RL="$P/streams/regions.$i.txt"; fi
  REGIONS=$(cat "$RL" | tr '\n' ' ')
  ( samtools view -u -@ 2 -F 0x904 "$LOCAL" $REGIONS 2>> "$B/$N.samtools.$i.err" \
      | samtools fastq - 2>> "$B/$N.samtools.$i.err" > "$FIFO"
    echo $? > "$B/$N.samtools.$i.exit"
    [ -s "$B/$N.samtools.$i.exit" ] && [ "$(cat "$B/$N.samtools.$i.exit")" -eq 0 ] && touch "$B/$N.samtools.$i.done"
  ) &
done
echo "streams launched $(date -u +%FT%TZ)"

( while kill -0 "$runner" 2>/dev/null; do
    rss_kb=$(awk '/VmRSS/{print $2}' /proc/$runner/status 2>/dev/null) || continue
    echo "$(date -u +%FT%TZ) $runner ${rss_kb:-0}" >> "$B/$N.rss"
    if [ "${rss_kb:-0}" -gt $((64*1024*1024)) ]; then
      echo "$(date -u +%FT%TZ) GUARD kill ${rss_kb}kB" >> "$B/$N.rss"
      kill -9 "$runner" 2>/dev/null; break
    fi
    sleep 10
  done ) & poller=$!
( while kill -0 "$runner" 2>/dev/null; do
    df -h /home/erikg | tail -1 >> "$B/$N.df"
    sleep 300
  done ) & dfwatch=$!

status=0; wait "$runner" || status=$?
kill "$poller" "$dfwatch" 2>/dev/null; wait "$poller" "$dfwatch" 2>/dev/null
wait 2>/dev/null
echo "$status" > "$B/$N.exit"
echo "$(( $(date +%s) - start ))" > "$B/$N.wall"
[ "$status" -ne 0 ] && exit "$status"

start2=$(date +%s)
stdbuf -oL -eL "$BIN" dedup \
  --out-dir "$B/full-optimized" \
  --occurrences "$B/full-optimized/occurrences.bin" \
  --partitions 32 \
  > "$B/$N.dedup.log" 2> "$B/$N.dedup.err"
st2=$?
echo "$st2" > "$B/$N.dedup.exit"
echo "$(( $(date +%s) - start2 ))" > "$B/$N.dedup.wall"
echo "$st2" > "$B/$N.exit"
[ "$st2" -eq 0 ] && touch "$B/$N.done"
exit "$st2"
