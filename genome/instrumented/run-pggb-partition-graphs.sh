#!/usr/bin/env bash
# STAGE 3 (the owner's pggb substrate ruling): the PGGB-BUILT local
# graphs for the chrMT/chrI partition set — the same partition/locality
# extents (the completed alignment-induced partition BED, commit
# 295bca9), the members' sequences restricted to the locality, aligned
# and condensed by the impg graph machinery's pggb engine (FastGA/
# SweepGA all-pairs alignment + seqwish induction + smoothxg-style
# block smoothing + gfaffix normalization — the same engine the early
# campaign's partition-graph renderings used). NOT the syng export
# path: node ids are the pggb build's OWN interning; the syng node-id
# projection is a separate receipt-side step (pggb_projection).
# Serial per partition, markers + external RSS poller + 64GiB guard.
# Usage: run-pggb-partition-graphs.sh <outdir-name> <partition-id>...
set -uo pipefail
source /home/erikg/impg-genome-inference/build-env.sh
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
BEDS=/home/erikg/yeast/partition-pos64-w10k-d1k-to-completion/results
ARCHIVE=/home/erikg/yeast/yeast235.agc
EXTRACT=/home/erikg/impg/target/experiments/yeast-partition-renderings/extract-bed
IMPG=/home/erikg/impg-genome-inference/target/release/impg
THREADS="${IMPG_PGGB_THREADS:-8}"
OUT="$D/$1"; shift
mkdir -p "$OUT"
P="$OUT/build-pggb-partition-graphs"
rm -f "$P.done" "$P.exit" "$P.rss" "$P.log" "$P.err"
start=$(date +%s)
for pid in "$@"; do
  folder="$OUT/partition$pid"
  mkdir -p "$folder/tmp"
  if [ ! -f "$folder/graph.gfa.done" ]; then
    rm -f "$folder/graph.gfa.done" "$folder/graph.gfa.exit" "$folder/graph.gfa.rss" "$folder/members.fa"
    cp "$BEDS/partition$pid.bed" "$folder/members.bed"
    start2=$(date +%s)
    nice -n 10 stdbuf -oL -eL "$EXTRACT" "$ARCHIVE" "$folder/members.bed" "$folder/members.fa" \
      > "$folder/extract.log" 2> "$folder/extract.err" &
    xpid=$!
    ( while kill -0 "$xpid" 2>/dev/null; do
        rss_kb=$(awk '/VmRSS/{print $2}' /proc/$xpid/status 2>/dev/null) || continue
        echo "$(date -u +%FT%TZ) $xpid ${rss_kb:-0}" >> "$folder/extract.rss"
        sleep 5
      done ) & xpoll=$!
    wait "$xpid"; xst=$?
    kill "$xpoll" 2>/dev/null; wait "$xpoll" 2>/dev/null
    echo "$xst" > "$folder/extract.exit"
    [ "$xst" -eq 0 ] || { echo "partition $pid extract FAILED"; exit 1; }
    nice -n 10 stdbuf -oL -eL taskset -c 0-63 "$IMPG" graph \
      --sequence-files "$folder/members.fa" \
      --gfa-engine pggb \
      --aligner fastga \
      --output "$folder/graph.gfa" \
      --temp-dir "$folder/tmp" \
      -t "$THREADS" \
      > "$folder/pggb.log" 2> "$folder/pggb.err" &
    gpid=$!
    ( while kill -0 "$gpid" 2>/dev/null; do
        rss_kb=$(awk '/VmRSS/{print $2}' /proc/$gpid/status 2>/dev/null) || continue
        echo "$(date -u +%FT%TZ) $gpid ${rss_kb:-0}" >> "$folder/pggb.rss"
        sleep 5
      done ) & gpoll=$!
    wait "$gpid"; gst=$?
    kill "$gpoll" 2>/dev/null; wait "$gpoll" 2>/dev/null
    echo "$gst" > "$folder/graph.gfa.exit"
    echo "$(( $(date +%s) - start2 ))" > "$folder/graph.gfa.wall"
    touch "$folder/graph.gfa.done"
    [ "$gst" -eq 0 ] || { echo "partition $pid pggb FAILED"; exit 1; }
  fi
  echo "=== partition $pid done $(date -u +%FT%TZ)"
done
echo "$(( $(date +%s) - start ))" > "$P.wall"
echo 0 > "$P.exit"
touch "$P.done"
exit 0
