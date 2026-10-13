#!/usr/bin/env bash
# THE SORTED-VARIANT IDENTITY GATES: the v2 --region-batching runs at the
# two decisive measurement points (s32 the best rate, s192 the committed
# thread count) vs the rebuilt binary's default-path 192-thread control
# (ctl192) - same gate.py runs mode, same script pattern as run-g23.sh.
# The optimized path must be record-identical to the default path; the
# default path is already record-identical to the committed binary's
# output (G2). Every run carries exactly equal record counts 9,557,434.
B=/home/erikg/hprcv2/proj-bandwidth
G="$B/gates"; mkdir -p "$G"
gate_pair () {  # RUNA RUNB TAG
  RUNA=$1; RUNB=$2; TAG=$3
  args=""
  for i in $(seq 0 18); do
    args="$args $RUNA/shard.fifo_${i}_fastq.bin.gz $RUNB/shard.fifo_${i}_fastq.bin.gz"
    args="$args $RUNA/shard.fifo_${i}_fastq.occ.bin $RUNB/shard.fifo_${i}_fastq.occ.bin"
  done
  if nice -n 19 python3 "$B/gate.py" runs $args > "$G/$TAG.runs.txt" 2>&1; then
    echo "$TAG RUNS-GATE PASS" >> "$G/$TAG.runs.txt"; echo "$TAG PASS"
  else
    echo "$TAG RUNS-GATE FAIL" >> "$G/$TAG.runs.txt"; echo "$TAG FAIL"
  fi
}
gate_pair "$B/optimized2/s32"  "$B/optimized/ctl192" g4-sorted-s32
gate_pair "$B/optimized2/s192" "$B/optimized/ctl192" g5-sorted-s192
echo DONE > "$G/svars.done"
