#!/usr/bin/env bash
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
# G2: the rebuilt binary's default path vs the committed-binary sweep run
gate_pair "$B/sweep/t192" "$B/optimized/ctl192" g2-default-path
# G3-strict: the v1 strict region-batched run vs the default-path control
gate_pair "$B/optimized/ctl192" "$B/optimized/o192" g3-strict
