#!/usr/bin/env bash
# THE H2 A/B LEG (the owner's 12-lane shape vs one shared index): 12
# processes, each holding its OWN private in-RAM copy of the panel syng
# (~46GB each, ~550GB total - the owner's observed shape), each fed a
# DISJOINT 1/12 of the same measured slice, 16 rayon threads per lane
# (12 x 16 = 192, the owner's "shouldn't even be running 192 cores"
# configuration). Measured: per-lane post-load walls, the aggregate
# reads/s, the SUM of per-lane peak RSS - against the single-process
# sweep points on the same slice. RAM-guarded (abort if the box's
# available memory falls below 120GiB during the load wave).
set -uo pipefail
B=/home/erikg/hprcv2/proj-bandwidth
BIN=/home/erikg/impg-genome-inference/target/release/examples/sample_projection
PANEL=/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng
L="$B/lanes12"
rm -rf "$L"; mkdir -p "$L/inputs" "$L/runs"
# 12 disjoint lane inputs over the same slice: record-aligned awk split
# of the concatenated 19 stream windows (4 lines per read; 7,782,400
# reads -> 12 disjoint lanes of ~648,533 reads each, none dropped).
awk -v n=12 -v total=7782400 -v base="$L/inputs" \
  'BEGIN{for(i=0;i<n;i++) f[i]=sprintf("%s/lane_%02d", base, i)}
   { rec=int((NR-1)/4); print >> f[int(rec*n/total)] }' \
  <(cat "$B"/slice/fifo_*.fastq)
wc -l "$L"/inputs/lane_* | head -13
( while :; do
    AVAIL=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)
    echo "$(date -u +%FT%TZ) avail_GiB $AVAIL" >> "$L/guard"
    if [ "$AVAIL" -lt 120 ]; then
      echo "$(date -u +%FT%TZ) GUARD TRIP avail ${AVAIL}GiB - killing lanes" >> "$L/guard"
      pkill -9 -f "lanes12/runs" || true
      exit 1
    fi
    sleep 10
  done ) & guard=$!
for j in $(seq 0 11); do
  JJ=$(printf "%02d" $j)
  ( RAYON_NUM_THREADS=16 nice -n 10 /usr/bin/time -v -o "$L/runs/lane$JJ.time" \
      "$BIN" project --panel "$PANEL" \
      --inputs "$L/inputs/lane_$JJ" \
      --out-dir "$L/runs/lane$JJ" \
      --map-cap 4000000 --progress 400000 \
      > "$L/runs/lane$JJ.log" 2> "$L/runs/lane$JJ.err" ) &
done
wait
kill "$guard" 2>/dev/null; wait "$guard" 2>/dev/null
echo DONE > "$L/done"
