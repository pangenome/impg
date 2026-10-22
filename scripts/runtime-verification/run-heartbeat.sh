#!/usr/bin/env bash
# THE HEARTBEAT: every 60s, the disk (the high-water mark of record) + RAM +
# the stage counters, into wave-logs/heartbeat.log. THE DISK HIGH-WATER is
# measured from these lines (the streaming reaper's own receipt).
RV=/home/erikg/hprcv2/wgs35x-fleet/runtime-verification
while true; do
  AV=$(df -BG --output=avail / | tail -1 | tr -dc '0-9')
  MA=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)
  QD=$(ls $RV/query/shard-*.done 2>/dev/null | wc -l)
  SUB=$(ls $RV/prep/*/substrate.done 2>/dev/null | wc -l)
  CUT=$(ls $RV/prep/*/cut.done 2>/dev/null | wc -l)
  POIS=$(ls $RV/prep/*/cut.poison 2>/dev/null | wc -l)
  CALL=$(ls $RV/prep/*/call.done 2>/dev/null | wc -l)
  FAIL=$(ls $RV/prep/*/call.exit 2>/dev/null | wc -l)
  SKIP=$(ls $RV/prep/*/call.skipped-no-evidence 2>/dev/null | wc -l)
  RVSZ=$(du -sBG $RV 2>/dev/null | awk '{print $1}')
  echo "$(date -u +%FT%TZ) df ${AV}G memavail ${MA}GiB query $QD/14 substrate $SUB cut $CUT poison $POIS call $CALL fail $FAIL noev $SKIP rvsize $RVSZ" >> $RV/wave-logs/heartbeat.log
  sleep 60
done
