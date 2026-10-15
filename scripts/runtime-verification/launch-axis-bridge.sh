#!/usr/bin/env bash
# launcher: the axis bridge (the drift-correct axis, the committed bridge
# recipe's second leg) - detached (setsid) so it survives the orchestrating
# shell. Writes axis-bridges.exit / axis-bridges.wall at completion.
cd /home/erikg/hprcv2/wgs35x-fleet/runtime-verification
setsid bash -c '
  s=$(date +%s)
  python3 rv-axis-bridge.py ../grid/windows.fa hg002-axis-target.fa axis-bridges.paf ${NPROC:-24} \
    > axis-bridge.log 2>&1
  echo $? > axis-bridges.exit
  echo $(( $(date +%s) - s )) > axis-bridges.wall
' >> axis-bridges.mm.err 2>&1 < /dev/null &
echo "axis bridge relaunched (mappy is_primary fix)"
