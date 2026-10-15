#!/usr/bin/env bash
cd /home/erikg/hprcv2/wgs35x-fleet/runtime-verification
setsid env POOL_W=32 bash ./run-prep-master.sh >> prep-master.nohup 2>&1 < /dev/null &
