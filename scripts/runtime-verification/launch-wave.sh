#!/usr/bin/env bash
# launcher: the wave at LANES lanes (env override), detached (setsid) so it
# survives the orchestrating shell
cd /home/erikg/hprcv2/wgs35x-fleet/runtime-verification
setsid env LANES="${LANES:-32}" bash ./run-wave.sh >> wave.nohup 2>&1 < /dev/null &
