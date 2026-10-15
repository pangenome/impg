#!/usr/bin/env bash
# launcher: the cut campaign at the guard-safe width (env override), detached
# (setsid) so it survives the orchestrating shell.
# THE WIDTH vs THE 64GiB GUARD: W=32 measured 65-69GB peak RSS and was
# guard-killed 4/4; W=16 is the guard-safe width (measured batch-0 W=1 peak
# 38.99GB; batch memory grows with the batch's aggregate routing-index
# nodes). The handoff's W=96 widen plan is IMPOSSIBLE under the guard - the
# concurrency C is the widening axis after queries.done (RAM-bound), never W.
# During the query stage: BATCH_W=16 N_CONCURRENT=2 RAYON=48 (~2x50GB RSS).
# After queries.done: N_CONCURRENT=4 (or 6 if RAM allows).
cd /home/erikg/hprcv2/wgs35x-fleet/runtime-verification
setsid env BATCH_W="${BATCH_W:-16}" N_CONCURRENT="${N_CONCURRENT:-2}" \
  RAYON="${RAYON:-48}" bash ./run-cuts.sh >> cuts-campaign.nohup 2>&1 < /dev/null &
