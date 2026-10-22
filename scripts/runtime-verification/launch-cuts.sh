#!/usr/bin/env bash
# launcher: the cut campaign at the guard-safe width (env override), detached
# (setsid) so it survives the orchestrating shell.
# THE WIDTH vs THE 64GiB GUARD: W=32 measured 65-69GB peak RSS and was
# guard-killed 4/4; W=16 measured UNSAFE at average density (both first-wave W=16
# batches guard-killed at 69-72.5GB); the ADAPTIVE WIDTH: guard-kill
# halves the width per killed batch (8/4/2/1, the guard is the only oracle). The
# concurrency C is the widening axis after queries.done (RAM-bound), never W.
# During the query stage: BATCH_W=8 N_CONCURRENT=2 RAYON=48 (~2x50GB RSS).
# After queries.done: N_CONCURRENT=8-12 (RAM-bound; the box has ~900GB free).
cd /home/erikg/hprcv2/wgs35x-fleet/runtime-verification
setsid env BATCH_W="${BATCH_W:-8}" N_CONCURRENT="${N_CONCURRENT:-2}" \
  RAYON="${RAYON:-48}" bash ./run-cuts.sh >> cuts-campaign.nohup 2>&1 < /dev/null &
