#!/usr/bin/env bash
set -euo pipefail
/usr/bin/time -v /usr/bin/timeout 120s /home/erikg/impg/.wg-worktrees/agent-609/target/debug/impg graph --sequence-files tests/test_data/local_compression/microtangle_repeat_motif/input.fa --paf-file docs/evaluations/local-compression-testbed-fast/fixtures/microtangle_repeat_motif/smoothxg_control/empty.paf --gfa-engine pggb --threads 1 --min-match-len 1 --min-map-length 1 --target-poa-length 64 --max-node-length 64 --poa-padding-fraction 0 -g docs/evaluations/local-compression-testbed-fast/fixtures/microtangle_repeat_motif/smoothxg_control/output.gfa

# expected_exit_code=0
# command_status=pass
