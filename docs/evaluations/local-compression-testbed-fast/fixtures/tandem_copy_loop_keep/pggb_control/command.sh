#!/usr/bin/env bash
set -euo pipefail
/usr/bin/time -v /usr/bin/timeout 120s /home/erikg/impg/.wg-worktrees/agent-606/target/debug/impg graph --sequence-files tests/test_data/local_compression/tandem_copy_loop_keep/input.fa --gfa-engine pggb --fastga --threads 1 --min-match-len 1 --min-map-length 1 --target-poa-length 64 --max-node-length 64 --poa-padding-fraction 0 -g docs/evaluations/local-compression-testbed-fast/fixtures/tandem_copy_loop_keep/pggb_control/output.gfa --scaffold-jump 0 --scaffold-mass 1

# expected_exit_code=0
# command_status=pass
