#!/usr/bin/env bash
set -euo pipefail
python3 scripts/local_compression_testbed.py run --profile fast --manifest tests/test_data/local_compression/manifest.json --fixtures nested_top_level_right --methods top_flubble_nonoverlap_sweepga --out-dir docs/evaluations/local-compression-testbed-fast

# expected_exit_code=0
# command_status=pass
