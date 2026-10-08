#!/usr/bin/env bash
set -uo pipefail
source /home/erikg/impg-genome-inference/build-env.sh
cd /home/erikg/impg-genome-inference
{
echo "=== suite started $(date) ==="
echo "--- examples ---"
taskset -c 252-255 nice -n 10 cargo test --release --offline --locked --examples -- --test-threads=1 2>&1 | grep -E "Running|^test result|^error"
echo "=== suite finished $(date) ==="
} > tests-suite-16-finalists.log 2>&1
echo done > tests-suite-16-finalists.done
