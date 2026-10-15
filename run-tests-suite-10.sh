#!/usr/bin/env bash
set -uo pipefail
source /home/erikg/impg-genome-inference/build-env.sh
cd /home/erikg/impg-genome-inference
{
echo "=== suite-10 started $(date) ==="
echo "--- lib ---"
taskset -c 252-255 nice -n 10 cargo test --release --offline --locked --lib -- --test-threads=1 2>&1 | grep -E "^test result"
echo "--- integration (all tests/) ---"
taskset -c 252-255 nice -n 10 cargo test --release --offline --locked --tests -- --test-threads=1 2>&1 | grep -E "Running|^test result"
echo "--- examples ---"
taskset -c 252-255 nice -n 10 cargo test --release --offline --locked --examples -- --test-threads=1 2>&1 | grep -E "Running|^test result"
echo "=== suite-10 finished $(date) ==="
} > tests-suite-10.log 2>&1
echo done > tests-suite-10.done
