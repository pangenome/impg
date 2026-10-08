#!/usr/bin/env bash
set -uo pipefail
source /home/erikg/impg-genome-inference/build-env.sh
cd /home/erikg/impg-genome-inference
{
taskset -c 252-255 nice -n 10 cargo test --release --offline --locked --examples -- --test-threads=1 2>&1 | grep -E "Running|^test result|^error"
} > tests-cosigt2.log 2>&1
echo $? > tests-cosigt2.exit
echo done > tests-cosigt2.done
