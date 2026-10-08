#!/usr/bin/env bash
# tests-item2a: the FULL example test suite on the ITEM-1 fixed tree
# (cores 252-255, --test-threads=1, release).
set -uo pipefail
source /home/erikg/impg-genome-inference/build-env.sh
cd /home/erikg/impg-genome-inference
taskset -c 252-255 nice -n 10 cargo test --release --offline --locked --examples -- --test-threads=1
echo $? > tests-item2a.exit
echo done > tests-item2a.done
