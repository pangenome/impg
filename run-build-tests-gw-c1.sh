#!/usr/bin/env bash
# build+tests for the genome-wide cross-support background (gw-c1):
# release build of the touched example, then the three test groups
# (lib / tests / examples) under the campaign convention (taskset 252-255,
# single-threaded, offline, locked).
set -uo pipefail
source /home/erikg/impg-genome-inference/build-env.sh
cd /home/erikg/impg-genome-inference
rm -f build-gw-c1.exit build-gw-c1.done
status=0
taskset -c 252-255 nice -n 10 \
  cargo build --release --offline --locked --examples \
  > build-gw-c1.log 2>&1 || status=$?
echo "$status" > build-gw-c1.exit
if [ "$status" -ne 0 ]; then touch build-gw-c1.done; exit "$status"; fi
rm -f tests-gw-c1.exit
taskset -c 252-255 nice -n 10 \
  cargo test --release --offline --locked --lib -- --test-threads=1 \
  > tests-gw-c1.log 2>&1 || status=$?
taskset -c 252-255 nice -n 10 \
  cargo test --release --offline --locked --tests -- --test-threads=1 \
  >> tests-gw-c1.log 2>&1 || status=$?
taskset -c 252-255 nice -n 10 \
  cargo test --release --offline --locked --examples -- --test-threads=1 \
  >> tests-gw-c1.log 2>&1 || status=$?
echo "$status" > tests-gw-c1.exit
touch build-gw-c1.done
