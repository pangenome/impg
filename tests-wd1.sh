#!/usr/bin/env bash
set -uo pipefail
source /home/erikg/impg-genome-inference/build-env.sh
cd /home/erikg/impg-genome-inference
rm -f tests-wd1.exit tests-wd1.done
status=0
taskset -c 252-255 nice -n 10 \
  cargo test --release --offline --locked --lib -- --test-threads=1 \
  > tests-wd1.log 2>&1 || status=$?
taskset -c 252-255 nice -n 10 \
  cargo test --release --offline --locked --test test_panel_route_diploid_search -- --test-threads=1 \
  >> tests-wd1.log 2>&1 || status=$?
taskset -c 252-255 nice -n 10 \
  cargo test --release --offline --locked --test test_panel_route_window_domain -- --test-threads=1 \
  >> tests-wd1.log 2>&1 || status=$?
taskset -c 252-255 nice -n 10 \
  cargo test --release --offline --locked --examples -- --test-threads=1 \
  >> tests-wd1.log 2>&1 || status=$?
echo "$status" > tests-wd1.exit
touch tests-wd1.done
