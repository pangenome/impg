#!/usr/bin/env bash
set -uo pipefail
source /home/erikg/impg-genome-inference/build-env.sh
cd /home/erikg/impg-genome-inference
taskset -c 240-243 nice -n 10 cargo build --release --offline --locked --example panel_route_mem_routed > build-finalists.log 2>&1
echo $? > build-finalists.exit
echo done > build-finalists.done
