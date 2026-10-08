#!/usr/bin/env bash
set -uo pipefail
source /home/erikg/impg-genome-inference/build-env.sh
cd /home/erikg/impg-genome-inference
taskset -c 240-243 nice -n 10 cargo check --release --offline --locked --example panel_route_mem_routed 2>&1 | tail -80
echo ${PIPESTATUS[0]} > build-item2a.exit
echo done > build-item2a.done
