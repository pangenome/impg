#!/usr/bin/env bash
set -uo pipefail
source /home/erikg/impg-genome-inference/build-env.sh
cd /home/erikg/impg-genome-inference
taskset -c 240-243 nice -n 10 cargo build --release --offline --locked --example panel_route_mem_routed
echo $? > build-cosigt1.exit
echo done > build-cosigt1.done
