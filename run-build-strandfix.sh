#!/usr/bin/env bash
cd /home/erikg/impg-genome-inference
source build-env.sh
taskset -c 252-255 cargo build --release --example panel_route_mem_routed > build-strandfix.log 2>&1
echo $? > build-strandfix.exit
taskset -c 252-255 cargo test --release --example panel_route_mem_routed > build-strandfix-tests.log 2>&1
echo $? > build-strandfix-tests.exit
touch build-strandfix.done
