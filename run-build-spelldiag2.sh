#!/usr/bin/env bash
cd /home/erikg/impg-genome-inference
source build-env.sh
taskset -c 252-255 cargo build --release --example panel_route_mem_routed > build-spelldiag2.log 2>&1
echo $? > build-spelldiag2.exit
touch build-spelldiag2.done
