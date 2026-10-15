#!/usr/bin/env bash
cd /home/erikg/impg-genome-inference
source build-env.sh
timeout 1800 taskset -c 252-255 cargo test --release > tests-collapse1.log 2>&1
echo $? > tests-collapse1.exit
touch tests-collapse1.done
