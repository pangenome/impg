#!/usr/bin/env bash
set -uo pipefail
source /home/erikg/impg-genome-inference/build-env.sh
cd /home/erikg/impg-genome-inference
rm -f build-wd1.exit build-wd1.done
status=0
taskset -c 240-243 nice -n 10 \
  cargo build --release --offline --locked --examples \
  > build-wd1.log 2>&1 || status=$?
echo "$status" > build-wd1.exit
touch build-wd1.done
