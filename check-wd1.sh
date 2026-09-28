#!/usr/bin/env bash
set -uo pipefail
source /home/erikg/impg-genome-inference/build-env.sh
cd /home/erikg/impg-genome-inference
rm -f check-wd1.exit check-wd1.done
taskset -c 240-243 nice -n 10 \
  cargo check --release --offline --locked --examples \
  > check-wd1.log 2>&1
echo "$?" > check-wd1.exit
touch check-wd1.done
