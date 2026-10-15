#!/usr/bin/env bash
cd /home/erikg/impg-genome-inference
source build-env.sh
taskset -c 252-255 cargo build --release --example rc_view_probe > build-rcprobe.log 2>&1
echo $? > build-rcprobe.exit
D=/home/erikg/impg/target/experiments/genome-mem-bwt-pipeline/chrIII-finalist-reranking-v1
taskset -c 240-243 target/release/examples/rc_view_probe \
  /home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng \
  /home/erikg/yeast/genome-panel-routes-rebuild-v1/routes > $D/rc-view-probe.log 2> $D/rc-view-probe.err
echo $? > $D/rc-view-probe.exit
touch rcprobe.done
