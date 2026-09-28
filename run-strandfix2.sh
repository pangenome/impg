#!/usr/bin/env bash
cd /home/erikg/impg-genome-inference
source build-env.sh
taskset -c 252-255 cargo build --release --example panel_route_mem_routed > build-strandfix2.log 2>&1
echo $? > build-strandfix2.exit
taskset -c 252-255 cargo test --release --example panel_route_mem_routed > build-strandfix2-tests.log 2>&1
echo $? > build-strandfix2-tests.exit
D=/home/erikg/impg/target/experiments/genome-mem-bwt-pipeline/chrIII-finalist-reranking-v1
cd $D
rm -f run-router-placement-diag2.exit run-router-placement-diag2.done
IMPG_ROUTER_PLACEMENT_DIAG=$PWD/router-placement-diag-spec.json \
taskset -c 240-243 nice -n 10 \
  /home/erikg/impg-genome-inference/target/release/examples/panel_route_mem_routed \
  --panel /home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng \
  --routes /home/erikg/yeast/genome-panel-routes-rebuild-v1/routes \
  --sample /home/erikg/yeast/genome-blind-mosaic-20260914T121834Z/sample/sample.membwt \
  --reads /home/erikg/yeast/genome-blind-mosaic-20260914T121834Z/reads.fastq.gz \
  --axis /home/erikg/yeast/genome-mem-bwt-bootstrap-20260911T055203Z/reference-axis-v1.json \
  --bed-directory /home/erikg/yeast/partition-pos64-w10k-d1k-to-completion/results \
  --component S288C#0#chrIII --routing-universe component \
  > run-router-placement-diag2.log 2> run-router-placement-diag2.err
echo $? > run-router-placement-diag2.exit
touch run-router-placement-diag2.done
