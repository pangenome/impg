#!/bin/sh
set -eu

# Primary impg command
/home/erikg/impg/target/release/impg query -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng -d 100k -r 'GRCh38#0#chr6:31982056-32035418' --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -o gfa:syng-local:blunt,k=163,s=22,seed=7:mask,top=0.00014658303,max-occ=377,min-run=8,sequence-k=147:crush,method=auto,max-rounds=3,min-traversal-len=2,max-traversal-len=3330,max-median-traversal-len=3330,max-traversals=874,polish-rounds=3,polish-max-traversal-len=11642,polish-max-median-traversal-len=2822,replacement-flank-bp=192,auto-poasta-max-len=1771:sort,pipeline=Ygs -O data/c4_cmaes_targetshape_runs/full_query_syng_local_auto/trial-000001/output -t 32

# External graph-report scoring command
/home/erikg/impg/target/release/impg graph-report -g data/c4_cmaes_targetshape_runs/full_query_syng_local_auto/trial-000001/output.gfa -o data/c4_cmaes_targetshape_runs/full_query_syng_local_auto/trial-000001/graph-report.tsv --format tsv --top 20
