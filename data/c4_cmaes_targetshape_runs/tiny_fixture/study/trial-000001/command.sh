#!/bin/sh
set -eu

# Primary impg command
data/c4_cmaes_targetshape_runs/tiny_fixture/fake-impg crush -g data/c4_cmaes_targetshape_runs/tiny_fixture/input.gfa -o data/c4_cmaes_targetshape_runs/tiny_fixture/study/trial-000001/output.gfa -t 1 --method auto --max-rounds 4 --min-traversal-len 37 --max-traversal-len 6280 --max-median-traversal-len 751 --max-traversals 3940 --polish-rounds 1 --polish-max-traversal-len 29289 --polish-max-median-traversal-len 462 --replacement-flank-bp 51 --auto-poasta-max-len 6188

# External graph-report scoring command
data/c4_cmaes_targetshape_runs/tiny_fixture/fake-impg graph-report -g data/c4_cmaes_targetshape_runs/tiny_fixture/study/trial-000001/output.gfa -o data/c4_cmaes_targetshape_runs/tiny_fixture/study/trial-000001/graph-report.tsv --format tsv --top 20
