#!/bin/sh
set -eu

# Primary impg command
/home/erikg/impg/target/release/impg crush -g /home/erikg/impg/data/calibrate_spectrum_driven_20260530T203220Z/c4_default_spectrum/c4-local-syng-k63-spectrum-default.gfa -o data/c4_cmaes_targetshape_runs/crush_local_k63_auto/trial-000001/output.gfa -t 32 --method auto --max-rounds 2 --min-traversal-len 59 --max-traversal-len 118267 --max-median-traversal-len 1100 --max-traversals 706 --polish-rounds 3 --polish-max-traversal-len 6240 --polish-max-median-traversal-len 146 --replacement-flank-bp 6 --auto-poasta-max-len 10064

# External graph-report scoring command
/home/erikg/impg/target/release/impg graph-report -g data/c4_cmaes_targetshape_runs/crush_local_k63_auto/trial-000001/output.gfa -o data/c4_cmaes_targetshape_runs/crush_local_k63_auto/trial-000001/graph-report.tsv --format tsv --top 20
