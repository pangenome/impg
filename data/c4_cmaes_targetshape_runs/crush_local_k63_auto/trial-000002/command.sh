#!/bin/sh
set -eu

# Primary impg command
/home/erikg/impg/target/release/impg crush -g /home/erikg/impg/data/calibrate_spectrum_driven_20260530T203220Z/c4_default_spectrum/c4-local-syng-k63-spectrum-default.gfa -o data/c4_cmaes_targetshape_runs/crush_local_k63_auto/trial-000002/output.gfa -t 32 --method auto --max-rounds 2 --min-traversal-len 698 --max-traversal-len 26373 --max-median-traversal-len 7572 --max-traversals 6004 --polish-rounds 1 --polish-max-traversal-len 38211 --polish-max-median-traversal-len 669 --replacement-flank-bp 12 --auto-poasta-max-len 44496

# External graph-report scoring command
/home/erikg/impg/target/release/impg graph-report -g data/c4_cmaes_targetshape_runs/crush_local_k63_auto/trial-000002/output.gfa -o data/c4_cmaes_targetshape_runs/crush_local_k63_auto/trial-000002/graph-report.tsv --format tsv --top 20
