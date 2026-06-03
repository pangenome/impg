#!/bin/sh
set -eu
/home/erikg/.cargo/bin/impg crush --gfa /home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_minmatch1_scaffold0.initial.gfa --output /home/erikg/impg/.wg-worktrees/agent-450/data/c4_threshold_sweep_20260603T154441Z/poasta_median5000/poasta_median5000.gfa --method poasta --max-iterations 20 --max-span 0 --max-traversal-len 100k --max-median-traversal-len 5000 --max-total-sequence 500m --max-traversals 100k --poa-scoring 1,4,6,2,26,1 --threads 32 -v 1
