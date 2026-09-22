# C4 CMA-ES Results

Date: 2026-05-31

Task: `run-c4-cmaes`

## Summary

I ran the real `scripts/c4-crush-cmaes.py` CMA-ES wrapper with the installed
Python `cma` package (`optimizer=cma-es`) on the C4 53 kb graph problem. All
trials used 32 impg threads and one trial at a time.

Initial crush-only studies used this existing same-locus syng-local blunt GFA:

```text
/home/erikg/impg/data/calibrate_spectrum_driven_20260530T203220Z/c4_default_spectrum/c4-local-syng-k63-spectrum-default.gfa
```

The PGGB target-shape artifacts were:

```text
/home/erikg/impg/data/c4_coverage_distribution_20260530T2342Z/summary.tsv
/home/erikg/impg/data/c4_coverage_distribution_20260530T2342Z/freq_bins.tsv
/home/erikg/impg/data/local_syng_parameter_sweep_20260529T193452Z/pggb_control/reports/C4_GRCh38_53kb.tsv
```

Exact trial audit rows, graph reports, commands, path validation JSON, and
compressed best-by-study output GFAs are under:

```text
data/c4_cmaes_runs/analysis/
data/c4_cmaes_runs/best_by_study/
```

## Answer

The currently best CMA-ES parameterization by PGGB shape distance is the
full-query `gfa:syng-local` auto trial:

```text
gfa:syng-local:blunt,k=173,s=31,seed=7:mask,top=0.0010159722,max-occ=296,min-run=4,sequence-k=125:crush,method=auto,max-rounds=2,min-traversal-len=55,max-traversal-len=8414,max-median-traversal-len=1512,max-traversals=529,polish-rounds=3,polish-max-traversal-len=8494,polish-max-median-traversal-len=2359,replacement-flank-bp=469,auto-poasta-max-len=4474:sort,pipeline=Ygs
```

It is not good enough to render/upload. It has a segment count close to PGGB
(`7463` vs `7170`), but the similarity is misleading: total segment sequence is
`9,041,867 bp`, which is `101.2x` the PGGB target `89,342 bp`; bp-weighted
depth is `5.88` vs PGGB `595.08`; target-bin frequency total-variation
distance is `0.947`; and path whitespace p99 is `1,376,322 bp` vs PGGB `3 bp`.
This is still a sparse, whitespace-heavy graph with most sequence stored as
singletons.

The best crush-only optimizer score was the `sweepga` trial (`498.463`), but it
is even farther from the PGGB shape (`13,959,198 bp`, `156.2x` PGGB segment
bp). It should not be promoted for render/upload.

For context, the pre-existing non-CMA `syng-local k=63,s=8:crush` run in the
target artifacts remains much closer to PGGB than any CMA trial here:
`259,877 bp` (`2.91x` PGGB), bp-weighted depth `194.86`, singleton bp
`26,484`, and frequency TV distance `0.670`. The CMA search space used by
`scripts/c4-crush-cmaes.py` did not include key current defaults such as
`max-iterations=until-done`, `min-traversal-len=0`, and
`polish-rounds=until-done`, which likely explains why these bounded samples did
not rediscover the current hand-tuned default behavior.

## Ranking Method

The table below is ranked by `PGGB shape distance`, not by raw CMA objective.
Raw CMA objective is still shown. The full-query and crush-only objectives use
different bp scales, so they are not directly comparable across modes.

Frequency distance is total variation distance over the PGGB target
coverage/frequency bins from `freq_bins.tsv`. PGGB shape distance is:

```text
frequency TV
+ abs(log2(segment_bp / PGGB_segment_bp))
+ abs(log2(bp_weighted_depth / PGGB_bp_weighted_depth))
+ 0.25 * abs(log2((path_ws_p99 + 1) / (PGGB_path_ws_p99 + 1)))
```

Lower is better. This scalar is only a ranking aid; the size, depth, and
whitespace columns are the actionable interpretation.

## Ranked Valid Trials

| Rank | Study | Trial | CMA score | PGGB shape distance | Freq TV | Segment bp vs PGGB | BP-weighted depth vs PGGB | Whitespace p99 vs PGGB | Path validation | Wall / RSS | Output GFA |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | --- |
| 1 | full-query syng-local auto | 1 | 81.013 | 18.869 | 0.947 | 101.2x | 0.0099x | 458,774x | pass, 466 paths | 302.3s / 52.5 GiB | `data/c4_cmaes_runs/full_query_syng_local/auto/trial-000001/output.gfa` |
| 2 | crush-only auto | 4 | 504.758 | 18.970 | 0.970 | 85.5x | 0.0117x | 2,199,849x | pass, exact input paths | 109.0s / 32.0 GiB | `data/c4_cmaes_runs/crush_local_k63/auto/trial-000004/output.gfa` |
| 3 | full-query syng auto | 1 | 65.281 | 19.352 | 0.948 | 118.7x | 0.0084x | 488,331x | pass, 466 paths | 353.2s / 54.0 GiB | `data/c4_cmaes_runs/full_query_syng/auto/trial-000001/output.gfa` |
| 4 | crush-only sweepga | 2 | 498.463 | 20.947 | 0.978 | 156.2x | 0.0064x | 4,153,814x | pass, exact input paths | 21.4s / 1.4 GiB | `data/c4_cmaes_runs/crush_local_k63/sweepga/trial-000002/output.gfa` |
| 5 | crush-only chain-povu | 1 | 556.703 | 24.275 | 0.990 | 430.7x | 0.0023x | 12,256,236x | pass, exact input paths | 54.1s / 1.2 GiB | `data/c4_cmaes_runs/crush_local_k63/chain-povu/trial-000001/output.gfa` |

Compressed copies of these best-by-study output GFAs are in
`data/c4_cmaes_runs/best_by_study/`.

## Raw Metrics

Key raw `impg graph-report --format tsv` fields for the ranked rows:

| Study | Segments | Links | Paths | Path steps | Segment bp | BP-weighted depth | Singleton bp | WS p99 | WS max | WS bridges >= threshold | Duplicate seq frac | Local repeat nodes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| PGGB target | 7,170 | 8,664 | 466 | 3,565,435 | 89,342 | 595.083 | 523 | 3 | 40,274 | 9,607 | 0.818550 | 356 |
| full-query syng-local auto | 7,463 | 9,599 | 466 | 143,359 | 9,041,867 | 5.876 | 4,754,103 | 1,376,322 | 6,705,560 | 30,437 | 0.324668 | 464 |
| crush-only auto | 280,047 | 354,199 | 466 | 1,937,609 | 7,637,777 | 6.957 | 7,400,077 | 6,599,546 | 7,595,130 | 227,150 | 0.985799 | 548 |
| full-query syng auto | 18,698 | 24,974 | 466 | 95,793 | 10,608,051 | 5.009 | 5,985,488 | 1,464,992 | 8,279,760 | 22,134 | 0.761793 | 230 |
| crush-only sweepga | 282,001 | 358,684 | 466 | 1,722,424 | 13,959,198 | 3.806 | 13,590,141 | 12,461,442 | 13,915,890 | 223,766 | 0.985135 | 491 |
| crush-only chain-povu | 332,388 | 442,496 | 466 | 911,053 | 38,475,316 | 1.381 | 38,415,958 | 36,768,708 | 38,405,991 | 249,402 | 0.993080 | 115 |

Full raw rows are in `data/c4_cmaes_runs/analysis/trial_summary.tsv`.

## Study Outcomes

| Study | Trials | Result | Best command / params | Score | Path validation | Output |
| --- | ---: | --- | --- | ---: | --- | --- |
| crush-only auto | 4 | 4 valid | `data/c4_cmaes_runs/crush_local_k63/auto/trial-000004/command.sh` | 504.758 | pass, exact path spelling | `data/c4_cmaes_runs/crush_local_k63/auto/trial-000004/output.gfa` |
| crush-only sweepga | 2 | 1 valid, 1 timeout | `data/c4_cmaes_runs/crush_local_k63/sweepga/trial-000002/command.sh` | 498.463 | pass, exact path spelling | `data/c4_cmaes_runs/crush_local_k63/sweepga/trial-000002/output.gfa` |
| crush-only chain-povu | 2 | 2 valid | `data/c4_cmaes_runs/crush_local_k63/chain-povu/trial-000001/command.sh` | 556.703 | pass, exact path spelling | `data/c4_cmaes_runs/crush_local_k63/chain-povu/trial-000001/output.gfa` |
| crush-only iterative-multi-level | 2 | pruned | trial 1 timed out at 900s; trial 2 was manually terminated after repeating the slow replacement-build path | invalid | fail, missing output GFA | none |
| crush-only coverage-multi-bubble | 2 | pruned | `--multi-level-window-mode sibling`; both candidates timed out at 300s | invalid | fail, missing output GFA | none |
| full-query syng-local auto | 1 | 1 valid | `data/c4_cmaes_runs/full_query_syng_local/auto/trial-000001/command.sh` | 81.013 | pass, 466 output/report paths | `data/c4_cmaes_runs/full_query_syng_local/auto/trial-000001/output.gfa` |
| full-query syng auto | 1 | 1 valid | `data/c4_cmaes_runs/full_query_syng/auto/trial-000001/command.sh` | 65.281 | pass, 466 output/report paths | `data/c4_cmaes_runs/full_query_syng/auto/trial-000001/output.gfa` |

## Notes

- Sweepga had the best crush-only objective, but one candidate hit the 900s
  timeout and the valid candidate was worse than auto by PGGB shape.
- Iterative-multi-level and coverage-multi-bubble repeatedly entered slow
  replacement-building paths and were pruned as requested by the task.
- Auto was the only family robust enough to promote to full-query mode in this
  run. Both `gfa:syng-local` and `gfa:syng` full-query checks completed with
  valid path counts.
- None of the CMA results should be rendered/uploaded as a current best C4
  graph. The current practical recommendation remains to keep the existing
  default syng-local/global crush outputs as the comparison baseline and expand
  the CMA search space before spending more compute.
