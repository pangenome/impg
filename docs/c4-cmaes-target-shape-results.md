# C4 CMA-ES Target-Shape Results

Date: 2026-05-31

Task: `pggb-target-shape-objective-v2`

## Summary

I changed `scripts/c4-crush-cmaes.py` so the PGGB target spectrum is part of
the objective sent back to CMA-ES. With `--target-gfa`, the wrapper now parses
each valid trial GFA and scores these target terms before assigning
`raw_objective` and `tell_objective`:

- bp-weighted node copy-frequency distribution, where node depth is occurrence
  count divided by path count, using the supplied PGGB frequency bins.
- Ordered frequency distance (`target_frequency_emd`) plus total variation
  (`target_frequency_tv`). This is the dominant default term:
  `1000 * (TV + EMD)`.
- bp-weighted node length distribution distance.
- excess total segment bp relative to the PGGB target.
- the existing path-whitespace penalties at a much smaller weight.

Path spelling/path-count validation remains a hard invalid constraint. The
target-shape code is only in the optimizer wrapper; no `impg` quality gate was
added.

Target data:

```text
/home/erikg/impg/data/local_syng_parameter_sweep_20260529T193452Z/pggb_control/graphs/C4_GRCh38_53kb.gfa
/home/erikg/impg/data/c4_coverage_distribution_20260530T2342Z/freq_bins.tsv
/home/erikg/impg/data/c4_coverage_distribution_20260530T2342Z/summary.tsv
```

Artifacts:

```text
data/c4_cmaes_targetshape_runs/analysis/trial_summary.tsv
data/c4_cmaes_targetshape_runs/analysis/comparison_summary.tsv
data/c4_cmaes_targetshape_runs/best_by_study/
```

## Reruns

This was a bounded real C4 rerun, not only a fixture:

| Study | Trials | Valid | Best trial | Target-shaped score | Path validation |
| --- | ---: | ---: | ---: | ---: | --- |
| crush-only auto from syng-local k63 | 2 | 2 | 1 | 4069.024 | exact 466 input paths, no spelling mismatches |
| full-query syng-local auto | 1 | 1 | 1 | 3688.987 | 466 expected, observed, and reported paths |

The full-query best configuration was:

```text
gfa:syng-local:blunt,k=163,s=22,seed=7:mask,top=0.00014658303,max-occ=377,min-run=8,sequence-k=147:crush,method=auto,max-rounds=3,min-traversal-len=2,max-traversal-len=3330,max-median-traversal-len=3330,max-traversals=874,polish-rounds=3,polish-max-traversal-len=11642,polish-max-median-traversal-len=2822,replacement-flank-bp=192,auto-poasta-max-len=1771:sort,pipeline=Ygs
```

The crush-only best configuration was:

```text
crush --method auto --max-rounds 2 --min-traversal-len 59 --max-traversal-len 118267 --max-median-traversal-len 1100 --max-traversals 706 --polish-rounds 3 --polish-max-traversal-len 6240 --polish-max-median-traversal-len 146 --replacement-flank-bp 6 --auto-poasta-max-len 10064
```

## Target Distances

Lower is better for the target-shaped score and all distance columns. Total bp
ratio is relative to the PGGB target (`89,342 bp`).

| Graph / run | Target score | Freq TV | Freq EMD | Length distance | Total bp | BP ratio | BP-weighted depth | Singleton bp | WS p99 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| PGGB target | 23.738 | 0.000000 | 0.000000 | 0.000000 | 89,342 | 1.00x | 595.083 | 523 | 3 |
| syng-local k63 baseline | 1090.028 | 0.670297 | 0.277371 | 0.655807 | 259,877 | 2.91x | 194.861 | 26,484 | 644 |
| previous CMA full-query syng-local auto | 4044.256 | 0.947302 | 0.394033 | 1.268603 | 9,041,867 | 101.21x | 5.876 | 4,754,103 | 1,376,322 |
| target-shaped full-query syng-local auto | 3688.987 | 0.947191 | 0.393471 | 1.243618 | 7,772,088 | 86.99x | 6.841 | 3,913,290 | 1,806,092 |
| target-shaped crush-only auto, trial 1 | 4069.024 | 0.970412 | 0.393640 | 0.847226 | 7,657,292 | 85.71x | 6.939 | 7,420,333 | 6,664,767 |
| target-shaped crush-only auto, trial 2 | 4109.992 | 0.972566 | 0.393718 | 0.851953 | 7,761,275 | 86.87x | 6.846 | 7,556,067 | 6,766,375 |

The PGGB target itself has a nonzero score because the objective still includes
small residual graph-report penalties, such as path whitespace and duplicate
sequence diagnostics. Its target-shape distances are zero.

## Comparison

Against the previous CMA result reported in `docs/c4-cmaes-results.md`, the
bounded target-shaped full-query run improved under the new target-shaped
objective:

- score improved from `4044.256` to `3688.987`;
- total segment bp improved from `9,041,867` to `7,772,088`;
- bp-weighted depth improved from `5.876` to `6.841`;
- singleton bp improved from `4,754,103` to `3,913,290`;
- frequency TV improved slightly from `0.947302` to `0.947191`;
- frequency EMD improved from `0.394033` to `0.393471`;
- length distance improved from `1.268603` to `1.243618`.

This is a real improvement over the previous CMA run, but it still does not beat
the existing syng-local k63 baseline. The baseline remains much closer to PGGB:
`259,877 bp` total (`2.91x` PGGB), bp-weighted depth `194.861`, singleton bp
`26,484`, frequency TV `0.670297`, and target-shaped score `1090.028`.

The result is therefore corrective for the optimizer objective, but not a new C4
graph recommendation. The next useful search should spend more samples around
the full-query syng-local region and expand the CMA dimensions to include the
current hand-tuned defaults that remain outside this bounded search space.

## Validation Notes

- `python -m py_compile scripts/c4-crush-cmaes.py` passed.
- `scripts/c4-crush-cmaes.py --help` documents `--target-gfa`,
  `--target-freq-bins`, `--target-summary`, and the four target weights.
- A tiny fixture under `data/c4_cmaes_targetshape_runs/tiny_fixture/` produced
  nonzero `target_frequency_weighted`, `target_length_weighted`, and
  `target_total_bp_weighted` terms in `trial.json`; `raw_objective` equaled the
  weighted total.
- Best-by-study gzip artifacts under
  `data/c4_cmaes_targetshape_runs/best_by_study/` passed `gzip -t`, along with
  the existing `data/c4_cmaes_runs/best_by_study/*.gfa.gz` artifacts.
