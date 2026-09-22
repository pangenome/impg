# C4 high-frequency run-aware mask plus aggressive crush sweep

Date: 2026-05-31

Task: `c4-highfreq-mask-crush-sweep`

Output directory:

```text
/home/erikg/impg/data/c4_highfreq_mask_crush_sweep
```

Reproducibility helper:

```text
scripts/c4-highfreq-mask-crush-sweep.py
```

## Summary

The new run-aware high-frequency occurrence mask works mechanically, but the
plain top-fraction sweeps do **not** fix C4. The important failure mode is that
the high-frequency policy can rescue the selected top-frequency occurrences
while the separate scaffold-context/spectrum-glue mask still privately splits
about 373k shared syncmer occurrences. That dominates graph size and leaves the
same underaligned, low-depth topology.

The only useful row in this sweep was the absolute cap:
`top=0,max-occ=1000,freq-run=10,freq-span=1k`. It is a partial lever:
it reduced the graph to 24,168 segments / 691,137 bp and target-shape score
1,580.890. This beats the target-shaped CMA result, but it still does **not**
beat the existing syng-local k63 baseline against the PGGB target distribution.

POA/POASTA pass-through is fixed for these runs: replacement logs show
`replacement build progress ... (accepted)` and compression diagnostics report
`threshold=disabled (diagnostic only); no alternate aligner attempted and no
replacement was rejected by ratio`.

## Commands

All primary metric runs used the requested C4 region, `-d 50k`, the HPRCv2
syng prefix and AGC, and `:nosort`:

```bash
python3 scripts/c4-highfreq-mask-crush-sweep.py \
  --threads 32 \
  --configs default_top0005_run10_span0,top001_run3_span0,top001_run5_span0,top001_run10_span0,top001_run16_span0 \
  --render-top 0

python3 scripts/c4-highfreq-mask-crush-sweep.py \
  --threads 32 \
  --configs top001_run10_span500,top001_run10_span1k,top002_run10_span1k,maxocc1000_run10_span1k \
  --render-top 0

python3 scripts/c4-highfreq-mask-crush-sweep.py \
  --threads 32 \
  --configs top001_run10_span1k_sharedrun10_seq1k,top001_run10_span1k_noskip_polish50k \
  --render-top 0

python3 scripts/c4-highfreq-mask-crush-sweep.py \
  --threads 32 --skip-runs --render-top 3
```

The runner wrote per-row `command.sh`, `engine.txt`, `query.stderr.log`,
`graph-report.tsv`, and `/usr/bin/time -v` output under the data directory.
Primary runs deliberately did not include `:sort,pipeline=Ygs`; only selected
renders used external `gfasort -p Ygs` and `gfalook -m`.

## Matrix

Target metrics are computed against:

```text
/home/erikg/impg/data/local_syng_parameter_sweep_20260529T193452Z/pggb_control/graphs/C4_GRCh38_53kb.gfa
```

Lower target score, frequency TV/EMD, total bp ratio, and singleton bp are
better. All rows completed with 465 paths.

| row | mask/crush change | S | bp | depth | singleton bp | HF split/rescued | scaffold split | TV | EMD | bp ratio | score |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `maxocc1000_run10_span1k` | `top=0,max-occ=1000,freq-run=10,freq-span=1k` | 24,168 | 691,137 | 154.368 | 205,569 | 127 / 1,083,946 | 5,977 | 0.677260 | 0.298795 | 7.74x | 1,580.890 |
| `top001_run10_span1k_sharedrun10_seq1k` | plus generic `min-run=10,sequence-k=1k` | 272,684 | 7,403,381 | 14.411 | 6,931,393 | 5 / 5,186 | 373,022 | 0.949836 | 0.388706 | 82.87x | 4,094.066 |
| `top001_run10_span500` | `top=0.001,freq-run=10,freq-span=500` | 310,063 | 8,494,510 | 12.560 | 7,988,298 | 3 / 5,188 | 372,996 | 0.950780 | 0.389820 | 95.08x | 4,431.161 |
| `top001_run10_span1k` | `top=0.001,freq-run=10,freq-span=1k` | 310,065 | 8,494,516 | 12.560 | 7,988,304 | 5 / 5,186 | 372,996 | 0.950780 | 0.389820 | 95.08x | 4,431.167 |
| `top001_run10_span1k_noskip_polish50k` | no candidate-size skip; polish 50k/10k | 310,065 | 8,494,516 | 12.560 | 7,988,304 | 5 / 5,186 | 372,996 | 0.950780 | 0.389820 | 95.08x | 4,431.167 |
| `top002_run10_span1k` | `top=0.002,freq-run=10,freq-span=1k` | 310,066 | 8,494,532 | 12.560 | 7,988,320 | 9 / 9,540 | 372,996 | 0.950780 | 0.389820 | 95.08x | 4,431.174 |
| `top001_run3_span0` | `top=0.001,freq-run=3,freq-span=0` | 313,092 | 8,532,439 | 12.504 | 8,026,373 | 3,116 / 2,075 | 372,996 | 0.950957 | 0.389855 | 95.50x | 4,449.942 |
| `default_top0005_run10_span0` | current default HF mask | 312,763 | 8,552,235 | 12.475 | 8,046,160 | 2,707 / 20 | 372,997 | 0.951041 | 0.389874 | 95.72x | 4,454.180 |
| `top001_run5_span0` | `top=0.001,freq-run=5,freq-span=0` | 314,877 | 8,575,605 | 12.441 | 8,069,542 | 4,904 / 287 | 372,996 | 0.951139 | 0.389895 | 95.99x | 4,464.477 |
| `top001_run10_span0` | `top=0.001,freq-run=10,freq-span=0` | 315,139 | 8,577,037 | 12.439 | 8,071,089 | 5,171 / 20 | 372,996 | 0.951158 | 0.389897 | 96.00x | 4,464.916 |
| `top001_run16_span0` | `top=0.001,freq-run=16,freq-span=0` | 315,156 | 8,577,073 | 12.439 | 8,071,133 | 5,191 / 0 | 372,996 | 0.951159 | 0.389897 | 96.00x | 4,464.927 |

Source table: `/home/erikg/impg/data/c4_highfreq_mask_crush_sweep/metrics.tsv`.

## Findings

`freq-run` alone is not enough. At `top=0.001`, lowering `freq-run` from 16 to
3 rescued 2,075 of 5,191 high-frequency occurrences, but the final graph was
still 313k segments / 8.53 Mb because scaffold-context splitting stayed at
372,996 occurrences.

`freq-span=500/1k` rescues the selected top-frequency nodes almost completely:
the `top=0.001,freq-run=10,freq-span=1k` row privately split only 5 of 5,191
high-frequency occurrences. The graph still had 310,065 segments because the
generic scaffold-context split remained at 372,996 occurrences.

The generic rescue follow-up (`min-run=10,sequence-k=1k`) proved that exact
sequence support can be detected at this layer (`sequence-supported=4,207,123`),
but it did not stop the spectrum-glue-selected nodes from being split:
scaffold-context split was still 373,022. It improved the top-fraction rows
from 310k to 273k segments, but it was still much worse than the max-occ row.

The no-skip crush follow-up removed the explicit oversized-root skip log and
used `max-traversal-len=0,max-median-traversal-len=0` with higher POASTA polish
caps. It selected 7 round-1 replacements and produced the same metrics as the
ordinary `top001_run10_span1k` row. Candidate-size skipping is therefore not
the primary blocker for this mask shape; the graph is already dominated by
private-copy explosion before crush can help.

The absolute `max-occ=1000` row is the only promising direction in this sweep.
It selected 1,695 high-occurrence nodes covering 1,084,073 local occurrences,
rescued 1,083,946 of them, and privately split only 127 high-frequency
occurrences. More importantly, generic scaffold-context splitting dropped to
5,977 occurrences. This is a partial degluing lever, but still far from PGGB
and still worse than the existing syng-local k63 baseline.

## Baseline Comparison

Against the existing syng-local k63 baseline from
`docs/c4-cmaes-target-shape-results.md`, the new best row does **not** improve:

| graph | score | TV | EMD | bp | depth | singleton bp |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| syng-local k63 baseline | 1,090.028 | 0.670297 | 0.277371 | 259,877 | 194.861 | 26,484 |
| best here: `maxocc1000_run10_span1k` | 1,580.890 | 0.677260 | 0.298795 | 691,137 | 154.368 | 205,569 |

Against the target-shaped CMA result, the new best row **does** improve:

| graph | score | TV | EMD | bp | depth | singleton bp |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| target-shaped full-query syng-local auto | 3,688.987 | 0.947191 | 0.393471 | 7,772,088 | 6.841 | 3,913,290 |
| best here: `maxocc1000_run10_span1k` | 1,580.890 | 0.677260 | 0.298795 | 691,137 | 154.368 | 205,569 |

So the answer is: run-aware high-frequency masking plus aggressive crush is a
real improvement over the failed target-shaped CMA search when configured with
an absolute occurrence cap, but it is not a replacement for the syng-local k63
baseline and it does not reach the PGGB-shaped distribution.

## PNG

Uploaded sorted Ygs render:

- https://hypervolu.me/~erik/impg/c4-highfreq-mask-crush-sweep-maxocc1000_run10_span1k.png

Only the best row was kept as a PNG artifact. The second-ranked known-bad
272k-segment row spent more than 12 minutes in `gfasort` topo-sort, and the
user guidance was to stop chasing worse-row renders. No raw GFA was uploaded.

## Validation

- The sweep used the run-aware high-frequency knobs from
  `syng-mask-highfreq-run-aware`: `freq-run`, `freq-span`,
  `freq-run-aware=true`, plus `max-occ` for the best row.
- 11 real C4 configurations were attempted; 11 produced path-count-valid
  graph-report metrics after the duplicate-parameter no-skip command was fixed.
- All completed metric GFAs have 465 paths and `path_count_ok=true`.
- The best sorted Ygs GFA was exact/path-count valid against its no-sort source:
  `/home/erikg/impg/data/c4_highfreq_mask_crush_sweep/maxocc1000_run10_span1k/path-validation.sorted-vs-nosort.json`
  reports `ok=true`, 465 vs 465 paths, and zero path sequence mismatches.
- `scp` upload of the best PNG was verified with
  `ssh erik@hypervolu.me 'ls -lh www/impg/c4-highfreq-mask-crush-sweep-*.png'`.
