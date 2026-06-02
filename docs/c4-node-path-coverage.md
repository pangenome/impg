# C4 Node Path-Coverage Comparison

Task: `add-node-path`

This note adds node path-coverage metrics to the C4 crush comparison. The table
uses real generated C4 GFAs from `/home/erikg/impg/data`; no mocked graph inputs
were used.

## Metric Definition

Node coverage here means **path-step visits per segment**:

- For each `P` or `W` line, every segment step increments that node's coverage.
- If one path visits the same node twice, that contributes two visits.
- `mean_cov = total path-step visits / segment_count`.
- `bp_weighted_cov = sum(segment_len * path_step_visits) / sum(segment_len)`.

The bp-weighted value is also the haplotype-base depth through stored segment
bases. For graph condensation comparisons it is the most useful single coverage
number: it rewards long shared sequence and discounts cases where many tiny
high-depth nodes inflate the unweighted mean.

The high-coverage shared-node threshold in this run is coverage `>=233`, because
all compared graphs have 465 paths and the report default is 50% of paths.

## Comparison Table

| Run | Segments | Links | Segment bp | Trivial stringy | Path steps | Mean cov | bp-weighted cov | p10 / median / p90 | Singleton nodes / bp | High-cov nodes / bp |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| iterative-multi baseline | 16,547 | 23,610 | 396,438 | 103 | 6,390,646 | **386.21** | 256.37 | 2 / 398 / 745 | 1,444 / 31,048 | 10,076 / 187,140 |
| neighbor-merge | 18,761 | 25,763 | 461,241 | **72** | 4,264,049 | 227.28 | 220.35 | 2 / 177 / 465 | 1,803 / 31,769 | 8,170 / 201,542 |
| iterative-multi-level final1 | 16,407 | 23,200 | 391,391 | 101 | 6,098,683 | 371.71 | 259.67 | 2 / 398 / 745 | 1,455 / 31,057 | 9,996 / 188,522 |
| expand-multi-bubble best | **16,227** | **22,882** | 392,671 | 100 | 5,696,957 | 351.08 | 258.83 | 2 / 398 / 745 | 1,461 / 31,056 | 9,803 / 189,778 |
| diagnose-residual-two best | 16,477 | 23,288 | **384,319** | 100 | 5,945,457 | 360.83 | **264.45** | 2 / 398 / 745 | 1,453 / 30,925 | 10,057 / 189,151 |
| low-min default-off + POASTA smoothing | 21,346 | 26,845 | 570,842 | 191 | 4,403,601 | 206.30 | 178.04 | 1 / 100 / 465 | 3,164 / 85,938 | 9,676 / 220,291 |

## Rankings

By unweighted mean node coverage:

1. iterative-multi baseline: 386.21
2. iterative-multi-level final1: 371.71
3. diagnose-residual-two best: 360.83
4. expand-multi-bubble best: 351.08
5. neighbor-merge: 227.28
6. low-min default-off + POASTA smoothing: 206.30

By bp-weighted mean path-depth:

1. diagnose-residual-two best: 264.45
2. iterative-multi-level final1: 259.67
3. expand-multi-bubble best: 258.83
4. iterative-multi baseline: 256.37
5. neighbor-merge: 220.35
6. low-min default-off + POASTA smoothing: 178.04

## Interpretation

The unweighted mean does **not** match the visual-quality story by itself. It
puts the iterative-multi baseline first because that graph has the most path
steps per segment, even though later C4 outputs are visually and structurally
better by segment count, segment bp, or residual-bubble inspection.

The bp-weighted mean is more aligned with graph condensation. It downranks
neighbor-merge and low-min smoothing, both of which have larger stored sequence
and weaker shared representation despite neighbor-merge's low trivial-stringy
count. It also lifts `diagnose-residual-two best`, whose accepted high-bp
candidate reduced segment bp enough to produce the strongest haplotype-base
depth through stored sequence.

This still does not perfectly reproduce visual quality. The
`expand-multi-bubble` graph remains the best segment-count output and was the
best prior visual among the multi-bubble runs, while `diagnose-residual-two`
wins bp-weighted coverage by reducing segment bp with a different accepted
candidate. Use bp-weighted coverage as a primary condensation metric alongside
segment bp, singleton bp, white-space/path-jump metrics, and visual inspection;
do not use unweighted mean coverage as a standalone quality score.

## Source GFAs

| Run | GFA |
| --- | --- |
| iterative-multi baseline | `/home/erikg/impg/data/c4_crush_iterative_multi_20260527T050056Z/run.Ygs.gfa` |
| neighbor-merge | `/home/erikg/impg/data/c4_iterative_multi_level_20260528T235244Z/baseline/neighbor-merge.Ygs.gfa` |
| iterative-multi-level final1 | `/home/erikg/impg/data/c4_iterative_multi_level_20260528T235244Z/sibling_final1_iterative/run.Ygs.gfa` |
| expand-multi-bubble best | `/home/erikg/impg/data/c4_expand_multi_bubble_20260529T020548Z/sibling_prior_shape_fb_r3/run.Ygs.gfa` |
| diagnose-residual-two best | `/home/erikg/impg/data/c4_diagnose_residual_two_20260529T141339Z/outward_highbp_r3/run.Ygs.gfa` |
| low-min default-off + POASTA smoothing | `/home/erikg/impg/data/c4_low_min_match_20260528T163541Z/default_off_poasta_smooth/run.Ygs.gfa` |

## Limitations

- Path-step coverage is not distinct haplotype support. A repeated visit by the
  same path increments coverage more than once.
- Orientation is ignored for coverage aggregation; `node+` and `node-` both
  contribute to the same segment's coverage.
- Zero-length or `*` segments contribute to node counts and path-step coverage
  but not to bp-weighted coverage.
- A high bp-weighted mean can still coexist with bad topology if shared sequence
  is connected by long path jumps or white-space bridges.
- Singleton bp is useful context: private alleles are real in C4, but high
  singleton bp after crushing can also mean underalignment or replacement
  fragmentation.
