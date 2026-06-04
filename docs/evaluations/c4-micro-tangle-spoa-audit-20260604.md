# C4 Micro-Tangle SPOA Audit

Date: 2026-06-04
Branch: `eg/c4-crush-resolution-controls`

## Question

Find micro-tangles on the order of hundreds of bp, check whether the current
POVU/crush path collects them as unique bubbles for resolution, and test whether
direct SPOA fixes them.

## Inputs

- Seed graph:
  `data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_scaffold0.initial.gfa`
- Existing full direct-SPOA series:
  `data/c4_one_many_initial_spoa_true_series_20260603T193935Z/`
- Audit output:
  `data/c4_micro_tangle_audit_20260604T095221/`

## Findings

POVU sees many micro-sites in the seed graph:

```text
POVU flubble decomposition on reference CHM13#0#chr6:31744284-31976975(+):
2849 sites, 2812 leaves
levels: 0:1945, 1:305, 2:50, 3:13, 4:82, 5:77, 6:11, 7:26, 8:74, 9:228, 10:38
```

But many visible micro-sites are not collected as clean unique-anchor bubbles.
The current candidate builder requires a unique entry/exit range per path
(`unique_anchor_range`). Several top POVU micro-sites have entry nodes that
occur thousands of times, so they produce zero usable unique-anchor path ranges.

Examples from the top POVU-site probe:

| site | level | leaf | ref span steps | start occurrences | end occurrences | unique-anchor paths | median bp |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| `>1271>1326` | 0 | true | 21 | 4373 | 465 | 0 | |
| `>7204>7207` | 0 | true | 15 | 9092 | 465 | 0 | |
| `>6361>6364` | 0 | true | 14 | 8740 | 465 | 0 | |
| `>5927>5932` | 0 | true | 10 | 8472 | 464 | 0 | |
| `>1643>1644` | 0 | true | 7 | 522 | 465 | 409 | 230 |
| `>536>544` | 1 | false | 7 | 459 | 459 | 459 | 172 |

So the answer is mixed:

- Some micro-tangles are collected as unique bubbles.
- Many high-depth micro-tangles are seen by POVU but are not resolvable by the
  current unique-anchor candidate model.

## Full-Graph SPOA Probe

A one-round SPOA probe restricted to 200-500 bp candidates selected and accepted
340 replacements:

```text
selected n=340, max-len median/max=293/499, median-len median/max=293/499
340 resolved, 0 bailed
replacement output-bp p50=294, max=500
```

This did not fix the graph. It worsened the same metrics that are visible in the
PNG:

| graph | segments | bp | bp-weighted coverage | singleton bp | path-jump p99 | white p99 | white max | long white bridges |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| seed | 7411 | 234828 | 454.330293 | 2922 | 5 | 66 | 33624 | 8790 |
| SPOA 200-500 r1 | 7725 | 239956 | 444.620989 | 2930 | 3260 | 112076 | 239579 | 215144 |
| full SPOA median100 | 8119 | 239264 | 445.906923 | 3895 | 3454 | 111591 | 238821 | 212876 |
| full SPOA median500 | 8317 | 243630 | 437.915996 | 4745 | 3428 | 111829 | 243072 | 214187 |

## Isolated Micro-Sites

Two concrete 100s-bp sites were extracted and crushed independently.

### `>1643>1644`

Input:

```text
409 paths, 4 nodes, 6 links
median traversal 230 bp, max 232 bp
```

SPOA result:

```text
1 candidate selected, 1 resolved, 0 bailed
segments 4 -> 232
links 6 -> 231
path steps 3017 -> 94224
```

This is not a fix. SPOA unfolds an already compact micro-site into many small
segments.

### `>536>544`

Input:

```text
459 paths, 9 nodes, 11 links
median traversal 172 bp, max 177 bp
```

SPOA result:

```text
2 candidates resolved, 0 bailed
segments 9 -> 8
links 11 -> 10
path steps 3802 -> 2818
white p99 0 -> 80
```

This is only a small local improvement and still introduces local white-space.

## Interpretation

Direct SPOA is not the missing C4 cleanup pass in its current form.

For micro-sites that are already compact in the syng/SweepGA seed, SPOA can
decompact them into per-base or near-per-base POA structure. For other small
sites, it can make a modest local improvement, but this does not translate into
better full-graph structure. On full C4, applying hundreds of these replacements
preserves path spellings but severely worsens path-order and white-space
metrics.

The sharper diagnosis is:

- POVU sees many micro-tangles.
- The current unique-anchor candidate model loses many high-depth micro-sites.
- The micro-sites it does collect are not necessarily the right replacement
  unit for full-graph cleanup.
- SPOA should not be applied blindly to every 100s-bp candidate. It needs either
  a stricter “actually tangled/underaligned” target definition or a different
  context-aware larger-block pass before any small-site polish.

## Next Step

The next useful experiment is not another median threshold sweep. It is a
candidate-classifier audit:

1. classify POVU sites into unique-anchor resolvable vs non-resolvable;
2. identify which class corresponds to visible white-space/tangle artifacts;
3. apply larger context-aware SweepGA/seqwish replacement to the surrounding
   block first;
4. run small-site POA only on residual sites that are not already compact and
   whose local replacement reduces path steps without increasing white-space.

