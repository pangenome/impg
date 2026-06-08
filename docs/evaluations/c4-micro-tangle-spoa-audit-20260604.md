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

Direct SPOA is not a no-op and should not be rejected because global
white-space/layout metrics increase. In these repetitive C4 micro-sites,
increased white-space is an expected rendering/layout effect when previously
compact repeated sequence is unfolded into lower-coverage local POA structure.
White-space is diagnostic only; it is not a correctness gate.

The hard correctness gate is path spelling. On this gate, direct SPOA is clean.
An accept-all C4 run with wide budgets preserved every path sequence while
accepting every replacement it could materialize.

## Accept-All SPOA Run

Command class:

```text
impg crush --method poa --max-iterations 20 \
  --max-traversal-len 50k --max-median-traversal-len 50k \
  --max-total-sequence 2g --max-traversals 100k
```

Output directory:

```text
data/c4_micro_tangle_accept_all_20260604T102704Z/
```

This run accepted all path-valid, non-overlapping frontiers until convergence:

```text
initial POVU: 2849 sites, 2550 polymorphic candidates, 1747 level-0 roots
round 1: 1786 candidates, 759 selected, 759/759 accepted
round 2: 897 candidates, 618 selected, 618/618 accepted
round 3: 254 candidates, 157 selected, 157/157 accepted
round 4: 97 candidates, 1 selected, 1/1 accepted
round 5: 2934 POVU sites, 0 eligible candidates
total: 1535 resolved, 0 bailed
runtime: 5:43.19
max RSS: 1.6 GB
```

Path validation against the seed graph:

```text
expected_paths    465
observed_paths    465
missing_paths     0
extra_paths       0
spelling_mismatches 0
```

Metrics, diagnostic only:

| graph | segments | bp | bp-weighted coverage | singleton bp | path-jump p99 | white p99 | white max | long white bridges |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| seed | 7411 | 234828 | 454.330293 | 2922 | 5 | 66 | 33624 | 8790 |
| accept-all SPOA | 8387 | 250010 | 426.740826 | 4609 | 3445 | 112207 | 249452 | 214719 |

The increase in white-space/path-jump should be interpreted as a rendering and
local-unfolding signal, not as a reason to reject the replacements.

The sharper diagnosis is:

- POVU sees many micro-tangles.
- The current direct-SPOA loop accepts all materializable, non-overlapping,
  path-valid candidates.
- Many high-depth micro-sites remain visible to POVU but do not become eligible
  candidates after the descent converges.
- Therefore the remaining gap is not a quality/rejection gate. It is the mapping
  from POVU sites to active materializable replacement units, especially around
  repeated boundary nodes and nested/reused local sites.

## Next Step

The next useful experiment is not another median threshold sweep. It is a
candidate-admission audit:

1. classify all final POVU sites into materializable vs non-materializable;
2. for non-materializable sites, record why: non-unique entry/exit, no active
   tree key, nested site already consumed, repeated boundary, or path range
   ambiguity;
3. extend the replacement-unit builder so repeated-boundary micro-sites can be
   grouped into a larger safe unit instead of disappearing from the active
   frontier;
4. keep path spelling as the only hard acceptance gate.
