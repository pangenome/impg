# Diagnose residual two-cluster C4 bubbles

Task: `diagnose-residual-two`

Input graph:
`/home/erikg/impg/data/c4_expand_multi_bubble_20260529T020548Z/sibling_prior_shape_fb_r3/run.nosort.gfa`

Input visual:
`/home/erikg/impg/data/c4_expand_multi_bubble_20260529T020548Z/sibling_prior_shape_fb_r3/c4-expand-multi-bubble-best.png`

Baseline metrics: 16227 segments, 22882 links, 465 paths, 392671 segment bp,
100 trivial stringy nodes.

Output directory:
`/home/erikg/impg/data/c4_diagnose_residual_two_20260529T141339Z`

## Summary

The largest residual A/B regions are exposed by POVU/flubble. The main misses
were window-generation and objective-selection issues, not a failure to find
the residuals.

Two important cases:

- `>6987..>7045` is exposed as a level-0 POVU site, with nested sites
  `>6933..>6934` at level 2 and `>7226..>6886` at level 3. The old outward
  windows did not force a candidate across this full residual because the
  30000 bp target was just below the 31218 bp root span and the generator did
  not admit large single residual sites directly. After adding direct residual
  admission and running a 42000 bp target, candidates spanning the region were
  generated and built. SweepGA alignments were present and replacements were
  valid, but all local objectives were negative.
- `>2394..>4224` is exposed as a level-0 POVU site with a 41428 bp root span
  and 16.46 Mb of traversal sequence. It was intentionally not built in the
  focused runs because it exceeds the tested 6 Mb total-traversal cap. This is
  currently a runtime/cap boundary, not a flubble detection failure.

One outward high-bp run did find an acceptable SweepGA condensation in a
different residual high-bp region (`>538..>1021` around `>544..>1004`), reducing
segment bp by 8352 bp while preserving all 465 paths. It increased segment count
by 250, so it improves the current objective but does not reduce the headline
segment count.

## Code Changes

Added diagnostics and an explicit outward mode for iterative multi-level crush:

- `--window-mode outward` / `window-mode=outward`
- candidate source: `outward-residual-window`
- round logs for largest discovered residual sites
- round logs for top generated candidate boundaries and traversal sizes
- build logs with SweepGA evidence:
  - raw PAF record count
  - aligned bp
  - PAF byte count
  - fastga frequency
  - replacement segment count
  - replacement shared-segment count
  - replacement bp
  - objective deltas and rejection reason
- regression coverage:
  - `iterative_multi_level_generates_outward_residual_windows`
  - `test_gfa_output_format_accepts_outward_window_mode`

## Residual Sites

POVU/flubble discovery on the current best C4 output found 4613 sites, including
the major residuals below. These entries are from the new diagnostic logs.

| Site | POVU level | Parent | Root span | Traversals | Total traversal bp | Note |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| `>2394..>4224` | 0 | `.` | 41428 bp | 464 | 16459883 | Large exposed level-0 residual; exceeds tested traversal cap. |
| `>6987..>7045` | 0 | `.` | 31218 bp | 372 | 4884981 | Large exposed residual; direct outward candidates were later generated and built. |
| `>6952..>6953` | 0 | `.` | 31187 bp | 57 | 1197 | Thin level-0 sibling near the large `>6987..>7045` residual. |
| `>6933..>6934` | 2 | `>6987>7045` | 29283 bp | 56 | 392 | Nested inside `>6987..>7045`. |
| `>7226..>6886` | 3 | `>6933>6934` | 28866 bp | 197 | 4872128 | Nested A/B-like high-bp residual inside the `>6987..>7045` tree. |
| `>544..>1004` | 0 | `.` | 5787 bp | 464 | 3528804 | High-bp residual area where outward_highbp accepted one candidate. |

## Run Results

| Run | Iterations requested | Runtime | Peak RSS | Output metrics | Accepted | R&D acceptability |
| --- | ---: | ---: | ---: | --- | ---: | --- |
| `outward_focus_r3` | 3 | 0:59.32 | 3.22 GB | 16227 S / 392671 bp / 100 stringy | 0 | Acceptable. Built 12 candidates quickly, but selected global objective failed. |
| `outward_highbp_r3` | 3 | 2:32.56 | 3.84 GB | 16477 S / 384319 bp / 100 stringy | 1 | Acceptable for iterative R&D. Round 1 improved bp; round 2 stopped after negative local objectives. |
| `outward_direct_residual_r1` | 1 | 0:38.78 | 3.01 GB | 16227 S / 392671 bp / 100 stringy | 0 | Acceptable focused diagnostic. Built direct `>6987..>7045` candidates and explained rejection. |

No run timed out or stalled. The 2.5 minute high-bp run is within the observed
interactive R&D envelope for 2-5 rounds on this C4 locus.

## Candidate Evidence

### `outward_focus_r3`

Purpose: SweepGA-only outward windows with a 30000 bp target, 10 max sites, 12
candidate limit, 3 iterations.

Observed problem: the generated top candidates were smaller path-order windows
near root steps 9660..10567. They did not include the large
`>6987..>7045` residual because its 31218 bp root span was above the 30000 bp
target and the previous outward generator did not admit one large residual site
as a direct candidate.

Evidence:

- built candidates: 12
- local objective passing: 8
- local objective rejected: 4
- cheap candidates selected before global check: 2
- global objective delta: -2565694, below floor 1
- rejection: global objective would regress the graph to 16334 segments,
  395353 segment bp, 101 stringy nodes
- final accepted: 0

Representative built candidate:

- boundary: `>13550..>17022`
- source: `outward-residual-window`
- sites: 10
- method: SweepGA
- root span: 1405 bp
- raw PAF records: 215296
- aligned bp: 302993674
- PAF bytes: 28303900
- replacement: 123 segments, 1434 bp
- objective score delta: 65189 locally
- acceptance: rejected by global objective when combined with the other selected
  non-overlapping candidate

### `outward_highbp_r3`

Purpose: prioritize high traversal-bp residual windows, 30000 bp target, 12 max
sites, 8 candidate limit, 3 iterations.

Accepted candidate:

- boundary: `>538..>1021`
- source: `outward-residual-window`
- sites: 12
- method: SweepGA
- root steps: 448..650
- root span: 6406 bp
- traversals: 465
- total traversal bp: 3825062
- raw PAF records: 223565
- aligned bp: 1763885766
- PAF bytes: 39411008
- fastga frequency: 1000000
- replacement: 620 segments, 594 shared segments, 8045 bp
- objective score delta: +8299948
- segment delta: -250 by the objective's signed convention
- segment bp delta: +8352
- global objective delta: +8299948
- output: 16477 segments, 384319 segment bp, 100 stringy nodes
- path validation: 465/465 paths, 0 spelling mismatches

Round 2 built 8 outward candidates after the round-1 rewrite, but all had
negative local objective deltas. No second candidate was accepted.

### `outward_direct_residual_r1`

Purpose: explicitly target the previously missed large residual with a 42000 bp
target and direct single-site residual admission.

Generated candidates all spanned `>6987..>7045` plus nearby context:

| Candidate | Sites | Root span | Traversals | Total traversal bp |
| --- | ---: | ---: | ---: | ---: |
| `>6919..>7045` | 12 | 31345 bp | 233 | 5139619 |
| `>6921..>7045` | 11 | 31342 bp | 233 | 5138920 |
| `>6929..>7045` | 10 | 31328 bp | 233 | 5135658 |
| `>6931..>7045` | 9 | 31314 bp | 233 | 5132396 |

All four candidates had SweepGA alignment evidence and non-empty replacements.
They were rejected by the local objective:

| Candidate | Raw PAF records | Aligned bp | PAF bytes | Replacement | Objective score delta | Segment bp delta | Reason |
| --- | ---: | ---: | ---: | --- | ---: | ---: | --- |
| `>6919..>7045` | 117159 | 877306254 | 25969121 | 1296 S / 31640 bp | -1093116 | -1404 | local objective rejected |
| `>6921..>7045` | 117159 | 877004648 | 25963956 | 1298 S / 31639 bp | -1096408 | -1406 | local objective rejected |
| `>6929..>7045` | 116923 | 875557859 | 25739916 | 1298 S / 31635 bp | -1099440 | -1408 | local objective rejected |
| `>6931..>7045` | 116545 | 873532898 | 25462081 | 1297 S / 31649 bp | -1126523 | -1436 | local objective rejected |

This answers the central question for the visible `>6987..>7045`/nested
`>7226..>6886` residual: it is not an alignment-absence problem and not an
invalid-replacement problem. The outward-expanded candidate can be generated and
SweepGA can build a replacement, but the induced graph is not smaller under the
current objective.

## Path Spelling Validation

The three real diagnostic outputs were compared against the prior-best input
graph by spelling every GFA path from segment sequence and orientation.

| Run | Paths | Missing | Extra | Spelling mismatches |
| --- | ---: | ---: | ---: | ---: |
| `outward_focus_r3` | 465/465 | 0 | 0 | 0 |
| `outward_highbp_r3` | 465/465 | 0 | 0 | 0 |
| `outward_direct_residual_r1` | 465/465 | 0 | 0 | 0 |

## Best PNG

Best diagnostic output by segment-bp reduction:

`/home/erikg/impg/data/c4_diagnose_residual_two_20260529T141339Z/outward_highbp_r3/c4-diagnose-residual-two-best.png`

Uploaded:

`https://hypervolu.me/impg/c4-diagnose-residual-two-best.png`

Remote confirmation:

`ssh erik@hypervolu.me 'ls -lh www/impg/c4-diagnose-residual-two-best.png'`

Result:

`-rw-r--r-- 1 erik erik 965K May 29 14:40 www/impg/c4-diagnose-residual-two-best.png`

## Conclusion

Residual A/B bubbles are not primarily a flubble decomposition problem. POVU
exposes the major residuals, including the large level-0 site and nested sites.

There are two distinct remaining blockers:

1. Window generation boundaries were too conservative for the `>6987..>7045`
   residual. A 30000 bp target did not include a 31218 bp level-0 residual, and
   the generator did not admit large single residual sites directly. The new
   outward residual mode and direct residual admission expose this path.
2. Once generated, the major `>6987..>7045` candidates have real SweepGA
   alignment evidence and valid replacements, but the replacements fail the
   current graph-size objective. This points to objective-selection or
   replacement induction quality as the next improvement area.

The next productive experiment is not more flubble detection. It is to improve
SweepGA/seqwish induction or post-induction cleanup for the directly generated
large residual candidates, and then rerun the same direct-residual windows under
the existing diagnostic logs.
