# Coverage-driven repeat-aware C4 crush

Task: `coverage-driven-repeat`

This pass adds an explicit coverage objective for multi-bubble crushing and
uses the node path-depth metrics from `graph-report` / `add-node-path` as the
acceptance signal. The implementation is intentionally transparent: it does
not rewrite PAF or hide alignments behind wrapper filters. Repeat handling is
only a logged candidate-selection decision.

## Lineage

- Parent graph: `c4-diagnose-residual-two-best`, commit `ed594fa`, input
  `/home/erikg/impg/data/c4_diagnose_residual_two_20260529T141339Z/outward_highbp_r3/run.nosort.gfa`.
- Coverage metrics came from `add-node-path`, commit `b8e8be3`.
- This pass extends the previous `iterative-multi-level` resolver from
  `88aa106` / `0f38df0` rather than adding a separate graph rewrite path.

## Implementation

New resolver mode:

```bash
--method coverage-multi-bubble
```

Aliases include `coverage-driven`, `coverage-window`, `coverage-crush`,
`node-coverage`, `path-coverage`, `coverage-repeat`, and
`repeat-aware-coverage`.

New multi-level objective:

```bash
--objective coverage
```

The existing `iterative-multi-level` default remains `--objective size`.
`coverage-multi-bubble` forces the coverage objective and turns on
repeat-aware candidate boundaries by default.

Coverage acceptance is two-stage:

1. Local candidate score compares the source window to its SweepGA/seqwish
   replacement using bp-weighted path-step coverage and singleton bp.
2. The rewritten whole graph must also improve the coverage objective.

A candidate passes the coverage objective if either bp-weighted node coverage
improves or singleton bp decreases. Both local and global logs include:

- bp-weighted coverage before/after;
- scaled coverage delta;
- singleton bp before/after;
- segment and segment-bp deltas;
- raw PAF record counts;
- collinear, off-diagonal, and dust-like PAF record counts;
- repeat-like and low-complexity node diagnostics.

Repeat-aware boundary selection is explicit. Candidate generation computes
node visit counts and treats a high-frequency node as repeat-like when it is
seen in at least half the paths and is either tiny (`<=64 bp`) or
low-complexity. A generated window is rejected only when both root anchors are
repeat-like. The round log reports this as `repeat-boundary-skipped=N`.
Accepted window alignments are still passed to SweepGA/seqwish unchanged.

## Code Surface

- `src/resolution.rs`
  - `ResolutionMethod::CoverageMultiBubble`
  - `MultiLevelObjectiveMode::{Size,Coverage}`
  - coverage objective scoring and whole-graph gate
  - outward residual windows included by default for coverage combined/sliding
    modes
  - repeat-aware boundary diagnostics
  - PAF alignment-profile diagnostics
- `src/main.rs`
  - `--objective` / `--window-objective`
  - `--repeat-aware-boundaries` / `--repeat-aware`
  - gfa-engine parser support for the same parameters
- `examples/compare_gfa_paths.rs`
  - explicit real-run path spelling validator used below

## Real C4 Runs

Output root:

```text
/home/erikg/impg/data/c4_coverage_driven_repeat_20260529T1812Z
```

Reference input:

```text
/home/erikg/impg/data/c4_diagnose_residual_two_20260529T141339Z/outward_highbp_r3/run.nosort.gfa
```

### Run 1: coverage objective, outward residual windows

```bash
target/release/impg crush \
  --gfa "$input" \
  --output "$run/run.nosort.gfa" \
  --method iterative-multi-level \
  --objective coverage \
  --window-mode outward \
  --window-target-bp 42000 \
  --max-window-sites 12 \
  --candidate-limit 8 \
  --max-total-sequence 6000000 \
  --min-objective-delta 1 \
  --max-iterations 3 \
  --auto-spoa-max-traversal-len 0 \
  --auto-poasta-max-traversal-len 0 \
  --min-match-length off \
  --max-pair-alignments 0 \
  --max-replacement-paf-bytes 0 \
  --sweepga-no-filter true \
  --polish-rounds 0 \
  --threads 32 \
  --verbose 1
```

Result: two 12-site outward windows accepted in rounds 1 and 2; round 3
was rejected by the global coverage gate.

### Run 2: repeat-aware combined multi-bubble windows

```bash
target/release/impg crush \
  --gfa "$input" \
  --output "$run/run.nosort.gfa" \
  --method coverage-multi-bubble \
  --window-mode combined \
  --window-target-bp 42000 \
  --max-window-sites 12 \
  --candidate-limit 8 \
  --max-total-sequence 6000000 \
  --min-objective-delta 1 \
  --max-iterations 3 \
  --auto-spoa-max-traversal-len 0 \
  --auto-poasta-max-traversal-len 0 \
  --min-match-length off \
  --max-pair-alignments 0 \
  --max-replacement-paf-bytes 0 \
  --sweepga-no-filter true \
  --polish-rounds 0 \
  --threads 32 \
  --verbose 1
```

Result: the generator produced top-level, sibling-run, parent-descendant,
sliding, level, and outward residual windows. Repeat-aware boundary selection
skipped 157,249 repeat-boundary candidates. One locally strong replacement
was rejected by the whole-graph coverage gate, so the graph remained unchanged.

## Metrics

| Graph | Segments | Links | Paths | Path steps | Segment bp | bp-weighted coverage | p10/median/p90 | Singleton bp | Wall | Accepted |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| PGGB control | 13,288 | 16,240 | 465 | 5,538,879 | 234,524 | 454.919 | 1/465/929 | 2,890 | 13:38 | n/a |
| `c4-diagnose-residual-two-best` input | 16,477 | 23,288 | 465 | 5,945,457 | 384,319 | 264.451 | 2/398/745 | 30,925 | n/a | n/a |
| Coverage outward r3 | 15,308 | 21,666 | 465 | 5,622,365 | 386,713 | 262.814 | 2/417/745 | 30,908 | 4:24 | 2 |
| Repeat-aware combined r3 | 16,477 | 23,288 | 465 | 5,945,457 | 384,319 | 264.451 | 2/398/745 | 30,925 | 1:40 | 0 |

The best accepted variant is `coverage_objective_outward_r3`: it reduces
singleton bp by 17 and reduces segment count, but it regresses total stored bp
and bp-weighted coverage. It is therefore a useful diagnostic and prototype
rather than a replacement for the previous visual best.

PNG artifact:

```text
https://hypervolu.me/~erik/impg/c4-coverage-driven-repeat-best.png
```

Remote upload was confirmed with:

```text
-rw-r--r-- 1 erik erik 967K May 29 18:31 www/impg/c4-coverage-driven-repeat-best.png
```

## Runtime Breakdown

No upstream query was run in these variants. Both started from the existing
`c4-diagnose-residual-two-best` GFA, so query/load was not timed separately;
the crush command wall time includes GFA input parsing and resolver setup.

| Variant | Phase | Wall time | Notes |
|---|---|---:|---|
| Coverage outward r3 | Total crush command | 4:24.16 | `/usr/bin/time`, max RSS 4,329,252 KB |
| Coverage outward r3 | Round 1 candidate generation | 11.99s | 15,044 outward candidates, 8 considered; detail render 2.20s, POVU parse 1.42s, POVU decompose 197.27ms |
| Coverage outward r3 | Round 1 condenser/graph induction | 52.70s | 8 SweepGA/seqwish candidates built |
| Coverage outward r3 | Round 1 lacing/rewrite/validate | 2.92s | 1 candidate accepted |
| Coverage outward r3 | Round 1 total | 68.52s | Global objective delta `8192` |
| Coverage outward r3 | Round 2 candidate generation | 9.66s | 15,061 outward candidates, 8 considered; detail render 1.87s, POVU parse 1.12s |
| Coverage outward r3 | Round 2 condenser/graph induction | 71.29s | 8 SweepGA/seqwish candidates built |
| Coverage outward r3 | Round 2 lacing/rewrite/validate | 4.17s | 1 candidate accepted |
| Coverage outward r3 | Round 2 total | 85.93s | Global objective delta `35718` |
| Coverage outward r3 | Round 3 candidate generation | 11.11s | 15,059 outward candidates, 8 considered; detail render 2.04s, POVU parse 1.11s |
| Coverage outward r3 | Round 3 condenser/global trial | 95.86s | Derived from 106.97s round total minus 11.11s generation; no separate build timer emitted on global rejection |
| Coverage outward r3 | Round 3 total | 106.97s | 1 locally accepted candidate rejected globally |
| Coverage outward r3 | Graph report | 0:12.69 | Rerun with `/usr/bin/time`, max RSS 1,661,452 KB |
| Coverage outward r3 | Sort for PNG | 0:27.82 | `gfasort -p Ygs -t 32`, max RSS 309,884 KB |
| Coverage outward r3 | Render PNG | 0:02.15 | `gfalook -m -x 2200 -y 1200`, max RSS 212,992 KB |
| Repeat-aware combined r3 | Total crush command | 1:40.42 | `/usr/bin/time`, max RSS 4,249,700 KB |
| Repeat-aware combined r3 | Round 1 candidate generation | 29.94s | 23,346 generated after 157,249 repeat-boundary skips; detail render 2.14s, POVU parse 1.34s |
| Repeat-aware combined r3 | Round 1 condenser/global trial | 69.70s | Derived from 99.64s round total minus 29.94s generation; no accepted lacing |
| Repeat-aware combined r3 | Round 1 total | 99.64s | 1 locally accepted candidate rejected globally |
| Repeat-aware combined r3 | Graph report | 0:13.57 | Rerun with `/usr/bin/time`, max RSS 1,766,128 KB |

## Candidate Diagnostics

### Coverage outward r3

| Round | Window | Raw PAF | Collinear | Off-diagonal | Dust-like | Repeat-like nodes | Local coverage | Singleton bp | Global result |
|---:|---|---:|---:|---:|---:|---|---|---|---|
| 1 | `>6919..>7045`, 12 sites, 31,345 bp span, 233 traversals | 117,159 | 51,083 | 66,076 | 66,076 | 802 nodes / 2,485 bp | 169.983 -> 162.441 | 309 -> 281 | accepted, delta `8192` |
| 2 | `>6917..>7130`, 12 sites, 31,814 bp span, 233 traversals | 122,475 | 51,083 | 71,392 | 71,392 | 182 nodes / 722 bp | 163.408 -> 166.097 | 282 -> 272 | accepted, delta `35718` |
| 3 | `>6978..>7384`, 12 sites, 32,007 bp span, 357 traversals | 238,631 | 80,071 | 158,560 | 158,560 | 297 nodes / 1,368 bp | 164.355 -> 168.863 | 276 -> 282 | rejected, delta `-19462481` |

The accepted round-1 candidate was singleton-driven, not coverage-driven: it
lost local bp-weighted coverage but reduced singleton bp enough to pass both
local and global gates. Round 2 improved both local coverage and singleton bp.
Round 3 looked locally attractive on coverage, but increased singleton bp
inside the candidate and expanded the whole graph enough that the global gate
rejected it.

### Repeat-aware combined r3

Round 1 summary:

```text
sources top-level=2034, sibling-run=9996, parent-descendants=43,
sliding-window=59, level-window=2036, outward-residual-window=9178,
repeat-boundary-skipped=157249
```

The strongest local candidate was:

| Window | Raw PAF | Collinear | Off-diagonal | Dust-like | Repeat-like nodes | Local coverage | Singleton bp | Global result |
|---|---:|---:|---:|---:|---|---|---|---|
| `>7226..>6886`, 1 site at level 3, 28,866 bp span, 197 traversals | 50,193 | 49,663 | 0 | 10 | 2,322 nodes / 13,996 bp | 81.541 -> 157.720 | 6,185 -> 202 | rejected, delta `-19826678` |

A short outward residual diagnostic window also showed the suspected repeat
signature:

| Window | Raw PAF | Collinear | Off-diagonal | Dust-like | Repeat-like nodes | Local score |
|---|---:|---:|---:|---:|---|---:|
| `>7262..>6916`, 882 bp root span, 364 traversals | 142,956 | 74,938 | 66,198 | 66,198 | 1,652 nodes / 12,217 bp | 5,728,811 |

That pattern is exactly what the task hypothesis predicted: high-frequency
short anchors plus many off-diagonal/dust-like PAF records in a residual
cluster. The global rejection is also important: the local replacement would
have lowered whole-graph bp-weighted coverage from 264.451 to 244.856 and
raised whole-graph singleton bp from 30,925 to 31,154, so it was correctly not
applied.

## Repeat and Coverage Findings

Residual separation is not only a window-generation problem. The combined
generator found broader adjacent, nested, and outward residual windows, and
the coverage logs show the candidate boundaries and alignment evidence.

The remaining failure appears to be repeat/dust-driven and condenser-limited:

- Highest-depth nodes in both the input and output are 1-bp sequence nodes
  visited thousands of times. Examples in the input include `17276` (`C`,
  total 17,417 visits across 465 paths), `17888` (`G`, total 16,112 visits
  across 459 paths), and `17885` (`T`, total 14,408 visits across 459 paths).
- Rare repeated local contexts remain high: 1,581 nodes / 2,198 occurrences
  in the input and 1,562 nodes / 2,176 occurrences after the best accepted
  coverage pass.
- The residual A/B windows around `>69xx..>73xx` repeatedly contain tens of
  thousands of off-diagonal PAF records. Round 3 had 158,560 off-diagonal and
  dust-like records out of 238,631 raw records.
- Some locally good replacements are rejected because they damage the global
  graph. This means the global objective is doing useful work, but the current
  SweepGA/seqwish condenser is still not PGGB-equivalent for this region.

Compared with PGGB, the gap is still large: PGGB has bp-weighted node coverage
454.919 and only 2,890 singleton bp. The best accepted coverage pass has
bp-weighted coverage 262.814 and 30,908 singleton bp.

## Path Preservation

Both real C4 outputs preserved all paths exactly against the input graph using
the committed `examples/compare_gfa_paths.rs` validator:

```text
coverage_objective_outward_r3:
expected_paths        465
observed_paths        465
missing_paths         0
extra_paths           0
spelling_mismatches   0

coverage_repeat_aware_combined_r3:
expected_paths        465
observed_paths        465
missing_paths         0
extra_paths           0
spelling_mismatches   0
```

The resolver also validates path spellings internally after every accepted
replacement, so a spelling-changing candidate cannot be applied.

## Validation

Focused tests run:

```text
source ./env.sh && cargo check --bin impg
source ./env.sh && cargo test coverage --lib -- --nocapture
source ./env.sh && cargo test repeat_aware --lib -- --nocapture
source ./env.sh && cargo test iterative_multi_level_generates_outward_residual_windows --lib -- --nocapture
source ./env.sh && cargo test test_gfa_output_format_accepts_coverage_multi_bubble_crush_params --bin impg -- --nocapture
source ./env.sh && cargo test graph_report::tests::reports_node_path_coverage_metrics_on_known_fixture --lib -- --nocapture
source ./env.sh && cargo test --test test_crush_integration c4_fragment -- --test-threads=1 --nocapture
source ./env.sh && cargo test --test test_crush_integration c4_top_flubble_seqwish_indexes_observed_exact_runs -- --test-threads=1 --nocapture
source ./env.sh && cargo build --release --bin impg
```

Full-suite attempt:

```text
source ./env.sh && cargo test --all -- --test-threads=1
```

This passed library tests, CLI tests, C4 integration tests, syng integration
tests, and then stalled in the pre-existing
`test_one_edge_rskip_side_survives_load` case from
`tests/test_syng_startcount.rs`. That case ran for roughly ten minutes at
100% CPU with no output before I interrupted it. The rest of that test binary
passes when that one case is skipped:

```text
source ./env.sh && cargo test --test test_syng_startcount -- --skip test_one_edge_rskip_side_survives_load --test-threads=1
```

`cargo fmt --check` still reports pre-existing formatting diffs under
`vendor/gfaffix`; the changed Rust files were formatted directly with
`rustfmt`.
