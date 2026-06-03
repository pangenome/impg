# C4 Blocker 04: Stop Repeated Broad Residual Churn

Task: `c4-blocker-04-stop`

Branch tested: `wg/agent-425/c4-blocker-04-stop`

Result: the repeated broad residual churn is stopped for the full eight-round
C4 command. The run preserves all 465 paths exactly. Graph-report metrics remain
diagnostics only and were not used as candidate-application gates.

## Implementation

`src/resolution.rs` now carries residual-novelty diagnostics and selection state
for iterative multi-level frontier selection. The guard is keyed by:

- semantic/root-coordinate interval from the root path range
- source ancestry from POVU site id, parent id, level, and leaf flag
- traversal count and traversal length profile
- root span and unique-step footprint

For broad top-level residuals, near-identical candidates are deferred before the
candidate cap if they overlap the retained root interval by at least 95 percent,
have similar root span and traversal structure, and are not materially expanded.
Structurally changed or expanded intervals are retained. The guard is retained
across intervening non-broad acceptances so a deferred broad residual cannot
return one round later after unrelated small work.

The change does not duplicate SweepGA filtering logic. Replacement building and
alignment evidence still come from the actual selected builder and SweepGA
semantics, and the logs report the raw evidence that was run.

## Targeted Tests

The regression tests added in `src/resolution.rs` are:

- `residual_novelty_defers_immediate_near_identical_broad_top_level`
- `residual_novelty_guard_survives_intervening_small_acceptance`
- `residual_novelty_retains_same_interval_when_structure_changes`

The second test reproduces the failure observed during the first full-run
attempt: round 5 deferred the broad site, accepted a small top-level site, and
round 6 was able to select the same broad site again until the guard was made
persistent across intervening small acceptances.

## Validation

Build and tests:

```bash
source ./env.sh
cargo build --bin impg
cargo test --lib resolution::
cargo install --path .
```

Results:

- `cargo build --bin impg`: passed
- `cargo test --lib resolution::`: passed, 95/95 tests
- `cargo install --path .`: passed after applying the existing local
  `wfmash-rs` `OUT_DIR` workaround for the vendored `rkmh.cpp` build issue

The validation worktree had the same local build prerequisites documented by
earlier C4 tasks:

- `vendor/syng` and `vendor/gfaffix` submodules initialized in the WG worktree
- existing `/home/erikg/.cargo/bin/wfmash` copied into the relevant Cargo
  `wfmash-rs` `OUT_DIR` so the build script did not rebuild vendored wfmash

## Full C4 Run

Input GFA:

`/home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/run.nosort.gfa`

Output directory:

`/home/erikg/impg/data/c4_blocker_04_stop_20260603T045904Z`

Command:

```bash
/usr/bin/time -v -o /home/erikg/impg/data/c4_blocker_04_stop_20260603T045904Z/time.txt \
  /home/erikg/.cargo/bin/impg crush \
  --gfa /home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/run.nosort.gfa \
  --output /home/erikg/impg/data/c4_blocker_04_stop_20260603T045904Z/full_c4_after_residual_novelty.gfa \
  --method iterative-multi-level \
  --window-mode largest \
  --window-target-bp 30k \
  --max-window-sites 8 \
  --candidate-limit 1 \
  --max-iterations 8 \
  --polish-method poasta \
  --polish-rounds 0 \
  --max-pair-alignments 0 \
  --max-replacement-paf-bytes 0 \
  --threads 32 \
  -v 2
```

Completion:

- Exit status: `0`
- Wall time: `1:52.81`
- Maximum resident set size: `2245372` KB
- Accepted candidates: `8/8`
- Per-run summary: `crush: 8 resolved, 0 bailed, 8 candidates seen across 8 rounds`

## Path Validation

External validation against the pre-crush C4 input:

```text
expected_paths	465
observed_paths	465
missing_paths	0
extra_paths	0
spelling_mismatches	0
```

Visible GFA path names remain semantic FASTA/PAF/GFA names. No sidecar name
recovery was used.

## Churn Guard Evidence

Blocker-03 repeatedly selected overlapping broad top-level intervals in rounds
5-8. In this run, the same broad root interval was discovered but deferred, and
the selected top-level candidate changed to smaller non-overlapping work:

| round | selected method | selected interval | selected root span | deferred broad interval(s) |
| --- | --- | --- | ---: | --- |
| 5 | Poasta | `>272214386..>272214802` | 7737 bp | `>272218972..>272218984`, `>272218972..>272218981` |
| 6 | Poasta | `>272214384..>272214803` | 7765 bp | `>272218972..>272218984`, `>272218972..>272218981` |
| 7 | Poasta | `>272218195..>272218464` | 6588 bp | `>272218972..>272218984`, `>272218972..>272218981` |
| 8 | Poasta | `>272235574..>272235721` | 7760 bp | `>272218972..>272218984`, `>272218972..>272218981` |

Representative round-8 log evidence:

```text
round 8 residual-novelty detail: #1 decision=defer reason=near-identical-immediate-broad-residual ... candidate=[round=8 source=top-level sites=1 ancestry=site=>272218972>272218984 ... root_steps=4702..5333 traversals=388 median=26752 total=7301649 root_span=33120 unique_steps=1506] previous=[round=4 source=top-level sites=1 ancestry=site=>272218973>272218980 ... root_steps=5102..5708 traversals=269 median=26631 total=7243526 root_span=32999 unique_steps=902]
round 8 generated candidate detail: #1 source=top-level sites=1 ancestry=site=>272235574>272235721 ... >272235574..>272235721 root_steps=511..658 root_span=7760bp traversals=465 max=7810 median=7761 total=3561277 unique_steps=224
```

## Artifacts

- GFA:
  `/home/erikg/impg/data/c4_blocker_04_stop_20260603T045904Z/full_c4_after_residual_novelty.gfa`
- Crush stderr log:
  `/home/erikg/impg/data/c4_blocker_04_stop_20260603T045904Z/crush.stderr.log`
- Crush stdout log:
  `/home/erikg/impg/data/c4_blocker_04_stop_20260603T045904Z/crush.stdout.log`
- Time report:
  `/home/erikg/impg/data/c4_blocker_04_stop_20260603T045904Z/time.txt`
- Graph-report TSV:
  `/home/erikg/impg/data/c4_blocker_04_stop_20260603T045904Z/full_c4_after_residual_novelty.graph-report.tsv`
- Path-validation TSV:
  `/home/erikg/impg/data/c4_blocker_04_stop_20260603T045904Z/full_c4_after_residual_novelty.path-validation.tsv`
- Scoreboard TSV:
  `/home/erikg/impg/data/c4_blocker_04_stop_20260603T045904Z/c4_blocker_04_scoreboard.tsv`
- Metric deltas TSV:
  `/home/erikg/impg/data/c4_blocker_04_stop_20260603T045904Z/c4_blocker_04_metric_deltas.tsv`

## Diagnostic Scoreboard

These metrics are diagnostic only. They were not candidate-application gates.

| label | segments | total segment bp | singleton bp | bp-weighted coverage | path whitespace p99 | path whitespace max |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| PGGB control, 465 paths | 13288 | 234524 | 2890 | 454.919215 | 14 | 219917 |
| Integrated C4 baseline | 19026 | 462940 | 78435 | 230.460695 | 160890 | 459430 |
| Blocker-03 full rerun | 20674 | 532484 | 101312 | 200.361840 | 196939 | 529093 |
| Blocker-04 residual novelty | 20256 | 537191 | 104683 | 198.606220 | 202069 | 533829 |

Current deltas:

| metric | vs integrated baseline | vs blocker-03 full rerun | current / PGGB |
| --- | ---: | ---: | ---: |
| segments | +6.46% | -2.02% | 1.524x |
| total segment bp | +16.04% | +0.88% | 2.291x |
| singleton bp | +33.46% | +3.33% | 36.222x |
| bp-weighted coverage | -13.82% | -0.88% | 0.437x |
| path whitespace p99 | +25.59% | +2.60% | 14433.500x |
| path whitespace max | +16.19% | +0.90% | 2.427x |

Interpretation: blocker-04 fixes the repeated broad residual frontier churn and
preserves exact paths, but it does not make the final graph PGGB-quality. The
next work should treat graph-shape metrics as diagnostics and investigate why
the alternative smaller replacements do not improve the final full-C4 graph.
