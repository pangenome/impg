# C4 Blocker 05b: Complete Homologous Traversal Aggregation

Task: `c4-blocker-05b-complete-traversal-aggregation`

Branch tested: `wg/agent-432/c4-blocker-05b-complete-traversal-aggregation`

Result: complete homologous traversal aggregation fixes the concrete C4 failure
mode where partial traversal sets were induced and laced while the homologous
old structure remained in non-candidate paths. The full C4 run preserves all
465 paths exactly and improves the final graph diagnostics substantially versus
the integrated baseline, blocker-03, and blocker-04. It is still not PGGB-sized.

## Implementation

`src/resolution.rs` now emits a `complete-homologous-window` candidate source
for iterative multi-level candidate generation.

The aggregation is structural rather than metric-gated:

- broad incomplete top-level residuals are treated as seeds
- nearby POVU sites are clustered by root-coordinate overlap and bounded root
  span
- path ranges are unioned across the cluster
- bounded homologous node expansion is used only when range-only union does not
  cover all paths
- candidates are emitted only when the aggregate reaches every graph path

Important code locations:

- candidate source and counters: `MultiLevelCandidateSource` and
  `MultiLevelGeneratedCandidates` at `src/resolution.rs:2613`
- generator hook before parent/sibling candidates are capped:
  `add_complete_homologous_windows` call at `src/resolution.rs:3239`
- largest-mode priority that prevents a partial broad residual from outranking
  the complete aggregate: `complete_homologous_window_priority` at
  `src/resolution.rs:3407`
- aggregation helpers: `add_complete_homologous_windows` and
  `candidate_from_complete_homologous_cluster` at `src/resolution.rs:4101`
- path coverage and length logging: `format_candidate_path_coverage` at
  `src/resolution.rs:5214`

The change does not duplicate SweepGA filtering. Alignment evidence, PAF counts,
transclosure evidence, replacement segment counts, and replacement bp are logged
from the actual SweepGA/wfmash/seqwish replacement builder.

## Targeted Tests

The regression tests added in `src/resolution.rs` are:

- `partial_homologous_lacing_can_improve_locally_but_grow_globally`
- `iterative_multi_level_aggregates_complete_homologous_traversal_set_before_cap`
- `largest_mode_prioritizes_complete_homologous_before_partial_top_level_scale`

The first test is the deterministic reproducer for the blocker-04 diagnosis: an
incomplete candidate improves local replacement size, but applying it leaves old
homologous structure in uncovered paths and makes the global graph grow. The
second test demonstrates the new aggregation behavior and would fail without
the complete-window aggregation. The third test covers the selection hazard
where a larger partial top-level residual could otherwise win largest mode ahead
of a complete aggregate.

The fixture also checks visible replacement input headers stay semantic and do
not use local synthetic names such as `__impg` or `candidate_`.

## Validation

Build and tests:

```bash
cargo test --lib homologous
cargo test --lib resolution::
cargo build --bin impg
cargo install --path .
```

Results:

- `cargo test --lib homologous`: passed, 3/3 tests
- `cargo test --lib resolution::`: passed, 98/98 tests
- `cargo build --bin impg`: passed
- `cargo install --path .`: passed

The native build required the local htslib/wfmash library include and library
paths already used by earlier C4 work.

## Full C4 Run

Input GFA:

`/home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/run.nosort.gfa`

Output directory:

`/home/erikg/impg/data/c4_blocker_05b_complete_traversal_wfmash_final_20260603T070646Z`

Command:

```bash
LD_LIBRARY_PATH=/home/erikg/wfmash/build/vendored_htslib/lib:/home/erikg/micromamba/lib:/home/erikg/htslib-local/lib:$LD_LIBRARY_PATH \
/usr/bin/time -v \
  -o /home/erikg/impg/data/c4_blocker_05b_complete_traversal_wfmash_final_20260603T070646Z/time.txt \
  /home/erikg/.cargo/bin/impg crush \
    --gfa /home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/run.nosort.gfa \
    --output /home/erikg/impg/data/c4_blocker_05b_complete_traversal_wfmash_final_20260603T070646Z/full_c4_after_complete_traversal.gfa \
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
    --sweepga-aligner wfmash \
    --threads 32 \
    -v 2
```

Completion:

- Exit status: `0`
- Wall time: `1:19:58`
- Maximum resident set size: `5901676` KB
- Accepted candidates: `8/8`
- Per-run summary: `crush: 8 resolved, 0 bailed, 8 candidates seen across 8 rounds`

The default FastGA replacement builder was also probed. It selected the same
complete 465-path aggregate in the target C4 region, but the final validation
artifact is the wfmash-backed run above because the earlier FastGA attempt hit
a local FAtoGDB buffer failure. Both paths use actual SweepGA/seqwish semantics.

## Path Validation

External validation against the pre-crush C4 input:

```text
expected_paths	465
observed_paths	465
missing_paths	0
extra_paths	0
spelling_mismatches	0
```

No sidecar path-name recovery was used.

## Aggregation Evidence

The old blocker-04 partial residual neighborhood was still discovered in round
1 as incomplete fragments:

```text
site=>272218634>272219318 ... root_span=32834bp traversals=13
site=>272218544>272219303 ... root_span=32809bp traversals=2
site=>272218959>272219339 ... root_span=32797bp traversals=34
site=>272218582>272219308 ... root_span=32793bp traversals=40
site=>272218963>272219340 ... root_span=32787bp traversals=117
site=>272218968>272219341 ... root_span=32776bp traversals=58
```

The selected replacement was the complete aggregate spanning the same region:

```text
round 1 generated candidate detail: #1 source=complete-homologous-window
  sites=441 >272217628..>272219165 root_span=37498bp
  traversals=465 path_coverage=465/465
  min=24343 median=63231 p90=69959 max=117858 total=29491050
```

The builder evidence for that aggregate was:

```text
sequences=465
raw_paf_records=374025
filtered_paf_records=221722
transclosure: Using 32 threads
replacement_segments=2546
replacement_shared_segments=2199
replacement_bp=89584
input_segments=5805, input_bp=142891
output_segments=2546, output_bp=89584
```

All eight accepted C4 candidates were `complete-homologous-window` candidates
with `path_coverage=465/465`; no selected candidate in the final run was one of
the old partial `117/143/256/269/372/388` traversal fragments.

Selected candidate coverage by round:

| round | interval | sites | path coverage | root span bp | raw PAF | replacement segments | replacement bp |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | `>272217628..>272219165` | 441 | 465/465 | 37498 | 374025 | 2546 | 89584 |
| 2 | `>272219455..>272220964` | 410 | 465/465 | 37369 | 215760 | 1789 | 70051 |
| 3 | `>272234337..>272235792` | 249 | 465/465 | 38158 | 211386 | 1564 | 40384 |
| 4 | `>272233081..>272237497` | 577 | 465/465 | 37487 | 330849 | 2115 | 37164 |
| 5 | `>272214168..>272216077` | 420 | 465/465 | 37414 | 215760 | 1321 | 37397 |
| 6 | `>272240125..>272235901` | 401 | 465/465 | 37463 | 199879 | 1714 | 37707 |
| 7 | `>272233058..>272242402` | 565 | 465/465 | 38181 | 275611 | 1990 | 37691 |
| 8 | `>272245186..>272235905` | 485 | 465/465 | 37463 | 199892 | 1710 | 37845 |

## Diagnostic Scoreboard

These metrics are diagnostic only. They were not candidate-application gates.

| label | segments | total segment bp | singleton bp | bp-weighted coverage | path whitespace p99 | path whitespace max |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| PGGB | 13288 | 234524 | 2890 | 454.919215 | 14 | 219917 |
| Integrated baseline | 19026 | 462940 | 78435 | 230.460695 | 160890 | 459430 |
| Blocker 03 | 20674 | 532484 | 101312 | 200.361840 | 196939 | 529093 |
| Blocker 04 | 20256 | 537191 | 104683 | 198.606220 | 202069 | 533829 |
| 05b complete traversal | 14601 | 372408 | 29744 | 286.485451 | 148832 | 370891 |

Current deltas:

| metric | vs integrated baseline | vs blocker-03 | vs blocker-04 | current / PGGB |
| --- | ---: | ---: | ---: | ---: |
| segments | -23.26% | -29.38% | -27.92% | 1.099x |
| total segment bp | -19.56% | -30.06% | -30.67% | 1.588x |
| singleton bp | -62.08% | -70.64% | -71.59% | 10.292x |
| bp-weighted coverage | +24.31% | +42.98% | +44.25% | 0.630x |
| path whitespace p99 | -7.49% | -24.43% | -26.35% | 10630.857x |
| path whitespace max | -19.27% | -29.90% | -30.52% | 1.687x |

Interpretation: complete traversal aggregation fixes the concrete blocker-04
partial-candidate application mode and improves the full C4 graph substantially
over the integrated baseline and blockers 03/04. It does not yet match PGGB:
PGGB still has fewer segments, less total segment bp, much lower singleton bp,
and much lower whitespace p99.

## Artifacts

- GFA:
  `/home/erikg/impg/data/c4_blocker_05b_complete_traversal_wfmash_final_20260603T070646Z/full_c4_after_complete_traversal.gfa`
- Crush stderr log:
  `/home/erikg/impg/data/c4_blocker_05b_complete_traversal_wfmash_final_20260603T070646Z/crush.stderr.log`
- Crush stdout log:
  `/home/erikg/impg/data/c4_blocker_05b_complete_traversal_wfmash_final_20260603T070646Z/crush.stdout.log`
- Time report:
  `/home/erikg/impg/data/c4_blocker_05b_complete_traversal_wfmash_final_20260603T070646Z/time.txt`
- Graph-report TSV:
  `/home/erikg/impg/data/c4_blocker_05b_complete_traversal_wfmash_final_20260603T070646Z/full_c4_after_complete_traversal.graph-report.tsv`
- Path-validation TSV:
  `/home/erikg/impg/data/c4_blocker_05b_complete_traversal_wfmash_final_20260603T070646Z/full_c4_after_complete_traversal.path-validation.tsv`
- Scoreboard TSV:
  `/home/erikg/impg/data/c4_blocker_05b_complete_traversal_wfmash_final_20260603T070646Z/c4_blocker_05b_scoreboard.tsv`
- Metric deltas TSV:
  `/home/erikg/impg/data/c4_blocker_05b_complete_traversal_wfmash_final_20260603T070646Z/c4_blocker_05b_metric_deltas.tsv`
- PNG:
  `/home/erikg/impg/data/c4_blocker_05b_complete_traversal_wfmash_final_20260603T070646Z/c4_blocker_05b_scoreboard.png`
