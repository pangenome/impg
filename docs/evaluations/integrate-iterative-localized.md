# Evaluation: integrate-iterative-localized

Task: `integrate-iterative-localized`
Evaluator: `agent-674`
Date: 2026-06-08

## Grade

Overall score: **0.24 / 1.00**

Confidence: **0.86**

Rubric underspecified: **false**. The task gives clear required behavior, artifact policy, and validation criteria.

## Dimension Scores

| Dimension | Score | Rationale |
| --- | ---: | --- |
| Integrated runnable iterative pipeline | 0.08 | The branch has separate local seed, dirty detector, and localized resolver APIs, but no production caller that wires seed induction -> dirty-region detection -> localized resolution -> re-analysis iterations. `src/lib.rs:1192-1214` returns the induced `syng-local` seed GFA directly, and `src/local_seed.rs:613-638` writes dirty reports only as debug artifacts. |
| Localized SmoothXG-like behavior over SYNG local sequence sets | 0.25 | SYNG-collected local seed induction exists and preserves names, and localized resolution can polish detector-selected chunks. However, the task required a localized polishing loop over SYNG local sequence sets; the resolver API only consumes already-supplied chunks and is not reached from the `syng-local` path. |
| Required controls and conservative defaults | 0.28 | Component-level knobs exist for dirty detector chunk merge distance/flank/budget (`src/graph_badness.rs:15-52`) and resolver flank/budget/method thresholds (`src/resolution.rs:1102-1254`). Missing are CLI/API controls for the integrated loop, including max localized-polish iterations and runtime/budget controls across iterations. |
| Exact path hard failure and path-name preservation | 0.60 | The localized resolver has exact path spelling as a hard gate and tests path preservation (`tests/test_localized_resolver.rs:69-220`). This is a strong component-level implementation, but not proven through the full requested loop. |
| Metrics reporting as diagnostics | 0.55 | Dirty-region metrics are diagnostic-only and localized resolver reports before/after metrics without quality gating. The reports are not surfaced through an integrated command path. |
| CLI/API fit | 0.18 | The existing graph/query output syntax is extended for `syng-local` seed induction, which is directionally compatible. There is no CLI/API entry point for iterative localized polishing or resolver thresholds in the graph/query pipeline. |
| Validation coverage | 0.20 | Unit/integration tests cover dirty detection and localized resolver components. The required synthetic end-to-end tests for clean unchanged, underaligned improvement, dirty chunks reported, path spellings preserved through the loop, exact corruption stopping replacement, and diagnostic metrics through the integrated path are absent. A focused cargo test attempt did not reach tests because the environment lacks `htslib/faidx.h` for the vendored `wfmash-rs` build. |
| Artifact and git hygiene | 0.20 | No task-specific generated artifacts were staged by this evaluator, but the task has no task-specific implementation commit or push. The branch already contains many prior commits and fixtures from dependency work, so artifact compliance for the intended integration task cannot be demonstrated. |

## Key Evidence

The current `syng-local` path builds and returns a seed graph:

- `src/lib.rs:1192-1214` collects local sequences, calls `local_seed::induce_seed_graph`, and returns `seed.gfa`.

Dirty-region analysis is currently debug-report-only from seed induction:

- `src/local_seed.rs:613-638` calls `graph_badness::analyze_gfa` and writes JSON/TSV reports when a debug directory exists.

The localized resolver is only an API taking precomputed chunks:

- `src/resolution.rs:1113-1254` defines `resolve_gfa_dirty_chunks(gfa, chunks, config)` and applies chunk replacements, but there is no production caller outside tests.
- `rg "resolve_gfa_dirty_chunks\\(" src tests` finds only the API definition and `tests/test_localized_resolver.rs` calls.

Component tests are meaningful but not end-to-end:

- `tests/test_localized_resolver.rs:69-220` covers flanked replacement, path preservation, diagnostic-only quality metrics, budget skips, and failed SweepGA without fallback.
- `tests/test_graph_badness.rs` covers dirty detector behavior independently.

The requested integrated loop is absent:

- There is no code path that re-runs `graph_badness::analyze_gfa` after each localized replacement.
- There is no convergence/budget loop over dirty chunks.
- There is no CLI/API surface for integrated localized-polish iterations, chunk merge distance, flank length, resolver thresholds, or runtime budget as a single runnable path.
- There are no synthetic end-to-end tests for the integrated behavior.

## Validation Notes

Attempted focused validation:

```text
cargo test --test test_graph_badness --test test_localized_resolver --test test_graph_output_crush query_syng_local_seed_driver_preserves_collected_path_spellings
```

Result: build failed before tests ran because `wfmash-rs` could not find `htslib/faidx.h` while compiling the vendored C++ dependency. This blocks local cargo validation in the current environment and is not evidence that the Rust tests themselves fail.

## Rationale

This task was primarily an integration task. The prerequisite components are present and individually useful, especially the localized resolver's exact path guard and diagnostic reporting. However, the required deliverable was a runnable iterative localized polishing pipeline over SYNG-collected local sequence sets. The current branch stops at seed graph construction and optional dirty diagnostics; it does not consume dirty chunks, apply localized reinduction/polishing, lace replacements back through repeated iterations, or expose the required controls through the graph/query pipeline.

Because the core requested integration is missing, this should not be accepted as satisfying `integrate-iterative-localized`. The score gives partial credit for the prerequisite pieces already in place and for component-level tests, but the task-level outcome is incomplete.
