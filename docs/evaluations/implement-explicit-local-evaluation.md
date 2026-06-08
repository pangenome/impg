# Evaluation: implement-explicit-local

Task: `implement-explicit-local`
Evaluator: `agent-660`
Date: 2026-06-08
Overall score: **0.24 / 1.00**
Confidence: **0.72**

## Rubric Status

The task rubric was sufficiently specified: it listed required behavior, validation commands, artifact policy, and commit/push expectations. There is, however, an evaluation-state mismatch: `wg evaluate run implement-explicit-local` could not run because the task record was still `InProgress`, not `done`, `failed`, or `pending-eval`. This manual evaluation therefore grades the submitted repository state and visible task evidence rather than a structured actor output bundle.

## Dimension Scores

| Dimension | Score | Rationale |
| --- | ---: | --- |
| Local sequence collection to seed graph induction path | 0.25 | Query-selected sequence extraction exists, and there is an anchor-seeded PAF to seqwish tail in `src/lib.rs`. The normal CLI `syng` path defaults to blunt SYNG GFA mode, so the anchor-seeded induction route is not cleanly exposed as the required local seed driver. |
| Current best SweepGA/FastGA + seqwish-style seed route | 0.20 | Existing graph and replacement code can run seqwish and SweepGA/FastGA-style replacement paths, but the submitted state does not make the whole-region SweepGA/FastGA + seqwish route the explicit initial local graph route for SYNG-collected regions. |
| Existing graph-output modes | 0.40 | Parser and CLI tests cover several `gfa:*:crush`, `syng`, `syng-local`, `seqwish`, and `pggb` combinations, but the relevant route is not proven end-to-end on a SYNG-collected local fixture. |
| Semantic path names and PanSN identifiers | 0.45 | Some path-name and path-sequence preservation tests exist, including PanSN-like names, and the sequence-to-FASTA helpers use `meta.path_name()`. Coverage is not focused on the selected seed induction route and does not prove the absence of synthetic local IDs throughout that route. |
| Runtime and reproducibility recording | 0.30 | The code logs timings and writes debug summaries for seqwish-tail builds when debug output is enabled. It does not record a concise command/config bundle for this explicit local driver, and no task validation logs show reproduced commands. |
| Diagnostic-only metrics and corruption checks | 0.45 | Existing seqwish-tail and crush code reports diagnostic counters and exact path preservation tests exist. Because there is no clean explicit local driver, I could not verify that graph-quality metrics remain purely diagnostic for the requested path. |
| Focused tests and validation | 0.15 | Tests cover adjacent parser/path-preservation behavior, but I did not find focused tests proving path names and path sequences on small SYNG-collected local sequence fixtures for the selected seed route. Full `cargo build` and targeted cargo tests were blocked by an external `wfmash-rs`/htslib header failure in this environment. |
| Artifact and commit hygiene | 0.15 | No task-specific implementation commit was visible before this evaluation. I did not identify a clean submitted implementation commit for `implement-explicit-local`; the branch mainly contains inherited prior work plus this evaluation artifact. |

## Evidence

- `src/lib.rs:722` collects query-selected sequences for `SyngLocal`, but then rebuilds a fresh sparse regional SYNG GFA via `write_syng_region_gfa_from_sequences_with_params`, not the required seed induction route.
- `src/lib.rs:1066` has a `GfaEngine::SyngNative` branch that would collect sequences, build an anchor-seeded PAF, and call `syng_graph::build_gfa_from_paf_and_sequences` (`src/lib.rs:1088-1144`), but parser defaults force `syng_gfa_mode = Some(Blunt)` for SYNG engines (`src/main.rs:3936-3937`). The normal `gfa:syng:crush` test confirms this mode is `Some(SyngGfaMode::Blunt)` (`src/main.rs:15054-15070`), so the anchor-seeded seed induction branch is not the clean exposed path.
- `src/graph.rs:1036-1077` already provides local interval sequence extraction to `graph --gfa-engine seqwish`, preserving `meta.path_name()` in FASTA headers. This is useful infrastructure but not specifically the SYNG-local explicit seed graph driver requested.
- `src/syng_graph.rs:955-1134` exposes a PAF plus sequence set to the shared seqwish induction tail and records debug summaries. This is a reusable component, but it is not tied into an externally clear local seed induction driver with task-specific tests.
- `tests/test_graph_output_crush.rs:38-125` proves a PAF-backed `gfa:poa:crush` query preserves source-coordinate path spellings. That is relevant but not the selected SYNG-collected seed route.
- `src/lib.rs:1913-1992` proves a `syng-local` blunt GFA round-trip before and after crush. That still starts from sparse SYNG topology, which the task explicitly wanted to avoid treating as final.

## Validation Notes

Attempted validation:

- `wg evaluate run implement-explicit-local` failed because the task was still `InProgress`.
- `git diff --check` and `git diff --cached --check` passed before adding this evaluation artifact.
- `cargo build` failed before crate validation due to missing system header `htslib/faidx.h` while building vendored `wfmash-rs`.
- Targeted cargo tests were consequently blocked by the same `wfmash-rs` build failure.

Key unmet task criteria:

- No clean, exposed local sequence collection to explicit seed graph induction driver was demonstrated for SYNG-collected local regions.
- The current best whole-region SweepGA/FastGA + seqwish-style route was not selected as the initial local seed graph route.
- Focused tests for path names and path sequences on small SYNG-collected local sequence fixtures were not added or demonstrated.
- Full `cargo build`, `cargo test`, targeted graph-output/local-compression tests, and push evidence were not available.

## Calibration

This is above zero because the repository contains substantial reusable pieces: sequence extraction, seqwish induction, SYNG region GFA output, `:crush` transforms, path-preservation helpers, and related tests. It remains a low score because the required behavior was to consolidate and expose those pieces as an explicit local seed graph induction driver, and the submitted state does not show that completed or validated.
