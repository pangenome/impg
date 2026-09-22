# Evaluation: implement-local-compression-fixtures

Task: `implement-local-compression-fixtures`
Evaluator: `agent-588`
Date: 2026-06-07

## Grade

Overall score: **0.00 / 1.00**
Confidence: **0.98**
Rubric underspecified: **false**

The task is explicit enough to grade. It requires a discoverable synthetic
fixture suite under a stable repository path, all thirteen required fixture
classes, per-fixture sequence/GFA inputs, machine-readable metadata with exact
path spellings and topology/compression expectations, CI versus local tiering,
and a validation script or cargo test that recovers expected path sequences
before the method matrix runs.

At evaluation time, the branch contained no task-specific implementation beyond
the inherited design document. There was no `tests/test_data/local_compression/`
directory, no manifest for local compression fixtures, no fixture metadata, no
new validation script or cargo test, no artifacts, and no implementation commit
for this task. Existing `tests/test_data/crush/` files predate this assignment
and do not satisfy the new fixture schema or required fixture classes.

## Dimension Scores

| Dimension | Score | Rationale |
| --- | ---: | --- |
| Fixture coverage | 0.00 | None of the thirteen required local-compression fixture classes were added to a manifest or stable fixture path. |
| Fixture inputs | 0.00 | No per-fixture sequence/GFA inputs were added for this task. |
| Metadata completeness | 0.00 | No metadata files exist with fixture ID, class, input paths, path names, exact expected path spellings, topology/compression assertion, metric ranges, known failure mode, tier, or render hints. |
| Path naming and semantics | 0.00 | No new fixture path names were submitted, PanSN-compatible or otherwise. |
| Manifest/discovery | 0.00 | No `manifest.json` or equivalent index exists for fixture discovery. |
| Exact path validation | 0.00 | No validation script or cargo test was added to load fixtures and verify expected path spellings are recoverable. |
| Topology assertions | 0.00 | No simple-fixture topology assertions or messy-fixture accepted ranges were added. |
| CI/local tiering | 0.00 | No CI-scaled fixture subset or local-only fixture metadata was provided. |
| Reuse of existing helpers | 0.00 | No implementation exists, so there is no helper reuse to credit. |
| Task hygiene and validation evidence | 0.00 | No task-specific commit, push, artifact, command log, or successful validation log was produced before evaluation. |

## Evidence

- `wg show implement-local-compression-fixtures` reported this task as
  `in-progress` for `agent-588`, with no actor progress logs before the
  evaluator start and no listed artifacts.
- `git show --stat HEAD --` showed branch head
  `7b7e0d1 feat: design-local-graph-compression-testbed (agent-585)`, changing
  only `docs/local-graph-compression-testbed-design.md`. That is the dependency
  design artifact, not an implementation of fixtures.
- `git status --short` showed only an untracked `.wg` entry before this
  evaluation artifact was written; there were no task-specific fixture files.
- `rg --files tests/test_data | rg 'local[_-]compression|compression|fixture'`
  returned no files.
- `ls tests/test_data` listed only pre-existing generic and `crush` test data:
  `a.fa`, `b.fa`, `c.fa`, `ref.fa`, `ref2.fa`, `test.agc`,
  `yeast.chrV.fa.gz`, and `crush/`.
- `find tests -maxdepth 3 -type f ...` found only older GFA-related tests and
  crush fixtures, such as `tests/test_data/crush/small_insertion.gfa` and
  `tests/test_data/crush/nested_bubbles_real.gfa`; it found no
  local-compression manifest, metadata, or validation fixtures.
- Repository search for schema terms such as `fixture_id`,
  `expected_path_spellings`, and `local_compression` found the design document
  and quality-pass document, but no implemented fixture metadata or test suite.

## Acceptance Gate Assessment

| Required validation item | Result | Notes |
| --- | ---: | --- |
| All required fixture classes present in a manifest/index | **fail** | No local-compression fixture manifest exists. |
| Readable stable path names | **fail** | No new fixture paths or path names were added. |
| Every fixture has required metadata | **fail** | No fixture metadata files exist. |
| Sequence/GFA inputs for every fixture | **fail** | No fixture input directories or files exist. |
| Exact expected path spellings checked before method matrix | **fail** | No script or cargo test was added. |
| Simple fixtures include explicit topology assertions | **fail** | No topology assertion metadata exists. |
| Messy fixtures include accepted metric ranges where relevant | **fail** | No messy/local-only fixture metadata exists. |
| CI subset is small and deterministic | **fail** | No CI fixture subset exists. |
| Heavier fixtures marked local-only and excluded from default CI/smoke | **fail** | No local-only fixture metadata exists. |
| Relevant cargo tests or fixture validation scripts pass and commands logged | **fail** | No validation command for this fixture set exists and no task log records such a run. |

## Calibration

A score of **0.00** is appropriate because the submitted state contains no
task-scoped deliverable for the requested fixture implementation. This is not a
case of partial fixture coverage or incomplete metadata; none of the required
fixture suite, discovery mechanism, or validation machinery was present.
Inherited design work and older crush fixtures are useful project context, but
they are not credited as performance on this implementation task.
