# Evaluation: `crush-smoothxg-on-2`

Task: `crush-smoothxg-on-syng-blunt`
Evaluator: `agent-188`
Date: 2026-05-27

## Rubric Status

Underspecification flag: **false**.

The task supplied explicit hard gates. It did not provide numeric weights, so this evaluation weights the outcome-critical real C4 run, metrics, publication, and required documentation more heavily than static parser support.

One status caveat matters: at evaluation time `wg show crush-smoothxg-on-2` reported the task as `in-progress`, with no task-specific completion logs or artifacts. I therefore grade the available branch state and workgraph record as the actor deliverable at grading time.

## Evidence Inspected

- `wg show crush-smoothxg-on-2`: no recorded task artifacts and no required validation logs from an implementation actor.
- `git log main..HEAD`: the branch contains prior project commits, with the visible head `1182ca7 feat: crush-smoothxg-on (agent-173)`, but no commit specific to `crush-smoothxg-on-2` or `crush-smoothxg-on-syng-blunt`.
- Workspace search: `docs/crush-smoothxg-on-syng-blunt.md` is absent, and there is no local reference to `c4-syng-blunt-smoothxg.png`.
- Static code inspection:
  - `src/main.rs` parses a generic `smooth` / `smoothxg` pipeline stage and does not require a preceding `crush` stage.
  - `src/lib.rs` applies `smooth_after_crush` after any optional crush step, so the current inherited implementation likely allows `gfa:syng:blunt:smooth` to run smoothing and then gfaffix.
  - Existing tests cover `gfa:syng:crush:smooth:nosort`, but I found no direct test for `gfa:syng:blunt:smooth` skipping crush.

## Hard Gate Assessment

| Gate | Result | Rationale |
| --- | ---: | --- |
| Pipeline accepts `gfa:syng:blunt:smooth` or equivalent, skipping crush | **partial** | The inherited generic `:smooth` parser and transform path appear to allow this shape. There is no task-specific commit or direct test proving the exact requested path. |
| `cargo test --all` passes | **fail / no evidence** | No actor log or artifact records this. I did not find a target-specific test run result. |
| Real C4 run completes, 465/465 paths preserved | **fail / absent** | No C4 output, stderr, summary, or path preservation report was recorded for this task. |
| Final metrics: segs, segment-bp, dup-extras, trivial-stringy, wall | **fail / absent** | The required metrics are not present in a task doc or artifact. |
| PNG uploaded as `c4-syng-blunt-smoothxg.png`, confirmed with `ssh ls` | **fail / absent** | No upload proof or matching filename was found. |
| `docs/crush-smoothxg-on-syng-blunt.md` committed | **fail** | The required document does not exist. |
| `wg artifact crush-smoothxg-on-syng-blunt docs/crush-smoothxg-on-syng-blunt.md` | **fail** | No matching artifact is visible in `wg show`; the referenced document is absent. |

## Dimension Scores

| Dimension | Score | Notes |
| --- | ---: | --- |
| Pipeline implementation and skip-crush routing | 0.40 | Some credit for inherited code that likely makes `syng:blunt:smooth` work, but the task did not add or verify the requested path directly. |
| Automated validation | 0.00 | No `cargo test --all` evidence and no direct parser regression for the requested pipeline. |
| Real C4 experiment | 0.00 | No completed run, no 465/465 path preservation report. |
| Metrics quality and comparison | 0.00 | Required metrics were not produced. |
| PNG publication proof | 0.00 | Required file and `ssh ls` confirmation absent. |
| Required documentation and workgraph artifact | 0.00 | Required task doc and artifact absent. |
| Workgraph hygiene | 0.00 | Target remained `in-progress` at evaluation time and lacked implementation progress logs. |

## Overall Grade

**Score: 0.10 / 1.00**

Confidence: **0.94**.

This is a near-failing result. The only meaningful credit is for pre-existing or inherited pipeline support that appears to make a direct `syng:blunt:smooth` route possible. The actor did not deliver the experiment that the task was primarily about: no real C4 run, no path preservation proof, no final metrics, no PNG upload confirmation, no required documentation, and no workgraph artifact. Because six of seven hard gates are absent and the remaining gate is only statically inferred, a score around 0.10 is calibrated.
