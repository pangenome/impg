# Evaluation: add-pggb-frequency

Task: `add-pggb-frequency`
Evaluator: `agent-348`
Date: 2026-05-31

## Grade

Overall score: **0.08 / 1.00**
Confidence: **0.90**
Rubric underspecified: **false**

The task description is specific enough to grade without inventing a rubric: it
names the required CLI options, target-frequency bin computation, objective
terms, hard invalid constraints, bounded reruns, validation commands, and final
results document. The available branch state does not contain a task-specific
implementation commit or the requested deliverables. The current
`scripts/c4-crush-cmaes.py` still uses the pre-existing scalar graph-report
objective and the previous docs describe PGGB frequency distance only as a
posthoc ranking metric.

Status caveat: at evaluation time, `wg show add-pggb-frequency` reported this
task as `in-progress` and assigned to this evaluator worktree. I found no
implementation commit or artifact for `add-pggb-frequency`; this grade is based
on the current branch contents and WG record available to the evaluator.

## Dimension Scores

| Dimension | Score | Rationale |
| --- | ---: | --- |
| Target-frequency CLI surface | 0.00 | `--help` exposes `--target-range` and `--target-bed` for query selection, but no `--target-freq-bins`, `--target-report`, `--target-frequency-weight`, `--target-shape-weight`, or equivalent target-objective controls. |
| Trial bp-weighted path-count frequency bins | 0.00 | No code computes or records trial frequency bins from `graph-report` or by parsing output GFA. There is no `data/c4_cmaes_targetfreq_runs` directory and no target-frequency trial JSON evidence. |
| Direct CMA objective integration | 0.00 | `raw_objective` and `tell_objective` are still derived only from existing graph-report terms. Frequency TV distance, log segment-bp ratio, log bp-weighted-depth ratio, and optional whitespace-p99 ratio are absent from `objective_terms`. |
| Hard invalid constraints / no impg quality gate | 0.70 | The inherited wrapper still keeps path spelling and full-query path count checks as invalid constraints outside impg. I found no evidence that a new impg quality gate was added. This is useful preservation of a constraint, not task-specific progress. |
| Search-space settings needed to rediscover defaults | 0.05 | Existing dimensions include `max_rounds` only in `[1, 4]`, `min_traversal_len` only in `[1, 20000]`, and `polish_rounds` only in `[1, 4]`, so the requested `min-traversal-len=0`/off and until-done/default behavior is not represented. |
| Required target-frequency CMA reruns | 0.00 | No short direct-target objective rerun is present for auto or full-query syng-local auto. The older `data/c4_cmaes_runs` studies predate the requested objective and do not satisfy this criterion. |
| Required target-frequency results document | 0.00 | `docs/c4-cmaes-targetfreq-results.md` is missing, so there is no explicit statement whether the direct frequency objective improved the C4 graph. |
| Validation evidence | 0.20 | The existing script passes `python -m py_compile`, `--help`, and a tiny dry-run fixture smoke. Existing gzip artifacts under `data/c4_cmaes_runs/best_by_study` pass `gzip -t`. These validate the prior wrapper, not the requested target-frequency work. |

## Evidence

- `scripts/c4-crush-cmaes.py:430-455` lists the current CLI options. Target
  selection options are present for query mode, but target-objective inputs and
  weights are absent.
- `scripts/c4-crush-cmaes.py:741-744` calls `compute_objective` only after the
  standard graph report and path validation pass; there is no target-frequency
  data path before objective calculation.
- `scripts/c4-crush-cmaes.py:892-921` defines the full objective as
  `total_segment_bp`, singleton fraction, whitespace p99/max, long bridges,
  duplicate sequence fraction, repeat-context fraction, and coverage reward.
  The requested PGGB frequency TV and target-ratio terms are not present.
- `docs/c4-cmaes-results.md:67-82` explicitly describes PGGB shape distance as
  a ranking aid rather than the raw CMA objective. That matches the task's
  stated "current state" rather than the desired new behavior.
- `docs/c4-cmaes-results.md:60-63` notes that the old CMA search space did not
  include `max-iterations=until-done`, `min-traversal-len=0`, or
  `polish-rounds=until-done`, which is the opposite of the requested rerun
  settings.
- `docs/c4-cmaes-targetfreq-results.md` and `data/c4_cmaes_targetfreq_runs/`
  are absent.

## Validation Attempted

The checks below passed:

```bash
python -m py_compile scripts/c4-crush-cmaes.py
scripts/c4-crush-cmaes.py --help
scripts/c4-crush-cmaes.py --mode crush-only \
  --input-gfa tests/test_data/crush/small_insertion.gfa \
  --method-family sweepga \
  --max-trials 1 \
  --population-size 2 \
  --study-dir /tmp/c4-cmaes-eval-dry-run \
  --dry-run
find data/c4_cmaes_runs/best_by_study -name '*.gfa.gz' -type f -print0 | \
  xargs -0 -r gzip -t
```

These checks only show that the pre-existing wrapper and old gzip artifacts are
not broken. They do not provide target-frequency objective validation because
the target-frequency feature and rerun artifacts are missing.

## Calibration

The score is close to zero because the central requested behavior is absent:
PGGB copy-frequency distance is still not part of `raw_objective` or
`tell_objective`, no trial frequency bins are recorded, no direct-objective
CMA-ES rerun was performed, and the required target-frequency results document
does not exist. The small amount of credit is for inherited infrastructure that
already preserves path constraints, records objective terms in trial JSON, runs
graph-report externally, and has a syntactically valid wrapper. That inherited
work would be useful starting material, but it does not complete the assigned
task.
