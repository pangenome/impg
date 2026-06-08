# Evaluation: c4-refine-one

Task: `c4-refine-one`
Evaluator: `agent-457`
Date: 2026-06-03

## Grade

Overall score: **0.00 / 1.00**
Confidence: **0.98**
Rubric underspecified: **false**

The task description is specific enough to grade. It requires two SPOA-only
`impg crush` refinements from the visually preferred
`one_many_scaffold0.initial.gfa` seed, with `max_median_traversal_len` set to
5k and 10k, using the installed `/home/erikg/.cargo/bin/impg` binary and no
POASTA. It also requires graph reports, exact 465-path validation and no path
spelling/name regressions versus the seed, sorted GFA output, `gfalook` PNGs,
upload of two named PNG files, a runtime/metric table, a recommendation versus
the seed and previous `one_many_scaffold0.spoa5`, timestamped data-directory
artifacts, commit/push, and WG completion.

No task-scoped implementation artifacts were produced before this evaluation.
The only completed predecessor, `.assign-c4-refine-one`, was an assignment node
that marked itself done after spawning. At evaluation start, `wg show
c4-refine-one` reported this task as `in-progress` under the evaluator branch,
with no recorded actor artifacts and no validation or completion logs for the
requested refinement workflow.

## Dimension Scores

| Dimension | Score | Rationale |
| --- | ---: | --- |
| Correct seed and command constraints | 0.00 | No run command, command log, or output directory exists showing use of the requested seed, `/home/erikg/.cargo/bin/impg`, `--method poa`, `--max-span 0`, the requested POA scoring, or the 5k/10k median thresholds. |
| Required SPOA refinements | 0.00 | Neither the 5k nor the 10k SPOA-only crush/refinement output exists. Searches found no `data/c4_one_many_initial_spoa_large_*` directory and no matching median5k/median10k files. |
| Graph reports and path validation | 0.00 | No `graph-report.tsv`, path-validation TSV/JSON, or equivalent validation evidence was produced for either output. The required 465-path and no-regression checks were not run. |
| Sorted GFAs, PNG rendering, and upload | 0.00 | No sorted GFA, `gfalook` PNG, gfasort/gfalook logs, or upload log exists. There is no evidence that `c4-one-many-initial-spoa-median5k.png` or `c4-one-many-initial-spoa-median10k.png` was uploaded. |
| Runtime/metric table and recommendation | 0.00 | No report compares 5k, 10k, the seed, and previous `one_many_scaffold0.spoa5`; no runtime table or recommendation was delivered. |
| Git and WG hygiene | 0.00 | Before this evaluation report, the branch contained inherited C4 commits and only an untracked `.wg` entry, with no commit attributable to `c4-refine-one` artifacts. The WG task had not been marked done by an actor and had no task-specific validation logs. |

## Evidence

- `wg show c4-refine-one` reported status `in-progress`, assigned to
  `agent-457`, with only task-pause/publish/coordinator-spawn logs and this
  evaluator's startup log.
- `wg show .assign-c4-refine-one` reported the assignment helper as `done`, but
  it had no artifacts and only logged spawn and task-marked-done events.
- `wg context c4-refine-one` reported: "No artifacts available from
  dependencies yet."
- `find data -maxdepth 2 -type d -name 'c4_one_many_initial_spoa_large_*'` and
  the same search rooted at `/home/erikg/impg/data` returned no output
  directories.
- `find /home/erikg/impg/data -maxdepth 2 -type f` for matching
  `*one_many*spoa*`, `*one-many*spoa*`, `*median5k*`, and `*median10k*` files
  returned no relevant artifacts.
- `wg list | rg 'c4-refine-one|flip-c4-refine|assign-c4-refine'` showed only
  `c4-refine-one` as active among these names; the downstream FLIP task remains
  blocked on it.
- `git status --short` showed only untracked `.wg` before this evaluation
  report was added.
- `git log --oneline main..HEAD` showed inherited prior C4 work but no
  task-specific `c4-refine-one` refinement-artifact commit.

## Acceptance Gate Assessment

| Required item | Result | Notes |
| --- | ---: | --- |
| Use visually preferred initial one-many seed, not minmatch1 | **fail** | No command or output proves any seed was used. |
| Run SPOA-only crush/refinement at 5k and 10k median traversal thresholds | **fail** | No runs or output GFAs exist. |
| Use installed `/home/erikg/.cargo/bin/impg`; do not rebuild locally; no POASTA | **fail** | No task command log exists to verify binary or method choices. |
| Run graph-report TSV for each output | **fail** | No output-specific graph reports exist. |
| Validate 465 paths and no path spelling/name regressions versus seed | **fail** | No path validation artifacts exist. |
| Sort with `gfasort -p Ygs` and render with `gfalook -m` | **fail** | No sorted GFAs, PNGs, or tool logs exist. |
| Upload two specifically named PNGs to `erik@hypervolu.me:www/impg/` | **fail** | No upload artifacts or logs exist. |
| Produce runtime/metric table and recommendation versus seed and previous `one_many_scaffold0.spoa5` | **fail** | No comparative report was produced. |
| Write outputs under timestamped `data/c4_one_many_initial_spoa_large_*` directory | **fail** | No such directory exists. |
| Commit/push artifacts and mark done | **fail** | No artifact commit or actor completion is visible. |

## Calibration

A score of **0.00** is warranted because the actor-visible task execution did
not produce any of the required files, validations, metrics, uploads, or
completion records. This is not a partial-credit case where one threshold
completed, validation failed, upload was blocked, or a recommendation was
missing; the requested refinement workflow was not attempted in the available
task record.

The confidence is high because the task's required artifact names and directory
prefix are concrete and direct repository/worktree searches found no matching
outputs. The only ambiguity is procedural: the WG prompt scheduled an evaluator
on the same task id rather than exposing a separate completed actor task. That
does not materially affect the grade, because the visible state contains no
deliverable to assess.
