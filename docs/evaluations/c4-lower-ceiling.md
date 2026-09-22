# Evaluation: c4-lower-ceiling

Task: `c4-lower-ceiling`
Evaluator: `agent-474`
Date: 2026-06-04

## Rubric Note

Rubric underspecified: **true**

The task did not provide an explicit numerical grading rubric or weights. It did,
however, give concrete acceptance criteria: run bounded iterative-multi-level
C4 crush variants from the specified seed, keep sPOA for median traversals under
2 kb, test auto-POASTA ceilings at 2.5 kb and 5 kb for small/medium bubbles,
route complete-homologous 10-15 kb windows to SweepGA, use modest candidate
limits and 32 threads, validate path preservation as the only fatal gate, run
`graph-report`, sort with `gfasort -p Ygs`, render/upload PNGs, compare against
the seed plus the two named prior variants, and deliver a concise markdown report
with runtime, max RSS, candidate counts, graph metrics, and PNG URLs.

I graded against those deliverables as the implicit rubric.

## Grade

Overall score: **0.00 / 1.00**
Confidence: **0.97**

No actor-visible task deliverable was produced. At evaluation start, `wg show
c4-lower-ceiling` reported this task as `in-progress` on the evaluator branch,
with no task artifacts, no task-specific validation logs, and no completion log
from an actor. The only completed predecessor, `.assign-c4-lower-ceiling`, was an
assignment helper that marked itself done after spawning and did not record any
outputs.

Repository and data-directory searches found prior C4 data, including the
pre-existing comparison directory
`/home/erikg/impg/data/c4_tiered_poasta_spoa_20260604T1420Z` containing
`direct_auto_spoa2k_poasta15k` and
`multilevel_cap64_spoa2k_poasta15k_r1` artifacts. Those artifacts predate this
task and are only the requested comparison baselines. I found no
`c4-lower-ceiling` report, no lower-ceiling sweep output directory, no 2.5 kb or
5 kb auto-POASTA ceiling variants from the specified seed, and no uploaded PNG
URLs for such variants.

## Dimension Scores

| Dimension | Score | Rationale |
| --- | ---: | --- |
| Required sweep execution | 0.00 | No command log or output directory shows any bounded iterative-multi-level lower-ceiling crush variants run from `/home/erikg/impg/data/c4_grouped_multi_bubble_seed_top1_20260604T121002Z/sorted/combined_compact_12000bp_4sites_top1_r1.Ygs.gfa`. |
| Routing and parameter compliance | 0.00 | No evidence shows sPOA retained for median <2 kb, auto-POASTA ceilings tested at 2.5 kb and 5 kb, complete-homologous 10-15 kb windows routed to SweepGA, modest candidate limits, or 32-thread execution. |
| Path preservation and fatal-gate behavior | 0.00 | No task-specific path-preservation validation exists. There is no evidence that path corruption, and only path corruption, was treated as fatal. |
| Graph reporting, sorting, rendering, and upload | 0.00 | No lower-ceiling `graph-report.tsv`, `gfasort -p Ygs` output, `gfalook -m` PNG, upload log, or PNG URL exists for the requested variants. |
| Comparative markdown report | 0.00 | No concise report compares the seed, `direct_auto_spoa2k_poasta15k`, `multilevel_cap64_spoa2k_poasta15k_r1`, and new lower-ceiling variants with runtimes, max RSS, accepted/rejected/failed counts, graph metrics, and PNG URLs. |
| Provenance and WG hygiene | 0.00 | Before this evaluation file was added, no commit or WG artifact attributable to `c4-lower-ceiling` task execution was visible. The task record contained no actor completion or validation trail. |

## Evidence Checked

- `wg show c4-lower-ceiling` showed status `in-progress`, no artifacts, and no
  actor logs beyond task publication/spawn plus this evaluator's logs.
- `wg show .assign-c4-lower-ceiling` showed the assignment node as `done`, but it
  had only spawn/done logs and no artifacts.
- `wg context c4-lower-ceiling` reported no dependency artifacts.
- `find . -maxdepth 4` and `rg` searches for `lower-ceiling`, `lower ceiling`,
  `auto-poasta ceiling`, `direct_auto_spoa2k_poasta15k`, and
  `multilevel_cap64_spoa2k_poasta15k` found no task-scoped report or
  lower-ceiling output in the repository.
- Targeted searches under `/home/erikg/impg/data` found the existing
  `/home/erikg/impg/data/c4_tiered_poasta_spoa_20260604T1420Z` comparison
  artifacts for the two named prior variants, but no matching lower-ceiling
  variant directories or files.
- The specified seed directory
  `/home/erikg/impg/data/c4_grouped_multi_bubble_seed_top1_20260604T121002Z`
  exists and contains prior grouped-multi seed outputs, reports, validation
  files, sorted GFAs, and renders. Those are inputs/baseline context rather than
  evidence that this task's lower-ceiling sweep was run.
- `git status --short` showed only untracked WG metadata before this evaluation
  report was created.

## Acceptance Gate Assessment

| Required item | Result | Notes |
| --- | ---: | --- |
| Use the specified `combined_compact_12000bp_4sites_top1_r1.Ygs.gfa` seed | **fail** | No task command or output proves the seed was used. |
| Run bounded iterative-multi-level variants | **fail** | No lower-ceiling variant output exists. |
| Keep sPOA for median <2 kb | **fail** | No routing logs or command config exist. |
| Try auto-POASTA ceilings of 2.5 kb and 5 kb | **fail** | No matching variant names, command logs, or outputs exist. |
| Route complete-homologous 10-15 kb windows to SweepGA | **fail** | No SweepGA routing evidence exists for this task. |
| Keep candidate limits modest and use 32 threads | **fail** | No task run command exists. |
| Validate path preservation and avoid metric-based gates | **fail** | No task validation artifacts or acceptance logs exist. |
| Run graph-report, `gfasort -p Ygs`, and `gfalook -m` | **fail** | No lower-ceiling graph reports, sorted GFAs, or PNGs exist. |
| Upload PNGs to `erik@hypervolu.me:www/impg/` | **fail** | No lower-ceiling upload log or URLs exist. |
| Compare against seed and the two named prior variants | **fail** | The two prior variants exist as baseline artifacts, but no comparison report was delivered. |
| Deliver concise markdown report with runtime/RSS/counts/metrics/URLs | **fail** | No actor report exists. |

## Calibration

A score of **0.00** is warranted because none of the requested task outputs are
present. This is not a case of partial execution with missing uploads or
incomplete metrics; there is no visible lower-ceiling sweep attempt, no
validation evidence, and no markdown deliverable to assess.

The confidence is high because the task's requested artifact names and comparison
targets are specific, the prerequisite seed and prior comparison artifacts are
visible, and targeted searches found no new lower-ceiling artifacts. The only
residual uncertainty is whether an external-only run or upload occurred without
any local WG artifact, commit, command log, or report. That would still fail the
task as specified because the deliverable was a concise markdown report with
reproducible local evidence.
