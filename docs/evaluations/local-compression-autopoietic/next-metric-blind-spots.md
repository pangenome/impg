# Next Local-Compression Metric Blind Spots

Date: 2026-06-07

Task: `analyze-next-local-2`

## Basis

This report follows iteration 1 of `autopoietic-local-compression-loop` and uses the current fast scoreboard only as evidence, not as a new benchmark run.

Primary artifacts:

- Iteration report: `docs/evaluations/local-compression-autopoietic/iter-1.md`
- Rolling summary: `docs/evaluations/local-compression-autopoietic/summary.md`
- Current fast scoreboard: `docs/evaluations/local-compression-testbed-fast/scoreboard.tsv`
- Scoreboard report: `docs/evaluations/local-compression-testbed-fast/report.md`
- Runner metric implementation, read-only reference: `scripts/local_compression_testbed.py`

Iteration 1 changed the scoreboard in the intended direction: `nested_top_level_wrong` moved from skipped to CI-fast coverage, exact paths still pass on graph-producing rows, and only `chunk_window_smooth_or_crush` passes that adversarial fixture. The evidence is recorded in `docs/evaluations/local-compression-autopoietic/iter-1.md:32-64` and the current row set is visible in `docs/evaluations/local-compression-testbed-fast/scoreboard.tsv:134-144`.

The blind spot is also explicit in the iteration conclusion: the current winner is not C4-ready, it is larger and deeper than the one-parent placeholder, and controls still fail topology despite exact path preservation (`docs/evaluations/local-compression-autopoietic/iter-1.md:66-79`). The summary repeats the same guardrail: `chunk_window_smooth_or_crush` is the current fast-profile winner only for the promoted synthetic fixture, while controls fail all included topology assertions (`docs/evaluations/local-compression-autopoietic/summary.md:13-18`).

## Current Scoreboard Signal

Fast-profile totals after iteration 1 are:

| Signal | Current value | Source |
| --- | ---: | --- |
| Rows | 143 | `docs/evaluations/local-compression-autopoietic/summary.md:15` |
| Command pass / skipped | 121 / 22 | `docs/evaluations/local-compression-autopoietic/summary.md:15` |
| Exact-path pass / hard corruption | 121 / 0 | `docs/evaluations/local-compression-autopoietic/summary.md:15` |
| Topology pass / fail / not_run | 71 / 50 / 22 | `docs/evaluations/local-compression-autopoietic/summary.md:15` |

The key current scoreboard rows are:

| Fixture and rows | Observed signal | Blind-spot implication |
| --- | --- | --- |
| `nested_top_level_wrong`, `scoreboard.tsv:134-144` | All 11 graph-producing rows preserve exact paths. The six one-parent local compact/window rows fail with 1 bubble/flubble, `chunk_window_smooth_or_crush` passes with 2 candidates, 2 bubbles, and 2 flubbles, and all three internal controls fail with 0 bubbles/flubbles. | Topology counts now catch the intended overmerge, but they still do not rank compression quality or control-like behavior. |
| `tandem_copy_loop_keep`, `scoreboard.tsv:90-100` | Raw and three internal controls fail because they emit no loop, while seven local compact/window rows pass with 1 self loop and 1 repeat loop. | Loop correctness is count-based and does not yet prove that the loop is the right motif, is traversed correctly by each path, or is not merely a synthetic shape that satisfies the fixture. |
| `dispersed_repeat_glue_break`, `scoreboard.tsv:101-111` | The only explicitly long-glue/dispersed-repeat fixture is still local-tier and skipped in the fast profile. | The fast profile currently has no positive or adversarial long-link row that exercises long-link detection. |

## Blind Spots

### White-Space Proxy

The scoreboard column is named `white_space_proxy_bp_*`, but the current implementation is not measuring base-pair white space. It computes `max(0, len(path_steps) - 3)` per path and then reports total, p95, p99, and max (`scripts/local_compression_testbed.py:1150-1183`).

That makes the proxy useful as a cheap path-step inflation smell, but not as a compression-quality metric:

- Of 121 graph-producing rows, 111 have `white_space_proxy_bp_p95=0`. The metric is mostly silent.
- The 10 nonzero rows all have topology pass, so the current topology gate is not using white-space to discriminate bad structure.
- Iteration 1's current winner has a higher mean white-space proxy than the one-parent placeholder family: `chunk_window_smooth_or_crush` has mean white-space proxy 2.545, while local compact placeholders have 0.545 and controls have 0.000 (`docs/evaluations/local-compression-autopoietic/iter-1.md:68-77`).
- On `nested_top_level_wrong`, the passing chunk row has `white_space_proxy_bp_total=8`, `p95=2`, and `max=2`, while one-parent rows have 0. The proxy therefore penalizes the method that fixed the fixture because it splits the path into more replay steps.

Conclusion: keep the proxy as a diagnostic, but do not let it stand in for compression or graph reuse. It is path-step count dressed as bp.

### Long-Link Detection

All graph-producing fast rows currently report `long_link_count=0` and `long_link_max_span_bp=0`. The runner reads these fields from parsed metrics and generated `# testbed_metrics` comments rather than deriving a genomic-distance span for every GFA link (`scripts/local_compression_testbed.py:1081-1105` and `scripts/local_compression_testbed.py:1184-1189`).

The fast profile therefore only proves that current generators and controls did not report long links. It does not prove that an output with a wrong distant glue edge would be detected. The fixture most directly named for this risk, `dispersed_repeat_glue_break`, is skipped by fast-profile policy in every method row (`docs/evaluations/local-compression-testbed-fast/scoreboard.tsv:101-111`).

Conclusion: long-link detection is still unexercised in the fast scoreboard. A future fixture should force one bad long glue edge and require the metric to detect it, but this is second priority behind adding a control-independent compression metric.

### Loop Correctness

The tandem fixture gives useful coverage, but the current pass/fail criterion is still narrow. Self loops are counted from same-source/same-target GFA links, while `repeat_loop_count` is carried through the testbed metric channel (`scripts/local_compression_testbed.py:1081-1105` and `scripts/local_compression_testbed.py:1184-1185`). Topology scoring then checks exact counts and allowed ranges (`scripts/local_compression_testbed.py:1215-1252`).

On `tandem_copy_loop_keep`, seven local compact/window rows pass with one self loop and one repeat loop, while raw and the three controls fail with zero loops (`docs/evaluations/local-compression-testbed-fast/scoreboard.tsv:90-100`). That is useful, but it leaves three correctness gaps:

- A row can satisfy the count target without proving the loop is traversed in the expected copy-number order.
- A row can satisfy the count target without proving the loop sequence is the expected tandem motif.
- A row can satisfy the count target without proving there is not an unrelated self-loop elsewhere in the graph.

Conclusion: loop count is a good regression guard, not a semantic loop proof. A later improvement should add path replay assertions for loop segment multiplicity and motif spelling.

### Topology Pass/Fail

The topology bit did its job in iteration 1: it exposed the one-parent overmerge on `nested_top_level_wrong` and made `chunk_window_smooth_or_crush` the only passing row for that fixture (`docs/evaluations/local-compression-autopoietic/iter-1.md:50-64`; `docs/evaluations/local-compression-testbed-fast/scoreboard.tsv:134-144`).

The problem is that binary pass/fail collapses distinct outcomes:

- Internal controls preserve exact paths and produce smaller graphs, but fail all included topology assertions under the synthetic bubble/flubble expectations (`docs/evaluations/local-compression-autopoietic/summary.md:15-18`).
- Local compact placeholders pass 10 of 11 included topology rows, even though iteration 1 showed they were not distinguishable until `nested_top_level_wrong` was promoted (`docs/evaluations/local-compression-autopoietic/iter-1.md:11` and `docs/evaluations/local-compression-autopoietic/iter-1.md:54-61`).
- The chunk-window row now passes 11 of 11 included rows, but it is also larger and has more path steps than the one-parent placeholder (`docs/evaluations/local-compression-autopoietic/iter-1.md:73-79`).

Conclusion: topology pass/fail is a necessary safety signal, but it is not a ranking objective. It can reward synthetic shape conformance while missing whether the graph is compact, shared, and control-like.

## Highest-Value Next Improvement

Add a comment-independent path replay compression metric to the scoreboard before another algorithm-tuning iteration.

Recommended metric:

```text
path_replay_compression_ratio = total_spelled_path_bp / total_segment_bp
```

Where `total_spelled_path_bp` is computed from the output GFA paths after exact path replay, and `total_segment_bp` is already in the scoreboard. Larger is better because it means more path sequence is represented by shared graph sequence rather than separate path-copy segments.

Why this is the highest-value improvement:

- It directly separates path-copy controls from structured sharing without trusting `# testbed_metrics` comments.
- It reframes the white-space issue: extra path steps are not automatically bad if they substantially improve sequence reuse.
- It provides a continuous ranking signal beside topology pass/fail, so future synthesis can distinguish "topology-shaped but over-expanded" from "compact but path-copy-like."
- It is cheap: the runner already parses segments and paths, already checks exact spellings, and already reports `total_segment_bp`.

Concrete expected effect on current rows:

| Fixture / row family | Existing evidence | Replay-compression interpretation |
| --- | --- | --- |
| `nested_top_level_wrong` chunk row | `scoreboard.tsv:140`: `total_segment_bp=37`, exact paths pass, topology pass | Total spelled path bp is 128, so the ratio would be about 3.46. This captures real sharing across the two independent event windows. |
| `nested_top_level_wrong` one-parent rows | `scoreboard.tsv:135-139` and `scoreboard.tsv:141`: `total_segment_bp=62`, exact paths pass, topology fail | Ratio would be about 2.06. Better than path-copy, but worse than chunk-window sharing and still topology-wrong. |
| `nested_top_level_wrong` raw/control rows | `scoreboard.tsv:134` and `scoreboard.tsv:142-144`: `total_segment_bp=128`, exact paths pass, topology fail | Ratio would be 1.00. These are exact path-copy-like outputs. |
| `tandem_copy_loop_keep` compact/window rows | `scoreboard.tsv:91-97`: `total_segment_bp=28`, exact paths pass, topology pass | Total spelled path bp is 108, so the ratio would be about 3.86 and would recognize useful loop reuse despite nonzero path-step white-space. |
| `tandem_copy_loop_keep` raw/control rows | `scoreboard.tsv:90` and `scoreboard.tsv:98-100`: `total_segment_bp=108`, exact paths pass, topology fail | Ratio would be 1.00, separating path-copy outputs from loop-shaped compact outputs. |

This should be introduced as a diagnostic metric first, not as a new hard gate. After one fast-profile rerun, add per-fixture expected ranges only where the metric has stable behavior.

## Secondary Follow-Ups

- Add one CI-fast long-link-positive fixture, or promote a minimal slice of `dispersed_repeat_glue_break`, so `long_link_count` and `long_link_max_span_bp` are exercised by at least one graph-producing row.
- Add tandem loop replay assertions: per-path copy count over the loop segment, loop motif spelling, and a check that any counted repeat loop is traversed by expected paths.
- Keep the current topology pass/fail signal, but report it beside replay compression, graph bytes, path steps, candidate count, and exact-path status rather than treating it as the single winner signal.

## Validation

- This report references the current scoreboard rows with file paths: `docs/evaluations/local-compression-testbed-fast/scoreboard.tsv:90-111` and `docs/evaluations/local-compression-testbed-fast/scoreboard.tsv:134-144`.
- This report references iteration artifacts with file paths: `docs/evaluations/local-compression-autopoietic/iter-1.md` and `docs/evaluations/local-compression-autopoietic/summary.md`.
- The highest-value next improvement is identified as `path_replay_compression_ratio`.
- No runner defaults, testbed code, fixtures, or generated scoreboards were modified by this task.
