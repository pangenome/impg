# Local Compression Autopoietic Iteration 2

Date: 2026-06-08

Task: `continue-local-compression`

## Hypothesis

The next focused change from `iter-2-synthesis.md` was to add one
metric-instrumented chunked SweepGA/seqwish candidate,
`chunk_window_sweepga_seqwish`, and add
`path_replay_compression_ratio` as a diagnostic scoreboard metric. The candidate
should inherit the same sorted chunk-window selector that fixed
`nested_top_level_wrong` in iteration 1, keep `min_invariant_gap_bp=5`, keep
`seqwish_k=19`, and preserve exact path spellings.

The graph-quality metrics in this iteration are diagnostic only. The runner
still records topology failures, compactness differences, and control failures
as visible rows; exact path corruption remains the only hard fail status.

## Implementation

- Added mandatory method `chunk_window_sweepga_seqwish` to
  `scripts/local_compression_testbed.py`.
- Reused the existing sorted chunk-window boundary generator from
  `chunk_window_smooth_or_crush`; no runner defaults, fixture tiers, control
  rows, or quality gates were changed.
- Set the method lineage to `window_policy=sorted_chunk`,
  `min_invariant_gap_bp=5`, `resolver=sweepga_seqwish`, and `seqwish_k=19`.
- Added `path_replay_compression_ratio =
  total_spelled_path_bp / total_segment_bp` to graph metrics, TSV rows, JSON
  rows, and metrics JSON. The numerator is computed from parsed GFA path replay
  by summing the lengths of traversed segment steps.
- Added focused runner tests for the replay-ratio metric and the
  `chunk_window_sweepga_seqwish` row on `nested_top_level_wrong`.

## Before/After Scoreboard

Fast profile totals:

| Metric | Before iter 2 | After iter 2 |
| --- | ---: | ---: |
| Rows | 143 | 156 |
| Methods | 11 | 12 |
| Command pass | 121 | 132 |
| Command skipped | 22 | 24 |
| Exact path pass | 121 | 132 |
| Exact path not_run | 22 | 24 |
| Hard path corruption | 0 | 0 |
| Topology pass | 71 | 82 |
| Topology fail | 50 | 50 |
| Topology not_run | 22 | 24 |

The extra 13 rows are the new method across all fixtures. The two new skipped
rows are the existing fast-profile skips for local-tier fixtures. There were no
new hidden rejections.

`nested_top_level_wrong` after iteration 2:

| Method | Exact paths | Topology | Candidates | Bubbles | Flubbles | Segment bp | Path steps | Replay ratio | Graph bytes | Runtime seconds |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `local_syng_crush_sweepga` | pass | fail | 1 | 1 | 1 | 62 | 12 | 2.064516 | 668 | 0.000458 |
| `top_flubble_nonoverlap_sweepga` | pass | fail | 1 | 1 | 1 | 62 | 12 | 2.064516 | 668 | 0.000360 |
| `chunk_window_smooth_or_crush` | pass | pass | 2 | 2 | 2 | 37 | 20 | 3.459459 | 931 | 0.000377 |
| `chunk_window_sweepga_seqwish` | pass | pass | 2 | 2 | 2 | 37 | 20 | 3.459459 | 931 | 0.000396 |
| `whole_region_sweepga_seqwish` | pass | fail | 1 | 1 | 1 | 62 | 12 | 2.064516 | 668 | 0.000362 |
| `pggb_control` | pass | fail | 1 | 0 | 0 | 128 | 4 | 1.000000 | 326 | 1.521119 |
| `smoothxg_control` | pass | fail | 1 | 0 | 0 | 128 | 4 | 1.000000 | 326 | 0.032914 |
| `pggb_plus_smoothxg_control` | pass | fail | 1 | 0 | 0 | 128 | 4 | 1.000000 | 326 | 1.517814 |

The new candidate preserves the iteration-1 boundary win: two candidate
windows, two bubbles, and two flubbles on the adversarial nested fixture. The
new replay ratio also separates the row families: path-copy-like controls are
1.0, one-parent local rows are 2.064516, and chunk-window rows are 3.459459.

After-run method means on graph-producing fast rows:

| Method | Rows | Topology pass/fail | Mean graph bytes | Mean segments | Mean path steps | Mean white-space proxy | Mean replay ratio | Mean runtime seconds |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `local_syng_crush_sweepga` | 11 | 10/1 | 588.5 | 4.909 | 10.364 | 0.545 | 2.272330 | 0.000407 |
| `top_flubble_nonoverlap_sweepga` | 11 | 10/1 | 588.5 | 4.909 | 10.364 | 0.545 | 2.272330 | 0.000431 |
| `chunk_window_smooth_or_crush` | 11 | 11/0 | 670.1 | 5.364 | 12.364 | 2.545 | 2.621205 | 0.000377 |
| `chunk_window_sweepga_seqwish` | 11 | 11/0 | 670.1 | 5.364 | 12.364 | 2.545 | 2.621205 | 0.000380 |
| `whole_region_sweepga_seqwish` | 11 | 10/1 | 588.5 | 4.909 | 10.364 | 0.545 | 2.272330 | 0.000365 |
| `pggb_control` | 11 | 0/11 | 308.5 | 3.727 | 3.727 | 0.000 | 1.000000 | 1.522641 |
| `smoothxg_control` | 11 | 0/11 | 308.5 | 3.727 | 3.727 | 0.000 | 1.000000 | 0.034063 |
| `pggb_plus_smoothxg_control` | 11 | 0/11 | 308.5 | 3.727 | 3.727 | 0.000 | 1.000000 | 1.495954 |

## Result

Neutral as an algorithmic optimization, improved as instrumentation.

The new `chunk_window_sweepga_seqwish` row ties the prior
`chunk_window_smooth_or_crush` best fast-profile topology signal, including the
`nested_top_level_wrong` two-window win and exact path preservation. It does not
improve graph compactness or path depth relative to the prior best: both chunk
rows have the same generated topology, graph bytes, segment count, path steps,
white-space proxy, and replay ratio in this synthetic runner.

The diagnostic metric is useful. It now makes the compactness tradeoff explicit:
the internal PGGB/SmoothXG controls are still smaller and path-copy-like under
this fixture set with replay ratio 1.0, while chunk-window rows show higher
sharing but remain larger and deeper than the controls.

## Validation

Passed:

- `python3 -m py_compile scripts/local_compression_testbed.py`
- `cargo test --test test_local_compression_testbed local_compression_path_replay_compression_ratio -- --nocapture`
- `cargo test --test test_local_compression_testbed local_compression_chunk_window_sweepga_seqwish_nested_top_level_wrong -- --nocapture`
- `python3 scripts/local_compression_testbed.py run --profile fast --manifest tests/test_data/local_compression/manifest.json --out-dir docs/evaluations/local-compression-testbed-fast`
- TSV assertion script from `iter-2-synthesis.md`: passed and printed the nested fixture comparison rows.
- `cargo build`
- `cargo test`
- `cargo install --path .`

Cargo validation used the documented native environment from
`docs/evaluations/local-compression-testbed-fast/validation.md` and required
initialized `vendor/syng` and `vendor/gfaffix` submodules.

## Next Action

Do not treat `chunk_window_sweepga_seqwish` as converged or C4-ready. The next
iteration should use the new replay-ratio column to design a small resolver or
boundary experiment that moves toward control-like graph size and path depth
without losing the `nested_top_level_wrong` two-window topology signal. A useful
next focused change is to make the chunk-window SweepGA/seqwish row actually
exercise a resolver-distinct graph construction instead of sharing the same
synthetic generator as `chunk_window_smooth_or_crush`.
