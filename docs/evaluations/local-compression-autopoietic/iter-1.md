# Local Compression Autopoietic Iteration 1

Date: 2026-06-07

Task: `autopoietic-local-compression-loop`

## Hypothesis

Bounded SmoothXG-style sorted chunk windows should split same-length independent variant clusters in the `nested_top_level_wrong` adversarial fixture, while the existing highest-parent compact placeholder should remain visible as an overmerge topology failure. Exact path spelling is a hard requirement.

This was the highest-value first iteration because the current fast scoreboard had real internal controls, but it still could not distinguish the seven local compact/window method IDs on a fixture designed to punish top-level overmerge. Before this change, `nested_top_level_wrong` was local-tier and all 11 rows were skipped in the fast profile.

## Primitive Modification

Lineage:

- Parent method primitive: `chunk_window_smooth_or_crush`, previously routed through the same `compact_bubble` placeholder as the other local compact methods.
- Offspring method behavior: same method id, now uses a same-length sorted chunk-window generator when variant runs are separable by at least 5 invariant bases; otherwise it falls back to `compact_bubble`.
- Parent fixture primitive: `nested_top_level_wrong`, previously tier `local` with `bubble_count` and `flubble_count` minimum 1.
- Offspring fixture primitive: tier `ci`; `bubble_count` and `flubble_count` minimum 2, preserving the same exact path spellings and assertion id.

Implementation:

- Added chunk-window run discovery and graph generation to `scripts/local_compression_testbed.py`.
- Recorded method parameter lineage with `min_invariant_gap_bp=5`.
- Recorded method-specific `candidate_count` so the chunk row reports two candidate windows on the adversarial fixture.
- Added `local_compression_chunk_window_exposes_nested_parent_overmerge` to `tests/test_local_compression_testbed.rs`.
- Regenerated `tests/test_data/local_compression/manifest.json` and `nested_top_level_wrong` metadata/notes.

## Before/After Scoreboard

Fast profile totals:

| Metric | Before | After |
| --- | ---: | ---: |
| Rows | 143 | 143 |
| Command pass | 110 | 121 |
| Command skipped | 33 | 22 |
| Exact path pass | 110 | 121 |
| Exact path not_run | 33 | 22 |
| Hard path corruption | 0 | 0 |
| Topology pass | 70 | 71 |
| Topology fail | 40 | 50 |
| Topology not_run | 33 | 22 |

The increase in topology failures is intentional: the promoted adversarial fixture now exposes overbroad one-parent outputs instead of skipping them.

`nested_top_level_wrong` before: all 11 rows skipped with `profile_excludes_local_fixture`.

`nested_top_level_wrong` after:

| Method | Exact paths | Topology | Candidates | Bubbles | Flubbles | Runtime seconds | Result |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| `local_syng_raw` | pass | fail | 0 | 0 | 0 | 0.000487 | raw path-copy graph remains topology-incomplete |
| `local_syng_crush_auto` | pass | fail | 1 | 1 | 1 | 0.000449 | one-parent placeholder exposed |
| `local_syng_crush_poa` | pass | fail | 1 | 1 | 1 | 0.000435 | one-parent placeholder exposed |
| `local_syng_crush_poasta` | pass | fail | 1 | 1 | 1 | 0.000448 | one-parent placeholder exposed |
| `local_syng_crush_sweepga` | pass | fail | 1 | 1 | 1 | 0.000459 | one-parent placeholder exposed |
| `top_flubble_nonoverlap_sweepga` | pass | fail | 1 | 1 | 1 | 0.000352 | overbroad top-level boundary exposed |
| `chunk_window_smooth_or_crush` | pass | pass | 2 | 2 | 2 | 0.000384 | hypothesis supported on this fixture |
| `whole_region_sweepga_seqwish` | pass | fail | 1 | 1 | 1 | 0.000385 | one-parent placeholder exposed |
| `pggb_control` | pass | fail | 1 | 0 | 0 | 1.440512 | internal control preserves paths but does not induce the expected two-bubble topology |
| `smoothxg_control` | pass | fail | 1 | 0 | 0 | 0.032247 | internal control preserves paths but does not induce the expected two-bubble topology |
| `pggb_plus_smoothxg_control` | pass | fail | 1 | 0 | 0 | 1.410196 | internal control preserves paths but does not induce the expected two-bubble topology |

After-run method means on graph-producing fast rows:

| Method | Included rows | Topology pass/fail | Mean graph bytes | Mean segments | Mean path steps | Mean white-space proxy | Mean runtime seconds |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| `local_syng_raw` | 11 | 0/11 | 432.0 | 3.364 | 3.364 | 0.000 | 0.000456 |
| `local_syng_crush_auto` | 11 | 10/1 | 588.5 | 4.909 | 10.36 | 0.545 | 0.000409 |
| `top_flubble_nonoverlap_sweepga` | 11 | 10/1 | 588.5 | 4.909 | 10.36 | 0.545 | 0.000347 |
| `chunk_window_smooth_or_crush` | 11 | 11/0 | 670.1 | 5.364 | 12.36 | 2.545 | 0.000364 |
| `whole_region_sweepga_seqwish` | 11 | 10/1 | 588.5 | 4.909 | 10.36 | 0.545 | 0.000351 |
| `pggb_control` | 11 | 0/11 | 308.5 | 3.727 | 3.727 | 0.000 | 1.433 |
| `smoothxg_control` | 11 | 0/11 | 308.5 | 3.727 | 3.727 | 0.000 | 0.03354 |
| `pggb_plus_smoothxg_control` | 11 | 0/11 | 308.5 | 3.727 | 3.727 | 0.000 | 1.423 |

The chunk-window row improves topology discrimination and remains SYNG-fast on this synthetic fixture, but it is larger and has more path steps than the one-parent placeholder. It does not prove PGGB/SmoothXG-like compression quality, and the internal controls still fail these topology assertions despite preserving exact paths.

## Validation

Passed:

- `python3 -m py_compile scripts/local_compression_testbed.py`
- `python3 scripts/local_compression_testbed.py write-fixtures --root tests/test_data/local_compression`
- `cargo test --test test_local_compression_testbed local_compression_chunk_window_exposes_nested_parent_overmerge -- --nocapture`
- `python3 scripts/local_compression_testbed.py run --profile fast --manifest tests/test_data/local_compression/manifest.json --out-dir docs/evaluations/local-compression-testbed-fast`
- `cargo build`
- `cargo test`
- `cargo install --path .`

Cargo validation required initialized `vendor/syng` and `vendor/gfaffix` submodules plus local htslib/jemalloc include/library paths, recorded in `docs/evaluations/local-compression-testbed-fast/validation.md`.

## Result

Supported, with limits. The iteration produced a real distinction in the fast scoreboard:

- Exact path preservation stayed intact for every graph-producing row.
- `chunk_window_smooth_or_crush` is now the only method that passes `nested_top_level_wrong`.
- Bad outputs remain visible as rows; no graph-quality gate or hidden filter was introduced.

This is an improvement in quality and clarity rather than C4 readiness. It shows that bounded sorted chunking can avoid one specific parent-overmerge pattern in the synthetic testbed, not that the local method is generally PGGB/SmoothXG-like.

## Next Hypothesis

Improve the graph-quality metrics before another broad tuning pass. The next iteration should add a metric that separates "compact but path-copy-like" controls from "topology-shaped but over-expanded" local outputs, such as path replay compression ratio, per-path shared-bp reuse, or branch-window count consistency. This matters because internal controls currently preserve exact paths and produce small graphs but fail bubble/flubble topology assertions, while local placeholders can pass topology with synthetic structure.
