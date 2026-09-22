# Nested bubbles where top-level boundary over-merges

- Fixture ID: `nested_top_level_wrong`
- Class: `nested_bubbles_top_level_wrong`
- Tier: `ci`
- Assertion: `nested_parent_boundary_wrong`
- Known failure mode: Top-level selection over-merges independent events and creates bad loops, long links, or excess path depth.
- Render status: skipped for fast profile rows; reason `render_tool_not_configured_for_fast_profile`.
- Method status counts: pass=12

| Method | Command status | Exact paths | Topology | Graph | Render |
| --- | --- | --- | --- | --- | --- |
| `local_syng_raw` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/nested_top_level_wrong/local_syng_raw/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_auto` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/nested_top_level_wrong/local_syng_crush_auto/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_poa` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/nested_top_level_wrong/local_syng_crush_poa/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_poasta` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/nested_top_level_wrong/local_syng_crush_poasta/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_sweepga` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/nested_top_level_wrong/local_syng_crush_sweepga/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `top_flubble_nonoverlap_sweepga` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/nested_top_level_wrong/top_flubble_nonoverlap_sweepga/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `chunk_window_smooth_or_crush` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/nested_top_level_wrong/chunk_window_smooth_or_crush/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `chunk_window_sweepga_seqwish` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/nested_top_level_wrong/chunk_window_sweepga_seqwish/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `whole_region_sweepga_seqwish` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/nested_top_level_wrong/whole_region_sweepga_seqwish/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `pggb_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/nested_top_level_wrong/pggb_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `smoothxg_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/nested_top_level_wrong/smoothxg_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `pggb_plus_smoothxg_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/nested_top_level_wrong/pggb_plus_smoothxg_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
