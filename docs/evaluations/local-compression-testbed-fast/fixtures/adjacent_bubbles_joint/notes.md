# Adjacent variants that should be joint

- Fixture ID: `adjacent_bubbles_joint`
- Class: `adjacent_bubbles_compress_together`
- Tier: `ci`
- Assertion: `adjacent_bubbles_joint_context`
- Known failure mode: Per-leaf flubble boundaries leave two stringy bubbles separated by artificial white-space.
- Render status: skipped for fast profile rows; reason `render_tool_not_configured_for_fast_profile`.
- Method status counts: pass=12

| Method | Command status | Exact paths | Topology | Graph | Render |
| --- | --- | --- | --- | --- | --- |
| `local_syng_raw` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/adjacent_bubbles_joint/local_syng_raw/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_auto` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/adjacent_bubbles_joint/local_syng_crush_auto/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_poa` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/adjacent_bubbles_joint/local_syng_crush_poa/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_poasta` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/adjacent_bubbles_joint/local_syng_crush_poasta/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_sweepga` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/adjacent_bubbles_joint/local_syng_crush_sweepga/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `top_flubble_nonoverlap_sweepga` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/adjacent_bubbles_joint/top_flubble_nonoverlap_sweepga/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `chunk_window_smooth_or_crush` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/adjacent_bubbles_joint/chunk_window_smooth_or_crush/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `chunk_window_sweepga_seqwish` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/adjacent_bubbles_joint/chunk_window_sweepga_seqwish/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `whole_region_sweepga_seqwish` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/adjacent_bubbles_joint/whole_region_sweepga_seqwish/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `pggb_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/adjacent_bubbles_joint/pggb_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `smoothxg_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/adjacent_bubbles_joint/smoothxg_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `pggb_plus_smoothxg_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/adjacent_bubbles_joint/pggb_plus_smoothxg_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
