# Short insertion/deletion bubble

- Fixture ID: `short_indel_3path`
- Class: `short_indel`
- Tier: `ci`
- Assertion: `short_indel_bubble`
- Known failure mode: Resolver clips inserted bases, creates a dangling tip, or collapses deletion into the wrong path.
- Render status: skipped for fast profile rows; reason `render_tool_not_configured_for_fast_profile`.
- Method status counts: pass=12

| Method | Command status | Exact paths | Topology | Graph | Render |
| --- | --- | --- | --- | --- | --- |
| `local_syng_raw` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/short_indel_3path/local_syng_raw/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_auto` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/short_indel_3path/local_syng_crush_auto/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_poa` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/short_indel_3path/local_syng_crush_poa/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_poasta` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/short_indel_3path/local_syng_crush_poasta/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_sweepga` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/short_indel_3path/local_syng_crush_sweepga/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `top_flubble_nonoverlap_sweepga` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/short_indel_3path/top_flubble_nonoverlap_sweepga/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `chunk_window_smooth_or_crush` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/short_indel_3path/chunk_window_smooth_or_crush/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `chunk_window_sweepga_seqwish` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/short_indel_3path/chunk_window_sweepga_seqwish/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `whole_region_sweepga_seqwish` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/short_indel_3path/whole_region_sweepga_seqwish/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `pggb_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/short_indel_3path/pggb_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `smoothxg_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/short_indel_3path/smoothxg_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `pggb_plus_smoothxg_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/short_indel_3path/pggb_plus_smoothxg_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
