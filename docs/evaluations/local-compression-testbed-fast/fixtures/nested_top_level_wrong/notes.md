# Nested bubbles where top-level boundary over-merges

- Fixture ID: `nested_top_level_wrong`
- Class: `nested_bubbles_top_level_wrong`
- Tier: `local`
- Assertion: `nested_parent_boundary_wrong`
- Known failure mode: Top-level selection over-merges independent events and creates bad loops, long links, or excess path depth.
- Render status: skipped for fast profile rows; reason `render_tool_not_configured_for_fast_profile`.
- Method status counts: skipped=11

| Method | Command status | Exact paths | Topology | Graph | Render |
| --- | --- | --- | --- | --- | --- |
| `local_syng_raw` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_auto` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_poa` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_poasta` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_sweepga` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `top_flubble_nonoverlap_sweepga` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `chunk_window_smooth_or_crush` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `whole_region_sweepga_seqwish` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `pggb_control` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `smoothxg_control` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `pggb_plus_smoothxg_control` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
