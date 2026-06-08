# Dispersed repeat glue should be broken or ignored

- Fixture ID: `dispersed_repeat_glue_break`
- Class: `dispersed_repeat_glue_break_or_ignore`
- Tier: `local`
- Assertion: `dispersed_repeat_no_long_glue`
- Known failure mode: seqwish/SYNG transitive closure glues distant repeat copies into long links.
- Render status: skipped for fast profile rows; reason `render_tool_not_configured_for_fast_profile`.
- Method status counts: skipped=12

| Method | Command status | Exact paths | Topology | Graph | Render |
| --- | --- | --- | --- | --- | --- |
| `local_syng_raw` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_auto` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_poa` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_poasta` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_sweepga` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `top_flubble_nonoverlap_sweepga` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `chunk_window_smooth_or_crush` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `chunk_window_sweepga_seqwish` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `whole_region_sweepga_seqwish` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `pggb_control` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `smoothxg_control` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
| `pggb_plus_smoothxg_control` | skipped | not_run | not_run | `profile_excludes_local_fixture` | skipped: render_tool_not_configured_for_fast_profile |
