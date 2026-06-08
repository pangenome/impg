# Duplicated flank requiring path context

- Fixture ID: `duplicated_flank_context`
- Class: `duplicated_flank_requires_path_context`
- Tier: `ci`
- Assertion: `duplicated_flank_context_required`
- Known failure mode: Replacement laces into the wrong duplicate flank or merges distinct occurrences.
- Render status: skipped for fast profile rows; reason `render_tool_not_configured_for_fast_profile`.
- Method status counts: pass=12

| Method | Command status | Exact paths | Topology | Graph | Render |
| --- | --- | --- | --- | --- | --- |
| `local_syng_raw` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/duplicated_flank_context/local_syng_raw/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_auto` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/duplicated_flank_context/local_syng_crush_auto/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_poa` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/duplicated_flank_context/local_syng_crush_poa/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_poasta` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/duplicated_flank_context/local_syng_crush_poasta/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_sweepga` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/duplicated_flank_context/local_syng_crush_sweepga/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `top_flubble_nonoverlap_sweepga` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/duplicated_flank_context/top_flubble_nonoverlap_sweepga/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `chunk_window_smooth_or_crush` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/duplicated_flank_context/chunk_window_smooth_or_crush/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `chunk_window_sweepga_seqwish` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/duplicated_flank_context/chunk_window_sweepga_seqwish/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `whole_region_sweepga_seqwish` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/duplicated_flank_context/whole_region_sweepga_seqwish/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `pggb_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/duplicated_flank_context/pggb_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `smoothxg_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/duplicated_flank_context/smoothxg_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `pggb_plus_smoothxg_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/duplicated_flank_context/pggb_plus_smoothxg_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
