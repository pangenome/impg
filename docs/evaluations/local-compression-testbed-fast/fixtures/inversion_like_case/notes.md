# Orientation-sensitive inversion-like branch

- Fixture ID: `inversion_like_case`
- Class: `inversion_like`
- Tier: `ci`
- Assertion: `inversion_like_orientation_preserved`
- Known failure mode: Resolver reverse-complements the wrong interval or creates long links between inverted anchors.
- Render status: skipped for fast profile rows; reason `render_tool_not_configured_for_fast_profile`.
- Method status counts: pass=12

| Method | Command status | Exact paths | Topology | Graph | Render |
| --- | --- | --- | --- | --- | --- |
| `local_syng_raw` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/inversion_like_case/local_syng_raw/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_auto` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/inversion_like_case/local_syng_crush_auto/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_poa` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/inversion_like_case/local_syng_crush_poa/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_poasta` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/inversion_like_case/local_syng_crush_poasta/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `local_syng_crush_sweepga` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/inversion_like_case/local_syng_crush_sweepga/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `top_flubble_nonoverlap_sweepga` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/inversion_like_case/top_flubble_nonoverlap_sweepga/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `chunk_window_smooth_or_crush` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/inversion_like_case/chunk_window_smooth_or_crush/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `chunk_window_sweepga_seqwish` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/inversion_like_case/chunk_window_sweepga_seqwish/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `whole_region_sweepga_seqwish` | pass | pass | pass | `docs/evaluations/local-compression-testbed-fast/fixtures/inversion_like_case/whole_region_sweepga_seqwish/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `pggb_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/inversion_like_case/pggb_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `smoothxg_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/inversion_like_case/smoothxg_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
| `pggb_plus_smoothxg_control` | pass | pass | fail | `docs/evaluations/local-compression-testbed-fast/fixtures/inversion_like_case/pggb_plus_smoothxg_control/output.normalized.gfa` | skipped: render_tool_not_configured_for_fast_profile |
