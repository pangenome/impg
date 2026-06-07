# Local Compression Autopoietic Summary

Task: `autopoietic-local-compression-loop`

## Rolling Status

| Iteration | Date | Hypothesis | Change | Result | Next recommendation |
| --- | --- | --- | --- | --- | --- |
| 1 | 2026-06-07 | Bounded sorted chunk windows can avoid `nested_top_level_wrong` parent overmerge while preserving exact paths. | Implemented same-length chunk-window generation under `chunk_window_smooth_or_crush`; promoted/tightened `nested_top_level_wrong` to CI fast profile; added targeted runner test. | Supported on the synthetic adversarial fixture: chunk window passes with two candidates/two bubbles, one-parent placeholders fail visibly, exact path corruption remains 0. | Improve compression/quality metrics so future tuning compares topology, graph size, path depth, white-space, and internal-control behavior without over-rewarding synthetic structure. |

## Current Best Fast-Profile Signal

`chunk_window_smooth_or_crush` is the current fast-profile winner only for the promoted `nested_top_level_wrong` adversarial fixture. It is not a C4-ready winner:

- Fast profile after iteration 1: 143 rows; 121 command pass, 22 skipped; 121 exact-path pass, 0 hard path corruption; 71 topology pass, 50 topology fail, 22 topology not_run.
- `chunk_window_smooth_or_crush`: 11 included rows, 11 topology pass, exact paths pass on every graph-producing row, mean runtime 0.000364 seconds.
- One-parent local compact placeholders: 11 included rows, 10 topology pass and 1 topology fail after `nested_top_level_wrong` promotion.
- Internal controls: 11 included rows per control, exact paths pass, topology fail on all included fixtures under current synthetic assertions.

## Lineage Notes

- `chunk_window_smooth_or_crush` descends from the prior `compact_bubble` placeholder route. Iteration 1 adds a bounded same-length sorted-chunk branch with `min_invariant_gap_bp=5`, falling back to the parent behavior when no multi-window structure is detected.
- `nested_top_level_wrong` descends from the original local-tier adversarial fixture. Iteration 1 changes tier to `ci` and requires at least two bubbles/flubbles, matching its stated role as an overbroad-parent boundary check.

## Guardrails

- Exact path spelling remains the hard requirement.
- No hidden filtering or graph-quality acceptance gate has been added.
- Do not claim C4 readiness from these synthetic fixture results.
- If the next two algorithm iterations do not improve quality, spend the following iteration on metrics/fixtures instead of further tuning.
