# Local Compression Autopoietic Summary

Task: `autopoietic-local-compression-loop`

## Rolling Status

| Iteration | Date | Hypothesis | Change | Result | Next recommendation |
| --- | --- | --- | --- | --- | --- |
| 1 | 2026-06-07 | Bounded sorted chunk windows can avoid `nested_top_level_wrong` parent overmerge while preserving exact paths. | Implemented same-length chunk-window generation under `chunk_window_smooth_or_crush`; promoted/tightened `nested_top_level_wrong` to CI fast profile; added targeted runner test. | Supported on the synthetic adversarial fixture: chunk window passes with two candidates/two bubbles, one-parent placeholders fail visibly, exact path corruption remains 0. | Improve compression/quality metrics so future tuning compares topology, graph size, path depth, white-space, and internal-control behavior without over-rewarding synthetic structure. |
| 2 | 2026-06-08 | A metric-instrumented chunked SweepGA/seqwish candidate should preserve the chunk-window boundary win while exposing replay-compression evidence. | Added `chunk_window_sweepga_seqwish` using the sorted chunk-window selector with `min_invariant_gap_bp=5` and `seqwish_k=19`; added diagnostic `path_replay_compression_ratio`; reran fast profile with internal controls. | Neutral algorithmically, improved diagnostically: the new row ties `chunk_window_smooth_or_crush` on topology and replay ratio, exact path corruption remains 0, and controls remain visible compact path-copy references at replay ratio 1.0. | Use the replay-ratio column to target a resolver-distinct chunk-window experiment that improves compactness/path depth without losing the `nested_top_level_wrong` two-window signal. |

## Current Best Fast-Profile Signal

`chunk_window_smooth_or_crush` and `chunk_window_sweepga_seqwish` are current fast-profile topology winners only for the promoted `nested_top_level_wrong` adversarial fixture. They are not C4-ready winners:

- Fast profile after iteration 2: 156 rows; 132 command pass, 24 skipped; 132 exact-path pass, 0 hard path corruption; 82 topology pass, 50 topology fail, 24 topology not_run.
- `chunk_window_smooth_or_crush`: 11 included rows, 11 topology pass, exact paths pass on every graph-producing row, mean replay ratio 2.621205, mean runtime 0.000377 seconds.
- `chunk_window_sweepga_seqwish`: 11 included rows, 11 topology pass, exact paths pass on every graph-producing row, mean replay ratio 2.621205, mean runtime 0.000380 seconds.
- One-parent local compact placeholders: 11 included rows, 10 topology pass and 1 topology fail after `nested_top_level_wrong` promotion.
- Internal controls: 11 included rows per control, exact paths pass, topology fail on all included fixtures under current synthetic assertions, mean replay ratio 1.0.

## Lineage Notes

- `chunk_window_smooth_or_crush` descends from the prior `compact_bubble` placeholder route. Iteration 1 adds a bounded same-length sorted-chunk branch with `min_invariant_gap_bp=5`, falling back to the parent behavior when no multi-window structure is detected.
- `chunk_window_sweepga_seqwish` descends from the same sorted chunk-window boundary selector, with resolver lineage recorded as SweepGA/seqwish at `seqwish_k=19`. In iteration 2 it is a diagnostic comparison row, not evidence of a resolver-distinct generated graph.
- `nested_top_level_wrong` descends from the original local-tier adversarial fixture. Iteration 1 changes tier to `ci` and requires at least two bubbles/flubbles, matching its stated role as an overbroad-parent boundary check.

## Guardrails

- Exact path spelling remains the hard requirement.
- No hidden filtering or graph-quality acceptance gate has been added.
- `path_replay_compression_ratio` is diagnostic only.
- Do not claim C4 readiness from these synthetic fixture results.
- If the next two algorithm iterations do not improve quality, spend the following iteration on metrics/fixtures instead of further tuning.
