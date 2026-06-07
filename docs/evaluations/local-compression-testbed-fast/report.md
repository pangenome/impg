# Local Compression Testbed Fast Profile Artifact Index

Profile: `fast`

## Reproduction

```bash
python3 scripts/local_compression_testbed.py write-fixtures --root tests/test_data/local_compression
python3 scripts/local_compression_testbed.py run --profile fast --manifest tests/test_data/local_compression/manifest.json --out-dir docs/evaluations/local-compression-testbed-fast
```

## Source Inputs

- Manifest: `tests/test_data/local_compression/manifest.json`
- Fixture metadata/input directories: `tests/test_data/local_compression/<fixture-id>/`

## Scoreboards

- TSV: `docs/evaluations/local-compression-testbed-fast/scoreboard.tsv`
- JSON: `docs/evaluations/local-compression-testbed-fast/scoreboard.json`
- Fixture validation log: `docs/evaluations/local-compression-testbed-fast/fixture-validation.log`
- Fixture validation JSON: `docs/evaluations/local-compression-testbed-fast/fixture-validation.json`

## Coverage

- Fixtures represented: 13 (`adjacent_bubbles_joint`, `alu_like_insertion`, `dispersed_repeat_glue_break`, `duplicated_flank_context`, `fake_repeat_anchor_split`, `inversion_like_case`, `microtangle_repeat_motif`, `mid_insertion_200bp`, `nested_top_level_right`, `nested_top_level_wrong`, `short_indel_3path`, `snp_bubble_3path`, `tandem_copy_loop_keep`)
- Methods represented: 11 (`local_syng_raw`, `local_syng_crush_auto`, `local_syng_crush_poa`, `local_syng_crush_poasta`, `local_syng_crush_sweepga`, `top_flubble_nonoverlap_sweepga`, `chunk_window_smooth_or_crush`, `whole_region_sweepga_seqwish`, `pggb_control`, `smoothxg_control`, `pggb_plus_smoothxg_control`)
- Rows: 143
- Command statuses: pass=80, path_corrupt=30, skipped=33
- Topology statuses: fail=10, not_run=63, pass=70

## Control Execution

- Control timeout: 120s per row, one thread, no hidden topology gate.
- Exact path corruption is the only hard rejection; graph-quality and topology failures remain visible rows.
- Standalone binary discovery is recorded in command logs. This bounded profile uses the repository `impg graph --gfa-engine pggb` path when available; `smoothxg_control` uses that local smoothxg-style stage with an explicit empty PAF because this repo does not expose a standalone smooth-only CLI path.
- `pggb_control`: command_status path_corrupt=10, skipped=3; exact_paths fail=10, not_run=3; topology not_run=13; tool_available true=13.
  Skips: `profile_excludes_local_fixture`.
  Availability detail: profile_excludes_local_fixture;provider=local_impg_graph_pggb;standalone_pggb=not_found;standalone_smoothxg=not_found;FAtoGDB=/home/erikg/impg/.wg-worktrees/agent-603/target/debug/FAtoGDB;FastGA=/home/erikg/impg/.wg-worktrees/agent-603/target/debug/FastGA;GIXmake=/home/erikg/impg/.wg-worktrees/agent-603/target/debug/GIXmake;gfaffix=/home/erikg/impg/.wg-worktrees/agent-603/target/debug/gfaffix;impg=/home/erikg/impg/.wg-worktrees/agent-603/target/debug/impg;time=/usr/bin/time;timeout=/usr/bin/timeout
- `smoothxg_control`: command_status path_corrupt=10, skipped=3; exact_paths fail=10, not_run=3; topology not_run=13; tool_available true=13.
  Skips: `profile_excludes_local_fixture`.
  Availability detail: profile_excludes_local_fixture;provider=local_impg_smoothxg_stage_via_pggb_engine;standalone_pggb=not_found;standalone_smoothxg=not_found;smoothxg_path=local impg smoothxg-style stage exposed through `impg graph --gfa-engine pggb` with an explicit empty PAF;gfaffix=/home/erikg/impg/.wg-worktrees/agent-603/target/debug/gfaffix;impg=/home/erikg/impg/.wg-worktrees/agent-603/target/debug/impg;time=/usr/bin/time;timeout=/usr/bin/timeout
- `pggb_plus_smoothxg_control`: command_status path_corrupt=10, skipped=3; exact_paths fail=10, not_run=3; topology not_run=13; tool_available true=13.
  Skips: `profile_excludes_local_fixture`.
  Availability detail: profile_excludes_local_fixture;provider=local_impg_graph_pggb;standalone_pggb=not_found;standalone_smoothxg=not_found;FAtoGDB=/home/erikg/impg/.wg-worktrees/agent-603/target/debug/FAtoGDB;FastGA=/home/erikg/impg/.wg-worktrees/agent-603/target/debug/FastGA;GIXmake=/home/erikg/impg/.wg-worktrees/agent-603/target/debug/GIXmake;gfaffix=/home/erikg/impg/.wg-worktrees/agent-603/target/debug/gfaffix;impg=/home/erikg/impg/.wg-worktrees/agent-603/target/debug/impg;time=/usr/bin/time;timeout=/usr/bin/timeout

### Control Comparison

- Raw SYNG graph-producing rows: 10; exact paths pass=10; topology fail=10.
- Compact/window graph-producing rows: 70; exact paths pass=70; topology pass=70.
- PGGB/SmoothXG control graph-producing rows: 30; exact paths fail=30; topology not_run=30.


## Output Roots

- Per-fixture notes, command logs, graph outputs, normalized graphs, and metrics JSON:
  `docs/evaluations/local-compression-testbed-fast/fixtures/<fixture-id>/<method-id>/`
- Renders are explicitly skipped in the fast profile with row-level `render_skip_reason` values.

## Fixture Notes

- `adjacent_bubbles_joint`: `docs/evaluations/local-compression-testbed-fast/fixtures/adjacent_bubbles_joint/notes.md`
- `alu_like_insertion`: `docs/evaluations/local-compression-testbed-fast/fixtures/alu_like_insertion/notes.md`
- `dispersed_repeat_glue_break`: `docs/evaluations/local-compression-testbed-fast/fixtures/dispersed_repeat_glue_break/notes.md`
- `duplicated_flank_context`: `docs/evaluations/local-compression-testbed-fast/fixtures/duplicated_flank_context/notes.md`
- `fake_repeat_anchor_split`: `docs/evaluations/local-compression-testbed-fast/fixtures/fake_repeat_anchor_split/notes.md`
- `inversion_like_case`: `docs/evaluations/local-compression-testbed-fast/fixtures/inversion_like_case/notes.md`
- `microtangle_repeat_motif`: `docs/evaluations/local-compression-testbed-fast/fixtures/microtangle_repeat_motif/notes.md`
- `mid_insertion_200bp`: `docs/evaluations/local-compression-testbed-fast/fixtures/mid_insertion_200bp/notes.md`
- `nested_top_level_right`: `docs/evaluations/local-compression-testbed-fast/fixtures/nested_top_level_right/notes.md`
- `nested_top_level_wrong`: `docs/evaluations/local-compression-testbed-fast/fixtures/nested_top_level_wrong/notes.md`
- `short_indel_3path`: `docs/evaluations/local-compression-testbed-fast/fixtures/short_indel_3path/notes.md`
- `snp_bubble_3path`: `docs/evaluations/local-compression-testbed-fast/fixtures/snp_bubble_3path/notes.md`
- `tandem_copy_loop_keep`: `docs/evaluations/local-compression-testbed-fast/fixtures/tandem_copy_loop_keep/notes.md`

## Method Matrix

| Fixture | `local_syng_raw` | `local_syng_crush_auto` | `local_syng_crush_poa` | `local_syng_crush_poasta` | `local_syng_crush_sweepga` | `top_flubble_nonoverlap_sweepga` | `chunk_window_smooth_or_crush` | `whole_region_sweepga_seqwish` | `pggb_control` | `smoothxg_control` | `pggb_plus_smoothxg_control` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `adjacent_bubbles_joint` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | path_corrupt | path_corrupt | path_corrupt |
| `alu_like_insertion` | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture |
| `dispersed_repeat_glue_break` | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture |
| `duplicated_flank_context` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | path_corrupt | path_corrupt | path_corrupt |
| `fake_repeat_anchor_split` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | path_corrupt | path_corrupt | path_corrupt |
| `inversion_like_case` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | path_corrupt | path_corrupt | path_corrupt |
| `microtangle_repeat_motif` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | path_corrupt | path_corrupt | path_corrupt |
| `mid_insertion_200bp` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | path_corrupt | path_corrupt | path_corrupt |
| `nested_top_level_right` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | path_corrupt | path_corrupt | path_corrupt |
| `nested_top_level_wrong` | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture |
| `short_indel_3path` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | path_corrupt | path_corrupt | path_corrupt |
| `snp_bubble_3path` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | path_corrupt | path_corrupt | path_corrupt |
| `tandem_copy_loop_keep` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | path_corrupt | path_corrupt | path_corrupt |
