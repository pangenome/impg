# Local Compression Testbed Fast Profile Artifact Index

Profile: `fast`

## Merge Payload Note

The `resolve-merge-for` merge keeps the validated runner changes, reports,
scoreboards, fixture notes, and row-level `metrics.json` files, but intentionally
omits bulky generated artifacts from the merge payload. The full validated
fast-profile artifact tree remains recoverable from
`wg/agent-626/continue-local-compression` at commit `9464fef`, including:

- `docs/evaluations/local-compression-testbed-fast/fixtures/**/output.gfa`
- `docs/evaluations/local-compression-testbed-fast/fixtures/**/output.normalized.gfa`
- `docs/evaluations/local-compression-testbed-fast/fixtures/**/{command.sh,stdout.log,stderr.log,empty.paf}`
- historical generated C4/local run payloads under `data/**`

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
- Validation note: `docs/evaluations/local-compression-testbed-fast/validation.md`

## Coverage

- Fixtures represented: 13 (`adjacent_bubbles_joint`, `alu_like_insertion`, `dispersed_repeat_glue_break`, `duplicated_flank_context`, `fake_repeat_anchor_split`, `inversion_like_case`, `microtangle_repeat_motif`, `mid_insertion_200bp`, `nested_top_level_right`, `nested_top_level_wrong`, `short_indel_3path`, `snp_bubble_3path`, `tandem_copy_loop_keep`)
- Methods represented: 12 (`local_syng_raw`, `local_syng_crush_auto`, `local_syng_crush_poa`, `local_syng_crush_poasta`, `local_syng_crush_sweepga`, `top_flubble_nonoverlap_sweepga`, `chunk_window_smooth_or_crush`, `chunk_window_sweepga_seqwish`, `whole_region_sweepga_seqwish`, `pggb_control`, `smoothxg_control`, `pggb_plus_smoothxg_control`)
- Rows: 156
- Command statuses: pass=132, skipped=24
- Topology statuses: fail=50, not_run=24, pass=82

## Control Execution

- Control timeout: 120s per row, one thread, no hidden topology gate.
- Exact path corruption is the only hard rejection; graph-quality and topology failures remain visible rows.
- Controls are generated through repository `impg graph` engines, not standalone `pggb` or `smoothxg` binaries. Standalone binary discovery is recorded only as diagnostics.
- Mapping: `pggb_control` uses the internal `--gfa-engine pggb` pipeline with one smoothing target; `pggb_plus_smoothxg_control` uses the same internal pggb engine with two smoothing targets; `smoothxg_control` uses the pggb engine's smoothxg-style smoothing stage over an explicit empty PAF source graph because this repo does not expose a separate smooth-only CLI.
- Path names: every included control graph-producing row now round-trips exact input FASTA path names and spellings; local-tier control rows remain skipped by fast-profile fixture policy only.
- `pggb_control`: command_status pass=11, skipped=2; exact_paths not_run=2, pass=11; topology fail=11, not_run=2; tool_available true=13.
  Skips: `profile_excludes_local_fixture`.
  Availability detail: profile_excludes_local_fixture;provider=internal_impg_graph_pggb;engine_invocation=impg graph --gfa-engine pggb;external_pggb_required=false;external_smoothxg_required=false;standalone_pggb=not_found;standalone_smoothxg=not_found;pggb_mapping=internal pggb engine: seqwish induction plus smoothxg-style smoothing and gfaffix;gfaffix=/home/erikg/impg/.wg-worktrees/agent-626/target/debug/gfaffix;impg=/home/erikg/impg/.wg-worktrees/agent-626/target/debug/impg;time=/usr/bin/time;timeout=/usr/bin/timeout
- `smoothxg_control`: command_status pass=11, skipped=2; exact_paths not_run=2, pass=11; topology fail=11, not_run=2; tool_available true=13.
  Skips: `profile_excludes_local_fixture`.
  Availability detail: profile_excludes_local_fixture;provider=internal_impg_smoothxg_stage_via_pggb_engine;engine_invocation=impg graph --gfa-engine pggb;external_pggb_required=false;external_smoothxg_required=false;standalone_pggb=not_found;standalone_smoothxg=not_found;smoothxg_mapping=internal pggb engine smoothing stage over an explicit empty PAF source graph;gfaffix=/home/erikg/impg/.wg-worktrees/agent-626/target/debug/gfaffix;impg=/home/erikg/impg/.wg-worktrees/agent-626/target/debug/impg;time=/usr/bin/time;timeout=/usr/bin/timeout
- `pggb_plus_smoothxg_control`: command_status pass=11, skipped=2; exact_paths not_run=2, pass=11; topology fail=11, not_run=2; tool_available true=13.
  Skips: `profile_excludes_local_fixture`.
  Availability detail: profile_excludes_local_fixture;provider=internal_impg_graph_pggb;engine_invocation=impg graph --gfa-engine pggb;external_pggb_required=false;external_smoothxg_required=false;standalone_pggb=not_found;standalone_smoothxg=not_found;pggb_mapping=internal pggb engine: seqwish induction plus smoothxg-style smoothing and gfaffix;gfaffix=/home/erikg/impg/.wg-worktrees/agent-626/target/debug/gfaffix;impg=/home/erikg/impg/.wg-worktrees/agent-626/target/debug/impg;time=/usr/bin/time;timeout=/usr/bin/timeout

### Control Comparison

- Raw SYNG graph-producing rows: 11; exact paths pass=11; topology fail=11.
- Compact/window graph-producing rows: 88; exact paths pass=88; topology fail=6, pass=82.
- PGGB/SmoothXG control graph-producing rows: 33; exact paths pass=33; topology fail=33.


## Output Roots

- Committed per-fixture payload is limited to notes and metrics JSON:
  `docs/evaluations/local-compression-testbed-fast/fixtures/<fixture-id>/<method-id>/`
- Per-fixture command logs, graph outputs, normalized graphs, and empty PAF files
  were validated on the source branch and are intentionally omitted here by the
  generated-artifact policy described above.
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

| Fixture | `local_syng_raw` | `local_syng_crush_auto` | `local_syng_crush_poa` | `local_syng_crush_poasta` | `local_syng_crush_sweepga` | `top_flubble_nonoverlap_sweepga` | `chunk_window_smooth_or_crush` | `chunk_window_sweepga_seqwish` | `whole_region_sweepga_seqwish` | `pggb_control` | `smoothxg_control` | `pggb_plus_smoothxg_control` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `adjacent_bubbles_joint` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/fail | pass/fail | pass/fail |
| `alu_like_insertion` | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture |
| `dispersed_repeat_glue_break` | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture | skipped:profile_excludes_local_fixture |
| `duplicated_flank_context` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/fail | pass/fail | pass/fail |
| `fake_repeat_anchor_split` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/fail | pass/fail | pass/fail |
| `inversion_like_case` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/fail | pass/fail | pass/fail |
| `microtangle_repeat_motif` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/fail | pass/fail | pass/fail |
| `mid_insertion_200bp` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/fail | pass/fail | pass/fail |
| `nested_top_level_right` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/fail | pass/fail | pass/fail |
| `nested_top_level_wrong` | pass/fail | pass/fail | pass/fail | pass/fail | pass/fail | pass/fail | pass/pass | pass/pass | pass/fail | pass/fail | pass/fail | pass/fail |
| `short_indel_3path` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/fail | pass/fail | pass/fail |
| `snp_bubble_3path` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/fail | pass/fail | pass/fail |
| `tandem_copy_loop_keep` | pass/fail | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/pass | pass/fail | pass/fail | pass/fail |
