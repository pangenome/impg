# C4 Smoke Report for Generic Graph-Output Crush

Task: `c4-smoke-report`

Branch: `wg/agent-551/c4-smoke-report`

Validated implementation commit: `b77bf99`

Artifact directory: `data/c4_smoke_report_20260606T075500Z`

## Summary

This is a bounded smoke/report validation, not a full C4 run. The primary
runtime comparison uses a 2-sequence, 226 bp-per-sequence slice derived from
`tests/test_data/crush/c4_fragments/easy_shared_flank.fa`, with all generated
validation output kept under this task directory.

The direct bounded seqwish graph-output attempts did not reach finalization:
they hit a tiny-PAF transclosure timeout before graph construction. I recorded
that as follow-up WG task `fix-seqwish-graph`. To avoid an unbounded full-C4
run while still validating the generic graph-output crush finalization path, I
generated a control-vs-graph-output-crush report using the same bounded C4
slice and the `gfa:poa` / `gfa:poa:crush` query path. The shared parser and
finalization tests for `gfa:pggb:crush` and `gfa:seqwish:crush` also passed.

## Bounds

- No full C4 run was started.
- Primary slice: 2 C4-derived sequences, 226 bp each.
- Primary runtime commands used `timeout 30s` and `-t 1`.
- Direct seqwish exploratory commands were also bounded with `timeout 30s`,
  `timeout 60s`, or `timeout 120s` and `-t 1`.
- Peak RSS for the successful graph-output crush query was 12,288 KB
  (`logs/c4_2seq_renamed.query-poa-crush.time.txt`).

## Primary Commands And Outputs

The exact single-line commands are recorded by `/usr/bin/time` in the matching
`*.time.txt` files. The primary commands were:

Control GFA:

```sh
timeout 30s target/release/impg query -d 0 -i data/c4_smoke_report_20260606T075500Z/outputs/c4_2seq_renamed.impg -a data/c4_smoke_report_20260606T075500Z/inputs/c4_easy_shared_flank_2seq.renamed.paf -r C4SMOKE_A#0#chr6:0-226 --min-transitive-len 0 -o gfa:poa --sequence-files data/c4_smoke_report_20260606T075500Z/inputs/c4_easy_shared_flank_2seq.renamed.fa --temp-dir data/c4_smoke_report_20260606T075500Z/tmp --debug-dir data/c4_smoke_report_20260606T075500Z/debug/c4_2seq_renamed.query-poa-control -O data/c4_smoke_report_20260606T075500Z/outputs/c4_2seq_renamed.query-poa-control -t 1
```

Graph-output crush GFA:

```sh
timeout 30s target/release/impg query -d 0 -i data/c4_smoke_report_20260606T075500Z/outputs/c4_2seq_renamed.impg -a data/c4_smoke_report_20260606T075500Z/inputs/c4_easy_shared_flank_2seq.renamed.paf -r C4SMOKE_A#0#chr6:0-226 --min-transitive-len 0 -o gfa:poa:crush,method=poa,max-rounds=1,polish-rounds=0 --sequence-files data/c4_smoke_report_20260606T075500Z/inputs/c4_easy_shared_flank_2seq.renamed.fa --temp-dir data/c4_smoke_report_20260606T075500Z/tmp --debug-dir data/c4_smoke_report_20260606T075500Z/debug/c4_2seq_renamed.query-poa-crush -O data/c4_smoke_report_20260606T075500Z/outputs/c4_2seq_renamed.query-poa-crush -t 1
```

Generated output GFA paths:

- `outputs/c4_2seq_renamed.query-poa-control.gfa`
- `outputs/c4_2seq_renamed.query-poa-crush.gfa`
- `outputs/c4_2seq.query-poa-crush.gfa` for the original C4-derived names,
  including embedded coordinate colons.

Generated report artifacts:

- `README.md` (this report)
- `validation/control_vs_graph_output_crush.metrics.tsv`
- `validation/c4_2seq_renamed.query-poa-control.graph-report.md`
- `validation/c4_2seq_renamed.query-poa-control.graph-report.tsv`
- `validation/c4_2seq_renamed.query-poa-crush.graph-report.md`
- `validation/c4_2seq_renamed.query-poa-crush.graph-report.tsv`
- `validation/c4_2seq_renamed.source_interval.graph-report.md`
- `validation/c4_2seq_renamed.source_interval.graph-report.tsv`

## Control Vs Graph-Output Crush

Key metrics from `validation/control_vs_graph_output_crush.metrics.tsv`:

| label | segments | links | paths | path_steps | total_segment_bp | direct_self_loop_edges | self_loop_repeat_runs | path_depth_p95 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| source | 2 | 0 | 2 | 2 | 452 | 0 | 0 | 1 |
| control_poa | 4 | 4 | 2 | 6 | 227 | 0 | 0 | 2 |
| graph_output_crush | 4 | 4 | 2 | 6 | 227 | 0 | 0 | 2 |

For this tiny slice, crush finalization preserved the control topology and path
spellings. The graph-output crush log shows the crush and final sort stages ran:

- `crush: 1 resolved, 0 bailed, 1 candidates seen across 1 rounds`
- `gfasort: pipeline 'Ygs' sorted graph in 0.000s`

Source: `logs/c4_2seq_renamed.query-poa-crush.stderr.log`.

## Path Preservation

The path preservation checks compare a FASTA-derived source GFA against the
observed GFA using `target/release/examples/compare_gfa_paths`, which checks
both path names and spelled path sequences.

Renamed C4-derived smoke names:

```text
expected_paths    2
observed_paths    2
missing_paths     0
extra_paths       0
spelling_mismatches       0
```

Files:

- `validation/c4_2seq_renamed.source_interval.gfa`
- `validation/compare_gfa_paths.query-poa-crush.stdout.tsv`
- `validation/path-name-diff.txt` (empty)

Original C4-derived names, preserving embedded coordinate colons:

```text
expected_paths    2
observed_paths    2
missing_paths     0
extra_paths       0
spelling_mismatches       0
```

Files:

- `validation/c4_2seq.source_interval.gfa`
- `validation/compare_gfa_paths.c4_2seq.query-poa-crush.stdout.tsv`
- `validation/path-name-diff.original-c4-query-poa-crush.txt` (empty)

The original C4-derived paths observed in the graph-output crush GFA were:

- `C4FIXTURE#0#easy_shared_flank:14-240:0-226`
- `C4FIXTURE#0#easy_shared_flank:330-556:0-226`

## Finalization Checks

Runtime finalization evidence:

- `logs/c4_2seq_renamed.query-poa-crush.stderr.log` records crush running and
  then `gfasort: pipeline 'Ygs' sorted graph in 0.000s`.
- `validation/c4_2seq_renamed.query-poa-crush.graph-report.tsv` reports
  `direct_self_loop_edges=0`, `adjacent_same_step_path_steps=0`, and
  `self_loop_repeat_runs=0`.

The C4 smoke slice has no direct self-loop edges, so the runtime log cannot show
a self-loop count changing. The final self-loop normalization path was verified
with the focused library test:

```sh
cargo test --release --lib graph_output_crush -- --nocapture
```

Result:

- `graph_output_crush_preserves_path_names_and_spellings_through_finalization`
  passed.
- `graph_output_crush_finalization_normalizes_self_loop_runs` passed.

The parser/finalization hook for the requested graph-output specs was verified
with:

```sh
cargo test --release --bin impg test_gfa_output_format_accepts_alignment_crush_stage_with_finalization -- --nocapture
```

Result:

- `gfa:pggb:crush` parses as the pggb engine with crush config and default
  graph sort pipeline `Ygs`.
- `gfa:seqwish:crush` parses as the seqwish engine with crush config and
  default graph sort pipeline `Ygs`.

Test logs:

- `logs/cargo-test-lib-graph_output_crush.stdout.log`
- `logs/cargo-test-lib-graph_output_crush.stderr.log`
- `logs/cargo-test-bin-gfa-output-crush-finalization.stdout.log`
- `logs/cargo-test-bin-gfa-output-crush-finalization.stderr.log`

## Direct Seqwish Attempts

These bounded attempts are recorded for transparency:

- `logs/query-seqwish-crush.time.txt`: bounded query using
  `-o gfa:seqwish:crush,method=poa,max-rounds=1,polish-rounds=0` on the
  original C4 fixture. It exited in 1.13 seconds because FastGA/FAtoGDB failed
  before graph finalization.
- `logs/easy_shared_flank.seqwish.control.time.txt`: bounded `impg graph`
  seqwish control on the six-sequence C4 fixture. It timed out at 120 seconds
  in transclosure before graph finalization.
- `logs/c4_2seq_renamed.seqwish.control.time.txt`: bounded `impg graph`
  seqwish control on the two-sequence C4-derived slice. It timed out at
  60 seconds in transclosure before graph finalization.
- `logs/c4_2seq_empty.seqwish.crush.time.txt`: bounded `impg graph`
  `seqwish:crush` with an empty PAF on the two-sequence C4-derived slice. It
  timed out at 30 seconds in transclosure before graph finalization.

Follow-up created: `fix-seqwish-graph`.

## Validation Checklist

- Bounded C4 slice or documented cap: satisfied. All runtime attempts were
  capped and used a two-sequence C4-derived smoke slice for primary validation.
- Exact command, output GFA path, and report artifact: satisfied. Exact
  commands are in `logs/*.time.txt`; GFA and report paths are listed above.
- Path names preserved against source FASTA/GFA names: satisfied. Both renamed
  and original C4-derived path-name/spelling comparisons passed with no missing,
  extra, or mismatched paths.
- Final self-loop normalization and Ygs sorting/finalization ran: satisfied by
  runtime `Ygs` evidence plus focused finalization tests that force and remove
  self-loop runs through the shared graph-output crush finalization path.
- `wg artifact` records generated report/PNG: satisfied. This report, the
  generated graph-report artifacts, primary GFAs, and path-comparison outputs
  were recorded with `wg artifact`.
