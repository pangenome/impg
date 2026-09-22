# Local Compression Testbed Scoreboard Report

Task: `evaluate-local-compression-testbed`
Date: 2026-06-07
Status: **blocked by missing runner outputs**
Follow-up created: `follow-up-produce`

## Source-of-truth preflight

This report did not run a comparative quality analysis because the required
runner outputs from `implement-local-compression-testbed-runner` are absent.
The task contract says to treat the runner scoreboard as the source of truth
and to stop rather than infer results from chat context when a required runner
output is missing or malformed.

Required runner artifacts were not available:

| Required input | Status |
| --- | --- |
| Scoreboard TSV path | **missing**; no task-scoped `scoreboard.tsv` exists |
| Scoreboard JSON path | **missing**; no task-scoped `scoreboard.json` exists |
| Fixture manifest | **missing**; no `tests/test_data/local_compression/manifest.json` exists |
| Per-fixture notes | **missing** |
| Command logs | **missing** |
| Graph output paths | **missing** |
| Optional PNG renders | **missing** |

Evidence checked:

```bash
wg msg read evaluate-local-compression-testbed --agent "$WG_AGENT_ID"
wg show evaluate-local-compression-testbed
wg show implement-local-compression-testbed-runner
rg --files docs/evaluations
rg -n "scoreboard|local compression|testbed|runner|fixture|render|PNG|json|tsv" docs/evaluations
rg --files | rg 'scoreboard|local.*compression|compression.*testbed|testbed|fixture|render|\.tsv$|\.json$|\.png$'
find . -path './.git' -prune -o -path './target' -prune -o -path './.wg' -prune -o \
  \( -name '*scoreboard*.tsv' -o -name '*scoreboard*.json' -o -name '*local*compression*' \
     -o -path '*local_compression*' -o -path '*compression_testbed*' \) -print
find tests/test_data -maxdepth 4 -type f | sort
```

Those checks found only design/evaluation documents:

- `docs/local-graph-compression-testbed-design.md`
- `docs/local-graph-compression-testbed-quality-pass.md`
- `docs/evaluations/implement-local-compression-fixtures.md`
- `docs/evaluations/implement-local-compression-testbed-runner.md`
- this report

The dependency evaluation at
`docs/evaluations/implement-local-compression-testbed-runner.md` also records
that no runner implementation, scoreboard TSV/JSON, per-fixture notes, command
logs, graph paths, or renderer outputs were produced.

## Reproduction detail

There is no fast-profile runner command to rerun because the runner artifact is
missing. The design document specifies the intended deterministic output shape:

```text
target/local-compression-testbed/<profile>/
  scoreboard.tsv
  scoreboard.json
  report.md
  fixtures/<fixture-id>/<method-id>/
    command.sh
    stdout.log
    stderr.log
    output.gfa
    output.normalized.gfa
    metrics.json
    render.png
```

The follow-up task `follow-up-produce` was created to implement the missing
fixtures/runner and produce real fast-profile scoreboards before this
comparison is attempted again.

## Fixture-by-method matrix

Legend:

- `NO ROW`: no scoreboard row exists, so the method cannot be scored.
- `NO SKIP ROW`: optional-control availability was not recorded by a runner;
  this is not evidence that the external tool is installed or missing.

| Fixture ID | Class | local_syng_raw | local_syng_crush_auto | local_syng_crush_poa | local_syng_crush_poasta | local_syng_crush_sweepga | top_flubble_nonoverlap_sweepga | chunk_window_smooth_or_crush | whole_region_sweepga_seqwish | pggb_control | smoothxg_control | pggb_plus_smoothxg_control |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `snp_bubble_3path` | `snp_bubble` | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO SKIP ROW | NO SKIP ROW | NO SKIP ROW |
| `short_indel_3path` | `short_indel` | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO SKIP ROW | NO SKIP ROW | NO SKIP ROW |
| `mid_insertion_200bp` | `insertion_50_500bp` | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO SKIP ROW | NO SKIP ROW | NO SKIP ROW |
| `alu_like_insertion` | `alu_like_insertion` | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO SKIP ROW | NO SKIP ROW | NO SKIP ROW |
| `adjacent_bubbles_joint` | `adjacent_bubbles_compress_together` | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO SKIP ROW | NO SKIP ROW | NO SKIP ROW |
| `fake_repeat_anchor_split` | `bubble_split_by_fake_repeat_anchor` | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO SKIP ROW | NO SKIP ROW | NO SKIP ROW |
| `microtangle_repeat_motif` | `repeated_motif_microtangle` | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO SKIP ROW | NO SKIP ROW | NO SKIP ROW |
| `duplicated_flank_context` | `duplicated_flank_requires_path_context` | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO SKIP ROW | NO SKIP ROW | NO SKIP ROW |
| `tandem_copy_loop_keep` | `tandem_copy_number_loop_cyclic` | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO SKIP ROW | NO SKIP ROW | NO SKIP ROW |
| `dispersed_repeat_glue_break` | `dispersed_repeat_glue_break_or_ignore` | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO SKIP ROW | NO SKIP ROW | NO SKIP ROW |
| `inversion_like_case` | `inversion_like` | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO SKIP ROW | NO SKIP ROW | NO SKIP ROW |
| `nested_top_level_right` | `nested_bubbles_top_level_right` | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO SKIP ROW | NO SKIP ROW | NO SKIP ROW |
| `nested_top_level_wrong` | `nested_bubbles_top_level_wrong` | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO ROW | NO SKIP ROW | NO SKIP ROW | NO SKIP ROW |

No fixture class has a working method in the source-of-truth scoreboard because
there is no source-of-truth scoreboard. This is not a quality failure by any
specific method; it is an artifact-generation failure in the upstream handoff.

## Diagnostic separation

The requested correctness and quality categories are separated below. None can
be collapsed into another category because no runner row exists to evaluate it.

| Category | Result | Interpretation |
| --- | --- | --- |
| Exact path preservation | Not run | No output graph or path-validation result exists. No hard path-corruption row exists. |
| Hard correctness failure | Not assessable | There may be no path corruption, but the absence of checks means correctness is unknown. |
| Expected topology | Not run | No fixture metadata assertion was consumed and no assertion status/message exists. |
| Excess white-space | Not measured | No `white_space_proxy_*` metrics exist. |
| Self-loops or repeat loops | Not measured | No loop-count metrics exist. |
| Long links | Not measured | No long-link count or max-span metrics exist. |
| Bubble/flubble anomalies | Not measured | No bubble/flubble counts or candidate-window sidecars exist. |
| Runtime/tool failures | Not recorded | No command scripts, stdout/stderr logs, exit codes, timing, or RSS metrics exist. |
| Optional PGGB/SmoothXG controls | Not recorded | No skipped optional-control rows exist, so no tool-availability conclusion is possible from the runner. |

## Required questions

Which failures are due to aligner/resolver behavior versus candidate-window
selection?

No source-of-truth rows exist with candidate counts, selected windows, resolver
choices, exit statuses, graph diagnostics, or path-preservation results. The
failure attribution cannot be made without inventing evidence.

When do exact flubble boundaries work, and when do SmoothXG-style chunks work
better?

No exact-flubble rows and no chunk-window rows exist. The intended fixtures for
this comparison are listed in the design, but this report cannot claim either
strategy works better until the runner emits rows for both method families.

Does highest-level non-overlapping flubble selection help on nested/messy
fixtures?

No `top_flubble_nonoverlap_sweepga` rows exist for `nested_top_level_right`,
`nested_top_level_wrong`, or any other fixture. This remains untested.

Is whole-region SweepGA/seqwish close enough to PGGB/SmoothXG quality on small
fixtures, and where does it fail?

No `whole_region_sweepga_seqwish`, `pggb_control`, `smoothxg_control`, or
`pggb_plus_smoothxg_control` rows exist. Optional controls were not recorded as
run or skipped, so there is no quality baseline.

Which method should be the default local graph path for users who want fast but
PGGB-like compression?

No default should be selected from this testbed yet. The only evidence-backed
recommendation is to keep SYNG as collection-only for this decision path until
the fast-profile scoreboard exists. Promoting `gfa:pggb`, chunk windows,
top-level flubble traversal, or whole-region SweepGA/seqwish would be arbitrary
without the fixture-by-method results.

## Renders

No representative PNG renders are available. This is not a render-tool failure
reported by a runner row; the graph outputs and render commands were never
produced. The follow-up runner task must either generate representative failure
renders or emit explicit render-skip reasons such as missing render tooling,
unparseable graph, profile exclusion, or time budget.

## Optional controls

No optional-control rows exist. Therefore this report does not state that PGGB,
SmoothXG, or PGGB+SmoothXG were skipped because tools were unavailable. The
runner must make that distinction in TSV/JSON with
`skipped_optional_tool_reason`, for example `tool_not_found`,
`profile_excludes_optional`, `input_format_unsupported`, or
`estimated_runtime_exceeds_limit`.

## Credible path and bottleneck

There may be a credible path toward "fast like SYNG, compressed like PGGB", but
the present bottleneck is earlier than method tuning: the fixture manifest,
runner, scoreboard, logs, graph outputs, and render records are missing. Until
those exist, the project cannot distinguish path corruption from topology
diagnostics, candidate-window failures from resolver failures, or local
repository methods from external optional controls.

## Recommendation

Do not productionize a default local compression path from this testbed yet.
The next action is to complete `follow-up-produce`: implement the fixture
manifest and runner, run the fast profile, and emit complete TSV/JSON rows for
all mandatory method families plus optional-control skip rows. After that,
rerun this report from the scoreboard source of truth and choose between
SmoothXG-style chunk windows, highest-level flubble traversal, whole-region
SweepGA/seqwish, or keeping SYNG as collection-only based on exact path
preservation first and graph-quality diagnostics second.
