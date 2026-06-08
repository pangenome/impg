# Local Graph Compression Testbed Design

Date: 2026-06-07
Task: `design-local-graph-compression-testbed`

## Purpose

C4-scale runs are too large and ambiguous for rapid debugging of local graph
compression. This testbed defines tiny, reproducible graph and sequence
fixtures where the input path spellings and expected graph behavior are known in
advance. The fixtures are intentionally adversarial: they should expose whether
candidate-window selection, flubble traversal, chunking, SYNG rebuilds,
SweepGA/seqwish induction, POA/abPOA/SPOA polishing, and PGGB/SmoothXG-style
controls preserve paths and improve topology.

The testbed is not a replacement for C4 validation. It is the first debugger
that should fail fast when a local method corrupts paths, chooses the wrong
window, glues repeats, over-compresses loops, leaves avoidable white-space, or
looks good only because hidden filters suppressed difficult candidates.

## Batch Contract

The dependency chain is:

```text
design -> synthetic fixtures -> method runner -> scoreboard/report
```

Each stage owns a different artifact and must not rely on chat-only context:

| Stage | Owner task | Required output | Consumed by |
| --- | --- | --- | --- |
| Design | `design-local-graph-compression-testbed` | This document | Fixture, runner, and report workers |
| Synthetic fixtures | `implement-local-compression-fixtures` | Fixture manifest, fixture inputs, fixture metadata, fixture validation script or test | Runner |
| Method runner | `implement-local-compression-testbed-runner` | Method execution outputs, TSV scoreboard, JSON scoreboard, logs, optional renders | Report worker |
| Scoreboard/report | `evaluate-local-compression-testbed` | Human-readable markdown report under `docs/evaluations/` | Future implementation and defaults work |

The fixture worker must not invent a different schema. The runner worker must
consume fixture metadata instead of hard-coding fixture expectations. The report
worker must treat the runner scoreboards as the source of truth.

## Responsibility Boundaries

The testbed separates seven responsibilities. Implementations may share helper
code, but the runner and report must attribute failures to these layers where
possible.

| Responsibility | Definition | Required outputs or checks |
| --- | --- | --- |
| Sequence collection | Build or load the exact per-fixture input sequences and names. Synthetic sequence generation belongs here. | Stable FASTA or sequence files, stable path names, expected path spelling map. |
| Graph induction | Convert sequence inputs into an initial graph for each method family. Examples include raw local SYNG-derived construction and whole-region SweepGA/seqwish. | Output GFA, command line, command status, graph path count. |
| Candidate-window selection | Choose local regions for compression, smoothing, or rebuild. Examples include POVU flubbles, highest-level non-overlapping flubbles, sorted/chunk windows, and whole-region windows. | Candidate/window sidecar or log with window IDs, coordinates, source path/context, parent/child relation where relevant, and skip reasons. |
| Crush/polish/rebuild | Replace selected regions using the selected resolver. Examples include SPOA/POA, abPOA/POASTA, SweepGA/seqwish, SmoothXG-style chunk smoothing, and local SYNG rebuilds. | Per-window method parameters, replacement graph path, accepted/rejected status, exact path-spelling check result. |
| Final normalization | Apply deterministic graph normalization only after method output exists. Examples include GFA parsing, path-name normalization, optional `gfasort`/`Ygs`, and loop normalization where explicitly requested. | Normalized output GFA, normalization command, pre/post path-spelling comparison. |
| Scoring | Compare output against fixture metadata. Scoring is diagnostic except for hard path corruption. | Per-fixture/per-method metrics, topology assertion result, hard path-corruption flag. |
| Reporting | Summarize scoreboards for humans without hiding rows. | Scoreboard TSV/JSON paths, per-fixture notes, markdown report, command logs, optional PNG render paths. |

## Repository Layout

The fixture implementation should use this layout unless an existing repository
convention is clearly better:

```text
tests/test_data/local_compression/
  manifest.json
  <fixture-id>/
    input.fa
    seed.gfa                 # optional if graph input is fixture-owned
    metadata.json
    expected_paths.tsv       # optional convenience file
    notes.md                 # optional human fixture note
```

Runner outputs should use a deterministic output root, for example:

```text
target/local-compression-testbed/<profile>/
  scoreboard.tsv
  scoreboard.json
  report.md                  # optional generated stub before evaluation task
  fixtures/<fixture-id>/<method-id>/
    command.sh
    stdout.log
    stderr.log
    output.gfa
    output.normalized.gfa
    metrics.json
    render.png               # optional
```

The fixture manifest must list every fixture ID, metadata path, tier, and class.
The runner must discover fixtures from this manifest.

## Fixture Metadata Schema

Each fixture metadata file must contain the fields below. JSON is preferred for
machine readability. YAML is acceptable only if the repository already has a
YAML parser in the fixture or runner path.

```json
{
  "fixture_id": "snp_bubble_3path",
  "fixture_class": "snp_bubble",
  "title": "Three-path SNP bubble",
  "tier": "ci",
  "input_sequences": [
    {
      "name": "ref#0#chrTest:1-33",
      "file": "input.fa",
      "expected_spelling": "ACGT..."
    }
  ],
  "input_sequence_names": [
    "ref#0#chrTest:1-33",
    "hapA#0#chrTest:1-33",
    "hapB#0#chrTest:1-33"
  ],
  "expected_path_spellings": {
    "ref#0#chrTest:1-33": "ACGT...",
    "hapA#0#chrTest:1-33": "ACGA...",
    "hapB#0#chrTest:1-33": "ACGG..."
  },
  "expected_topology": {
    "assertion_id": "single_ordered_snp_bubble",
    "description": "One source-to-sink SNP bubble, no loop, no long links.",
    "exact": {
      "bubble_count": 1,
      "flubble_count_min": 1,
      "self_loops": 0,
      "repeat_loops": 0,
      "max_long_link_span_bp": 0,
      "branch_order": ["ref", "hapA", "hapB"]
    }
  },
  "allowed_ranges": {
    "segment_count": {"min": 3, "max": 8},
    "link_count": {"min": 4, "max": 12},
    "path_depth_median": {"min": 1, "max": 3},
    "node_depth_p95": {"min": 1, "max": 3},
    "white_space_proxy_bp_p95": {"min": 0, "max": 0},
    "bubble_count": {"min": 1, "max": 1},
    "flubble_count": {"min": 1, "max": 2}
  },
  "long_link_policy": {
    "max_count": 0,
    "max_span_bp": 0
  },
  "known_failure_mode": "Path-corrupting crush can swap the SNP branch spelling.",
  "render_hints": {
    "reference_path": "ref#0#chrTest:1-33",
    "highlight_paths": ["hapA#0#chrTest:1-33", "hapB#0#chrTest:1-33"],
    "expected_shape": "single small bubble"
  }
}
```

Required top-level fields:

| Field | Required | Meaning |
| --- | --- | --- |
| `fixture_id` | Yes | Stable lowercase ID used in file paths and scoreboard rows. |
| `fixture_class` | Yes | One of the required fixture classes listed in this design. |
| `title` | Yes | Short human-readable name. |
| `tier` | Yes | CI/local-only tier. Use `ci` for the fast subset, `local` for heavier local-only benchmark fixtures. |
| `input_sequences` | Yes | Input sequence records, files, and expected spelling per record. |
| `input_sequence_names` | Yes | Stable path/sequence names. Prefer PanSN-compatible names or document any deviation. |
| `expected_path_spellings` | Yes | Exact expected spelled sequence for every path that every method must preserve. |
| `expected_topology` | Yes | Topology/compression assertion. Simple fixtures should use exact motifs; messy fixtures may combine motifs with explicit ranges. |
| `allowed_ranges` | Yes when exact topology is brittle | Explicit accepted ranges for size/depth/topology proxies. Empty or omitted only when exact topology is complete. |
| `known_failure_mode` | Yes | The algorithmic failure this fixture is meant to expose. |
| `render_hints` | Optional | Reference path, highlighted paths, expected visual shape, and whether failures should render by default. |

`expected_topology.assertion_id` must be stable because runner rows should emit
the assertion ID and pass/fail message. The runner may add derived assertion
messages, but it must not reinterpret a fixture without the metadata field.

## Required Fixture Set

The fixture worker must implement at least the fixture classes below. More than
one fixture may be added per class when needed, but every class must appear in
the manifest. The IDs below are recommended stable IDs.

| Fixture ID | Class | Tier | Expected behavior | Known failure mode |
| --- | --- | --- | --- | --- |
| `snp_bubble_3path` | `snp_bubble` | CI | One small ordered bubble; exact path spelling; no loop; no long links. | Branch path spellings get swapped or the SNP remains as path-specific singleton nodes. |
| `short_indel_3path` | `short_indel` | CI | One indel bubble with insertion/deletion branch; no long links; compact shared flanks. | Resolver clips inserted bases, creates a dangling tip, or collapses deletion into the wrong path. |
| `mid_insertion_200bp` | `insertion_50_500bp` | CI | 50-500 bp inserted branch remains spelled exactly and shares both flanks. | POA/SYNG crush over-fragments the insertion or drops the insertion branch. |
| `alu_like_insertion` | `alu_like_insertion` | local | Alu-sized, repeat-rich inserted branch with exact spelling; topology may use size/depth ranges. | Repeat content glues unrelated anchors or fails under short exact-match floors. |
| `adjacent_bubbles_joint` | `adjacent_bubbles_compress_together` | CI | Two adjacent variants should be selected and compressed as one local problem or as a topology-equivalent joint result. | Per-leaf flubble boundaries leave two stringy bubbles separated by artificial white-space. |
| `fake_repeat_anchor_split` | `bubble_split_by_fake_repeat_anchor` | CI | A bubble containing a repeated internal k-mer should not split into two unrelated events. | Candidate discovery treats the repeat as a hard anchor and under-compresses the true event. |
| `microtangle_repeat_motif` | `repeated_motif_microtangle` | CI | Short repeated motif creates a bounded microtangle; exact spelling mandatory; accepted topology ranges explicit. | Local graph explodes, self-loops appear, or all motif copies collapse into one unsupported node. |
| `duplicated_flank_context` | `duplicated_flank_requires_path_context` | CI | Duplicate flank sequence appears in two positions; selection must use path context to choose the right occurrence. | Replacement laces into the wrong duplicate flank or merges distinct occurrences. |
| `tandem_copy_loop_keep` | `tandem_copy_number_loop_cyclic` | CI | Tandem copy-number variation should preserve a cyclic or loop-capable representation when the fixture declares a loop required. | Normalization or crush linearizes the copy-number loop and loses cyclic topology. |
| `dispersed_repeat_glue_break` | `dispersed_repeat_glue_break_or_ignore` | local | Dispersed repeat similarity should be broken or ignored; accepted graph should avoid long repeat-glue links. | seqwish/SYNG transitive closure glues distant repeat copies into long links. |
| `inversion_like_case` | `inversion_like` | CI | Inversion-like branch should preserve orientation-sensitive path spellings and avoid illegal cross-orientation lacing. | Resolver reverse-complements the wrong interval or creates long links between inverted anchors. |
| `nested_top_level_right` | `nested_bubbles_top_level_right` | CI | Nested bubbles should be compressed by top-level non-overlapping boundaries; descendant details included in one replacement. | Leaf-only selection leaves avoidable nested white-space or repeated replacement conflicts. |
| `nested_top_level_wrong` | `nested_bubbles_top_level_wrong` | local | Top-level compression is intentionally wrong; smaller/chunked boundaries should outperform it. | Top-level selection over-merges independent events and creates bad loops, long links, or excess path depth. |

Minimum CI subset:

```text
snp_bubble_3path
short_indel_3path
mid_insertion_200bp
adjacent_bubbles_joint
fake_repeat_anchor_split
microtangle_repeat_motif
duplicated_flank_context
tandem_copy_loop_keep
inversion_like_case
nested_top_level_right
```

Minimum local-only benchmark subset:

```text
alu_like_insertion
dispersed_repeat_glue_break
nested_top_level_wrong
```

The local subset may be heavier but must remain small enough for local
iteration. It should not require C4 or HPRC indexes. A local-only fixture should
target seconds to a few minutes, not hours.

## Expected Topology Guidance

Exact path spelling is mandatory for every method and every fixture result.
This is the only correctness condition allowed to abort or reject a produced
method output.

Simple fixtures should require exact topology motifs whenever feasible:

- expected bubble/flubble count;
- loop presence or absence;
- branch ordering when order is deterministic after normalization;
- no long links beyond fixture-defined limits;
- expected shared left and right flanks;
- expected path count and path-name stability.

Messy fixtures may use ranges when exact topology is too brittle. Ranges must be
fixture-specific and recorded in metadata. Allowed range examples:

- `segment_count.min/max`;
- `link_count.min/max`;
- `path_depth_median.min/max`;
- `path_depth_p95.min/max`;
- `node_depth_p95.min/max`;
- `white_space_proxy_bp_p95.min/max`;
- `white_space_proxy_bp_max.min/max`;
- `bubble_count.min/max`;
- `flubble_count.min/max`;
- `repeat_loop_count.min/max`;
- `long_link_count.min/max`;
- `max_long_link_span_bp.min/max`.

If a fixture permits a loop, the metadata must say whether the loop is required,
optional, or forbidden. For example, `tandem_copy_loop_keep` should mark a
repeat/copy-number loop as required, while `snp_bubble_3path` should forbid
self-loops and repeat loops.

## Method Matrix

The runner must emit a row for every fixture and every method family below. A
method row may be `pass`, `fail`, `path_corrupt`, `error`, or `skipped`, but it
must not disappear.

Mandatory methods are repository-local capabilities. A missing mandatory method
is a runner/testbed failure unless this design or a later committed design
revision explicitly marks the method deferred. Optional controls are external
or heavy tools and may be skipped with precise reasons.

### Mandatory Methods

| Method ID | Family | Required behavior | Typical command source | Notes |
| --- | --- | --- | --- | --- |
| `local_syng_raw` | Raw/local SYNG-derived graph construction | Build the fixture using the local SYNG-derived graph path with no crush/polish replacement. | `impg query -o gfa:syng...` or an equivalent fixture-local SYNG/GFA helper. | Baseline for path collection and native topology. |
| `local_syng_crush_auto` | SYNG plus existing flubble/local crush variant | Run the current default local crush behavior. | `gfa:syng:...:crush,method=auto` or equivalent. | Must record selected windows and resolver choices. |
| `local_syng_crush_poa` | SYNG plus POA/SPOA local crush | Force POA/SPOA where supported by existing crush options. | `crush,method=poa` or `polish-method=poa`. | Useful for tiny SNP/indel/insertion fixtures. |
| `local_syng_crush_poasta` | SYNG plus abPOA/POASTA local crush | Force abPOA/POASTA where supported. | `crush,method=poasta` or equivalent. | Should expose insertion and adjacent-bubble behavior. |
| `local_syng_crush_sweepga` | SYNG plus SweepGA/seqwish local crush | Force local SweepGA/seqwish replacement on eligible windows with bounded parameters. | `crush,method=sweepga,seqwish-k=<preset>` or equivalent. | Must keep exact path preservation as hard guard only. |
| `top_flubble_nonoverlap_sweepga` | Highest-level non-overlapping flubble-window crush | Select only highest-level non-overlapping POVU flubbles and rebuild each selected top-level window. | `method=top-flubble-sweepga` / `top-level-flubble-sweepga` / `chain-povu-sweepga`. | Required to test nested fixtures where top-level is right or wrong. |
| `chunk_window_smooth_or_crush` | SmoothXG-style sorted/chunk-window smoothing or crush where locally available | Partition by sorted/chunk windows rather than exact flubble boundaries, then smooth/crush each window with local repo capabilities. | Existing `smooth`, `crush`, or graph pipeline commands if present. | This is a local method family. If the repo has no local implementation, emit `skipped` with `mandatory_unavailable` and fail validation until implemented or deferred. |
| `whole_region_sweepga_seqwish` | Whole-region SweepGA/seqwish induction | Ignore local candidate boundaries and induce a graph over all fixture sequences using bounded SweepGA/seqwish. | `impg graph --gfa-engine seqwish` or equivalent. | Baseline for whether the whole tiny region can be reconstructed cleanly. |

Each mandatory method row must include the exact parameter string used. For
bounded tiny-fixture presets, prefer explicit values for threads, `seqwish-k`,
filter/no-filter mode, maximum PAF bytes if applicable, sort pipeline, and
polish rounds. Defaults are not enough for reproducibility.

### Optional Controls

| Method ID | Family | Run condition | Skip reason examples |
| --- | --- | --- | --- |
| `pggb_control` | PGGB-style whole-region control | Run when PGGB or repository `gfa:pggb` wrapper is installed and bounded enough for the fixture profile. | `tool_not_found`, `profile_excludes_optional`, `estimated_runtime_exceeds_limit`, `unsupported_fixture_input`. |
| `smoothxg_control` | External SmoothXG control | Run when SmoothXG is installed, accepts the fixture graph/input, and fits the profile budget. | `tool_not_found`, `profile_excludes_optional`, `input_format_unsupported`, `estimated_runtime_exceeds_limit`. |
| `pggb_plus_smoothxg_control` | Full PGGB/SmoothXG-style external control | Run only in local benchmark profile when both tools are available and bounded. | Same as above. |

Missing optional controls are not failures. They must still produce scoreboard
rows with `command_status=skipped` and `skipped_optional_tool_reason`.

## SmoothXG-Style Chunking Tradeoff

Exact flubble boundaries can be too literal on sparse or messy graphs. A POVU
flubble boundary follows the current graph topology, which may contain fake
anchors, nested artifacts, duplicated flanks, or repeat-glue nodes. If the
current graph is under-aligned, the "exact" flubble may be a symptom of the
problem instead of the right replacement unit.

SmoothXG-style chunking can outperform exact flubble boundaries because it uses
sorted or bounded windows that include nearby ambiguous context. That wider
context can:

- align adjacent bubbles together instead of preserving artificial separation;
- ignore a fake internal repeat anchor that split the event;
- give duplicated flanks enough path context to avoid wrong lacing;
- smooth sparse local graphs where a flubble boundary contains too little
  alignment evidence;
- avoid top-level over-merging by using bounded chunks rather than a single
  parent flubble.

Fixtures intended to expose this tradeoff:

| Fixture | Expected tradeoff |
| --- | --- |
| `adjacent_bubbles_joint` | Chunk windows should be able to compress adjacent bubbles together; exact leaf flubbles may leave two artifacts. |
| `fake_repeat_anchor_split` | Chunk windows should bridge a fake repeat anchor; exact flubble boundaries may split the true event. |
| `microtangle_repeat_motif` | Chunk windows may reduce sparse motif white-space, but too-wide chunks can over-collapse copies. |
| `duplicated_flank_context` | Chunk windows with path context should avoid wrong duplicate-flank lacing better than context-free local boundaries. |
| `dispersed_repeat_glue_break` | Chunking should help only if it breaks/ignores repeat glue; whole-region transitive closure may do worse. |
| `nested_top_level_right` | Highest-level flubble should win or tie because the parent event is the correct unit. |
| `nested_top_level_wrong` | Chunk windows or smaller local selection should beat top-level flubble because the parent boundary over-merges independent events. |

The report must answer where exact flubble boundaries work and where
SmoothXG-style chunks work better using scoreboard evidence, not only visual
preference.

## Metric Schema

The runner must emit the fields below in both TSV and JSON. TSV column names
should match these names where possible. JSON may nest structures, but the same
information must be present.

### Identity and Provenance

| Field | Required | Meaning |
| --- | --- | --- |
| `fixture_id` | Yes | Fixture ID from metadata. |
| `fixture_class` | Yes | Fixture class from metadata. |
| `tier` | Yes | `ci` or `local`. |
| `method_id` | Yes | Method ID from the matrix. |
| `method_family` | Yes | Method family from the matrix. |
| `method_parameters` | Yes | Reproducible parameter string or JSON object. |
| `profile` | Yes | `ci`, `local`, or another explicit runner profile. |
| `input_manifest_path` | Yes | Fixture manifest path used by runner. |
| `metadata_path` | Yes | Fixture metadata path. |
| `output_gfa_path` | Yes for produced graphs | Final graph path, blank only for skipped/error rows without graph output. |
| `normalized_gfa_path` | Yes when normalization ran | Normalized graph path. |
| `command_log_path` | Yes | Path to raw command script/log. |
| `stdout_log_path` | Yes | Captured stdout. |
| `stderr_log_path` | Yes | Captured stderr. |

### Correctness and Topology

| Field | Required | Meaning |
| --- | --- | --- |
| `exact_path_preservation` | Yes | `pass`, `fail`, or `not_run`. |
| `hard_path_corruption` | Yes | Boolean. True only when exact path spelling differs, path names are missing/extra, or paths cannot be parsed. |
| `path_corruption_detail` | Yes when corrupted | Path names and mismatch summary. |
| `expected_topology_status` | Yes | `pass`, `fail`, `not_run`, or `not_applicable`. |
| `expected_topology_assertion_id` | Yes | Assertion ID from metadata. |
| `expected_topology_message` | Yes | Human-readable result message. |
| `missing_path_count` | Yes | Number of expected paths absent from output. |
| `extra_path_count` | Yes | Number of unexpected paths in output. |
| `path_name_stable` | Yes | Boolean or `not_run`. |

### Graph Size and Shape

| Field | Required | Meaning |
| --- | --- | --- |
| `graph_size_bytes` | Yes | Output graph file size. |
| `segment_count` | Yes | GFA segment/node count. |
| `link_count` | Yes | GFA link count. |
| `path_count` | Yes | GFA path count. |
| `total_segment_bp` | Yes | Sum of segment sequence lengths where available. |
| `total_path_steps` | Yes | Total path steps. |
| `node_depth_min` | Yes | Minimum node/path occupancy depth. |
| `node_depth_p05` | Yes | Node-depth p05. |
| `node_depth_median` | Yes | Node-depth median. |
| `node_depth_p95` | Yes | Node-depth p95. |
| `node_depth_max` | Yes | Node-depth max. |
| `path_depth_min` | Yes | Per-path depth or occupancy min under the runner's depth definition. |
| `path_depth_p05` | Yes | Path-depth p05. |
| `path_depth_median` | Yes | Path-depth median. |
| `path_depth_p95` | Yes | Path-depth p95. |
| `path_depth_max` | Yes | Path-depth max. |
| `white_space_proxy_bp_total` | Yes | Total path white-space proxy or equivalent sorted-order gap proxy. |
| `white_space_proxy_bp_p95` | Yes | p95 white-space proxy. |
| `white_space_proxy_bp_p99` | Yes | p99 white-space proxy. |
| `white_space_proxy_bp_max` | Yes | Maximum white-space proxy. |
| `self_loop_count` | Yes | Direct self-loop edge count. |
| `repeat_loop_count` | Yes | Repeat/tandem-loop count under runner heuristic. |
| `bubble_count` | Yes | Bubble count where measurable. |
| `flubble_count` | Yes | POVU flubble/site count where measurable. |
| `long_link_count` | Yes | Count of links exceeding fixture or profile threshold. |
| `long_link_max_span_bp` | Yes | Maximum long-link span if measurable. |

### Runtime and Status

| Field | Required | Meaning |
| --- | --- | --- |
| `runtime_seconds` | Yes | Wall-clock runtime for method command. |
| `max_rss_kb` | Preferred | Peak RSS when available from `/usr/bin/time -v` or equivalent. |
| `command_status` | Yes | `pass`, `fail`, `path_corrupt`, `error`, or `skipped`. |
| `exit_code` | Yes for commands | Process exit code. |
| `skip_reason` | Yes for skipped rows | General skip reason. |
| `skipped_optional_tool_reason` | Yes for skipped optional controls | Precise optional-tool reason. |
| `mandatory_unavailable` | Yes | Boolean. True when a mandatory method family cannot run. |
| `candidate_count` | Preferred | Number of selected candidate windows. |
| `candidate_skipped_count` | Preferred | Number of explicit candidate skips. Must not hide the row. |
| `candidate_skip_reasons` | Preferred | Aggregated explicit reasons. |

Every metric that cannot be measured must be present with `null`, `NA`, or a
documented sentinel plus a reason field. Missing columns are not acceptable.

## Output Artifacts

The runner and report stages must produce these human-readable and
machine-readable artifacts:

| Artifact | Required | Contents |
| --- | --- | --- |
| `scoreboard.tsv` | Yes | One row per fixture/method with every required metric column. |
| `scoreboard.json` | Yes | One object per fixture/method with equivalent data and nested details. |
| Per-fixture notes | Yes | Expected behavior, observed pass/fail state, notable regressions, and method-specific comments. |
| Markdown summary/report | Yes | Fixture-by-method matrix, links/paths to scoreboards, command reproduction detail, recommendations. |
| Raw command logs | Yes | Command script, stdout, stderr, exit status, timing. |
| Output graph paths | Yes for produced graphs | Raw and normalized GFA paths. |
| Optional PNG renders | Optional but recommended | Required for failures/regressions when render tooling is available; representative controls may also render. |

PNG renders should be generated for:

- hard path-corruption cases if a graph was produced and parseable enough to
  render;
- expected-topology failures;
- regressions against a previous committed scoreboard;
- representative successes for optional PGGB/SmoothXG controls when run.

If renders are skipped, the notes/report must say whether the reason was
missing render tooling, unparseable graph, profile exclusion, or time budget.

## No-Hidden-Filtering Rule

No hidden quality filters, path guards, candidate suppression, or silent
acceptance gates are allowed in this testbed.

The only permitted hard guard is aborting or rejecting a method output when
exact path spelling is corrupted. A hard path-corruption row must still be
recorded in TSV/JSON and in the report.

All other quality problems must be scored and reported, not hidden:

- too many segments or links;
- excess white-space proxy;
- self-loops or repeat loops when forbidden;
- missing expected loops when required;
- too many/few bubbles or flubbles;
- long links;
- bad branch ordering;
- poor candidate selection;
- resolver timeout or nonzero exit;
- optional tool unavailable;
- no backend alignments;
- graph much larger than baseline.

Candidate suppression is allowed only when it is an explicit method parameter
or an explicit cost/profile bound recorded in `method_parameters`,
`candidate_skip_reasons`, and the command logs. It must not silently make a
method look better by removing difficult windows from the denominator.

## Profiles and Budgets

The testbed has two standard profiles.

### CI Profile

Purpose: deterministic fast smoke coverage for every mandatory method that is
safe to run in CI.

Required properties:

- includes every `tier=ci` fixture;
- excludes `tier=local` fixtures by default;
- runs mandatory methods with tiny bounded settings;
- emits skipped rows for optional controls unless explicitly enabled;
- should target minutes, not hours;
- must validate exact path spelling for every produced graph.

Recommended CI budget targets:

| Budget | Target |
| --- | ---: |
| Per fixture/method wall time | <= 30 seconds when practical |
| Whole CI profile | <= 10 minutes when practical |
| Sequence length per CI fixture | Usually <= 2 kbp, except explicitly justified microtangle/context fixtures |
| Path count per CI fixture | Usually 3-8 paths |

### Local Benchmark Profile

Purpose: richer comparison for repeat-heavy or nested cases and optional
external controls.

Required properties:

- includes CI fixtures and local-only fixtures unless a narrower fixture filter
  is specified;
- may run optional PGGB/SmoothXG controls when installed and bounded;
- records time/RSS and skip reasons;
- should remain local-debuggable without C4 or HPRC indexes.

Recommended local budget targets:

| Budget | Target |
| --- | ---: |
| Per fixture/method wall time | <= 5 minutes unless explicitly marked heavier |
| Whole local profile | <= 60 minutes when optional controls are enabled |
| Sequence length per local fixture | Usually <= 20 kbp |
| Path count per local fixture | Usually <= 16 paths |

## Scoring Semantics

The runner should compute row-level status as follows:

| Status | Meaning |
| --- | --- |
| `pass` | Exact path spelling passes and fixture topology/ranges pass. |
| `fail` | Exact path spelling passes, but topology/ranges or graph-quality assertions fail. |
| `path_corrupt` | Exact path spelling, path names, or path parsing failed. This is the only hard rejection. |
| `error` | Command failed before producing a scoreable graph, or metrics could not be computed for a non-path-corruption reason. |
| `skipped` | Method intentionally not run with an explicit reason. |

The report must separate:

- hard correctness failures (`path_corrupt`);
- topology assertion failures;
- graph-quality diagnostics such as white-space, loops, long links, and
  depth/size anomalies;
- runtime/tool failures;
- optional-tool availability.

A method can be useful even when it fails one fixture. The report should name at
least one working method per fixture class or state that none works.

## Implementation Notes for Fixture Workers

Fixture sequences should be readable by eye. Prefer short synthetic alphabets
over opaque random strings, but avoid repetitive sequences so simple fixtures do
not accidentally create repeat glue. For repeat fixtures, make the repeated
motif intentional and document it in `known_failure_mode`.

Path names should be stable and meaningful. PanSN-compatible names are
preferred:

```text
sample#hap#contig:start-end
```

If a fixture uses a graph input, the graph must spell exactly the same
`expected_path_spellings` as the sequence input. A fixture validation script or
cargo test must prove this before the runner executes the method matrix.

## Implementation Notes for Runner Workers

Reuse existing graph and path validation helpers where possible. Do not invent
parallel ad hoc GFA parsing if repository helpers can parse paths, compute graph
metrics, decompose POVU flubbles, render reports, or validate path spellings.

Every command should be reproducible from its command log. Environment details
that affect output, such as `TMPDIR`, thread counts, sort pipeline, and external
tool paths, should be captured.

The runner should run exact path comparison in two places when practical:

1. immediately after graph induction or method output;
2. after final normalization.

If normalization corrupts paths, mark the row `path_corrupt` and include whether
the pre-normalized graph was still valid.

## Implementation Notes for Report Workers

The report should include:

- scoreboard TSV and JSON paths;
- command/profile summary;
- fixture-by-method matrix;
- one short section per fixture class;
- optional control availability and skip reasons;
- representative PNG paths or a render-skip explanation;
- answer to whether "fast like SYNG, compressed like PGGB" is credible from the
  local evidence;
- actionable next-step recommendation tied to scoreboard rows.

Do not drop inconvenient rows. If a method corrupts paths, keep it in the matrix
as a correctness failure and exclude it only from graph-quality comparisons that
require valid path spellings.

## Acceptance Checklist

This design is complete when downstream workers can answer all of the following
without chat context:

- where fixtures live and how they are discovered;
- what metadata every fixture must define;
- which fixture classes are required and which are CI versus local-only;
- which method rows must appear in the scoreboard;
- which methods are mandatory versus optional controls;
- which metrics must be emitted;
- which artifacts must be written;
- how exact path preservation, topology assertions, and graph-quality metrics
  are separated;
- which hard guard is allowed;
- why SmoothXG-style chunk windows are being compared with exact flubble
  boundaries.
