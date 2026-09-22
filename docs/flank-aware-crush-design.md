# Flank-Aware Crush Design

Task: `design-flank-aware-crush`
Date: 2026-06-06
Status: implementation design for `implement-flank-aware-crush`

## Purpose

Bubble-local and motif-local crush sometimes aligns only the variable-looking
bytes inside a local region. For short indels, tandem motifs, homopolymers, and
microtangles, that interior can lack enough common context for global
SPOA/abPOA/POASTA/SweepGA/seqwish induction to keep homologous bases together.

The canonical model is:

1. Define the exact replaceable target interval first.
2. Add occurrence-local left and right flanks only as resolver context.
3. Feed all resolver inputs in one explicit target orientation.
4. Run the chosen resolver in global/end-to-end mode on the flanked sequences.
5. Trim the resolved graph back to the original target interval.
6. Lace only the trimmed replacement over the original interval.
7. Accept the replacement only if full-graph path names and path spellings are
   unchanged.

The flank sequence itself is never part of the replacement interval.

## Current State To Refactor

The repository already has a first-generation `replacement_flank_bp` mechanism:

- `ResolutionConfig::replacement_flank_bp` defaults to `0`.
- `PathRange` carries `sequence`, `extended_sequence`, `left_flank_bp`, and
  `right_flank_bp`.
- `materialize_flanked_sequences` walks outward from `begin_step` and
  `end_step`.
- `clip_replacement_to_interior` removes leading and trailing flank bases after
  replacement construction.
- `apply_replacement_frontier` remains the lacing step and
  `path_sequences_equal` remains the full-graph preservation gate.

The implementation task should preserve the useful shape but promote it into a
single explicit extraction/resolution/trimming/lacing contract shared by every
resolver. The important missing semantics are occurrence identity for repeated
contexts, target-orientation metadata, resolver-independent trim mapping,
diagnostic accounting, and motif-local/nested-region scope handling.

## Terms

Target interval:
The exact per-path step and bp interval that may be replaced. In the current
code this is each `PathRange { path_idx, begin_step, end_step,
source_begin_bp, source_end_bp }`. The target may be a simple POVU bubble, a
flubble, a nested site found inside a local replacement graph, a generated
multi-site window, or a motif-local anchor-pair window. The implementation must
not silently expand it when flanks are requested.

Occurrence:
One path-local traversal of the target interval. A repeated flank sequence in a
different position is a different occurrence and must not be substituted as
context.

Resolver input:
The sequence handed to SPOA, abPOA, POASTA, SweepGA/seqwish, POASTA block
polish, or a future graph-output resolver. With flanks enabled it is
`canonical_left_flank + canonical_target + canonical_right_flank`; otherwise it
is just `canonical_target`.

Trim plan:
Per-occurrence metadata that maps resolver-output path coordinates back to the
target interval. It is based on path coordinate offsets, not on searching for a
flank string.

Lacing sequence:
The replacement path spelling that must be inserted into the original path
orientation over the original target interval.

## Data Contract

Introduce an internal data object equivalent to this shape. Names can differ,
but the fields and invariants should be present.

```rust
struct FlankAwareCandidate {
    candidate_id: String,
    source_kind: CandidateSourceKind,
    target_orientation: TargetOrientation,
    replacement_scope: ReplacementScope,
    occurrences: Vec<TargetOccurrence>,
    stats: TraversalStats,
    lineage: CandidateLineage,
}

struct TargetOccurrence {
    path_idx: usize,
    source_path_name: String,
    original_path_name: String,
    target_begin_step: usize,
    target_end_step: usize,
    target_begin_bp: usize,
    target_end_bp: usize,
    source_orientation: OccurrenceOrientation,
    target_steps: Vec<Step>,
    target_sequence_path_orientation: Vec<u8>,
    target_sequence_canonical: Vec<u8>,
    left_flank: FlankContext,
    right_flank: FlankContext,
    resolver_sequence: Vec<u8>,
    trim_plan: TrimPlan,
}

struct FlankContext {
    requested_bp: usize,
    actual_bp_path_orientation: usize,
    actual_bp_canonical: usize,
    path_side: FlankSide,
    canonical_side: FlankSide,
    path_step_range: Option<(usize, usize)>,
    path_bp_range: Option<(usize, usize)>,
    truncation: FlankTruncationReason,
    sequence_path_orientation: Vec<u8>,
    sequence_canonical: Vec<u8>,
}

struct TrimPlan {
    canonical_left_trim_bp: usize,
    canonical_right_trim_bp: usize,
    expected_flanked_sequence_canonical: Vec<u8>,
    expected_target_sequence_path_orientation: Vec<u8>,
    restore_orientation: OccurrenceOrientation,
}
```

Required invariants:

- `target_begin_step..target_end_step` is the only laced interval.
- `target_sequence_path_orientation` equals the current graph path spelling over
  the target interval.
- `resolver_sequence` equals
  `left_flank.sequence_canonical + target_sequence_canonical +
  right_flank.sequence_canonical`.
- `trim_plan.canonical_left_trim_bp + target_sequence_canonical.len() +
  trim_plan.canonical_right_trim_bp == resolver_sequence.len()`.
- Path names are carried from the source graph. Resolver-local names may be
  sanitized for FASTA or external tools, but the mapping back to source path
  names must be one-to-one and lossless.

## Algorithm

### 1. Freeze The Target Interval

Candidate discovery must first emit a fully specified target interval for every
occurrence. This happens before any flank collection.

For POVU bubble/flubble/nested candidates:

- Use the existing discovered `PathRange` step interval for each path
  occurrence.
- Preserve POVU lineage: site id, parent id, level, leaf flag, and any local
  replacement parent signature used for nested descent.
- Preserve the root-path begin/end and root span for scheduling only. These are
  not a license to extend the replaceable interval.

For multi-bubble or complete-homologous generated windows:

- Store the generated merged step range per path after all child ranges have
  been merged.
- Store source ancestry: source site ids, window mode, source-site count, and
  any homologous-node expansion flag.
- The target interval is the merged per-path interval, not each child site.

For motif-local regions:

- Store the exact anchor-pair interval returned by motif-local discovery. If
  the current motif-local builder includes anchor steps in the target interval,
  that is the target. Flanks start outside that interval.
- Store motif source ancestry, such as sparse-offshoot path/core steps, repeat
  run, anchor labels, support counts, order jump, and white-space evidence.
- A motif-local candidate that cannot be represented as a deterministic
  per-path interval is not eligible for lacing and must be skipped with a
  diagnostic.

### 2. Collect Occurrence-Local Flanks

For each occurrence, collect up to `replacement_flank_bp` bases immediately
adjacent to the target interval on the same path occurrence.

Left flank:

- Walk path steps before `target_begin_step` toward the path start.
- Stop at the requested bp length, path start, graph boundary, containing
  replacement-scope boundary, or an explicit N-run boundary if N-cutting is
  configured to split contig context.
- If the last included step overshoots the requested length, trim bases at the
  far end so the flank remains adjacent to the target boundary.

Right flank:

- Walk path steps from `target_end_step` toward the path end.
- Stop for the same reasons as the left flank.
- If it overshoots, trim bases at the far end.

The collector must never locate flanks by sequence lookup. Repeated `ACGT...`
context at another graph position is irrelevant unless it is adjacent to this
same path occurrence. Store path coordinates and step ranges with every flank
so diagnostics can prove provenance.

Actual flank length is per occurrence and per side. It may be shorter than the
request because of:

- path or contig end;
- graph boundary;
- containing local-replacement scope boundary;
- N-cut boundary when configured;
- fragmented haplotype path;
- missing traversal on a path;
- inconsistent entry/exit that prevents a deterministic occurrence interval.

Default behavior accepts short, zero-length, and one-sided flanks. A path-end
target with only a right flank is still a valid resolver input.

### 3. Normalize Orientation

Choose one canonical target orientation per candidate:

- For POVU candidates, use the root/reference occurrence orientation from the
  POVU site start to end.
- For generated multi-site windows, use the selected root-path occurrence.
- For motif-local anchor-pair candidates, use the orientation of the anchor pair
  that created the candidate.

For every occurrence, determine whether its path spelling is forward or reverse
relative to the canonical target orientation. The current implementation often
only admits forward anchor-pair occurrences; the flank-aware design must make
this explicit so future reverse traversals do not silently corrupt paths.

Forward occurrence:

- canonical left flank = path-left flank;
- canonical target = path target sequence;
- canonical right flank = path-right flank;
- trim coordinates are `(left_len, right_len)`;
- lacing uses resolver path steps as emitted after trim.

Reverse occurrence:

- canonical left flank = reverse-complement(path-right flank);
- canonical target = reverse-complement(path target sequence);
- canonical right flank = reverse-complement(path-left flank);
- trim coordinates are `(path_right_len, path_left_len)`;
- after trimming, restore the lacing path to original path orientation by
  reversing the trimmed path walk and toggling step orientation, or by emitting
  an equivalent graph path whose `path_sequence` equals the original
  `target_sequence_path_orientation`.

The entire resolver input must be in canonical orientation. Do not feed a mix
of forward and reverse-complemented target sequences and expect the resolver to
infer strand. Strand inference may still be used inside a pairwise aligner, but
the wrapper must validate the output paths against the canonical resolver
inputs before trimming.

### 4. Resolve In Global Mode

Every resolver must consume the same `resolver_sequence` list and must emit a
graph/path representation with one path per input sequence, in original input
order after name mapping.

Resolver requirements:

| Resolver | Required global mode | Output contract |
|---|---|---|
| SPOA/direct POA | Use existing global/end-to-end SPOA engine. No local clipping by the resolver. | Replacement GFA paths spell the flanked input sequences exactly before trim. |
| abPOA | Invoke abPOA in global mode with GFA output. Reject empty traversals that would require sentinel bases. | Normalize abPOA path names back to source headers and validate exact flanked spellings. |
| POASTA | Use POASTA global alignment with the configured two-piece affine scoring. | Convert clipped W/P output into exact graph paths; output paths must spell flanked inputs before trim. |
| SweepGA/seqwish | Run local all-vs-all replacement alignment on the flanked inputs, then seqwish induction using the shared replacement filter settings. | The induced graph must contain exact walks for every flanked input. Zero PAF is a build failure/diagnostic skip, not a shape-quality decision. |
| POASTA block helpers | Use the same flanked input and trim contract when called as a crush replacement builder. | Do not keep separate flank logic in block-specific code. |
| Future graph-output resolver | Implement `build_replacement_from_flanked_inputs(inputs, config) -> ResolverGraphOutput`. | Must provide one exact path/walk per input and stable path names. Otherwise it is not eligible for lacing. |

Global/end-to-end mode is the default resolver mode for crush replacement.
Local/semi-global resolver modes may be added only behind an explicit different
method name because they change the trim contract.

### 5. Trim By Coordinates, Not By Flank String

After the resolver emits a graph, validate each resolver-output path before
trimming:

- path count equals occurrence count;
- path names can be mapped one-to-one back to input headers;
- each path spells `expected_flanked_sequence_canonical` exactly;
- no duplicate path name loses source identity;
- path walk references valid segments with known sequence.

Then trim each resolver-output path by its `TrimPlan`:

- `interior_start = canonical_left_trim_bp`;
- `interior_end = path_len - canonical_right_trim_bp`;
- `interior_end - interior_start` must equal
  `target_sequence_canonical.len()`;
- map these offsets through the resolver-output path walk;
- if a boundary lands inside a segment, split that oriented segment;
- drop segments that are used only by flank sequence;
- preserve graph sharing inside the trimmed target interval when the same slice
  appears on multiple paths.

Repeated flank sequence is not ambiguous because trimming is by path offset.
Ambiguity is only about the resolver-output graph/path mapping.

Hard trim failures:

- missing resolver path for an input;
- duplicate resolver path name without a lossless input mapping;
- resolver path spelling differs from the expected flanked sequence;
- resolver path length is shorter than recorded flanks;
- trim boundaries invert or do not leave the expected target length;
- a path offset cannot be mapped to a unique segment slice;
- reverse-orientation restoration cannot produce the original target spelling;
- the trimmed replacement has the wrong path count or wrong per-path target
  spelling.

These failures are path-corruption risks and must produce a diagnostic
skip/fail. They are not quality guards.

### 6. Restore Path Orientation And Names

The trimmed replacement graph must be converted from canonical orientation to
the lacing orientation expected by the original graph:

- For forward occurrences, keep the trimmed path walk orientation.
- For reverse occurrences, reverse the path walk and toggle step orientation, or
  otherwise emit a path that spells the original target sequence.
- Path names in the trimmed replacement remain resolver-local only; lacing uses
  path index and the original path name from `TargetOccurrence`.
- Replacement/source names written to debug artifacts must retain the semantic
  source token, including sample/haplotype/contig decorations and visible
  source coordinate ranges. Synthetic duplicate suffixes are allowed only as
  internal resolver ids and must be mapped back before graph lacing.

Before lacing, validate:

- `path_sequence(trimmed_replacement, replacement.paths[i]) ==
  occurrence.target_sequence_path_orientation` for every occurrence;
- replacement path count equals occurrence count;
- each occurrence still points to an existing path and the original
  `target_begin_step..target_end_step` range in the current graph.

### 7. Lace Only The Original Target Interval

Lacing must reuse the existing `apply_replacement_frontier` semantics:

- Replace `target_begin_step..target_end_step` on each source path.
- Insert only the trimmed replacement path walk.
- Do not insert left or right flank bases/nodes.
- Preserve all original path names exactly.
- Preserve full path spellings exactly.
- Keep non-overlap checks based on original target intervals, not flanked
  intervals.

The post-lace hard gate is exact full-graph path preservation by spelling:

```text
path_sequence_map(before_graph) == path_sequence_map(after_graph)
```

Graph-shape proxies are diagnostics only. Segment count, link count, singleton
bp, path-step count, coverage, trivial-stringy counts, white-space p99, and
compression ratio must not reject or roll back a path-valid replacement.

## CLI And Config Defaults

Keep the existing user-facing flank knob:

- `--replacement-flank-bp <N>` on `impg crush`;
- aliases `--flank-bp` and `--flank`;
- engine-stage aliases `replacement-flank-bp=<N>`, `replacement-flank=<N>`,
  `flank-bp=<N>`, and `flank=<N>`;
- parse size suffixes through the existing size parser;
- default `0`, preserving current behavior unless the user opts in.

Recommended initial C4 experiment value is `500` bp per side, matching the
earlier wider-context experiment. This is not a new default until validation
proves it should become one.

Minimum usable flank behavior:

- default minimum is `0`;
- zero-length and one-sided actual flanks are allowed;
- requested length and actual length are recorded independently per occurrence
  and side;
- candidates are not skipped merely because path ends, N boundaries, fragmented
  haplotypes, or missing adjacent context shorten the flanks;
- an optional future `replacement-min-flank-bp` must be documented as an
  explicit extraction policy, not as a graph-quality guard. If added, it should
  fail before resolver invocation with a clear diagnostic such as
  `flank-too-short`, never after inspecting output graph shape.

Diagnostics to add:

- requested flank bp;
- per-candidate min/median/max actual left/right flank bp;
- counts of zero-left, zero-right, and zero-both occurrences;
- truncation reasons by side;
- reverse-orientation occurrence count;
- resolver input bp including flanks versus target bp;
- trim failures by hard-failure reason;
- candidate source kind and lineage.

Resolver global-mode selection:

- global/end-to-end is the crush replacement default for SPOA, abPOA, POASTA,
  SweepGA/seqwish induction, and generic graph-output resolvers;
- no new CLI switch is needed for normal flank-aware crush;
- a non-global mode would need a separate method name or explicit experimental
  parameter and a different trim contract.

## Edge-Case Semantics

Repeated contexts:

- Store flank provenance as path index, step range, and bp range.
- Never find or validate flanks by searching for their sequence in the graph.
- Repeated identical flanks are safe because trim boundaries are path offsets
  inside each resolver-output path.

Requested flank longer than available context:

- Use the available adjacent context.
- Record actual length and truncation reason.
- Allow zero-length and one-sided flanks by default.
- Do not compensate by borrowing sequence from another path or another repeat
  occurrence.

Nested bubbles and flubbles:

- Target interval is the nested candidate range in the current working graph.
- Flanks may be collected from adjacent context in that same current graph path.
- If the nested candidate has a known parent/local replacement scope, the flank
  collector must not leave that scope unless the implementation explicitly
  records that whole-path flanking is intended for that candidate source.
- Child discovery and non-overlap scheduling continue to use target intervals,
  not flanked intervals.

Motif-local regions:

- The anchor-pair window is the target interval if that is what the candidate
  builder emits.
- Flanks are outside the anchor-pair window.
- Sparse/offshoot support metrics decide discovery only. They must not become
  replacement acceptance guards.
- Motif-local candidates with inconsistent path entry/exit are skipped before
  resolver invocation with a diagnostic.

Ns and fragmented haplotypes:

- N-cut policies may stop flank collection at configured N-run boundaries.
- Ns inside the target interval remain part of the target and must be preserved
  by path spelling.
- Missing traversals or paths that do not enter and exit consistently are not
  resolver inputs for that candidate. If fewer than two deterministic
  occurrences remain, skip the candidate.

Reverse orientation:

- Store orientation metadata before sequence extraction.
- Reverse-complement the full flanked occurrence into canonical orientation for
  resolver input.
- Swap left/right trim lengths under reverse orientation.
- Restore path orientation before lacing and validate per-path target spelling.

Path-name stability:

- The graph after lacing must contain exactly the same path-name set as before.
- Path names include PanSN/sample/haplotype/contig decorations and any visible
  coordinate ranges.
- Resolver-local FASTA ids may be sanitized, deduplicated, or coordinate-trimmed
  only if there is a lossless map back to source path names.
- Debug replacement paths should prefer semantic source names and must not leak
  `__impg`, `local_`, or duplicate-copy names into the final graph.

## Implementation Plan

1. Add a shared flank-aware extraction module inside `resolution.rs` or a
   nearby private module. The module owns target interval freezing, flank
   collection, orientation normalization, and `TrimPlan` construction.
2. Replace direct use of `PathRange.extended_sequence` with a shared
   `FlankedResolverInput` builder. Existing fields may remain temporarily, but
   resolver builders should depend on the new object.
3. Make `build_replacement_with_method_report` call one shared wrapper:
   `build_flanked_replacement(candidate, config, method)`.
4. Inside the wrapper:
   - build resolver inputs;
   - call the method-specific builder with those inputs;
   - validate exact flanked output paths;
   - trim through `TrimPlan`;
   - restore orientation;
   - validate target spellings;
   - return an ordinary trimmed replacement graph to lacing.
5. Convert SPOA, abPOA, POASTA, SweepGA/seqwish, AllWave/seqwish,
   WFMash/seqwish, chain/POASTA helpers, and motif-local paths to use the same
   wrapper.
6. Keep `apply_replacement_frontier` operating on target intervals only.
7. Extend logs and debug artifacts with flank/trim diagnostics.
8. Add tests before changing defaults.

Do not duplicate flank extraction in resolver-specific code. Resolver builders
should not know how to walk graph context; they should only know how to turn
named flanked sequences into a graph with exact input paths.

## Test Matrix

Synthetic fixtures should live under `tests/test_data/crush/flank_aware/` or be
generated inline next to existing `resolution.rs` unit tests when smaller.

| Fixture | Purpose | Expected assertions |
|---|---|---|
| `flank_required_indel_alignment.gfa` | Two or more traversals with a local indel/homopolymer where interior-only global alignment can split homologous sequence but flanks anchor it. | With `replacement_flank_bp > 0`, resolver input includes flanks; trimmed replacement excludes flanks; all path spellings are preserved; diagnostic target bp excludes flank bp. |
| `repeated_context_microtangle.gfa` | Same flank sequence appears at multiple graph/path positions around similar motif copies. | Flank provenance points to the adjacent occurrence; only the chosen target interval is laced; an identical flank elsewhere is unchanged; no sequence-search ambiguity. |
| `path_end_one_sided_flank.gfa` | Candidate starts at path/contig beginning or ends at path/contig end. | Actual flank lengths record `0` on the missing side; resolver still runs by default; lacing preserves paths. |
| `reverse_orientation_traversal.gfa` | One occurrence traverses the target in reverse orientation or through reverse-complemented W/P steps. | Resolver inputs are canonical; reverse occurrence has swapped trim lengths; post-trim path is restored to original orientation; path names/spellings unchanged. |
| `trim_boundary_ambiguity_fake_resolver.gfa` | Fake resolver emits missing path, duplicate path, wrong flanked spelling, or path length shorter than recorded flanks. | Candidate is diagnostic `path-invalid`/hard skip; no partial lacing occurs; error message names the trim failure. |
| `path_name_stability_pansn.gfa` | Paths have sample/haplotype/contig decorations, coordinate ranges, whitespace tokens, and local replacement source names. | Final graph path names exactly equal input names; debug names remain semantic; no `__impg`, `local_`, or duplicate-copy token leaks into final paths. |
| `motif_local_anchor_window.gfa` | Sparse motif offshoot or self-loop run with high-support flanking anchors. | Target interval is the motif-local anchor-pair range; flanks are outside it; replacement laces only the target; support metrics are discovery diagnostics only. |
| `fragmented_or_missing_traversals.gfa` | Some paths enter without a matching exit, contain N-boundaries, or have fragmented haplotype context. | Inconsistent paths are excluded from that candidate; actual flank truncation reasons are logged; candidate skips only if fewer than two deterministic occurrences remain. |

Focused unit tests:

- `flank_extraction_records_requested_and_actual_lengths`
- `flank_extraction_uses_occurrence_coordinates_not_sequence_search`
- `orientation_normalization_swaps_reverse_flanks`
- `trim_plan_splits_segments_at_flank_boundaries`
- `trim_plan_rejects_missing_or_misspelled_resolver_path`
- `lacing_uses_original_target_interval_not_flanked_interval`
- `graph_shape_metrics_are_diagnostic_not_acceptance_gates`

Integration tests:

- Existing `c4_slice_auto_crush_with_flank_preserves_path_sequences` should be
  strengthened to assert diagnostic flank counts and no flank nodes in the
  laced interval.
- Add method-specific small tests for `poa`, `poasta`, `abpoa` when binaries or
  fake resolvers are available.
- Add one SweepGA/seqwish smoke test guarded by existing external-tool
  availability conventions.

## C4 Measurement Plan

The C4 validation task should write a durable report under:

```text
docs/evaluations/flank-aware-crush-c4.md
```

Use an output directory like:

```text
/home/erikg/impg/data/flank_aware_crush_c4_YYYYMMDDTHHMMSSZ/
```

Required local artifacts:

- `commands.sh` with every command executed;
- `tool_versions.tsv`;
- input GFA path and checksum;
- one output GFA per variant, both unsorted and Ygs-sorted when sorted;
- `graph-report` markdown/stdout per output;
- `compare_gfa_paths` stdout/stderr per output;
- `candidate-accounting.tsv` with candidate source, method, requested flank,
  actual flank distribution, resolver status, trim status, and applied status;
- `flank-diagnostics.tsv` with per-candidate actual flank/truncation summaries;
- rendered PNG path per output;
- uploaded PNG URL per output.

Minimum command shape:

```bash
OUT=/home/erikg/impg/data/flank_aware_crush_c4_YYYYMMDDTHHMMSSZ
SEED=/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_minmatch1_scaffold0.initial.gfa
COMPARE=/home/erikg/impg/target/release/examples/compare_gfa_paths

mkdir -p "$OUT"/{graphs,reports,validation,renders,logs}

/usr/bin/time -v target/release/impg crush \
  --gfa "$SEED" \
  --output "$OUT/graphs/c4.flank500.nosort.gfa" \
  --method auto \
  --max-rounds until-done \
  --replacement-flank-bp 500 \
  --max-pair-alignments 0 \
  --max-replacement-paf-bytes 0 \
  --polish-rounds until-done \
  --polish-method poasta \
  --polish-max-traversal-len 10k \
  --polish-max-median-traversal-len 1k \
  -t 32 -v 1 \
  >"$OUT/logs/c4.flank500.stdout.log" \
  2>"$OUT/logs/c4.flank500.stderr.log"

gfasort -i "$OUT/graphs/c4.flank500.nosort.gfa" \
  -o "$OUT/graphs/c4.flank500.Ygs.gfa" -p Ygs -t 32

target/release/impg graph-report \
  -g "$OUT/graphs/c4.flank500.Ygs.gfa" \
  >"$OUT/reports/c4.flank500.graph-report.md"

"$COMPARE" "$SEED" "$OUT/graphs/c4.flank500.Ygs.gfa" \
  >"$OUT/validation/c4.flank500.compare_gfa_paths.stdout.log" \
  2>"$OUT/validation/c4.flank500.compare_gfa_paths.stderr.log"

gfalook -i "$OUT/graphs/c4.flank500.Ygs.gfa" \
  -o "$OUT/renders/c4.flank500.Ygs.png" -m -x 2200 -y 1200

scp "$OUT/renders/c4.flank500.Ygs.png" \
  hypervolu.me:/var/www/erik/impg/flank-aware-crush-c4-flank500.png
printf '%s\t%s\n' c4.flank500 \
  'https://hypervolu.me/~erik/impg/flank-aware-crush-c4-flank500.png' \
  >>"$OUT/renders/png_urls.tsv"
```

Also run at least:

- `flank0` control with the same resolver settings;
- `flank100`;
- `flank500`;
- `flank1000` if runtime allows;
- one synthetic/small-locus sanity run that exercises reverse or one-sided
  flanks.

Required report tables:

Candidate accounting:

| variant | selected | built | applied | hard_trim_fail | path_invalid | zero_paf | skipped_inconsistent | zero_left | zero_right | reverse_occurrences |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|

Path preservation:

| variant | input_paths | output_paths | spelling_mismatches | missing_paths | extra_paths | result | compare_stdout | compare_stderr |
|---|---:|---:|---:|---:|---:|---|---|---|

Graph/report metrics:

| variant | S | L | path_steps | segment_bp | singleton_bp | self_loops | adjacent_same_node_repeats | duplicate_seq_frac | trivial_stringy | ws_p99 | ws_max | wall | max_rss_kb | PNG |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|

The C4 report may recommend a flank length based on these measurements, but
these measurements must not be implemented as acceptance gates. A replacement
is allowed to make graph-shape metrics worse if exact path spelling is
preserved and trim/lacing are unambiguous.

## Acceptance Policy

Allowed hard blockers for one replacement:

- unmappable or ambiguous trim boundary;
- missing or duplicate resolver path without a lossless source mapping;
- resolver-output path spelling differs from expected flanked input;
- orientation metadata missing or inconsistent;
- fewer than two deterministic occurrences after excluding inconsistent paths;
- post-trim per-path target spelling differs from original target spelling;
- post-lace full-graph path-name or path-spelling corruption.

Rejected as replacement blockers:

- segment count regression;
- link count regression;
- path-step count regression;
- singleton-bp increase;
- duplicate-sequence fraction;
- trivial-stringy count;
- white-space p99/max;
- compression ratio;
- heuristic graph quality score;
- speculative visual quality.

Those graph-shape values belong in diagnostics, measurements, and C4 reports
only. They must not reject, roll back, or reroute a replacement that satisfies
the hard path-corruption checks.

## Lineage

This design is downstream of `quality-pass-flank-aware-crush` and incorporates
the canonical model recorded in `docs/flank-aware-crush-quality-pass.md`. It
also supersedes the underspecified parts of the earlier wider-context
experiment in `docs/crush-wider-context-bubbles.md`: the earlier implementation
proved the basic `add flanks, align globally, clip before lacing` shape, while
this document defines the stricter shared contract needed before expanding the
feature to motif-local and graph-output resolver paths.
