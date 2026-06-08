# Audit: SmoothXG/PGGB Mechanics for Localized Polishing

Task: `audit-smoothxg-pggb`
Date: 2026-06-08

## Summary

The repository already has most of the mechanics needed for localized
PGGB/SmoothXG-style polishing, but they are wired as whole-graph controls or
post-crush transforms, not as a first-class localized dirty-region driver.

The next implementation pass should not restart broad C4-only algorithm churn.
It should reuse the existing mechanics in this order:

1. Collect a region-local sequence set with the existing SYNG/query sequence
   extraction path.
2. Induce a region-level seed graph with the existing seqwish/SweepGA tail.
3. Run the existing SmoothXG-style path-overlap smoother on only selected dirty
   region seeds.
4. Lace the replacement back through the existing path-preserving replacement
   and flank-trimming invariants.
5. Run one final normalization/sort pass.

The main missing work is adapter code: region/dirty selection, local seed object
construction, a public replacement/lacing API for non-POVU dirty regions, and
guarded integration with shared graph-output paths.

## Invocation Surfaces Inspected

- `scripts/local_compression_testbed.py` defines mandatory local methods and
  optional PGGB/SmoothXG controls (`METHODS`, lines 93-178). The local methods
  include raw SYNG, flubble/local crush variants, sorted chunk-window rows, and
  whole-region SweepGA/seqwish rows (lines 45-55, 94-156). Optional controls are
  `pggb_control`, `smoothxg_control`, and `pggb_plus_smoothxg_control` (lines
  58-70, 157-177).
- The controls are internal IMPG invocations, not standalone PGGB/SmoothXG
  binaries. `control_availability` records `engine_invocation=impg graph
  --gfa-engine pggb` and marks standalone `pggb`/`smoothxg` discovery as
  diagnostic only (`scripts/local_compression_testbed.py:1565-1620`).
- `control_command` invokes:

```bash
impg graph --sequence-files <fixture.fa> --gfa-engine pggb --threads 1 \
  --min-match-len 1 --min-map-length 1 \
  --target-poa-length <64 or 32,64> --max-node-length 64 \
  --poa-padding-fraction 0 -g <output.gfa>
```

  For `pggb_control` and `pggb_plus_smoothxg_control`, the runner inserts
  `--fastga --scaffold-jump 0 --scaffold-mass 1`; for `smoothxg_control`, it
  inserts an explicit empty PAF before `--gfa-engine` so the pggb engine's
  smoothing stage runs over a source-sequence graph
  (`scripts/local_compression_testbed.py:1623-1670`).
- `tests/test_local_compression_testbed.rs` asserts that every `*_control`
  row uses `impg graph`, uses `--gfa-engine pggb`, does not invoke standalone
  `smoothxg`, and preserves exact path spellings on graph-producing rows
  (lines 252-287).
- The fast-profile artifact index confirms 13 fixtures, 12 methods, internal
  control execution, exact-path hard gating, and row-visible topology failures
  rather than hidden rejection (`docs/evaluations/local-compression-testbed-fast/report.md:25-55`).

Important limitation: the local chunk-window rows in the testbed are synthetic
GFA generators, not real SmoothXG/seqwish invocations. `sorted_chunk_window_gfa`
and `chunk_window_sweepga_seqwish_gfa` split equal-length expected spellings by
variant runs and shared invariant gaps (`scripts/local_compression_testbed.py:936-1129`).
They are useful fixtures for expected topology, but they are not adapter code for
real local graph polishing.

## Current Mechanics and Parameters

| Area | Current implementation | Parameters currently used |
| --- | --- | --- |
| Whole-graph PGGB path | `src/commands/graph.rs::run_graph_build_pggb` runs seqwish induction, sorts for 1D layout, calls `smooth::smooth_gfa`, then `graph::normalize_and_sort` and restores direct FASTA path names (lines 606-676). `src/lib.rs` mirrors this for query/partition `GfaEngine::Pggb` (lines 1013-1064). | CLI defaults: target POA lengths `700,1100`, max node length `100`, padding fraction `0.001` (`src/main.rs:2085-2101`). Testbed controls override target POA to `64` or `32,64`, max node `64`, padding `0` (`scripts/local_compression_testbed.py:65-71`, `1633-1655`). |
| Seqwish induction | `src/commands/graph.rs::induce_graph_from_alignment` builds the seqwish index, indexes PAF, computes transitive closures, compacts, derives links, emits GFA, and unchops (`src/commands/graph.rs:173-445`). `src/syng_graph.rs::build_gfa_from_paf_and_sequences` writes local FASTA/PAF, filters PAF, optionally lowers adaptive min-match, then reuses `induce_graph_from_alignment` (`src/syng_graph.rs:955-1088`). | `GraphBuildConfig` defaults: `min_match_len=23`, `repeat_max=0`, `min_repeat_dist=0`, `sparse_factor=0.0`, `transclose_batch=10_000_000`, `disk_backed=false`, `aligner=wfmash`, `num_mappings=many:many`, `scaffold_jump=50000`, `scaffold_mass=10000`, `scaffold_filter=many:many`, `min_map_length=0`, `map_pct_identity=90` (`src/commands/graph.rs:43-145`). Testbed controls use `--min-match-len 1 --min-map-length 1`. Local replacement tail normally resolves seqwish min to `1` via `replacement_min_match_len_policy=Fixed(1)` (`src/resolution.rs:570-613`, `9097-9155`). |
| SmoothXG-style smoothing | `src/smooth.rs` implements smoothxg-style block decomposition and POA smoothing over GFA strings. `SmoothConfig::new` defaults to two passes `700,1100`, `max_node_length=100`, `poa_padding_fraction=0.001`, scoring `(1,4,6,2,26,1)`, `block_source=PathOverlap` (`src/smooth.rs:19-78`). `smooth_gfa` runs one pass per target length, re-decomposing each pass (`src/smooth.rs:199-254`). | Per pass: `max_block_weight = target_poa_length * n_haps`, `max_poa_length = 2 * target_poa_length`; unchop+sort unless `pre_sorted`; chop; block placement; long-block break; per-block SPOA/POASTA; lace (`src/smooth.rs:299-508`). Padding has a 311 bp base for blocks at depth <= 100 when padding is enabled (`src/smooth.rs:2103-2120`). |
| Smooth block placement | The default `SmoothBlockSource::PathOverlap` is the closest SmoothXG-style behavior. It walks sorted node order, accumulates unseen path steps, finalizes when block weight or max path length exceeds the target, sorts path ranges longest-first, then topologically splits disconnected components (`src/smooth.rs:47-56`, `715-984`). | This differs from flubble/bubble-only methods because it is driven by sorted graph layout and path overlap, not only POVU bubble boundaries. It can span adjacent bubbles and decompose non-bubble dirty graph regions. |
| Flubble and neighbor placement | `SmoothBlockSource::Flubble` uses POVU flubble extents and gap passthrough; `NeighborMergePoasta` greedily merges reference-path neighbor sites and uses POASTA (`src/smooth.rs:348-401`, `1036-1410`). Parser aliases include `smooth`, `smoothxg`, `flubble-smooth`, and `neighbor-merge-poasta` (`src/main.rs:3213-3382`, `3873-3888`). | Neighbor-merge defaults to `target_poa_lengths=[10000,10000,10000]` when parsed as `neighbor-merge-poasta` (`src/main.rs:3223-3229`). Chain-POVU replacement uses one 10 kbp smoothing pass, `poa_padding_fraction=0.0`, then two POASTA cleanup iterations (`src/resolution.rs:521-522`, `9332-9388`). |
| Local replacement and polishing | `src/resolution.rs` implements path-preserving replacement over POVU/flubble candidates. It has direct POA, POASTA, allwave/seqwish, SweepGA/seqwish, wfmash/seqwish, top-flubble-sweepga, iterative multi-level, and motif-local variants (`src/resolution.rs:272-337`). `ResolutionPolishMethod::Smooth` runs the SmoothXG-style smoother on a replacement graph (`src/resolution.rs:406-424`, `10637-10710`). | Defaults route by median traversal length: `<1 kb` -> SPOA, `<10 kb` -> POASTA, otherwise SweepGA/seqwish (`src/resolution.rs:6461-6496`). Pair defaults: k-nearest `3`, k-farthest `1`, tree count `1`, random fraction `0.01`, mash-k `15` (`src/resolution.rs:530-535`, `592-596`). SweepGA replacement uses `fastga`, auto kmer frequency `max(1000, 10 * traversals)`, min aln `0`, max pair alignments `10000`, max PAF bytes `64 MiB` (`src/resolution.rs:542-547`, `9158-9183`). |
| Flank gathering | `replacement_flank_bp` is implemented but defaults to `0`; the canonical wider-context run uses `500` (`src/resolution.rs:554-557`). Flanks are materialized as aligner-visible context only, normalized for orientation, then clipped from the replacement before substitution (`src/resolution.rs:8076-8338`, `10746-10822`). | `tests/test_crush_integration.rs::c4_slice_auto_crush_with_flank_preserves_path_sequences` sets `replacement_flank_bp=500` and checks every path sequence is unchanged (`tests/test_crush_integration.rs:1215-1258`). |
| Normalization and sorting | Shared final step is `graph::normalize_and_sort`: run sibling `gfaffix`, normalize repeat self-loop runs, then `sort_gfa` with pipeline `Ygs` (`src/graph.rs:858-895`). `sort_gfa_pipeline` supports `Y`, `g`, `s`, `S`, and `u`; `Ygs` uses gfasort's combined YGS path (`src/graph.rs:897-962`). | `DEFAULT_SYNG_GFA_SORT_PIPELINE` is `Ygs` (`src/lib.rs:65`). Post-crush smooth runs `gfaffix` after smoothing and leaves final self-loop normalization/sort to the optional graph sort stage (`src/lib.rs:855-920`). |
| SYNG local sequence collection | `graph::prepare_sequences` fetches query intervals, strand-normalizes them, and creates coordinate path names (`src/graph.rs:14-58`, `485-560`). `build_syng_local_region_gfa_from_intervals` collects those sequences and calls `write_syng_region_gfa_from_sequences_with_params` (`src/lib.rs:722-752`). | `write_syng_region_gfa_from_sequences_with_params` builds a fresh regional SYNG GBWT, temporary FASTA/index, and writes the selected ranges as a local SYNG GFA (`src/lib.rs:598-705`). Local SYNG default frequency mask: top fraction `0.0005`, min shared run `5`, high-frequency run `10`, high-frequency sequence span `1000`, local repeat max minor `2`, dominant fraction `0.80`, N-cut default min run `1` when enabled (`src/commands/syng2gfa.rs:24-61`, `70-118`; `src/syng.rs:165-166`). |

## SmoothXG-Style vs Flubble/Bubble-Only

Flubble/bubble-only approaches depend on existing graph decomposition boundaries:
POVU flubbles, top-level roots, or exact local bubbles. This is visible in
`ResolutionMethod::TopFlubbleSweepga`, which selects level-0/root flubbles and
runs local SweepGA/seqwish for each complete bounded region without merging
adjacent roots (`src/resolution.rs:313-317`, `2182-2406`).

SmoothXG-style path-overlap smoothing is different:

- It first needs a stable 1D graph order. The PGGB path sorts the raw seqwish
  GFA before smoothing (`src/commands/graph.rs:642-663`; `src/lib.rs:1022-1035`).
- It places blocks by sorted graph order, path overlap, and target path weight,
  then splits disconnected components (`src/smooth.rs:715-984`).
- It runs multiple smoothing passes, re-decomposing after each pass, so later
  passes can address longer variation exposed by earlier passes (`src/smooth.rs:199-254`).
- It retains path coverage explicitly. Gap blocks are passthrough by design,
  smoothable blocks must produce a non-empty path-preserving GFA, and padded
  trimming falls back to no-padding POA if it would shrink path-coordinate
  coverage (`src/smooth.rs:425-508`, `1857-1937`, `2257-2290`).

That makes SmoothXG-style smoothing the right mechanism for localized dirty
regions that are not clean flubbles, such as adjacent bubbles, residual stringy
regions, dispersed repeat glue, and local graph-layout badness. Flubble-only
methods remain useful as candidate discovery and bounded replacement paths, but
they should not be the only localized polishing strategy.

## Mapping to the Target Pipeline

### 1. SYNG-collected local sequence sets

Reuse:

- `graph::prepare_sequences` for query interval sequence extraction and canonical
  path naming (`src/graph.rs:485-560`).
- `write_syng_region_gfa_from_sequences_with_params` when the desired seed is a
  local SYNG graph over the exact sequence set (`src/lib.rs:598-705`).
- `SyngGfaFrequencyMask::local_default` for local repeat-aware syncmer sharing
  unless a task explicitly wants raw unmasked SYNG (`src/commands/syng2gfa.rs:104-118`).

Adapter needed:

- A first-class `LocalSequenceSet`/region object that carries source intervals,
  path names, sequence bytes, orientation, optional flank extents, and a stable
  mapping back to the parent graph.
- Dirty-region selection that feeds local sequence sets from graph badness or
  dirty-region evidence, not only POVU flubble candidates.

### 2. Region-level seed induction

Reuse:

- `syng_graph::build_gfa_from_paf_and_sequences` as the seed-induction tail
  when a local aligner or SYNG-anchor PAF is available (`src/syng_graph.rs:955-1088`).
- `commands::graph::induce_graph_from_alignment` for the actual seqwish graph
  induction path (`src/commands/graph.rs:173-445`).
- The local replacement seqwish config pattern: min-match floor effectively `1`,
  optional adaptive lowering, and no hidden scaffold-mass floor for local
  replacements (`src/resolution.rs:9097-9155`).

Adapter needed:

- A region-level seed induction driver that chooses between local SYNG seed,
  SweepGA/seqwish seed, or anchor-seeded PAF seed for the same sequence set.
- A clear policy for local PAF filtering. The next pass should start with the
  local replacement semantics (`min_match_len=1`, local no-filter/many-to-many
  when explicitly requested) rather than whole-graph PGGB defaults.

### 3. Localized dirty-region smoothing

Reuse:

- `smooth::smooth_gfa` with `SmoothBlockSource::PathOverlap` as the first
  SmoothXG-style smoother over region-level seed GFAs (`src/smooth.rs:199-254`,
  `299-508`, `715-984`).
- `SmoothPipelineConfig`/`:smooth` parsing for command-line exposure after the
  local adapter exists (`src/lib.rs:67-90`; `src/main.rs:3213-3382`).
- The `ResolutionPolishMethod::Smooth` path when smoothing is a replacement
  polish stage rather than a primary seed smoother (`src/resolution.rs:406-424`,
  `10687-10710`).

Adapter needed:

- A localized smoother entry point that takes a dirty-region object or
  region-level seed GFA and returns a replacement graph plus path-range mapping.
  Today `smooth_gfa` accepts a whole GFA string and relies on path-coordinate
  names for lacing; it does not directly expose a dirty-region replacement plan.
- Tests that compare exact path spellings before and after localized smoothing
  on non-C4 local fixtures, because the current fast-profile synthetic chunk
  rows do not prove the real smoother has run.

### 4. Replacement and final graph output

Reuse:

- Flank-context mechanics and clipping from `materialize_flanked_sequences`,
  `normalize_range_for_resolution`, `clip_replacement_to_interior`, and
  `validate_interior_replacement_paths` (`src/resolution.rs:8076-8338`,
  `10746-10822`).
- Full path spelling validation via `path_sequences`/`path_sequence_map`
  (`src/resolution.rs:6254-6260`, `11281-11290`).
- `graph::normalize_and_sort` for final normalization/sort once the localized
  replacements have been laced (`src/graph.rs:858-895`).

Adapter needed:

- A safe public replacement API, or a new internal module, so dirty-region
  smoothing can lace region replacements without duplicating private
  `resolution.rs` internals.
- A guard that applies normalization/sort once at the end of a batch of localized
  replacements, not after every tiny dirty region.

## Implemented Now vs Adapter Work

| Mechanic | Status |
| --- | --- |
| Internal PGGB control invocation through `impg graph --gfa-engine pggb` | Implemented and tested in the local compression runner. Reuse as control evidence, not as the production local adapter. |
| Whole-graph PGGB seqwish -> sort -> smooth -> gfaffix/sort path | Implemented in `run_graph_build_pggb` and query/partition dispatch. |
| SmoothXG-style path-overlap block decomposition, multi-pass smoothing, long-block breaking, block lacing | Implemented in `src/smooth.rs`. Reuse directly for region seed GFAs. |
| Flubble-guided and neighbor-merge smoothing variants | Implemented. Keep as secondary block sources or diagnostics, not as the primary localized SmoothXG substitute. |
| Local replacement, path validation, flank context, clipping, and exact-spelling guards | Implemented in `src/resolution.rs`, but many useful functions are private and coupled to POVU candidates. |
| SYNG local sequence extraction and local SYNG GFA rendering | Implemented through `graph::prepare_sequences`, `build_syng_local_region_gfa_from_intervals`, and `write_syng_region_gfa_from_sequences_with_params`. |
| Region-level explicit seed induction from a SYNG-collected sequence set | Partially implemented as reusable pieces. Needs a driver that owns the sequence set and chooses local SYNG, SweepGA/seqwish, or anchor-seeded PAF induction. |
| Dirty-region detection independent of flubbles | Not implemented as the target adapter. Current graph-quality reports diagnose badness, and iterative/motif replacement policies generate residual windows, but there is no clean localized SmoothXG dirty-region API yet. |
| Localized SmoothXG replacement plan and parent-graph lacing | Not implemented. This is the highest-risk adapter because it touches replacement, graph output, and exact path preservation. |
| Testbed real SmoothXG chunk execution | Not implemented. Current chunk-window testbed rows are synthetic expected-topology generators. |

## Recommended Next Implementation Pass

1. Implement an explicit local seed induction driver, not a new C4 search loop.
   File scope should stay close to `src/lib.rs`, `src/graph.rs`,
   `src/syng_graph.rs`, and a small new adapter module. Inputs should be
   SYNG/query-collected local sequence sets; outputs should be seed GFA plus
   path-name/source-range metadata.

2. Start seed induction with existing seqwish tail semantics. Use
   `build_gfa_from_paf_and_sequences` for local PAF -> seqwish GFA, and
   `write_syng_region_gfa_from_sequences_with_params` for local SYNG seed GFA.
   Do not call the whole `impg graph --gfa-engine pggb` CLI from inside the
   adapter.

3. Apply `smooth::smooth_gfa` to region-level seed GFAs with
   `SmoothBlockSource::PathOverlap`. Begin with pggb defaults
   `target_poa_lengths=[700,1100]`, `max_node_length=100`, padding `0.001` for
   realistic regional data; keep the testbed's `64` and `32,64` controls only
   as bounded fixture controls.

4. Add dirty-region lacing through existing path-preserving replacement
   mechanics. Prefer exposing a narrow helper around existing validation,
   flank-trim, and replacement application rather than duplicating lacing logic.

5. Normalize after the localized batch with `gfaffix`, self-loop normalization,
   and `Ygs` sort. Avoid normalizing each dirty-region subgraph unless the
   smoother or lacer requires it for correctness.

## Implementation-Order Risks

- Shared graph-output dispatch is already dense. `dispatch_gfa_engine`,
  `apply_graph_transforms`, partitioned lacing, `:crush`, `:smooth`, and
  `:sort/:nosort` interact in `src/lib.rs:830-920` and `1201-1332`. Add the
  local adapter behind a narrow stage or helper before changing CLI defaults.
- Replacement application is path-preservation-critical. `resolution.rs`
  currently validates by exact path spelling after replacement and flank
  trimming. Any new localized smoothing plan must preserve this invariant before
  it can touch graph-output code.
- The current `smooth_gfa` API consumes and emits whole GFA strings. If the next
  task tries to use it directly on the parent graph, it risks smoothing unrelated
  clean regions. Build the dirty-region seed GFA first, smooth locally, then
  lace.
- Sorting/normalization ordering matters. PGGB sorts before smoothing for block
  layout and normalizes after smoothing. SYNG/crush graph-output can also sort
  through `graph_sort_pipeline`. A double sort is mostly cosmetic, but a double
  normalize/lace sequence can obscure path-name or coordinate regressions.
- Flank context is implemented but off by default. Enabling it globally can
  increase replacement graph size and must be opt-in or dirty-region scoped until
  duplicated-flank and reverse-orientation fixtures cover the new adapter.
- The local compression chunk rows should not be treated as proof of real
  SmoothXG execution. They show what topology the target can represent; they do
  not exercise `src/smooth.rs`.
- C4 remains diagnostic. The next implementation should include non-C4 local
  compression fixtures and exact path checks before comparing C4 metrics against
  the PGGB/SmoothXG-style control.

## Validation Notes

This audit only added this report. It intentionally did not add generated GFA,
GFA.gz, GFA.zst, PNG, PDF, stdout/stderr logs, or per-fixture command artifacts.
