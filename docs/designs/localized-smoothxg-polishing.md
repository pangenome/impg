# Localized SmoothXG-Like Polishing For SYNG Local Graphs

Date: 2026-06-08
Task: `design-localized-smoothxg`
Status: design proposal over existing implementation pieces; integrated
production caller available as `gfa:syng-local:localized`

## Scope

This design defines a localized graph-polishing architecture for local haplotype
sequence sets collected by SYNG. The intended flow is:

```text
SYNG query / selected FASTA records / optional input GFA
  -> collect exact local path sequences
  -> build a seed graph for the whole local region
  -> detect dirty subregions in that seed graph
  -> smooth or rebuild only those subregions
  -> accept only replacements that preserve path spellings exactly
  -> emit polished GFA and diagnostic report
```

The goal is not another broad C4-only parameter search. C4 remains a required
stress case, but the implementation must first be driven by small synthetic
fixtures and at least a small non-C4 SVR panel when local data is present.

Graph-quality metrics are diagnostic and prioritization signals. The only hard
acceptance gate for replacing graph content is exact path corruption: a
replacement that changes any required path spelling or final path name is
rejected or left unchanged.

## Current Implementation

The integrated production path is `--gfa-engine syng-local:localized` (or the
compact output form `-o gfa:syng-local:localized,...`). It keeps the existing
query/partition graph-output conventions: SYNG collects the local sequence set,
plain `syng-local` builds the explicit whole-region SweepGA/FastGA + seqwish
seed graph, and `src/localized_polish.rs` iterates dirty-region detection plus
localized resolver tiers over selected flanked chunks.

The defaults are conservative and reproducible: three iterations, one chunk per
iteration, sixteen total chunks, a 500 kbp total chunk budget, detector defaults
from `DirtyRegionOptions`, and resolver defaults from `LocalizedResolverConfig`.
CLI parameters on the `:localized` stage expose iteration count, merge distance,
flank length, detector thresholds, resolver method/thresholds, chunk budgets,
wall-clock budget, and debug report directory.

Path preservation is the hard gate. The loop validates seed paths before work,
rejects any resolver `PathInvalid` report, revalidates after each applied batch,
and revalidates the final graph against the original local sequence set. Metrics
and dirty-region reports are logged and optionally written under the configured
debug directory; they are diagnostic only and never acceptance gates.

## Implemented Pieces

The repository already has most of the primitives needed for this design:

| Area | Existing references | Notes |
| --- | --- | --- |
| SYNG index, query, sidecars, and path naming | `src/syng.rs`, `docs/designs/syng-integration.md` | Provides `.1gbwt`, `.1khash`, `.names`, `.pstep`, `.spos`, `.meta`, path-position lookup, regional query, and sequence-name metadata. |
| SYNG GFA extraction | `src/commands/syng2gfa.rs`, `src/lib.rs` (`build_syng_region_gfa_from_intervals`, `build_syng_local_region_gfa_from_intervals`) | Supports whole-index and selected-region GFA output, blunt/native modes, local SYNG rebuild, and graph transform application. |
| Pipeline parsing | `src/graph_pipeline.rs`, `src/main.rs` (`parse_crush_stage`, `parse_smooth_stage`, `GraphEngineCli`) | Existing `gfa:<engine>:crush:smooth:nosort` syntax can host a new localized polishing stage. |
| Graph induction | `src/syng_graph.rs` (`build_gfa_syng_native_from_sequences`, `build_gfa_from_paf_and_sequences`), `src/commands/graph.rs` | Provides internal seqwish induction from sequence sets and PAF, plus PGGB-like graph build/smoothing through the graph command path. |
| SmoothXG-like smoothing | `src/smooth.rs` (`SmoothConfig`, `SmoothBlockSource`, `smooth_gfa`) | Current smoothing operates at whole-graph or existing block-source level. This design proposes dirty-region driven blocks. |
| Local replacement and path validation | `src/resolution.rs` (`ResolutionConfig`, `resolve_gfa_bubbles`, `build_*_replacement`, `validate_replacement_paths`, `render_rewritten_graph`) | Existing crush replacement machinery already validates replacement path spellings and has SPOA/abPOA/POASTA/SweepGA/seqwish builders. |
| Graph diagnostics | `src/graph_report.rs`, `src/gfa_self_loops.rs`, `examples/compare_gfa_paths.rs`, `examples/validate_gfa_path_sources.rs` | Existing reports include singleton bp, path white-space, self-loop diagnostics, POVU embedding, TSV/JSON/markdown output, and exact GFA path comparison helpers. |
| Synthetic local compression fixtures | `tests/test_data/local_compression/`, `scripts/local_compression_testbed.py`, `tests/test_local_compression_testbed.rs` | Provides fixture classes, expected path spellings, local methods, PGGB/SmoothXG-style controls, and diagnostic scoreboards. |

## Proposed Pieces

Add a localized polishing layer instead of broadening the existing C4-only crush
experiments. The layer should be implemented as either:

- a graph pipeline stage, for example `gfa:syng-local:<seed>:local-polish`,
  parsed beside `:crush` and `:smooth`; or
- a standalone internal entry point invoked by `query`, `partition`, and
  `graph` when the engine options request localized polishing.

The concrete API should be small and reusable:

| Proposed type | Likely owner | Purpose |
| --- | --- | --- |
| `LocalizedPolishConfig` | new `src/local_polish.rs` or `src/resolution.rs` if reuse is simpler | Seed choice, detector thresholds, chunk merge/flank/budget limits, resolver tier thresholds, iteration limits, report settings. |
| `SeedGraphSource` | same | `SyngDerived`, `WholeRegionSweepgaSeqwish`, `PggbLike`, `ExistingGfa`, and `ChunkSweepgaSeqwish`. |
| `DirtySite` | same or `src/graph_report.rs` if shared | A detector hit with metric source, graph node/path coordinates, severity, and explanatory fields. |
| `DirtyChunk` | same | Merged dirty sites plus flanks, selected path slices, source names, sequence spellings, and budget accounting. |
| `ResolverTier` | same | `DirectPoa`, `AbpoaPoasta`, `SweepgaSeqwish`, and `SkipUnchanged` behavior. |
| `LocalizedPolishReport` | same plus `src/graph_report.rs` serialization helpers | Per-round metrics, dirty sites, chunks, resolver decisions, accepted/rejected replacements, validation status, and runtime/resource counters. |

The new layer should call existing graph builders and replacement validators
rather than duplicating GFA parsing, path replay, or path serialization.

## Inputs And Outputs

### Inputs

The local sequence set can come from any of these sources:

| Input | Required fields | Notes |
| --- | --- | --- |
| SYNG query | `.syng` prefix, query interval or BED/BEDPE row, path/range selector, optional transitive query settings | Primary production path. Uses SYNG to collect local homologous haplotype intervals quickly. |
| FASTA or AGC sequence files | exact source sequences matching SYNG `.names` or direct graph input path names | Required whenever GFA induction needs DNA sequence spelling for selected intervals. |
| Existing seed GFA | GFA text/path, optional source FASTA for exact validation | Allows polishing a seed graph produced by another step without recollecting sequences. |
| Optional seed method config | `syng-derived`, `whole-region-sweepga-seqwish`, `pggb-like`, `existing-gfa` | Selects how the initial graph is built before dirty-region smoothing. |
| Budget config | max region bp, max haplotypes, max chunks, max total chunk bp, max pairwise records/bytes, max wall time, max RSS if enforceable | Budgets control work attempted, not replacement correctness. |

For SYNG-driven operation, the collection step must emit a `LocalSequenceSet`:

```text
region_id
source_query_name
source_query_interval
records[]:
  stable_path_name
  pansn_sample / haplotype / contig when parseable
  source_interval
  strand
  sequence
  original_name_token
```

Names should be the same user-visible names that would appear in final GFA
paths. If a temporary aligner needs unique chunk-local aliases, the alias map
must be private to the chunk and rewritten before validation or final output.

### Outputs

Every run should produce:

| Output | Contract |
| --- | --- |
| Seed graph GFA | The graph after seed induction and before localized polishing. For CLI output this may be a debug/report artifact, not necessarily committed in docs. |
| Polished GFA | Final GFA with exactly the same required path-name set and path spellings as the seed/local sequence set, except for user-requested filtering. |
| Metrics/report | JSON plus optional markdown/TSV with seed metrics, final metrics, per-round deltas, detector hits, selected chunks, resolver tier used, accepted/rejected status, path-validation status, and skip/budget reasons. |
| Path preservation record | Machine-readable result keyed by path name: present/missing/extra/spelling-match. This can reuse `examples/compare_gfa_paths.rs` behavior or an internal equivalent. |

### Path-Name And PanSN Preservation Contract

1. Final GFA path names must exactly match the seed graph path names for all
   retained paths. No `local_`, `chunk_`, `__impg`, temp FASTA aliases, or
   replacement IDs may leak into final `P`/`W` path names.
2. PanSN names remain PanSN names. The sample, haplotype, and contig tokens are
   not rewritten by smoothing. Coordinate suffixes from source collection are
   preserved unless the user explicitly requested coordinate normalization.
3. For GFA 1.1 `W`-line output, PanSN fields must be regenerated from the
   original name metadata, not inferred from temporary chunk aliases.
4. Every replacement stores the original chunk path slices as
   `(path_name, source_interval, orientation, sequence)`. The replacement is
   accepted only if replaying the replacement graph recovers the same sequence
   for each slice.
5. Whole-graph validation runs after all accepted replacements. If the final
   graph fails exact path comparison, the run fails or writes the unmodified
   seed graph with a hard failure in the report. It must not silently emit a
   corrupt polished graph.

## Seed Graph Options

Seed graph construction is separate from localized polishing. The seed can be
high quality or intentionally cheap, but its provenance must be reported.

| Seed option | Existing implementation reference | Use when | Avoid when |
| --- | --- | --- | --- |
| Whole-region SweepGA/FastGA + seqwish | `src/commands/align.rs`, `src/syng_graph.rs::build_gfa_from_paf_and_sequences`, `src/commands/graph.rs` | The local region has manageable haplotype count and bp size, and a good region-wide seed is likely to reduce later dirty work. This is the default for moderate SVR panels and C4-sized local loci when resource budgets allow. | Pairwise alignment would exceed configured pair/byte/time budgets, or repeat-rich regions are known to glue distant copies without chunk context. |
| PGGB-like internal engine | `src/commands/graph.rs`, `src/lib.rs` PGGB engine path and `smooth::smooth_gfa` | A control-grade seed is needed, or the validation task wants a PGGB/SmoothXG comparison without relying on external binaries. Use as a bounded control and optionally as a seed for small panels. | It would make the local method indistinguishable from a whole-region PGGB/SmoothXG run, or it hides which localized smoother fixed the graph. |
| SYNG-derived seed | `src/lib.rs::build_syng_region_gfa_from_intervals`, `src/lib.rs::build_syng_local_region_gfa_from_intervals`, `src/commands/syng2gfa.rs` | Fastest collection path, large regions, exploratory diagnostics, and cases where the whole-region aligner is too expensive. This keeps the architecture grounded in SYNG sequence collection. | The seed is so sparse or white-space dominated that all selected chunks immediately exceed smoothing budgets; in that case promote to whole-region or chunk-level induction. |
| Existing input GFA | `src/resolution.rs::resolve_gfa_bubbles`, `src/graph_report.rs` | Downstream debugging, reproduction from saved seed graphs, or comparison against PGGB/SmoothXG controls. | The GFA lacks path names/spellings needed for exact validation and no source sequence set is available. |

The implementation should default to SYNG sequence collection plus
whole-region SweepGA/FastGA + seqwish seed when the collected region is within
budget. If the whole-region seed exceeds budget, it should fall back to a
SYNG-derived seed and rely on chunk-level induction only for dirty chunks.

## Dirty-Region Detectors

Dirty-region detection runs on the current graph at the start of every round.
Each detector emits `DirtySite` records. These records drive selection order
and report diagnostics; they do not by themselves force acceptance or rejection
of a replacement.

| Detector | Signal | Existing references | Output |
| --- | --- | --- | --- |
| Low path depth / white space | High `path_white_space_bp_p95/p99/max`, long top white-space jumps, sparse coverage runs, low bp-weighted depth in local interval | `src/graph_report.rs::describe_gfa`, `GraphMetrics`, `top_white_space_jumps`, `top_white_space_regions` | Path-coordinate intervals spanning long jumps or sparse path-depth regions. |
| Singleton bp | High singleton bp or singleton nodes in a local interval | `src/graph_report.rs` node coverage metrics | Candidate intervals where private singletons should be aligned against flanks or confirmed as true private sequence. |
| Path jumps after sort | Large adjacent path jumps visible only after `gfasort`/`Ygs` or sort-aware graph report, especially when final order creates long edges between nearby path steps | `src/lib.rs::apply_graph_transforms`, `src/graph_report.rs`, existing sort stages in `src/main.rs` | Intervals around the jump endpoints, reported with pre-sort and post-sort coordinates when both are available. |
| Residual loops and self-loops | Direct self-loop edges, adjacent same-node path steps, repeated self-loop runs | `src/gfa_self_loops.rs`, `src/graph_report.rs` self-loop fields | Loop nodes plus the smallest path ranges containing the loop traversals. |
| Bubble/path replay failures | POVU decomposition failure, candidate extraction failure, replacement path replay mismatch, missing paths in a local candidate graph | `src/resolution.rs::validate_replacement_paths`, POVU use in `src/graph_report.rs`, tests under `tests/test_crush_integration.rs` | Hard dirty sites when the graph cannot be reliably replayed; these should be prioritized or skipped with explicit failure reason. |
| Underaligned chunks | Whole-region or chunk PAF has too few records, no indexable exact run, excessive dropped alignments, low local alignment density, or seqwish emits mostly disconnected path-specific pieces | `src/syng_graph.rs::build_gfa_from_paf_and_sequences`, adaptive seqwish diagnostics, `scripts/local_compression_testbed.py` metrics | Chunk candidates marked for smaller flanks, lower-tier POA, or fallback skip depending on budget. |

Detector severity should be transparent. A recommended scoring model is:

```text
severity = weighted sum of normalized detector values
priority = severity, then path-span bp descending, then stable region id
```

The weights are for ordering only. They are not acceptance gates.

## Chunk Selection

Chunk selection converts dirty sites into non-overlapping `DirtyChunk` records.
This is where localized smoothing is made explicit.

1. Map every dirty site to one or more path-coordinate intervals. Prefer a
   reference path or user-specified anchor path when present; otherwise use the
   dominant path embedding reported by POVU or the path with the clearest
   source interval.
2. Sort intervals by `(source contig, path name, start, end)` and merge with an
   interval-union sweep. This keeps selection near `O(n log n)` and avoids
   pairwise dirty-site comparisons.
3. Merge nearby sites when the gap between them is below `merge_gap_bp`, when
   they share a loop/bubble identifier, when there is no clean invariant node
   between them, or when smoothing them separately would produce overlapping
   replacement flanks.
4. Include flanks outside the dirty core. Defaults should be configurable, but
   a practical starting policy is `max(100 bp, min(500 bp, 10% of core span))`
   with extension to the nearest stable shared node when that extension stays
   within budget.
5. Preserve path context by extracting the same chunk interval from every path
   that traverses the selected local region, plus any path needed to keep
   homology unambiguous in duplicated flanks.
6. Avoid quadratic blowups before calling aligners. Enforce chunk budgets on
   haplotype count, total input bp, max path length, pairwise combination count,
   expected PAF bytes, and per-round total chunk bp. If a chunk exceeds budget,
   try one split along clean invariant gaps; otherwise skip it unchanged with a
   report reason.
7. Avoid over-fragmentation. Do not create chunks smaller than the resolver's
   minimum useful flank/core size unless the detector is a hard replay failure.
   Prefer one merged chunk over several adjacent chunks when the gap is made of
   singleton bp, a sort jump, or path-specific white space.
8. Select a maximal non-overlapping set per round. When chunks overlap, choose
   the higher-severity merged chunk; ties prefer the chunk with more dirty
   detectors and then the smaller total bp. Deferred overlaps can be retried in
   the next round after accepted replacements are applied.

Every chunk must carry:

```text
chunk_id
round
dirty_site_ids[]
core_intervals[]
flanked_intervals[]
selected_path_names[]
source_slice_records[]  # path name, source interval, strand, sequence
budgets
resolver_hint
```

## Resolver Tiers

Resolution is tiered by chunk size, haplotype count, and evidence quality. The
same path-correctness gate applies to every tier.

| Tier | Resolver | Existing references | Use when |
| --- | --- | --- | --- |
| 0 | Skip unchanged | report-only | Chunk exceeds hard budget, lacks sequence, has no safe boundaries, or all resolver attempts failed path validation. |
| 1 | SPOA/POA or abPOA/POASTA over flanked chunk paths | `src/resolution.rs::build_poa_replacement`, `build_abpoa_replacement`, `build_poasta_replacement`; `src/smooth.rs::poasta_block` | Small flanked bubbles, short indels/SNPs, low to moderate haplotype count, or underaligned chunks where pairwise induction has insufficient anchors. |
| 2 | SweepGA/FastGA + seqwish over flanked chunk paths | `src/commands/align.rs`, `src/syng_graph.rs::build_gfa_from_paf_and_sequences`, `src/resolution.rs::build_sweepga_seqwish_replacement` | Larger chunks where pairwise alignment is affordable and expected to recover shared structure better than POA. |
| 3 | Chunk split and retry | chunk selector plus tier 1/2 | A chunk is biologically dirty but exceeds pairwise or memory limits; split on clean invariant gaps and retry children in the same or next round. |

Suggested initial thresholds should be conservative and configurable:

```text
tier 1:
  median path length <= 2 kb, max path length <= 10 kb, total input bp <= 2 Mb
tier 2:
  total input bp <= region budget, pair count <= max_pair_alignments,
  expected PAF bytes <= max_paf_bytes
fallback:
  split once or skip unchanged; never accept a path-corrupt replacement
```

The resolver should record why it chose a tier. If a higher tier fails due to
zero alignments or seqwish underalignment, it may retry a lower tier when the
lower tier is within budget. If a lower tier preserves paths but produces an
ugly graph, it is still accepted unless exact path validation fails; the graph
quality movement is reported as diagnostic evidence.

## Iteration Rule

Each run proceeds in rounds:

1. Build or load the seed graph.
2. Compute seed graph report and exact path baseline.
3. Detect dirty sites on the current graph.
4. Merge dirty sites into non-overlapping chunks with flanks.
5. Resolve chunks in deterministic order subject to per-round budgets.
6. Validate every replacement against chunk source spellings.
7. Apply only path-correct replacements.
8. Validate whole-graph path names and spellings.
9. Recompute graph-quality diagnostics and continue if dirty sites remain and
   budgets allow.

Acceptance rule:

```text
accept replacement iff replacement path replay is exact
and whole-graph path comparison remains exact after applying the round
```

Stop when any of these is true:

- no dirty sites are detected;
- no chunks are selectable after merging and budget checks;
- a round accepts zero replacements;
- the same dirty-region signature repeats after an accepted round;
- `max_rounds`, wall time, memory, chunk count, total chunk bp, or pairwise
  alignment budget is exhausted.

Do not add hidden metric gates. Segment count, singleton bp, white-space p99,
self-loop counts, replay ratio, and topology scores are logged before/after and
used for human diagnosis and next-round prioritization only. Exact path
corruption is the hard failure.

## Report Schema

The report should be stable enough for validation tasks and downstream
implementation agents:

```json
{
  "region_id": "sample#0#chr:start-end",
  "seed_source": "whole-region-sweepga-seqwish",
  "path_contract": {
    "expected_paths": 0,
    "observed_paths": 0,
    "missing_paths": [],
    "extra_paths": [],
    "spelling_mismatches": [],
    "status": "pass"
  },
  "rounds": [
    {
      "round": 1,
      "metrics_before": {},
      "dirty_sites": [],
      "chunks": [],
      "accepted_replacements": 0,
      "rejected_replacements": 0,
      "metrics_after": {},
      "diagnostic_delta": {}
    }
  ],
  "budget": {},
  "status": "pass"
}
```

For markdown reports, include a short top-level conclusion, then tables for
seed/final metrics, per-round chunks, rejected replacements, and validation
controls. Avoid embedding bulky logs or generated GFA contents in docs.

## Validation Plan

Validation must prove that this is localized SYNG sequence-set smoothing, not
C4-only parameter churn.

1. Synthetic fixtures:
   - Use `tests/test_data/local_compression/manifest.json` and
     `scripts/local_compression_testbed.py` as the first validation target.
   - Require exact path names and spellings for every fixture row.
   - Exercise dirty detectors with fixtures that already target singleton bp,
     fake repeat anchors, nested bubbles, duplicated flanks, loops, inversion,
     and adjacent bubbles.
   - Add new fixtures only as small text/FASTA/JSON artifacts; do not commit
     generated GFA outputs.
2. Unit and integration tests:
   - Add targeted tests for dirty-site extraction from `GraphReport`.
   - Add chunk merge/flank tests with overlapping, adjacent, and duplicated
     flank cases.
   - Add resolver-routing tests that distinguish SPOA/abPOA from
     SweepGA/seqwish by chunk size and budget.
   - Reuse exact path comparison helpers for replacement and final graph
     validation.
3. C4 diagnostic validation:
   - Run the current C4 local region path with a bounded number of rounds and
     an explicit seed choice.
   - Compare against the existing whole-region SweepGA/seqwish seed, PGGB-like
     control, and current crush/smooth variants.
   - Report singleton bp, path white-space p99/max, self-loops, path steps,
     segment bp, runtime, and RSS. These are diagnostics, not hidden gates.
4. Non-C4 SVR panel:
   - When local data is available, use a small panel of non-C4 structural
     variant regions collected through SYNG.
   - Include insertions, repeat-rich alleles, inversion-like paths, and nested
     or adjacent variants.
   - Keep panel size small enough for repeatable local validation before
     expanding.
5. PGGB/SmoothXG controls:
   - Use internal PGGB-like and SmoothXG-like controls from `src/commands/graph.rs`
     and `src/smooth.rs` where external `pggb`/`smoothxg` binaries are not
     available.
   - If external tools are locally available, record exact versions and command
     lines in the report, but commit only markdown/JSON/TSV summaries.
6. Runtime targets:
   - Synthetic CI subset should complete in seconds to low minutes.
   - Small non-C4 SVR panel should complete in a local developer run without
     exceeding configured budgets.
   - C4 diagnostic run should show bounded localized work: report selected
     chunk count and total chunk bp so reviewers can see that runtime is not a
     disguised whole-region exhaustive sweep.

## Implementation Ordering

Shared graph-output and replacement code should remain sequential. These files
are high collision risk and should not be edited by parallel subtasks:

- `src/resolution.rs`
- `src/lib.rs`
- `src/main.rs`
- `src/graph_report.rs`
- `src/smooth.rs`
- `src/syng_graph.rs`
- `src/commands/graph.rs`

Recommended sequence:

1. Explicit local seed induction driver:
   - Add the config/API that collects a SYNG local sequence set and builds a
     seed graph with a chosen seed option.
   - Wire only enough CLI/pipeline parsing to exercise the seed choices.
2. Dirty-region detector layer:
   - Add `DirtySite` extraction from existing graph reports, self-loop reports,
     replay failures, and underalignment summaries.
3. Chunk selection:
   - Implement interval merge/flank/budget logic and tests before resolver
     mutation code.
4. Resolver tiers:
   - Reuse existing POA/abPOA/POASTA and SweepGA/seqwish builders through a
     narrow local-polish interface.
   - Keep validation and application in the existing replacement code path
     where possible.
5. Iterative integration:
   - Add round loop, per-round report, exact final path validation, and stop
     conditions.
6. Validation and synthesis:
   - Run synthetic fixtures, C4 diagnostic, non-C4 SVR panel when available,
     and PGGB/SmoothXG controls.

Docs and fixtures can proceed independently because they do not need to edit
the shared graph-output code:

- fixture metadata and tiny FASTA additions under `tests/test_data/local_compression/`;
- documentation under `docs/designs/` and `docs/evaluations/`;
- report schema examples using small JSON snippets;
- validation notes that summarize generated artifacts without committing them.

## Artifact Policy

This design expects implementation and validation tasks to commit only source,
tests, docs, and tiny metrics summaries. Do not commit generated GFA, GFA.gz,
GFA.zst, PNG, PDF, bulky stdout/stderr logs, per-fixture command artifacts, or
any staged blob over 1 MiB. Generated graph outputs should stay under `target/`,
`/tmp`, or an explicitly ignored local data directory and be summarized in
markdown/JSON/TSV only.
