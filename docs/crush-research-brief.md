# Crush Research Brief

**Branch:** `eg/c4-crush-resolution-controls`
**Handoff read:** `docs/c4-crush-handoff.md` (full)
**Date:** 2026-05-24

---

## 1. What "Crush" Means in This Project

"Crush" is the name for **exact path-preserving bubble resolution**: detecting bounded single-entry/single-exit bubbles (flubbles) in a blunt GFA via POVU decomposition, replacing each bubble with a locally re-induced graph using a selected alignment method, then validating that all path sequences are unchanged before accepting.

### SyngGraph Crush (syng-specific path)

Triggered via: `impg query -o gfa:syng:crush[,options]` or `gfa:syng:blunt:crush`.

Pipeline:
1. Build blunt syng GFA from the syng index: `build_syng_region_gfa_from_intervals()` → `write_gfa_with_mode(..., SyngGfaMode::Blunt)` (`src/commands/syng2gfa.rs`)
2. Apply `apply_graph_transforms()` → `resolution::resolve_gfa_bubbles()` (`src/lib.rs:745-756`)

Enforced constraint: crush requires blunt mode. Raw syng GFA is explicitly rejected at `src/main.rs:3277-3281`:
```
"Invalid --gfa-engine '{}': crush requires blunt syng GFA; use syng:crush or syng:blunt:crush"
```

This path only works when `impg` is given a `.syng` index (`-a <prefix>.syng`). It builds the graph in-process from syng data structures before crushing.

### General Bubble Crush (any-graph path)

Two entry points:
- **Standalone CLI:** `impg crush -g input.gfa -o output.gfa` — reads any existing blunt GFA, calls `resolution::resolve_gfa_bubbles()` (`src/main.rs:7875`). CLI is defined at `src/main.rs:4692-4693`.
- **Inline pipeline stage:** Any GFA engine + `:crush` stage: `impg query -o gfa:<engine>:crush`. The crush stage is parsed at `src/main.rs:3218-3225` and applied via `apply_graph_transforms()`.

Both entry points use the identical `ResolutionConfig` → `resolve_gfa_bubbles()` path in `src/resolution.rs:413`.

### Resolution Methods (both paths)

All methods are in `src/resolution.rs`:

| Method | Enum | Description |
|--------|------|-------------|
| `auto` | `ResolutionMethod::Auto` | Small bubbles → SPOA; large → AllWave/SweepGA |
| `poa` | `ResolutionMethod::Poa` | Direct SPOA replacement |
| `poasta` | `ResolutionMethod::Poasta` | POASTA replacement, falls back to SPOA |
| `star-biwfa` | `ResolutionMethod::StarBiwfa` | Star-column BiWFA graph |
| `allwave` | `ResolutionMethod::Allwave` | Sparse many-to-many BiWFA + seqwish + SPOA polish |
| `sweepga` | `ResolutionMethod::Sweepga` | SweepGA/FastGA pair selection + seqwish + SPOA polish |

---

## 2. Current State Inventory

### What's Done

- **Core resolver** (`src/resolution.rs`, 4009 lines): complete, stable, heavily tested.
  - `resolve_gfa_bubbles()` at line 413: public entry point from any GFA string.
  - `resolve_graph_bubbles()` at line 429: internal multi-round frontier loop.
  - `find_candidate_frontier()` at line 1107: POVU decomposition, flubble selection, budget gating.
  - `apply_replacement_frontier()` at line 1452: rewrite graph, validate, or error.
  - `validate_replacement_paths()` at line 2436: exact path-sequence validation per replacement.
  - `path_sequences_equal()` at line 2796: round-level path-preservation guard.

- **All resolution methods**: fully implemented — `build_poa_replacement()`, `build_poasta_replacement()`, `build_biwfa_inmemory_replacement()`, `build_allwave_seqwish_replacement()`, `build_sweepga_seqwish_replacement()`.

- **Visual-tail quality scoring** (`graph_quality()` at `src/resolution.rs:994`): computes a score weighted toward long whitespace bridges (p99, max). **Currently diagnostic/log-only** — not used for acceptance decisions (changed in commit `0af1a4c`).

- **Standalone `impg crush` command** (`src/main.rs:4692`): full CLI with all resolution parameters. Dispatches to `resolve_gfa_bubbles()` at line 7875.

- **SweepGA auto-frequency** (`replacement_sweepga_kmer_frequency()` at `src/resolution.rs:1628`, commit `98fd538`): default is now `0` (auto), which scales to `max(1000, traversal_count × 10)` instead of hard-coded `10`. Prevents seeding starvation in bubble-local repetitive contexts.

- **Scaffold filter** (`replacement_scaffold_filter` field in `ResolutionConfig`, commit `f453983`): now configurable via CLI. Default `"1:1"`.

- **Path validation** (commit `0af1a4c`): Quality-score-based round rejection was **removed**. Replacements are now accepted if and only if exact path-sequence validation passes. The old `DEFAULT_MAX_ROUND_SCORE_GROWTH = 0.02` field is gone.

- **Tests**:
  - 40+ unit tests in `src/resolution.rs` (lines 2967–4009): SNP, insertion, deletion, SV, long-gap, allwave, pairwise-induced, frontier selection, scaffold filter parsing, frequency auto-scale.
  - CLI integration test: `tests/test_syng_integration.rs:207` (`test_crush_cli_resolves_blunt_gfa`).
  - Parameter parsing tests: `src/main.rs:12536–13059`.
  - Path-preservation test: `src/lib.rs:1307` (`test_crush_transform_is_applied_path_preserving`).

### What's Broken / Unproven

- **Current HEAD is behaviorally different from the known-good C4 artifact in two ways**:
  1. `98fd538` changed the default SweepGA frequency from `10` to auto-scaled. The known-good run used `10` (implicit default at the time). CLI string `gfa:syng:crush,method=sweepga,aligner=fastga,...` is no longer semantically identical.
  2. `0af1a4c` removed round-level quality-score gating. The known-good artifact was produced under old code that rejected rounds if visual-tail score grew >2%. Current code never rolls back a round on quality grounds.

- **Unverified**: Does current HEAD reproduce the known-good GRCh38 C4 metrics (465 paths, ~422 kb segments, p99 link jumps ~45, 25414 white-space bridges)? This is the critical open question from the handoff.

- **No integration regression test**: No test freezes the known-good C4 graph-stat output. Handoff explicitly requests one (`docs/c4-crush-handoff.md` lines 244-249).

- **Remaining graph quality issues** in the known-good artifact:
  - 25414 path white-space bridges ≥1000 bp
  - link jumps p99=45, max=3270
  - The handoff identifies this as still underaligned in places.

### What's Missing

- **Integration/regression test** at graph-stat level for GRCh38 C4. The handoff specifies concrete acceptance criteria (paths=465, segment bp ~422 kb, link p99 bounded, no white-space explosion, exact path-sequence match).
- **No benchmarking infrastructure**: no `benches/` directory, no criterion harness. Performance is tracked only via adhoc `/usr/bin/time -v` output in `data/` artifact dirs.
- **General bubble crush via non-syng engines is untested at integration level** beyond the one small synthetic test. The `impg crush` CLI test covers only a 2-traversal SNP bubble.

---

## 3. Goal Restatement

The user wants three things:

### A. SyngGraphs crushed (syng-specific path)

`impg query -a <prefix>.syng -o gfa:syng:crush[,options]`

This builds the full population graph from the syng index and crushes it. The known-good command is documented in `docs/c4-crush-handoff.md:43-54`. The syng path is the primary production path.

### B. Bubbles crushed from any graph (general path)

`impg crush -g input.gfa` or `impg query -o gfa:<engine>:crush`

This applies the same `resolve_gfa_bubbles()` to any blunt GFA produced by any engine. Useful for post-processing existing graphs (e.g., frozen CHM13 debug graphs) or non-syng construction pipelines.

### C. FAST — performance is first-class

The known-good C4 run consumed:
- Wall: 4:39 (mostly query + alignment)
- Max RSS: 55 GB

The handoff does not explicitly state regression budgets for time or memory, but "do not increase" is implied. There is no formal benchmark harness. Performance work needs one.

---

## 4. Constraints Discovered

### Invariants That Must Hold

1. **Exact path-sequence preservation**: every replacement must pass `validate_replacement_paths()` (`src/resolution.rs:2436`) — sequence lengths per path must be identical before and after. This is the only hard gate on acceptance in current code.

2. **Round-level consistency**: `path_sequences_equal()` (`src/resolution.rs:2796`) is called after the full frontier is applied. If it fails, the round is an error.

3. **Blunt GFA only**: crush parses GFA with 0M links (`parse_gfa()` internal to `resolution.rs`). Raw syng (variable overlap) will silently mis-parse or be explicitly rejected (`src/main.rs:3277`).

4. **POVU root = first path**: `find_candidate_frontier()` uses the first GFA path as the POVU root (`src/resolution.rs:1107`). Bubble spans are measured in this root's coordinate space. Changing path order changes which candidates are eligible under `max_bubble_span`.

5. **"Validated" in "Accept validated crush replacements" (commit `0af1a4c`)**: means path-sequence-preserving validation only — not quality-score validation. Quality score is informational. The old quality gating (`DEFAULT_MAX_ROUND_SCORE_GROWTH`) is gone. Current code cannot roll back a round based on graph quality.

6. **SweepGA FastGA backend is always all-vs-all**: sparse pair dispatch is not supported with FastGA (`src/resolution.rs` — "FastGA does not have a cheap pairs-file mode"). Sparse pairs only work with wfmash. This means the known-good C4 run (which used FastGA) ran full all-vs-all.

### Edge Cases

- `max_pair_alignments=0` / `max_paf_bytes=0`: disables those guards entirely (used in known-good command).
- `sweepga_kmer_frequency=10` pins the old behavior; `kmer_frequency=0` (new default) auto-scales.
- `min_traversal_len` and `max_traversal_len` can be set such that `min > max`: this forces the sweepga/allwave path for large bubbles while skipping small ones (tested in `test_crush_stage_allows_large_pairwise_pass_above_direct_budget` at `src/main.rs:12588`).

### API Surfaces That Must Not Break

- `resolve_gfa_bubbles(gfa: &str, config: &ResolutionConfig) -> io::Result<ResolvedGfa>` (`src/resolution.rs:413`) — public API used by both `lib.rs` and `main.rs`.
- `ResolutionConfig` struct fields — used directly in `src/main.rs` CLI dispatch and tests.
- `impg crush` CLI flags — external users invoke these; renaming breaks scripts.

---

## 5. Open Questions

1. **Does current HEAD reproduce the known-good GRCh38 C4 metrics?** The handoff mandates two runs (current auto-frequency vs legacy `kmer-frequency=10`) before any further algorithm changes. This has not been done. If HEAD regresses, which of the two behavioral changes is responsible?

2. **Is the new auto-frequency policy (commit `98fd538`) actually better or just different?** For small bubbles (2–10 traversals), the auto formula gives 1000–100 vs the old flat `10`. For large bubbles (463 traversals), it gives 4630. The CHM13 debug graph explosion was fixed by raising frequency, but it is unknown whether raising it for C4's large repetitive bubbles is helpful or whether it causes over-filtering.

3. **With round-level quality gating removed (commit `0af1a4c`), can rounds now regress graph quality indefinitely?** The known-good run rejected a second round that added three replacements because quality grew >2%. Under current code, that round would be accepted. Is the result better or worse?

4. **What are the remaining 25414 long white-space bridges from?** The handoff lists four possible sources (syng/blunt glue, SweepGA alignment, polish/POA, sorting). Identifying the dominant source is a prerequisite for designing the fix.

5. **What performance budget is acceptable?** The known-good run is 4:39 wall / 55 GB RSS. Is this the target to match, or the ceiling to improve upon? No formal regression budget exists.

6. **Should there be a quality score guard as an option (not default)?** Now that validation is path-only, users have no guard against a crush round that produces a valid-but-ugly graph. Adding an opt-in `--max-round-score-growth` could reintroduce the old behavior for cases where quality feedback is wanted.

7. **Integration test fixture**: the handoff requests a regression test asserting broad graph-stat metrics for GRCh38 C4. What is the infrastructure for this? No integration test currently runs against real population data. It would need the HPRC assemblies and syng index, which are too large for CI. Where does this live?

---

## Source File Reference

| File | Purpose | Key Items |
|------|---------|-----------|
| `src/resolution.rs` | Core bubble resolver | `resolve_gfa_bubbles()` (L413), `ResolutionConfig` (L29), `validate_replacement_paths()` (L2436), `graph_quality()` (L994), `find_candidate_frontier()` (L1107), `replacement_sweepga_kmer_frequency()` (L1628) |
| `src/lib.rs` | Entry point for GFA transforms | `apply_graph_transforms()` (L745), `EngineOpts.crush_config` (L81), `build_syng_region_gfa_from_intervals()` (L644) |
| `src/main.rs` | CLI definition and dispatch | `Args::Crush` (L4693), `parse_crush_stage()` (L2536), syng+crush blunt enforcement (L3277) |
| `src/commands/syng2gfa.rs` | Syng GFA serialization | `SyngGfaMode::Blunt` (L122), `write_gfa_with_mode()` |
| `src/graph_pipeline.rs` | Pipeline stage parser | `GraphPipelineSpec::parse()` (L33), crush stage name handling |
| `tests/test_syng_integration.rs` | CLI integration tests | `test_crush_cli_resolves_blunt_gfa()` (L207) |
| `docs/c4-crush-handoff.md` | C4 state snapshot | Known-good artifact paths, baseline reproduction procedure, open questions |
