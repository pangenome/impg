# Crush Design

**Branch:** `wg/agent-7/crush-design` (parent: `eg/c4-crush-resolution-controls`)
**Predecessor:** `docs/crush-research-brief.md`, `docs/c4-crush-handoff.md`
**Date:** 2026-05-24
**Goal:** Specify how the fast unified crush is structured and which speed wins
the implementer should ship first.

The crush algorithm is already correct on real data — the C4 known-good
artifact validates that. The work this design covers is to make crush
**faster without changing the path-preservation contract or the C4 baseline
output**, exposing one unified entry point that serves both the SyngGraph
pipeline and the standalone bubble-crush CLI.

---

## 1. Unified API

The unified API already exists. The function does *not* need to be renamed
or wrapped — it needs an in-memory sibling so callers that already hold a
parsed graph can skip the GFA string round-trip.

### 1.1 Public surface (current)

```rust
pub struct ResolutionConfig { /* all knobs, unchanged */ }
pub struct ResolvedGfa     { pub gfa: String, pub stats: ResolutionStats }

/// Parse a blunt GFA string, resolve bubbles, return the rewritten GFA.
/// Used by both the syng pipeline (`apply_graph_transforms`) and the
/// standalone `impg crush` CLI.
pub fn resolve_gfa_bubbles(
    gfa: &str,
    config: &ResolutionConfig,
) -> io::Result<ResolvedGfa>;
```

This is the surface external scripts and `src/main.rs` depend on. It MUST NOT
change shape or break existing behavior.

### 1.2 New in-memory entry point (additive)

```rust
/// Resolved in-memory graph: same internal type the resolver mutates between
/// rounds. Avoids the render -> parse round-trip when the caller already
/// holds a parsed graph (e.g. blunt syng GFA built in-process).
pub struct ResolvedGraph {
    pub graph: Graph,            // currently private; promote to pub(crate)
    pub stats: ResolutionStats,
}

/// Same algorithm, but takes an already-parsed graph and returns the resolved
/// graph without a final render. Callers that need GFA text wrap with
/// `render_graph(&resolved.graph)`; callers that hand off to the next pipeline
/// stage (sort, vcf, etc.) operate on the in-memory graph.
pub fn resolve_graph_bubbles_inmemory(
    graph: Graph,
    config: &ResolutionConfig,
) -> io::Result<ResolvedGraph>;

/// Convenience: render the resolved graph back to a GFA 1.0 string.
pub fn render_resolved_graph(resolved: &ResolvedGraph) -> String;
```

`resolve_gfa_bubbles` becomes a thin wrapper:

```rust
pub fn resolve_gfa_bubbles(gfa: &str, config: &ResolutionConfig)
    -> io::Result<ResolvedGfa>
{
    let graph = parse_gfa(gfa)?;
    let resolved = resolve_graph_bubbles_inmemory(graph, config)?;
    Ok(ResolvedGfa {
        gfa: if resolved.stats.resolved == 0 && resolved.stats.bailed == 0 {
            gfa.to_string()         // unchanged-graph fast path (current behavior)
        } else {
            render_resolved_graph(&resolved)
        },
        stats: resolved.stats,
    })
}
```

### 1.3 Caller integration

| Caller                                 | Today                                                        | After                                                                   |
|----------------------------------------|--------------------------------------------------------------|-------------------------------------------------------------------------|
| `apply_graph_transforms` (`src/lib.rs:745`) | `resolve_gfa_bubbles(&gfa_string, …)`                        | unchanged (still has GFA string in hand from `build_syng_region_gfa_*`) |
| `impg crush` CLI (`src/main.rs:7875`)  | `resolve_gfa_bubbles(&fs::read_to_string(input)?, …)`        | unchanged                                                               |
| Future in-process syng pipeline        | could build `Graph` directly and call `_inmemory`            | new, opt-in                                                             |

Because `apply_graph_transforms` already has the GFA as a string (the
serializer in `build_syng_region_gfa_from_intervals` returns a `String`), the
in-memory entry point is most valuable for tests, benchmarks, and a future
syng-side change that emits the in-memory `Graph` directly. We add it now
because the rest of this design uses it as a measurement vehicle.

---

## 2. Algorithm Choice

**Keep the current algorithm.** Do not change resolution methods, scoring,
or selection ordering. The C4 handoff is explicit: reproduce the known-good
artifact first; only then change behavior. The design's job is to make the
*same* algorithm run faster.

### 2.1 Where the time goes today

From the C4 known-good run (4:39 wall, 55 GB max RSS) and from reading
`resolve_graph_bubbles`, the time decomposes roughly as:

1. **Alignment** (SweepGA/FastGA + seqwish per bubble): dominant — minutes.
2. **POVU decomposition per round**: rebuilds POVU's `NativeGfa` from a
   freshly-rendered GFA string every round.
3. **Per-round full-graph rewrite**: `render_rewritten_graph` then
   `parse_gfa` of the result for the next round.
4. **Path-sequence validation**: `path_sequences_equal` walks every path,
   reconstructs every sequence as `Vec<u8>`, allocates an `FxHashMap`.

(2)-(4) is the part this design speeds up. (1) is the algorithm and is
out of scope here — `crush-perf` and `crush-impl-syng/bubble` may revisit
it after the baseline is locked.

### 2.2 Speed targets (quantitative)

Reference workload: GRCh38 C4 known-good command from
`docs/c4-crush-handoff.md:43-54` on 32 threads.

| Metric                              | Today (known-good)     | Target after design     | Stretch |
|-------------------------------------|------------------------|-------------------------|---------|
| Wall                                | 4:39                   | ≤ 4:00 (≥ 15% faster)   | ≤ 3:30  |
| Max RSS                             | 55 GB                  | ≤ 55 GB (no regression) | ≤ 50 GB |
| Path-validation share of wall       | not measured           | ≤ 2 % after change      | ≤ 1 %   |
| Per-round POVU re-parse share       | not measured           | ≤ 5 % after change      | ≤ 2 %   |

Reference micro-workload (synthetic, repeatable, runs in CI):

| Input                                                    | Today  | Target (single thread) |
|----------------------------------------------------------|--------|------------------------|
| `tests/test_data/yeast.chrV.fa.gz`-derived 32-path GFA   | TBD    | crush in < 10 s        |
| 1k-path × 100 kb synthetic blunt graph (gen at test time)| TBD    | crush in < 60 s        |
| Inline 3-segment SNP bubble (existing CLI test)          | < 1 s  | unchanged              |

Targets are *first* measured (current HEAD) and *then* the design's changes
must move them in the right direction. No target is asserted as a regression
gate before the baseline is recorded — see § 5.

### 2.3 Time and memory complexity

Let `S` = segment count, `L` = total segment bp, `P` = path count,
`T` = total path-step count, `B` = bubbles selected per round,
`R` = rounds, `K` = mean traversal count per bubble.

|                              | Current                          | After                            | Notes |
|------------------------------|----------------------------------|----------------------------------|-------|
| GFA parse per round          | `O(L + T)`                       | `O(L + T)` once, then `O(1)`     | parse only on entry |
| POVU decompose per round     | `O(S + T)` (POVU native)         | `O(S + T)` (POVU native)         | unchanged |
| Render full graph per round  | `O(L + T)`                       | removed (in-memory)              | saved per round |
| Re-parse rewritten graph     | `O(L + T)`                       | removed (in-memory)              | saved per round |
| `path_sequences_equal`       | `O(P · L)` + alloc               | `O(T)` hashes, no alloc          | streaming xxhash compare |
| Per-replacement alignment    | algorithm-dominated              | algorithm-dominated              | unchanged |

Memory complexity is bounded by `O(L + T)` for the working graph plus
`O(K · L_bubble)` per active bubble. The design preserves both bounds.

### 2.4 Parallelism plan

| Loop                                                  | Today                       | After                                       |
|-------------------------------------------------------|-----------------------------|---------------------------------------------|
| `path_step_index` per path                            | `rayon::par_iter`           | unchanged                                   |
| candidate enumeration over POVU sites                 | `rayon::par_iter`           | unchanged                                   |
| `materialize_candidate_sequences` per candidate       | `rayon::par_iter`           | unchanged                                   |
| `build_replacement_with_method` per candidate         | `rayon::par_iter`           | unchanged                                   |
| `path_sequences_equal` (path-by-path hash)            | sequential                  | `rayon::par_iter` over paths, streaming hash |
| `render_resolved_graph` segment write                 | sequential `String::push`   | sequential (cheap; <2 % wall after change)  |

Crate choice: `rayon` — already in the codebase and already used in this
module. No SIMD work in this design; the per-round savings come from
eliminating string round-trips, not from inner-loop vectorization.

Expected speedup from parallel validation alone: the C4 graph has 465 paths;
parallel xxhash over them across 32 threads turns an O(P · L) sequential
walk into an O(L) parallel walk, removing path validation from the wall.

---

## 3. Data Structures

### 3.1 In-memory graph (already present)

`Graph` in `src/resolution.rs` is already the right shape:

```rust
pub(crate) struct Graph {
    pub segments: Vec<Segment>,   // id + Vec<u8> seq, dense ids
    pub paths:    Vec<Path>,      // name + Vec<Step>
}
pub(crate) struct Step { node: usize, rev: bool }   // 9 bytes packed
```

Cache locality is good: segment sequences are contiguous `Vec<u8>`, paths
are contiguous `Vec<Step>`, and the hot inner loops (path sequence
materialization, POVU adjacency) walk these in order. Promote `Graph`,
`Segment`, `Path`, `Step` to `pub(crate)` so the new in-memory entry point
can expose them without leaking via a `String` boundary, but keep them
private to the crate to preserve the API surface.

### 3.2 Reused parser state

`parse_gfa` builds an `FxHashMap<String, usize>` id→idx table; this is
already on the hot path once per entry. After the design, the table is
built **once** per `resolve_gfa_bubbles` call instead of once per round.

### 3.3 Per-round POVU input

POVU consumes its own `NativeGfa`, currently produced by rendering the
working graph and re-parsing inside POVU. The render is the largest
non-alignment per-round cost on big inputs. There are two options:

| Option | Source                              | Approach                                                                                          | Risk                                                                       |
|--------|-------------------------------------|---------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------|
| A      | crate `povu::NativeGfa::parse(...)` | keep the render → parse path; reduce allocation via a reusable `String` buffer                    | low — same code path, just less churn                                       |
| B      | `povu::NativeGfa::from_graph(...)`  | add a constructor that ingests our in-memory `Graph` directly; no GFA text intermediate           | **medium** — requires a POVU API change or a `pub(crate)` shim in `povu-rs` |

**Pick A first.** Land a reusable buffer (a `String` owned by the resolver
that is `clear()`-and-write-into between rounds) so we stop reallocating
megabytes per round. This is a non-invasive win and unblocks Option B as a
later `crush-perf` task. Option B requires touching `povu-rs` and should
be specified separately if profiling shows POVU parse is still a measurable
share after Option A.

### 3.4 Streaming path-sequence validation

Today `path_sequences_equal` builds two `FxHashMap<String, Vec<u8>>` and
compares them. For the C4 graph that allocates ~hundreds of MB per round.

Replace with a streaming compare:

```rust
fn path_sequences_equal_streaming(before: &Graph, after: &Graph) -> io::Result<bool> {
    if before.paths.len() != after.paths.len() { return Ok(false); }
    let before_by_name: FxHashMap<&str, &Path> =
        before.paths.iter().map(|p| (p.name.as_str(), p)).collect();
    after.paths.par_iter().try_fold(
        || true,
        |acc, after_path| {
            if !acc { return Ok::<bool, io::Error>(false); }
            let Some(before_path) = before_by_name.get(after_path.name.as_str())
                else { return Ok(false); };
            Ok(path_streams_eq(before, before_path, after, after_path)?)
        },
    ).try_reduce(|| true, |a, b| Ok(a && b))
}
```

`path_streams_eq` walks both paths' step sequences byte-by-byte (computing
segments lazily) and short-circuits on the first mismatch. No materialized
`Vec<u8>` of full path sequence is ever held.

This is the single highest-leverage memory and time win in the design, and
it is algorithmically a refactor — no behavior change — so it can be
shipped behind a feature flag and benchmarked head-to-head against the
current implementation for byte equivalence.

---

## 4. Integration With Validated-Replacement Logic

The recent commits this design must respect:

- **`0af1a4c`** (Accept validated crush replacements): removed
  `DEFAULT_MAX_ROUND_SCORE_GROWTH` quality gating. Acceptance is now
  exact path-sequence preservation only.
- **`98fd538`** (Fix SweepGA crush local seed frequency): SweepGA default
  `kmer_frequency=0` → auto-scaled per-bubble (see
  `replacement_sweepga_kmer_frequency`, `src/resolution.rs:1628`).
- **`f453983`** (scaffold filter configurable): `replacement_scaffold_filter`
  default `"1:1"`.

### 4.1 What the design preserves

1. **Two validation gates remain, both as today:**
   - `validate_replacement_paths` (`src/resolution.rs:2436`) — per-bubble,
     before accepting a replacement plan. Every replacement path's sequence
     length must equal the original traversal's.
   - `path_sequences_equal` (`src/resolution.rs:2796`) — per-round, after
     applying the frontier. This is the design's streaming replacement
     target (§ 3.4), but the *contract* is unchanged: round fails iff any
     emitted path's full sequence differs from before.
2. **Auto kmer-frequency for SweepGA stays the default.** No behavior
   change. The C4 reproduction in § 5 must run *both* the current default
   and the pinned `kmer-frequency=10` semantics so the baseline lock-in
   step can confirm which one matches the known-good artifact.
3. **Scaffold filter and num-mappings stay configurable through
   `ResolutionConfig`.** The streaming validator and in-memory entry point
   do not look at these fields — they pass through unchanged.

### 4.2 What the design does NOT do

- It does not reintroduce quality-score gating, even as an opt-in flag.
  That belongs to a separate task (`crush-research` open question 6).
- It does not change which method `auto` picks. Routing
  (`candidate_selection_method`, `src/resolution.rs:838`) is untouched.
- It does not change POVU root-path selection (first path). Bubble span
  measurement stays in root coordinates.

### 4.3 Validation invariants asserted by tests

```rust
// New unit test in src/resolution.rs:
#[test]
fn inmemory_and_string_entry_points_produce_identical_resolved_graph() {
    let gfa = include_str!("../tests/test_data/crush/small_insertion.gfa");
    let config = ResolutionConfig { ..ResolutionConfig::default() };
    let via_string = resolve_gfa_bubbles(gfa, &config).unwrap();
    let parsed     = parse_gfa(gfa).unwrap();
    let via_memory = resolve_graph_bubbles_inmemory(parsed, &config).unwrap();
    assert_eq!(via_string.gfa, render_resolved_graph(&via_memory));
    assert_eq!(via_string.stats, via_memory.stats);
}

#[test]
fn streaming_validator_matches_materializing_validator_on_known_inputs() {
    for gfa_path in crush_test_fixtures() {           // real files
        let before = parse_gfa(&fs::read_to_string(&gfa_path).unwrap()).unwrap();
        let after  = round_trip_rewrite(&before);
        assert_eq!(
            path_sequences_equal(&before, &after).unwrap(),
            path_sequences_equal_streaming(&before, &after).unwrap(),
            "validator disagreement on {}", gfa_path.display(),
        );
    }
}
```

---

## 5. Test Strategy (Real Inputs Only)

Per project rules: no synthetic mocks of GFA structure, no fabricated
`Graph` literals invented in tests. Every input is either a checked-in
real GFA fixture or one generated at test time from real sequence in
`tests/test_data/`.

### 5.1 Correctness inputs (committed or generated from committed FASTA)

| Input                                                              | Source                                                      | Purpose                                                                |
|--------------------------------------------------------------------|-------------------------------------------------------------|------------------------------------------------------------------------|
| Inline 3-segment SNP bubble                                        | `tests/test_syng_integration.rs:218-230` (existing literal) | Smoke test the new in-memory API hits the same code path as the CLI    |
| `tests/test_data/yeast.chrV.fa.gz`                                 | already committed real sequence                             | Build a small blunt syng GFA at test time; crush it; assert paths preserved |
| `tests/test_data/{a,b,c,ref,ref2}.fa`                              | already committed                                           | Build several multi-haplotype blunt graphs; cross-validate string vs in-memory paths |
| `tests/test_data/test.agc` + `tests/test_data/ref.fa.fai`          | already committed                                           | Drive `impg query -o gfa:syng:crush` end-to-end against the AGC bundle |
| `tests/test_syng_render_bundle_preserves_source_namespace`'s graph | generated from `numeric_to_ascii(make_sequence_numeric(...))` in the existing test | Reused as a fixture to assert deterministic crush output on a synthetic-but-real-sequence multi-hap graph |

### 5.2 Performance inputs (real, large-enough-to-measure)

| Input                                                               | Where it lives                                  | Notes                                                                  |
|---------------------------------------------------------------------|-------------------------------------------------|------------------------------------------------------------------------|
| GRCh38 C4 syng + AGC bundle                                         | `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.{syng,agc}` (developer machine, not CI) | The handoff baseline. **The single authoritative perf workload.** Runs are recorded under `data/` with `/usr/bin/time -v` as today. |
| `tests/test_data/yeast.chrV.fa.gz`-derived 32-path GFA              | generated at test/bench time                    | Repeatable, CI-friendly micro-benchmark                                |
| Frozen CHM13 C4 debug GFA (if needed for narrow regression)         | `data/c4_crush_rnd_20260523T235304Z/...mask_min3.gfa` (developer machine) | Mentioned in the handoff as the debugging fixture; not for primary perf measurement (different range, different shape) |
| `q.gfa` (root-of-tree scratch from main worktree)                   | not present in this worktree; mentioned in handoff as untracked scratch | If the implementer wants a quick crush input, generate one fresh from the yeast fixture instead — don't rely on a scratch file that may move |

`vendor/gfaffix` is referenced in the task description but is not checked
out in this worktree (the submodules under `vendor/` are empty
directories — see `vendor/gfaffix/` and `vendor/syng/` listings). The
implementer should `git submodule update --init vendor/gfaffix` if they
want gfaffix-shipped examples, and use those only as supplementary
sanity-check inputs. The yeast and AGC fixtures above are the
authoritative real inputs.

### 5.3 Benchmark harness (new, minimal)

Add a `benches/crush.rs` driven by `criterion` (already a common dep
choice; if it is not in `Cargo.toml`, vendor it under `[dev-dependencies]`
in the same PR as the harness). The harness:

1. Generates the yeast-derived blunt GFA at bench startup
   (`impg syng` library calls, not subprocess), so the input is
   reproducible without checked-in binary blobs.
2. Runs `resolve_gfa_bubbles` and `resolve_graph_bubbles_inmemory` on
   the same input, reports wall time and (via `cap_max_rss`) RSS for both.
3. Runs the streaming validator and the current validator on the same
   pre/post pair and reports their wall.

Output is `criterion`'s JSON in `target/criterion/`. The C4 perf run
remains script-driven (`docs/c4-crush-handoff.md:156-213`) and is not
something CI runs; it is a developer-machine procedure documented in this
design's § 5.4.

### 5.4 C4 baseline lock-in procedure (precondition for any algo change)

Before the implementer ships any code change beyond the in-memory entry
point and the streaming validator, they MUST run the two-baseline
reproduction from `docs/c4-crush-handoff.md:156-213` (current auto-frequency
and pinned `kmer-frequency=10`) and record the resulting graph-stat reports
under `data/c4_design_baseline_<UTC>/`. Then:

1. Pick whichever of the two matches the known-good
   `C4A.parent5k.sweepga_allvsall_fastga.k311.done.Ygs.report.md`
   metrics within ±2% on (paths, segment bp, link p99, white-space bridges).
2. Set the resolver defaults so that command, unchanged, hits that
   semantically-equivalent run.
3. Land the design's speed work and rerun the same command end-to-end. If
   the graph-stat report drifts more than ±2% on those four metrics, the
   speed change is rejected and rolled back.

This is the binding regression test for the design. There is no looser
contract.

### 5.5 What the implementer MUST add as a smoke scenario

A new `tests/smoke/scenarios/crush_inmemory_matches_string.sh` that:

1. Generates a small blunt GFA from `tests/test_data/ref.fa` and friends
   via `impg syng` + `gfa:syng:blunt`.
2. Runs `impg crush -g <gfa> -o /tmp/out_string.gfa`.
3. Runs the new in-memory path through a tiny Rust helper bin
   (`tests/smoke/bin/crush_inmemory.rs`) on the same input.
4. Asserts byte-equal output.

Listed under `owners` for `crush-impl-syng` / `crush-impl-bubble` /
`crush-integrate` so the smoke gate catches any divergence.

---

## 6. Risk List

| # | Risk                                                                                          | Likelihood | Impact   | Mitigation                                                                                                          |
|---|-----------------------------------------------------------------------------------------------|------------|----------|---------------------------------------------------------------------------------------------------------------------|
| 1 | Streaming validator disagrees with materializing validator on an edge case (reverse strands, duplicate path names) | medium     | critical | § 4.3 cross-validates the two on every committed crush fixture before the materializing version is removed. The materializing version is kept behind `cfg(test)` as the oracle. |
| 2 | In-memory entry point and string entry point diverge silently (e.g. one keeps an unused-segment ID the other drops) | medium     | high     | § 4.3 asserts byte-equal output across every fixture. Smoke scenario in § 5.5 enforces this in CI.                  |
| 3 | C4 baseline lock-in shows current HEAD does *not* reproduce the known-good artifact          | high       | high     | The lock-in (§ 5.4) is the first design milestone. If HEAD regresses, the design's first PR pins the default (and only the default) of `sweepga_kmer_frequency` to `10` for the syng path, restoring pre-`98fd538` semantics, and `crush-perf` later proposes a separate default if data justifies it. |
| 4 | Promoting `Graph` to `pub(crate)` accidentally exposes internals through `pub use`            | low        | medium   | No `pub use Graph` — only the new typed entry points (`resolve_graph_bubbles_inmemory`, `render_resolved_graph`) cross the crate boundary, and they take `Graph` by value. External callers stay on `&str`-typed `resolve_gfa_bubbles`. |
| 5 | POVU re-parse stays expensive even with a reusable buffer                                    | medium     | medium   | Document Option B in § 3.3 as a follow-up `crush-perf` task; gather a profile after the buffer change to confirm. |
| 6 | Parallel streaming validator runs out of bandwidth on small inputs (rayon overhead > work)   | low        | low      | Use `par_iter().with_min_len(8)` so paths shorter than ~8 don't fan out across threads; sequential fallback is `O(T)` and still fine. |
| 7 | Benchmark numbers regress on machines without 32 cores                                       | medium     | low      | The targets in § 2.2 are stated for 32 threads; smaller machines are expected to scale roughly linearly down. The C4 perf run is the binding number, not the micro-benchmarks. |
| 8 | The implementer fans out the design into too many PRs and the changes drift                  | medium     | medium   | The design is one PR worth of work: in-memory entry point + streaming validator + benchmark harness + baseline lock-in script. `crush-impl-syng` and `crush-impl-bubble` should consume this as one batch, not split. |
| 9 | Submodule `vendor/gfaffix` not checked out, gfaffix-derived inputs unavailable               | observed   | low      | § 5 calls this out and routes test inputs through `tests/test_data/` instead. The implementer can `git submodule update --init` if they want gfaffix examples; the design does not depend on them. |
| 10 | Quality-tail score logged but never gated; user-visible "ugly" graph slips through         | low        | medium   | Out of scope for this design. The handoff (open question 6) tracks this; it is `crush-perf`'s call whether to add an opt-in flag. The design does not foreclose it. |

---

## 7. Out of Scope (Explicit)

- New resolution methods.
- Changes to POVU itself, including the proposed `NativeGfa::from_graph`
  (called out as a follow-up under § 3.3 Option B, not part of this design).
- Reintroducing round-level quality gating.
- Changing default `sweepga_kmer_frequency` semantics (the baseline
  lock-in step decides the default; this design does not pick it).
- Replacing rayon with another scheduler.
- SIMD work in inner alignment loops.

---

## 8. Provenance & Absorptive-Capacity Notes

The pieces this design imports from outside the current resolver:

| Import                                          | Source                                                                                  | Absorptive capacity                                                          |
|-------------------------------------------------|-----------------------------------------------------------------------------------------|------------------------------------------------------------------------------|
| Streaming pre/post hash compare for paths       | standard pattern from `git`, `seqwish`, and `gfaffix` itself                            | high — the codebase already uses `FxHashMap`, `rayon`, and `Vec<u8>` walks   |
| In-memory pipeline entry point alongside string entry point | mirror of how `apply_graph_transforms` already toggles between in-memory and string for sort | high — same conditional pattern already present                              |
| Reusable scratch `String` buffer between rounds | mirror of `seqwish` and `wfmash` round-buffer patterns                                  | high — purely a Rust allocation pattern                                      |
| `criterion`-based benchmark harness             | standard Rust pattern, used in adjacent rust-bio crates                                 | medium — repo has no `benches/` today; adding one is mechanical              |
| POVU `from_graph` constructor (Option B)        | needs an upstream change to `vendor/povu` or a `pub(crate)` shim                        | **low** — flagged as risk #5; treated as a separate follow-up task           |

Each high-capacity import is implementable with the patterns already in
`src/resolution.rs`. The single low-capacity import (Option B) is
deliberately deferred to a separate `crush-perf` task and is not on the
critical path for the speed targets in § 2.2.

---

## 9. Open Questions Carried Forward

These are not blockers for the design but the implementer should keep
them in mind:

1. After the baseline lock-in, which `sweepga_kmer_frequency` default is
   correct? (Handoff open question 1.)
2. Is the C4 second-round rejection by old quality gating
   biologically/visually correct, or noise the user prefers gone?
   (Handoff open question 5.)
3. Should there be an opt-in `--max-round-score-growth` flag? (Research
   brief open question 6.) The design leaves this to `crush-perf`.

---

## Source Reference (changes this design proposes touch only these files)

| File                                          | Change                                                                                         |
|-----------------------------------------------|------------------------------------------------------------------------------------------------|
| `src/resolution.rs`                           | Promote `Graph` etc. to `pub(crate)`; add `resolve_graph_bubbles_inmemory`, `render_resolved_graph`; rewrite `path_sequences_equal` as streaming + parallel; add reusable POVU render buffer |
| `src/lib.rs`                                  | No change (still calls `resolve_gfa_bubbles` with a `String`)                                  |
| `src/main.rs`                                 | No change (CLI still calls `resolve_gfa_bubbles`)                                              |
| `benches/crush.rs` (new)                      | Criterion harness over the yeast-derived GFA                                                   |
| `tests/smoke/scenarios/crush_inmemory_matches_string.sh` (new) | Smoke scenario asserting byte-equality between the two entry points              |
| `tests/smoke/manifest.toml`                   | Add the new scenario, listing the relevant crush impl tasks as `owners`                        |
| `Cargo.toml`                                  | Add `criterion` to `[dev-dependencies]` if not already present                                 |
