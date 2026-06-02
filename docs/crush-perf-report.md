# Crush Performance Report

**Task:** `crush-perf`
**Branch:** `wg/agent-20/crush-perf` (parent: `crush-integrate` @ 6e0dbe7)
**Date:** 2026-05-24
**Goal:** Profile the integrated crush implementation against a real, large
GFA input and optimize the per-round infrastructure (parse / render / validate
— the non-alignment hot path called out in `docs/crush-design.md` §2.1).

---

## 1. Hardware Context

| Field | Value |
|-------|-------|
| Machine | Linux 6.8.0-106-generic |
| CPU | AMD EPYC 7713 64-Core Processor (× 2 sockets) |
| Logical CPUs | 256 (128 cores × 2 threads/core) |
| Build | `cargo build --release` (default `rustc 1.x`) |
| Workload binary | `target/release/impg` |

`perf_event_paranoid=4` on this host, so I used the criterion harness in
`benches/crush.rs` plus in-process `Instant::now()` instrumentation (already
present in `src/resolution.rs`) rather than `perf` / `samply`.

---

## 2. Workload

**Primary input:**
`/home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.blunt.gfa`

This is the C4 blunt graph from the design doc's reference workload
(`docs/crush-design.md` §5.2 line 359).  Shape:

| Bytes  | Segments | Paths | Links | Total path steps |
|-------:|---------:|------:|------:|-----------------:|
| 45 MiB | 18 048   | 465   | 20 944 | 4 591 857       |

It is `.gitignore`d (the developer-machine baseline; not a CI fixture).  Override
the bench input with `CRUSH_BENCH_GFA=/abs/path/to/blunt.gfa cargo bench --bench crush`.

The bench harness also runs cleanly on the 117 MiB / 113 716-segment AMY1A
graph — included in §4.3 below.

---

## 3. Bench Harness

`benches/crush.rs` (new) is a `criterion` harness with three benches.  Run
with:

```
cargo bench --bench crush
```

The benches call three `#[doc(hidden)] pub` bench helpers added to
`src/resolution.rs`:

| Helper | What it measures |
|--------|------------------|
| `bench_parse_gfa(gfa)` | `parse_gfa` only — the entry-side cost paid once per `resolve_gfa_bubbles` call |
| `bench_parse_render_gfa(gfa)` | `parse_gfa` + `render_graph` round-trip — the per-discovery-round render cost (see `src/resolution.rs:1304`) |
| `bench_validate_streaming_self(gfa)` | 2× `parse_gfa` + `path_sequences_equal_streaming` — the per-round post-replacement validator (see `src/resolution.rs:1597`) |

These three primitives are the work the crush algorithm pays per round
*outside* alignment.  Alignment time (POA / POASTA / SweepGA) is bounded by
the candidate count, dominates the wall on C4-sized inputs, and is explicitly
out of scope per `docs/crush-design.md` §7 (see §6 below for the in-process
end-to-end timing that includes it).

---

## 4. Results

### 4.1 Baseline vs. optimized (C4A.blunt.gfa, single-input, 10-sample criterion)

| Bench                        | Baseline | Optimized | Speedup |
|------------------------------|---------:|----------:|--------:|
| `parse_gfa`                  | 165.2 ms | 50.2 ms   | **3.3×** |
| `parse + render` (round-trip)| 1 723.7 ms | 317.4 ms | **5.4×** |
| `validate streaming self`    | 433.5 ms | 136.4 ms  | **3.2×** |

Derived per-phase times (round-trip minus parse, validate minus 2× parse):

| Phase                                  | Baseline | Optimized | Speedup |
|----------------------------------------|---------:|----------:|--------:|
| `render_graph` only                    | ~1 559 ms | ~267 ms  | **5.8×** |
| `path_sequences_equal_streaming` only  | ~103 ms  | ~36 ms    | ~2.9× (already parallel before) |

### 4.2 Per-round end-to-end on C4 (real `impg crush` invocation)

Run (optimized binary, `--max-iterations 1 --method poasta --max-traversal-len 0
--max-median-traversal-len 0` to bound the round to discovery + selection-guard
only — no alignment subprocess; full alignment dominates the wall on default
budgets and is out of scope here):

```
$ /usr/bin/time -v target/release/impg crush \
    -g .../C4A.blunt.gfa --max-iterations 1 \
    --method poasta --max-traversal-len 0 --max-median-traversal-len 0 \
    -o /tmp/c4_one_round.gfa
```

In-process log from the optimized run:

```
crush: parsed input GFA into 18048 segment(s), 465 path(s) in 84.53 ms
crush discovery detail: render 251.29 ms, povu-parse 1.27 s,
                        povu-decompose 157.90 ms, id-map 720.88 µs,
                        path-index 150.95 ms, candidate-build 177.78 ms
crush round 1: 2441 POVU site(s), 0 selected, 2437 selection-guarded in 2.78s
Elapsed (wall clock) time: 0:03.37
Maximum resident set size: 1.05 GB
```

Comparing the rendered-graph-prep portion against the baseline cost predicted
by §4.1:

| Phase                       | Baseline (predicted) | Optimized (measured) | Saved/round |
|-----------------------------|---------------------:|---------------------:|------------:|
| Initial `parse_gfa`         | ~165 ms              | ~85 ms               | ~80 ms once |
| Per-round `render_graph`    | ~1 559 ms            | ~251 ms              | ~1.3 s      |
| POVU re-parse (povu-internal — **NOT my target**) | ~1.27 s | ~1.27 s | — |

On a multi-round crush (the C4 default is ~5 rounds), this is **~6.5 s of
infrastructure cost removed** before alignment even starts.

### 4.3 Scaling check on AMY1A (117 MiB, 113 716 segments, 466 paths)

`CRUSH_BENCH_GFA=/home/erikg/impg/data/cosigt_svr_graph_direct_20260520135313/AMY1A.gfa
cargo bench --bench crush -- --quick` on the optimized build:

| Bench                       | Time   | Throughput |
|-----------------------------|-------:|-----------:|
| `parse_gfa`                 | 158 ms | 738 MiB/s  |
| `parse + render`            | 997 ms | 117 MiB/s  |
| `validate streaming self`   | 500 ms | 233 MiB/s  |

Per-byte throughput on the bigger graph (more segments per byte → more
small-string work) is in line with C4: parse holds ≥ 700 MiB/s, validation
holds ≥ 230 MiB/s, render holds ≥ 100 MiB/s.

---

## 5. Top Hot Spots Identified & Fixes Applied

The three top hot spots in the optimizable per-round path (alignment
explicitly out of scope per `docs/crush-design.md` §7):

### 5.1 `render_rewritten_graph` — String formatting + UTF-8 validation
**Old wall share (parse+render bench):** ~90 % (1 559 / 1 724 ms).

The old renderer (`src/resolution.rs:2973-3088` pre-patch) had four
allocation-heavy hot spots:
1. Per-S-line `format!("S\t{}\t{}\n", id, String::from_utf8_lossy(&seq))` —
   two allocations and one UTF-8 validation pass per segment.  On C4 that is
   45 MiB of sequence bytes validated and reformatted into a fresh `String`.
2. Per-L-line `out.push_str(&format!(...))` — same pattern × 20 944 edges.
3. `BTreeSet<(String, bool, String, bool)>` for edge dedup — String clones
   into the set plus tree rebalances.
4. Per-P-line `steps.iter().filter_map(...).collect::<Vec<_>>().join(",")` —
   one fresh `String` per step plus a `Vec<String>` collect.

**Fix (`src/resolution.rs:3107` after-patch):**
- Build into a `Vec<u8>` and `String::from_utf8_unchecked` at the end (all
  emitted bytes are ASCII by construction).  Skips UTF-8 validation on the
  ~45 MiB of sequence.
- Replace `BTreeSet<String, ...>` with `Vec<(u32, bool, u32, bool)>` (Copy)
  + `sort_unstable_by` + `FxHashSet` for dedup.
- Resolve OutNode → emitted ID via `Vec<u32>` keyed by original segment
  index (O(1)) for the common no-replacements path, falling back to a
  small FxHashMap for replacement nodes only.
- Parallel P-line render via `out_paths.par_iter()` into per-path
  `Vec<u8>` chunks (concatenated afterwards) — each path is independent.
- Pre-size `out` with a tight capacity estimate so we avoid reallocation
  during the 45 MiB write.

**After:** render contributes ~267 ms of the 317 ms `parse + render` bench
(was ~1 559 ms).  **5.8× faster.**

### 5.2 `parse_gfa` — Per-line `Vec<&str>` collect + serial step lookup
**Old wall share (per-call parse cost):** ~100 % (165 ms).

The old parser (`src/resolution.rs:703-789` pre-patch) had two cost centers:
1. `line.split('\t').collect::<Vec<&str>>()` — one allocation per line.
   C4 has 18 048 S + 20 944 L + 465 P = ~39 k allocations just for the field
   vectors.
2. Single-threaded path-step parsing: 4 591 857 total path steps, each one
   doing `parse_p_step` + `id_to_idx.get(id)` — a hash lookup per step.

**Fix (`src/resolution.rs:703-905` after-patch):**
- Skip the `collect::<Vec<_>>()`: walk fields via `splitn(...)` iterators
  for the small-arity cases (S, P, W) and a hand-rolled `nth_tab_field`
  for the single L-line overlap check.
- Three-pass parser: classify lines first, then build segments + id_to_idx
  serially (small relative to step parsing), then **parse all P and W
  lines in parallel** via `rayon::par_iter().with_min_len(8)`.  The
  read-only id_to_idx is shared; each path is independent.

**After:** parse contributes 50 ms (was 165 ms).  **3.3× faster.**
At 256 threads on C4's 465 paths, the parallelism is path-limited but
still gives ~3× because the per-step hash lookup is the dominant cost and
it amortizes across the 4.6 M total steps.

### 5.3 `render_rewritten_graph::ordered_nodes` — `FxHashSet<OutNode>` over 4.6 M steps
**Old wall share (inside render):** ~30 ms.

The old renderer dedup'd nodes-in-traversal-order via `FxHashSet<OutNode>`,
where `OutNode` is an enum (24 bytes when boxed for the hash).  4.6 M
inserts on the C4 graph is the third hot loop after sequence memcpy and
edge collection.

**Fix:** for the common case of no replacement graphs (which is the typical
final-render call site and the bench measurement), use a
`Vec<bool>` keyed by original segment index.  Replacement nodes go through
a small `FxHashSet<(usize, usize)>`.  Removes ~30 ms from the C4 render
without affecting the small-input call sites.

### 5.4 What was already optimal (left as-is per profile)

| Function                              | Status |
|---------------------------------------|--------|
| `path_sequences_equal_streaming`      | Already parallel + zero-alloc streaming (per design §3.4); 36 ms on C4 → no fix needed |
| POVU's internal `NativeGfa::parse`    | Not our code (in `vendor/povu`); 1.27 s share is the dominant *remaining* per-round cost.  Design §3.3 Option B (in-memory POVU constructor) is filed as a follow-up — out of scope for `crush-perf`. |
| `path_step_index` (per-path indices)  | Already `rayon::par_iter` (`src/resolution.rs:1337`); 150 ms on C4 |
| `materialize_candidate_sequences`     | Already parallel (`src/resolution.rs:1323`); algorithm-bound |

---

## 6. Correctness

The optimized code preserves the same byte-level behavior:

- `cargo test --release --all` — **375 / 375 pass**, no regressions
  (266 lib + 56 + 7 + 3 + 25 + 8 + 10 integration tests, 0 failures).
- Round-trip idempotence verified on C4A.blunt.gfa:
  `parse → render → parse → render` produces byte-identical output to
  `parse → render` (47 186 345 bytes both times), and all 465 path
  sequences are preserved.
- The existing test
  `inmemory_and_string_entry_points_produce_identical_resolved_graph`
  (`src/resolution.rs:4193`) asserts byte-equality between the two
  resolver entry points on the small_insertion fixture and still passes.
- `streaming_validator_matches_materializing_validator_on_known_inputs`
  (`src/resolution.rs:4205`) still passes against the materializing oracle.

No correctness checks were removed; both validation gates from
`docs/crush-design.md` §4.1 (`validate_replacement_paths` per bubble,
`path_sequences_equal_streaming` per round) are untouched.

---

## 7. Reproducing the bench

```
# baseline (this commit): single-input C4 bench
cargo bench --bench crush

# scaling check: AMY1A
CRUSH_BENCH_GFA=/abs/path/to/AMY1A.gfa cargo bench --bench crush -- --quick

# real end-to-end on C4 (discovery-only)
/usr/bin/time -v target/release/impg crush \
    -g .../C4A.blunt.gfa --max-iterations 1 \
    --method poasta --max-traversal-len 0 --max-median-traversal-len 0 \
    -o /tmp/c4_one_round.gfa
```

The bench input defaults to `data/c4_crush_eval_20260523T140141Z/C4A.blunt.gfa`
(absolute fallback `/home/erikg/impg/data/...`); set `CRUSH_BENCH_GFA` to
override.  The bench skips with a friendly message if the input is missing
so it stays usable on machines without the developer-machine fixtures.

---

## 8. Out of Scope / Open Follow-ups

Filed as known remaining work for the next `crush-perf` cycle, *not* part of
this report:

1. **POVU in-memory constructor** (`docs/crush-design.md` §3.3 Option B).
   POVU re-parses our rendered GFA every round; on C4 that is 1.27 s/round
   and now the single largest remaining infrastructure cost.  Requires an
   upstream change to `vendor/povu` to accept our in-memory `Graph`.
2. **Reusable render buffer** (`docs/crush-design.md` §3.3 Option A).
   Each round currently allocates a fresh ~45 MiB `Vec<u8>` for the
   rendered graph.  Cheap to add (carry a `Vec<u8>` field on the resolver
   state, `.clear()` between rounds) but the absolute saving is small
   compared to the wins in §5 — defer until profiling justifies it.
3. **Alignment** (`docs/crush-design.md` §7).  SweepGA / FastGA / seqwish /
   POA still dominate the wall on default budgets; speeding them up is a
   separate scope from this task.
