# Diagnose residual underaligned C4 chunks

Task: `diagnose-residual-underaligned`

Date: 2026-06-05 UTC

Input graph:

```text
/home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/graphs/c4.k311.poa2kb.gfa
```

Experiment root:

```text
/home/erikg/impg/data/c4_residual_motif_polish_20260605T150000Z
```

Driver:

```text
scripts/c4-motif-local-polish.py
```

## Summary

The residual offshoots are real path-step/local-coverage events, not just a
PNG impression. The scan found 1,431 sparse offshoot candidates where <=5
paths leave high-coverage flanks, plus 53 direct self-loop repeat candidates.
The top offshoots line up with the graph-report long path jumps, for example
`10858+ -> 10861+ -> 10860+`, `10927+ -> 7245+ -> 10929+`,
`10921+ -> 10924+ -> 10923+`, and `7170+ -> 7173+ -> 10870+`.

For the five representative extracted chunks, lowering local seqwish
min-match to 79, 127, or 191 did not improve underalignment. It did not even
change the induced result versus 311: FastGA emitted a 0-byte raw PAF for
every tested chunk and every k value, so seqwish received no alignments and
emitted one disconnected segment per path. This makes `k=311` not the
proximate cause for these residual chunks. The proximate blocker is candidate
and evidence selection: the local windows are too short/tight for the current
FastGA -> seqwish path, and the existing POVU-boundary/direct POA paths do not
produce a path-preserving improvement.

Direct self-loop motifs are a separate, already-addressable class. The existing
`impg normalize-self-loops` pass is motif-local, independent of POVU flubble
boundaries, and preserves all path spellings. On the full C4 graph it removes
all 57 direct self-loop edges and all 252,715 adjacent same-node repeat steps,
but it is not a singleton-offshoot fix: full-graph singleton bp increases from
6,350 to 6,748 and path white-space p99 increases from 120,462 to 186,283.

## Artifacts

Machine-readable outputs:

```text
/home/erikg/impg/data/c4_residual_motif_polish_20260605T150000Z/discovered_chunks.tsv
/home/erikg/impg/data/c4_residual_motif_polish_20260605T150000Z/selected_chunks.tsv
/home/erikg/impg/data/c4_residual_motif_polish_20260605T150000Z/windows.tsv
/home/erikg/impg/data/c4_residual_motif_polish_20260605T150000Z/variant_results.tsv
```

Each selected chunk has an extracted local input GFA and FASTA under:

```text
/home/erikg/impg/data/c4_residual_motif_polish_20260605T150000Z/chunks/<chunk_id>/
```

The FASTA and local GFA records preserve the original path names. Full path
lists and path-step windows are in `windows.tsv`; only representative/example
names are shown below because several chunks have 459-465 local windows.

Uploaded full-graph self-loop prototype PNGs:

```text
https://hypervolu.me/~erik/impg/c4.k311.poa2kb.before-selfloop-normalize.Ygs.gfalook-m.png
https://hypervolu.me/~erik/impg/c4.k311.poa2kb.selfloop-normalized.Ygs.gfalook-m.png
```

## Selected chunks

Path-step coordinates are 0-based indices into the original GFA `P` record.
`core_steps` is the sparse/offending motif; anchors are the high-coverage
flanks used to extract same-anchor local windows across paths.

| chunk id | kind | example / involved path | path-step coord | anchors | core steps | local windows | traversal bp min/median/max | input segment bp | input singleton bp | input self-loops |
| --- | --- | --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `off_001_10861_10` | singleton offshoot | `NA20806#2#CM102496.1:31859718-32092404(+)` | 10-10 | `10858+` -> `10860+` | `10861+` | 465 | 449/449/449 | 450 | 1 | 0 |
| `off_002_7245_104` | singleton offshoot | `NA21110#2#CM089683.1:31841000-32073722(+)` | 104-104 | `10927+` -> `10929+` | `7245+` | 465 | 287/287/287 | 288 | 1 | 0 |
| `off_003_10924_99` | singleton offshoot | `NA20282#1#JBIRDT010000005.1:5167979-5453555(+)` | 99-99 | `10921+` -> `10923+` | `10924+` | 465 | 153/153/153 | 154 | 1 | 0 |
| `off_004_7173_21` | singleton offshoot | `NA19338#2#CM087775.1:31924004-32150386(+)` | 21-21 | `7170+` -> `10870+` | `7173+` | 465 | 123/123/123 | 124 | 1 | 0 |
| `loop_001_7068_464` | self-loop repeat | 459 paths; max run in `HG01358#2#CM089068.1:31846088-32072428(+)` | 464-492 | `7067+` -> `7069+` | 29 copies of `7068+` | 459 | 108/110/115 | 87 | 0 | 1 |

## Local variant results

All path-preserving variant GFAs were checked with:

```bash
/home/erikg/impg/target/release/examples/compare_gfa_paths <local-input.gfa> <variant.gfa>
```

For the 25 path-preserving non-input variant outputs, `compare_gfa_paths`
reported 0 missing paths, 0 extra paths, and 0 spelling mismatches. The five
`impg graph --gfa-engine poa` outputs were rejected as rewrite candidates:
they produced the right number of paths but not the original path names
(`missing_paths == observed_paths`, `extra_paths == observed_paths`).

| chunk id | chosen method | k / min-match | segment bp before -> after | singleton bp before -> after | self-loop count before -> after | visual underalignment improved? |
| --- | --- | --- | ---: | ---: | ---: | --- |
| `off_001_10861_10` | no accepted improving local rewrite; least-bad pass-through was `abpoa_crush` | n/a | 450 -> 450 | 1 -> 1 | 0 -> 0 | no |
| `off_002_7245_104` | no accepted improving local rewrite; least-bad pass-through was `abpoa_crush` | n/a | 288 -> 288 | 1 -> 1 | 0 -> 0 | no |
| `off_003_10924_99` | no accepted improving local rewrite; least-bad pass-through was `abpoa_crush` | n/a | 154 -> 154 | 1 -> 1 | 0 -> 0 | no |
| `off_004_7173_21` | no accepted improving local rewrite; least-bad pass-through was `abpoa_crush` | n/a | 124 -> 124 | 1 -> 1 | 0 -> 0 | no |
| `loop_001_7068_464` | `normalize-self-loops` local motif cleanup | n/a | 87 -> 291 | 0 -> 27 | 1 -> 0 | repeat loop yes; singleton/offshoot no |

Local seqwish k matrix details:

| chunk class | k=79 | k=127 | k=191 | k=311 | interpretation |
| --- | --- | --- | --- | --- | --- |
| four singleton offshoots | raw PAF 0 bytes, filtered 0 mappings, 0 links | same | same | same | Lower k cannot help because no alignments reach seqwish. |
| self-loop repeat | raw PAF 0 bytes, filtered 0 mappings, 0 links | same | same | same | Seqwish removes the loop only by exploding to one singleton segment per path. |

Representative log line pattern, repeated for every tested chunk and k:

```text
FastGA completed, converting .1aln to PAF
Conversion completed successfully, output size: 0 bytes
[sweepga] Summary: 0 -> 0 mappings (0.0% kept)
Compaction complete (465 nodes)
Links derived (0 links)
```

## Whole-graph self-loop prototype

Command:

```bash
impg normalize-self-loops \
  -g /home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/graphs/c4.k311.poa2kb.gfa \
  -o /home/erikg/impg/data/c4_residual_motif_polish_20260605T150000Z/whole_graph_selfloop_norm/c4.k311.poa2kb.selfloop-normalized.gfa \
  --report /home/erikg/impg/data/c4_residual_motif_polish_20260605T150000Z/whole_graph_selfloop_norm/selfloop-normalize-report.json
```

Validation:

```text
expected_paths      465
observed_paths      465
missing_paths       0
extra_paths         0
spelling_mismatches 0
```

Key graph-report deltas:

| metric | before | after |
| --- | ---: | ---: |
| segments | 8,160 | 8,350 |
| links | 10,928 | 11,289 |
| path steps | 3,374,388 | 3,121,673 |
| total segment bp | 251,692 | 254,167 |
| singleton bp | 6,350 | 6,748 |
| direct self-loop edges | 57 | 0 |
| adjacent same-node path steps | 252,715 | 0 |
| adjacent same-step path steps | 252,715 | 0 |
| self-loop repeat runs | 28,031 | 0 |
| local repeat-context nodes | 667 | 651 |
| POVU sites / leaves | 2,886 / 2,701 | 2,829 / 2,646 |
| path white-space p99 | 120,462 | 186,283 |

The self-loop prototype is therefore useful as a motif-local repeat cleanup,
but should not be treated as solving the residual singleton/offshoot rendering
issue.

## Cause assessment

For the concrete residual chunks tested here:

- `k=311` is not the direct blocker. Lowering the seqwish floor to 79, 127,
  or 191 produced identical failures because the local FastGA step emitted no
  PAF records.
- Candidate/evidence selection is the direct blocker. The current local windows
  around sparse offshoots are real, but the selected sequences and tight
  boundaries do not feed usable alignments into the existing FastGA -> seqwish
  induction path. Direct POA graph induction is not path-name preserving, and
  abPOA crush preserves path spellings but makes no useful change on these
  windows.
- Both still matter architecturally: k is relevant once the local candidate
  has alignable evidence, but these residual offshoots first need a candidate
  source that is independent of clean POVU boundaries and can use wider/flank
  aware sequence evidence or a direct path-preserving POASTA/abPOA builder.

## Patch plan

A safe implementation should add a motif-local candidate source rather than
lowering the global seqwish default.

1. Add a new candidate source and CLI mode.

   Files/functions:

   - `src/resolution.rs`: extend `ResolutionMethod` with a name such as
     `MotifLocal` or extend `CoverageMultiBubble` with `window-mode=motif`.
   - `src/resolution.rs`: update `ResolutionMethod::parse_name()` and
     `ResolutionMethod::method_name()`.
   - `src/main.rs`: update `parse_crush_stage()` help/validation strings and
     add options for `motif-max-sparse-paths`, `motif-min-flank-paths`,
     `motif-min-order-jump`, and `motif-max-window-bp`.

2. Discover candidates independent of POVU.

   Files/functions:

   - Add `discover_motif_local_candidates(graph: &Graph, config: &ResolutionConfig)`
     near the existing discovery path used by `discover_all_candidates()`.
   - Reuse the path-index and coverage data already built for
     `materialize_candidate_sequences()` and multi-level windows.
   - Detect:
     - direct self-loop repeat runs, delegating to
       `src/gfa_self_loops.rs::normalize_repeat_self_loops` when the target is
       only a repeat-unit loop;
     - sparse offshoot runs with coverage path count <=5 between flanks with
       high path support;
     - long path-step order/white-space jumps with support <=5.

3. Materialize anchor-pair windows across all spanning paths.

   Files/functions:

   - Add an anchor-pair materializer beside
     `build_multi_level_window_candidate()` and
     `materialize_candidate_sequences()` in `src/resolution.rs`.
   - Use original `source_path_name` / path names from the input graph, not
     synthetic FASTA IDs.
   - Use high-coverage left/right anchors and collect every path traversal
     between the anchors, with bounded flank bp and `max_window_bp`.

4. Route replacement by motif shape, not only by median length.

   Files/functions:

   - In `multi_level_window_replacement_method()` or the new motif builder,
     route sub-2 kb motifs to existing path-preserving POASTA/abPOA replacement
     builders, not `impg graph --gfa-engine poa`.
   - For offshoot motifs, run direct POASTA/abPOA first; use FastGA/seqwish only
     if raw PAF count is nonzero. If raw PAF is zero, reject that variant
     immediately and log the evidence gap.
   - Keep the acceptance gate exact-path preserving and add a local objective
     that rejects singleton-bp explosions even if self-loop count decreases.

5. Tests.

   Add focused tests:

   - `src/gfa_self_loops.rs`: existing normalizer tests already cover exact
     path preservation; add a graph-report regression asserting direct-loop
     metrics go to zero for a C4-like repeat run.
   - `src/resolution.rs`: synthetic sparse-offshoot graph where one path has
     `left -> private -> right` and the majority has `left -> shared -> right`;
     assert motif discovery returns one candidate with original path names.
   - `tests/test_crush_integration.rs`: a small extracted C4 offshoot fixture
     derived from `off_004_7173_21` or a synthetic equivalent; assert local
     FastGA/seqwish zero-PAF variants are rejected and do not replace the
     graph.
   - CLI parser test in `src/main.rs` for the new motif mode/options.

Validation commands after implementation:

```bash
source ./env.sh && cargo test --lib gfa_self_loops
source ./env.sh && cargo test --lib motif_local
source ./env.sh && cargo test --test test_crush_integration c4_motif
source ./env.sh && cargo test --bin impg motif
source ./env.sh && cargo build
```

This task did not implement the new Rust candidate source because the
experiments show the correct behavior depends on a new discovery/materialize
path and an acceptance objective, not a small threshold edit. The safe code
landed here is the reusable experiment driver plus this patch plan.

## Validation notes

- `compare_gfa_paths` passed for every path-preserving local variant and for
  the whole-graph self-loop-normalized prototype.
- Direct `poa` graph outputs were not accepted because `compare_gfa_paths`
  failed due path-name replacement.
- Full PNG renders were produced with `gfasort -p Ygs` and `gfalook -m`, then
  uploaded only to `erik@hypervolu.me:www/impg/`.
- A fresh local `cargo build --release --example compare_gfa_paths` failed in
  the vendored `wfmash-rs` C++ build because `htslib/faidx.h` was not visible.
  The experiment therefore used the existing documented binary:
  `/home/erikg/impg/target/release/examples/compare_gfa_paths`.
- The documented environment was also tried:
  `source ./env.sh && cargo build` and
  `source ./env.sh && cargo test --lib gfa_self_loops`. Both failed before
  task-specific Rust code ran because this worktree is missing the vendored
  syng C sources required by the build script, for example
  `vendor/syng/syngbwt3.c`, `vendor/syng/rskip.c`, and
  `vendor/syng/impg_syng_helpers.c`.
- `python3 -m py_compile scripts/c4-motif-local-polish.py` passed.
