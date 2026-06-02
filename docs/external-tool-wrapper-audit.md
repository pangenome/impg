# External Tool Wrapper Semantics Audit

Task: `audit-external-tool`

Date: 2026-05-28

## Invariant

External tool wrappers in `impg` must not silently change utility semantics.

Allowed wrapper behavior:

- marshal inputs and outputs;
- invoke external tools or library utilities with documented parameters;
- parse tool outputs;
- validate correctness invariants, such as path preservation and parse validity;
- log diagnostics.

Disallowed hidden behavior:

- filtering alignments or records beyond explicit user-selected tool parameters;
- substituting fallback outputs when a tool or wrapper stage fails;
- rewriting tool output in a way that changes graph, path, or alignment semantics;
- imposing hidden minimum length, exact-run, identity, path, or graph complexity thresholds;
- dropping or mutating sequence/path names except for reversible, documented marshalling.

## Decisions

This audit found and removed hidden semantic changes in the alignment and local graph-induction wrappers:

- `src/syng_graph.rs:796` no longer re-adds raw PAF records after the explicit SweepGA PAF filter drops them. Regression coverage is in `src/syng_graph.rs:1409`.
- `src/resolution.rs:4841` no longer applies an `impg`-side 200 bp wfmash input floor. wfmash/MashMap now owns minimum-length semantics and should fail explicitly if a non-empty input is unsupported.
- `src/resolution.rs:4821` no longer rejects otherwise valid SweepGA/seqwish replacements solely because the seqwish output has no shared replacement path segments. It logs this as a diagnostic.
- `src/resolution.rs:4402` no longer substitutes an unpolished replacement graph when POA/POASTA/smooth polishing fails.
- `src/resolution.rs:4453` and `src/resolution.rs:4553` no longer fall back to direct POASTA/SPOA when chain-POVU smoothing or AllWave pair scheduling fails.
- `src/smooth.rs:424` no longer substitutes an unsmoothed subgraph when SPOA/POASTA block smoothing fails. Gap blocks are still passthrough by construction because they are the coverage complement between selected smoothable blocks.
- `src/smooth.rs:1051` and `src/smooth.rs:1193` no longer fall back from explicit POVU reference hints to the first path when the hint does not match.
- `src/commands/align.rs:1028` no longer returns a partial PAF if a generic per-pair alignment invocation fails.

## Wrapper Table

| Integration point | Code | External/local utility | Classification | Decision |
| --- | --- | --- | --- | --- |
| `graph` FASTA combine, aligner invocation, PAF filtering | `src/commands/graph.rs:150`, `src/commands/graph.rs:731` | SweepGA, wfmash, FastGA | marshal/invoke/parse plus explicit user parameters | Retained. FASTA input is combined, optional `input_paf` is copied, aligner parameters come from `GraphBuildConfig`, and filtering is delegated to the shared SweepGA PAF filter. |
| Shared graph PAF filter | `src/commands/graph.rs:675` | SweepGA `PafFilter` | explicit user parameter | Retained. `no_filter`, mapping mode, scaffold mode, identity, and minimum map length are user-visible config. No post-filter rescue remains in local seqwish induction. |
| seqwish graph induction tail | `src/commands/graph.rs:167` | seqwish crate | marshal/invoke/parse plus correctness validation | Retained. `min_match_len`, sparse factor, transitive closure settings, disk-backed trees, and threads come from `GraphBuildConfig`. The global lock is a correctness guard for seqwish temp-file bookkeeping, not a semantic transform. |
| In-memory SweepGA alignment API | `src/commands/align.rs:801` | SweepGA, wfmash, FastGA | marshal/invoke/parse plus explicit user parameters | Retained with one fix. The generic per-pair path now errors if an alignment invocation fails instead of silently returning a partial PAF. |
| wfmash job list emission | `src/commands/align.rs:283` | SweepGA joblist/wfmash | marshal/invoke | Retained. Emits documented wfmash PanSN job commands through SweepGA's joblist helper. |
| FastGA job list emission | `src/commands/align.rs:355` | FastGA | marshal/invoke | Retained. Emits one command per unique FASTA file pair. File-stem output names are command marshalling, not sequence/path name mutation. |
| Aligner factory | `src/commands/mod.rs:49` | SweepGA wfmash/FastGA integrations | explicit user parameter | Retained. It passes documented backend options through to SweepGA integrations. Adaptive wfmash parameters are owned by SweepGA/wfmash integration rather than hidden downstream filtering. |
| Local PAF to seqwish helper | `src/syng_graph.rs:796` | seqwish tail, SweepGA PAF filter | marshal/invoke/parse plus explicit parameter | Fixed. It writes exact FASTA names, applies only the explicit shared PAF filter, optionally lowers `min_match_len` only when `adaptive_min_match_len` is true, and then invokes the shared seqwish tail. |
| Syng-native full-pair BiWFA PAF | `src/syng_graph.rs:335`, `src/syng_graph.rs:410` | WFA2/BiWFA | correctness validation | Retained. Invalid CIGAR consumption is rejected because it would break seqwish coordinates. The output is a PAF representation of the BiWFA result. |
| Syng-native sparse BiWFA PAF | `src/syng_graph.rs:430` | SweepGA kNN pair selection, WFA2/BiWFA | explicit user parameter | Retained. Pair sampling parameters are explicit. Pairwise CIGAR validation is correctness validation. |
| Anchor-seeded syng BiWFA PAF | `src/syng_graph.rs:973` | syng query, SweepGA anchor scaffold filter, WFA2/BiWFA | marshal/invoke/parse plus correctness validation | Retained. Falling back from anchor-seeded gap alignment to full-pair BiWFA affects performance and coverage of the same requested pair, not filtering or quality-based substitution. It is logged. |
| AllWave replacement PAF to seqwish | `src/resolution.rs:4336`, `src/resolution.rs:4553` | AllWave, seqwish | explicit user parameter | Fixed. Pair schedule parameters remain explicit. Zero non-empty traversal or zero selected pair cases now error instead of substituting direct SPOA. |
| SweepGA/seqwish replacement | `src/resolution.rs:4696` | SweepGA, wfmash/FastGA, seqwish | explicit user parameter | Fixed. Top-flubble mode trusts the explicit SweepGA filter and disables the second seqwish-tail PAF filter. The hidden wfmash input floor was removed. No-shared replacement output is diagnostic-only. |
| Replacement seqwish filter config | `src/resolution.rs:4218` | seqwish tail, SweepGA PAF filter | explicit user parameter | Retained. Local replacement defaults now keep `min_match_len=1`; `off`, `adaptive`, and fixed floors are explicit via CLI and config. `no_filter=true` disables downstream exact-run, identity, and PAF filters. |
| Pairwise replacement polish | `src/resolution.rs:4402` | SPOA, POASTA, smooth | correctness validation | Fixed. A polish failure now fails the replacement path rather than substituting the unpolished graph. |
| Chain-POVU smooth to POASTA replacement | `src/resolution.rs:4453` | smoothxg-style smoothing, POASTA | correctness validation | Fixed. Empty traversal or smooth failure now errors instead of falling back to direct POASTA. |
| POASTA replacement/block wrappers | `src/resolution.rs:4943`, `src/resolution.rs:5021`, `src/smooth.rs:2016` | POASTA | marshal/invoke/parse plus correctness validation | Retained. POASTA output must parse and preserve replacement paths. Smooth-block callers now error on failure instead of substituting a passthrough graph. |
| SPOA graph generation | `src/graph.rs:310`, `src/resolution.rs:5956`, `src/smooth.rs:2043` | SPOA | marshal/invoke/parse | Retained with smoothing fix. SPOA is used to build the requested POA graph. Smoothing no longer hides SPOA failures behind unsmoothed output. |
| smoothxg-style block smoothing | `src/smooth.rs:198`, `src/smooth.rs:424` | SPOA/POASTA plus lacing | marshal/invoke/parse plus correctness validation | Fixed. Gap passthrough remains explicit structural marshalling. Smoothable-block failures now error. |
| POVU flubble and neighbor block placement | `src/smooth.rs:1051`, `src/smooth.rs:1193` | POVU native Rust API | explicit user parameter | Fixed. Explicit reference hints no longer fall back to the first graph path if POVU rejects the hint. |
| gfaffix normalization | `src/graph.rs:838` | `gfaffix` binary | marshal/invoke/parse | Retained. The wrapper searches only beside the running executable to avoid version drift and errors if the binary is missing or exits non-zero. |
| gfasort and unchop helpers | `src/graph.rs:732`, `src/graph.rs:763` | gfasort crate | explicit user parameter plus correctness validation | Retained. The requested pipeline is validated, then gfasort is invoked in-process. Empty pipeline and trivial graph passthrough are identity cases. |
| Graph engine dispatch transforms | `src/lib.rs:772` | crush, smooth, gfaffix, gfasort | explicit user parameter | Retained. Transforms are applied only when selected by `EngineOpts` and return errors from their wrappers. |
| Query local graph render bundle | `src/commands/render.rs:133`, `src/commands/render.rs:321` | syng, seqwish, pggb, POA | marshal/invoke/parse | Retained. Builds a rendered FASTA, invokes the selected local graph engine, then records translation metadata. Unsupported engines error. |
| gfalook image rendering | `src/main.rs:10067` | `gfalook` binary | marshal/invoke/parse | Retained. Skips only empty graphs with no segment records. Missing binary and non-zero status are errors. |
| Syng direct region GFA | `src/lib.rs:593`, `src/lib.rs:970` | syng index, syng2gfa, bluntg mode | explicit user parameter | Retained. `syng:raw` preserves syng overlaps, `syng:blunt` requests blunt graph materialization, and syncmer parameter assertions are explicit. |
| syng2gfa materialization | `src/commands/syng2gfa.rs:49`, `src/commands/syng2gfa.rs:1535` | syng, pangenome/bluntg-equivalent native path | explicit user parameter | Retained. Raw/blunt mode, frequency mask, shared-run mask, local-repeat clone thresholds, and N-run cutting are explicit stages or defaults. `nomask` disables masking. |
| syng2gfa frequency/shared-run masks | `src/commands/syng2gfa.rs:992`, `src/commands/syng2gfa.rs:1183`, `src/main.rs:3466` | syng2gfa local graph materializer | explicit user parameter | Retained. Query syng GFA defaults to `local_default()` masking, and `mask`, `nomask`, `cut-ns`, and parameterized stages are parsed in the GFA engine string. This is documented graph-generation policy, not a wrapper around seqwish induction. |
| syng render bundle | `src/commands/render.rs:58` | syng index and syng2gfa | marshal/invoke/parse | Retained. Fetches source intervals, builds a regional syng index, optionally emits GFA, and writes reversible translation tables. |
| Syng-backed query wrapper | `src/lib.rs:282` | syng index | explicit user parameter plus marshal/parse | Retained. Raw syncmer intervals flow through unless `with_chain_filter` is selected. Chain count and span fraction are explicit syng query settings. |
| Syng transitive anchor chaining | `src/syng_transitive.rs:173` | syng query, SweepGA scaffold filter, BiWFA refinement | explicit user parameter plus correctness validation | Retained. Anchor count, span fraction, extension budget, depth, and seed filters are explicit syng query settings. |

## Regression Coverage

Added or updated regression coverage:

- `src/syng_graph.rs:1409` verifies that `build_gfa_from_paf_and_sequences` does not add raw PAF records back after an explicit `min_map_length` filter removes them.
- `src/smooth.rs:2606` verifies explicit POVU flubble reference hints now error when the hint misses graph paths.
- `src/resolution.rs:6738` verifies wfmash and FastGA local replacement input floors do not hide backend minimum-length semantics.
- `tests/test_crush_integration.rs:465` now exercises local C4 fragment seqwish induction with default local graph-induction semantics: `no_filter=true`, `min_match_len=1`, and no adaptive hidden exact-run floor.

Focused validation run before the full test pass:

- `cargo test syng_graph::tests::seqwish_tail_does_not_rescue_records_dropped_by_explicit_filter`
- `cargo test smooth::tests::test_flubble_guided_blocks_errors_when_reference_hint_misses`
- `cargo test resolution::tests::wfmash_aligner_input_floor_does_not_hide_tool_length_semantics`
- `cargo test --test test_crush_integration c4_fragment_seqwish_regressions_induce_shared_graphs -- --test-threads=1 --nocapture`
- `cargo test --test test_crush_integration c4_fragment -- --test-threads=1 --nocapture`
- `cargo test --test test_crush_integration c4_top_flubble_seqwish_indexes_observed_exact_runs -- --test-threads=1 --nocapture`
- `cargo build`
- `cargo test --all -- --test-threads=1`
- `cargo install --path .`

The representative C4 fragment validation preserved all fixture paths and still induced shared graph structure without reintroducing an unfolded replacement regression:

| Fixture | Paths | Segments | Segment bp | Shared segments | Duplicate segment sequences |
| --- | ---: | ---: | ---: | ---: | ---: |
| `easy_shared_flank` | 6 | 4 | 227 | 4 | 0 |
| `bounded_multi_bubble` | 12 | 9 | 1768 | 1 | 6 |
| `unfolded_minrun` | 5 | 22 | 593 | 16 | 6 |
| `short_floor` | 9 | 10 | 450 | 9 | 1 |
| `duplicated_repeat` | 12 | 10 | 225 | 9 | 2 |

## Remaining Guardrails

- Historical evaluation documents may still mention removed behavior as historical context. Current wrapper policy is this document's invariant.
- The shared `ImpgIndex` trait still returns vectors rather than `Result`, so syng query wrappers that implement the trait cannot propagate every query-layer error without a larger trait migration. This audit did not change that trait because the failing C4 graph-induction path is below the graph wrapper and now has explicit error propagation.
