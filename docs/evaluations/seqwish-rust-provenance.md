# Rust seqwish provenance in IMPG

Task: `trace-rust-seqwish`

## Verdict

IMPG uses Rust seqwish as a direct Cargo git dependency. The seqwish core is
not implemented inside IMPG, not vendored as an IMPG crate/submodule, and not
invoked by shelling out to an external `seqwish` binary.

The dependency is declared in the root manifest:

- `Cargo.toml:74`: `seqwish = { git = "https://github.com/pangenome/seqwish", branch = "rust-2" }`
- `Cargo.lock:2897-2899`: package `seqwish` version `0.1.3`, source
  `git+https://github.com/pangenome/seqwish?branch=rust-2#15aee55987115bf82b1e0769e52543a265ecfb1a`

Resolved provenance:

- Repository URL: `https://github.com/pangenome/seqwish`
- Branch: `rust-2`
- Commit/revision used by this checkout:
  `15aee55987115bf82b1e0769e52543a265ecfb1a`
- Cargo source: git dependency, not a path dependency.
- Local Cargo cache observed during this audit:
  `/home/erikg/.cargo/git/checkouts/seqwish-be7416204a0fe48d/15aee55`

Checks:

- `cargo tree -i seqwish --locked` reports
  `seqwish v0.1.3 (https://github.com/pangenome/seqwish?branch=rust-2#15aee559) -> impg v0.4.1`.
- `git ls-remote https://github.com/pangenome/seqwish refs/heads/rust-2`
  returned `15aee55987115bf82b1e0769e52543a265ecfb1a` for `refs/heads/rust-2`
  on 2026-06-08.
- `git submodule status` lists `vendor/gfaffix` and `vendor/syng`, not
  `seqwish`.
- A Rust-source search found `seqwish::` calls only in
  `src/commands/graph.rs` and `src/main.rs`; there is no
  `Command::new("seqwish")` shell-out in the Rust implementation.

## What Code Is Local

Because IMPG does use an upstream seqwish branch, the "not upstream" case does
not apply. The local code is orchestration around the upstream crate:

- Alignment generation and filtering before seqwish: SweepGA/wfmash/FastGA for
  the normal graph path, syng/BiWFA for syng-native regional paths, and
  AllWave/SweepGA/wfmash for crush replacement paths.
- Temporary-file setup, debug artifact writing, and a process-wide lock around
  seqwish induction because the seqwish crate uses process-global tempfile
  bookkeeping.
- Empty-PAF and tiny-single-PAF fallbacks that emit source-sequence GFA instead
  of running transitive closure.
- Post-processing after seqwish: unchop, gfaffix normalization, gfasort, and
  optional smoothxg-style smoothing for the `pggb` engine.

The seqwish algorithmic core itself comes from the dependency. IMPG imports and
calls the upstream crate's library modules in `src/commands/graph.rs:29-35`:
`alignments::unpack_paf_alignments`, `seqindex::SeqIndex`,
`transclosure::compute_transitive_closures`, `compact::compact_nodes`,
`links::derive_links`, and `gfa::emit_gfa`. The upstream crate's local cached
`/home/erikg/.cargo/git/checkouts/seqwish-be7416204a0fe48d/15aee55/src/lib.rs:1-16`
describes the Rust seqwish stages as sequence indexing, PAF alignment
processing, transitive closure, node compaction, link derivation, and GFA
emission; IMPG's `induce_graph_from_alignment` follows those same stages at
`src/commands/graph.rs:194-445`.

## IMPG Files And Commands Using The Rust seqwish Path

Cargo/provenance:

- `Cargo.toml:69-74` declares graph-building dependencies and the seqwish git
  dependency.
- `Cargo.lock:2897-2917` pins `seqwish` version/source and its dependency set.

Direct seqwish crate integration:

- `src/commands/graph.rs:1-4` defines the graph-building module as
  SweepGA plus seqwish integration.
- `src/commands/graph.rs:29-35` imports seqwish crate APIs directly.
- `src/commands/graph.rs:156-164` runs alignment/filtering and then calls
  `induce_graph_from_alignment`.
- `src/commands/graph.rs:166-178` documents `induce_graph_from_alignment` as
  the post-alignment seqwish tail.
- `src/commands/graph.rs:178-184` serializes in-process seqwish inductions with
  `SEQWISH_INDUCTION_LOCK`.
- `src/commands/graph.rs:202-253` builds `SeqIndex`, creates seqwish temp
  files, and calls `unpack_paf_alignments`.
- `src/commands/graph.rs:314-325` calls
  `compute_transitive_closures`.
- `src/commands/graph.rs:343-413` calls `compact_nodes`, `derive_links`, and
  `emit_gfa`.
- `src/commands/graph.rs:443-445` cleans up seqwish temp files.
- `src/commands/graph.rs:513-532` implements `run_graph_build`, the plain
  `seqwish` graph engine entry point.
- `src/commands/graph.rs:606-640` implements `run_graph_build_pggb`, which
  first runs the seqwish pipeline and then smooths.
- `src/commands/graph.rs:941-1005` defines `AlignmentResult` and the shared
  generated-PAF filter used before seqwish.
- `src/commands/graph.rs:1306-1330` starts the partitioned graph pipeline,
  which aligns, partitions, then dispatches per-partition engines.

CLI wiring:

- `src/main.rs:1980-2003` configures seqwish's internal tempfile directory via
  `seqwish::tempfile::set_dir`.
- `src/main.rs:2056-2083` defines `SeqwishOpts` shared by `graph`, `query`,
  and `partition` commands.
- `src/main.rs:2123-2148` includes `SeqwishOpts` in `EngineCliOpts`.
- `src/main.rs:3798-3801` parses `--gfa-engine seqwish` to
  `GfaEngine::Seqwish`; `pggb` is parsed to `GfaEngine::Pggb`.
- `src/main.rs:4024-4050` builds `EngineOpts` and forwards seqwish options.
- `src/main.rs:5426-5443` declares `impg graph`, described as
  SweepGA plus seqwish.
- `src/main.rs:8732-8753` routes partitioned `impg graph --gfa-engine
  seqwish:WINDOW`/`pggb:WINDOW` through `run_graph_build_partitioned`.
- `src/main.rs:8789-8790` routes `impg graph --gfa-engine seqwish` through
  `graph::run_graph_build`.
- `src/main.rs:8801-8817` routes `impg graph --gfa-engine pggb` through
  `graph::run_graph_build_pggb`.
- `src/main.rs:12102-12132` routes `query -o gfa` through
  `impg::dispatch_gfa_engine`; `src/main.rs:12147-12176` does the same for
  VCF output, which first builds GFA.

Library dispatch used by query/partition:

- `src/lib.rs:36-58` defines `GfaEngine::Pggb`, `Seqwish`, and `SyngNative`;
  `SyngNative` is documented as emitting PAF then feeding seqwish.
- `src/lib.rs:92-120` stores the shared `GraphBuildConfig` in `EngineOpts`.
- `src/lib.rs:828-845` defines `dispatch_gfa_engine`, shared by `query -o gfa`
  and `partition -o gfa`.
- `src/lib.rs:984-1015` dispatches `GfaEngine::Seqwish` and `GfaEngine::Pggb`
  through `graph::generate_gfa_seqwish_from_intervals`.
- `src/lib.rs:1066-1160` dispatches `GfaEngine::SyngNative`; anchor-seeded or
  pairwise BiWFA PAF is passed to `syng_graph` and then to the shared seqwish
  tail.
- `src/lib.rs:1204-1305` runs per-partition `dispatch_gfa_engine_inner`, so
  partitioned seqwish/pggb/syng-native paths still share the same tail.

Interval and syng wrappers:

- `src/graph.rs:1030-1078` implements
  `generate_gfa_seqwish_from_intervals`, which extracts intervals to a temp
  FASTA and calls `crate::commands::graph::build_graph`.
- `src/syng_graph.rs:1-15` documents syng-native as pairwise BiWFA -> PAF ->
  existing seqwish induction.
- `src/syng_graph.rs:924-950` builds PAF from named sequences and calls
  `build_gfa_from_paf_and_sequences`.
- `src/syng_graph.rs:952-1088` writes FASTA/PAF, filters generated PAF, adapts
  `min_match_len` when requested, wraps the data as `AlignmentResult`, and
  calls `crate::commands::graph::induce_graph_from_alignment`.

Crush/local replacement callers:

- `src/resolution.rs:9097-9148` builds the `GraphBuildConfig` used for local
  replacement seqwish induction.
- `src/resolution.rs:8679` routes `TopFlubbleSweepga` to the SweepGA/seqwish
  replacement path.
- `src/resolution.rs:8685-8687` routes `Allwave`, `Sweepga`, and `Wfmash` to
  seqwish-backed replacement builders.
- `src/resolution.rs:9431-9506` builds AllWave PAF and calls
  `syng_graph::build_gfa_from_paf_and_sequences`.
- `src/resolution.rs:9571-9714` builds SweepGA/wfmash PAF, configures the
  seqwish tail, and calls `syng_graph::build_gfa_from_paf_and_sequences`.
- `src/resolution.rs:9897-9907` implements `method=wfmash` as the SweepGA
  seqwish replacement path with the aligner pinned to wfmash.

Render/local graph helper:

- `src/commands/render.rs:45-49` accepts render engines `seqwish` and `pggb`
  as local graph engines.
- `src/commands/render.rs:321-345` maps render `seqwish` to
  `graph::run_graph_build` and render `pggb` to `graph::run_graph_build_pggb`.

Parser/tests that exercise or preserve the path:

- `src/graph_pipeline.rs:148-157` preserves legacy pipeline specs such as
  `wfmash:seqwish:10k`.
- `tests/test_graph_seqwish.rs` exercises the CLI seqwish graph path.
- `tests/test_crush_integration.rs` exercises seqwish-backed crush replacement
  behavior.
- `tests/test_syng_integration.rs` exercises `query -o gfa --gfa-engine
  seqwish`/`pggb` parsing and syng-query graph output paths.
