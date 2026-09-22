# SYNG Integration Architecture

Status: implementation-grounded design note. This is not a redesign proposal.
It describes the current IMPG/SYNG integration as implemented in the Rust code.

Primary implementation references:

- `src/syng.rs`: `SyngIndex`, syncmer extraction, GBWT/KmerHash ownership,
  sidecar formats, path/position lookup, query, mapping, and regional GBWT build.
- `src/syng_ffi.rs`: FFI declarations for the vendored SYNG C library.
- `src/syng_parallel.rs`: packed-syncmer dictionary extraction and deterministic
  sorted/deduplicated dictionary construction for `--parallel-dictionary`.
- `src/main.rs`: CLI parsing and dispatch for `syng`, `syng-repair`, `syng2gfa`,
  `query`, `partition`, `map`, `read-index`, `genotype`, and `infer`.
- `src/syng_transitive.rs`: chain filtering, SweepGA scaffold chaining, ends-only
  boundary refinement, transitive query hops, and chain thresholds.
- `src/lib.rs`: `SyngImpgWrapper`, GFA engine dispatch, `syng` and `syng-local`
  graph output, crush/smooth/sort graph transforms.
- `src/commands/syng2gfa.rs`: whole-index and regional SYNG-to-GFA conversion,
  raw/blunt modes, frequency-aware private splitting, gap fill, and N cutting.
- `src/pack.rs`, `src/projection.rs`, `src/commands/genotype.rs`,
  `src/commands/infer.rs`: read projection, pack/proj evidence, cosine genotype,
  and infer/stitching interaction.

Existing user-facing and validation notes that overlap this design:
`README.md`, `docs/syng-gfa-query.md`, `docs/syng-position-query-index.md`,
`docs/syng-parallel-construction.md`, `docs/genotype-evidence-audit.md`,
`docs/graph-quality-validation.md`, `docs/local-syng-parameter-sweep.md`, and
`docs/evaluations/local-compression-testbed-fast-synthesis.md`.

## High-Level Approach

IMPG embeds SYNG as an alignment-free backend. The SYNG contribution is a
syncmer dictionary plus a GBWT-like path index over syncmer nodes. IMPG uses
that implicit graph to discover homologous intervals, project reads onto
syncmer nodes, materialize regional graphs, and build feature vectors for
genotype/infer. It does not require a precomputed all-vs-all PAF or an explicit
variation graph for the initial query step.

The representation is intentionally succinct:

- A SYNG node is a canonical syncmer stored in SYNG's `KmerHash`.
- A SYNG path is a sequence of signed syncmer node IDs plus inter-syncmer
  offsets stored in SYNG's `SyngBWT`.
- For each input sequence, IMPG inserts the forward path and a reverse-complement
  path into the GBWT. The IMPG name map records the forward path metadata used
  by queries and graph output.
- Non-syncmer sequence between adjacent syncmers is not stored as full DNA in
  the SYNG index. It is recovered from FASTA/AGC sequence files when graph or
  FASTA output needs source spelling. If sequence files are not supplied to
  GFA conversion, gap DNA is emitted as `N` with warnings.
- Coordinates are made practical by IMPG sidecars (`.syng.pstep` and
  `.syng.spos`) layered on top of SYNG's files. These are IMPG indexes, not
  native SYNG files.

In IMPG, SYNG participates in several larger workflows:

- `impg query -a <syng-prefix>` discovers homologous target intervals from a
  path/range, optionally refines boundaries with BiWFA, and emits BED/BEDPE,
  FASTA, regional GBWT, GFA, or VCF.
- `impg partition -a <syng-prefix>` wraps `SyngIndex` in `SyngImpgWrapper` so
  the existing partition code can use SYNG queries instead of alignment-backed
  interval trees.
- `impg map` projects reads or probes onto syncmer nodes and can emit GAF,
  PAF-like coordinate projections, pack vectors, or projection bundles.
- `query -o gfa:syng...` converts query-selected SYNG walks into local GFA.
  `gfa:syng-local...` extracts query-selected sequence, builds a fresh local
  SYNG index, and then materializes GFA.
- `genotype cos` and `infer` use SYNG query/map evidence as graph-feature
  vectors. The implemented scorer is cosine over syncmer-node or graph-node
  features; it is not a probabilistic imputation engine.

## Build Integration

The build command is `impg syng` (`Args::Syng` in `src/main.rs`). It accepts
either `--agc <file>` or `-f/--fasta <file>` and writes a prefix:

```bash
impg syng -f panel.fa -o panel.syng \
  --syncmer-length 63 --smer-length 8 --syncmer-seed 7 \
  --position-sample-rate 256 -t 16
```

For AGC input, IMPG uses `ragc_core`. Path names are created by
`syng_sequence_name`: the primary contig token is used directly when it already
contains `#` PanSN separators; otherwise the path name is `contig@sample`.
For FASTA input, the sequence name is the first whitespace-delimited token after
`>`.

Syncmer CLI parameters follow the paper-facing names:

- `--smer-length` is the inner `s` length. Default: `8`.
- `--syncmer-length` is the total syncmer length `k`. Default: `63`.
- The total syncmer length must be odd and greater than `s`.
- `--syncmer-seed` defaults to `7`.
- Legacy `--syncmer-k` maps to `--smer-length`, and legacy `--syncmer-w` is the
  internal window component. Internally `SyncmerParams { k, w, seed }` means
  `k = s` and total length `k + w`.

Build behavior:

1. Resolve syncmer parameters with `resolve_syng_syncmer_params`.
2. Create `SyngIndex::new(params)`, which owns a SYNG `KmerHash`, `SyngBWT`,
   and `Seqhash`.
3. Enable online sampled path-step collection with
   `enable_online_sampled_positions(position_sample_rate)`. The default sample
   rate is `256`, and `0` is rejected.
4. For each sequence, convert bases to SYNG numeric encoding. Non-ACGT bases
   are mapped to `A`, matching the current SYNG convention in `src/syng.rs`.
5. Extract syncmers with the SYNG seqhash iterator. Add missing canonical
   syncmers to the `KmerHash` unless the two-pass dictionary path is active.
6. Insert a forward GBWT path. Store the first signed syncmer node, its
   occurrence count/rank, number of syncmers, and first syncmer bp position in
   `.names` metadata.
7. Insert the reverse-complement GBWT path by reversing the syncmer order and
   negating node IDs.
8. Record per-path sampled path steps on the regular syncmer-step grid:
   step `0`, every `N`th step, and the final step on each path.
9. Finalize sampled positions. Because GBWT occurrence ranks are not stable
   until all paths are inserted, the final `.pstep` checkpoint index is rebuilt
   from the completed GBWT rather than trusting ranks captured during insertion.
10. Save core SYNG files, position sidecars, and metadata.

The `--parallel-dictionary` option changes only dictionary construction. It
performs a deterministic prepass that extracts packed canonical syncmers from
all sequences, sorts and deduplicates them, preloads the `KmerHash` in that
global order, and then replays paths with `RequireExisting` lookup. This avoids
thread-order-dependent node ID assignment for large FASTA/AGC inputs. Path
insertion itself remains an ordered replay through the SYNG GBWT.

Sequences shorter than one syncmer, or sequences that yield no syncmers, are
kept in the name map with length but are not indexed in the GBWT. They therefore
cannot contribute syncmer anchors.

### Build Outputs

For a prefix `panel.syng`, the current file set is:

| File | Owner | Contents |
| --- | --- | --- |
| `panel.syng.1gbwt` | SYNG | GBWT/path index over signed syncmer nodes and offsets. |
| `panel.syng.1khash` | SYNG | Canonical syncmer dictionary; node IDs are `1..=N`. |
| `panel.syng.names` | IMPG | Path number, name, sequence length, GBWT forward-path start info, number of syncmers, first syncmer bp position. |
| `panel.syng.pstep` | IMPG | Per-path sampled path-step checkpoints for forward path-coordinate walking. |
| `panel.syng.spos` | IMPG | Occurrence-major sampled syncmer-position checkpoints for reverse lookup from syncmer occurrence/rank to path position. |
| `panel.syng.meta` | IMPG | TSV metadata: format version, internal `syncmer_k`, `syncmer_w`, seed, total syncmer length, and `smer_length`. |

If the prefix does not end in `.syng`, sidecars are named
`<prefix>.syng.names`, `<prefix>.syng.pstep`, `<prefix>.syng.spos`, and
`<prefix>.syng.meta`. If the prefix already ends in `.syng`, IMPG avoids
duplicating `.syng` and uses `<prefix>.names`, `<prefix>.pstep`, etc. The core
SYNG files are always `<prefix>.1gbwt` and `<prefix>.1khash`.

Load-time prefix detection accepts a bare prefix or paths ending in `.1khash`,
`.1gbwt`, `.syng.names`, `.names`, `.syng.spos`, `.spos`, `.syng.pstep`,
`.pstep`, `.syng.meta`, or `.meta` where the sibling `.1khash` is available.

Metadata is authoritative for syncmer parameters on load. Callers often pass
`SyncmerParams::default()` to `SyngIndex::load` or `SyngMatcher::load`; the
actual parameters are read from `.meta`.

`impg syng-repair -a <prefix> --position-sample-rate <N> --force` can rebuild
or resample `.pstep` and `.spos` from existing `.1gbwt`, `.1khash`, `.names`,
and `.meta` files. Normal `SyngIndex::load` requires the sidecars and rejects
missing/stale `.spos` or `.pstep`; `load_for_repair` is the tolerant path used
by repair.

## Index And Data Structures

### `.1gbwt`

`SyngBWT` stores signed syncmer paths and inter-syncmer offsets. IMPG uses it
for:

- path traversal from a recorded start node/count,
- exact walk matching (`syngBWTmatchStart`, `syngBWTmatchNext`),
- rank advancement from a syncmer occurrence to a sampled checkpoint
  (`syngBWTadvanceRank`),
- reconstructing path steps during sidecar repair and regional graph output.

The GBWT has both forward and reverse-complement path spellings, but IMPG's
name map points at the forward path start for each input sequence.

### `.1khash`

`KmerHash` stores canonical syncmer sequences and maps extracted syncmers to
node IDs. It is required for build, query, GFA segment sequence extraction, and
read matching. `SyngMatcher` intentionally loads only `.meta` and `.1khash` for
fast read-to-syncmer projection modes that do not need coordinates.

### `.meta`

`.meta` is a small TSV written last during `save`. It records the IMPG metadata
format version and the syncmer scheme. Query, map, pack, and read-index commands
load parameters from this file so reads are scanned with the same scheme as the
index.

### `.names`

`.names` maps internal path numbers to user path names and lengths. New files
carry seven columns:

```text
path_number  name  length  start_node  start_count  num_syncmers  first_syncmer_pos
```

The loader also accepts older three- and six-column formats. In six-column
files, `first_syncmer_pos` defaults to `0`, which can preserve old behavior but
does not carry the exact first-syncmer coordinate for freshly rebuilt indexes.

### `.pstep`

`.pstep` is the forward positional index: per-path sampled GBWT checkpoints.
Each checkpoint records the path index, syncmer step index, bp position, signed
node, previous node, previous offset, and traversal rank needed to resume a GBWT
walk. Data is varint/delta encoded per path and memory-mapped on load.

This index answers "what syncmer nodes overlap path interval `[start, end)`?"
without scanning the whole path:

1. Find the checkpoint at or before `start - (syncmer_length - 1)`.
2. Recreate the GBWT path state from the checkpoint.
3. Walk forward with `syngBWTpathNext`.
4. Emit every syncmer whose span overlaps `[start, end)`.
5. Stop after the interval end or the path's final syncmer.

`walk_path_range` and graph/candidate feature extraction depend on `.pstep`.

### `.spos`

`.spos` is the reverse positional index: occurrence-major sampled checkpoints
keyed by `(signed_node, absolute incoming rank)`. It is built from `.pstep`
after the GBWT is complete. It stores checkpoint path index and bp position, and
is memory-mapped/binary-searched on load.

This index answers "where does this signed syncmer occurrence/rank land in
path coordinates?"

1. Start from a signed node and a GBWT occurrence rank.
2. Look for a matching checkpoint in `.spos`.
3. If no checkpoint is present at the current occurrence, advance one GBWT step
   with `syngBWTadvanceRank`, accumulating walked bp.
4. Repeat up to the sample rate.
5. When a checkpoint is found, subtract walked bp from checkpoint bp position to
   recover the original occurrence's exact path coordinate.

This is the hot path for locating query syncmer anchors across all target paths.

### Replacement And Legacy Notes

Older code still contains a `SampledPositions` representation, but normal load
and query now use `SampledPathSteps` (`.pstep`) plus `SampledCheckpointIndex`
(`.spos`). If a persistent `.spos` is missing in an in-memory/test context, IMPG
can build a fallback hash table from `.pstep`; normal on-disk index loads
require the persistent `.spos`.

## Query Integration

`impg query` routes to SYNG when `-a` names a SYNG prefix or sidecar. The SYNG
query path supports `bed`, `bedpe`, `gfa`, `vcf`, `fasta`, and `gbwt` output.
`auto` resolves to `bed`. PAF is not a SYNG query output; coordinate projection
for reads is under `impg map -o paf`.

Important CLI defaults in `src/main.rs`:

- `--syng-padding 120`
- `--syng-extension 0`
- `--syng-extend-budget 1000`
- `--syng-min-chain-anchors 20` as a cap on the adaptive threshold; `0`
  disables anchor-count filtering
- `--syng-min-chain-fraction 0.0`
- `--syng-seed-drop-top-fraction 0.0005`
- `--syng-seed-max-occurrences 0` (disabled)
- `--syng-seed-walk-anchors 5`
- boundary realignment is enabled by default and requires `--sequence-files` or
  `--sequence-list`; `--syng-raw` disables it and emits raw syncmer-resolution
  intervals.

### Raw Region Lookup

The core path is `SyngIndex::query_region_with_anchors_ext_visited_seed_filtered`:

1. Resolve the source path name through `.names`.
2. Expand the source interval by `--syng-extension` for syncmer discovery.
3. Use `.pstep` to walk only syncmers overlapping that expanded path interval.
4. If seed filtering is disabled, group query syncmers by unsigned node ID and
   query orientation.
5. Locate every occurrence of each seed node in both GBWT orientations via
   `.spos`.
6. For each located occurrence and each query occurrence, assign target strand:
   `+` when query and target orientations agree, `-` otherwise.
7. Create a small padded interval around the target syncmer occurrence. Use a
   diagonal-like signature (`target_pos - query_pos` on `+`, `target_pos +
   query_pos` on `-`) to deduplicate repeated anchors.
8. Group by target path and strand. Merge only overlapping/adjacent intervals
   within each group. Different strands are kept separate.
9. Sort by genome, strand, start for deterministic output.

The interval merge inside `SyngIndex` is strict overlap/adjacency
(`iv.start <= last.end`). User-facing `-d/--merge-distance` is applied later by
the shared IMPG output pipeline and by the transitive chaining limits.

### Frequency Masking And Bounded Walk Seeds

The default query seed filter is enabled because `drop_top_fraction` defaults
to `0.0005`. In that mode IMPG does not locate every single shared syncmer
individually. Instead:

1. Count occurrence totals for unique query-local seed nodes.
2. Drop the highest-frequency fraction and any node above
   `--syng-seed-max-occurrences` if the cap is set.
3. Split the ordered query walk into runs that exclude dropped nodes.
4. From each retained run, create bounded exact GBWT walk seeds of
   `--syng-seed-walk-anchors` consecutive syncmers, stepping by that same
   number.
5. Add equivalent reverse-orientation seed tasks.
6. Match each seed exactly in the GBWT with `syngBWTmatchStart` /
   `syngBWTmatchNext`.
7. Locate the final seed node occurrence range through `.spos`, reconstruct the
   seed target start from relative seed offsets, emit anchors, then merge.

If too few seed steps remain for bounded walk seeding, the code falls back to
the node-position path for the retained nodes. Dropped high-frequency syncmers
are skipped only for seeding candidate ranges. Once a target range survives,
later path walking, graph materialization, pack/genotype features, and GFA
conversion can still see all syncmers in that path interval.

### Chaining, Boundary Refinement, And Transitive Query

Default non-raw `impg query` calls
`syng_transitive::query_transitive_ext_with_seed_filter`. One hop:

1. Runs the SYNG anchor query above.
2. Chains anchors with `chain_anchors_with_sweepga_scaffold_mass`, which converts
   anchors to PAF-like records and reuses SweepGA's `PafFilter` in scaffold mode.
   There is no separate local plane-sweep chainer in the SYNG module.
3. Computes an effective minimum anchor count from query length, syncmer
   density, an assumed 95% identity, and a 0.10 expected-syncmer fraction, then
   caps it at `--syng-min-chain-anchors`.
4. Filters chains below the effective anchor count and chains whose query
   extent is below `--syng-min-chain-fraction * query_range_len`.
5. Removes overlapping opposite-strand duplicate chains on the same target path,
   keeping the majority-anchor strand.
6. Refines only chain boundaries. Multi-anchor chains use two ends-free BiWFA
   alignments, one from the first anchor backward and one from the last anchor
   forward, with windows bounded by `--syng-extend-budget`, neighboring chains,
   and source/target span caps. Singleton chains use linear projection.

The `-d/--merge-distance` value participates in chain gap/drift/span limits and
the shared output merge. Large `-d` values can absorb larger SVs into a chain;
small values keep intervals more fragmented.

When `--transitive` is set, `impg query` raises SYNG max depth from one hop to
the requested transitive depth. The transitive implementation keeps a visited
seed-node set so repeated hops do not keep re-expanding through the same
syncmer nodes. The `SyngImpgWrapper` used by partition exposes
`query_transitive_dfs` and `query_transitive_bfs` through the common
`ImpgIndex` trait, but for partition those trait methods delegate to SYNG
query behavior rather than to an alignment interval tree.

### Query Outputs

For `bed` and `bedpe`, SYNG intervals are converted to IMPG
`AdjustedInterval` records and flow through the common BED/BEDPE output and
merge code.

For `fasta`, IMPG fetches each SYNG interval from the supplied sequence index.

For `gbwt`, IMPG fetches each interval's sequence and calls
`SyngIndex::build_region_gbwt` with the source index's syncmer parameters,
writing a regional `.1gbwt`, `.1khash`, `.names`, `.pstep`, `.spos`, and
`.meta` set under `-O`.

For `gfa` and `vcf`, SYNG first discovers intervals, then dispatches graph
generation. `vcf` converts the generated GFA through the native POVU path.

## Graph Output Integration

SYNG graph output is reachable through:

- `impg query -a <syng-prefix> ... -o gfa:syng...`
- `impg query -a <syng-prefix> ... -o vcf:syng...`
- `impg syng2gfa -a <syng-prefix> ...`

The compact output form `-o gfa:<engine-spec>` is parsed by
`apply_gfa_output_engine_shorthand` and the same engine parser as
`--gfa-engine`. For SYNG indexes, `gfa:syng` defaults to `syng:blunt` and the
default final sort pipeline is `Ygs`. `gfa:syng:raw` preserves the native
overlap graph and cannot be combined with `:crush`. `gfa:syng:crush` means:

1. query-selected SYNG intervals,
2. direct SYNG regional GFA extraction,
3. blunt exact path materialization,
4. optional frequency/private-split policy,
5. crush graph-resolution transform,
6. optional smooth stage if requested,
7. self-loop normalization and `gfasort` pipeline `Ygs` unless `:nosort` is
   supplied.

`gfa:syng-local` is a different engine. It extracts the query-selected
sequences, builds a fresh temporary regional SYNG index with either requested
local `k/s/seed` or the source index parameters, writes a temporary FASTA, loads
the temporary regional index, and then calls the same `syng2gfa` range writer.
For `syng-local`, `k/s/seed` are rebuild parameters. For `syng`, they are
assertions against the loaded global index.

`impg syng2gfa` dumps a whole SYNG index instead of query-selected intervals.
It emits GFA 1.0 with `P` lines by default or GFA 1.1 with `W` lines when
requested. It supports `--gfa-mode raw` and `--gfa-mode blunt`; whole-index
conversion does not apply the query-local frequency mask unless called through
the range writer with an engine mask.

### Raw And Blunt Materialization

The raw writer emits:

- `H` header,
- one `S` per used syncmer node,
- additional `S` lines for inter-syncmer, prefix, and suffix gap segments,
- `L` edges with nonzero overlaps where syncmers overlap in the source path,
- `P` or `W` paths depending on GFA version.

Gap segments are interned by exact sequence plus local signed-syncmer context,
so identical gap DNA can be shared only when its local graph context also
matches. If source sequence files are missing, gap segments are filled with
`N`s.

The blunt writer currently uses an exact zero-overlap materializer
(`write_exact_blunt_path_work_gfa`) rather than requiring consumers to run an
overlap-aware parser. It emits full untrimmed syncmer anchor segments, trimmed
syncmer-slice segments when needed, gap segments, zero-overlap links, and exact
source-spelling paths. Blunt range output currently requires GFA 1.0.

### Frequency Mask, Private Splitting, And N Cutting

For query-local `gfa:syng` and `gfa:syng-local`, the default frequency policy is
`SyngGfaFrequencyMask::local_default()`:

- drop/private-split the top `0.0005` local syncmer fraction,
- use occurrence-aware high-frequency masking,
- rescue high-frequency occurrences that lie in supported runs of at least
  `10` syncmers or exact spans of at least `1000` bp,
- require local shared context run length `5` by default,
- split rare repeated-copy local contexts with max minor count `2` and dominant
  fraction `0.80`,
- split dispersed scaffold-glue candidates selected by path-copy/scaffold
  context rules.

Use `:nomask` / `:nofilter` to disable this policy, or `:mask,...` to tune
parameters such as `top`, `max-occ`, `freq-run`, `freq-span`, `min-run`,
`sequence-k`, `repeat-minor`, and `repeat-dominance`.

`:cut-ns` / `:cut-n` in a SYNG graph engine does not clip input intervals. It
cuts N-runs out of fetched gap DNA during SYNG GFA materialization and splits
emitted paths at those breaks. `cut-n-min-run` controls the minimum N-run
length, defaulting to `1` when N cutting is enabled.

There is also a separate graph pipeline stage `cut-n=<bp>` before non-SYNG or
sequence-extracted engines, for example `-o gfa:cut-n=100:pggb`. That clips
terminal N-runs from extracted query sequences before graph construction and is
not the same as SYNG gap-DNA `:cut-ns`.

## Mapping Integration

`impg map` is the read/probe projection command for SYNG indexes:

```bash
impg map -a panel.syng -q reads.fq -o gaf
impg map -a panel.syng -q reads.fq -o paf --min-anchors 3 --max-hits 10
impg map -a panel.syng -q reads.fq -o pack -O sample.pack
impg map -a panel.syng -q reads.fq -o pack-tsv -O sample.pack.tsv.zst
impg map -a panel.syng -q reads.fq -o proj -O sample.proj
```

The default output is `gaf`. Supported formats are `gaf`, `paf`, binary
`pack` aliases (`pack`, `packbin`, `pack-bin`, `bpack`), text pack aliases
(`pack-tsv`, `pack-text`, `packtsv`), and `proj`.

Read syncmers are always extracted using parameters loaded from the index
metadata. `SyngMatcher` scans both the input read and its reverse complement and
keeps the orientation with more dictionary hits. Reverse-complement hits are
reported with negated signed nodes and query positions transformed back into
the original read coordinate system.

### GAF

GAF output uses the lightweight `SyngMatcher` and loads only `.meta` and
`.1khash`. It streams FASTA/FASTQ in chunks, uses worker-local SYNG seqhash
handles across Rayon/thread workers, and writes output chunks in input order.

Each retained read emits one GAF record with:

- query name/length and covered query span,
- strand `+`,
- path string as `>node` / `<node` over syncmer node IDs,
- synthetic path length equal to `syncmer_count * syncmer_length`,
- tags `an:i:<anchor-count>`, `sk:i:<syncmer-length>`, and `qp:B:I,...` query
  positions.

`--min-anchors` for GAF is a raw matched-syncmer occurrence count.

### PAF

PAF output loads the full `SyngIndex`, including `.pstep` and `.spos`, reads the
query FASTA/FASTQ into memory, projects query syncmers to indexed genome
coordinates with `SyngIndex::map_sequence`, chains anchors with the same
SweepGA scaffold-mass chainer, filters by `--min-anchors`, sorts by anchor
count, and truncates by `--max-hits` when nonzero. It emits PAF-like rows with
`an:i` and `sk:i` tags. It is a syncmer-anchor projection, not a base-level
alignment with CIGAR.

### Pack And Projection Bundles

Pack output uses the lightweight `.1khash` matcher and streams reads. For each
read:

1. Match syncmers in the best read orientation.
2. Collect unsigned node IDs.
3. Sort/deduplicate node IDs within the read.
4. Apply `--min-anchors` to the number of distinct node IDs.
5. Increment each distinct retained node by one.

Consequences: pack counts are per-read distinct-node support. Repeated
occurrences of the same node inside one read do not increase that node's count,
orientation is ignored, and counts are not weighted by read length, node length,
base coverage, mapping quality, or syncmer occurrence count.

Text pack is a TSV `#node_id count`. Binary pack uses `src/pack.rs`: a dense
u8 vector over node IDs, independently zstd-compressed blocks, and an overflow
table for counts above `255`. Binary `--pack-compression-level` defaults to
`12`, must be `1..=22`, and `--pack-block-size` defaults to
`DEFAULT_BINARY_BLOCK_SIZE`.

`-o proj -O sample.proj` writes:

```text
sample.proj/
  manifest.json
  sample.pack
  reads.gaf.zst
```

The manifest declares feature space `syng-syncmer-node`, records the SYNG
prefix, and points to the pack and optional GAF read-walk file. `proj` requires
`-O`.

`impg read-index` is adjacent but separate: it builds a sampled node-major
read-to-syncmer inverted index (`.r2s.meta`, `.r2s.sample`, `.r2s.post`) from
reads and a SYNG dictionary. It is not required by `query`, `map`, or
`genotype cos`.

## Genotype And Infer Interaction

Implemented SYNG genotype is the `genotype cos` / `cosigt` path in
`src/commands/genotype.rs`.

For a SYNG backend:

1. Load the SYNG index and sample pack.
2. Parse `--target-range`.
3. Discover haplotype candidates with `query_region_with_anchors_ext`.
4. In `spanning` mode, group all hits by `(path, strand)`, merge to the
   candidate hull, deduplicate anchors, then require `--min-anchors` and
   `--min-span-fraction`.
5. In `overlapping` mode, score each gathered hit independently and require
   nonempty interval plus `--min-anchors`; span fraction is reported but not
   used as a filter in that branch.
6. For each candidate interval, walk the indexed path range with `.pstep` and
   count unsigned syncmer node IDs. Repeated node occurrences in a candidate
   interval increase the candidate feature count.
7. Compute one-haplotype cosine scores against the sample pack, optionally
   truncate candidates by `--candidate-top-k`, enumerate ploidy-sized
   combinations with replacement up to `--max-combinations`, and report top
   cosine genotypes.

The current cosine score uses unoriented syncmer-node feature vectors. It does
not use mapping quality, read likelihoods, node length normalization in the
SYNG-node backend, orientation-aware cosine components, or a haplotype-copy
HMM.

`genotype cos` also has a graph backend that can consume a GFA, render bundle,
or dynamically built query graph plus graph-projected pack evidence. That is a
different feature space (`gfa-segment` or `variation-graph-node`) and is
documented in the genotype architecture notes.

`impg infer` reuses the SYNG local cosine scorer. It can type one range,
target BED rows, partition BED rows, or internally discovered SYNG partitions.
It then optionally stitches local calls into a mosaic with a beam search.

Implemented infer relationships:

- Local call rows use the same pack coverage, candidate discovery, candidate
  top-k, cosine combination search, and QV behavior as `genotype cos`.
- If `--proj` or `--gaf` provides GAF read walks, infer can build read-walk
  evidence for stitching and optional mosaic/FASTA/GFA emission.
- GAF read-walk evidence does not change the local cosine rows unless stitching
  or mosaic emission is active.
- With `qp:B:I` query positions, infer can compute SYNG GBWT MEMs over signed
  read walks and candidate walks; without positions, it falls back to
  whole-walk orientation/LIS-style support.

Roadmap/speculative behavior to keep separate from current implementation:

- full probabilistic imputation from read likelihoods,
- base-level read realignment into candidate haplotypes,
- orientation-aware or length-normalized SYNG-node cosine,
- cohort-scale imputation using read-index postings,
- production-quality truth-calibrated genotype confidence.

## Algorithm Sketches

### Build

```text
build_syng(input, prefix, params, sample_rate):
  index = SyngIndex::new(params)
  index.enable_online_sampled_positions(sample_rate)

  if --parallel-dictionary:
    packed = parallel extract canonical syncmers from all sequences
    dictionary = sort_dedup(packed)
    index = SyngIndex::new_with_packed_syncmer_dictionary(params, dictionary)
    index.enable_online_sampled_positions(sample_rate)
    replay sequences with RequireExisting lookup
  else:
    stream sequences and add missing syncmers during replay

  for each sequence:
    numeric = encode A/C/G/T as 0/1/2/3; other bases as A
    syncmers = seqhash iterator over numeric sequence
    if no syncmers:
      record name and length only
      continue

    for each syncmer:
      node_id = khash add/find canonical syncmer

    start forward GBWT path at first signed node
    for each next syncmer:
      add signed node with offset from previous syncmer
      collect path-step record
    finish forward path

    start reverse-complement GBWT path at -last node
    add remaining nodes in reverse order with negated signs
    finish reverse path

    record name, length, start_node, start_count, num_syncmers,
      first_syncmer_pos
    record sampled path steps: step 0, every Nth step, final step

  rebuild sampled path steps from completed GBWT
  derive occurrence-major checkpoints from path steps
  save .1gbwt, .1khash, .names, .pstep, .spos, .meta
```

### Path-Position To Syncmer Traversal

```text
walk_path_range(path_idx, start, end):
  scan_start = max(0, start - (syncmer_len - 1))
  checkpoint = .pstep.checkpoint_at_or_before(path_idx, scan_start)
  if no checkpoint:
    checkpoint = path start from .names

  state = gbwt_state_from(checkpoint)
  current = checkpoint
  while current exists:
    if current.bp_pos < end and current.bp_pos + syncmer_len > start:
      emit (current.signed_node, current.bp_pos)
    if current.bp_pos >= end or current is final path step:
      break
    current = advance one GBWT path step from state
```

### Syncmer Occurrence To Path Position

```text
locate_signed_node_occurrence_range(signed_node, low, high):
  for rank in low..high:
    node = signed_node
    abs_rank = rank
    walked_bp = 0

    repeat at most sample_rate + 1 times:
      if .spos has checkpoint(node, abs_rank):
        emit path_idx and checkpoint.bp_pos - walked_bp
        break

      next = syngBWTadvanceRank(node, abs_rank)
      if no next:
        break
      walked_bp += next.offset
      node = next.node
      abs_rank = next.rank
```

### Query Interval Discovery And Merge

```text
query_region(path, start, end, padding, extension, seed_filter):
  query_walk = walk_path_range(path, start - extension, end + extension)

  if seed_filter.enabled:
    counts = occurrence counts for unique query nodes
    dropped = top_fraction(counts) union nodes with count > max_occ
    runs = split query_walk at dropped nodes
    seed_tasks = consecutive walk_anchors syncmers from each run and RC run

    for task in parallel:
      (last_node, low, high) = exact GBWT match for task
      hits = locate_signed_node_occurrence_range(last_node, low, high)
      for hit:
        reconstruct target_start from seed relative offsets
        emit padded interval plus seed anchors
  else:
    query_node_positions = node -> [(query_pos, query_orientation)]
    for node in parallel:
      hits = locate occurrences of node in both target orientations
      for hit and query position:
        strand = '+' if orientations agree else '-'
        signature = target_pos - query_pos on '+'
                  = target_pos + query_pos on '-'
        deduplicate anchors by signature
        emit padded interval

  group intervals by (target_path, strand)
  sort each group by start
  merge only when next.start <= current.end
  sort/deduplicate anchors inside merged intervals
  return deterministic genome/strand/start order
```

### Chain And Boundary Refinement

```text
refined_one_hop(query_range):
  raw_hits = query_region_with_anchors(...)
  chain_gap, drift, span_cap = query_scaled_chain_limits(query_len,
    extend_budget, syncmer_len, merge_distance)
  effective_min = min(user_cap,
    ceil(query_len * syncmer_density * 0.95^syncmer_len * 0.10))

  chains = SweepGA scaffold-mass chain(raw_hits,
    scaffold_gap=chain_gap,
    min_scaffold_length=effective_min * syncmer_len)

  keep chains with anchor_count >= effective_min
  keep chains with query_anchor_extent >= min_chain_fraction * query_len
  dedupe overlapping opposite-strand chains on same target

  for each chain in parallel:
    if chain has <2 anchors:
      linearly project query endpoints from anchor
    else:
      fetch source and target edge windows
      run ends-free BiWFA on left and right ends only
      fall back to linear projection on failure
    emit refined HomologousInterval
```

### Graph Extraction And Conversion

```text
syng_region_gfa(intervals, mode, sequence_index, frequency_mask):
  ranges = map IMPG query intervals to SYNG path_idx/start/end/strand

  if mode == syng-local:
    sequences = fetch interval DNA
    temp_index = build regional SYNG index from sequences
    ranges = full paths in temp_index
    index = temp_index

  if frequency_mask enabled:
    local_walks = collect syncmer walks for ranges
    selected_nodes = top-frequency and context/glue candidates
    split_occurrences = occurrence-level private splits,
      with run/span/scaffold rescue

  for range in parallel:
    walk syncmers with .pstep
    fetch gap DNA from sequence index or fill with N
    if cut-ns enabled:
      split gap DNA at N-runs and insert path breaks
    collect raw path steps and links

  if raw:
    emit syncmer S lines, gap S lines, overlap L lines, P/W paths
  if blunt:
    emit exact zero-overlap syncmer/gap/clone segments,
      zero-overlap links, and exact source-spelling P paths

  apply optional crush transform
  apply optional smooth transform
  apply optional self-loop normalization + gfasort pipeline
```

### Syncmer Read Mapping

```text
map_read(read):
  params = load from .meta
  fwd = syncmers in read that exist in .1khash
  rev = syncmers in reverse_complement(read) that exist in .1khash
  matches = orientation with more hits

  if output == gaf:
    require raw match count >= min_anchors
    sort matches by query position
    emit >node/<node walk and qp:B:I positions

  if output == pack:
    distinct_nodes = sorted unique unsigned node IDs
    require distinct_nodes.len >= min_anchors
    for node in distinct_nodes:
      counts[node] += 1

  if output == paf:
    group matches by unsigned node/query positions
    locate target occurrences through .spos
    chain anchors with SweepGA scaffold-mass chainer
    emit coordinate hits sorted by anchor count
```

## CLI Examples

Environment-specific paths below are examples. The HPRCv2 paths match existing
repo notes and scripts but are developer-machine local, not CI fixtures.

C4 from ODGI's test graph:

```bash
curl -sL -O https://raw.githubusercontent.com/pangenome/odgi/master/test/chr6.C4.gfa
odgi paths -i chr6.C4.gfa -f > chr6.C4.fa
samtools faidx chr6.C4.fa

impg syng -f chr6.C4.fa -o c4.syng \
  --syncmer-length 63 --smer-length 8 --syncmer-seed 7 \
  --position-sample-rate 256 -t 4

impg query -a c4.syng \
  -r 'grch38#chr6:31972046-32055647:0-10000' \
  --sequence-files chr6.C4.fa -d 150 -o bed

samtools faidx chr6.C4.fa 'grch38#chr6:31972046-32055647:5000-7000' > probe.fa
impg map -a c4.syng -q probe.fa -o gaf
impg map -a c4.syng -q probe.fa -o paf
impg map -a c4.syng -q probe.fa -o pack-tsv
impg map -a c4.syng -q probe.fa -o pack -O probe.pack
impg map -a c4.syng -q probe.fa -o proj -O probe.proj
```

HPRCv2 AGC build pattern:

```bash
impg syng \
  --agc /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -o /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  --syncmer-length 63 --smer-length 8 --syncmer-seed 7 \
  --position-sample-rate 256 \
  --parallel-dictionary \
  -t 32
```

HPRCv2 C4/C4-style local graph:

```bash
impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o gfa:syng:crush \
  -O c4.syng.crush \
  -v 1
```

HPRCv2 local SYNG parameter sweep style:

```bash
impg query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -b loci.bed \
  -d 100k \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -o 'gfa:syng-local:blunt,k=127,s=16,seed=7:crush' \
  -O local-k127-s16-graphs \
  --render-graph --render-graph-output local-k127-s16-renders \
  --describe-graph --graph-report-output local-k127-s16-reports \
  --graph-report-format tsv \
  -t 16 -v 1
```

Pack/genotype/infer:

```bash
impg map \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -q sample.reads.fq.zst \
  -o proj -O sample.proj \
  --pack-compression-level 12 \
  -t 32

impg genotype cos \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  --proj sample.proj \
  -r GRCh38#0#chr6:31891045-32123783 \
  --ploidy 2 --top-n 20

impg infer \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  --proj sample.proj \
  --partitions partitions.bed \
  --top-n 1 \
  -O sample.infer.tsv.zst
```

Whole-index GFA dump:

```bash
impg syng2gfa \
  -a c4.syng \
  --sequence-files chr6.C4.fa \
  --gfa-mode blunt \
  --gfa-version 1.0 \
  -o c4.syncmers.gfa
```

## Limitations And Open Issues

Be explicit about validation scope:

- The core build/query/map paths have extensive integration tests in
  `tests/test_syng_integration.rs` and targeted tests for start-count/rank
  behavior, graph output, pack/genotype semantics, and parser shorthands.
- HPRCv2 and C4 examples in docs/scripts are real developer-machine workloads,
  not CI-reproducible fixtures. They should be treated as operational evidence,
  not portable acceptance tests.
- `gfa:syng:crush` is useful and actively validated for path preservation, but
  it is not PGGB-quality by default on every hard locus. Current work is about
  reducing fragmentation and repeat-glue artifacts while preserving exact path
  spellings.
- Raw SYNG GFA is topology-incomplete for ordinary variation-graph consumers:
  it preserves the syncmer overlap topology/path spellings but does not induce
  bubbles/flubbles the way seqwish/PGGB-style graph induction does.
- The fast local-compression testbed shows path preservation and topology
  assertions passing for compact/window methods on small CI fixtures, but many
  method IDs currently route through the same compact-bubble placeholder. That
  testbed does not yet prove resolver-specific superiority.
- C4/local graph quality is still iterative. Existing notes identify C4 as the
  canonical hard workload and track syng frequency masks, local rebuild
  parameters, crush methods, SmoothXG-like smoothing, and SweepGA/seqwish
  replacement experiments. This is current engineering evidence, not a settled
  production guarantee.
- Boundary refinement is ends-only. Interior CIGAR/base-level alignment through
  every chain is not implemented in the SYNG query path.
- `impg map -o paf` is syncmer-anchor projection, not base alignment. Default
  `gaf` and pack outputs do not load coordinate sidecars and intentionally do
  not assign reads to genomic intervals.
- Pack evidence is distinct-node-per-read support. It is not base depth, not
  mapping-quality weighted, and not occurrence-count depth within a read.
- `genotype cos` and `infer` are implemented scoring/inference tools over
  graph-feature vectors, but they are not full pangenome imputation. Read-walk
  GAF evidence currently affects stitching/mosaic paths, not local cosine rows.
- Ambiguous bases are mapped to `A` during syncmer extraction. Graph GFA gap
  materialization can preserve or cut `N` in fetched gap DNA, but syncmer
  selection itself follows the SYNG numeric-base convention.
- Old `.names` files without `first_syncmer_pos` can load, but fresh rebuilds
  are preferred for exact coordinate behavior.
