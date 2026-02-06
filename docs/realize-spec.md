# Realize Engine Specification

## Overview

The realize engine converts implicit pangenome graph regions into explicit GFA variation graphs through recursive decomposition. Given a set of homologous sequence intervals (from `impg query`), it produces a high-quality pangenome graph by:

1. **Large regions (>threshold)**: Running sweepga all-vs-all alignment, then recursively partitioning into smaller chunks
2. **Small regions (<=threshold)**: Running SPOA partial order alignment with boundary padding
3. **Lacing**: Stitching sub-graphs back together into a single coherent GFA

This replaces the current single-shot approaches (`gfa-poa` and `gfa-seqwish`) with a multi-resolution strategy that produces better graphs for regions of any size.

## Motivation

The current graph generation has two modes, each with limitations:

- **gfa-poa** (`generate_gfa_from_intervals`): Uses SPOA directly. Works well for small regions (~1kb) but degrades on larger regions because POA is O(n*L^2) and struggles with structural variation.
- **gfa-seqwish** (`generate_gfa_seqwish_from_intervals`): Uses sweepga+seqwish. Handles large regions but produces coarse graphs that miss fine-grained variation in short segments.

The realize engine combines both: sweepga for global structure, SPOA for local consensus, with recursive decomposition bridging the two scales.

## Architecture

```
realize(intervals, sequences, depth=0)
    |
    |-- if total_span <= poa_threshold:
    |       return padded_poa(intervals, sequences)
    |
    |-- alignments = sweepga_align(sequences, sparsification[depth])
    |-- chunks = partition(alignments, chunk_size)
    |
    |-- for each chunk:
    |       sub_intervals = extract_chunk_intervals(chunk)
    |       sub_gfa = realize(sub_intervals, sequences, depth+1)
    |       collect sub_gfa
    |
    |-- return lace(sub_gfas)
```

### Module Layout

```
src/realize/
    mod.rs          -- Public API: realize() entry point and RealizeConfig
    orchestrator.rs -- Recursive decomposition logic
    poa.rs          -- Padded POA with boundary handling
    align.rs        -- sweepga alignment execution (pairs -> PAF)
    partition.rs    -- Chunk partitioning from alignment results
    lace.rs         -- Sub-graph stitching (thin wrapper around existing lace)
```

## Core Types

```rust
/// Configuration for the realize engine
pub struct RealizeConfig {
    /// Regions at or below this size (bp) go directly to POA.
    /// Above this, sweepga + recursive decomposition is used.
    /// Default: 1000
    pub poa_threshold: usize,

    /// Target chunk size (bp) when partitioning large regions.
    /// Chunks will be approximately this size after sweepga alignment.
    /// Default: 5000
    pub chunk_size: usize,

    /// Boundary padding (bp) added on each side of POA regions.
    /// Overlapping padded regions are trimmed during lacing.
    /// Default: 100
    pub padding: usize,

    /// Maximum recursion depth. Safety limit to prevent runaway recursion.
    /// Default: 10
    pub max_depth: usize,

    /// Number of threads for parallel operations.
    pub num_threads: usize,

    /// Sparsification strategy per recursion depth level.
    /// Index 0 = depth 0 (coarsest), higher indices = deeper levels.
    /// If depth exceeds the length, the last entry is reused.
    /// Default: [Connectivity(0.99)]
    pub sparsification: Vec<SparsificationStrategy>,

    /// SPOA scoring parameters: (match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2)
    /// Default: (5, 4, 6, 2, 24, 1)
    pub scoring_params: (u8, u8, u8, u8, u8, u8),

    /// SeqwishConfig for sweepga+seqwish alignment steps
    pub seqwish_config: SeqwishConfig,

    /// Optional temp directory for intermediate files
    pub temp_dir: Option<String>,

    /// Whether to sort the final GFA using gfasort
    pub sort_output: bool,
}
```

```rust
/// A region to be realized: a set of homologous sequence intervals
pub struct RealizeRegion {
    /// The intervals (from impg query results)
    pub intervals: Vec<Interval<u32>>,

    /// Genomic span of the anchor region (the query target)
    pub anchor_start: i32,
    pub anchor_end: i32,

    /// Recursion depth (0 = top level)
    pub depth: usize,
}
```

```rust
/// Result of realizing a region
pub struct RealizeResult {
    /// The GFA graph as a string
    pub gfa: String,

    /// Statistics about the realization
    pub stats: RealizeStats,
}

pub struct RealizeStats {
    /// Total sequences in the region
    pub num_sequences: usize,
    /// Maximum recursion depth reached
    pub max_depth_reached: usize,
    /// Number of POA leaf calls
    pub poa_calls: usize,
    /// Number of sweepga alignment calls
    pub sweepga_calls: usize,
    /// Total time in milliseconds
    pub total_ms: u64,
}
```

## Algorithm Detail

### 1. Entry Point: `realize()`

```rust
pub fn realize(
    impg: &impl ImpgIndex,
    intervals: &[Interval<u32>],
    sequence_index: &UnifiedSequenceIndex,
    config: &RealizeConfig,
) -> io::Result<RealizeResult>
```

The public entry point. Extracts sequences from intervals using `prepare_sequences()` (already in `src/graph.rs`), then delegates to the recursive orchestrator.

**Sequence extraction**: Uses the existing `prepare_sequences()` function which handles:
- Parallel sequence fetching via rayon
- Strand detection and reverse-complement
- MAF-style coordinate normalization
- Longest-first sorting

### 2. Recursive Orchestrator

```rust
fn realize_recursive(
    sequences: &[(String, SequenceMetadata)],
    config: &RealizeConfig,
    depth: usize,
    stats: &mut RealizeStats,
) -> io::Result<String>
```

**Decision logic**:

```
let total_span = sequences.iter().map(|(s, _)| s.len()).max().unwrap_or(0);

if total_span <= config.poa_threshold || depth >= config.max_depth {
    // Base case: POA
    return padded_poa(sequences, config);
}

// Recursive case: align, partition, recurse, lace
let strategy = config.sparsification
    .get(depth)
    .or(config.sparsification.last())
    .cloned()
    .unwrap_or(SparsificationStrategy::Connectivity(0.99));

let alignments = sweepga_align(sequences, &strategy, config)?;
let chunks = partition_into_chunks(sequences, &alignments, config.chunk_size)?;

let sub_gfas: Vec<String> = chunks
    .into_iter()
    .map(|chunk| realize_recursive(&chunk.sequences, config, depth + 1, stats))
    .collect::<Result<Vec<_>, _>>()?;

lace_subgraphs(&sub_gfas, sequences, config)
```

**Why `max().unwrap_or(0)` for span**: The longest sequence determines whether the region is "large" since it dominates alignment complexity.

**Depth-varying sparsification**: At the top level (depth 0), we want sparser alignments for speed. At deeper levels with smaller regions, we can afford denser alignments. Example configuration:

```rust
vec![
    SparsificationStrategy::Connectivity(0.99),  // depth 0: sparse
    SparsificationStrategy::Connectivity(0.999), // depth 1: denser
    SparsificationStrategy::None,                // depth 2+: all pairs
]
```

### 3. Padded POA (Base Case)

```rust
fn padded_poa(
    sequences: &[(String, SequenceMetadata)],
    config: &RealizeConfig,
) -> io::Result<String>
```

For small regions, run SPOA with boundary padding to produce clean graph edges.

**Padding strategy**:

```
Original region:     [==========]
With padding:    [pp|==========|pp]
                  ^                ^
             padding bp       padding bp

POA alignment:   [pp|==========|pp]  <- all sequences padded symmetrically
Trim result:         [==========]    <- remove padding columns from GFA
```

**Algorithm**:

1. For each sequence, extend coordinates by `config.padding` bp on each side (clamped to sequence boundaries).
2. Fetch extended sequences.
3. Run SPOA alignment (using existing `prepare_poa_graph_and_sequences` logic).
4. Generate GFA from the SPOA graph.
5. Trim padding: identify and remove graph nodes that are exclusively within padding regions.
   - Walk each path from the start; nodes visited before reaching the non-padded region are "left-padding" nodes.
   - Walk each path from the end; nodes visited before reaching the non-padded region are "right-padding" nodes.
   - If a node is padding for ALL paths, remove it. Otherwise keep it (it's part of the core alignment).
6. Update paths and edges after trimming.

**Why padding matters**: Without padding, POA graph edges at region boundaries are artifacts of sequence truncation rather than true variation. Padding lets SPOA "see" flanking context, producing correct graph topology at boundaries. The padding is then trimmed from the final output.

**Edge case**: When a sequence is near a contig boundary (start < padding or end + padding > contig_length), clamp the extension. Short asymmetric padding is acceptable.

### 4. Sweepga Alignment Execution

```rust
fn sweepga_align(
    sequences: &[(String, SequenceMetadata)],
    strategy: &SparsificationStrategy,
    config: &RealizeConfig,
) -> io::Result<Vec<PafRecord>>
```

Runs sweepga all-vs-all alignment on the input sequences, filtered by the sparsification strategy.

**Implementation approach**: This wraps the existing sweepga integration from `src/commands/graph.rs`:

1. Write sequences to a temporary FASTA file.
2. Generate alignment pairs using `generate_pairs()` from `src/commands/align.rs`.
3. For each pair, run `FastGAIntegration::align_to_temp_paf()`.
4. Apply sweepga filtering (plane-sweep with scaffolding) unless `config.no_filter`.
5. Parse resulting PAF records.

**Key difference from `impg align`**: The existing `align` command generates job lists but doesn't execute alignments directly (the `TODO` at line 662 of `align.rs`). The realize engine needs to actually run the alignments. This is the `complete-impg-align` dependency task.

**Integration with sparsification**: At each recursion level, pair selection is controlled by the depth-indexed sparsification strategy. This means:
- Depth 0 (large regions, many sequences): Use connectivity-based sparsification to keep O(n log n) pairs
- Deeper levels (smaller regions, fewer sequences): Can afford denser or all-pairs alignment

### 5. Chunk Partitioning

```rust
fn partition_into_chunks(
    sequences: &[(String, SequenceMetadata)],
    alignments: &[PafRecord],
    chunk_size: usize,
) -> io::Result<Vec<Chunk>>
```

After sweepga alignment, partition the aligned region into smaller chunks for recursive processing.

**Partitioning approach**: Use the longest (anchor) sequence as a coordinate system:

1. Identify the anchor sequence (longest, or the original query target).
2. Divide the anchor into windows of `chunk_size` bp.
3. For each window, project through the sweepga alignments to find the corresponding intervals in all other sequences.
4. Each chunk is a set of subsequences covering one window.

**Chunk overlap**: Adjacent chunks overlap by `config.padding` bp to ensure clean boundaries when laced together.

```
Anchor:  [----chunk0----|--pad--|----chunk1----|--pad--|----chunk2----]
                        ^^^^^^^^              ^^^^^^^^
                        overlap               overlap
```

**Output**:

```rust
struct Chunk {
    /// Subsequences for this chunk (one per input sequence, possibly fewer if
    /// some sequences don't span this region)
    sequences: Vec<(String, SequenceMetadata)>,

    /// Position in the anchor coordinate system
    anchor_start: i32,
    anchor_end: i32,

    /// Index for ordering during lacing
    index: usize,
}
```

### 6. Sub-graph Lacing

```rust
fn lace_subgraphs(
    sub_gfas: &[String],
    sequences: &[(String, SequenceMetadata)],
    config: &RealizeConfig,
) -> io::Result<String>
```

Stitches chunk sub-graphs back into a single GFA.

**This wraps the existing lace infrastructure** from `src/commands/lace.rs`. The key operations are:

1. Parse each sub-GFA into segments, links, and paths.
2. Remap node IDs to a unified namespace (existing lace logic handles this).
3. Trim overlapping padding regions between adjacent chunks.
4. Link contiguous ranges (where chunk N ends and chunk N+1 begins).
5. Emit combined GFA.

**Overlap trimming**: Because chunks overlap by `padding` bp, the overlapping regions must be reconciled. The lace module already handles this via its overlap detection and node-splitting logic (`trim_overlaps` in `lace.rs`).

**Path continuity**: Paths must be continuous across chunk boundaries. The lace module's gap-filling logic (`fill_gaps`) inserts connecting nodes where needed.

**Final sort**: If `config.sort_output` is true, sort the combined GFA using `sort_gfa()` (gfasort Ygs pipeline).

## Integration Points

### Wire into `impg query --format gfa`

Replace the current single-shot GFA generation in `perform_query()`:

```rust
// Current (single-shot):
"gfa" | "gfa-seqwish" => generate_gfa_seqwish_from_intervals(...)
"gfa-poa" => generate_gfa_from_intervals(...)

// New (realize engine):
"gfa" | "gfa-realize" => realize(impg, &results, &sequence_index, &realize_config)
"gfa-seqwish" => generate_gfa_seqwish_from_intervals(...)  // keep as fallback
"gfa-poa" => generate_gfa_from_intervals(...)              // keep as fallback
```

The default `gfa` format switches from seqwish to the realize engine. The old modes remain accessible as `gfa-seqwish` and `gfa-poa`.

### Wire into `impg graph`

The `graph` command (`src/commands/graph.rs`) currently runs the sweepga+seqwish pipeline directly. With the realize engine, it can optionally use recursive decomposition:

```rust
// In build_graph():
if config.use_realize {
    // Use realize engine for graph construction
    realize(impg, &intervals, &sequence_index, &realize_config)
} else {
    // Existing sweepga+seqwish pipeline (default for now)
    generate_gfa_seqwish_from_intervals(...)
}
```

### Wire into `impg similarity`

The similarity command could use realize-generated graphs for more accurate pairwise similarity computation, comparing graph-based similarity against alignment-based similarity.

## CLI Parameters

New flags for commands that use the realize engine:

```
--poa-threshold <BP>       POA threshold in bp [default: 1000]
--chunk-size <BP>          Chunk size for recursive decomposition [default: 5000]
--padding <BP>             Boundary padding in bp [default: 100]
--max-realize-depth <N>    Maximum recursion depth [default: 10]
--sparsification <SPEC>    Sparsification strategy per depth level.
                           Can be specified multiple times for depth-varying.
                           [default: giant:0.99]
--no-realize               Disable realize engine, use direct seqwish/POA
```

## Error Handling

- **Empty regions**: If a chunk has no sequences (all sequences are shorter than the anchor at that position), skip it and let lacing fill the gap.
- **Single-sequence chunks**: A chunk with one sequence produces a trivial GFA (one segment, one path). This is correct and handled by both POA and seqwish.
- **Recursion limit**: At `max_depth`, always use POA regardless of region size. Log a warning if the region is still large.
- **sweepga failure**: If alignment fails for a chunk, fall back to POA for that chunk. Log the error.
- **Lacing failure**: If sub-graphs cannot be laced (e.g., no shared paths), concatenate them with gap nodes.

## Performance Considerations

- **Parallelism**: sweepga alignment and sequence fetching are already parallel (rayon). Chunk processing within a recursion level can also be parallelized, but nested parallelism should use rayon's scope to avoid thread explosion.
- **Memory**: Each recursion level materializes sequences and sub-GFAs. For very deep recursion, this could be significant. The `max_depth` limit and `poa_threshold` prevent pathological cases.
- **Temp files**: sweepga requires temp files for FASTA and PAF. Use `config.temp_dir` (e.g., `/dev/shm`) for better I/O performance.
- **Complexity**: For a region of size L with n sequences:
  - Single-level (POA only): O(n * L^2) - acceptable for L <= 1kb
  - Single-level (seqwish): O(n^2 * L) for alignment + O(n * L) for graph induction
  - Realize (recursive): O(n^2 * L / log(L)) amortized, with much better constants due to smaller sub-problems

## Testing Strategy

1. **Unit tests**: Test padded_poa with known small sequences. Verify padding is correctly trimmed.
2. **Integration tests**: Run realize on yeast chromosome regions at various sizes (100bp, 1kb, 10kb, 100kb). Compare path content against known haplotypes.
3. **Regression tests**: Verify that realize produces graphs containing all input sequences as paths (no sequence loss).
4. **Boundary tests**: Test with regions at contig boundaries, single-sequence regions, and regions with only two sequences.
5. **Round-trip test**: Extract sequences from realize GFA paths, align back to originals, verify identity.

## Dependency Tasks

This spec unblocks:

1. **complete-impg-align**: Implement actual sweepga execution in `align.rs` (currently only generates job lists). The realize engine needs `sweepga_align()` to return PAF records directly.
2. **implement-padded-poa**: Implement the padded POA with boundary padding and trimming logic in `src/realize/poa.rs`.
3. **implement-recursive-realize**: Implement the recursive orchestrator, chunk partitioning, and lacing in `src/realize/mod.rs` and `src/realize/orchestrator.rs`.
4. **wire-realize-into** (x3): Wire the realize engine into query, graph, and similarity commands.
