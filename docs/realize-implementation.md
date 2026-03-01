# Realize Engine — Implementation Reference

Deep code review and documentation of the recursive realize engine as implemented on branch `smooth-it`.

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Module Map](#module-map)
3. [Data Flow](#data-flow)
4. [Entry Points](#entry-points)
5. [Recursive Orchestrator](#recursive-orchestrator)
6. [Padded POA (Base Case)](#padded-poa-base-case)
7. [Sweepga Alignment](#sweepga-alignment)
8. [Chunk Partitioning](#chunk-partitioning)
9. [Sub-graph Lacing](#sub-graph-lacing)
10. [Configuration Reference](#configuration-reference)
11. [Known Issues and Edge Cases](#known-issues-and-edge-cases)
12. [Deviations from Spec](#deviations-from-spec)

---

## Architecture Overview

The realize engine converts a set of homologous sequence intervals into a GFA variation graph through recursive multi-resolution decomposition:

```
realize(intervals)
├── prepare_sequences()          # Extract, strand-normalize, sort longest-first
└── realize_recursive(seqs, depth=0)
    ├── [base case] padded_poa_from_sequences()   # depth ≥ max OR span ≤ threshold
    └── [recursive case]
        ├── sweepga_align()                       # all-vs-all alignment → PAF
        ├── parse_paf_file()                      # PAF → PafRecord structs
        ├── partition_into_chunks()               # project anchor windows → chunks
        ├── realize_recursive(chunk, depth+1)     # recurse on each chunk
        └── lace_subgraphs()                      # ID remap + path merge + linking
```

Key design decisions:
- **Anchor-centric partitioning**: The longest sequence serves as the coordinate reference for windowing.
- **Linear interpolation projection**: PAF alignments are used to project anchor windows onto other sequences via proportional coordinate mapping.
- **Padding for boundary quality**: POA base cases include flanking context that is trimmed via MSA column analysis.
- **Fallback cascade**: Empty alignments → POA; single chunk → POA; single sequence → trivial GFA.

## Module Map

| File | Lines | Purpose |
|------|-------|---------|
| `src/realize/mod.rs` | 974 | Public API, recursive orchestrator, chunk partitioning, lacing |
| `src/realize/poa.rs` | 1098 | Padded POA: SPOA alignment with boundary padding and MSA-based trimming |
| `src/commands/align.rs` | 1195 | Sweepga alignment execution, sparsification strategies, PAF filtering |
| `src/graph.rs` | 1086 | Shared utilities: `prepare_sequences()`, `sort_gfa()`, `SequenceMetadata`, GFA-to-MSA |

## Data Flow

### Type Pipeline

```
Interval<u32>                      # impg query result (seq_id, start, end)
    ↓ prepare_sequences()
(String, SequenceMetadata)         # extracted sequence + coordinate metadata
    ↓ realize_recursive()
    ├─ padded_poa_from_sequences() → PaddedPoaResult { gfa, metadata }
    │  via SPOA → MSA → trim → re-POA → GFA
    └─ sweepga_align()             → NamedTempFile (PAF)
       ↓ parse_paf_file()
       PafRecord[]                 # parsed PAF alignment records
       ↓ partition_into_chunks()
       Chunk[]                     # windowed subsequences with anchor coords
       ↓ realize_recursive() per chunk
       String[] (sub-GFAs)
       ↓ lace_subgraphs()
       String (merged GFA)
```

### Key Type: `SequenceMetadata`
```rust
// src/graph.rs:24
pub struct SequenceMetadata {
    pub name: String,         // Contig name (e.g. "chr1")
    pub start: i32,           // MAF-style start (RC-adjusted for '-' strand)
    pub size: i32,            // Span in bp
    pub strand: char,         // '+' or '-'
    pub total_length: usize,  // Full contig length
}
```

`start` uses MAF convention: for reverse strand, `start = total_length - end` (forward coordinates).

### Key Type: `PafRecord`
```rust
// src/realize/mod.rs:119
struct PafRecord {
    query_name: String,
    query_len: usize,       // #[allow(dead_code)]
    query_start: usize,
    query_end: usize,
    target_name: String,
    target_len: usize,       // #[allow(dead_code)]
    target_start: usize,
    target_end: usize,
    strand: char,            // #[allow(dead_code)]
}
```

Note: `query_len`, `target_len`, and `strand` are parsed but currently unused (marked `dead_code`).

## Entry Points

### `realize()` — Primary API
**File**: `src/realize/mod.rs:142`

```rust
pub fn realize(
    impg: &impl ImpgIndex,
    intervals: &[Interval<u32>],
    sequence_index: &UnifiedSequenceIndex,
    config: &RealizeConfig,
) -> io::Result<RealizeResult>
```

Steps:
1. Calls `prepare_sequences()` (from `graph.rs`) to extract, RC-normalize, and sort sequences longest-first.
2. If empty, returns a header-only GFA.
3. Creates atomic counters for `poa_calls`, `sweepga_calls`, `max_depth_reached`.
4. Calls `realize_recursive()` at depth 0.
5. Optionally sorts the final GFA via `sort_gfa()` (gfasort Ygs pipeline).
6. Returns `RealizeResult` with GFA string and stats.

### `realize_from_sequences()` — Sequences-in-Memory API
**File**: `src/realize/mod.rs:195`

Same as `realize()` but takes pre-extracted `&[(String, SequenceMetadata)]` directly. Used by modules that already have sequences in memory (e.g., testing, or chained pipeline stages).

## Recursive Orchestrator

**File**: `src/realize/mod.rs:244`

```rust
fn realize_recursive(
    sequences: &[(String, SequenceMetadata)],
    config: &RealizeConfig,
    depth: usize,
    poa_calls: &AtomicUsize,
    sweepga_calls: &AtomicUsize,
    max_depth_reached: &AtomicUsize,
) -> io::Result<String>
```

### Decision Logic

1. **Max depth tracking**: Uses `compare_exchange_weak` CAS loop to atomically track the maximum recursion depth reached (lines 253-264). Thread-safe for potential future parallel chunk processing.

2. **Empty input**: Returns `"H\tVN:Z:1.0\n"`.

3. **Single sequence**: Returns a trivial GFA (one segment `S`, one path `P`) via `make_trivial_gfa()`.

4. **Base case** (`total_span <= poa_threshold` OR `depth >= max_depth`):
   - Logs a warning if max depth forced POA on a large region.
   - Calls `padded_poa_from_sequences()`.

5. **Recursive case** (`total_span > poa_threshold` AND `depth < max_depth`):
   - Selects sparsification strategy from `config.sparsification[depth]` (or last entry if depth exceeds array length).
   - Builds `SweepgaAlignConfig` via `build_sweepga_config()`.
   - Formats sequence names as `"{name}:{start}-{end}({strand})"` for PAF naming.
   - Calls `sweepga_align()` → PAF temp file.
   - Parses PAF records via `parse_paf_file()`.
   - **Fallback**: If no alignments produced, falls back to POA with a warning.
   - Calls `partition_into_chunks()` to divide the region.
   - **Fallback**: If partitioning yields ≤1 chunk, falls back to POA.
   - Recurses on each chunk at `depth + 1`.
   - Laces sub-GFAs via `lace_subgraphs()`.

### Fallback Chain

```
sweepga produces 0 alignments  →  POA (warn)
partition yields ≤1 chunk       →  POA
chunk has 0 sequences          →  skip
sub_gfas is empty              →  header-only GFA
sub_gfas has 1 entry           →  return it directly (no lacing)
```

**Important**: The sparsification strategy is currently fetched but stored in `_strategy` (unused variable, line 295). The `build_sweepga_config()` always uses hardcoded defaults. This means depth-varying sparsification is **not yet implemented** — only the alignment config defaults are used.

## Padded POA (Base Case)

**File**: `src/realize/poa.rs`

### Two Entry Points

1. **`padded_poa()`** (line 40): Takes impg + intervals + sequence_index. Fetches sequences with extended coordinates (clamped to contig boundaries), tracks actual padding amounts, handles RC padding swap.

2. **`padded_poa_from_sequences()`** (line 119): Takes pre-extracted sequences. Estimates padding as `min(padding, seq_len/2)` for left and `min(padding, seq_len - left_pad)` for right. This is the entry point used by the realize engine.

### Algorithm

```
Input sequences (with padding baked in)
    ↓
Build SPOA graph (kSW local alignment, convex gap model)
    ↓
Check if any_padding is true
    ├── false: generate GFA directly, post-process strand orientations
    └── true:
        ↓ generate MSA
        ↓ find_core_column_range() — identify non-padding MSA columns
        ↓ build_trimmed_gfa() — extract core subsequences, re-run SPOA, emit GFA
```

### SPOA Configuration

- Alignment type: `kSW` (Smith-Waterman local alignment)
- Scoring: Convex gap model with two gap penalty tiers
- Default: match=5, mismatch=4, gap_open1=6, gap_ext1=2, gap_open2=24, gap_ext2=1
- Sequences fed longest-first for best SPOA quality
- Weights: uniform (`vec![1u32; seq.len()]`)

### Core Column Range Detection

**`find_core_column_range()`** (line 306):

For each MSA row, counts non-gap characters to determine where padding ends and core begins:

- **`core_start`** = max over all sequences of the column where left padding is exhausted
- **`core_end`** = min over all sequences of the column where right padding would begin

This ensures the core region contains only "real" aligned bases from all sequences.

**Degenerate case**: If `core_start >= core_end` (padding consumed everything), returns `(mid, mid)` where `mid = ncols / 2`. This produces an empty trimmed result.

### Trimming Strategy

**`build_trimmed_gfa()`** (line 387):

Rather than surgically editing the SPOA graph, this function:
1. Extracts non-gap characters from MSA columns `[core_start, core_end)` for each sequence
2. Builds a **fresh** SPOA graph from only the trimmed core sequences
3. Generates GFA from the new graph

This is cleaner than trying to remove nodes from an existing graph but means the POA is run **twice** when padding is present.

### Strand Post-Processing

**`post_process_gfa_for_strands()`** (line 472 in poa.rs, also line 70 in graph.rs):

For reverse-strand paths, reverses the path step order and flips orientations (`+` ↔ `-`). The implementation in `poa.rs` is a near-duplicate of the one in `graph.rs`, with one difference: `graph.rs` panics on missing orientation, `poa.rs` defaults to `seg.to_string()`.

### Path Naming Convention

Headers use forward-strand coordinates: `"{name}:{fwd_start}-{fwd_end}"`. For reverse strand metadata, the conversion is:
```
fwd_start = total_length - start - size
fwd_end   = total_length - start
```

## Sweepga Alignment

**File**: `src/commands/align.rs:924`

```rust
pub fn sweepga_align(
    sequences: &[(String, &[u8])],
    config: &SweepgaAlignConfig,
) -> io::Result<tempfile::NamedTempFile>
```

### Steps

1. If fewer than 2 sequences, returns an empty PAF temp file.
2. Writes all sequences to a temporary FASTA file (`.fa` suffix required by FastGA).
3. Sets `TMPDIR` environment variable if `config.temp_dir` is specified.
4. Runs `FastGAIntegration::align_to_temp_paf()` (all-vs-all: same file as both query and target).
5. Applies PAF filtering via `apply_paf_filter()` unless `config.no_filter`.
6. Returns the (filtered) PAF as a `NamedTempFile`.

### PAF Filtering Pipeline

**`apply_paf_filter()`** (line 815):

Uses sweepga's `PafFilter` with:
- Filter mode parsed from `"1:1"` string → `FilterMode` + per-query/per-target limits
- Scaffold chaining: jump=50kb, mass=10kb, filter=1:1
- Overlap threshold: 0.95
- Scoring function: `LogLengthIdentity`
- Self-alignments excluded (`.with_keep_self(false)`)

### Default Alignment Parameters (from realize engine)

```rust
// build_sweepga_config() at realize/mod.rs:386
kmer_frequency: 10,
min_alignment_length: 100,
num_mappings: "1:1",
scaffold_jump: 50_000,
scaffold_mass: 10_000,
scaffold_filter: "1:1",
overlap: 0.95,
min_identity: 0.0,
```

## Chunk Partitioning

**File**: `src/realize/mod.rs:451`

```rust
fn partition_into_chunks(
    sequences: &[(String, SequenceMetadata)],
    paf_records: &[PafRecord],
    chunk_size: usize,
    padding: usize,
) -> Vec<Chunk>
```

### Algorithm

1. **Anchor selection**: `sequences[0]` (pre-sorted longest-first).

2. **Short-circuit**: If anchor length ≤ `chunk_size + padding`, returns a single chunk with all sequences.

3. **Name formatting**: Sequences are named `"{name}:{start}-{end}({strand})"` to match PAF record names from sweepga.

4. **Projection collection**: For each non-anchor sequence, collects all PAF records where the anchor is either query or target. Stores as `(anchor_start, anchor_end, seq_start, seq_end)` tuples, sorted by `anchor_start`.

5. **Window generation**: Creates overlapping windows along the anchor:
   ```
   pos=0: [0-padding .. 0+chunk_size+padding]
   pos=chunk_size: [chunk_size-padding .. 2*chunk_size+padding]
   ...
   ```
   Note: First window starts at `max(0, pos - padding)`, so effectively position 0 for the first window.

6. **Sequence extraction per window**: For each window `[win_start, win_end)`:
   - Always includes the anchor's slice for this window.
   - For each non-anchor sequence, calls `project_window_to_sequence()`.
   - Clamps projected coordinates to sequence length.
   - Sorts chunk sequences by length descending.

### Coordinate Projection

**`project_window_to_sequence()`** (line 585):

For each alignment block overlapping the anchor window, performs **linear interpolation**:

```
fraction_start = (overlap_start - anchor_start) / anchor_span
fraction_end   = (overlap_end   - anchor_start) / anchor_span
projected_start = seq_start + fraction_start * seq_span
projected_end   = seq_start + fraction_end   * seq_span
```

Across multiple overlapping alignment blocks, takes the union: `min(all projected_starts)` to `max(all projected_ends)`.

Returns `None` if no alignment covers the window (sequence is absent from that chunk).

**Limitation**: Linear interpolation assumes uniform alignment density within each PAF block. Insertions/deletions within a block cause the projected coordinates to be approximate. For the purpose of defining chunk boundaries (which include padding), this is usually acceptable.

**Limitation**: Strand information from PAF records is currently ignored during projection. The `strand` field is parsed but not used. This means reverse-strand alignments are projected as if they were forward-strand, which could produce incorrect coordinate projections for inverted regions.

## Sub-graph Lacing

**File**: `src/realize/mod.rs:635`

```rust
fn lace_subgraphs(sub_gfas: &[String], _num_threads: usize) -> io::Result<String>
```

### Algorithm

1. **ID Remapping**: Each sub-GFA gets fresh node IDs in a unified namespace. `next_id` starts at 1 and increments for each segment.

2. **GFA Parsing**: Parses each sub-GFA line-by-line:
   - `S` lines: Remap segment IDs, store in `id_remap` HashMap.
   - `L` lines: Remap both endpoints using `id_remap`.
   - `P` lines: Remap each step's node ID, collect per-path ordered steps.

3. **Path Merging**: For each path name (across all sub-GFAs):
   - Concatenates the step lists from each sub-GFA where the path appears.
   - Between consecutive sub-GFA contributions, inserts a linking edge: `L {last_node} {last_orient} {first_node} {first_orient} 0M`.

4. **Link Deduplication**: Uses a `HashSet<&String>` to remove duplicate links.

5. **Assembly**: Emits `H` header, all segments, deduplicated links, and merged paths.

### Path Naming and Ordering

Path order is preserved by tracking first-appearance order across sub-GFAs. A path that appears in sub-GFA 0 and sub-GFA 2 (but not 1) gets its steps from both concatenated in order.

### Linking Edges

The `0M` overlap in linking edges indicates zero-base overlap between adjacent chunks. This is correct when chunks are non-overlapping in the final graph (padding has been trimmed by the POA base case).

**Note**: The `_num_threads` parameter is accepted but unused. Lacing is currently single-threaded.

## Configuration Reference

### `RealizeConfig`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `poa_threshold` | `usize` | 1000 | Max region size (bp) for direct POA |
| `chunk_size` | `usize` | 5000 | Target window size for partitioning |
| `padding` | `usize` | 100 | Boundary padding for POA regions |
| `max_depth` | `usize` | 10 | Maximum recursion depth |
| `num_threads` | `usize` | 4 | Thread count for parallel operations |
| `sparsification` | `Vec<SparsificationStrategy>` | `[Connectivity(0.99)]` | Per-depth sparsification (last entry reused) |
| `scoring_params` | `(u8,u8,u8,u8,u8,u8)` | `(5,4,6,2,24,1)` | SPOA scoring: match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2 |
| `seqwish_config` | `SeqwishConfig` | `default()` | Seqwish parameters (currently unused by realize) |
| `temp_dir` | `Option<String>` | `None` | Temp directory for intermediate files |
| `sort_output` | `bool` | `true` | Whether to sort final GFA via gfasort |

### Sparsification Strategies

| Strategy | Syntax | Description |
|----------|--------|-------------|
| `None` | `"none"` / `"all"` | All pairs — O(n^2) |
| `Random(f)` | `"random:0.5"` | Random fraction of pairs |
| `Connectivity(p)` | `"giant:0.99"` | Erdos-Renyi giant component guarantee |
| `TreeSampling{..}` | `"tree:2:1:0.1:15"` | k-nearest + k-farthest + random (Mash-distance based) |

## Known Issues and Edge Cases

### 1. Sparsification strategy is fetched but not applied

**Location**: `src/realize/mod.rs:295-300`

The depth-indexed sparsification strategy is retrieved into `_strategy` (note the underscore — unused). The `build_sweepga_config()` function always uses hardcoded defaults. This means all recursion depths use the same alignment parameters regardless of the `sparsification` config vector.

### 2. Strand not used in coordinate projection

**Location**: `src/realize/mod.rs:499-519`

When building projections from PAF records, the strand field is parsed (line 421) but not used in `project_window_to_sequence()`. For reverse-complement alignments, the query coordinates in PAF are on the reverse strand, but the projection treats them as forward-strand coordinates. This could produce incorrect chunk boundaries for inverted alignments.

### 3. Duplicate `post_process_gfa_for_strands` implementation

**Locations**: `src/graph.rs:70` and `src/realize/poa.rs:472`

Near-identical implementations exist in both files. The `graph.rs` version panics on missing segment orientation; the `poa.rs` version silently defaults to the raw step string. Should be consolidated.

### 4. Double POA execution with padding

**Location**: `src/realize/poa.rs:119-204`

When `padding > 0`, `padded_poa_from_sequences()` builds the SPOA graph once (to get the MSA for trimming column detection), then discards it and builds a **second** SPOA graph from the trimmed core sequences. This doubles the POA computational cost for padded regions.

### 5. Padding estimation in `padded_poa_from_sequences`

**Location**: `src/realize/poa.rs:141-154`

The realize engine calls `padded_poa_from_sequences()` with `padding = config.padding` (default 100). The function estimates left/right padding as `min(padding, seq_len/2)` and `min(padding, seq_len - left_pad)`. For sequences much larger than `2*padding`, this works correctly. But the padding estimate assumes the caller already embedded padding into the sequences — it doesn't verify this.

### 6. Large-padding degenerate case

**Location**: `src/realize/poa.rs:372-377` and MEMORY.md gotcha

When `padding > seq_len/2`, both `left_pad` and `right_pad` consume the entire sequence, producing `core_start >= core_end`. The function returns `(mid, mid)` — an empty core. The test at line 907 (`test_padded_poa_large_padding_clamped`) confirms this is handled without panicking, but the result is essentially empty.

The MEMORY.md documents this: "padded_poa_from_sequences with padding > seq_len/2 produces empty core region (use padding=0 for tiny seqs)."

### 7. `SequenceMetadata.start` is `i32`, coordinates may overflow for very large genomes

**Location**: `src/graph.rs:25`

`start` and `size` are `i32`, which limits them to ~2.1 Gbp. Human chromosomes fit, but some plant genomes could overflow. The `total_length` field is `usize` (no overflow).

### 8. Path name format mismatch between realize and POA

The realize engine formats names as `"{name}:{start}-{end}({strand})"` (with strand in parentheses) for sweepga alignment naming (line 308). But `make_trivial_gfa()` and `make_headers()` format path names as `"{name}:{fwd_start}-{fwd_end}"` (without strand). These two naming conventions coexist in the same GFA when lacing combines trivial GFAs from single-sequence chunks with POA-generated GFAs.

### 9. Link deduplication uses `HashSet<&String>` (order not deterministic)

**Location**: `src/realize/mod.rs:785`

The link deduplication in `lace_subgraphs()` uses a `HashSet`, which means link output order is non-deterministic across runs. This doesn't affect correctness but makes GFA output non-reproducible for diffing/testing.

### 10. `seqwish_config` field is unused

**Location**: `src/realize/mod.rs:53`

`RealizeConfig` carries a `seqwish_config` field but the realize engine never uses seqwish — it uses sweepga for alignment and POA for graph construction. This field is a leftover from the spec.

### 11. Chunk processing is sequential

**Location**: `src/realize/mod.rs:353-367`

Chunks are processed in a `for` loop, not in parallel. The atomic counters suggest parallel processing was anticipated (they'd be unnecessary for single-threaded execution), but rayon parallelism hasn't been added. For large regions with many chunks, this could be a performance bottleneck.

## Deviations from Spec

The original spec is in `docs/realize-spec.md`. Notable deviations in the implementation:

| Spec | Implementation | Notes |
|------|---------------|-------|
| Separate modules: `orchestrator.rs`, `partition.rs`, `lace.rs`, `align.rs` | All in `realize/mod.rs` except POA | Simpler structure; no separate lace/partition modules |
| `RealizeRegion` struct | Not implemented | Sequences passed directly, not wrapped in a region struct |
| Wraps existing lace infrastructure from `commands/lace.rs` | Reimplements lacing from scratch | Self-contained implementation in `lace_subgraphs()` |
| Depth-varying sparsification active | Sparsification fetched but unused (`_strategy`) | All depths use same alignment config |
| `sweepga_align` returns `Vec<PafRecord>` | Returns `NamedTempFile` (PAF on disk) | Parsed separately via `parse_paf_file()` |
| Overlap trimming via lace module | `0M` linking edges, no overlap trimming | Padding is trimmed at POA level, not during lacing |
| Parallel chunk processing mentioned | Sequential `for` loop | Atomic counters suggest parallelism was planned |
| Wired into `impg query --format gfa` | Not yet wired | Exists as standalone `realize()` / `realize_from_sequences()` API |

## Test Coverage

### `realize/mod.rs` tests

| Test | What it covers |
|------|---------------|
| `test_make_trivial_gfa` | Single-sequence GFA generation with correct path naming |
| `test_realize_from_sequences_empty` | Empty input → header-only GFA |
| `test_realize_from_sequences_single` | Single sequence → trivial GFA, 1 POA call |
| `test_realize_from_sequences_small_region_uses_poa` | Two short sequences below threshold → direct POA, 0 sweepga calls |
| `test_project_window_to_sequence_no_overlap` | Window before/after alignment → None |
| `test_project_window_to_sequence_full_overlap` | Window exactly covering alignment → exact projection |
| `test_project_window_to_sequence_partial_overlap` | Window covering half → proportional projection |
| `test_partition_single_chunk` | Short sequences → single chunk returned |
| `test_lace_subgraphs_single` | Single GFA → returned unchanged |
| `test_lace_subgraphs_two` | Two GFAs → ID remapping, linking edge, merged path |
| `test_lace_subgraphs_empty` | Empty input → header-only GFA |

### `realize/poa.rs` tests

| Test | What it covers |
|------|---------------|
| `test_find_core_column_range_no_padding` | Zero padding → full MSA is core |
| `test_find_core_column_range_symmetric_padding` | Equal padding on both sides → correct core range |
| `test_find_core_column_range_with_gaps` | MSA gaps interact correctly with padding counting |
| `test_find_core_column_range_asymmetric_padding` | Different padding amounts per sequence → conservative core |
| `test_find_core_column_range_empty_msa` | Empty input → (0, 0) |
| `test_find_core_column_range_all_padding` | Padding consumes everything → degenerate midpoint |
| `test_padded_poa_from_sequences_no_padding` | Two identical sequences, padding=0 → valid GFA |
| `test_padded_poa_from_sequences_with_padding` | Core region with flanking padding → valid trimmed GFA |
| `test_padded_poa_from_sequences_with_variation` | SNP in core → multiple segments (bubble) |
| `test_padded_poa_empty_input` | Empty → empty GFA |
| `test_padded_poa_single_sequence` | One sequence → valid single-path GFA |
| `test_padded_poa_reverse_strand_metadata` | Strand preserved in metadata |
| `test_padded_poa_large_padding_clamped` | Huge padding → no panic (graceful clamp) |
| `test_padded_poa_three_sequences_with_indel` | Indel variation → bubble graph |
| `test_post_process_gfa_for_strands_forward` | Forward strand path unchanged |
| `test_post_process_gfa_for_strands_reverse` | Reverse strand → path reversed, orientations flipped |
| `test_make_headers_forward_and_reverse` | Header coordinate conversion for both strands |

### Missing test coverage

- **No integration test with sweepga**: All tests use sequences below `poa_threshold`, so the recursive case (sweepga + partition + lace) is never exercised in tests.
- **No test for `parse_paf_file()`**: PAF parsing is only tested indirectly.
- **No test for multi-chunk lacing with different path names**: Only tests single shared path name.
- **No test for `project_window_to_sequence` with multiple overlapping alignment blocks**.
- **No test for chunk extraction with sequences absent from some chunks**.
