# impg depth

## What it does

`impg depth` computes pangenome coverage depth: for every locus in the pangenome, it counts how many distinct haplotypes (samples) cover that position. The output tells you which regions are shared across many genomes (high depth) versus unique to a single assembly (depth=1).

This is useful for:

- **Identifying conserved regions** shared by all or most assemblies
- **Detecting sample-specific sequence** (depth=1 regions, e.g., novel insertions)
- **Quality-checking pangenome alignments** (unexpected depth drops may indicate alignment gaps)
- **Guiding downstream analyses** like variant calling or structural variant detection in regions with sufficient sample coverage

## Core idea

Traditional depth tools (like `samtools depth`) require reads mapped to a single reference. In a pangenome, there is no single reference -- you have N assemblies with pairwise alignments between them. `impg depth` solves this by treating pairwise alignments as an implicit graph and traversing it to find all haplotypes covering each locus, without building an explicit graph structure.

The key insight: if sample A aligns to sample B, and sample B aligns to sample C, then A, B, and C all cover the same locus. `impg depth` discovers these transitive relationships through BFS traversal of the alignment graph.

## Basic usage

By default, `impg depth` computes depth across the entire pangenome. Every base is assigned to exactly one output row -- no double-counting. Hub sequences (those with the most alignments) are auto-detected via a degree pre-scan and processed first as anchors.

```bash
impg depth -a alignments.paf -O output_prefix -t 64
```

Output format (`output_prefix.depth.tsv`):

```
#id  length  depth  SampleA          SampleB           SampleC
1    50000   3      chrA:0-50000     chrB:100-50100    chrC:200-50200
2    30000   1      chrA:50000-80000 NA                NA
```

Each row is an interval. The `depth` column is the number of distinct samples covering it. Each sample column shows `seq:start-end` coordinates or `NA`. Gaps between alignment intervals are filled with depth=1 (the anchor sequence only).

## Modes

### 1. Ref-anchored mode (`--ref SAMPLE`)

Prioritizes the specified sample's sequences as anchors. This guarantees that all regions covered by the ref sample are reported in the ref sample's coordinate system.

```bash
impg depth -a alignments.paf --ref CHM13 -O output_prefix
```

The ref sample's sequences are processed in Phase 1 (before all others), so they claim all their homologous regions first. Non-ref sequences only anchor regions the ref sample doesn't cover. All sequences from all samples still appear in the output.

### 2. Ref-only mode (`--ref SAMPLE --ref-only`)

Only outputs depth for regions anchored on the ref sample's sequences. Regions not covered by the ref sample are omitted from the output entirely.

```bash
impg depth -a alignments.paf --ref CHM13 --ref-only -O output_prefix
```

This is useful when you want a reference-centric view: "for each position on CHM13, how many other haplotypes align here?"

### 3. Region query mode (`-r` or `-b`)

Query depth for specific genomic regions rather than the whole pangenome. This mode queries both forward and reverse alignment directions (via `query_reverse_for_depth`) to ensure complete coverage.

```bash
# Single region
impg depth -a alignments.paf -r chr1:1000000-2000000

# Multiple regions from BED file
impg depth -a alignments.paf -b regions.bed
```

Output format:

```
#ref_seq   ref_start  ref_end  depth  Sample1            Sample2
chr1       1000000    1050000  3      sample2_seq:100-200  sample3_seq:50-150
```

Note: `--samples` and `--samples-file` filtering only applies in region query mode. In global mode, these options are loaded but not passed to the depth computation pipeline.

### Add-on: Statistics output (`--stats`)

The `--stats` flag can be combined with the global or ref-anchored modes to produce summary statistics and per-depth BED files instead of the default tabular output. Cannot be combined with region query mode (`-r`/`-b`).

```bash
# Basic statistics
impg depth -a alignments.paf --stats -O depth_stats

# Statistics with per-interval sample tracking
impg depth -a alignments.paf --stats --combined-output -O depth_stats

# Statistics anchored on a reference
impg depth -a alignments.paf --stats --ref CHM13 -O depth_stats
```

Outputs:
- `depth_stats.summary.txt`: depth distribution (depth N covers X bp, Y%)
- `depth_stats.depthN.bed`: BED file for each depth level
- `depth_stats.combined.bed` (with `--combined-output`): single sorted BED with depth and sample columns

Combined BED format:

```
#seq_name  length  depth  samples
chrA       50000   3      sample1,sample2,sample3
chrA       30000   1      sample1
```

Note: the `length` column is `end - start` (the interval width), not the full sequence length. When `--merge-tolerance` is set, adjacent intervals with similar depth are merged (depth = max of merged group, samples = union).

## How the algorithm works

### Step 1: Degree pre-scan

Before computing depth, `impg depth` scans all sequences in parallel to compute their **alignment degree** -- the number of distinct other samples each sequence has direct alignments to. Only sequences passing the `--min-seq-length` filter are considered.

```
Example (star topology: B,C,D,E all aligned to A):
  A: degree=4  (aligned to B, C, D, E)
  B: degree=1  (aligned to A)
  C: degree=1  (aligned to A)
  D: degree=1  (aligned to A)
  E: degree=1  (aligned to A)
```

This pre-scan identifies **hub sequences** (high-degree nodes) that should be processed first. The pre-scan runs in parallel and is fast (seconds for thousands of sequences).

When `--ref` is specified, the ref sample's sequences are used as Phase 1 directly, but the degree pre-scan still runs to sort the remaining sequences by connectivity.

### Step 2: Sequence ordering

Sequences are sorted for processing by:

1. `--ref` sample sequences first (if specified)
2. Alignment degree descending (hubs before leaves)
3. Sequence length descending (longer sequences first, as tiebreaker)

### Step 3: Two-phase parallel processing

**Phase 1 -- Hub sequences**: High-connectivity sequences are processed first, guaranteeing they become the coordinate-system anchors for their homologous regions.

- With `--ref`: the ref sample's sequences are Phase 1
- Without `--ref`: sequences with degree >= ceil(max_degree/2) are auto-detected as hubs (computed as `(max_degree + 1) / 2`). If max_degree <= 1, Phase 1 is skipped entirely.

For transitive mode (`-x`), Phase 1 uses **chunk-level parallelism**: each hub sequence is split into 5MB chunks, giving ~1200 parallel tasks for a human-sized genome -- fully utilizing all threads. For non-transitive mode, Phase 1 uses sequence-level parallelism (each hub sequence is one task).

**Phase 2 -- Remaining sequences**: Processes the rest. Most regions are already claimed by Phase 1 hubs, so Phase 2 only anchors regions the hubs didn't cover.

### Step 4: Per-anchor depth computation

For each anchor sequence (or unprocessed region within it):

**Non-transitive (default):**
Uses a fast path with pre-scanned raw alignment intervals. All raw intervals for the anchor sequence are loaded via `query_raw_intervals()`, sorted by target start, and searched with binary search -- O(log n) to find the relevant range, then linear scan over overlapping intervals. The bidirectional index captures both alignment directions, so a single 1-hop query finds all directly aligned samples.

**Transitive (`-x`) — default raw-interval BFS:**
Uses a depth-specific raw-interval BFS (`depth_transitive_bfs`) that operates on raw alignment extents via `query_raw_overlapping()` at each hop -- an O(n+k) coitrees range query that fetches only intervals overlapping the current region. Query coordinates are computed via linear interpolation (same as the non-transitive path). This avoids CIGAR walk overhead and chunk boundary gaps from CIGAR-precise projection, making it faster and more robust for depth computation where base-level precision is not needed. Starting from the anchor region, follows alignment chains up to `--max-depth` hops (default: 2). Large regions are split into 5MB chunks to prevent BFS from exploring the entire pangenome. Visited ranges are lazily allocated per sequence (not pre-allocated for all sequences), keeping memory footprint minimal.

**Transitive with `--use-BFS` — CIGAR-precise BFS:**
Uses the same CIGAR-based BFS as the `query` command (`query_transitive_bfs`/`query_transitive_dfs`), providing base-level precision at the cost of CIGAR walk overhead. Use this when exact coordinate projection is needed rather than region-level depth.

```
Example with -x and max_depth=2:
  Anchor: A (hub)
  Hop 1: A -> B, A -> C, A -> D (direct alignments)
  Hop 2: B -> E (E aligns to B but not to A)
  Result: depth=5 (A + B + C + D + E)
```

All discovered regions on other sequences are projected back to the anchor's coordinate system via chained linear interpolation.

**Self-alignment filtering**: Alignments where the query sample matches the anchor sample but the query sequence ID differs (intra-genome duplications) are filtered out. Same-sample, same-sequence alignments (the anchor covering itself) are explicitly added as depth=1.

**BFS visited tracking**: A region is considered "visited" if >80% of it is already covered by previously visited ranges on the same sequence. This prevents infinite loops and redundant work.

### Step 5: Sweep-line depth computation

Once all overlapping alignments are gathered for an anchor region, a sweep-line algorithm computes per-position depth:

1. Create events: (position, START/END, sample_id) for each alignment
2. Sort events by position, with starts before ends at the same position
3. Sweep left-to-right, maintaining a `SampleBitmap` of active samples
4. At each position change, emit a `SparseDepthInterval` with the current depth and sample set
5. Fill gaps between intervals with depth=1 (anchor sample only)

The `SampleBitmap` structure tracks active samples using 1 bit per sample plus a u32 overlap count per sample. Multiple overlapping alignments from the same sample are correctly tracked (counted as depth=1 for that sample, not double-counted). When multiple alignments from the same sample overlap, the one with the largest overlap with the current interval is chosen for coordinate reporting.

**Sparse storage**: `SparseDepthInterval` only stores samples that are present (as `Vec<(sample_id, query_seq_id, query_start, query_end)>`), sorted by sample_id for O(log n) binary search lookup. This is far more memory-efficient than storing all samples per interval.

### Step 6: Output formatting

After sweep-line produces intervals:

1. **Window splitting** (only when `--window-size` is explicitly provided): Intervals crossing window boundaries are cut into pieces, with sample query coordinates proportionally adjusted.
2. **Merge tolerance** (when `--merge-tolerance > 0`): Adjacent intervals with similar depth are merged. An interval can merge with the current group if `(max_depth - min_depth) / max_depth <= tolerance`. When merged: depth = max of group, samples = union.
3. **Output**: In normal mode, each interval becomes a TSV row with `#id length depth Sample1 Sample2 ...`. In stats mode, intervals are accumulated into depth distribution histograms and per-depth BED files.

### Step 7: Processed region tracking

A `ConcurrentProcessedTracker` ensures every base in the pangenome is counted exactly once:

- When an anchor discovers alignments on other sequences, those regions are marked as "processed"
- When a sequence comes up for processing, only its unprocessed regions are computed
- Per-sequence mutexes provide zero-contention parallelism (different sequences never share a lock)
- After Phase 1, all hub sequences are explicitly marked as fully processed

## Why hub-first ordering matters

Consider a star topology where samples B, C, D, E are all aligned to sample A:

```
Without hub-first ordering:
  B processed first -> finds A -> depth=2 (B+A)
  A is now "processed" -> skipped
  C processed -> finds A -> depth=2 (C+A)
  D, E same -> depth=2 each
  WRONG: true depth should be 5

With hub-first ordering:
  A processed first -> finds B,C,D,E -> depth=5 (A+B+C+D+E)
  B,C,D,E marked as processed -> skipped
  CORRECT: depth=5 reported in A's coordinate system
```

When `--ref` is not specified, the degree pre-scan automatically detects A as the hub (degree=4 >> others' degree=1) and processes it first. When `--ref A` is specified, the same result is guaranteed explicitly.

## Key options

### Transitive queries

| Option | Description |
|--------|-------------|
| `-x` | Enable transitive BFS traversal (richer depth for sparse alignments) |
| `--transitive-dfs` | Use DFS instead of BFS (conflicts with `-x`; slower, fewer overlapping results) |
| `-m, --max-depth N` | Maximum BFS/DFS hops (default: 2, 0 = unlimited) |
| `--min-transitive-len N` | Minimum alignment length for transitive hops (default: 100) |
| `--min-distance-between-ranges N` | Minimum distance between ranges on same sequence (default: 10) |

Without `-x`, only direct (1-hop) alignments are considered. This is fast but may underestimate depth when the anchor is not the hub. With `-x`, BFS discovers indirect connections at the cost of more computation.

### Reference selection

| Option | Description |
|--------|-------------|
| `--ref SAMPLE` | Prioritize this sample's sequences as Phase 1 anchors |
| `--ref-only` | Only output rows anchored on the ref sample (requires `--ref`) |

`--ref` ensures the named sample's sequences are processed first, so they become coordinate-system anchors for all regions they cover. All sequences are still included in the output. Adding `--ref-only` filters the output to only include the ref sample's anchored rows.

When `--ref` is not specified, hub sequences are auto-detected via the degree pre-scan.

### Sequence filtering

| Option | Description |
|--------|-------------|
| `--min-seq-length SIZE` | Exclude sequences shorter than SIZE (accepts k/m/g suffixes) |
| `--fai-list FILE` | FAI index files; sequences not in alignments get depth=1 |

### Sample filtering (region query mode only)

| Option | Description |
|--------|-------------|
| `--samples LIST` | Comma-separated sample names to include |
| `--samples-file FILE` | File with sample names to include (one per line) |

Note: these options only take effect in region query mode (`-r`/`-b`). In global mode, the filter is loaded but not passed to the depth computation.

### Output control

| Option | Description |
|--------|-------------|
| `-O, --output-prefix PREFIX` | Output file prefix (stdout if omitted) |
| `--window-size N` | Split output into fixed-size windows. Only applied when explicitly set; no windowing by default. |
| `--merge-tolerance F` | Merge adjacent intervals with similar depth (default: 0.05 = 5%). Formula: `(max - min) / max <= tolerance`. When merged: depth = max, samples = union. |
| `--merge-adjacent` | Merge adjacent results with identical depth (region query mode) |
| `--separator CHAR` | PanSN format separator for sample name extraction (default: `#`) |

### Statistics add-on

| Option | Description |
|--------|-------------|
| `--stats` | Produce summary statistics and per-depth BED files (conflicts with `-r`/`-b`) |
| `--combined-output` | Single BED file with depth and sample columns (requires `--stats`) |

### Performance

| Option | Description |
|--------|-------------|
| `-t N` | Number of threads |
| `--use-BFS` | Use CIGAR-precise BFS for transitive depth (slower but base-level precision). Default uses raw-interval BFS with linear interpolation, which is faster and sufficient for region-level depth. |
| `--index-mode per-file` | Per-file indexing for many alignment files |

### Memory management

In transitive mode, interval tree caching is enabled (`set_tree_cache_enabled(true)`) to reuse trees across BFS nodes within the same sequence. The sub-index cache is cleared after each sequence (`clear_sub_index_cache()`) to bound memory usage. With jemalloc, this avoids RSS fragmentation.

## Data structures

Key data structures used for memory-efficient depth computation:

| Structure | Purpose | Memory |
|-----------|---------|--------|
| `SampleIndex` | Maps sample names to compact u16 IDs (max 65535 samples) | O(num_samples) |
| `CompactSequenceLengths` | Sequence lengths, sample-to-seq and seq-to-sample mappings using numeric IDs | O(num_sequences) |
| `CompactAlignmentInfo` | Alignment record using numeric IDs instead of strings | ~60% less than string-based |
| `SampleBitmap` | Tracks active samples in sweep-line: 1 bit per sample + u32 overlap count | O(num_samples) per instance |
| `SparseDepthInterval` | Depth interval storing only present samples | 24 bytes base + 16 bytes per sample |
| `IntervalSet` | Sorted non-overlapping intervals with binary search operations | O(num_intervals) |
| `ConcurrentProcessedTracker` | Per-sequence mutexes wrapping IntervalSets | O(num_sequences) |

## Limitations

- **Non-transitive mode with star topologies**: Without `-x`, leaf nodes only see the hub (depth=2 instead of N). Hub-first ordering mitigates this for the hub itself, but leaf-anchored regions still report incomplete depth. Use `-x` for full coverage.
- **Transitive depth limit**: With `-x -m 2` (default), chains longer than 2 hops may miss distant samples. Increase `-m` for deeper traversal at the cost of more computation.
- **Coordinate projection is approximate by default**: Default raw-interval BFS uses linear interpolation (not CIGAR-level mapping). Coordinates in non-anchor sample columns are approximate. Use `--use-BFS` for CIGAR-precise mapping, though this may produce more discontinuities at chunk boundaries due to revealing actual alignment structure differences.
- **Parallel non-determinism within Phase 2**: Within Phase 2, rayon's work-stealing means the exact assignment of anchors can vary slightly between runs. Phase 1 (hub) results are deterministic.
- **Sample filtering scope**: `--samples` and `--samples-file` only work in region query mode. Global mode does not support sample-level filtering.
- **BFS visited heuristic**: The 80% coverage threshold for "visited" is a heuristic. In rare cases with highly fragmented overlaps, a region may be revisited or skipped prematurely.

## Examples

```bash
# Global depth with auto hub detection
impg depth -a all_vs_all.paf -O pangenome -t 64

# Transitive depth with reference anchoring
impg depth -a all_vs_all.paf -x --ref CHM13 -O pangenome -t 128

# Ref-only: depth from CHM13's perspective only
impg depth -a all_vs_all.paf -x --ref CHM13 --ref-only -O pangenome -t 128

# Large-scale with per-file indexing and length filter
impg depth --alignment-list files.list --index-mode per-file \
    --ref CHM13 -x --min-seq-length 1M -O result -t 128

# Global depth with statistics output
impg depth -a all_vs_all.paf --stats --combined-output -O stats

# Query specific region
impg depth -a all_vs_all.paf -x -r chr1:10000000-20000000

# Filter to specific samples (region query mode only)
impg depth -a all_vs_all.paf --samples HG002,HG003,CHM13 -r chr1:10000000-20000000

# DFS traversal instead of BFS
impg depth -a all_vs_all.paf --transitive-dfs --ref CHM13 -O pangenome -t 128

# With explicit window splitting
impg depth -a all_vs_all.paf --ref CHM13 --window-size 100000 -O pangenome -t 64

# CIGAR-precise BFS mode (slower, base-level precision)
impg depth -a all_vs_all.paf -x --use-BFS --ref CHM13 -O pangenome -t 128
```
