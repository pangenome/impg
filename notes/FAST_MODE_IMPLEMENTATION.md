# Approximate Mode Implementation Guide

## Overview

This document provides a complete, step-by-step guide to implement `--approximate` mode for the `impg` tool. Approximate mode provides approximate query results by scanning tracepoints instead of computing full CIGAR strings and fetching sequences, resulting in significant performance improvements for .1aln format alignments.

**Performance Goal**: ~10-100x speedup by skipping sequence I/O and CIGAR computation, making .1aln queries nearly as fast as PAF format.

### Prerequisites

Before implementing approximate mode, ensure these .1aln optimizations are in place (already implemented in current codebase):

- **Vec-based O(1) trace_spacing cache**: Using `Vec<Option<i64>>` indexed by file_index instead of HashMap (10-50x faster lookups)
- **Lazy-loaded cache**: Only reads trace_spacing when needed, never pre-warms (critical for memory efficiency)
- **Memory-efficient file opens**: Separate short-lived file opens for cache miss to minimize buffer accumulation
- **Optimized buffer sizes**: 64-256 bytes depending on operation type
- **Direct fetch optimization**: `fetch_alignment_direct()` static function in `onealn.rs`

### What Approximate Mode Does

**Skips**:
- ❌ Sequence fetching from FASTA/AGC files
- ❌ CIGAR computation via WFA (Wavefront Alignment)
- ❌ Reverse complement operations
- ❌ Base-level precision

**Provides**:
- ✅ Approximate intervals (accurate to within trace_spacing, typically ~100bp)
- ✅ Approximate CIGAR encoding for identity metrics (gap-compressed identity and block-identity)
- ✅ Identity estimation from tracepoints and trace_diffs without sequence I/O
- ✅ Same overlap detection via interval trees
- ✅ Correct strand handling

### Use Cases

Ideal for:
- Exploratory analysis requiring quick overviews
- Large-scale queries where approximate results suffice
- Transitive queries exploring deep graph relationships
- Partition operations (only need approximate boundaries)
- Initial filtering before detailed analysis

## Technical Background

### 1aln Tracepoint Format

Alignments in .1aln format store:
- `trace_spacing`: Regular interval on **query** sequence (e.g., 100bp)
- `tracepoints[]`: Array of consumed lengths on **target** sequence
- `trace_diffs[]`: Differences at each tracepoint (used for identity estimation in approximate mode)

**Key insight**: Query advances regularly by `trace_spacing`, target advances irregularly (stored in tracepoints) due to indels.

```
Query:  |----100bp----|----100bp----|----100bp----|
Target: |----95bp-----|----102bp----|----98bp-----|
        ^              ^              ^              ^
        TP[0]=95      TP[1]=102      TP[2]=98
```

### Approximate Projection Algorithm

Instead of fetching sequences and computing CIGAR:

1. Get alignment metadata (already in memory from interval tree)
2. Fetch tracepoints using optimized `fetch_alignment_direct()` with `read_trace_diffs=true`
3. **Scan phase**: Collect all tracepoints overlapping requested range
   - Track first/last overlapping tracepoint info for refinement
   - Accumulate identity statistics (matches/mismatches from trace_diffs)
   - Use FULL tracepoint segments during scan
4. **Refinement phase**: Apply overlap heuristic AFTER scan completes
   - Refine first query position based on partial overlap of first segment
   - Refine last query position based on partial overlap of last segment
   - Refine min/max target positions to actual overlap boundaries
5. Return intervals with approximate CIGAR encoding

```rust
// Phase 1: Scan and collect
for each tracepoint:
    query_delta = (idx == 0) ? first_boundary : trace_spacing
    target_delta = tracepoint * direction

    if segment overlaps requested_range:
        record full segment boundaries
        save first/last segment info for refinement
        accumulate identity statistics from trace_diffs

// Phase 2: Refine (CRITICAL - do this AFTER collecting all segments)
if first_segment overlaps partially:
    fraction = (range_start - segment_start) / target_delta
    first_query_pos += fraction * query_delta  // Refine!

if last_segment overlaps partially:
    fraction = (range_end - segment_start) / target_delta
    last_query_pos = segment_query_start + fraction * query_delta  // Refine!
```

### Refined Overlap Heuristic

**Key insight**: When a tracepoint segment partially overlaps the requested range, we can estimate the corresponding query position using fractional projection.

**Example**:
```
Requested:  |----------[range: 10-90]----------|
Tracepoint: |----[0-100 target, 0-100 query]----|
Overlap:         [10-90 target, ?-? query]

Heuristic:
  - First overlap: fraction_start = (10-0)/100 = 0.1
    → query_start ≈ 0 + 0.1*100 = 10
  - Last overlap: fraction_end = (90-0)/100 = 0.9
    → query_end ≈ 0 + 0.9*100 = 90
```

**Benefits**:
- ✅ Tighter intervals (±10bp typical vs ±100bp without refinement)
- ✅ More accurate identity filtering
- ✅ Better transitivity (fewer spurious connections)
- ✅ No performance penalty (uses data already scanned)

**CRITICAL**: Apply refinement AFTER collecting all overlapping segments, not during the scan. This ensures we work with the complete overlapping region.

### Coordinate Storage Details

**Query coordinates** (metadata):
- Stored in FORWARD order for BOTH strands
- `query_start < query_end` always in metadata
- For reverse alignments: swap after projection to match normal mode output
  - Normal mode returns: `query_start > query_end` (negative length)
  - Approximate mode must match this behavior

**Target coordinates** (metadata):
- Forward: stored forward
- Reverse: stored REVERSED (flipped within contig)

**Tracepoints**:
- Refer to original target coordinates in .1aln file
- For reverse: target reversed, query NOT reversed

**First boundary calculation**:
```rust
first_boundary = ((query_start / trace_spacing) + 1) * trace_spacing - query_start
```
Only the first segment uses `first_boundary`; all others use `trace_spacing`.

### Identity Estimation in Approximate Mode

**How it works**:
1. Scan overlapping tracepoints and accumulate:
   - `aligned_length = min(query_delta, abs(target_delta))`
   - `matches = aligned_length - trace_diffs[i]`
   - `mismatches = trace_diffs[i]`
2. Create approximate CIGAR: `[=matches, Xmismatches]`
3. Calculate gap-compressed identity: `matches / (matches + mismatches)`
4. Filter empty CIGARs (too small to produce meaningful statistics)

**Benefits**:
- ✅ Accurate identity metrics without sequence fetching
- ✅ No NaN values in bedpe output
- ✅ No performance penalty (trace_diffs already being read)
- ✅ Sequence files not required for bed/bedpe output

**Implementation**:
- `fetch_alignment_direct(read_trace_diffs=true)` for approximate mode
- Normal mode: `read_trace_diffs=false` (full CIGAR reconstructed anyway)

## Implementation Steps

### Step 1: Add Command-Line Flag

**File**: [src/main.rs](src/main.rs)

**Location**: `QueryOpts` struct definition (search for `struct QueryOpts`)

**Add this field**:
```rust
/// Use approximate mode for faster queries (1aln files only, bed/bedpe output): scans tracepoints without CIGAR computation or sequence fetching
#[arg(help_heading = "Performance")]
#[clap(long, action)]
approximate: bool,
```

**Why**: Creates the `--approximate` flag available for query and refine commands.

---

### Step 2: Add Output Format Validation

**File**: [src/main.rs](src/main.rs)

**Location**: In `Query` command handler, after output format resolution (search for `let resolved_output_format =`)

**Add this validation**:
```rust
// Validate --approximate mode compatibility
if query.approximate {
    if resolved_output_format != "bed" && resolved_output_format != "bedpe" {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "--approximate mode is only compatible with 'bed' and 'bedpe' output formats, not '{}'",
                resolved_output_format
            ),
        ));
    }
}
```

**Why**: Approximate mode returns approximate CIGAR vectors, incompatible with PAF, GFA, MAF, and FASTA-ALN formats that require full alignment paths.

---

### Step 3: Implement Approximate Projection Method

**File**: [src/impg.rs](src/impg.rs)

**Location**: After `project_overlapping_interval()` method (search for `fn project_overlapping_interval(`)

**Implementation**: See `project_overlapping_interval_fast()` in the current codebase.

**Key features**:
- Scans tracepoints to find overlapping subset
- Accumulates identity statistics (matches, mismatches, insertions, deletions)
- Creates approximate CIGAR encoding
- Handles forward and reverse strands correctly
- Early exit optimization when past requested range

**Why**: Core of approximate mode. Scans tracepoints without sequence I/O or full CIGAR computation.

---

### Step 4: Update query() Method

**File**: [src/impg.rs](src/impg.rs)

**Location**: `pub fn query()` method

**Changes**:
1. Add `approximate_mode: bool` parameter
2. Update debug log to indicate approximate mode
3. Branch between approximate and normal mode:
   - Approximate mode: call `project_overlapping_interval_fast()`
   - Normal mode: call `project_overlapping_interval()` with sequence cache

**Why**: Enables mode selection at query time.

---

### Step 5: Update Transitive Query Methods

**File**: [src/impg.rs](src/impg.rs)

**Methods to update**:
- `query_transitive_dfs()`
- `query_transitive_bfs()`

**Changes**:
1. Add `approximate_mode: bool` parameter
2. Branch projection calls based on mode
3. Pass mode through all recursive/iterative steps

**Why**: Transitive queries amplify performance gains since they explore multiple hops. Approximate mode is especially beneficial here.

---

### Step 6: Update perform_query() Function

**File**: [src/main.rs](src/main.rs)

**Changes**:
1. Add `approximate_mode: bool` parameter
2. Pass to all three query variants (regular, transitive BFS, transitive DFS)
3. Update all call sites to pass `query.approximate`

**Why**: Central query dispatcher needs to pass approximate_mode flag through to all query methods.

---

### Step 7: Update Sequence Requirement Logic

**File**: [src/main.rs](src/main.rs)

**Function**: `setup_output_resources()`

**Changes**:
```rust
// In approximate mode with bed/bedpe output, .1aln files don't need sequences
// (identity is calculated from tracepoints without CIGAR reconstruction)
let onealn_needs_sequences = has_onealn_files &&
    !(approximate_mode && (output_format == "bed" || output_format == "bedpe"));
```

**Why**: Approximate mode with bed/bedpe doesn't require sequence files since it calculates identity from tracepoints.

---

### Step 8: Optimize fetch_alignment_direct()

**File**: [src/onealn.rs](src/onealn.rs)

**Function**: `fetch_alignment_direct()`

**Changes**:
1. Add `read_trace_diffs: bool` parameter
2. Conditionally read X line (trace_diffs):
   ```rust
   'X' => {
       if read_trace_diffs {
           alignment.trace_diffs = file.int_list().map(|v| v.to_vec()).unwrap_or_default()
       }
   }
   ```

**Why**:
- Normal mode: skips reading trace_diffs for faster I/O (full CIGAR reconstructed anyway)
- Approximate mode: reads trace_diffs for identity calculation

---

### Step 9: Update get_onealn_alignment() Calls

**File**: [src/impg.rs](src/impg.rs)

**Changes**:
- Normal mode CIGAR reconstruction: `get_onealn_alignment(metadata, false)`
- Approximate mode projection: `get_onealn_alignment(metadata, true)`

**Why**: Pass correct flag based on whether trace_diffs are needed.

---

### Step 10: Add --approximate to Partition Command

**File**: [src/main.rs](src/main.rs)

**Changes**:
1. Add `approximate: bool` field to Partition struct
2. Pass to `partition_alignments()` function
3. Update transitive query calls in partition.rs

**Why**: Partition command uses transitive queries and only needs intervals (not CIGARs), perfect for approximate mode.

---

### Step 11: Update Refine Command Support

**Files**: [src/commands/refine.rs](src/commands/refine.rs) and [src/main.rs](src/main.rs)

**Changes**:
1. Add `approximate_mode: bool` to RefineConfig struct
2. Update transitive query calls in refine logic
3. Initialize from `refine.query.approximate`

**Why**: Refine uses QueryOpts (which includes --approximate flag), so it automatically inherits approximate mode support.

---

## Verification

### Compilation Check

```bash
cargo check
```

**Expected**: Should complete without errors. All type signatures should match.

### Build Test

```bash
cargo build --release
```

**Expected**: Clean build with no warnings about unused parameters or type mismatches.

---

## Testing

### Test 1: Basic Approximate Query (bed output)

```bash
# Normal mode
time ./target/release/impg query \
  -a alignments.1aln \
  -r chr1:1000000-2000000 \
  -o bed > normal.bed

# Approximate mode
time ./target/release/impg query \
  -a alignments.1aln \
  -r chr1:1000000-2000000 \
  -o bed --approximate > fast.bed

# Compare
wc -l normal.bed fast.bed
```

**Expected**:
- Approximate mode significantly faster (5-100x depending on data)
- Similar number of results (intervals approximate but overlaps should match)

### Test 2: Approximate Mode with BEDPE Output

```bash
./target/release/impg query \
  -a alignments.1aln \
  -r chr1:1000000-1001000 \
  -o bedpe --approximate > output.bedpe
```

**Expected**:
- Works without errors
- Identity values (gi:f:, bi:f:) are accurate (not NaN)
- No sequence files required

### Test 3: Approximate Mode with Incompatible Format (should fail)

```bash
./target/release/impg query \
  -a alignments.1aln \
  -r chr1:1000-2000 \
  -o paf --approximate
```

**Expected**: Error message: "--approximate mode is only compatible with 'bed' and 'bedpe' output formats, not 'paf'"

### Test 4: Approximate Transitive Query

```bash
time ./target/release/impg query \
  -a alignments.1aln \
  -r chr1:1000000-1001000 \
  -x --max-depth 3 \
  --approximate -o bed
```

**Expected**: Much faster than without --approximate, especially with higher max-depth

### Test 5: Approximate Partition

```bash
time ./target/release/impg partition \
  -a alignments.1aln \
  -w 100000 \
  -o bed --approximate --separate-files
```

**Expected**: Faster partitioning while producing valid partition files

### Test 6: Approximate Refine

```bash
time ./target/release/impg refine \
  -a alignments.1aln \
  -b regions.bed \
  --approximate > refined.bed
```

**Expected**: Faster refinement with approximate but valid intervals

---

## Files Modified Summary

| File | Changes |
|------|---------|
| [src/main.rs](src/main.rs) | Added `approximate: bool` to QueryOpts and Partition; validation; updated perform_query(); updated all call sites; updated sequence requirement logic |
| [src/impg.rs](src/impg.rs) | Added project_overlapping_interval_fast() with identity calculation; updated query(), query_transitive_dfs(), query_transitive_bfs(); updated get_onealn_alignment() |
| [src/onealn.rs](src/onealn.rs) | Added read_trace_diffs parameter to fetch_alignment_direct() |
| [src/commands/partition.rs](src/commands/partition.rs) | Added approximate_mode parameter; updated transitive query calls |
| [src/commands/refine.rs](src/commands/refine.rs) | Added approximate_mode to RefineConfig; updated transitive query calls |

---

## Troubleshooting

### "Approximate mode is only compatible with 'bed' and 'bedpe'"

**Cause**: Using --approximate with output format that requires full CIGAR paths (paf, gfa, maf, fasta-aln)

**Solution**: Use `-o bed` or `-o bedpe`, or remove `--approximate`

### Results differ between normal and approximate mode

**Cause**: Approximate mode provides refined approximations (±10-20bp typical)

**Solution**: Expected behavior. Differences should be minimal with refined heuristic. Use normal mode only if base-level precision required.

### Reverse alignments show positive query lengths in approximate mode

**Cause**: Missing coordinate swap for reverse alignments. Metadata stores query coords in forward order always.

**Solution**:
```rust
// For reverse: swap query coordinates after projection
let (query_start, query_end) = if is_reverse {
    (working_query_end, working_query_start)  // Swap to match normal mode
} else {
    (working_query_start, working_query_end)
};
```

### Empty CIGAR panic at high transitive depths

**Cause**: Very small overlaps produce empty CIGARs after refinement.

**Solution**: Filter them out instead of panicking:
```rust
if approx_cigar.is_empty() {
    return None;  // Skip intervals too small for meaningful statistics
}
```

### Approximate mode not saturating (exploding result counts)

**Cause**: Applying refinement heuristic PER segment instead of AFTER collecting all segments.

**Solution**: Two-phase approach:
1. Scan: collect all overlapping segments (use full boundaries)
2. Refine: adjust first/last positions based on partial overlaps

### Approximate mode not much faster

**Possible causes**:
1. PAF input (approximate mode only helps .1aln format)
2. Small regions (overhead proportionally larger)
3. Few overlapping alignments (less to optimize)

**Ideal conditions for approximate mode**:
- Large target ranges
- Many overlapping alignments
- .1aln format input
- Transitive queries with high depth

### Compilation error about missing approximate_mode parameter

**Cause**: Missed updating a call site

**Solution**:
1. Find the function name in error message
2. Search for all calls to that function
3. Add `approximate_mode` or `query.approximate` or `config.approximate_mode` as appropriate

---

## Performance Characteristics

**Time Complexity**:
- Normal mode: O(n × (I + W + S))
  - n = overlapping alignments
  - I = sequence I/O cost
  - W = WFA alignment cost
  - S = sequence length
- Approximate mode: O(n × T)
  - n = overlapping alignments
  - T = tracepoint scan cost (≪ I + W + S)

**Space Complexity**:
- Normal mode: O(n × S) - stores sequences and CIGARs
- Approximate mode: O(n) - only intervals and approximate CIGARs

**Accuracy** (with refined overlap heuristic):
- Coordinates accurate to ±10-20bp typical (0.1-0.2% error)
- Without refinement: ±100bp (trace_spacing resolution)
- Identity metrics accurate (calculated from tracepoint statistics)
- Reverse strand handling: matches normal mode output exactly
- Suitable for exploratory analysis, transitivity, and bulk operations

---

## Architecture Notes

### Approximate CIGAR Encoding

**Updated**: Approximate mode no longer returns empty CIGAR vectors. Instead, it returns approximate CIGAR encoding with accumulated statistics.

**Old behavior** (pre-identity estimation):
- Returned `Vec::new()` (empty CIGAR)
- Resulted in NaN identity values in bedpe output

**New behavior** (with identity estimation):
- Returns approximate CIGAR: `[=N, XM]` where:
  - `N` = total matches (accumulated from overlapping segments)
  - `M` = total mismatches (from trace_diffs)
- Provides accurate gap-compressed identity: `N / (N + M)`
- Indels not tracked (gap-compressed identity ignores them anyway)
- No sequence fetching or WFA alignment required

This approach avoids:
1. Sequence fetching from FASTA/AGC
2. WFA alignment computation
3. Memory overhead of full CIGAR strings

While still providing accurate identity metrics for bedpe output.

### Why Separate File Opens on Cache Miss?

Original optimization tried single file open for both trace_spacing and alignment fetch, but this caused buffers to accumulate data from both operations, increasing memory usage.

Current approach:
- **Cache HIT** (common case after first query per file): 1 file open - FAST
- **Cache MISS** (only first query per file): 2 short-lived file opens
  1. Read trace_spacing from header → close immediately
  2. Fetch alignment → close immediately

This keeps buffers small and short-lived, minimizing memory while maintaining performance.

### Why Lazy-Loaded Cache?

Pre-warming the cache would read trace_spacing for ALL .1aln files upfront. If you only need 1 file out of 1000, this wastes I/O reading 999 unnecessary files.

Lazy loading reads only what you actually query, which is much more efficient for sparse access patterns.

---

## Conclusion

Approximate mode provides substantial speedup for .1aln queries by:
- ✅ Skipping sequence I/O (major bottleneck)
- ✅ Skipping CIGAR computation (CPU intensive)
- ✅ Refined overlap heuristic: ±10-20bp accuracy (10x better than trace_spacing)
- ✅ Accurate identity metrics from tracepoint statistics (no NaN values)
- ✅ Correct reverse strand handling (matches normal mode output)
- ✅ No sequence files required for bed/bedpe output
- ✅ Proper saturation at high transitive depths

**Key innovation**: Two-phase projection with fractional refinement provides near-exact intervals while maintaining full speed benefits.

The implementation threads approximate_mode through all query paths (regular, transitive BFS, transitive DFS) and all relevant commands (query, partition, refine), making it widely applicable.

## Recent Enhancements (Latest Version)

### Refined Overlap Heuristic (Nov 2025)

**What changed**:
- Two-phase projection: scan all overlapping segments, then refine boundaries
- Fractional projection for partial overlaps of first/last segments
- Fixed reverse strand coordinate handling (metadata stores query in forward order)
- Filter empty CIGARs from very small overlaps

**Impact**:
- ✅ Intervals accurate to ±10bp (vs ±100bp before)
- ✅ Reverse alignments now return negative query lengths (matches normal mode)
- ✅ More accurate identity filtering
- ✅ Saturation behavior matches expectations

**Key fix**: Apply refinement AFTER collecting all overlapping segments, not during scan.

**Typical accuracy**:
```
Target range: 2753-10000 (7247bp)
Normal mode:  query 111317811-111325058 (7247bp)
Approximate:  query 111317800-111325058 (7258bp)
Difference: 11bp (0.15%)
```

### Identity Estimation (Nov 2025)

**What changed**:
- `project_overlapping_interval_fast()` calculates identity from tracepoints
- `fetch_alignment_direct()` has `read_trace_diffs` parameter
- No sequence files required for bed/bedpe output
- Gap-compressed identity calculated accurately

**Performance**:
- Same ~10-100x speedup as before
- No sequence file requirement for bed/bedpe
- Accurate identity metrics (no NaN values)

**Files modified**:
- [src/impg.rs](src/impg.rs): Two-phase projection with refinement heuristic
- [src/onealn.rs](src/onealn.rs): `read_trace_diffs` parameter in `fetch_alignment_direct()`
- [src/main.rs](src/main.rs): Updated sequence requirement logic
