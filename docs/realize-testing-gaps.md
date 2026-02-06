# Realize Engine: Test Coverage Evaluation

## Summary

The realize engine has **28 unit tests** across two modules (`realize/mod.rs`: 11 tests, `realize/poa.rs`: 17 tests). All 61 library tests pass. Coverage is strongest for the POA subsystem and pure-function helpers, but significant gaps exist in the recursive orchestrator, lacing logic, PAF parsing, and integration-level behavior.

---

## Current Test Inventory

### `realize/mod.rs` (11 tests)

| Test | What It Covers |
|------|---------------|
| `test_make_trivial_gfa` | Single-sequence GFA generation |
| `test_realize_from_sequences_empty` | Empty input → header-only GFA |
| `test_realize_from_sequences_single` | Single sequence → trivial GFA, stats |
| `test_realize_from_sequences_small_region_uses_poa` | Two short seqs → POA path (below threshold) |
| `test_project_window_to_sequence_no_overlap` | Window projection: no-overlap returns None |
| `test_project_window_to_sequence_full_overlap` | Window projection: exact alignment match |
| `test_project_window_to_sequence_partial_overlap` | Window projection: partial window coverage |
| `test_partition_single_chunk` | Short seqs → single chunk returned |
| `test_lace_subgraphs_single` | Lacing a single GFA is identity |
| `test_lace_subgraphs_two` | Lacing two GFAs: ID remapping, linking, path merge |
| `test_lace_subgraphs_empty` | Lacing empty input → header-only |

### `realize/poa.rs` (17 tests)

| Test | What It Covers |
|------|---------------|
| `test_find_core_column_range_no_padding` | Core range = full MSA when padding=0 |
| `test_find_core_column_range_symmetric_padding` | Symmetric padding trimming |
| `test_find_core_column_range_with_gaps` | Padding trimming with MSA gaps |
| `test_find_core_column_range_asymmetric_padding` | Different left/right padding per sequence |
| `test_find_core_column_range_empty_msa` | Empty MSA → (0,0) |
| `test_find_core_column_range_all_padding` | All-padding → degenerate midpoint range |
| `test_padded_poa_from_sequences_no_padding` | Two identical seqs, padding=0 → valid GFA |
| `test_padded_poa_from_sequences_with_padding` | Two seqs with 2bp padding → trimmed GFA |
| `test_padded_poa_from_sequences_with_variation` | SNP between two seqs → bubble graph |
| `test_padded_poa_empty_input` | Empty input → empty GFA |
| `test_padded_poa_single_sequence` | Single seq → trivial GFA |
| `test_padded_poa_reverse_strand_metadata` | Strand metadata preserved through POA |
| `test_padded_poa_large_padding_clamped` | Padding > seq_len/2 → no panic |
| `test_padded_poa_three_sequences_with_indel` | 3 seqs with indels → multi-segment graph |
| `test_post_process_gfa_for_strands_forward` | Forward strand path unchanged |
| `test_post_process_gfa_for_strands_reverse` | Reverse strand path: segments reversed, orientations flipped |
| `test_make_headers_forward_and_reverse` | Header generation for +/- strands |

---

## Gap Analysis

### CRITICAL: No tests for recursive path through sweepga

**Functions with zero test coverage:**
- `realize_recursive()` — the core orchestrator — is only tested via `realize_from_sequences()` at the POA base case. No test exercises the sweepga→partition→recurse→lace recursive path.
- `build_sweepga_config()` — never tested (not even indirectly).
- `parse_paf_file()` — never tested directly; only exercised if sweepga is invoked.

**Why this matters:** The recursive case is the main value proposition of the realize engine. A bug in partitioning, PAF parsing, or the sweepga→recurse→lace flow would be completely invisible to the current test suite.

**Proposed tests:**
1. **Mock-based recursive test**: Create a test that provides pre-made sequences large enough to exceed `poa_threshold`, with a mock or stub for `sweepga_align` that returns a known PAF file. This isolates the partition→recurse→lace logic from external dependencies.
2. **`parse_paf_file` unit tests**: Create temp PAF files with known content and verify parsed records. Include edge cases: empty file, malformed lines, lines with <12 fields, non-numeric fields.
3. **`build_sweepga_config` test**: Verify that RealizeConfig fields map correctly to SweepgaAlignConfig.

### HIGH: Partition logic under-tested

**`partition_into_chunks()`** has one test (`test_partition_single_chunk`) which only tests the trivial case where everything fits in one chunk.

**Missing scenarios:**
1. **Multi-chunk partitioning**: Sequences long enough to produce 2+ chunks. Verify:
   - Correct number of chunks created
   - Anchor slice coordinates are correct per chunk
   - Padding overlap between adjacent chunks exists
   - Non-anchor sequences appear in correct chunks based on projection
2. **Projection-based inclusion**: Non-anchor sequences that only partially span the anchor (should appear in some chunks but not others)
3. **Empty projections**: No PAF records for a non-anchor sequence → it should be absent from all chunks
4. **Edge case: chunk_size=0 or chunk_size=1**: Should not panic

**Proposed tests:**
```
test_partition_multiple_chunks
test_partition_non_anchor_partial_coverage
test_partition_no_projections_for_sequence
test_partition_small_chunk_size
```

### HIGH: Lacing correctness for complex graphs

**`lace_subgraphs()`** has 3 tests, but they only cover trivial cases (0, 1, or 2 single-segment GFAs). Missing:

1. **Multi-segment sub-GFAs**: Lacing two GFAs that each have multiple segments and links. Verify that:
   - All segments get unique remapped IDs
   - Intra-GFA links are correctly remapped
   - Cross-GFA linking edges connect correct nodes
2. **Multiple paths**: Sub-GFAs with 2+ paths. Verify each path is independently merged across chunks.
3. **Paths present in some chunks but not others**: A path that exists in GFA1 but not GFA2 should still appear in the output (with only GFA1's steps).
4. **Reverse-orientation steps in paths**: Steps like `3-` should survive ID remapping.
5. **Link deduplication**: Verify that duplicate links are correctly deduplicated.
6. **GFAs with optional fields**: Segments with 4+ tab-separated fields (e.g., tags) should be preserved.

**Proposed tests:**
```
test_lace_multi_segment_subgraphs
test_lace_multiple_paths
test_lace_partial_path_coverage
test_lace_reverse_oriented_steps
test_lace_link_deduplication
test_lace_segment_tags_preserved
```

### MEDIUM: `project_window_to_sequence()` edge cases

Has 3 tests (no-overlap, full, partial) but missing:

1. **Multiple overlapping projections**: Two alignment blocks that both overlap the window → should take union of projected ranges
2. **Window spans multiple alignment blocks with gap**: Alignment at [0,50) and [80,130), window [40,90) → should merge projections
3. **Zero-length alignment span**: `a_start == a_end` → should handle gracefully (division by zero risk in `a_span`)
4. **Very large coordinates**: Near `usize::MAX` values → test saturating arithmetic

### MEDIUM: `make_trivial_gfa()` reverse strand

The current test only covers forward strand. For reverse strand, the coordinate calculation `(total_length - start - size)` to `(total_length - start)` should be verified.

### MEDIUM: `realize()` entry point (with ImpgIndex)

`realize()` is never tested because it requires an `ImpgIndex` and `UnifiedSequenceIndex`. While `realize_from_sequences()` bypasses these, the `realize()` function contains the `prepare_sequences()` call + sort_gfa integration.

**Proposed approach:** Create a test with a small in-memory ImpgIndex (the crate has mock/test helpers for this based on existing impg.rs tests), call `realize()`, and verify GFA output.

### LOW: Stats tracking

`RealizeStats` fields (`poa_calls`, `sweepga_calls`, `max_depth_reached`) are tested for simple cases but not for:
- Multi-level recursion (verify `max_depth_reached` > 0)
- Mixed POA + sweepga paths (verify both counters increment)

### LOW: `sort_gfa` integration

`realize_from_sequences` and `realize` call `sort_gfa` when `config.sort_output = true` (default). All current tests set `sort_output = false`. No test verifies that GFA sorting works correctly end-to-end through the realize pipeline.

### LOW: Error paths

No tests for error conditions:
- `parse_paf_file` on non-existent file
- `sweepga_align` failure (e.g., temp dir not writable)
- GFA sorting failure

---

## Integration Test Coverage

### Existing integration tests

| Test File | Realize Coverage |
|-----------|-----------------|
| `test_pipeline_integration.rs` | Tests full pipeline (index→partition→graph→lace) but **not** the realize engine. Uses the older `graph` command, not `realize`. Requires wfmash + samtools (marked `#[ignore]`). |
| `test_transitive_integrity.rs` | Tests query/partition correctness. Does **not** exercise the realize engine at all. |
| `test_agc_integration.rs` | Tests AGC index integration. Not related to realize. |

**Gap:** There is no integration test that exercises the realize engine end-to-end. The test_pipeline_integration test uses the `graph` command (sweepga+seqwish), not the realize recursive engine.

**Proposed integration test:**
A test that:
1. Creates synthetic sequences with known variation
2. Calls `realize_from_sequences()` with config that will trigger recursion
3. Validates the output GFA has correct paths, segments, and links
4. Checks that all input sequences are represented as paths in the GFA
5. Optionally validates that the GFA MSA (via `gfa_to_msa`) contains the expected variation

---

## Priority-Ordered Recommendations

| Priority | Gap | Effort | Impact |
|----------|-----|--------|--------|
| **P0** | Recursive sweepga→partition→lace path untested | High (needs sweepga dependency or mock) | Critical — main engine logic |
| **P0** | `parse_paf_file()` unit tests | Low | Catches silent data corruption |
| **P1** | Multi-chunk `partition_into_chunks()` | Medium | Partitioning bugs = wrong graphs |
| **P1** | Complex `lace_subgraphs()` scenarios | Medium | Lacing bugs = broken GFA output |
| **P1** | Integration test for realize end-to-end | High | Confidence in real-world behavior |
| **P2** | `project_window_to_sequence()` edge cases | Low | Defensive testing |
| **P2** | Reverse-strand `make_trivial_gfa()` | Low | Simple addition |
| **P2** | `realize()` with ImpgIndex | Medium | Entry-point coverage |
| **P3** | Stats tracking in recursive cases | Low | Nice to have |
| **P3** | `sort_gfa` through realize pipeline | Low | Nice to have |
| **P3** | Error path testing | Low | Robustness |

---

## Approach Recommendations

### For P0 (recursive path testing)

The biggest challenge is that `sweepga_align` is an external dependency (calls FastGA). Options:

1. **Synthetic test with real sweepga**: Create sequences ~2-5kb that will trigger recursion. This is a true integration test but requires sweepga to be compiled (which it is, as a dependency).
2. **Extract and mock**: Factor out the alignment step behind a trait so tests can inject known PAF results. This is cleaner but requires refactoring.
3. **Test from PAF**: Write a helper that takes pre-computed PAF + sequences and runs just the partition→recurse→lace flow. This tests everything except the alignment step.

Recommendation: Option 3 first (isolate partition→lace), then option 1 for integration confidence.

### For P1 (partition + lace testing)

These are pure functions with no external dependencies — easy to unit test with synthetic data. The lacing function takes `&[String]` GFA strings and returns a merged GFA string. Easy to construct test cases.

### Testing infrastructure suggestions

- Add a `test_helpers` module in `realize/` with builders for common test fixtures (sequences, PAF records, GFA strings)
- Consider property-based testing (proptest) for `project_window_to_sequence` — the linear interpolation math has many corner cases
