# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`impg` (implicit pangenome graph) is a Rust tool for rapid querying of pangenome alignments. It treats all-vs-all pairwise alignments as an implicit pangenome graph, enabling fast projection of target ranges through alignment networks to extract homologous sequences without building explicit graph structures.

## Build and Test Commands

```bash
# Install (recommended for production use)
cargo install --force --path .

# Build for development
cargo build                    # Debug build
cargo build --release          # Release build (use for benchmarking)
cargo build --release --no-default-features  # Without AGC support

# Run without installing
cargo run --release -- query -a alignments.paf -r chr1:1000-2000

# Testing
cargo test                     # All tests
cargo test test_name           # Specific test
cargo test -- --nocapture      # With output

# Code quality
cargo check                    # Check compilation
cargo clippy                   # Linter
cargo fmt                      # Format code
```

## Core Architecture

### Alignment Index System
- **coitrees** (Cache Oblivious Interval Trees): Fast range lookups over alignments
- **Serialized indices**: Bincode-serialized indices cached on disk for fast startup
- **Two input formats**:
  - **PAF files**: Standard pairwise alignment format with CIGAR (from wfmash/minimap2 --eqx)
  - **.1aln files**: Compressed format using tracepoints instead of full CIGAR strings

### Key Modules

- `src/main.rs`: CLI argument parsing, command dispatch, output format handling
- `src/impg.rs`: Core interval projection logic, CIGAR delta encoding, transitive queries (BFS/DFS)
- `src/alignment_record.rs`: Unified representation for PAF and .1aln alignments
- `src/onealn.rs`: Parser for .1aln format with tracepoint-to-CIGAR reconstruction
- `src/paf.rs`: PAF file parsing and CIGAR handling
- `src/sequence_index.rs`: Unified sequence access for FASTA and AGC archives
- `src/forest_map.rs`: Interval tree forest mapping sequences to their alignment trees
- `src/commands/`: Subcommands (partition, refine, similarity, lace)

### .1aln Format Details

The .1aln format stores alignments compactly using tracepoints:
- `trace_spacing`: Regular interval on query sequence (typically 100bp)
- `tracepoints[]`: Array of consumed lengths on target sequence at each interval
- `trace_diffs[]`: Number of differences at each tracepoint (for identity calculation)

**Two execution modes for .1aln**:
1. **Normal mode**: Reconstructs full CIGAR from tracepoints using WFA (Wavefront Alignment)
2. **Approximate mode** (`--approximate`): Skips CIGAR reconstruction, uses refined tracepoint statistics for ~10-100x speedup

See `notes/FAST_MODE_IMPLEMENTATION.md` for detailed implementation notes on tracepoint subsetting and approximate mode.

### Transitive Queries

The tool supports finding alignments connected through multiple hops:
- **BFS mode** (default `-x`): Breadth-first exploration, finds more overlapping results
- **DFS mode** (`--transitive-dfs`): Depth-first exploration, fewer overlaps but deeper paths
- **Max depth** (`-m`): Control transitive search depth (0 = unlimited)

## Commands

- **query**: Query overlapping alignments for a region or BED file
- **partition**: Partition alignments into windows across sequences
- **refine**: Refine loci to maximize sample support
- **similarity**: Compute pairwise similarity matrices with optional PCA
- **lace**: Combine multiple GFA or VCF files
- **stats**: Print alignment statistics
- **index**: Create IMPG index from alignment files

## Development Notes

### Adding New Output Format
1. Update `OutputFormat` enum in `src/main.rs`
2. Implement format writer in `write_*_output()` function
3. Update format validation in command handlers

### Modifying Alignment Processing
1. Core projection logic: `src/impg.rs::project_overlapping_interval()`
2. For .1aln approximate mode: `project_overlapping_interval_fast()`
3. CIGAR operations use `CigarOp` enum with delta encoding

### Working with Sequences
- `UnifiedSequenceIndex` handles both FASTA and AGC formats
- AGC support is optional (controlled by `agc` feature flag)
- Use `#[cfg(feature = "agc")]` for AGC-specific code

## Important Implementation Notes

### CIGAR Delta Encoding
`impg` uses compact delta encoding for CIGAR operations rather than storing full strings. When working with CIGAR:
- Operations are `Vec<CigarOp>` where each op is an enum variant (Match, Insertion, Deletion, etc.)
- Convert to/from string representation as needed for output formats

### Coordinate Systems
- **PAF coordinates**: 0-based, half-open intervals [start, end)
- **BED coordinates**: 0-based, half-open intervals [start, end)
- **GFA coordinates**: 1-based, closed intervals [start, end]
- Target reverse strand: coordinates are flipped within the contig
- Query reverse strand: coordinates stored forward in metadata, swap after projection for output

### Memory and Performance
- Use `rayon` for parallelization (configured via `-t` threads parameter)
- Interval trees are memory-intensive for large alignments
- Index caching avoids re-parsing alignment files
- For .1aln approximate mode: skip sequence fetching for bed/bedpe output

### Approximate Mode Requirements
When using `--approximate` with transitive queries (`-x` or partition/refine commands):
- Must set `--min-transitive-len` > trace_spacing
- Typical: `--min-transitive-len 101` for trace_spacing=100bp
- Non-transitive queries don't need this restriction

## Dependencies of Note

- `rust-htslib`: Fixed at version 0.46.0 (see Cargo.toml comment about issue #434)
- `lib_wfa2`, `lib_tracepoints`, `onecode`: For .1aln format support
- `handlegraph`: For GFA graph operations in lace command

## Testing

Test data is in `tests/test_data/` with sample FASTA and AGC files. Key testing considerations:
- Test with both PAF and .1aln alignment formats
- Verify output formats: bed, bedpe, paf, gfa, maf, fasta, fasta-aln
- Test transitive queries at various depths (-m 1, -m 2, -m 10)
- Compare approximate mode results to normal mode for accuracy
- Test edge cases: reverse strand alignments, overlapping intervals

## Debugging

- Verbosity: `-v 2` for debug logs, or `RUST_LOG=debug`
- For .1aln issues: check trace_spacing cache and tracepoint reconstruction
- Performance profiling: always use `--release` build
