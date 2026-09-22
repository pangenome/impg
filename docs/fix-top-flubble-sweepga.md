# Fix top-flubble SweepGA graph induction bug

## Summary

The broken top-level SweepGA C4 run was not just a weak alignment result. It
had two concrete pipeline failures that allowed local replacements to become
unfolded path concatenations:

1. Nonempty, coordinate-consistent wfmash PAF could reach seqwish with a
   `min_match_len` larger than every exact CIGAR run in the block. Seqwish then
   consumed no usable matches and emitted one isolated path-specific segment per
   traversal.
2. Top-level blocks with zero emitted PAF records were still accepted as valid
   replacement graphs. Those replacements preserved path spellings but carried
   the original traversals forward as unfolded local graphs.

The fix makes seqwish's effective minimum match length respect the observed
exact-match runs in the filtered PAF, records per-block debug summaries, fails
when nonempty PAF produces no shared replacement nodes, and skips zero-PAF
top-level blocks instead of integrating an unfolded replacement.

## Reproducer Block

Small real C4-derived block:

```text
Before-fix artifact:
/home/erikg/impg/data/fix_top_flubble_repro_20260527T182100Z

After-fix artifact:
/home/erikg/impg/data/fix_top_flubble_repro_fixed2_20260527T184000Z

Committed fixture FASTA:
tests/test_data/crush/top_flubble_seqwish_minrun.fa

Committed fixture PAF:
tests/test_data/crush/top_flubble_seqwish_minrun.paf
```

The failing block had 5 real haplotype traversal records totaling 1922 bp and
16 raw/filtered wfmash PAF records. The PAF sequence names and lengths matched
the FASTA records. The broken run passed those alignments to seqwish with an
effective `min_match_len` of 299 even though the maximum exact CIGAR run was
218. Seqwish therefore emitted:

```text
segments: 5
segment bp: 1922
links: 0
paths: 5
shared replacement segments: 0
```

After the fix, the same block clamps the effective seqwish threshold to the
observed exact run:

```text
max_exact_match_run: 218
effective_min_match_len: 218
seqwish_segments: 6
seqwish_segment_bp: 1704
seqwish_links: 2
seqwish_paths: 5
shared replacement segments: 1
```

That is the first exact point where aligned traversals stopped compressing:
PAF survived filtering and names were consistent, but seqwish's match threshold
made those alignments unusable for graph induction.

## Code Changes

- `src/syng_graph.rs`
  - Adds filtered-PAF CIGAR inspection and exact-match run statistics.
  - Lowers the seqwish `min_match_len` only when filtered PAF records exist and
    the best observed exact run is below the configured threshold.
  - Writes `summary.tsv` under `IMPG_CRUSH_DEBUG_DIR/graph_build_*` with input,
    PAF, name/length consistency, seqwish, and effective threshold metrics.
- `src/resolution.rs`
  - Adds replacement shared-segment counting.
  - Records `replacement_shared_segments` in top-flubble evidence logs.
  - Fails loudly when nonempty PAF yields a replacement graph with no shared
    segment used by more than one path.
  - Returns a distinct zero-PAF replacement error and top-level mode treats it
    as a skipped local block. The original graph interval is left intact.
- `tests/test_crush_integration.rs`
  - Adds `c4_top_flubble_seqwish_indexes_observed_exact_runs`, a focused
    regression using the committed C4 FASTA/PAF fixture. It requires seqwish to
    induce shared nodes and reduce replacement segment bp below input bp.

## Real C4 Validation

Run directory:

```text
/home/erikg/impg/data/c4_top_flubble_fix_skipzero_20260527T185505Z
```

Command:

```bash
TMPDIR=/home/erikg/impg/data/c4_top_flubble_fix_skipzero_20260527T185505Z/tmp \
RAYON_NUM_THREADS=8 \
IMPG_CRUSH_DEBUG_DIR=/home/erikg/impg/data/c4_top_flubble_fix_skipzero_20260527T185505Z/debug \
RUST_LOG=impg=info \
target/release/impg query \
  -t 32 \
  --temp-dir /home/erikg/impg/data/c4_top_flubble_fix_skipzero_20260527T185505Z/tmp \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r 'GRCh38#0#chr6:31891045-32123783' \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=top-flubble-sweepga,aligner=wfmash,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=0:nosort' \
  -O /home/erikg/impg/data/c4_top_flubble_fix_skipzero_20260527T185505Z/run.nosort \
  -v 1
```

The command exited 0.

## Per-block Evidence

The real C4 run wrote 25 aligned block debug summaries. Each aligned block has:

```text
combined.fa
raw.paf
filtered.paf
seqwish.gfa
summary.tsv
```

Example:

```text
/home/erikg/impg/data/c4_top_flubble_fix_skipzero_20260527T185505Z/debug/graph_build_0000/summary.tsv

sequences: 48
total_sequence_bp: 10848
raw_paf_records: 2256
filtered_paf_records: 2256
parsed_filtered_paf_records: 2256
filtered_paf_records_with_cigar: 2256
paf_name_mismatches: 0
paf_length_mismatches: 0
max_exact_match_run: 226
effective_min_match_len: 226
seqwish_segments: 2
seqwish_segment_bp: 452
seqwish_paths: 48
```

The corresponding top-flubble evidence log proves the laced replacement used
the induced shared graph:

```text
region 11/484 evidence:
input_traversals=48
input_bp=10848
emitted_alignment_records=2256
aligned_bp=509856
paf_bytes=333309
replacement_segments=2
replacement_shared_segments=2
replacement_bp=452
applied=true
```

The final C4 run applied 25 aligned replacements and skipped 459 zero-PAF
blocks:

```text
crush: 25 resolved, 459 bailed, 484 candidates seen across 1 rounds
```

Skipping is deliberate here. A zero-PAF block does not contain alignment
evidence for seqwish to consume, so building a replacement graph from it would
only reintroduce unfolded local traversals.

## Before/After C4 Metrics

Path preservation against the no-crush reference:

```text
Reference paths: 465
Output paths: 465
Missing path names: 0
Extra path names: 0
Path spelling mismatches: 0
```

Graph metrics:

| Run | Applied replacements | Skipped zero-PAF blocks | Segments | Segment bp | Trivial-stringy |
| --- | ---: | ---: | ---: | ---: | ---: |
| Broken top-flubble SweepGA | 484 | 0 | 221342 | 26037441 | 398 |
| Fixed top-flubble SweepGA | 25 | 459 | 20118 | 823380 | 468 |

The fixed graph is still not a competitive crush result, but it is no longer
the 26 Mbp unfolded-concatenation failure. The segment-bp drops by about 31x
and segment count by about 11x relative to the broken top-flubble result while
preserving all 465 path spellings.

## PNG

```text
Local:
/home/erikg/impg/data/c4_top_flubble_fix_skipzero_20260527T185505Z/c4-top-flubble-fix.png

Remote:
www/impg/c4-top-flubble-fix.png

Confirmed:
ssh erik@hypervolu.me 'ls -lh www/impg/c4-top-flubble-fix.png'
-rw-r--r-- 1 erik erik 1.1M May 27 19:07 www/impg/c4-top-flubble-fix.png
```

## Validation Commands

Focused regression:

```bash
cargo test --test test_crush_integration c4_top_flubble_seqwish_indexes_observed_exact_runs -- --nocapture
```

Full test suite:

```bash
cargo test --all
```

Both were run with the native dependency environment from
`/home/erikg/impg/env.sh`.
