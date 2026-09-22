# C4 fragment regression suite

The C4 crush regression suite now uses small, committed fragments in
`tests/test_data/crush/c4_fragments/` plus the existing
`top_flubble_seqwish_minrun` reproducer. The goal is to keep the failure modes
from the full C4 top-flubble run under `cargo test` without rerunning the full
locus.

## Source

The new fixtures were extracted from the corrected top-flubble SweepGA debug
run:

```text
/home/erikg/impg/data/c4_top_flubble_fix_skipzero_20260527T185505Z/debug
```

For each selected debug block, the fixture generator kept a small set of FASTA
records and retained only raw PAF lines whose query and target names were still
present. This preserves the real C4 traversal spellings, sequence IDs, PAF
coordinates, and PAF lengths while making each fixture only a few kilobytes.

## Covered Fragments

| Test name | Fixture | Source block | Protected behavior |
| --- | --- | --- | --- |
| `easy_shared_flank` | `c4_fragments/easy_shared_flank.{fa,paf,gfa}` | `graph_build_0000` | Simple shared-flank/top-flubble case with clear allelic sharing. The lacing fixture must replace duplicate traversal interiors with induced shared graph nodes. |
| `bounded_multi_bubble` | `c4_fragments/bounded_multi_bubble.{fa,paf}` | `graph_build_0001` | Locally complex but bounded region with three compact allelic classes. The induced graph must stay small and shared. |
| `unfolded_minrun` | `top_flubble_seqwish_minrun.{fa,paf}` | prior top-flubble reproducer block | Previously broken class where nonempty PAF was present but seqwish emitted isolated path-specific segments because the effective min-match floor was too high. |
| `short_floor` | `c4_fragments/short_floor.{fa,paf}` | `graph_build_0002` | Short 159-300 bp traversals around the 311 bp filtering/seqwish floor. The PAF filter and seqwish min-match rescue must leave usable alignments. |
| `duplicated_repeat` | `c4_fragments/duplicated_repeat.{fa,paf}` | `graph_build_0004` | Repeated/duplicated C4 sequence where alignments must survive filtering and seqwish must consume them into shared nodes rather than duplicate unfolded segments. |

## Test Invariants

`c4_fragment_seqwish_regressions_induce_shared_graphs` feeds each committed
FASTA/PAF pair through the real `build_gfa_from_paf_and_sequences` pipeline. It
asserts:

- PAF records are nonempty for homologous fragments.
- PAF query/target names and lengths match the FASTA records handed to seqwish.
- Output path names and spellings equal the input FASTA records.
- The induced graph has shared nodes with fixture-specific depth.
- Segment count and segment bp remain under fixture-specific bounds.
- Segment bp is smaller than total input traversal bp, catching unfolded
  path-concatenation output.
- Duplicate segment-sequence signatures remain under fixture-specific bounds.

`c4_fragment_lacing_uses_pairwise_induced_graphs` wraps each same traversal set
in a tiny shared-flank GFA and runs `resolve_gfa_bubbles` with
`method=TopFlubbleSweepga` and the FastGA backend on a bounded thread pool. It
asserts:

- Path spellings are preserved after replacement/lacing.
- Exactly one top-level replacement is applied.
- No bailout/quality-guard fallback hides a bad replacement.
- The laced graph contains shared path-supported nodes and reduces segment bp
  relative to the unfolded input bubble.

These checks specifically would have caught the broken top-flubble
unfolded/concatenated output class: the old behavior had nonempty PAF but zero
shared replacement nodes and segment bp equal to the sum of all traversal
spellings.
