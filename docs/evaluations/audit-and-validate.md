# Audit and Validate: C4 SPOA Replacement Lacing

Task: `audit-and-validate`
Date: 2026-06-04

## Verdict

The first hypothesis was only partly right. SPOA replacement graphs are built
and laced into the graph; the failure is not that the aligner never runs. The
full C4 reruns show that direct SPOA is resolving hundreds of small flubbles,
but those replacements do not improve the visible graph and substantially worsen
path-order/white-space metrics.

A global `sequence -> segment id` projection for replacement nodes was tested
and rejected as a default. It reduced segment count on small slices, but on full
C4 it over-collapsed short identical replacement segments across distant local
contexts and made long white-space bridges worse. That patch was removed from
the working tree.

The retained code change is narrower: replacement nodes that occur before the
first surviving original node on a path are inserted before that original node
instead of being appended after all originals. This is a real renderer/lacing
corner-case fix, but it does not change the current C4 SPOA result. The C4
problem remains algorithmic: direct SPOA is polishing tiny local flubbles rather
than resolving the larger underaligned structures visible in the SweepGA seed
render.

Hard gate status: exact path spelling is preserved. No quality or metric gate
was added.

## Evidence

Baseline seed:

```text
segments=7411 links=10072 bp=234828 bp_weighted_cov=454.330293
singleton_bp=2922 path_jump_p99=5 white_p99=66 bridges>=1kb=8790
duplicate_sequence_frac=0.669950
```

Original true SPOA median100:

```text
segments=8119 links=10904 bp=239264 bp_weighted_cov=445.906923
singleton_bp=3895 path_jump_p99=3454 white_p99=111591 bridges>=1kb=212876
duplicate_sequence_frac=0.714127
```

Anchor-only patch rerun:

```text
segments=8119 links=10904 bp=239264 bp_weighted_cov=445.906923
singleton_bp=3895 path_jump_p99=3454 white_p99=111591 bridges>=1kb=212876
duplicate_sequence_frac=0.714127
```

The anchor-only rerun is metric-identical to the original true SPOA median100
output. Therefore the path-start insertion fix is correct but not the limiting
C4 issue.

The full C4 median100 run confirms the direct SPOA budget is applied correctly
when the binary is rebuilt:

```text
crush discovery detail: skipped 1112 explicit SPOA candidate(s) outside direct length budgets
crush round 1 traversal stats: selected n=353, max-len median/max=61/100,
median-len median/max=61/100
```

Across the full run, direct SPOA resolved 523 candidates across three rounds and
preserved all path spellings, but the result was worse than the seed.

## Validation

Validated commands:

```bash
cargo test replacement_at_path_start_is_ordered_before_first_surviving_original --lib
cargo test direct_poa_respects_median_traversal_budget --lib
cargo build
cargo test --lib
target/debug/impg crush --gfa data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_scaffold0.initial.gfa \
  --output data/c4_spoa_anchor_probe_20260604T090500Z/graphs/spoa_median100.anchor.gfa \
  --method poa --max-iterations 20 --max-span 0 --max-traversal-len 50k \
  --max-median-traversal-len 100 --max-total-sequence 2g --max-traversals 100k \
  --poa-scoring 1,4,6,2,26,1 --threads 32 -v 1
target/debug/impg graph-report \
  -g data/c4_spoa_anchor_probe_20260604T090500Z/graphs/spoa_median100.anchor.gfa \
  -o data/c4_spoa_anchor_probe_20260604T090500Z/reports/spoa_median100.anchor.graph-report.tsv \
  --format tsv
```

The full library test suite passed with `376` tests and `0` failures.

## Recommendation

Do not use direct SPOA as a polishing pass over the current C4 SweepGA seed.
Also do not add global sequence-deduplication as a default renderer behavior.

Next useful work should target the actual missing operation: build larger,
context-aware replacement blocks from the residual flubble/bubble structure and
induce those blocks with SweepGA/seqwish or another many-to-many aligner before
optionally applying bounded POA to small residual tangles.
