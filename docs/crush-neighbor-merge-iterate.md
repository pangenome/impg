# Crush experiment: neighbor-merge POASTA iteration

**Task:** `crush-neighbor-merge`
**Date:** 2026-05-27
**Output dir:** `/home/erikg/impg/data/c4_crush_neighbor_merge_iterate_20260527T124218Z/`
**PNG (local):** `/home/erikg/impg/data/c4_crush_neighbor_merge_iterate_20260527T124218Z/c4-crush-neighbor-merge-iterate.png`
**PNG (uploaded):** `https://hypervolu.me/~erik/impg/c4-crush-neighbor-merge-iterate.png`

## What Changed

This pass adds a local neighborhood smoother for post-crush graphs:

| Stage syntax | Behavior |
|---|---|
| `:neighbor-merge-poasta` | Decompose all POVU bubble sites on the current graph, select non-overlapping local sites, greedily merge path-adjacent sites into reference-span groups under a configurable cap, run POASTA per local neighborhood, lace all blocks back together, then repeat for the configured number of iterations. |
| `:neighbor-merge`, `:bubble-neighbor-merge` | Aliases for the same stage. |

The stage defaults to three 10 kb iterations. The cap and count can be set with:

```bash
gfa:syng:mask,min-run=3:crush,...:neighbor-merge-poasta,iterations=3,target-length=10k:nosort
```

Implementation notes:

- `src/main.rs` parses the new stage aliases and parameters.
- `src/smooth.rs` adds `SmoothBlockSource::NeighborMergePoasta`, POVU site selection, path-adjacent greedy grouping, per-iteration metrics, and per-block POASTA dispatch.
- `src/resolution.rs` exposes `poasta_sequences_to_gfa`, a path-preserving POASTA helper reused by the smoother.
- `examples/neighbor_merge_existing_gfa.rs` is a validation helper for applying the same smoother to an existing real GFA, followed by gfaffix and `Ygs` sort.

## Tests

Targeted tests:

| Test | Coverage |
|---|---|
| `smooth::tests::test_neighbor_site_groups_merge_adjacent_sites_up_to_cap` | Greedy adjacent-site grouping respects the cap. |
| `smooth::tests::test_neighbor_merge_blocks_select_povu_sites_on_nested_fixture` | POVU sites from a nested fixture become neighbor-merge blocks plus passthrough gaps. |
| `smooth::tests::test_neighbor_merge_poasta_preserves_path_sequences_on_nested_fixture` | Local POASTA smoothing preserves path-spelled sequences. |
| `main::tests::test_gfa_engine_neighbor_merge_poasta_stage_parses_defaults` | Stage defaults to `NeighborMergePoasta` with three 10 kb passes. |
| `main::tests::test_gfa_engine_neighbor_merge_poasta_stage_parses_iterations_and_cap` | `iterations=` and `target-length=` parse through the real CLI path. |

Commands run:

```bash
cargo test --lib neighbor_ -- --nocapture
cargo test --bin impg test_gfa_engine_neighbor_merge -- --nocapture
cargo build --release --lib --bin impg
cargo build --release --bin gfaffix
cargo build --release --example neighbor_merge_existing_gfa
```

The local build required the same validation workaround used by earlier C4 tasks: `wfmash-rs` tries to build vendored wfmash, but the vendored C++ source is missing `<limits>` in `rkmh.cpp`. I copied the existing `/home/erikg/bin/wfmash` binary into Cargo's `wfmash-rs` `OUT_DIR` before compiling, so repository files were not changed.

## C4 Validation

Starting point was the allowed iterative-multi graph:

```text
/home/erikg/impg/data/c4_crush_iterative_multi_20260527T050056Z/run.Ygs.gfa
```

Validation command:

```bash
out=/home/erikg/impg/data/c4_crush_neighbor_merge_iterate_20260527T124218Z
input=/home/erikg/impg/data/c4_crush_iterative_multi_20260527T050056Z/run.Ygs.gfa
RUST_LOG=info /usr/bin/time -v -o "$out/time.txt" \
  target/release/examples/neighbor_merge_existing_gfa \
  "$input" "$out/run" 3 10000 32 \
  > "$out/stdout.txt" 2> "$out/stderr.txt"
target/release/gfaffix "$out/run.raw-smooth.gfa" \
  -o "$out/run.gfaffix.gfa" -p 32
gfasort -i "$out/run.gfaffix.gfa" \
  -o "$out/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/run.Ygs.gfa" \
  -o "$out/c4-crush-neighbor-merge-iterate.png" -m -x 2200 -y 1200
scp "$out/c4-crush-neighbor-merge-iterate.png" \
  erik@hypervolu.me:www/impg/c4-crush-neighbor-merge-iterate.png
ssh erik@hypervolu.me 'ls -lh www/impg/c4-crush-neighbor-merge-iterate.png'
```

Upload confirmation:

```text
-rw-r--r-- 1 erik erik 1.4M May 27 13:10 www/impg/c4-crush-neighbor-merge-iterate.png
```

## Per-Iteration Metrics

From `stderr.txt`:

| Iteration | Cap | POVU sites | Candidate sites | Selected sites | Selected groups | Neighbor blocks | Gap blocks | Blocks merged | Sites in merged blocks | Overlap skipped | Over cap skipped | Segs before | Segs after | Segment-bp before | Segment-bp after | Wall |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 10,000 | 4,590 | 4,320 | 3,925 | 23 | 16 | 7 | 23 | 3,925 | 1,682 | 5 | 16,547 | 20,176 | 396,438 | 444,983 | 1,136.754 s |
| 2 | 10,000 | 4,309 | 4,276 | 3,895 | 23 | 16 | 7 | 23 | 3,895 | 1,662 | 12 | 20,176 | 19,022 | 444,983 | 462,907 | 171.902 s |
| 3 | 10,000 | 4,308 | 4,275 | 3,895 | 23 | 16 | 7 | 23 | 3,895 | 1,661 | 12 | 19,022 | 18,780 | 462,907 | 462,860 | 126.505 s |

Post-processing:

| Step | Wall | Result |
|---|---:|---|
| gfaffix | 0.43 s | 18,780 S / 25,787 L / 462,860 bp -> 18,761 S / 25,763 L / 461,241 bp |
| `gfasort -p Ygs` | 32.17 s | 18,761 S / 25,763 L / 461,241 bp |

## Path Preservation

Path count and byte-spelled path sequence preservation:

```text
before=465 after=465 missing=0 extra=0 mismatched=0
```

The comparison keyed paths by the PanSN prefix before `:` because lacing rewrites coordinate suffixes.

## Final Graph Metrics

| Graph | Segments | Links | Paths | Segment bp | Trivial-stringy |
|---|---:|---:|---:|---:|---:|
| Iterative-multi input | 16,547 | 23,610 | 465 | 396,438 | 103 |
| Neighbor-merge raw smooth | 18,780 | 25,787 | 465 | 462,860 | not rerun |
| Neighbor-merge gfaffix + Ygs | 18,761 | 25,763 | 465 | 461,241 | 72 |
| PGGB control target | 13,288 | 16,240 | 465 | 234,524 | 12 |

The result does not meet the victory threshold (`segments <= 14,500` or `stringy <= 30`). It does reduce the stringy heuristic versus the iterative-multi input, but it regresses segment count and segment-bp.

## Dup-Sequence Extras

Computed with `/tmp/gfa_seqdup.py`, which canonicalizes forward/reverse-complement segment sequences and reports extras as `segments_in_duplicate_group - 1`.

| Graph | <=4 bp | 5-10 bp | 11-50 bp | 51-200 bp | >200 bp |
|---|---:|---:|---:|---:|---:|
| Iterative-multi input | 11,415 | 282 | 83 | 6 | 1 |
| Neighbor-merge gfaffix + Ygs | 11,183 | 315 | 189 | 68 | 0 |

The high-value long duplicate band (`>200 bp`) is eliminated, but the 51-200 bp band regresses from 6 to 68 extras.

## Interpretation

This implementation confirms that local neighborhood merging is feasible and path-preserving on real C4, and it can reduce the trivial-stringy heuristic from 103 to 72. The current grouping is too aggressive or too coarse for compactness: first iteration expands the graph to 20,176 segments, and later iterations only partially recover to 18,761.

The largest slow block in iteration 1 was:

```text
block 7 sites=234 ref_steps=3904-4463 ref_span_bp=9907 ranges=330 input_bp=3345540 output_bp=29497
```

That one block dominated the first pass wall time. A likely next refinement is to keep path-adjacent grouping but add a secondary cap on observed block sequence mass or max traversal length, so a 10 kb reference-span group with high-depth divergent insertions is split before POASTA.
