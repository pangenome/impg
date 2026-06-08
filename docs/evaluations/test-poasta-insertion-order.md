# Test: POASTA Insertion Order On C4 Crush04

Date: 2026-06-05

## Scope

This is the focused insertion-order experiment for
`Poasta_crush04_top_span6460bp_med6461bp_cov47of465_sites1_steps3320-3523`.
The input was the applied replacement GFA from the prior audit artifact:

`/home/erikg/impg/.wg-worktrees/agent-478/data/audit_poasta_replacement_20260604T182011Z/debug/applied_frontier_0000/Poasta_crush04_top_span6460bp_med6461bp_cov47of465_sites1_steps3320-3523/replacement.gfa`

The prior audit identified this replacement at
`docs/evaluations/audit-poasta-replacement.md:91` and recorded the same compact
shape at `docs/evaluations/audit-poasta-replacement.md:102`: 72 segments, 96
links, 47 paths, and 1900 path steps. It also found no path corruption for the
matched POASTA debug replacement at
`docs/evaluations/audit-poasta-replacement.md:117`.

## Method

I added explicit-order replacement helpers in `src/resolution.rs:9217` and
`src/resolution.rs:9286`. These keep the output path order fixed to the input
headers while varying only the progressive insertion order. The baseline
`poasta_sequences_to_gfa` default remains longest sequence first, then path
name, at `src/resolution.rs:9196`.

The experiment driver is `examples/poasta_order_driver.rs`. It:

- reconstructs the 47 exact path sequences from the input replacement GFA;
- runs POASTA with current longest-then-name, reverse, medoid-first,
  nearest-neighbor guide-tree, and deterministic random orders;
- tries a reference/CHM13-first order when such a path exists;
- optionally runs abPOA with the same sequence set and insertion orders;
- writes per-order GFAs, path-labeled SVG renders, insertion orders, skipped
  orders, and diagnostic summary metrics.

Relevant driver entry points:

- order construction: `examples/poasta_order_driver.rs:76`
- longest-then-name: `examples/poasta_order_driver.rs:432`
- medoid selection: `examples/poasta_order_driver.rs:461`
- nearest-neighbor guide-tree order: `examples/poasta_order_driver.rs:475`
- deterministic random shuffle: `examples/poasta_order_driver.rs:534`
- exact path comparison: `examples/poasta_order_driver.rs:553`
- SCC/self-loop summary: `examples/poasta_order_driver.rs:582`
- SVG rendering through `gfalook`: `examples/poasta_order_driver.rs:729`
- TSV metric schema: `examples/poasta_order_driver.rs:776`

Command run:

```bash
cargo run --release --example poasta_order_driver -- \
  --input /home/erikg/impg/.wg-worktrees/agent-478/data/audit_poasta_replacement_20260604T182011Z/debug/applied_frontier_0000/Poasta_crush04_top_span6460bp_med6461bp_cov47of465_sites1_steps3320-3523/replacement.gfa \
  --output-dir data/test_poasta_insertion_order_20260605T064500Z \
  --abpoa-bin /tmp/abPOA-evaluate-abpoa-as/bin/abpoa \
  --render-svg \
  --random-orders 5
```

abPOA was available and reported version `1.5.6`.

## Artifacts

The generated experiment artifact directory is:

`data/test_poasta_insertion_order_20260605T064500Z`

Key files:

- `summary.tsv`: one row per replacement, with exact path preservation,
  segment/link/path/path-step/segment-bp counts, SCC/self-loop counts,
  path-depth metrics, and white-space metrics.
- `orders.tsv`: the exact insertion order used by each replacement.
- `skipped.tsv`: skipped requested order variants.
- `gfas/*.gfa`: generated replacement GFAs.
- `renders/*.svg`: path-labeled SVG render for every generated replacement.
- `input_sequences.fa`: reconstructed 47 candidate sequences used as aligner
  input.

Runtime for the release experiment, including SVG rendering, was 1:12 wall time
with 1.46 GiB maximum RSS.

## Skipped Variant

`reference_chm13_first` was not run because no CHM13 or GRCh38 path was present
among the 47 focused paths.

## Exact Path Preservation

Exact path preservation is the only acceptance gate for this experiment. Every
tested replacement preserved all 47 input path sequences exactly.

| Method | Order | Missing | Extra | Changed |
|---|---:|---:|---:|---:|
| POASTA | longest_then_name | 0 | 0 | 0 |
| POASTA | reverse_longest_then_name | 0 | 0 | 0 |
| POASTA | medoid_first_then_longest | 0 | 0 | 0 |
| POASTA | nearest_neighbor_guidetree | 0 | 0 | 0 |
| POASTA | random_seed_00 | 0 | 0 | 0 |
| POASTA | random_seed_01 | 0 | 0 | 0 |
| POASTA | random_seed_02 | 0 | 0 | 0 |
| POASTA | random_seed_03 | 0 | 0 | 0 |
| POASTA | random_seed_04 | 0 | 0 | 0 |
| abPOA | longest_then_name | 0 | 0 | 0 |
| abPOA | nearest_neighbor_guidetree | 0 | 0 | 0 |

## Graph Shape Metrics

The POASTA graph shape is invariant across every tested order. Segment IDs and
path-step encodings are not always byte-identical, but the measured topology and
path-depth/white-space diagnostics are unchanged.

| Method | Order | Segments | Links | Paths | Path steps | Segment bp |
|---|---:|---:|---:|---:|---:|---:|
| POASTA | longest_then_name | 72 | 96 | 47 | 1900 | 6484 |
| POASTA | reverse_longest_then_name | 72 | 96 | 47 | 1900 | 6484 |
| POASTA | medoid_first_then_longest | 72 | 96 | 47 | 1900 | 6484 |
| POASTA | nearest_neighbor_guidetree | 72 | 96 | 47 | 1900 | 6484 |
| POASTA | random_seed_00 | 72 | 96 | 47 | 1900 | 6484 |
| POASTA | random_seed_01 | 72 | 96 | 47 | 1900 | 6484 |
| POASTA | random_seed_02 | 72 | 96 | 47 | 1900 | 6484 |
| POASTA | random_seed_03 | 72 | 96 | 47 | 1900 | 6484 |
| POASTA | random_seed_04 | 72 | 96 | 47 | 1900 | 6484 |
| abPOA | longest_then_name | 6484 | 6508 | 47 | 239974 | 6484 |
| abPOA | nearest_neighbor_guidetree | 6484 | 6508 | 47 | 239974 | 6484 |

abPOA preserves paths exactly but produces a base-resolution exact graph after
normalization for this input. It is therefore not a compact replacement for the
current POASTA result on this case.

## SCC And Self-Loop Metrics

No tested replacement produced directed SCCs larger than one oriented node,
directed self-loop edges, segment self-loop links, or cyclic oriented nodes.

| Method | Order set | directed_sccs_gt1 | directed_self_loop_edges | segment_self_loop_links | largest_directed_scc | cyclic_oriented_nodes |
|---|---:|---:|---:|---:|---:|---:|
| POASTA | all 9 orders | 0 | 0 | 0 | 0 | 0 |
| abPOA | both orders | 0 | 0 | 0 | 0 | 0 |

## Path-Depth And White-Space Metrics

The order-invariant POASTA values are:

- bp-weighted path depth: 37.010179
- segment white-space bp fraction: 0.212549
- segment white-space bp total: 64774
- path depth median/p95/max: 36/47/47
- path white-space bridges >=1000 bp: 37

POASTA path-level white-space percentiles vary only with which path appears
first in the insertion order:

| POASTA order | First path bp | path_white_space_bp_p95 | path_white_space_bp_p99 | path_white_space_bp_max |
|---|---:|---:|---:|---:|
| longest_then_name | 6461 | 317 | 5950 | 5950 |
| reverse_longest_then_name | 93 | 51 | 5948 | 5948 |
| medoid_first_then_longest | 6461 | 317 | 5950 | 5950 |
| nearest_neighbor_guidetree | 6461 | 317 | 5950 | 5950 |
| random_seed_00 | 93 | 318 | 5948 | 5948 |
| random_seed_01 | 6461 | 317 | 5950 | 5950 |
| random_seed_02 | 6461 | 317 | 5950 | 5950 |
| random_seed_03 | 6460 | 317 | 5950 | 5950 |
| random_seed_04 | 6460 | 317 | 5950 | 5950 |

The segment-level white-space total, segment-level white-space fraction,
path-depth values, and graph size are stable across POASTA orders, so the small
path-level percentile changes are diagnostic artifacts rather than evidence of
meaningfully different alignment behavior.

## SVG Renders

Path-labeled renders were produced for every generated replacement:

- `renders/poasta_longest_then_name.svg`
- `renders/poasta_reverse_longest_then_name.svg`
- `renders/poasta_medoid_first_then_longest.svg`
- `renders/poasta_nearest_neighbor_guidetree.svg`
- `renders/poasta_random_seed_00.svg`
- `renders/poasta_random_seed_01.svg`
- `renders/poasta_random_seed_02.svg`
- `renders/poasta_random_seed_03.svg`
- `renders/poasta_random_seed_04.svg`
- `renders/abpoa_longest_then_name.svg`
- `renders/abpoa_nearest_neighbor_guidetree.svg`

The POASTA renders differ in label/layout details where segment IDs differ, but
the graph shape metrics above show that the visible splitting is not explained
by the path inclusion order on this focused 47-sequence case.

## Conclusion

Path inclusion order does not explain the visible POASTA underalignment or
splitting for this C4 crush04 replacement. Across current longest-then-name,
reverse, medoid-first, nearest-neighbor/guide-tree, and five deterministic
random orders, POASTA always emitted the same compact diagnostic shape: 72
segments, 96 links, 1900 path steps, 6484 segment bp, no SCCs or self-loops, and
identical segment-level white-space/depth metrics.

Recommendation: keep the current deterministic longest-then-name default. It
is stable, simple, and produced the same diagnostics as guide ordering here.
The explicit guide-order helper should remain available for experiments, but
this focused case does not justify changing the production default.
