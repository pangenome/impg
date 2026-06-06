# Graph Characterization

- Source: `data/c4_smoke_report_20260606T075500Z/validation/c4_2seq_renamed.source_interval.gfa`
- Status: `REVIEW`
- Size: 2 segments, 0 links, 2 paths, 2 path steps, 452 segment bp
- Review flags: `components>1, largest_component_frac, internal_tips>0, common_start_frac, common_end_frac`

## Large-Scale Architecture

- Path-depth profile across segment order: min 1, median 1, p95 1, max 1; reused nodes above path count: 0
- Node path coverage: 2 total path-step visits over 2 nodes; mean 1.00, bp-weighted mean 1.00; p10/median/p90 1/1/1
- Singleton nodes (coverage=1): 2 nodes, 452 bp; high-coverage shared nodes (coverage>=2): 0 nodes, 0 bp
- Node coverage histogram: 1:2 nodes/452 bp
- Path-by-segment occupancy: 0.500 covered, 0.500 white-space cells (452 missing path-bp); node occupancy fractions p05 0.500, median 0.500, p95 0.500
- Highest path-depth nodes:
  - `src1` order 0: total 1, paths 1, max/path 1 in `C4SMOKE_A#0#chr6:0-226`, seq_len 226, seq `CCTCGGTCTCGGTGTTTGTGGACCATCACCTGGCACCCTCCTTCTCTT...`
  - `src2` order 1: total 1, paths 1, max/path 1 in `C4SMOKE_B#0#chr6:0-226`, seq_len 226, seq `CCTCGGTCTCGGTGTTTGTGGACCATCACCTGGCACCCTCCTTCTCTT...`

## Topology

- Components: 2 (largest: 1 nodes, 0.500)
- Tips: 2 total, 2 internal/non-endpoint
- Most common path start: `src1+` (1/2, 0.500)
- Most common path end: `src1+` (1/2, 0.500)

## Layout And Adjacency

- Link jumps in current segment order: p95 0, p99 0, max 0
- Consecutive path-step jumps: p95 0, p99 0, max 0
- Path white-space bridges in current segment order: 0 nonzero, 0 >= 1000 bp; total 0 bp, mean 0.0, p95 0, p99 0, max 0

## Repeat And Glue Signals

- Duplicate segment sequences: 0 groups, 0 nodes, max copy count 1, node fraction 0.000
- Rare repeated local contexts: 0 nodes, 0 occurrences
- Direct self-loop L edges: 0 edge(s) on 0 node(s); adjacent same-node path steps 0, same-signed repeat steps 0, repeat runs 0, max run 0

## Interpretation

- The graph has 2 connected components; disconnected debris or over-split intervals may be present.
- Path starts or ends are dispersed; this can indicate contig breaks, incomplete interval merging, or uncollapsed duplicated boundary nodes.
- Current segment-order link jumps are bounded relative to graph size.
- The path-by-segment matrix has substantial white space; inspect the sparse-coverage runs to distinguish real insertion alleles from repeat-driven underalignment.
- Mean node coverage is 1.00 path-step visits per node; bp-weighted mean path-depth is 1.00 haplotype bases per stored segment base. Prefer the bp-weighted value when comparing graph condensation because it discounts many tiny high-depth nodes and emphasizes shared representation across longer sequence.
- 2 singleton node(s) covering 452 bp have coverage=1; persistent singleton bp usually marks private sequence, underalignment, or unshared replacement fragments.
- No rare repeated local-context signal was detected under the configured dominance threshold.
