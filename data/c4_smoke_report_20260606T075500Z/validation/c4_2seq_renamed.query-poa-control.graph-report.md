# Graph Characterization

- Source: `data/c4_smoke_report_20260606T075500Z/outputs/c4_2seq_renamed.query-poa-control.gfa`
- Status: `REVIEW`
- Size: 4 segments, 4 links, 2 paths, 6 path steps, 227 segment bp
- Review flags: `link_jump_max_frac`

## Large-Scale Architecture

- Path-depth profile across segment order: min 1, median 1, p95 2, max 2; reused nodes above path count: 0
- Node path coverage: 6 total path-step visits over 4 nodes; mean 1.50, bp-weighted mean 1.99; p10/median/p90 1/1/2
- Singleton nodes (coverage=1): 2 nodes, 2 bp; high-coverage shared nodes (coverage>=2): 2 nodes, 225 bp
- Node coverage histogram: 1:2 nodes/2 bp, 2:2 nodes/225 bp
- Path-by-segment occupancy: 0.996 covered, 0.004 white-space cells (2 missing path-bp); node occupancy fractions p05 0.500, median 0.500, p95 1.000
- Highest path-depth nodes:
  - `1` order 0: total 2, paths 2, max/path 1 in `C4SMOKE_A#0#chr6:0-226`, seq_len 65, seq `CCTCGGTCTCGGTGTTTGTGGACCATCACCTGGCACCCTCCTTCTCTT...`
  - `4` order 3: total 2, paths 2, max/path 1 in `C4SMOKE_A#0#chr6:0-226`, seq_len 160, seq `TGGAGACCACCAGTGGCCAACTCCCTGCGAGTGGATGTCCAGGCTGGG...`
  - `3` order 2: total 1, paths 1, max/path 1 in `C4SMOKE_B#0#chr6:0-226`, seq_len 1, seq `C`
  - `2` order 1: total 1, paths 1, max/path 1 in `C4SMOKE_A#0#chr6:0-226`, seq_len 1, seq `A`
- Highest-depth order runs:
  - orders 0-0 (1 nodes, `1` to `1`): depth min 2, max 2, mean 2.00
  - orders 3-3 (1 nodes, `4` to `4`): depth min 2, max 2, mean 2.00

## Topology

- Components: 1 (largest: 4 nodes, 1.000)
- Tips: 0 total, 0 internal/non-endpoint
- Most common path start: `1+` (2/2, 1.000)
- Most common path end: `4+` (2/2, 1.000)

## Layout And Adjacency

- Link jumps in current segment order: p95 2, p99 2, max 2
- Consecutive path-step jumps: p95 2, p99 2, max 2
- Path white-space bridges in current segment order: 2 nonzero, 0 >= 1000 bp; total 2 bp, mean 0.5, p95 1, p99 1, max 1
- Top long links:
  - `1` -> `3`: jump 2, path_support 1
  - `2` -> `4`: jump 2, path_support 1
  - `3` -> `4`: jump 1, path_support 1
  - `1` -> `2`: jump 1, path_support 1
- Top path white-space bridges:
  - `C4SMOKE_A#0#chr6:0-226` step 1: `2+` -> `4+` bridges 1 bp (orders 1 -> 3, gap bp 66-67)
  - `C4SMOKE_B#0#chr6:0-226` step 0: `1+` -> `3+` bridges 1 bp (orders 0 -> 2, gap bp 65-66)
- Top path jumps:
  - `C4SMOKE_A#0#chr6:0-226` step 1: `2+` -> `4+` jump 2
  - `C4SMOKE_B#0#chr6:0-226` step 0: `1+` -> `3+` jump 2
  - `C4SMOKE_A#0#chr6:0-226` step 0: `1+` -> `2+` jump 1
  - `C4SMOKE_B#0#chr6:0-226` step 1: `3+` -> `4+` jump 1

## Repeat And Glue Signals

- Duplicate segment sequences: 0 groups, 0 nodes, max copy count 1, node fraction 0.000
- Rare repeated local contexts: 0 nodes, 0 occurrences
- Direct self-loop L edges: 0 edge(s) on 0 node(s); adjacent same-node path steps 0, same-signed repeat steps 0, repeat runs 0, max run 0

## Interpretation

- The graph is connected in the undirected S/L topology.
- Most paths share common entry and exit nodes, consistent with a focused local locus.
- At least one long sorted-order link is also supported by path adjacency, so it is not just a renderer artifact.
- Mean node coverage is 1.50 path-step visits per node; bp-weighted mean path-depth is 1.99 haplotype bases per stored segment base. Prefer the bp-weighted value when comparing graph condensation because it discounts many tiny high-depth nodes and emphasizes shared representation across longer sequence.
- 2 singleton node(s) covering 2 bp have coverage=1; persistent singleton bp usually marks private sequence, underalignment, or unshared replacement fragments.
- No rare repeated local-context signal was detected under the configured dominance threshold.
