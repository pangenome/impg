# Crush Chain-Greedy Walk

Task: `crush-chain-greedy`

Date: 2026-05-27

## Summary

This adds `method=chain-greedy`, a greedy reference-path walk for `impg crush`.
The mode still uses POVU to discover candidate bubble boundaries, but it does
not use the POVU parent/child tree to decide replacement blocks. Instead it:

1. Sorts discovered candidates by the root path interval.
2. Drops strict root-path containers so the walk sees the smallest adjacent
   intervals.
3. Walks the root path and groups consecutive non-overlapping candidates into
   chains until the configured span cap is reached.
4. Extracts all path traversals across the outer chain anchors.
5. Builds one POASTA replacement graph for the whole chain and applies normal
   path-preserving replacement validation.

The CLI accepts `--method chain-greedy` for standalone `impg crush` and
`method=chain-greedy` inside `gfa:syng:crush,...`. Chain size is controlled by
`--chain-target-bp` / `chain-target-bp`; default is `10k`.

The C4 result is path-preserving and fast, but it is not a graph-quality win.
It forms the intended adjacent-bubble chains, but the POASTA whole-chain
replacements increase segment count, stored sequence bp, duplicate-sequence
extras, and the internal quality score relative to the syng+mask input. It is
therefore useful as a tested experimental knob, not a new recommended default.

## Implementation

Changed files:

- `src/resolution.rs`
  - Added `ResolutionMethod::ChainGreedy`.
  - Added `ResolutionConfig::chain_greedy_target_bp`.
  - Added `resolve_graph_bubble_chains`, which bypasses the normal tree
    descent loop and drives selection from `find_path_walk_chain_frontier`.
  - Added root-path interval filtering, path-walk ordering, greedy chain
    grouping, chain materialization, per-round chain logging, and chain-size
    distribution logging.
  - Added `chain_greedy_groups_adjacent_reference_bubbles`, a synthetic test
    proving two adjacent reference bubbles are grouped into one replacement
    block while path sequences are preserved.
- `src/main.rs`
  - Added parser support for `method=chain-greedy`.
  - Added `--chain-target-bp` and engine-stage aliases:
    `chain-target-len`, `chain-target-size`, `target-chain-bp`,
    `chain-size`, `chain-bp`, `greedy-chain-target-bp`.
  - Added CLI parser coverage for `gfa:syng:crush,method=chain-greedy`.

## Focused Tests

```bash
CFLAGS='-I/home/erikg/wfmash/build/vendored_htslib/include -I/home/erikg/micromamba/include' \
LDFLAGS='-L/home/erikg/wfmash/build/vendored_htslib/lib -L/home/erikg/micromamba/lib' \
CMAKE_PREFIX_PATH='/home/erikg/wfmash/build/vendored_htslib;/home/erikg/micromamba' \
cargo test --lib chain_greedy_groups_adjacent_reference_bubbles -- --nocapture
```

Result: passed.

```bash
CFLAGS='-I/home/erikg/wfmash/build/vendored_htslib/include -I/home/erikg/micromamba/include' \
LDFLAGS='-L/home/erikg/wfmash/build/vendored_htslib/lib -L/home/erikg/micromamba/lib' \
CMAKE_PREFIX_PATH='/home/erikg/wfmash/build/vendored_htslib;/home/erikg/micromamba' \
cargo test --bin impg test_gfa_output_format_accepts_chain_greedy_crush_params -- --nocapture
```

Result: passed.

Full gate:

```bash
CFLAGS='-I/home/erikg/wfmash/build/vendored_htslib/include -I/home/erikg/micromamba/include' \
LDFLAGS='-L/home/erikg/wfmash/build/vendored_htslib/lib -L/home/erikg/micromamba/lib' \
CMAKE_PREFIX_PATH='/home/erikg/wfmash/build/vendored_htslib;/home/erikg/micromamba' \
cargo test --all
```

Result: passed.

## C4 Validation Command

Input graph:

```text
/home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.gfa
```

Output directory:

```text
/home/erikg/impg/data/c4_crush_chain_greedy_walk_20260527T1240Z
```

Command:

```bash
/usr/bin/time -v target/release/impg crush \
  --gfa /home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.gfa \
  --output /home/erikg/impg/data/c4_crush_chain_greedy_walk_20260527T1240Z/c4-crush-chain-greedy-walk.gfa \
  --method chain-greedy \
  --chain-target-bp 5k \
  --max-iterations 2 \
  -v 2 \
  > /home/erikg/impg/data/c4_crush_chain_greedy_walk_20260527T1240Z/c4-crush-chain-greedy-walk.stdout \
  2> /home/erikg/impg/data/c4_crush_chain_greedy_walk_20260527T1240Z/c4-crush-chain-greedy-walk.stderr
```

Runtime from `/usr/bin/time -v`:

| Metric | Value |
|---|---:|
| Exit status | 0 |
| Wall time | 2:28.89 |
| User CPU | 329.47 s |
| System CPU | 4.60 s |
| CPU | 224% |
| Max RSS | 2,438,228 KiB |

## Per-Iteration Chain Logs

The run resolved 52 chain blocks in two rounds, with no bailed replacements.

```text
crush chain-greedy per-round chains: [r1=37, r2=15]; total resolved=52
crush: 52 resolved, 0 bailed, 52 candidates seen across 2 rounds
```

Round 1:

| Field | Value |
|---|---:|
| POVU sites | 2,441 |
| Candidate bubbles | 2,437 |
| Path-walk bubbles | 1,309 |
| Chains formed | 45 |
| Chains selected/resolved | 37 |

Round 1 selected chain distribution:

| Metric | Median | P90 | Max |
|---|---:|---:|---:|
| Source bubbles per chain | 29 | 34 | 37 |
| Root span bp | 4,908 | 4,985 | 4,993 |
| Total traversal bp | 2,236,695 | 2,314,608 | 2,912,101 |
| Max traversal bp | 4,945 | 5,193 | 36,139 |
| Median traversal bp | 4,902 | 4,974 | 5,047 |

Round 2:

| Field | Value |
|---|---:|
| POVU sites | 3,799 |
| Candidate bubbles | 3,799 |
| Path-walk bubbles | 3,221 |
| Chains formed | 23 |
| Chains selected/resolved | 15 |

Round 2 selected chain distribution:

| Metric | Median | P90 | Max |
|---|---:|---:|---:|
| Source bubbles per chain | 39 | 116 | 193 |
| Root span bp | 4,947 | 4,995 | 4,998 |
| Total traversal bp | 1,500,565 | 2,307,926 | 2,311,596 |
| Max traversal bp | 4,967 | 5,054 | 5,133 |
| Median traversal bp | 4,924 | 4,997 | 4,998 |

The observed chain distributions confirm that path adjacency, not POVU tree
membership, controlled replacement block formation. For example, round 2
formed one selected chain containing 193 adjacent source bubbles under a
5 kb root-span target.

## Path Preservation

Independent oriented-path sequence comparison:

```text
input  {'segments': 18048, 'links': 20944, 'paths': 465, 'segment_bp': 389316, 'unique_seq': 12721, 'dup_extra': 5327, 'path_steps': 4591857}
output {'segments': 21684, 'links': 27415, 'paths': 465, 'segment_bp': 540986, 'unique_seq': 10790, 'dup_extra': 10894, 'path_steps': 3625741}
path_preservation {'before': 465, 'after': 465, 'missing': 0, 'extra': 0, 'mismatched': 0}
```

Gate result: 465/465 paths preserved by name and oriented sequence.

## Final C4 Metrics

| Metric | Input syng+mask | Chain-greedy output | Delta |
|---|---:|---:|---:|
| Segments | 18,048 | 21,684 | +3,636 |
| Links | 20,944 | 27,415 | +6,471 |
| Paths | 465 | 465 | 0 |
| Segment bp | 389,316 | 540,986 | +151,670 |
| Unique segment sequences | 12,721 | 10,790 | -1,931 |
| Duplicate-sequence extras | 5,327 | 10,894 | +5,567 |
| Path steps | 4,591,857 | 3,625,741 | -966,116 |

Internal quality log:

| Round | Quality score | Segments | Segment bp | Links | Path steps | ws-total | ws-p99 | ws-max | ws-long >=10 kb |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Input | 217,455,065 | 18,048 | 389,316 | 20,934 | 4,591,857 | 24,770,698,526 | 177,792 | 388,883 | 208,494 |
| After round 1 | 224,958,710 | 19,331 | 479,878 | 24,418 | 3,673,000 | 24,279,368,823 | 187,501 | 475,595 | 178,424 |
| After round 2 | 265,712,116 | 21,684 | 540,986 | 27,415 | 3,625,741 | 29,736,176,652 | 223,551 | 536,727 | 195,046 |

The output reduces total path steps, but graph topology and stored sequence
quality regress. The second iteration especially worsens the quality score
because white-space tail metrics and segment bp grow.

## Comparison To Existing Baselines

The closest existing iterative/multi baselines are the prior C4 auto-routing
runs documented in `docs/crush-experiment-synthesis.md` and
`docs/crush-smoothxg-on-output.md`. PGGB control metrics are from
`docs/crush-vs-pggb-comparison.md`.

| Graph | Wall | Segments | Links | Segment bp | Paths | Notes |
|---|---:|---:|---:|---:|---:|---|
| Input syng+mask | - | 18,048 | 20,944 | 389,316 | 465 | Baseline graph, already path-preserving |
| Auto k=311 iterative/multi | 32:55 | 19,968 | - | 570,180 | 465 | Best synthesis row by quality score |
| `method=auto` + no-filter prior best | 36:53 | 19,836 | 23,384 | 553,585 | 465 | Prior best in PGGB comparison docs |
| Chain-greedy walk, this task | 2:28.89 | 21,684 | 27,415 | 540,986 | 465 | Faster, path-preserving, worse topology |
| PGGB control | 13:38 | 13,288 | 16,240 | 234,524 | 465 | Best graph-quality target |

Interpretation:

- Compared with the auto/no-filter prior best, chain-greedy is much faster and
  stores slightly less segment bp (540,986 vs 553,585), but it emits more
  segments and links.
- Compared with the synthesis auto k=311 iterative/multi run, chain-greedy has
  fewer segment bp (540,986 vs 570,180) but more segments (21,684 vs 19,968)
  and a worse final quality score (265,712,116 vs 193,762,774).
- Compared with PGGB, chain-greedy is still far behind: 1.63x as many segments,
  1.69x as many links, and 2.31x as much stored segment sequence.
- The greedy walk hypothesis is partially confirmed operationally: it does
  combine long runs of adjacent small bubbles. The quality hypothesis is not
  confirmed on C4 with POASTA as the whole-chain aligner.

## Visualization

Sorted GFA:

```bash
gfasort \
  -i /home/erikg/impg/data/c4_crush_chain_greedy_walk_20260527T1240Z/c4-crush-chain-greedy-walk.gfa \
  -o /home/erikg/impg/data/c4_crush_chain_greedy_walk_20260527T1240Z/c4-crush-chain-greedy-walk.Ygs.gfa \
  -p Ygs -t 32
```

PNG:

```bash
gfalook \
  -i /home/erikg/impg/data/c4_crush_chain_greedy_walk_20260527T1240Z/c4-crush-chain-greedy-walk.Ygs.gfa \
  -o /home/erikg/impg/data/c4_crush_chain_greedy_walk_20260527T1240Z/c4-crush-chain-greedy-walk.png \
  -m -x 2200 -y 1200
```

Upload:

```bash
scp /home/erikg/impg/data/c4_crush_chain_greedy_walk_20260527T1240Z/c4-crush-chain-greedy-walk.png \
  erik@hypervolu.me:www/impg/c4-crush-chain-greedy-walk.png
```

Confirmed with `ssh ls`:

```text
-rw-r--r-- 1 erik erik 2.2M May 27 12:49 www/impg/c4-crush-chain-greedy-walk.png
```

Public URL:

```text
https://hypervolu.me/~erik/impg/c4-crush-chain-greedy-walk.png
```

## Hard-Gate Checklist

- [x] Code change in `src/resolution.rs`.
- [x] CLI parser change in `src/main.rs`.
- [x] `cargo test --all` passes.
- [x] Real C4 run completed with exit status 0.
- [x] Real C4 run preserved 465/465 paths by oriented sequence.
- [x] Per-iteration chains-formed count logged.
- [x] Per-iteration chain-size distribution logged.
- [x] Final metrics recorded.
- [x] PNG uploaded as `c4-crush-chain-greedy-walk.png`.
- [x] PNG upload confirmed with `ssh ls`.
- [x] This document committed as `docs/crush-chain-greedy-walk.md`.
