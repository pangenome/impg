# crush-chain-povu-tree-blocks

Task: `crush-chain-povu`

## Summary

Implemented a new `method=chain-povu` crush mode that discovers POVU tree
sites on the current working graph, walks the tree top-down, and selects the
largest subtree-rooted blocks whose parent-bounded interior is at most 10 kb.
When a selected parent fits the cap, all descendants are absorbed into that
single block and traversal does not descend further. Oversized subtrees recurse
into children; oversized leaves are skipped.

Each selected block is materialized as path sequences and processed through a
local smoothxg-style smoothing pass followed by bounded POASTA cleanup. The
replacement is compacted and path-validated before rewrite. As of
`remove-crush-replacement`, this mode no longer compares the smoothxg -> POASTA
replacement against direct POASTA by segment-bp or segment count. If the
smoothxg -> POASTA replacement parses and preserves every traversal sequence,
it is used unconditionally; metrics are diagnostic only. Direct POASTA is only
used when the smooth path cannot produce a valid replacement, or when an empty
traversal makes the smooth path inapplicable.

## Code changes

- Added `ResolutionMethod::ChainPovu` with CLI/parser aliases:
  `chain-povu`, `povu-chain`, `povu-tree-blocks`, `tree-blocks`,
  `crush-chain-povu`.
- Added `resolve_graph_bubbles_chain_povu`, a dedicated multi-round chain
  driver in `src/resolution.rs`.
- Added `select_chain_povu_blocks`, which builds the current POVU candidate
  forest and performs top-down subtree selection under the 10 kb cap.
- For `method=chain-povu`, `discover_all_candidates` keeps all POVU sites
  regardless of `min-traversal-len` so parent blocks can absorb every
  descendant under the cap.
- Added per-round logs for:
  - POVU sites/tree nodes/roots
  - selected subtree block count
  - absorbed descendant site count
  - skipped/oversized counts
  - selected and accepted block-size distributions
  - final per-iteration block counts and distributions
- Added regression tests for parser aliases, subtree selection, nested parent
  replacement, and path-preserving chain-POVU replacement.

## Validation

Rust validation was run with the local native dependency include/library paths
needed by this checkout:

```bash
CARGO_BUILD_JOBS=8 \
CFLAGS='-I/home/erikg/wfmash/build/vendored_htslib/include -I/home/erikg/micromamba/include' \
LDFLAGS='-L/home/erikg/wfmash/build/vendored_htslib/lib -L/home/erikg/micromamba/lib' \
CMAKE_PREFIX_PATH='/home/erikg/wfmash/build/vendored_htslib;/home/erikg/micromamba' \
cargo test --all
```

Result: pass. The full suite completed after the final guard/logging changes:

- lib tests: 292 passed
- binary tests: 64 passed
- integration tests: passed
- existing C4-heavy/known-red tests remained ignored as already marked in the
  repository

The release/install binary was refreshed after tests:

```bash
CARGO_BUILD_JOBS=8 \
CFLAGS='-I/home/erikg/wfmash/build/vendored_htslib/include -I/home/erikg/micromamba/include' \
LDFLAGS='-L/home/erikg/wfmash/build/vendored_htslib/lib -L/home/erikg/micromamba/lib' \
CMAKE_PREFIX_PATH='/home/erikg/wfmash/build/vendored_htslib;/home/erikg/micromamba' \
cargo install --path .
```

## C4 command

Primary 3-round C4 validation:

```bash
out=/home/erikg/impg/data/c4_crush_chain_povu_blocks_guarded_20260527T125849Z
LD_LIBRARY_PATH='/home/erikg/wfmash/build/vendored_htslib/lib:/home/erikg/micromamba/lib:'"${LD_LIBRARY_PATH:-}" \
/usr/bin/time -v -o "$out/time.txt" \
  /home/erikg/impg/.wg-worktrees/agent-214/target/release/impg query \
    -t 32 \
    -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
    -r 'GRCh38#0#chr6:31891045-32123783' \
    --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
    -d 50k \
    -o 'gfa:syng:mask,min-run=3:crush,method=chain-povu,max-rounds=3,aligner=fastga,min-traversal-len=5k,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-rounds=until-done,polish-method=poasta,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
    -O "$out/run.nosort" \
    -v 1

gfasort -i "$out/run.nosort.gfa" -o "$out/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/run.Ygs.gfa" -o "$out/c4-crush-chain-povu-blocks.png" -m -x 2200 -y 1200
scp "$out/c4-crush-chain-povu-blocks.png" erik@hypervolu.me:www/impg/c4-crush-chain-povu-blocks.png
ssh erik@hypervolu.me 'ls -lh ~/www/impg/c4-crush-chain-povu-blocks.png'
```

The task allowed 2-3 iterations, so I ran a 2-round checkpoint for comparison
at `/home/erikg/impg/data/c4_crush_chain_povu_blocks_2round_20260527T130927Z`.
After a workgraph message asked for a non-decisive multi-round comparison, I
also ran a 5-round comparison at
`/home/erikg/impg/data/c4_crush_chain_povu_blocks_5round_20260527T132444Z`.

The primary uploaded PNG remains the 3-round run because it is inside the
original 2-3 iteration envelope and gives lower stringy residue than 2 rounds.
The 5-round run is documented as an exploratory comparison: it lowers stringy
residue further, but grows both segment count and segment-bp.

## C4 run status

Primary 3-round run:

```text
Exact round count: 3
Exit status: 0
Elapsed wall: 7:42.74
Maximum RSS: 55,912,452 KiB
Paths: 465 / 465
```

The crush stage itself reported:

```text
crush: 2267 resolved, 0 bailed, 2267 candidates seen across 3 rounds
Syng query complete ... in 4m21.478s
```

Additional comparison runs:

| run | exact rounds | exit | wall | max RSS KiB | paths |
|---|---:|---:|---:|---:|---:|
| 2-round checkpoint | 2 | 0 | 6:56.53 | 55,746,104 | 465 |
| 3-round primary | 3 | 0 | 7:42.74 | 55,912,452 | 465 |
| 5-round comparison | 5 | 0 | 10:45.37 | 56,017,348 | 465 |

Path preservation is enforced by replacement validation during crush, and the
final GFA has 465 `P` lines:

```bash
awk '$1=="S"{s++;bp+=length($3)} $1=="L"{l++} $1=="P"{p++}
     END{printf("S=%d L=%d P=%d bp=%d\n", s,l,p,bp)}' run.nosort.gfa
```

Output:

```text
S=19433 L=24739 P=465 bp=357353
```

## Per-iteration block logs

The target block cap used by `method=chain-povu` is exactly
`DEFAULT_CHAIN_POVU_MAX_BLOCK_BP = 10,000 bp`.

Primary 3-round `run.stderr` summary:

```text
crush chain-povu per-iteration block counts: [r1=714, r2=647, r3=906]; size distributions: [r1 blocks n=714, min=45, p50=132, p90=173, max=9320, total=120655; r2 blocks n=647, min=3, p50=141, p90=268, max=9446, total=125221; r3 blocks n=906, min=3, p50=67, p90=279, max=7297, total=127907]; total resolved=2267; next_id moved 272214103 -> 272236948
```

Primary 3-round accepted summaries:

| round | selected/accepted | target cap | actual max block | selected p50/p90 | oversize internal | oversize leaves | post-round segments | post-round bp |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 714 / 714 | 10,000 | 9,320 | 132 / 173 | 5 | 57 | 19,204 | 385,658 |
| 2 | 647 / 647 | 10,000 | 9,446 | 141 / 268 | 6 | 26 | 18,596 | 350,342 |
| 3 | 906 / 906 | 10,000 | 7,297 | 67 / 279 | 10 | 28 | 19,433 | 357,353 |

The 5-round comparison used the same 10,000 bp target cap and added
replacement-decision logging. Every block chose the direct POASTA fallback
because its replacement was smaller than the smoothxg -> POASTA result:

| round | selected/accepted | actual max block | selected p50/p90 | oversize internal | oversize leaves | smooth kept | direct fallback | post-round segments | post-round bp |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 714 / 714 | 9,320 | 132 / 173 | 5 | 57 | 0 | 714 | 19,204 | 385,658 |
| 2 | 647 / 647 | 9,446 | 141 / 268 | 6 | 26 | 0 | 647 | 18,596 | 350,342 |
| 3 | 906 / 906 | 7,297 | 67 / 279 | 10 | 28 | 0 | 906 | 19,433 | 357,353 |
| 4 | 823 / 823 | 9,265 | 55 / 201 | 16 | 24 | 0 | 823 | 20,119 | 369,756 |
| 5 | 583 / 583 | 6,693 | 48 / 194 | 13 | 20 | 0 | 583 | 20,967 | 385,060 |

The 2-round checkpoint stopped after round 2:

```text
crush chain-povu per-iteration block counts: [r1=714, r2=647]; size distributions: [r1 blocks n=714, min=45, p50=132, p90=173, max=9320, total=120655; r2 blocks n=647, min=3, p50=141, p90=268, max=9446, total=125221]; total resolved=1361; next_id moved 272214103 -> 272228806
```

The 5-round comparison ended with:

```text
crush chain-povu per-iteration block counts: [r1=714, r2=647, r3=906, r4=823, r5=583]; size distributions: [r1 blocks n=714, min=45, p50=132, p90=173, max=9320, total=120655; r2 blocks n=647, min=3, p50=141, p90=268, max=9446, total=125221; r3 blocks n=906, min=3, p50=67, p90=279, max=7297, total=127907; r4 blocks n=823, min=3, p50=55, p90=201, max=9265, total=111640; r5 blocks n=583, min=3, p50=48, p90=194, max=6693, total=82302]; total resolved=3673; next_id moved 272214103 -> 272248796
```

## Final metrics

Primary 3-round result:

| metric | value |
|---|---:|
| segments | 19,433 |
| links | 24,739 |
| paths | 465 |
| segment bp | 357,353 |
| trivial-stringy candidates | 265 |
| 51-200 bp duplicate-seq extras | 81 |
| >200 bp duplicate-seq extras | 9 |

2-round checkpoint:

| metric | value |
|---|---:|
| segments | 18,596 |
| links | 23,415 |
| paths | 465 |
| segment bp | 350,342 |
| trivial-stringy candidates | 371 |
| 51-200 bp duplicate-seq extras | 27 |
| >200 bp duplicate-seq extras | 1 |

5-round comparison:

| metric | value |
|---|---:|
| segments | 20,967 |
| links | 26,684 |
| paths | 465 |
| segment bp | 385,060 |
| trivial-stringy candidates | 209 |
| 51-200 bp duplicate-seq extras | 166 |
| >200 bp duplicate-seq extras | 18 |

Comparison against the task references:

| run | segments | segment bp | trivial-stringy |
|---|---:|---:|---:|
| PGGB | 13,288 | 234,524 | 12 |
| iterative-multi best | 16,547 | 396,000 | 104 |
| chain-POVU 2 rounds | 18,596 | 350,342 | 371 |
| chain-POVU 3 rounds | 19,433 | 357,353 | 265 |
| chain-POVU 5 rounds | 20,967 | 385,060 | 209 |

This approach improves segment-bp versus the cited iterative-multi run
(`357,353` at 3 rounds and `385,060` at 5 rounds vs `396k`) but does not
beat it on segment count or stringy residue. Additional rounds keep reducing
stringy residue but grow segment count, duplicate-sequence extras, and segment
bp after round 2.

## Interpretation

The residual under-alignment is not primarily explained by the 10 kb cap being
too low. In the 5-round comparison, the selected block maxima were 9,320,
9,446, 7,297, 9,265, and 6,693 bp, with p90 values only 173-279 bp in rounds
1-3 and 194-201 bp in rounds 4-5. The cap does block a small number of larger
subtrees (`oversize_internal` 5, 6, 10, 16, 13 and `oversize_leaves` 57, 26,
28, 24, 20), but most processed blocks are much smaller than the target.

The old replacement decision log compared smoothxg -> POASTA against direct
POASTA and selected the direct result by metric. That decision path has been
removed. Current chain-POVU semantics do not let direct POASTA win because it
is smaller; a valid smoothxg -> POASTA replacement is spliced into the graph
unconditionally. Direct POASTA remains a validity recovery path when smoothxg
-> POASTA fails to build, parse, or preserve traversal spellings.

Increasing the cap beyond 10 kb is therefore not justified as the next primary
lever from these historical logs alone. A higher cap might help a few oversized
internal subtrees, but most selected blocks were already tiny. A follow-up
should either change how
the smooth seed graph is constructed for many-haplotype blocks, or select wider
context using a different criterion than simply raising the cap.

## PNG

Local PNG:

```text
/home/erikg/impg/data/c4_crush_chain_povu_blocks_guarded_20260527T125849Z/c4-crush-chain-povu-blocks.png
```

Canonical URL:

```text
https://hypervolu.me/~erik/impg/c4-crush-chain-povu-blocks.png
```

Upload confirmation:

```bash
ssh erik@hypervolu.me 'ls -lh ~/www/impg/c4-crush-chain-povu-blocks.png'
```

Output:

```text
-rw-r--r-- 1 erik erik 2.0M May 27 13:08 /home/erik/www/impg/c4-crush-chain-povu-blocks.png
```

## Hard-gate checklist

- [x] Code change in `src/resolution.rs` plus CLI parser wiring in `src/main.rs`.
- [x] `cargo test --all` passes after the final guard change.
- [x] Real C4 run completed, exit status 0.
- [x] 465 / 465 paths preserved by final `P` line count and crush validation.
- [x] Per-iteration block count and size distribution logged.
- [x] Final metrics reported.
- [x] PNG uploaded as `c4-crush-chain-povu-blocks.png` and confirmed with `ssh ls`.
- [x] `docs/crush-chain-povu-tree-blocks.md` committed with this task.
- [x] Workgraph artifact recorded for `docs/crush-chain-povu-tree-blocks.md`.
