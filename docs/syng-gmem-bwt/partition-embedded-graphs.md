# The partition-embedded graph build and the expressibility census

Date: 2026-11-06. Branch `work/genome-mem-bwt-pipeline`. This note records
the owner-ruled partition-embedded graph build (BUILD AND MEASURE only — no
inference rewiring), the exact construction path, the receipts, and the
expressibility census headline. The per-arc working details live in the
running report (`genome/` artifacts,
`context-and-mem-likelihood.md`); this note preserves the construction path
and the census semantics against session compaction.

## 1. The ruling

Instead of windowed content-homed universes, build per-partition embedded
pangenome graphs — the pggb-equivalent induction via the impg graph
machinery — exported as GFA for the partitions covering chrMT and chrI.
KEY PROPERTY: the partition graph must be the panel's own alignment-induced
structure restricted to the partition's material, with syncmer node ids
CONTINUOUS with the global syng interning (frame-independent identity) so
the proven likelihood instrument and truth index work unchanged. Material
is in the partition graph because it GENUINELY ALIGNS into the partition
locality — no content-homing, no window territories.

## 2. The construction path (stated exactly)

- **Material**: the panel's own alignment-induced partitioning — the
  COMPLETED `impg partition` run
  (`/home/erikg/yeast/partition-pos64-w10k-d1k-to-completion`, source
  commit 295bca9, `-w 10000 -d 1000 --selection-mode haplotype
  --syng-min-chain-anchors 5 --syng-min-chain-fraction 0.5
  --no-rehome-singletons`, exit 0, wall 10,544 s, 19,421 partitions,
  3,336,976,856/3,336,976,856 bp = 100.0000%). Membership is purely
  alignment-induced (syng chain criteria). The run is NOT re-executed.
- **Graph**: `examples/partition_graph_export.rs` renders each partition's
  member rows as `SyngGfaPathRange`s through
  `impg::commands::syng2gfa::write_range_gfa_with_mode` — the same writer
  `impg syng2gfa` and `impg query -o gfa --gfa-engine syng` use — in RAW
  mode with `SyngGfaFrequencyMask::disabled()`: no private splits, no
  clones, no frequency removal, so every syncmer node keeps its GLOBAL
  interning abs id (the key property decides the mode; any masking
  re-interns nodes and breaks instrument/truth-index continuity).
  Inter-syncmer gap segments are spliced from real panel DNA (the AGC)
  and interned by the writer strictly after `num_syncmer_nodes()`, so S
  names never alias syncmer ids.
- **Scope**: 50 partitions — the 35 axis partitions covering the chrMT/chrI
  windows (chrI partitions 0–20; chrMT 14900, 6382, 14901, 14902, 15582,
  15968, 14903–14906, 15586, 16122, 16242, 16360), the truth-side
  partitions holding SK1 material for the previously non-expressible
  loci (chrMT 16239/15966/15972; chrI 14590/14609/16771/14610/541/3865/
  14482/14606/3093/1899 plus 102/3134/18304/14346), and the two
  near-twin-holder partitions 15970/16121.
- **Receipts** (`/home/erikg/yeast/genome-balanced-diploid-validation-20260930/partition-graphs/`):
  `partition<N>.gfa` + `partition<N>.gfa.map.json` (the node-id mapping:
  S names ≤ `n_syncmer_nodes` ARE the global interning ids; gap segments
  are the ids beyond it), `build-partition-graphs.{start,done,log,wall}`
  (wall 13 s), the alignment receipts
  (`partition-graph-homology-requests.tsv`, `partition-graph-homology.jsonl`
  — the syng's own `query_region_with_anchors` from every window extent,
  padding 120, the probe's own existing choice), the walk receipts
  (`partition-graph-walk-requests.tsv`, `partition-graph-walks.jsonl` —
  the panel's own path walks for every truth/twin member row), the census
  (`partition-graph-census.jsonl`, `partition-graph-census.tables.txt`),
  and the independent checker log `check-partition-graphs.log` — ALL
  PHASES PASS.

## 3. The census semantics (no thresholds)

- **(a) the 805-class test**: per previously non-expressible locus (the
  committed likelihood receipts' `truth_piece_presence`): the SK1
  POSITIONAL region (the window's coordinates on the SK1 truth contig — a
  coordinate statement) tiled by the partition structure's SK1 rows
  (which partitions hold it, with intervals), plus the syng's own
  window-query hits (every hit reported with anchors; the best hit named).
  Verdicts: IN-AXIS-PARTITION / TILED-ELSEWHERE / PARTIAL-ELSEWHERE /
  ABSENT.
- **(b) the pockets**: per chrI tail locus, every off-truth panel path
  COALESCED with the truth pair's material at the window (shared covered
  graph nodes, from the multi-census receipts) on contigs outside the
  chrI family — the pocket carriers; verdict by the syng's window-query
  hit plus the partition placement of the carrier's rows: IN-BY-ALIGNMENT
  (the axis partition holds a row overlapping the alignment extent),
  OUT-ALIGNED-ELSEWHERE, OUT-UNPLACED, NO-ALIGNMENT-INTO-WINDOW.
- **(c) the near-twins**: membership rows across the partition structure;
  per common BUILT partition the GFA spells give shared/left-only/
  right-only node sets; variant pockets as bp runs from the panel's own
  walks.
- **(d) the truth spells**: per locus, S288C rows and SK1 positional
  tiling (bp covered of the window), the same verdict classes, plus the
  window-query positional hit bp.
- **Id continuity is MEASURED, not assumed**: the checker's spell-equality
  phase re-derives every truth/twin member row's panel walk and demands
  the GFA P line spells exactly the same signed global node ids in order
  (109/109 walks).

## 4. The headline (measured; the running report carries the full tables)

- **(a) 8/8** previously non-expressible loci have SK1 material present by
  alignment in the partition structure (7 fully tiled, chrMT L13 tiled to
  SK1's contig end — the remaining 1,155 bp is SK1's shorter chrMT, a
  length polymorphism, not a placement gap). None are in the window's own
  axis partition — the frame gap closes in the partition universe, in
  NAMED partitions.
- **(b) the pockets GENUINELY ALIGN IN**: the big conserved scaffold
  pockets (BTE/ANL/ANM/BBT/CRE/ATV/CGH/BAD block*_contig1 and the
  chrVIII/chrVII/chrX/chrV/chrXI copies at L16) are IN-BY-ALIGNMENT with
  hundreds of anchors and ~10 kb extents tiling ALONG the chrI axis
  partitions; the windowed frame's content-homing was homing real aligned
  material.
- **(c)** the near-twin routes run on shared nodes (chrMT CBK/CBM share
  3,552 nodes across their 14 common partitions with 3 vs 5 variant bp
  runs; chrI BTE#3/#4 block28 share 100% — identical haplotypes through
  the graph).
- **(d)** at every expressible locus BOTH homologs are in the SAME axis
  partition (the truth pair is jointly spellable in one partition graph);
  the only untiled bp anywhere is chrMT L13's contig-end tail.

No selection change, no inference rewiring, no thresholds; the scoreboard
machinery is untouched; assessment-side only.
