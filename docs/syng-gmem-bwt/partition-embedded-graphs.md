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

## 5. The anchor-projection layer (slice A of the local-realignment model)

Date: 2026-11-06 (slice A session). The owner's architecture ruling: the
local-realignment evidence model EXTENDS the syng GBWT MEM-to-read
projection — the MEMs are the anchor layer and place each read exactly on
interned nodes; the realignment projects only THE BASES THE MEMs SKIP
through the graph along the anchored placement. No new aligner, no
parallel machinery; the scoring is a later slice.

The layer (`examples/partition_anchor_projection.rs`, runner
`genome/instrumented/run-anchor-projection.sh`, checker
`genome/instrumented/check-anchor-projection.py`):

- **Per record** (of the committed multi-matching census — the
  read-matched stage's verified occurrence set): the canonical anchor
  chain (signed node, canonical position; span = last + k), the read
  variants (mirror flag, flank bases, instance count — reads of one
  record differ in flanks AND mirror state, so the PHYSICAL strand is
  per (occurrence x variant): `(orientation == 0) != mirrored`), and the
  occurrences census-verbatim plus the SEQUENCE-VERIFIED origin shift
  (measured 0 at every occurrence on both components — the committed
  census starts are exact placements).
- **Per occurrence**: the re-derived merged read-matched intervals and
  covered-node lists (the census's own recipe) reproduce the committed
  census EXACTLY at every occurrence (234,472 chrMT + 1,415,203 chrI);
  the read hull equals the placed segment at every (occurrence, variant)
  strand-aware (2.0M chrMT + 8.6M chrI checks).
- **The binding lesson** (the likely death site of the timed-out
  attempt): the census occurrence starts live in the routing's
  CANONICAL-SCHEME territory step convention — the syng's stored path
  walk (`walk_path_range`) carries only ONE frame's syncmer selection,
  so anchor lookups must reproduce the territory step extraction (the
  forward frame from the stored walk, the rc frame from the raw
  extraction on the fetched range's reverse complement, per position the
  frame that spells the k-mer's canonical form forward).
- **The per-path graph context** (the lookups the scoring slice
  consumes): per candidate row (the member rows of the axis partitions
  of every window the occurrence touches, sharing at least one anchor
  node with the required sign), the anchor correspondence (the
  unique-monotone rule over the row's committed walk; repeated nodes
  unresolved are named ambiguous) and per skipped base (hull-relative,
  shared by all variants of one mirror state): the row-side coordinate,
  the covering walk-index range (the partition GFAs' own step
  sequences; empty = an inter-window gap pocket), and the row's base or
  the out-of-extent flag. Emitted for an evenly-spaced stated sample of
  1,000 occurrences per component (the full-fidelity per-path context
  over every occurrence x every candidate row measures in the tens of
  GB and is NOT emitted; the scoring slice computes contexts in-process
  with this same machinery).

Receipts (in the data directory, beside the earlier stage artifacts):
`anchor-projection-{chrMT,chrI}.jsonl` (the full placement layer),
`anchor-projection-context-{chrMT,chrI}.jsonl` (+ `.rows.jsonl` — the
sampled candidate rows' walks and sequences, validated by the checker
against the committed GFAs' P/L/S lines), run markers
`run-anchorproj-{chrMT,chrI}.*` (exit 0, walls 153s/583s, RSS peaks
2.3/3.5 GB), and the checker log `check-anchor-projection.log` — ALL
PHASES PASS (41,435,430 checks, zero failures).

Measured shape: chrMT 4,891 records / 34,088 read instances / 234,472
occurrences / 33,036 read variants / 2,727,452 skipped bp (~80 bp per
instance); chrI 13,905 / 231,210 / 1,415,203 / 183,844 / 18,844,823
(~81.5 bp per instance); interior hull gaps ZERO on both components
(the committed zero-gap measurement reproduced); multi-placement
records 4,803 chrMT / 13,264 chrI. Context samples: 998 chrMT
occurrences with 51,373 traversing rows and 10.0M per-base lookups;
1,000 chrI occurrences with 98,305 traversing rows and 22.1M lookups.
Assessment-side only; no thresholds; no scoring built; the scoreboard
machinery is untouched.

## 6. The chrIV scale-up (slice 1: survey, partition graphs, expressibility)

Date: 2026-11-06 (chrIV end-to-end session). The same construction
path and census semantics, extended to the first 1.5 Mb chromosome —
the balanced-diploid validation's chrIV component (167 loci, truth
S288C+SK1, `truth_pair_in_candidate_domain` 0/167 in the windowed
frame; 98 assessable / 69 bracketed loci in the committed balanced
receipt, which also supplies each assessable locus's SK1 ORTHOLOG
interval — the alignment-truth statement, reported beside the
POSITIONAL coordinate statement the chrMT/chrI census used).

- **Survey** (receipt-side, from the partition BEDs + the axis file +
  the names file): 19,421 partitions genome-wide carry 29,458
  chrIV-family member rows (181 paths, 162 strains, 12 multi-copy)
  across 1,822 partitions. The 167 windows sit over 158 axis
  partitions; S288C#0#chrIV = 167 rows / 1,566,853 bp; SK1#0#chrIV =
  163 rows / 1,486,921 bp, 145 rows in axis partitions, 17 SK1 holder
  partitions (117, 292, 544, 2073, 2086, 2088, 2312, 2322, 2323,
  2331, 2899, 2908, 3483, 3484, 3567, 3715, 14624).
- **The near-twin flag rule** (stated, structural, pre-scoring): a
  pair of chrIV-family paths is flagged iff EVERY member row of one
  path has an interval-close counterpart (both ends within 50 bp) in
  the same partition of the other. This flags exactly two pairs —
  `AAA#0#chrIV`/`SGDref#0#chrIV` (100% EXACT-identical rows, the
  identical-through-graph class) and `CLL#0#chrIV`/`CLL#1#chrIV`
  (100% interval-close, 0% exact, offsets within ±16 bp, the
  near-identical same-strain copy class); the runner-up `ALI#0`/`ALI#1`
  sits at 72.5% and the next pair at 48.8% (full spectrum in the
  receipt). 10 twin-only holder partitions (224, 351, 872, 929, 2332,
  12263, 13908, 13916, 17761, 19193).
- **Build**: 185 partitions (axis ∪ SK1 holders ∪ twin holders),
  the same export tool and mode; wall 20 s, RSS peak 5.46 GB, 64 GiB
  guard clean. Partition 16 (shared with the chrI build) re-rendered
  BYTE-IDENTICAL — a free determinism proof. Homology dump: 167 window
  queries (padding 120), wall 19 s; walks: 996 truth/twin member rows,
  wall 5 s.
- **The 805-class test at chrIV scale (all 167 loci)**: IN-AXIS 96 /
  SPLIT-ELSEWHERE 69 (NEIGHBOR-AXIS 61 — every holder is an axis
  partition, a one-window coordinate-offset class — / FOREIGN-REPEAT
  8, the partition-spine disease class: L0+L1 subtelomeric 292/544, L39
  2073, L52 2086/3483/14624, L53 2088/3484/14624, L54 2088, L123 117,
  L124 3715) / CONTIG-END-LENGTH-POLYMORPHISM 2 (L163, L166: SK1's
  chrIV ends 80 kb before S288C's — though both windows still carry
  window-query hits onto SK1, e.g. L163's 278-anchor hit tiled in axis
  partitions 285-287). By the ORTHOLOG statement the split is even
  milder: 96/98 assessable loci IN-AXIS, the 2 exceptions (L66, L93)
  in the repeat-locality axis partitions 110/234.
- **The near-twin condensation**: AAA/SGDref share 58,048/58,048
  nodes (100.0%) across all 163 common partitions, ZERO variant
  pockets — identical haplotypes through the graph; CLL#0/CLL#1 share
  54,121 nodes (99.6% of both sides) across 157 common partitions
  with 50/49 left/right variant-pocket bp runs.
- **Joint spellability**: JOINT-IN-AXIS 96/167, NEIGHBOR-AXIS-SPLIT 61,
  FOREIGN-SPLIT 8, CONTIG-END 2 — the partition universe holds BOTH
  homologs' material along the whole aligned length of the chromosome;
  the windowed frame's 0/167 is a frame artifact, not missing material.
- **Id continuity**: all 996 truth/twin member walks spell EXACTLY by
  their GFA P lines (checker phase 5).
- **Receipts** (beside the chrMT/chrI set, same directory):
  `partition-graph-chrIV-{build-list.txt,homology-requests.tsv,
  walk-requests.tsv,homology.jsonl,walks.jsonl,census.jsonl,
  census.tables.txt}`, build/dump markers + RSS pollers
  (`build-chrIV-partition-graphs.*`, `dumps-chrIV-*`), and the
  independent checker log `check-partition-graphs-chrIV.log` — ALL
  PHASES PASS. Scripts: `genome/instrumented/
  partition-expressibility-census-chrIV.py` (requests + census),
  `run-chrIV-partition-graphs.sh`, `check-partition-graphs-chrIV.py`.

Assessment-side only; no thresholds (the 50 bp interval-close window
and the 0.3 spectrum floor are stated classification/reporting
conventions, like the 120 bp homology padding); scoreboard machinery
untouched; no 502; no PR push. THE EXHAUSTIVE chrIV SCORING RUN IS
SLICE 2 — not started here.
