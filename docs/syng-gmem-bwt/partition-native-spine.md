# The partition/window duality — the architecture note of record

Date: 2026-09-29. Branch `work/genome-mem-bwt-pipeline`. This note records the
partition-vs-window perspective as the architecture's remaining structural
debt, the seam's measured costs from the 2026-09-28/29 arc, the path to a
partition-native spine, and the timing recommendations the owner has marked
planning-relevant. It is a docs note of record — no code changes ride with it.
The working details live in the per-arc reports (`genome/` artifacts, the
`context-and-mem-likelihood.md` running report); this note preserves the
PERSPECTIVE against session compaction and worker turnover.

## 1. The two frames

The phasing DP's spine is REFERENCE WINDOWS: ~10 kb intervals laid on the
reference axis in reference order — the bootstrap's scaffold. The material
the DP assigns to those windows is PARTITION-TERRITORY ROWS: per-source
frames, grouped by the graph's own partition territory (each source's own
contiguous DNA, the panel's physical structure). A window is a chunk of the
AXIS; a row is a chunk of a SOURCE. The window carries the observed evidence
(the reads' observed universe is pooled per window); the territory rows carry
the expressible truth (the vocabulary the sample's material can be spelled
from). The DP stitches rows across window boundaries — every such stitch is a
FRAME SEAM, and the seam is where this arc's defects lived.

The duality is the architecture's remaining structural debt: the evidence
frame (windows) and the expression frame (territories) are not the same
geometry, and every coupling between them pays the seam.

## 2. The seam's measured costs (this arc's evidence)

Every ITEM-2 defect class lived in the frame seam:

- **Frame-alignment failures on novel junctions.** A junction joins two
  INDEPENDENT frames — co-linearity applies only to same-source pairs. The
  truth's own 9564→9602 junction failed the old (pre-fix) co-linearity test
  by 10,959 bp; the fix (same-source pairs only) is the seam's first measured
  repair.
- **The coordinate-overlap staging rejecting attested stranded rows.** The
  loop-closure fix (locus-range scoping of the junction span index) repaired
  the window-anchor projection that stranded attested rows outside their
  windows.
- **The truth-piece vocabulary mismatch.** The truth's local material read
  through the window's row vocabulary required the traversal-locality filter,
  the existential in-domain predicate, and the segment-major derivation
  order fixes — all seam corrections at the checker side.
- **The capstone: the window's observed universe pools cross-source mass.**
  At chrIII locus 10 the window's observed universe measured 7,077 distinct
  feature keys / mass 28,289 while the TRUTH's own row explains 1.03% of it
  (the best rows explain 17.4%) — the window pools both copies' reads plus
  background at a diploid window, and no single row can explain it. The
  omission term dominates the truth row's loss 72x (self 1,315 vs omission
  94,179; rank 708/1210 at its own window). THIS IS ITEM 3'S GATE EVIDENCE:
  the window-mass instrument's objective-vs-truth anti-correlation, measured
  at full-chromosome scale.

The same instrument limitation was then measured on the DIPLOID axis (the
2026-09-29 smoke): a true 6x/4x S288C+SK1 sample selects HAPLOID on both
chrMT and chrI (the true pair priced 30k-165k nats worse than the single),
0 heterozygous pairs called, dosage columns reading 2.0-on-one-class against
the true 1.2/0.8 — the pooled evaluator preferring the truth on both. Two
independent gate evidences for the same structural fact: the window's pooled
mass cannot distinguish material it cannot frame.

## 3. The path

- **ITEM 3 de-weights the window.** The MEM-projection local likelihood:
  per-record positional evidence against the pair's spelled MEM spans.
  Source-specific records place on their true frames and credit them
  directly; the window becomes the DP's CHUNKING GRAIN, not the evidence's
  frame. No constants; Poisson bones; the S2 record-level machinery
  (InstanceStructure) as foundation.
- **The endgame is the owner's standing directive: de novo out of the
  graph.** The partition/block graph's OWN ADJACENCY as the DP's spine —
  territory boundaries, not reference windows, order the chain; the
  reference axis is demoted to diagnostics and scaffolding. The seam
  disappears when the spine and the material share the frame.

## 4. What this arc already built for it

The components of a partition-native spine are already in the tree:

- **The territory tables** — frame-agnostic: every partition's rows for
  every source, the row vocabulary the DP assigns.
- **The attested census** — sample-driven, not reference-driven: per-read
  chain-derived frame attribution, the frame-offset diagonal (rp = lp + F),
  the aggregate AttestedComposition form.
- **Maximal rows and partial chains** — the stitch machinery's maximal-row
  filter and the partial-chain form (a deletion is a real allele).
- **The honest-untypable brackets** — the copy-choice seam brackets with
  read-length falsifiable predictions (≥299 bp / ≥1,359 bp at 150 bp reads).
- **ITEM 4's END markers** — the continuation-exhaustion census against the
  panel's path termini: the intra-window analog that anchors the spine's
  own termini.

Each was gated and measured in the window frame; each carries over as the
spine's vocabulary when the frame changes.

## 5. Timing (planning-relevant)

- **The window-spine machinery is working and fast.** The whole-genome
  rerun (2026-09-29): 17/17 components, aggregate acc_H 0.9611 / 1 switch,
  the full genome in ~2.6 h serial (the largest component 63 min; the old
  4-core era's chrVII alone was 3.1 h). Iterate on it.
- **ITEM 3 next.** The window stops being the evidence frame; the two gate
  evidences (haploid chrIII, diploid smoke) name the defect it repairs.
- **Then the partition-spine migration as a planned iteration** with the
  window-frame baseline in hand.
- **The diploid validation runs BEFORE the migration unless the owner
  reorders** — measuring the scaffold first is the default; this note
  records the owner's standing option to reorder. (The diploid smoke
  completed 2026-09-29; the full validation is the owner's iteration
  campaign with real numbers in hand.)
