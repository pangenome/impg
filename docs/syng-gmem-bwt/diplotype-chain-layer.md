# The chain/phasing layer: two whole-chromosome molecules per component

The owner's product ruling from the campaign's start: **the diploid
product is two molecules**; per-locus genotyping is the product, and
the chain is a **thin phasing layer over it**. The CLI product
(`impg genome-infer call-diplotypes`) emits the per-locus diplotype
calls (the two called paths, the class identity, the model-internal
QUAL, E, the alternative log-gap). The chain layer
(`impg genome-infer chain-molecules`) links those CLOSED calls into
the two molecules. This document is the layer's design of record.

## The thin-layer property

The chain consumes `calls.jsonl` read-only and NEVER writes it: the
per-locus calls are bit-identical pre/post chain by construction.
The gate proves it three ways: the consumed-calls FNV fingerprint
emitted inside the product record equals the fingerprint of
`calls.jsonl` after the run; every molecule locus entry is EXACTLY
the emitted called pair (fold indices and member rows bit-equal),
reoriented only; and the per-locus accuracy columns of the switch
table are carried UNCHANGED from the product run's own truth-QV
artifact. No locus call is changed anywhere — the chain only links.

## The chain: a two-state orientation DP

The called pair per locus is UNORDERED. A chain state
`s_j ∈ {0,1}` says which emitted fold of locus j rides molecule 1.
Each boundary between consecutive called loci carries two hypotheses —
**same** (`s_{j+1} = s_j`) and **swap** (`s_{j+1} = 1 − s_j`) — with
measured evidence for both, in a stated lexicographic order
(no thresholds, no tuning constants):

1. **The edge layer** — the junction-spanning reads' votes (the
   layer COSIGT never had). A census occurrence touching BOTH
   consecutive windows is a verified-placement candidate of a read
   pattern spanning the junction. Its covered panel material (the
   census's contained segment-node windows on its own path) is
   attributed to the called folds' member rows; a record votes ONCE
   per boundary (the record-once discipline: its census multiplicity
   read instances are its mass), split uniformly over its crossing
   occurrences at the boundary and over the attributed fold pairs.
   Both-side attribution is required: one-sided material is **orphan
   crossing mass** — counted and reported, never a vote, because the
   called pairs do not express that molecule's continuity there.
   Occurrences whose census node lists are empty are the placements
   the scoring instrument's own binding would reject (no panel step
   at the anchor); the chain excludes them consistently.
2. **The adjacency layer** — the called folds' shared-segment
   continuation: member rows on the SAME panel path that continue
   across the junction (exact touches `row_f.end == row_g.start`
   dominate, then shared-overlap bp).
3. **The unidentifiable tie** — both channels exactly tied (a
   homozygous emitted pair, or symmetric evidence): the boundary is
   flagged `unidentifiable` and the stated deterministic rule carries
   the chain through ("same"). The chain may not invent certainty.

With only relative orientation evidence, the 2-state Viterbi
maximizes each boundary independently; the accumulation
`s_{j+1} = s_j XOR (boundary = swap)` builds the two molecules.

### How ambiguous loci enter the linkage

A locus with tied classes (`called_class_count > 1`) enters through
its EMITTED winner pair; every boundary at it is flagged
`through_ambiguous_locus` (an alternative tied pair could reorient
it). A homozygous emitted pair makes both adjacent boundaries
structurally unidentifiable. Every link confidence carries the two
endpoint QUALs (the model-internal per-locus call certainty
weighting the uncertain links) and the endpoint flags — never a
synthesized posterior.

## The emission (molecules.jsonl)

One record per component: the two molecules — the per-locus
orientation bits, the per-locus fold material with junction
accounting (the axis-window extent vs the fold's own material;
negative gaps are honest overlaps, positive gaps honest seams), and
the spelled sequences — plus every boundary's full evidence
(vote masses per link, adjacency both ways, orphan and unattributed
crossing counts, the window gap), the orientation, the decision
channel, and the link confidence.

## The test mode (assessment-side only)

`--truth-qv-file` consumes the product run's OWN
`calls.jsonl.truth-qv.jsonl` (never the default output) and emits
`molecules.jsonl.truth-qv.jsonl`: per boundary, the chain
orientation vs the truth orientation under the committed assignment
yardstick, with **recomputed rate-and-edit ties bracketing the
boundary** (the truth side's own unidentifiability — a homozygous
call, an inexpressible truth pair — never silently resolved), and
per-locus accuracy under best-case assignment carried UNCHANGED in
separate columns. **Switch errors and per-locus accuracy are never
conflated** — the owner's evaluation ruling since the haploid era.

## The gate (chain-molecules-gate.py)

Per component: the thin-layer proof (the calls fingerprint, the
molecule-loci bit-equality with the emitted pairs, the QUAL
carriage), the orientation self-consistency (the bits accumulate
the boundary orientations; every boundary is the stated
lexicographic decision of its own emitted evidence), the switch
summary arithmetic, and the rendered switch table.

## The measured gates

- **chrMT**: 14 loci, 13 boundaries, 8 assessable, **0 switches**
  (5 non-expressible brackets; decided crossing 8 / adjacency 4 /
  tie 1 — the homozygous locus 0). Gate 176/176 checks.
- **chrI**: 21 loci, 20 boundaries, 15 assessable, **2 switches**
  (3 non-expressible + 2 truth-tied brackets; crossing 14 /
  adjacency 1 / tie 5). Gate 260/260 checks. Both switch errors sit
  at UNIDENTIFIABLE boundaries adjacent to locus 16, a non-rank-1
  call (error 0.29) whose emitted pair order flips vs truth — the
  chain had no evidence there, carried the stated tie rule, and the
  switch table says so with the numbers.
