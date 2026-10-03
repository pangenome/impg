# The realignment scoring layer (slice B of the local-realignment evidence model) — note of record

Owner direction 2026-11-06; the anchor-projection layer (slice A,
`partition_anchor_projection.rs`) stands as the placement machinery of
record; this layer is the scoring built on it. Assessment-side only:
`examples/partition_realign_score.rs` + `genome/instrumented/run-realign-score.sh`
+ `genome/instrumented/check-realign-scoring.py`. No `src/` file touched,
no threshold, no tuning constant, scoreboard machinery unmodified.

## The model, derived and stated

**The scoring rate** comes from the reads' own base qualities: the
sample's FASTQ carries one uniform quality byte (`I` = Phred 40,
asserted by streaming the whole file), so the per-base error
probability is `eps = 10^(-40/10) = 1e-4` and the substitution model
gives `P(match) = 1 - eps`, `P(mismatch) = eps/3`. The per-base log
weights are `A = ln(1 - eps)`, `B = ln(eps/3)`. The measured error
profile is stated beside them: slice A measured zero unmatched
placement bp over 1,063,522 verified occurrences.

**Per (record, candidate row) skipped-base votes**, computed in-process
by slice A's committed machinery (the unique-monotone anchor
correspondence, the serving-anchor projection, the mirror equations —
the same functions the slice A receipts and checker validated; the
tens-of-GB full-fidelity per-path context is computed, never
materialized). A candidate row's log-likelihood for one read instance
is

```
LL(read | row) = -ln(2 * (len_row - 149))          # the uniform placement
                                                  # prior, both strands
                  + logsumexp over anchored placements of
                    [ backbone + skipped-base votes ]
                  with the derived all-mismatch FLOOR term
```

* **The backbone** is the merged bp of the PINNED anchor windows (the
  anchors the row's contained walk shares with the required sign):
  matches by node identity and orientation sign, identical across
  near-identical paths (rows sharing the anchor chain share the
  backbone), computed once per (record, placement), asserted per base
  against the row's sequence by the direct check over the whole scored
  domain.
* **The skipped-base votes** are the read's bases outside its merged
  anchor windows (the flanks; interior hull gaps measured zero, emitted
  honestly per mirror state when present): each skipped base is
  projected onto the row through its serving anchor (the bracketing
  pinned anchor, slice A's committed rule) and votes `A` on a match,
  `B` on a mismatch against the row's spelled sequence, complemented at
  physically-reverse placements.
* **The placement prior** is uniform over the row's full-read
  placements: `2 * (len_row - 149)` of them — the derived count.
* **The floor**: a read the row cannot anchor scores at the model's
  per-read MINIMUM, `150 * B` (every base mismatching — the rigorous
  lower bound of the read's likelihood under any placement), inside the
  same logsumexp. As a term it underflows to +0.0 exactly wherever an
  anchored placement exists, so it is a bit-exact no-op there. The
  floor's magnitude is NOT the true unanchored likelihood (the best
  unanchored alignment would sit near `-950`, not `-1546`); ranks are
  floor-magnitude-robust (the floor enters every class through the
  same constant), gap magnitudes are floor-dominated where unplaced
  counts differ — the `log_gap_decomposition` field carries the
  attribution (unplaced mass on each side and the both-placed
  per-base evidence gap).
* **Validity**: a placement exists only when every voted base
  (backbone windows and skipped votes) projects inside the row's
  spelled extent — the read lies fully inside the row, the uniform
  prior's own support. Edge-discarded placements are counted per locus.
* **The named remainder**: bases inside anchor windows the row does
  not pin (absent or ambiguous anchors, and edge-overlapping windows —
  a step whose 63bp window is not fully inside the row extent cannot
  be a backbone anchor) abstain. Counted per locus as
  `abstained_hull_bp`; measured negligible at the pilot (3.1 bp per
  scored placement at chrI L7).

**Record-once, no share-splitting**: each candidate row is an
independent hypothesis; the record votes once per row — multiple
occurrences of one record pin the SAME placement on a row (the pin is
a property of the pattern, the orientation and the row) and count
once. The record's instances are its evidence: every distinct read
variant (mirror state, hull start, sequence) votes with its own
instance count, and the class mixture is per instance.

**The locality ruling** (measured before implementation): the locus's
evidence set is the records with at least one occurrence TOUCHING the
locus's window (the routing's own territory-touch rule — the same
convention that defines the current instrument's per-locus record
sets), and the occurrences that vote at the locus are exactly those
touching occurrences. The measured outside class: ALL 2,142 chrMT +
9,779 chrI occurrences outside built axis-partition rows are
own-row-edge overhangs — the read's anchor hull crosses the
alignment-induced membership boundary of its own path (the partition
BED's chain end falls inside the read's hull; e.g. ADQ#0#chrI's
partition-7 row ends at 71,751 while a 63bp read sits at
71,744..71,807); ZERO occurrences are remote or foreign. Ruling: IN —
structurally, the membership row boundary is the window-seam artifact
the partition spine was built to dissolve, the likelihood a read
supports a row is a property of (read, row sequence), and the
placement-validity rule already restricts such reads to testify only
where the row spells their bases. At the pilot windows the class is
156/9,530 (chrI L2), 1,184/75,419 (L4), 688/90,206 (L7) of touching
occurrences.

**Candidates** fold by identical sequence AND contained walk —
coalesced copies are ONE hypothesis (the identical-through-graph pair
BTE#3/#4#block28_contig1 folds by construction; the checker verifies
the two member rows spell identical sequence and walk from the
committed GFAs). Classes are all unordered fold pairs `i <= j`; the
class log-likelihood is the sum over the locus's evidence units
(record, variant) of `count * ln(0.5 e^{LL_i} + 0.5 e^{LL_j})` — the
equal 1/2 mixture is the balanced 15x/15x diploid prior. QUAL is the
existing cluster machinery verbatim (`spectrum_knee`,
`cluster_form_qual`, `qual_p`, `qual_from_p` and the material-distance
spectrum over the likelihood ratios — mirrored pure functions from
`panel_route_spine/cosine_probe.rs`), with the material distance the
differing observed mass on merged node+edge usage multisets, observed
mass = record multiplicity (no share split).

## The exactness proof

* **In-process, the whole scored domain**: for EVERY scored (unit,
  fold, placement) the factorized `(m, c)` is asserted equal to the
  direct per-base comparison — the backbone verified base-by-base
  against the row's sequence (not trusted from node identity), over
  the MERGED pinned-window coverage (syncmer windows overlap heavily;
  counting per anchor would double-count the overlaps), and the
  skipped votes recomputed. Integer equality, loud failure. The pilot
  measured 768,827/768,827 EXACT. Two real defects were found and
  fixed on the way, both unit-proven: the direct backbone double-count
  (overlapping pinned windows), and the direct comparison direction at
  physically-reverse placements (the read's window is the row window
  REVERSED — offset `t` corresponds to row offset `k-1-t`).
* **The independent checker** (`check-realign-scoring.py`, ALL PHASES
  PASS, 6,637 checks, 0 failures) re-derives everything from the
  committed artifacts: the folds' sequences and contained walks from
  the GFAs (P/L/S spelling at the front-overhang offset, gap segments
  excluded), the locality counts, and the exactness sample — 773
  (record, fold) placements re-derived EXACTLY (114,382 per-base
  votes and backbone bases; every sampled read verified present in
  the FASTQ; the correspondence re-enumerated from the GFA walk), plus
  the bounded-offset dominance scan (the anchored placement attains
  the best single-offset substitution score within +/-150: 646/770
  collinear placements dominate; the 124 named exceptions carry FLANK
  indels — no single offset can collect both sides of a flank indel,
  the piecewise serving projection does; the stated
  shifted-mismatch rule) and the bounded edit-distance alignment
  (646 indel-free, 124 with flank indels).
* The class table is re-derived in full from the per-unit matrix
  (all 40,169 classes across the three loci; max deviation 1.2e-6),
  the truth ranks reproduced, the floor entries verified at
  `prior + 150B` exactly, the QUAL cluster states reproduced.

## The pilot receipts (chrI, exit 0, wall 61s, RSS peak 2.6 GB)

| locus | folds | units | classes | truth rank | log-gap | QUAL | before (read-matched instrument) |
|---|---|---|---|---|---|---|---|
| 2 (identical-pair locus) | 31 | 6,987 | 496 | 4 | 151,917.18 | unbounded | truth NOT EXPRESSIBLE (FAILED-AS-BRACKETED) |
| 4 (truth-rank1 control) | 190 | 4,576 | 18,145 | **1** | 0.0 | unbounded | rank 1, QUAL 109.51 |
| 7 (the tail locus) | 207 | 5,851 | 21,528 | 1,644 | 1,247,492.61 | unbounded | rank 1,086, gap 3,105.9 |

* **chrI L4 — the control holds**: the truth pair (S288C#0#chrI,
  SK1#0#chrI folds) is the bit-exact winner; rank 1, gap 0.0.
* **chrI L2 — the bracket opens**: the truth pair is now EXPRESSIBLE
  (rank 4; the windowed instrument could not express it at all). The
  winner is (the BTE#3/#4 block28_contig1 identical fold, the
  S288C/W303 identical fold): the BTE scaffold pair replaces SK1's
  homolog. The gap decomposes as 105 read-mass of unplaced
  seam-straddlers (the BTE row extends 556 bp past SK1's row end and
  places them; the truth pair floors them) plus 3,977 logs of per-base
  evidence on jointly placed reads.
* **chrI L7 — the tail widens**: truth rank 1,644 (was 1,086); the
  winner is (ATV#3/#4 block43_contig1, BAD#3/#4 block96_contig1) —
  scaffold pockets. The gap decomposes as 849 read-mass of floor (the
  truth rows cannot place the pocket reads at all) plus 54,015 logs
  of per-base evidence (the pocket rows spell the pocket reads
  exactly — the multi-census stage's foreign-scaffold-mass finding
  CONFIRMED in likelihood units: the tail is not an evidence-model
  artifact; the lever remains the universe/row structure).
* **The identical-through-graph pair stays tied — correctly**: 5397/
  5545 fold to ONE candidate at L2 by identical sequence+walk
  (verified from the GFAs); there are no separate candidates to tie.
* **QUAL saturates unbounded at all three loci**: the cluster p-form
  consumes RELATIVE likelihoods (`e^{ll - best}`), and at
  realignment-scale separations (10^4+ logs between the winner and the
  nearest divergent rival) they underflow to +0.0, giving `p = 1`. The
  machinery is unmodified per the ruling; the saturation is the
  honest behavior of likelihood ratios at these magnitudes (with the
  true unanchored floor the separations would still underflow: ~600
  logs per unplaced read).

## What remains (next slice)

The full exhaustive chrMT/chrI rerun under the realignment likelihood
(all 35 loci), the before/after table against the committed
read-matched receipts, and the prediction table (the tail loci, the
controls, the near-twins). The named refinements the pilot exposed for
the owner: (1) the unanchored FLOOR's magnitude (the derived minimum
vs the true unanchored likelihood — ranks are robust, gap magnitudes
are not); (2) the QUAL p-form's dynamic range under likelihood ratios;
(3) the within-anchor variant positions (the abstained hull
remainder, measured negligible at the pilot); (4) the row-extent seam
effect (rows extending past a neighbor's end collect the
straddler reads — a row-geometry lever, not an evidence lever).
