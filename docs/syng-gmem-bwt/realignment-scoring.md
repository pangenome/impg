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

## Slice C: the generated-here/generated-elsewhere marginalization (2026-11-06)

**The defect being fixed.** A locus's record set is defined by the
routing's territory touch, so conserved pockets bring reads into the
locale's evidence that a candidate spelling the pocket collects at full
match score while a candidate that does not pays the all-mismatch floor
(`150*B` ~ -1,546) — even though the read's own genome generated it
somewhere the candidate models only as "outside this locale". The
owner's rho-squared principle: conserved material cancels;
locally-variable material discriminates.

**The anatomy first (the form-deciding gate, measured at chrI L7).**
For every unit in L7's record set, classified by whether the slice-B
winner (the ATV/BAD scaffold pocket folds) and the truth pair place it,
the placement anatomy on the winner and truth folds (the
`--anatomy-out` emission on `partition_realign_score`): the
anchored/skipped/abstained decomposition, the collinearity verdict
(the pinned anchors at ONE offset), the slice-B piecewise (m, c) beside
the one-continuous-path whole-read score, the per-skipped-base votes
with their projected coordinates, and the read's FULL-READ match at
its verified occurrences on the sample's own donor paths
(S288C/SK1, origin-shift-corrected, strand-aware).

* **THE VERDICT: CONTINUOUS.** The winner's pocket placements are
  100% collinear — one continuous path through the pocket row's spell
  — with m mean 148.7/150 (c = 0 on 75.5% of the 1,603 winner-only
  placements, c = 1 on 17%). The flanks vote AGAINST the scaffold
  rows' own context and MATCH it: the pocket copies extend identical
  through the reads' full spans, and no read spans into unique flank
  material within its 150 bp. No placement-continuity rule is
  imposed; the marginalization carries the whole weight.
* **THE ANATOMY'S HEADLINE FINDING: the pocket reads are NOT
  foreign/remote-generated.** Every winner-only-placed unit (1,603
  units, mass 1,677) has a FULL-READ 150/150 match on the sample's own
  donor paths S288C#0#chrI / SK1#0#chrI: 1,051 units (66%) INSIDE
  window 7's own truth rows, 426 (27%) in the neighboring ±2 kb seam,
  only 126 far. Yet the truth folds place none of the in-window ones:
  2,059 units (mass 2,140) at L7 have an in-window full-match donor
  occurrence with NO valid truth placement.
* **THE CAUSE, proven with a named example (record 113 / unit 84 /
  node 92071): the pin skeleton's stored-walk frame blindness.** The
  fold pin skeletons are built from `walk_path_range` — ONE frame's
  syncmer selection — and reads whose anchor k-mers are
  rc-frame-qualified on the truth path are ABSENT from that walk
  (slice A's binding lesson, missed in the slice-B pin skeleton),
  while the same k-mers appear canonical-forward on the
  reverse-complement scaffold copies, which pin and collect the full
  match. The truth's `150*B` floor for these reads is a
  placement-machinery frame artifact, not an inability to spell them.
  The slice-B "54,015 per-base logs" also decompose into
  one-homolog floor asymmetries (±1,546-magnitude mixture terms), not
  per-base evidence — on jointly-placed reads the truth's per-base
  evidence is BETTER (m 146.8 / c 1.08 vs the winner's 142.9 / 2.93).
  The frame repair (canonical-scheme pin skeletons, the territory
  extraction slice A already validated) is the named next lever; the
  marginalization bounds its cost at E in the meantime.

**The marginalization (derived, no tuning constants).** Each read's
likelihood under a candidate marginalizes the two generation branches:

    LL(read | candidate) = logsumexp( local-spell branch , E(read) )

The local-spell branch is rule 2's per-(record, row) likelihood (the
uniform placement prior + the anchored placements). `E(read)` is the
CANDIDATE-INDEPENDENT elsewhere branch: the read explained by the
genome outside this locale at the derived background rate — the read
matches its origin at the MEASURED per-base rate (the census measured
zero unmatched placement bp over 1,063,522 verified occurrences:
reads are exact segments of their origins, at the reads' own
Phred-derived rate A), under the max-entropy uniform origin over the
donor's diploid genome:

    E(read) = len(read) * A - ln(2 * (G_diploid - len(read) + 1))

The same max-entropy spirit as the old Poisson model's beta = M/|U|:
the unit of match likelihood spread uniformly over the universe of
placements. The genome size is measured truth-free from the panel —
the instrument's own genome universe: the mean per-(strain,
haplotype) total path length (14,199,901.5 bp over 235 haplotype
path sets), doubled for the diploid generation model
(G_diploid = 28,399,803 bp; E = -17.870035 on this sample).

* **THE FLOOR for unexplained reads is E**, not the all-mismatch
  150*B (which survives only as the per-placement lower bound inside
  the local branch's logsumexp — a bit-exact no-op wherever an anchor
  pins, since `e^{150*B}` underflows).
* **A read whose local placements score no better than E does not
  vote beyond E** (the mixture absorbs every local branch below E);
  **the logsumexp is monotone in the local branch**, so a candidate
  that spells a read at least as well as every rival under the local
  branch keeps its per-unit ordering after the marginalization — the
  L4-control preservation property (unit-proven).

**The pilot, before (slice B) → after (slice C), all three loci
checker-verified (check-realign-scoring.py --marginal, ALL PHASES
PASS, 6,589 checks, 0 failures; the slice-B mode still passes its own
receipts, 6,757 checks, 0 failures):**

| locus | truth rank before → after | log gap before → after | winner before → after |
|---:|---|---:|---|
| 4 (control) | 1 → 1 | 0.0 → 0.0 | unchanged (bit-exact, checker-enforced) |
| 2 | 4 → 2 | 151,917.18 → 60.74 | unchanged (BTE identical fold + S288C/W303 fold) |
| 7 | 1,644 → 8 | 1,247,492.61 → 982.04 | (ATV+BAD scaffold pockets) → (S288C/AAA/SGDref fold + ATV pocket fold) |

The pocket harvest shrinks 99.92% at L7; the winner stops being the
pure scaffold-pocket pair; the both-placed evidence gap flips negative
(-2,900.7: the truth's evidence is better on jointly-placed reads);
the identical-through-graph fold stays one candidate by construction
(checker-verified from the GFAs); QUAL stays unbounded (the cluster
p-form's relative likelihoods underflow at the surviving separations).

**What remains (next slices).** (1) THE FRAME REPAIR — the
canonical-scheme pin skeletons (slice A's territory extraction) so
that rc-frame-qualified anchors pin on the rows that spell them: at
L7 alone this affects 2,059 in-window full-match units the truth
currently cannot place; it is a slice-B machinery repair, needs its
own exactness gates, and should land before or with the exhaustive
rerun. (2) The exhaustive chrMT/chrI rerun under the marginal model
(all 35 loci) with the before/after table. (3) The QUAL p-form's
dynamic range (unchanged per the ruling; the saturation is named).

## Slice D — the canonical-scheme pin skeleton (the frame repair)

Slice C's anatomy proved the pocket reads are not foreign: 66% have a
FULL-READ 150/150 match inside the window's own truth rows, and the
truth's failure to place them is a PLACEMENT-MACHINERY FRAME ARTIFACT —
the fold pin skeletons were built from the syng's stored path walk
(`walk_path_range`), which carries only ONE frame's syncmer selection,
so anchor k-mers that are rc-frame-qualified on a path are absent from
the walk (slice A's binding lesson, missed in slice B's pin skeleton).
Slice D repairs the skeletons.

**The repair (the routing's own convention, ported from
`build_territory_index_rows`):** a fold's pin skeleton is the
CANONICAL-SCHEME selection over its row extent — per position, the
frame that spells the k-mer's canonical form min(K, rc(K)) forward
decides; the rc frame's qualifying k-mers come from the raw extraction
on the fetched range's reverse complement (mapped back to path
coordinates, the sign carrying the frame's orientation; node identity
is frame-independent); a position whose canonical frame did not qualify
carries NO step. The census records are derived under this same scheme,
so every anchor of every record is present at every true occurrence of
either orientation. Every emitted step is verified BY SEQUENCE at
extraction time against the AGC-fetched panel sequence (`syncmer_seq`:
the interned window at the claimed position, orientation-aware per the
sign — this verification caught the kmerHashSeq lowercase convention on
its first run), the frame decision is verified against the window's
canonical form, and the forward part is verified complete against the
stored walk. The before-record stays reproducible under
`IMPG_REALIGN_STORED_WALK_SKELETON` (the identity gate).

**The diagnosis (all axis-partition rows, both components):** EVERY row
is affected — chrI 4,472 rows (added 678,789 rc-frame-only anchors,
dropped 735,805 stored positions whose canonical frame did not
qualify), chrMT 1,252 rows (added 97,371, dropped 174,719); the
canonical selection is ~50/50 by frame on chrI. The named proof case
(record 113/unit 84, node 92071): the truth row now carries (-92071, rc
frame) at 73,906 where the stored walk has no step (its stored-frame
neighbors 7446960/9705993 are what blinded it).

**The gates (check-realign-scoring.py --frame, ALL PHASES PASS, 579,969
checks, 0 failures):** the identity run reproduces the committed
slice-C receipts on every semantic field with byte-identical sidecars;
the factorized-vs-direct exactness gate re-ran in full (1,405,438/
1,405,438 placements exact, up from 768,827 — the repaired skeletons
place more); the no-regression sweep measured ZERO pin losses at all
three loci (the dropped stored positions never carried a true pin), pin
gains 39,114/297,196/300,301 (L2/L4/L7), pinned-both pairs with
changed ll 3,096/16,967/36,330 (added alternative monotone assignments;
max |delta| 0.0069/7.9709/0.0076, every changed placement per-base
exact); the blinded-unit census (units with a full-match donor
occurrence inside a truth-fold member row and no truth placement, the
paired anatomy receipts): L7 2,003/mass 2,083 → 14/14, L4 2,115 → 33,
L2 2,215 → 280 (the L7 remainders are edge-discard/ambiguity cases —
anchors present and pinned but the read overhangs the row).

**The pilot (chrI 2,4,7, before = the slice-C marginal receipts = the
identity gate):**

| locus | truth rank before → after | log gap before → after | winner before → after |
|---:|---|---:|---|
| 4 (control) | 1 → 1 | 0.0 → 0.0 | unchanged (bit-exact verdict; the winner IS the truth pair, checker-enforced) |
| 2 | 2 → 4 | 60.74 → 273.96 | unchanged material (BTE identical fold + S288C/W303 fold) |
| 7 | 8 → **1** | 982.04 → **0.0** | pocket pair → **the truth pair itself** (bit-exact winner) |

At L7 the winner-only class collapses 1,198 units → 0 (the 1,051-unit
winner-exclusive fuel fully re-absorbed: every unit the winner places,
the truth also places) and the tail locus CLOSES. At L2 the repair is
honest UNBIASED machinery — the rival BTE scaffold row gains its own
rc-frame anchors and collects 105 → 225 units of seam-straddler mass
(the row-geometry lever named since slice B: the BTE row extends 556 bp
past SK1's row end), so the truth rank returns to 4 with a
truth-favoring both-placed gap (-207.3). The L4 control holds bit-exact
(rank 1, gap 0.0, the winner the truth pair, winner material unchanged).

**What remains (the next slice):** THE EXHAUSTIVE chrMT/chrI RERUN
under the marginal model with the repaired skeletons (all 35 loci) plus
the before/after table — the instrument is now frame-unblinded and the
rerun is the delivery. The L2 seam story (row extents) and the QUAL
p-form's dynamic range stay named as before.
