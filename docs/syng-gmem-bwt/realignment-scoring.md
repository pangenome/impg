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

## Slice E — the exhaustive chrMT/chrI rerun under the marginal model with repaired skeletons (2026-11-06)

The delivery slice: every window of both components — all 14 chrMT + 21
chrI = **35 loci** — scored exhaustively under the frame-unblinded
instrument (`anchor-realign-v2-marginal-frame`: the marginal realignment
likelihood, the canonical-scheme pin skeletons, the existing cluster
QUAL machinery verbatim). No Rust source changed this slice; the
instrument is the committed slice-D binary semantics, re-run over the
full domains. Receipts `realign-exhaustive-{chrMT,chrI}.jsonl` (+ exactness,
ingredients, records, skeleton sidecars) and run markers
`run-realignexhaustive-{chrMT,chrI}-*` (exit 0; **walls 130 s + 472 s =
602 s serial**; **RSS peaks 2,168,796 kB + 5,949,056 kB = 2.07 + 5.66
GiB, both < 11× under the 64 GiB guard**). The in-process
factorized-vs-direct gate re-verified EVERY placement over both full
domains: **1,929,497/1,929,497 (chrMT) + 21,490,054/21,490,054 (chrI) =
23,419,551/23,419,551 EXACT**.

### The full per-locus table (before = the committed Poisson-era read-matched receipts, the instrument of record)

chrMT (expressible = L2..L10; rank/gap n/a = truth pair not expressible
under either instrument — L0, L1, L11, L12, L13 on both sides):

| locus | rank before → after | log gap after | QUAL before → after | in called set | winner (after) |
|---:|---|---:|---|---|---|
| 0 | n/a | n/a | 9.50 → unbounded | n/a | DG1768+S288C [0,9041) pair |
| 1 | n/a | n/a | 21.79 → 23.21 | n/a | BAL_1a [6913,12888) + DG1768/S288C/W303 fold |
| 2 | 3 → **48** | 3,213.3 | 60.76 → unbounded | no | CBK+CBM [10519,17115) + SK1 [8858,15598) |
| 3 | 1 → **1** | 0.0 | unbounded → unbounded | **yes** | **the truth pair** (DG1768/S288C + SK1) |
| 4 | 15 → **1** | 0.0 | 14.11 → unbounded | **yes** | **the truth pair** (AAA/DG1768/S288C/SGDref + CBK/CBM/SK1 fold) |
| 5 | 1 → **1** | 0.0 | unbounded → unbounded | **yes** | **the truth pair** (AAA/S288C/SGDref + SK1) |
| 6 | 173 → **117** | 8,290.8 | unbounded → unbounded | no | DBVPG6044 [25902,36857) + S288C [33428,39044) |
| 7 | 4 → **1** | 0.0 | unbounded → unbounded | **yes** | **the truth pair** (S288C + SK1) |
| 8 | 2 → **1** | 0.0 | unbounded → unbounded | **yes** | **the truth pair** (DBVPG6044+SK1 fold + S288C) |
| 9 | 3 → 3 | 2,072.1 | 156.54 → unbounded | no | AAA+SGDref [0,10902) + SK1 [58861,69265) |
| 10 | 3 → **119** | 8,030.3 | unbounded → unbounded | no | DBVPG6044+SK1 fold + W303 [58915,65826) |
| 11–13 | n/a | n/a | 140.41/unbounded/127.81 → unbounded | n/a | foreign/scaffold rows (non-expressible) |

chrI (expressible = L2..L19 after — **L2 newly expressible**; L0, L1,
L20 non-expressible on both sides):

| locus | rank before → after | log gap after | QUAL before → after | in called set | winner (after) |
|---:|---|---:|---|---|---|
| 2 | n/a → **4** | 274.0 | 153.53 → unbounded | no | BTE#3/#4 block28 identical fold + S288C/W303 fold |
| 3 | 13 → **2** | 12.9 | unbounded → unbounded | no | AAA+SGDref fold + SK1 fold |
| 4 | 1 → **1** | 0.0 | 109.51 → unbounded | **yes** | **the truth pair** (S288C + SK1) — bit-exact control |
| 5 | 1 → **1** | 0.0 | 57.10 → unbounded | **yes** | **the truth pair** (6-member S288C fold + SK1) |
| 6 | 1 → **4** | 2,712.2 | unbounded → unbounded | no | CBM#1#chrI_2 [0,11799) + SK1 |
| 7 | 1,086 → **1** | 0.0 | unbounded → unbounded | **yes** | **the truth pair** (AAA/S288C/SGDref + SK1) |
| 8 | 1,123 → **1** | 0.0 | unbounded → unbounded | **yes** | **the truth pair** (7-member S288C fold + SK1) |
| 9 | 259 → 4 | 882.2 | 151.76 → unbounded | no | S288C 6-member fold + CFF#1#chrI_1 [82727,95172) |
| 10 | 47 → **1** | 0.0 | 15.50 → unbounded | **yes** | **the truth pair** (S288C + SK1) |
| 11 | 43 → **1** | 0.0 | 18.87 → unbounded | **yes** | **the truth pair** (6-member S288C fold + SK1) |
| 12 | 1 → **2** | 313.2 | 76.37 → unbounded | no | CFF#2#chrI_1 [105350,116305) + SK1 |
| 13 | 3,487 → 1,082 | 35,914.9 | unbounded → unbounded | no | AGK [122220,132514) + CBM#2#chrI_3 [20643,31291) |
| 14 | 61 → 12 | 1,984.8 | unbounded → unbounded | no | BTE#3/#4 block23 identical fold + SK1 |
| 15 | 8,975 → 396 | 26,374.3 | unbounded → unbounded | no | ATV#3/#4 block146 identical fold + SK1 |
| 16 | 18 → **6,613** | 62,913.7 | unbounded → unbounded | no | AIG#0#chrIV [1162863,1168638) + ASB#2#chrIV — cross-chromosome |
| 17 | 3 → **2** | 488.5 | unbounded → unbounded | no | BEM#3/#4 block92 identical fold + SK1 |
| 18 | 1,307 → 621 | 35,958.8 | unbounded → unbounded | no | AIS#1#chrI [167305,177537) + BTE#3/#4 block1 identical fold |
| 19 | 108 → 182 | 19,779.7 | unbounded → unbounded | no | S288C 4-member fold + BFH#1#chrI [149862,160593) |

### The aggregate

Truth rank-1 under the realignment instrument vs the Poisson-era
read-matched receipts: **chrMT 5 of 9 expressible (L3, L4, L5, L7, L8)
vs 2 of 9 (L3, L5)** — the "2/8" convention used the older mask in
which chrMT L4 was bracketed; the receipts measure L4 expressible at
Poisson rank 15 and it CLOSES to rank 1 here. **chrI 6 of 18 expressible
(L4, L5, L7, L8, L10, L11) vs 4 of 17 (L4, L5, L6, L12)** — chrI L2 the
frame-gap-era bracket locus is NEWLY EXPRESSIBLE (rank 4). No locus
expressible under the Poisson instrument is inexpressible under the
realignment instrument (checker-enforced).

### The prediction table, scored honestly

- **The tail loci:** chrI L7 CLOSES to rank 1 (gap 0.0; expected from
  the slice-D pilot); **L8 CLOSES to rank 1** (1,123 → 1); **L10 and
  L11 CLOSE to rank 1** (47/43 → 1); L3 and L17 reach rank 2; L9 rank
  4; L14 61 → 12. **L13 (1,082), L15 (396) and L18 (621) do NOT close** —
  their gaps are unplaced-mass dominated (+19,298 / +12,095 / +17,316
  truth-unplaced read mass: scaffold/conserved rows the winners place
  and the truth rows cannot) with truth-favoring both-placed evidence
  (−8,063.6 / −3,567.4 / −6,529.6). **L16 WORSENS 18 → 6,613** and is
  the one locus where the both-placed evidence itself is winner-favoring
  (+17,119): the winner is CROSS-CHROMOSOME chrIV repeat-family material
  (AIG + ASB rows), 53,746 truth-unplaced mass — genuinely-matched
  repeat material, not an evidence-model artifact. chrMT L2 WORSENS
  3 → 48 and L10 WORSENS 3 → 119 (the winners' rows are 3×/2.7× the
  truth row extents); L6 improves 173 → 117 but stays deep.
- **The old wins:** chrMT L3 and L5 HOLD (bit-exact truth-pair
  winners); chrI L4 and L5 HOLD (L4 the enforced bit-exact control).
  **chrI L6 and L12 are LOST — the honest regressions of this rerun**
  (L6 rank 1 → 4, gap 2,712.2; L12 rank 1 → 2, gap 313.2). Both are the
  row-geometry/seam class: the winner replaces one truth homolog with a
  LONGER row that places more of the window's reads (L6: the
  CBM#1#chrI_2 scaffold row, len 11,799 vs the truth fold's 10,039,
  unplaced-mass asymmetry +812, both-placed evidence −2,758.9
  TRUTH-favoring; L12: the CFF#2#chrI_1 row len 10,955 vs 10,000,
  asymmetry +295, both-placed −1,465.2 truth-favoring). The Poisson
  instrument's rank-1 there was mass-share luck of the seed-matched
  universe; under read-level likelihood the longer row wins on placed
  mass while losing the per-base evidence — the lever is row extent /
  the partition boundary geometry, exactly the L2 class named since
  slice B.
- **The identical-through-graph pairs stay tied — correctly, by
  construction** (checker phase 7 re-proves from the GFAs: chrI L2's
  BTE#3/#4 block28_contig1 fold members spell identical sequence and
  walk; the L14/L15/L17/L18 winners carry their own identical
  scaffold folds as single hypotheses).
- **The frame-gap-era non-expressible loci now score:** chrI L2 — the
  named bracket locus — is expressible and scores at rank 4. Every
  locus the partition-graph expressibility census called jointly
  spellable is expressible here; the remaining non-expressible sets are
  unchanged (chrMT L0/L1/L11/L12/L13, chrI L0/L1/L20).
- **The QUAL story:** under likelihood-ratio separations the cluster
  p-form's relative likelihoods underflow (e^{ll−best} → 0) so p = 1 and
  QUAL saturates UNBOUNDED at 34 of the 35 loci (the honest behavior of
  the unmodified machinery, named since slice B); the single finite
  value is chrMT L1's 23.21 (a 4-fold, 10-class domain). The Poisson
  era's finite QUALs (0.37–156.5) all saturate; the QUAL signal must
  come from a form with more dynamic range if it is to rank confidence
  at realignment separations — the machinery itself stays untouched per
  the ruling.

### The residuals, named per locus with numbers

Three classes cover every non-rank-1 expressible locus:
1. **Row-geometry/seam (the L2 class):** the rival row is longer or
   shifted and places reads the truth rows cannot; the truth WINS the
   jointly-placed per-base evidence. chrI L2 (+225 unplaced mass,
   −207.3 both-placed), L6 (+812, −2,758.9), L12 (+295, −1,465.2),
   L17 (+856, −3,726.0), L9 (+554, −1,585.4), L3 (+2, −1.7), chrMT
   L9 (+457, −172.4). Levers: row extents and window-seam geometry.
2. **Truth-row fragmentation:** the alignment-induced partition
   boundary cuts the truth homolog's material into short/multi-member
   folds and reads over the cut cannot place. chrMT L2 (the truth
   S288C co-fold spans only 2,231 bp beside the winner's 6,596 bp row),
   chrMT L10 (truth S288C co-fold 2,527 bp vs the winner's 6,911 bp
   W303 row), chrI L19 (truth SK1 row cut to 6,634 bp; the winner's
   second fold is a foreign-window row BFH#1#chrI [149862,160593)),
   chrI L13/L18 partially.
3. **Genuinely-matched foreign/scaffold pocket mass (the multi-census
   finding, now in likelihood units):** chrI L13 (+19,298 unplaced
   mass), L15 (+12,095), L18 (+17,316), L16 (+12,511 and the only
   winner-favoring both-placed evidence, +17,119 — cross-chromosome
   chrIV repeat material), chrMT L6 (+1,857), L10 (+1,312). The pocket
   rows spell the pocket reads exactly; no evidence-model repair can
   close these — the lever is the candidate universe/row structure.
   chrI L14's winner (the BTE block23 identical fold) also collects
   +886 unplaced mass over a 12,580 bp row vs the truth fold's 10,027.

### The independent checker, all modes

`check-realign-scoring.py` gains the `--exhaustive chrMT|chrI` mode
(phases 1–8 over every locus: markers/exit/wall/64 GiB guard; every
fold's canonical skeleton re-verified against the partition GFAs with
every frame decision checked — chrMT 1,099 rows / chrI 3,545 rows,
added 89,318 + 536,927 and dropped 160,465 + 584,971 steps, 228,645 +
830,785 byte-level sequence-verified against the S segments; the
locality counts re-derived exactly; the sampled exactness re-derivation
— 3,360 + 5,178 placements, every vote and backbone base; the E
derivation audit and the ≥E bounding property over the whole matrices;
ALL class log-likelihoods recomputed from the matrices — 1,169,912 +
4,341,818 checks total, max deviations ≤ 1.2e-7; the QUAL states
reproduced; the generalized identical-fold proof from the GFAs; the
before/after verdicts with the old-wins PREDICTION VERDICTS — a lost
old-win is reported as a NAMED MEASURED FINDING, never a receipt
failure — and the chrI L4 control hard-gated and passing).
**ALL PHASES PASS in every mode:** the new `--exhaustive chrMT` (1,169,912
checks) and `--exhaustive chrI` (4,341,818 checks; the two prediction
violations at L6/L12 named in the log), and the three prior modes
unchanged on their committed receipts — default 6,760, `--marginal`
6,592, `--frame` 579,972 checks, all 0 failures. Checker walls: 118 s
(chrMT) + 366 s (chrI) = 484 s; the three prior modes ≈ 2 min each
(poll-bounded). Checker RSS peaks ≈ 3 GB (chrMT) / 5 GB (chrI),
observed well under the guard.

### What remains (the levers, unchanged in kind)

The candidate-universe/row-structure levers (row extents at window
seams; the alignment-induced partition boundary that fragments truth
rows; the scaffold/repeat pocket material and the cross-chromosome
windows), and the QUAL p-form's dynamic range. No selection swap, no
threshold, no scoreboard change; truth assessment-side only.

## Phase 0 of the runtime plan — the dominance measurement (2026-11-06)

The instrument's first wall-attribution slice (measure-before-lever;
timers and counters only, emitted to stderr — no receipt field
changes, no restructuring, no behavior change). `partition_realign_score`
gains phase-local timers at its natural phase boundaries — inputs
(panel+routes+census+partition maps), the FASTQ quality scan, the E
derivation, the derive-cache load, the pilot-record BINDING phase (with
sub-timers for its three cost drivers: the per-record canonical-scheme
step extraction, the per-candidate key-shape re-derivation, and the
per-(occurrence × shift) `verify_key_at` AGC fetch probes, plus measured
counts — records, occurrences, candidate keys, verify probes, fetch
calls, fetch bytes), the per-record read variants, the RowStore
partition build (context assembly: per-row fetch + stored-walk
extraction + canonical-scheme skeleton extraction, sequence-verified),
and per locus: the touching/locality classification, the scoring
matrix, the class enumeration, QUAL, the exactness sample, and the
receipts/IO writes, with a `locus_N_receipts` RSS probe per locus for
the memory-peak attribution. A global AGC-fetch counter attributes
fetch volume per phase.

### The timered exhaustive reruns (serial, exit 0, 64 GiB guard clean)

chrMT all 14 loci: wall 150 s, poller RSS peak 2,195,864 kB; chrI all
21 loci: wall 597 s, poller RSS peak 5,957,684 kB (the clean rerun;
the committed slice-E run measured 5,949,056 kB — the same peak).
Factorized-vs-direct 1,929,497 + 21,490,054 = 23,419,551/23,419,551
EXACT in both. **THE ANSWER-PRESERVATION GATE passes on both
components:** the timered receipts equal the committed slice-E
receipts on every semantic field (only the walls/rss_kb timing fields
differ; the walls field set unchanged), and all four sidecars
(exactness/ingredients/records/skeleton) are BYTE-IDENTICAL — the
instrumentation provably does not perturb the instrument.

### THE DOMINANCE TABLE (the measured phase walls; chrMT 147.6 s and
chrI 588.9 s in-process totals)

| phase | chrMT wall | chrMT share | chrI wall | chrI share | the measured cost driver |
|---|---:|---:|---:|---:|---|
| inputs (panel+routes+census+maps) | 5.4 s | 3.7% | 6.6 s | 1.1% | the 10M-vertex/24M-edge GBWT + 4,891/13,905 census records |
| quality scan | 2.7 s | 1.8% | 2.6 s | 0.4% | 2,439,103 reads streamed |
| E derivation | 0.002 s | ~0% | 0.002 s | ~0% | 235 haplotype sets |
| derive cache load | 1.4 s | 0.9% | 1.4 s | 0.2% | 2,439,103 reads / 666,327 keys |
| **binding (anchor/pin placement + verify)** | **134.5 s** | **91.1%** | **551.1 s** | **93.6%** | 4,891/13,905 records; 234,472/1,415,203 occurrences; 4,954/14,629 candidate keys; **703,605/4,247,781 verify probes; 708,496/4,261,686 AGC fetches (51.4/308.5 MB)** |
| read variants | 0.0 s | ~0% | 0.4 s | 0.1% | per-record variant census |
| rows (context assembly) | 1.5 s | 1.0% | 6.3 s | 1.1% | 1,252/4,472 rows; 2,504/8,944 fetches (16.5/78.2 MB) |
| loci: fold enum + units/locality + skeleton sidecar | 0.6 s | 0.4% | 2.0 s | 0.3% | Σ 1,099/3,545 folds; Σ 47,997/266,095 units |
| loci: scoring matrix | 0.2 s | 0.1% | 1.9 s | 0.3% | 1,929,497/21,490,054 placements (rayon-parallel) |
| loci: class LLs | 0.1 s | 0.1% | 3.0 s | 0.5% | Σ 59,225/378,306 classes × units mix_logsumexp (≈ 4.9 × 10^9 terms at chrI) |
| loci: QUAL + exactness | 0.3 s | 0.2% | 1.4 s | 0.2% | spectra + 3,360/5,178 sample pairs |
| loci: receipts/IO | 0.9 s | 0.6% | 12.1 s | 2.1% | the ingredients ll_matrix (Σ 47,997/266,095 units × folds; chrI L16 34.5 M entries) |

Within the binding phase the verify probes are **99.1% (chrMT) /
99.5% (chrI)** of the phase wall (133.3 s / 547.9 s); the step
extraction is 0.8/2.3 s and the key-shape derivations 0.0/0.1 s.

### THE INTERPRETATION — the lever mapping, measured

**The dominant phase is the binding phase's per-(occurrence × shift ×
candidate) sequence-fetch verification — 90% (chrMT) / 93% (chrI) of
the whole instrument wall.** The natural unit cost is per-CALL, not
per-byte: the average fetched verify range is ~72 bp (308.5 MB over
4,261,686 fetches at chrI) and the measured cost is ~129 µs/probe —
the random-access overhead of `sources.fetch` (per-call seek/decompress
into the route sources), not the compared bytes. The mapping to the
plan's levers:

- **(a) Candidate-independent recomputation (the hoist-to-once class) —
  THE lever this measurement names.** The verify loop fetches the SAME
  per-occurrence range once per (candidate key × shift): measured
  4,247,781 probes against 1,415,203 occurrences = a **3.01× pure
  re-fetch redundancy** (the candidate multiplicity is small here —
  14,629 candidate keys over 13,905 records — so the 3-shift loop is
  the redundancy). Hoisting the fetch to once per (occurrence,
  candidate) — or once per occurrence with the ±1 shift handled by
  fetching one base of margin — bounds the dominant phase at ~1/3 of
  its current wall (≈183 s at chrI) with NO model change. A per-range
  fetch cache (the binding phase's ranges and the rows phase's ranges
  overlap the already-fetched row sequences) would collapse the
  per-call overhead further. Second-order members of the same class,
  measured: the rows phase fetches each row's range twice (seq fetch +
  the canonical-skeleton fetch of [start−k, end+k)) — 8,944 fetches
  for 4,472 rows — at 6.3 s total; the per-(fold, unit)
  `decode_record_tokens` re-derivation inside scoring is inside the
  1.6 s scoring wall.
- **(b) Per-class re-derivation (the Phase-2 vote-vector lever):**
  the class enumeration re-walks all units per class — Σ classes ×
  units ≈ 4.9 × 10^9 mix_logsumexp terms at chrI — for 3.0 s of
  wall (0.5%; rayon-parallel; the CPU-seconds are larger but the wall
  is the measure). Real CPU, but not where the wall lives.
- **(c) Irreducible work:** the semantic per-occurrence verification
  (one fetch + per-anchor compare per occurrence — 1,415,203
  occurrences at chrI, each fetched range ~72 bp) and the scoring of
  21,490,054 placements. Even fully hoisted, the binding phase
  bounds at ~1/3 of its current wall; the remaining cost is the fetch
  mechanism itself.

**The memory-peak attribution (the 5.66 GiB):** NOT the binding phase
(RSS 3.21 GB at chrI) and not scoring (4.27 GB at L16) — the peak is
held in the per-locus **receipts/IO phase at the large-matrix loci**:
the ingredients sidecar materializes the full ll_matrix as a serde
JSON value tree before writing (chrI L16: 448 folds × 77,117 units =
34.5 M entries ≈ +1.15 GB live over the 4.27 GB baseline), and the
external poller catches 5.95 GB there (probe after the write: 4.80 GB
retained). A streaming serializer (to_writerf over the matrix without
the value-tree materialization) would remove the spike without
touching the receipts' content.

### The checker in the timered mode

`check-realign-scoring.py --exhaustive chrMT|chrI --timered-base
<base> --run-tag <tag>` re-runs every slice-E phase over the timered
receipts and adds **phase 8t, the answer-preservation gate**: per
locus, every semantic field must equal the committed slice-E receipt
(only walls/rss_kb excluded, the walls field set checked unchanged
against schema drift) and the four sidecars must be byte-identical.
**ALL PHASES PASS:** chrMT 1,170,642 checks, chrI 4,342,914 checks, 0
failures (the +730/+896 checks over slice E are phase 8t's own). The
default paths and all prior modes are untouched.

No thresholds, no tuning constants, no selection swap, no scoreboard
change; assessment-side only; serial runs; the 64 GiB guard clean
everywhere.

## Phase 1 of the runtime plan — THE FETCH-PATH REBUILD (2026-11-06)

The three levers the phase-0 dominance table named, all
answer-preserving (the receipts must reproduce slice E field for
field, the sidecars byte-identically — the identity-gate pattern from
slice D, enforced by the checker's phase 8t). Phase 2's vote-vector
lever stays unbuilt: the same measurement that named the levers
disproved it (0.1%/0.5% of wall).

### Lever 1 — the hoisted re-fetch

The binding verify loop re-fetched the same per-occurrence range once
per origin shift (the measured 3.01x redundancy: 4,247,781 probes over
1,415,203 occurrence-ranges at chrI). The three shift ranges differ
only in the low bound (`max(0, start + shift - 1)`, shifts −1/0/1) and
share the high bound `start + span + 1`, so one union fetch
`[max(0, start − 2), start + span + 1)` per (occurrence, candidate)
covers all three. `verify_key_at` (now top-level) verifies each shift
against the union buffer's tail from the shift's clamped low bound —
byte-for-byte the crop the per-shift fetch returned, including the
short/empty-crop edge cases (start near 0; range overhanging the path
end). Unit-proven against a copy of the pre-rebuild per-shift-fetch
code across interior/near-zero/path-end occurrence starts, both
orientations, both mirror states, both verdicts.

### Lever 2 — the per-lane in-memory AGC cache (the call-bound fix)

Phase 0 measured the fetch CALL-bound: ~129–190 µs per ~72 bp crop
(one random-access seek/decompress per `Sources::fetch` into the
AGC). Measured first, as ordered: the AGC is one 93.3 MB file (9,901
lanes, 3.34 GB uncompressed); the touched subset is 840 occurrence
lanes / 384.5 MB at chrI (141 / 11.5 MB at chrMT), 1,286 lanes /
514.7 MB with the partition-map member rows. **The derived choice:**
each TOUCHED lane loads once, whole, through the same validated
`Sources::fetch` path (one sequential decompression per lane — the
per-locality batch-prefetch option taken at its natural locality, the
contig; a window/block cache would still pay a decompression call per
window), then every crop is a memory slice of the cached buffer.
Byte-identity by construction: the lane buffer IS the 0..len crop
`Sources::fetch` returns (length- and alphabet-validated,
uppercased), so a crop of it is a crop of the answer. The cache holds
only touched lanes (≤ ~515 MB at chrI), never the 3.34 GB panel.

### Lever 3 — the streaming ingredients serializer

The ll_matrix now serializes directly from the typed matrix, not
through a serde_json value tree (phase 0 measured the chrI L16 tree
at +1.15 GB live over the 4.27 GB scoring baseline). Object keys are
written in the `json!` macro's map order (serde_json without
`preserve_order`: alphabetical), each field through the same
serializer — the sidecar stays byte-identical.

### THE WALL TABLE (before = phase 0, after = the rebuilt runs)

| phase | chrMT before | chrMT after | chrI before | chrI after |
|---|---:|---:|---:|---:|
| inputs | 5.4 s | 7.5 s | 6.6 s | 6.7 s |
| quality scan | 2.7 s | 3.2 s | 2.6 s | 2.7 s |
| derive cache | 1.4 s | 1.8 s | 1.4 s | 1.3 s |
| **BINDING** | **134.5 s (91.1%)** | **1.7 s (79x)** | **551.1 s (93.6%)** | **10.0 s (55x)** |
| read variants | 0.0 s | 0.0 s | 0.4 s | 0.3 s |
| rows/context assembly | 1.5 s | 0.7 s | 6.3 s | 2.3 s |
| the loci loop | 3.5 s | 3.2 s | 37.7 s | 11.4 s |
| **in-process total** | **147.6 s** | **18.1 s** | **588.9 s** | **35.4 s (16.6x)** |
| external wall | 150 s | **22 s** | 597 s | **37 s** |

Inside the rebuilt binding phase: the verify probes are unchanged as
a semantic count (703,605 / 4,247,781) but now run over 234,535 /
1,415,927 union range fetches (the measured 3.01x hoist), served from
141 / 840 lane loads (10.9 / 366.7 MB, 0.4 / 7.9 s of load time inside
the phase); the phase's fetch calls fell 708,496→239,426 (chrMT) and
4,261,686→1,429,832 (chrI). The factorized-vs-direct gate re-verified
every placement: 1,929,497 + 21,490,054 = 23,419,551/23,419,551 EXACT.

### THE MEMORY SPIKE (gate e)

The external poller peak at chrI fell 5,957,684 kB → **4,203,688 kB**
(the slice-E 5.66 GiB peak was held at the L16 receipts phase; the
lane cache itself adds 366.7 MB). The in-process receipts probe at L16
reads 4,422,448 kB after the write over the 3,846,188 kB scored
baseline; the residual ~576 MB is the records-sidecar named-classes
value tree — outside this lever's scope, named as the remaining
receipts-phase allocation. chrMT peak 2,077,092 kB (was 2,195,864).
The 64 GiB guard is clean everywhere.

### The gates

**(b) The identity gate:** the checker's `--timered-base` mode over
the rebuilt receipts (`--exhaustive chrMT --timered-base
realign-rebuilt-chrMT --run-tag realignrebuilt-chrMT-chrMT`, likewise
chrI): **ALL PHASES PASS — chrMT 1,170,642 checks, chrI 4,342,914
checks, 0 failures**, both including phase 8t (every semantic field
equal to the committed slice-E receipts, only walls/rss_kb differing,
walls field set unchanged; all four sidecars byte-identical). The
rebuilt runs are a proven no-op on the answers. **(c) The prior
modes** over their committed receipts: default 6,760, --marginal,
--frame, all 0 failures (the phase-0 no-regression sweep repeated).
**(d) Unit tests:** `partition_realign_score` 17/17 (16 prior + the
hoisted-verify equivalence test), `panel_route_mem_routed` 80/80,
`partition_anchor_projection` 8/8, `partition_graph_export` 1/1.

Receipts: `realign-rebuilt-{chrMT,chrI}.{jsonl,exactness.jsonl,skeleton.jsonl}`
+ sidecars + run markers (`run-realignrebuilt-*`, external RSS poller
files) + the checker logs (`check-rebuilt-*.log`) at the validation
dir; the committed slice-E and timered receipts untouched on disk.
Runner: `genome/instrumented/run-realign-rebuilt.sh`. Assessment-side
only (no src/ file touched); no thresholds, no tuning constants, no
selection swap, no scoreboard change, no 502, no PR push; serial; the
64 GiB guard clean everywhere.

What remains for the runtime plan: the parallelism pilot (the next
stage, not started here) — with binding collapsed to 1.7 s/10.0 s, the
loci loop and inputs are the next-largest walls, and the
per-locality work is now dominated by in-memory computation.

## Phase 3 of the runtime plan — the parallelism pilot, measured not assumed (2026-11-06)

The serial-only rule comes from a MEASURED race in the OLD machinery's
reference-rescore (7-wide and 2-wide both failed, solo passed; the
race file parked, one line). This instrument is receipt-side and
partition-local, but the discipline is measure-don't-assume: the race
detector here is BYTE-IDENTITY — every rung of the width ladder must
reproduce the committed phase-1 serial receipts on every semantic
field with all four sidecars byte-identical.

**The seams.** The two natural parallel seams in the phase-1 profile,
both in `examples/partition_realign_score.rs` (assessment-side; no
src/ file touched): **seam 1** the per-record binding loop (10.0 s of
35.4 s at chrI); **seam 2** the per-locus loop (11.4 s at chrI). The
shared-state proof is compiler-enforced, not asserted: each seam body
is a pure per-record / per-locus function — every outer capture is an
immutable shared reference (the panel, the census lines, the key
indexes, the derive cache, the reads, the partition maps, the
assembled row store, the bound keys/shifts/variants) or the
thread-safe fetch path (the `AtomicU64` fetch counters; the
`LaneCache`'s mutex-per-lane slots, whose loads hold the lane's mutex
through `Sources::fetch`, so loads are serialized per lane and
idempotent — no lane can load twice). The `&dyn Fn` row-store fetch
is widened to `&(dyn Fn + Sync)` and the seam map's bound enforces
`Sync` at compile time. Every seam result lands in its own pre-sized
slot read out in item order, and every sidecar line is buffered per
locus and written in LOCUS ORDER after the loop — the receipts are
thread-interleaving-independent BY CONSTRUCTION.

**The width semantics, measured twice.** The first attempt put the
seams on the rayon pool with `RAYON_NUM_THREADS=2`: the chrMT loci
phase went 1.8 s → 18.5 s. The control experiment isolated two
stacked effects: (1) the committed "serial" runs were never
single-threaded in the loci loop — the per-locus
scoring/classes/spectrum par_iters ride rayon's DEFAULT pool (the
whole box), and bounding the pool to the rung width silently
serialized those phases too (serial code at width 2: loci 8.2 s);
(2) nesting the seam par_iter in the same bounded pool adds
contention on top (18.5 s). NO RACE — even that pathological run was
byte-identical. The corrected mechanism (`seam_map`): the seams run
on DEDICATED scoped std::threads at width N
(`IMPG_REALIGN_PARALLEL_SEAMS=1 IMPG_REALIGN_SEAM_WIDTH=N`), and from
a non-rayon thread every inner par_iter installs rayon's global pool —
the committed inner behavior, unchanged. The committed default (no
env) is the serial path, byte-for-byte.

**The identity gate, every rung, both components: PASS.** All four
sidecars byte-identical to the committed phase-1 rebuilt receipts
(exactness/ingredients/records/skeleton), every semantic field equal
(only walls/rss_kb differ), binding counters identical (chrI 13,905
records / 1,415,203 occurrences / 14,629 candidates / 4,247,781
verify probes over 1,415,927 range fetches / 840 lane loads). No race
found at any width — the seam bodies touch no mutable shared state.

**The wall ladder (external wall, house runner, 64 GiB guard clean
everywhere; the poller RSS peak is the honest meter):**

| phase (chrI) | serial (re-baseline) | 2-wide | 4-wide | 8-wide |
|---|---:|---:|---:|---:|
| inputs | 7.0 s | 7.1 s | 6.1 s | 6.9 s |
| quality scan | 2.5 s | 2.6 s | 2.5 s | 2.6 s |
| derive cache | 1.4 s | 1.3 s | 1.3 s | 1.5 s |
| **binding (seam 1)** | **9.7 s** | **6.8 s** | **5.5 s** | **5.6 s** |
| rows/context | 2.5 s | 2.4 s | 2.4 s | 2.3 s |
| **the loci loop (seam 2)** | **13.1 s** | **14.7 s** | **16.1 s** | **24.2 s** |
| external wall | **41 s** | **34 s** | **30 s** | **30 s** |
| poller RSS peak | 5.71 GB | 5.52 GB | 6.26 GB | 6.66 GB |

(chrMT: serial 18 s / 2w 14 s / 4w 14 s / 8w 14 s; binding 1.1 → 0.8
s; loci 1.8 → 3.2 s at 8w; RSS peaks 2.31 / 2.33 / 1.73 / 1.79 GB —
the peak's locus alignment shifts with the stride schedule, poller
sampling noise included. The committed phase-1 reference walls under
their day's load: chrI 37 s / chrMT 22 s, chrI peak 4.20 GB — the
re-baseline serial peak reads higher because the emission buffering
holds the largest ingredients line (~706 MB at L16) in memory before
the ordered write; the same buffering is in every rung.)

**The honest reading, measured:** the binding seam pays (9.7 s →
5.5 s, plateauing at the serialized 840-lane load floor — the lane
loads hold their per-lane mutexes through the decompression); the
LOCI SEAM IS A WALL LOSS AT EVERY WIDTH (13.1 s → 14.7/16.1/24.2 s)
because the inner phases already saturate the box's default pool and
outer concurrency only adds live-loci memory pressure and scheduling
churn. The whole-instrument gain is 41 s → 30 s at chrI (~1.35x),
entirely from seam 1. The optimal measured rung is 4-wide; 8-wide
adds nothing on the wall and costs RSS.

**The scaling projection (EXTRAPOLATION, labeled as such).** A
chrIV-class chromosome (~1.53 Mb, ~6.7x chrI's window count) under
the measured per-phase rates projects: binding ~9.7 s x 6.7 ≈ 65 s
serial, ~37 s at 4-wide (the lane-load floor scales with the touched
lane set); the loci loop ~13.1 s x 6.7 ≈ 88 s serial (per-locus
matrix sizes assumed similar; the L16-class large loci dominate).
The whole genome (16 chromosomes, ~53x chrI) projects ~20 minutes
serial per full-genome pass and ~17 minutes at 4-wide — the seam
gain does not compound, because the loci seam is a measured loss and
the remaining wall (inputs ~6.7 s + quality 2.6 s + derive 1.3 s
constants + the serialized lane-load floor) is seam-invisible. The
levers for chrIV-class walls remain algorithmic (the named
records-sidecar value tree, the per-locus receipts/IO), not seam
width.

**The gates:** the checker `--timered-base` sweep over every rung's
receipts (all phases + 8t), plus the three prior modes over their
committed receipts; unit tests 17/17. Receipts `realign-par-{serial,2,4,8}-{chrMT,chrI}.*`
+ run markers (`run-realignpar-*`) + checker logs (`check-par-*.log`)
at the validation dir; the committed receipts untouched on disk.
Runner: `genome/instrumented/run-realign-parallel.sh`. Assessment-side
only; no thresholds, no tuning constants, no selection swap, no
scoreboard change, no 502, no PR push; the 64 GiB guard clean
everywhere.

## chrIV end-to-end, slice 2 — THE EXHAUSTIVE SCORING RUN (2026-11-06)

The measurement that replaces the whole-genome extrapolation with
arithmetic: all 167 chrIV loci under the marginal realignment model
with repaired skeletons, the anchor-projection receipts at 1.5 Mb
chromosome scale, and the domain/truth/wall tables against the chrI
profile.

**The input chain (the balanced-diploid validation's chrIV
component):** the multi-census run (`run-cosine-multicensus.sh chrIV`,
the committed pilots' recipe re-runnable per component) produced
92,883 census records (6.68x chrI) / 12,753,283 verified occurrences
(9.01x) / window span 0..166; the same run's Poisson-era read-matched
receipt measured the windowed frame at 0/167 truth-pair-expressible
(the slice-1 census's frame-gap prediction). The anchor-projection
receipts (`anchor-projection-chrIV.*`, 2.6 GB + context + rows
sidecars, exit 0, wall 704 s, poller peak 9.13 GB): 12,753,283/12,753,283
occurrences interval-count and node-list matched against the census,
shift 0 at every occurrence, 232 M hull checks; the sampled checker
`check-anchor-projection.py chrIV` ALL PHASES PASS (297,745,711
checks).

**THE REPEAT-DOMAIN REPAIR (the exact-DP anchor correspondence).**
The committed monotone-assignment enumeration is exponential in
anchors-per-repeat-node; the first chrIV serial attempt measured
locus 39 alone at 760.7 s in scoring (partition 175's repeat rows
carry a syncmer node at 34 positions; 5,788 units x 179 folds), every
other completed locus <= 0.6 s. The multiplicity map over the 158
axis partitions names the hazard set (L71 a scaffold row with a node
at 1,356 positions — benign there only because few units pin on it;
L147 338; L78 103; L81 101; L166 96). The repair replaces the
enumeration with an exact forward/backward feasibility DP: every
complete monotone chain contains every candidate-bearing anchor, so a
position is pinned iff it is the anchor's unique DISTINCT feasible
candidate; >= 2 feasible = ambiguous None; none = the zero-leaf
all-None outcome. Equivalence gates: unit test 18/18 (3,000
randomized synthetic walks vs a verbatim copy of the committed
enumeration — it caught the duplicate-position counting divergence on
the way) and the chrMT/chrI serial byte-identity gates (all four
sidecars byte-identical to the committed phase-1 rebuilt receipts, 0
semantic diffs). With the DP, locus 39 scores in 0.0 s with
IDENTICAL numbers (593,870/593,870 factorized placements).

**THE WINDOW GENERALIZATION.** chrIV is the first component where
windows are not 1:1 with axis partitions (167 windows over 158
partitions; six partitions host 2-4 windows). Windows derive from
the AXIS ROWS; the truth folds are window-aware (truth-first = the
fold carrying the window's own axis row; truth-second = the committed
path-name rule at single-window partitions, None at multi-window ones
where an SK1 row cannot be attributed to a window by placement —
chrIV partition 110's two SK1 rows are other windows' orthologs).
Byte-identity-gated on chrMT/chrI (all sidecars byte-identical).

**THE EXHAUSTIVE RUN OF RECORD (4-wide, the phase-3 clean rung; exit
0; 64 GiB guard clean):** external wall 148 s (serial identity base
236 s), poller RSS peak 28.58 GB (serial 26.79 GB); in-process phases
(serial): inputs 12.0 s, quality 2.2 s, derive 1.1 s, binding 39.9 s,
variants 0.9 s, rows/context 16.0 s, loci 138.0 s; at 4-wide the
per-locus summed walls inflate (loci 260.4 s summed — the phase-3
contention finding at chromosome scale) while the external wall
falls. THE IDENTITY GATE: the 4-wide receipts equal the serial
receipts on every semantic field, all FOUR sidecars byte-identical —
no race at chrIV scale. Receipts `realign-exhaustive-chrIV.*` +
`realign-par-serial-chrIV.*` + run markers at the validation dir.

**THE DOMAIN-SCALING TABLE (chrIV vs chrI, measured):** windows 167
vs 21 (7.95x); class pairs post-fold 4,279,081 vs 378,306 (11.31x);
units 2,502,395 vs 266,095 (9.40x); records 128,215 vs 16,020
(8.00x); factorized placements 202,871,144 vs 21,490,054 (9.44x);
per-locus means 25,623 class pairs vs 18,015 (1.42x) and 14,984 units
vs 12,671 (1.18x) — SUPERLINEAR in the class domain, dominated by the
repeat-locality partitions: the largest loci are L92/L106 (partition
24: 1,072 folds, 575,128 class pairs, 72,105 units, 18.4 M
placements), L51 (partition 187: 956 folds, 457,446 pairs), L166
(partition 289: 700 folds, 245,350 pairs — the chromosome-end
subtelomeric locality), L66/L93/L107 (partition 16: the chrI-L16
twin, 448 folds, 100,576 pairs), L121, L120/L134.

**THE TRUTH-RANK TABLE (all 167 loci; classes from the slice-1
expressibility census):** 143/167 truth-pair-expressible under the
committed path-domain convention; truth rank-1 at 70 loci. By class:
IN-AXIS-PARTITION 96 expressible, 48 rank-1; TILED-ELSEWHERE 46
expressible (the one-window coordinate-offset class — the window's
partition carries the NEIGHBOR window's aligned SK1 row; the paired
class is named per locus in the checker's phase 8c, never silently
passed as the ortholog), 21 of them rank-1; TILED-ELSEWHERE 23
inexpressible (13 multi-window partitions, 8 no-SK1-row, L120/L134
the ambiguous multi-window pair); ABSENT (the contig-end length
polymorphism) — L163 EXPRESSIBLE and truth rank 1 (its axis
partition 286 carries SK1's [1454080,1464091) row, the terminal
material the slice-1 census's window query found), L166 inexpressible
(its 289 partition carries no SK1#0#chrIV row at all; the winner
harvests the cross-chromosome subtelomeric repeat family — AAA/SGDref
chrXV + SK1 chrV folds — and is the only bounded-QUAL locus, 92.89).
The identical pair AAA#0#chrIV/SGDref#0#chrIV folds to ONE candidate
at every locus where either holds a row (checker phase 7: co-membership,
identical intervals, identical sequence+walk from the GFAs).

**THE WALL/MEMORY VERDICT vs THE EXTRAPOLATION:** the phase-3
projection guessed chrIV-class binding ~65 s serial -> ~37 s at
4-wide and a whole-instrument wall in the low minutes with an RSS
guess of 10-20 GB. Measured: binding 39.9 s serial -> 22.5 s at
4-wide (better than projected — the lane-cache floor scales with the
touched lane set, 3,399 lane loads / 1.7 GB at the anchor layer and
12,928,943 scoring fetches served from cache), the serial
whole-instrument wall 236 s and 148 s at 4-wide — the honest
arithmetic for the whole genome: chrIV is 7.95x chrI in windows but
the wall is 3.6-4.9x chrI's (30-41 s), sublinear because the
component constants (inputs, quality scan, derive cache) amortize;
the whole genome (16 components, ~53x chrI) projects to ~15-20
minutes at 4-wide per full-genome pass with the DP in place. The RSS
peak 28.6 GB exceeds the 10-20 GB guess — the repeat-locality loci's
matrices (L92's 72,105 units x 1,072 folds) are the cost; still well
under the 64 GiB guard.

**The gates:** the checker in all modes — `--exhaustive chrIV
--timered-base realign-exhaustive-chrIV --run-tag
realignexhaustive-chrIV-chrIV` over three --loci-slice invocations
(every check per sliced locus; the slicing is a wall measure, not a
sample), the anchor checker, the three prior modes unchanged on their
committed receipts, unit tests 18/18 + 8/8; tables receipt-side via
`realign-chrIV-tables.py`. The checker repairs the slice completion
itself exposed, all checker-side, the receipts unaffected at every
step (each failing check was one the receipts never claimed to
satisfy): phase 8c assumed the instrument orders truth_folds
(axis-row fold first — the instrument emits the pair sorted by fold
index, so the SK1 fold can carry the lower index; 10 loci affected,
fixed to the order-independent membership proof), and phase 7's
twin and generic identical-pair branches compared members' FULL
P-line spellings (and, in one repair round, the contained walk
against the fold's CLAIMED walk — the claimed walk is the REPAIRED
skeleton in this mode) — the fold key is (row sequence, contained
STORED walk), and full spellings legitimately extend past the row
extent differently (partition 110's 74bp repeat fold: 102 members,
AMP_1a#0#chrIII_chrX/ANL/AVN members spell 135bp P lines;
partition 111's 54bp fold at L91/L105: 137 members incl. fused-path
rows) — fixed to the fold key itself, cross-member equality of
the contained stored walk. The three prior modes re-run under the
final checker: 6,760 / 6,592 / 579,972 checks, 0 failures each. Assessment-side only; no thresholds; the
truth classes come from the committed slice-1 census receipt; no
selection swap, no scoreboard change, no 502, no PR push; the 64 GiB
guard clean everywhere.

## The whole-genome fleet — all 17 components under the realignment instrument (2026-11-06, in flight)

The owner's go: the remaining 14 balanced-diploid components
(chrII, chrIII, chrV, chrVI, chrVII, chrVIII, chrIX, chrX, chrXI,
chrXII, chrXIII, chrXIV, chrXV, chrXVI) through the chrIV recipe
end-to-end — survey, partition-graph build, expressibility census,
anchor-projection receipts, the serial identity base + the 4-wide
exhaustive run of record (the marginal model with repaired
skeletons, the 64GiB guard, external RSS polling), the corrected
checker in its component-sliced modes, the per-component tables.
chrMT, chrI and chrIV stand closed. The machinery is the committed
chrIV-stage instrument, component-parameterized with every
generalization gated byte-identical on chrIV (commit 87c0ded); the
census's chrIV-only windowed-frame assert became a reported count
(the fleet's balanced receipts carry windowed-frame truth-pair-
expressible loci at most chromosomes — the phase-8 old-wins
machinery is live genome-wide, unlike chrIV's vacuous 0/167 guard).

Per-component sections land below as each chain closes; the
whole-genome table is emitted by
`genome/instrumented/genome-truth-rank-tables.py` (re-runnable
mid-fleet; pending components reported as pending). Assessment-side
only; no thresholds; no scoreboard change; serial within component
where the recipe requires it.

### chrVI — the first fleet component closed end-to-end (2026-11-06)

The chrIV recipe, all gates: survey 47-partition build set / 29 loci;
twins AAA#0/SGDref#chrVI (100% exact) + ALH_1b/ALH_1c (interval-close,
the CLL class); census IN-AXIS 16 / NEIGHBOR 8 / FOREIGN 5 / CONTIG-END
0; anchor 191s (peak 3.9GB); serial base 67s (8.4GB); 4-wide run of
record 62s (9.6GB), byte-identity gate passed (phase 8t vs the serial
receipts); anchor checker ALL PHASES PASS (59,482,753 checks); realign
checker ALL PHASES PASS; graph checker ALL PHASES PASS. THE NUMBERS:
25 of 29 loci truth-pair-expressible under the path-domain convention,
TRUTH RANK-1 AT 15 (in-axis 11/16, tiled-elsewhere-neighbor 4/9); the
windowed frame's 2 old wins BOTH HOLD (no prediction violations), 23
loci newly expressible; the largest locus L28 (partition 376: 30,496
units, 8,778 classes) ranks the truth 525th — the repeat-domain
residual class. Census wall 937s (peak 23.1GB). The per-locus table:
`realign-chrVI-tables.txt` beside the receipts.

### chrIII — closed end-to-end (2026-11-06)

38 loci / 48-partition build; twins AAA/SGDref only. Census IN-AXIS 19
/ NEIGHBOR 11 / FOREIGN 8 / CONTIG-END 0. Walls: census 825s (24.8GB),
anchor 354s (4.5GB), serial 105s (11.7GB), 4-wide 65s (13.6GB), the
8t byte-identity gate passed. All three checkers ALL PHASES PASS
(anchor 111,671,505 checks; realign 8,932,744/0; graph). TRUTH RANK-1
AT 13 OF 28 EXPRESSIBLE (in-axis 10/19, tiled-elsewhere-neighbor 3/9);
the windowed frame's old wins: 7 windowed-expressible before (see the
phase-8 aggregate in the checker log); table at
`realign-chrIII-tables.txt`.

### chrVIII — closed end-to-end (2026-11-06)

60 loci / 74-partition build; twins AAA/SGDref only. Census IN-AXIS 44
/ NEIGHBOR 6 / FOREIGN 5 / CONTIG-END 5 (SK1's chrVIII ends short —
the chrIV L163/L166 contig-end class). Walls: census 566s (26.6GB),
anchor 197s (4.3GB), serial 96s (8.0GB), 4-wide 55s (8.4GB), the 8t
gate passed. All three checkers ALL PHASES PASS (anchor 81,539,693;
realign 11,431,917/0; graph). TRUTH RANK-1 AT 25 OF 50 EXPRESSIBLE
(in-axis 23/44, tiled-elsewhere-neighbor 2/6); table at
`realign-chrVIII-tables.txt`.
