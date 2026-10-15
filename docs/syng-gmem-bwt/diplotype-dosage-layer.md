# The dosage surface — the copy-number layer of the two-molecule emission

The product chain: `impg genome-infer call-diplotypes` (the per-locus
diplotype calls with quality) -> `chain-molecules` (the two
whole-chromosome molecules per component) -> **`emit-dosage`** (the
per-panel-segment copy count the two molecules imply). The dosage is
DERIVED from the chain's molecules.jsonl and from nothing else on the
product side: the same material, counted — never re-selected, never
re-inferred. This document is the layer's design of record.

## The dosage emission (truth-free, product-side)

Per panel path, the two molecules' member rows decompose into
**elementary segments**: the breakpoints are every row's start and end,
and the segments are the *covered* consecutive breakpoint pairs — every
row CONTAINS every elementary segment it touches (no partial overlaps
exist by construction), and the gaps between disjoint rows are not
material (0 copies by construction; the per-path summary states
`covered_bp` against the path length so the complement is visible).

A segment's **copy count** is the number of traversals of it by
molecule 1 plus the number by molecule 2 — one traversal per
(molecule, locus, member row) containing the segment:

- both molecules traverse -> 2; one -> 1; absent -> 0 (the complement);
- **repeat-visit multiplicity is carried**: a molecule traversing a
  segment twice (two loci, or two rows of one fold) counts its true
  multiplicity, so counts above 2 are emitted where the walk revisits.

The **copy-bp identity** is exact and gated: the segments tile the
rows, so `SUM(dosage * length)` over the emission equals the summed
member-row lengths (`row_bp`); the gate proves it.

Per segment the record states the **territory/axis extents**: every
visit carries its locus's axis window (the partition-graphs' own
derivation, the same window spine the chain consumes), and the segment
carries the merged axis extent of its visits. Per locus the record
states the restricted dosage view (segment count, covered bp, the
copy-count histogram, and the **flat-class-of-2 flag** — every covered
segment carried by both folds, the old single-class-of-2 signature
measured at the coverage level; a homozygous emitted pair is flat by
construction).

## The depth-consistency QC (truth-free, a derived field — the ratio stated, never binned, no thresholds)

The observed per-segment census mass vs the expected under the emitted
dosage. The likelihood's own depth model is the balanced diploid
prior: the class likelihood's **equal 1/2 mixture weights** over the
two homologs are the statement that each molecule's covered material
carries the same expected read mass per covered bp, per covered copy.
The per-covered-copy expectation is therefore uniform per bp:

    mu = attributed mass / emitted copy-bp     (derived per component
                                                from the census and
                                                the emission itself,
                                                never a tuned constant)
    expected_mass(s) = mu * dosage(s) * length(s)
    mass_ratio(s)    = observed_mass(s) / expected_mass(s)   (stated,
                                                null when expected = 0)

The observed mass is the census's own record-once material
attribution: a record's multiplicity-weighted covered bp — the union
over its occurrences of the contained-node material windows on the
occurrence's path (the same `(bp, bp+k)` windows the chain's edge layer
walks, k the syncmer length) — attributed to a segment by overlap bp.
The aggregate states the totals: the census mass, the material mass,
the attributed mass, and the **unattributed material mass** (the
census's multi-matching spread on panel material outside the emission,
counted, never hidden). No segment is binned, judged, or filtered —
the ratio is the statement.

## The test mode (assessment-side ONLY, the --truth-qv-file convention)

The product run's own `calls.jsonl.truth-qv.jsonl` names the truth
pair per locus (fold indices into the product run's own instrument
receipt, whose `fold_identities` carry the truth folds' member rows —
the receipt is a truth-free product artifact; the truth REFERENCE is
the test mode's own input, `--truth-qv-file` + `--folds` required
together). The truth pair's per-segment dosage is derived by the same
rule as the emission, and the agreement table per component states:

- per locus, the emitted pair's dosage vs the truth pair's over their
  union universe, with the honest classes the data names:
  - `correct_dosage` — the copy vectors agree;
  - `dosage_inherited_from_non_rank1_call` — the call itself is not
    rank-1; the dosage error is the call's, not the emission's
    (mechanism: `territory_split` vs `multiplicity_mismatch`);
  - `dosage_specific_emission_defect` — the locus IS rank-1 and the
    dosage still differs (necessarily a winner-pair mismatch at a
    tie-degenerate locus; same mechanisms);
  - `non_rank1_call_dosage_agrees` — a wrong call whose traversal
    structure coincides with the truth's (stated, never counted
    correct);
  - `bracketed_truth_pair_not_expressible` — no truth dosage exists;
- per segment, the component-wide agreement over the union universe
  with the same decomposition: `agree` / `inherited_from_non_rank1_calls`
  / `dosage_specific_at_rank1_locus` / `bracketed_locus_involved`
  (a segment a bracketed locus touches — no truth total can judge it);
- the per-locus product accuracy columns (`rank1`, `perfect`,
  `truth_pair_expressible`) carried UNCHANGED from the product
  artifact in separate columns, never conflated — the owner's
  evaluation ruling since the haploid era.

## The gates

`dosage-gate.py <component> <molecules.jsonl> <dosage-dir>
<partition-graphs-dir> <truth-qv.jsonl>`:

1. **The thin-layer proof** — the consumed molecules fingerprint
   (FNV) equals molecules.jsonl now (the dosage never wrote it), and
   every segment, visit, copy count, per-locus dosage view, and
   per-path summary is RE-DERIVED in independence from molecules.jsonl
   (+ the axis windows) and compared field-identically with
   dosage.jsonl — the dosage is bit-reproducible from the molecules.
   (The census-derived observed mass is the layer's own new
   measurement — the chain gate's own discipline: it does not
   re-derive the census walk; its arithmetic is proven instead.)
2. **The QC arithmetic** — `mu` is the derived expectation;
   `expected_s == mu * dosage_s * length_s`; `ratio_s ==
   observed_s / expected_s`; the sums equal the emitted totals; the
   copy-bp identity.
3. **Zero truth keys** in the default emission (no key naming
   truth/rank1/perfect/assignment anywhere in dosage.jsonl).
4. **The agreement table** — the carried product columns match the
   product run's artifact; the class counts, the correct-dosage
   fraction, the rank-1 agreement, and the per-segment classes are
   re-tallied from the per-locus lines and summed to the union
   universe.

chrMT: ALL CHECKS PASS 422/0 — 9/9 expressible agree, 9/9 rank-1
agree, zero differing segments (23 of 62 union segments involve
bracketed loci, honestly carried; the one flat-class-of-2 locus is
one of chrMT's 5 bracketed loci). chrI: ALL CHECKS PASS 689/0 — 14/18
expressible agree, **14/14 rank-1 agree**, the 4 disagreeing loci
exactly the committed non-rank-1 residuals (L2/L13/L16/L18), 0
flat-class-of-2 loci of 21.

## The old defect's measure

The old era's `genotype_calls.diplotype.dosage` emitted ONE class of
2 at 1,041 of 1,307 loci (only 266 two-1-copy classes) versus truth
1.0/1.0 everywhere — the named surface defect, separate from
selection. The new surface measures the disease at the coverage level
(the per-locus `flat_class_of_2` flag, truth-free) and judges it at
rank-1 loci (the test mode): at every rank-1 locus the emitted pair IS
the truth pair, so the emitted dosage equals the truth dosage
structurally — the fleet aggregate states the numbers.

## The fleet aggregate (17/17, `dosage-aggregate.py`)

Every component gated ALL CHECKS PASS (the thin-layer proof, the QC
arithmetic, the agreement table). The fleet totals:

- **Correct dosage 880 / 1,132 expressible = 77.74%** — and
  **rank-1 dosage agreement 880 / 880 = 100%**: the dosage surface
  is a faithful function of the call (correct exactly at the
  rank-1 loci, inherited disagreement exactly at the 252 non-rank-1
  expressible loci — the residue's own material, no dosage error
  that selection did not already make).
- **THE HONEST DECOMPOSITION**: per locus — 880 correct, 252
  `dosage_inherited_from_non_rank1_call`, 175 bracketed (truth pair
  not expressible), **0 dosage-specific emission defects**, **0
  non-rank-1 calls whose dosage coincides with the truth's**; every
  differing locus's mechanism is `territory_split` (252),
  `multiplicity_mismatch` 0. Per segment over the union universes
  (9,527 segments): 6,954 agree, 1,922 inherited from non-rank-1
  calls, 651 bracketed-locus-involved, **0 dosage-specific at
  rank-1**.
- **The emission**: 8,136 segments, copy-bp 77,959,457 (the row-bp
  identity gated at every component); the depth-consistency QC
  states the ratio per segment with NO binning — the fleet's stated
  distribution: median 0.686, p10 0.549, p90 1.511 (the census's
  own multi-matching spread: much read mass lands on panel material
  outside the emission — the unattributed material mass is stated
  per component, never hidden).
- **THE OLD DEFECT'S VERDICT**: the 1,041/1,307 single-class-of-2
  era is GONE — the emission states 45 flat-class-of-2 loci of
  1,307 (3.44%), and the 45 decompose honestly: 3 at rank-1 loci
  where the TRUTH pair itself is coverage-flat (chrXII loci whose
  two homologs the panel represents by one coverage; the emission
  matches truth exactly there), 42 at bracketed loci (truth
  unknowable). At every expressible locus where the truth is mixed
  (1,129 of 1,132), the emission is mixed — the real heterozygosity
  at rank-1 loci is carried in the per-segment copy counts, 880/880.
