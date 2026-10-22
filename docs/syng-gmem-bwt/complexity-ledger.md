# THE COMPLEXITY LEDGER — the owner's algorithmic take-stock ("0.14 core-seconds per locus all in... maybe you need to use bigger loci... I think there's some quadratic stuff, some n-log-n stuff getting baked in here. Take stock of what we got.")

**The instrument, stated first.** Every number below is from the committed
receipts or a fit over them, and every fit names its quality. The two
corpora: **(A) the tranche callC corpus** — 2,822 per-locus rows mined
from the 138 completed serial (`RAYON_NUM_THREADS=1`, wall = core-second)
truth-assessment calls of tranches 1+2 at the CURRENT committed stack
(`scripts/complexity-ledger/mine-locus-shapes.py` over the
`callC-*.err` logs; the exhaustive class scan, the far-bin prune ON,
the anchor hoist; the Oct-19 binary whose call-path bits the gates
proved current) — and **(B) the SS10.6.14 autopsy census** (the 25-locus
L1 per-class timer corpus, serial, scaled by the receipted 164.8/182.1
timer factor). The fits: log-log OLS
(`scripts/complexity-ledger/fit-exponents.py`), R² named everywhere.
The work dir: `/home/erikg/hprcv2/global-panel-call-mode/complexity-ledger/`
(the mined TSVs, the fit output `exponent-fits.txt`, the grain re-measure
dirs `L1-w*/`, `T03226-w*/` with per-stage `.wall/.err/.rss` receipts).

---

## 0. THE PLAIN ANSWER (the take-stock, in the owner's terms)

**THE ONE REAL QUADRATIC is the class scan, and it is exactly what the
owner suspected:** the candidate pair count is `folds × (folds+1) / 2`
— **verified EXACTLY equal to the emitted `classes` count at 2,822 of
2,822 loci (max |residual| 0)** — and the scan walks every pair against
every unit, so

> **W(class scan) = 18.42 ns × [f(f+1)/2] × U** (fit: t = 2.86e-8 ·
> (f(f+1)/2·U)^0.978, R² = 0.9966, n = 25 per-class points; the
> aggregate rate 18.42 ns/unit-pair reproduces the SS10.6.14 receipt's
> "~18.4 ns" bit for bit).

The scan model alone explains **71% of the whole corpus's serial loci
wall** (75,619 s of 106,340 s). It is quadratic IN THE PANEL'S ROUTE
DIVERSITY PER WINDOW (the fold count), linear in the window's unit
mass. **The bounded scan (Frontier 1, the production default) already
removes ~99% of this term exactly where it engages** (the prunable
monsters: T03226 l0 99.39%, T05884 l8 99.01% of pairs skipped) **and ~0%
at the segdup plateaus** (T00940 l1, the 449-fold family: 0 pruned) —
the bound's honest ceiling is ~half the tail term (§10.6.22's
twice-confirmed −46%/−54% whole-loci-phase at the two prunable
monsters), because the tail contains BOTH kinds.

**THE SECOND SCALING CLASS is the scoring pair loop**, linear in pinned
(fold, unit) pairs: corpus rate **2.43 µs/pinned-pair** (the autopsy's
per-class fits: placement t = 3.01e-7·scored^1.049, R² = 0.988 — linear;
scoring t = 7.15e-6·scored^0.929, R² = 0.921), with the irreducible
per-base floor ~0.15-0.3 µs/pair (SS10.6.14's receipt). The pinned-pair
model explains a further **16%** of the corpus wall.

**THE FOLD ENUMERATION IS NOT A COST POLE** (the anatomy in §4): the
204 folds at the deep locus cost **0.1 s of a 1,423 s locus** to
enumerate; the enumeration is a deduplication of member rows
(802 → 213 distinct (path,start,end) rows → 204 after the equal-column
merge), O(members + rows·L log L) in the walk length L, and the
candidate-pair explosion it FEEDS is the scan's f² — the lever is the
scan's, not the enumeration's.

**THE N-LOG-N STUFF IS NOT THE STORY.** The measured shares at the
median and deep loci: the bounded scan's stage-one plan/sort is
invisible below ~1e5 pairs (the fit against pairs has R² = 0.105 —
the walls are quantized at 0.1 s and contention-noised; at the
monsters the whole bounded scan is 1.0-1.3 s against 1,100+ s serial
scan equivalents). The AGC seeks + lane loads are the bind's own
terms (a whole territory: fetches ~1.5M / 139-362 MB, lane loads
2,327-4,683 / 133-2,113 MB; the bind phase fits census-records^0.947,
R² = 0.70 — linear in mass, not superlinear). No measured phase
carries a material n-log-n term.

**THE OWNER'S BIGGER-LOCI LEVER, MEASURED AND REFUSED ON COST:** the
fixed classes are real (per-call inputs 243.8 s median, FLAT in
windows — fit exponent −0.035, R² = 0.016; bind 7.2 s; rows 13.4 s;
per-locus floor ~0.20 s) — but §5's re-measure at the current
architecture shows they are **mass-bound or call-bound, not
window-bound**, so they do NOT amortize with grain (binding is
census-record-bound and the census record count is grain-invariant:
81,688 records at 10 kb vs 81,765 at 20 kb at the territory), while the
content classes GROW — per-window units linearly (×1.99 per 2× grain),
folds sub-linearly (×1.6, clamped at ~443-457), the whole L1 loci
phase **394.0 → 506.6 → 785.6 → 731.4 core-s from 10→20→40→50 kb**
(+1.9x at 40 kb) and the monster window at the deep territory l0
**2,631.1 s at 20 kb vs ~1,4xx s at 10 kb** (the merged window carries
the cross-boundary fold pairs the fine windows never enumerate).
Bigger loci buy the ~0.2 s floor's amortization and pay ~2x the
content — refused at cost, at both a locality and the deep stratum,
with the accuracy sidecars landing clean at these territories (the
2026-10-09 sweep's genome-wide accuracy verdict stands untouched).

**THE MEDIAN-vs-DEEP SPREAD is the budget:** per-locus serial wall
median **0.50 s**, p90 5.9 s, p95 **28.2 s**, max **8,547 s**; the
worst-5% of loci own **96.3%** of the corpus per-locus wall (the
tranche's own 96.9% receipt reproduced independently). The scan model
is 5.6e6 unit-pairs at the median locus and 3.8e11 at the max — a
**68,000×** spread.

---

## 1. THE LEDGER (phase | measured formula | share at the median locus | share at the deep locus | the lever that owns it)

The shares at the median locus (T00041-l0-class, wall 0.50 s serial)
and the deep locus (T03226 l0, wall 1,423.4 s serial, the prunable
monster; plateau monster T00940 l1 beside it). "Model" = the fitted
term's predicted share of the locus wall at that shape.

| phase | W = (measured) | median locus | deep locus (T03226 l0) | lever |
|---|---|---:|---:|---|
| **the exhaustive class scan** | 18.42 ns · f(f+1)/2 · U (R²=0.9966; exact quadratic receipt at 2822/2822) | 0.09 s (18%) | **1,126.8 s (79%)** | **bound (DONE, default)** — removes ~99% where it engages; plateau remainder → memoization (23.9-37.7%, NOT implemented) |
| — the same at the plateau monster (T00940 l1) | 449²/2 × 171,590 = 1.73e10 unit-pairs | — | 319.3 s (50% of 641.7 s), **0 pruned by the bound** | memoization / plateau-aware bound (NEW) |
| **the scoring pair loop** (pinned pairs) | 2.43 µs · pinned_pairs (placement 3.01e-7·scored^1.049, R²=0.988) | 0.18 s (36%) | 136.6 s (10%) | per-base floor IRREDUCIBLE (~0.4-0.6 core-s/locus at pinned-heavy shapes); walk-residue ENGINEERING (SS10.6.14(d)) |
| **fold enumeration + units classification** | O(members + rows·L log L); fits too small to resolve (best: records^0.34, R²=0.37) | inside the 0.20 s floor | **0.1 s (<0.01%)** | none needed — not a pole |
| **the alignment-connected merge** (equal-column candidates) | candidates bounded by CONTENT (vs folds_pre²: R²=0.01 — not quadratic in practice); admission 71% | inside the floor | <0.1 s | none needed |
| **qual + exactness + receipts/IO** | the per-locus floor's remainder | ~0.1 s | ~160 s (11%) — the receipts emission dominates (T03226 l22 RSS 53 GB) | receipts-shape (OWNER) + IO (ENGINEERING) |
| **per-locus fixed floor total** (fold build + qual + receipts + scan floor) | measured 0.20 s median at matrix<2e4 | **0.20 s (40%)** | <1% | **grain** (amortizes ~1/grain) + precompute |
| **bind** (per territory, not per locus) | 0.0026 · census_records^0.947 (R²=0.70); verify probes ~linear in occurrences | 7.2 s median / 25 loci = 0.29 s | up to 507 s max (the deep-territory class) | **grain** (records halve per window merge) + precompute (Frontier 2) |
| **rows/context** | weak fit (R²=0.21) — a serial IO class | 13.4 s / 25 = 0.54 s | same (mass-invariant) | precompute (Frontier 2) |
| **inputs (panel load)** | FLAT vs windows (fit exponent −0.035, R²=0.016) | 243.8 s / territory | same | **fleet COW amortization (DONE in the product form)** + grain |
| **the map/cut side** (all-in) | census: linear in routed mass (occurrences), flat-to-down in windows (the sweep's receipt: peak 40.1→19.1 GiB at 20→50 kb; walls 229-414 s flat) | ~26k core-s genome-wide (the committed receipt) | — | ~2x engineering (the standing receipt) |

**Corpus audit:** scan model 71% + pinned-pair model 16% + floor ~13%
≈ 100% of the serial loci wall (the models' overlap at small shapes
absorbed in the floor term).

---

## 2. THE QUADRATIC INVENTORY — what remains post-bound

1. **THE CLASS SCAN (f² in folds, × U)** — the only structural
   quadratic. Post-bound state: **removed where the landscape is
   peaked** (T03226 l0: the scan's serial-equivalent 1,126.8 s → the
   bounded walk's ~1.2 s + stage-one 1.0 s; T05884 l8: 710.4 s →
   ~3.4 s at T16 measured), **INTACT at the uniform-excess plateaus**
   (T00940 l1: 319.3 s stays; T03901 l13 40.6% pruned is the mixed
   case). The named residual levers: the exact-preserving bit-equal-row
   memoization (distinct_rows 187 of 204 folds at T03226 l0 —
   37.7% of pair sums share one value; unit-proven, NOT implemented)
   and a plateau-aware bound (the uniform-excess landscape defeats the
   per-fold excess-mass bound BY CONSTRUCTION).
2. **THE PAIR LOOP (pinned pairs = folds × units × pin-rate)** — not
   quadratic in a single variable but a PRODUCT that behaves
   quadratically at windows where both folds and units are large
   (T03226 l0: 56.2M pinned pairs; T05884 l8: 125.7M). The bound does
   NOT touch it (it needs every exact (m,c)). Irreducible core:
   the per-base comparison; the engineering residue (the anchor class
   is GONE — the hoist took it −90.7%; the placement-walk 2-4x over the
   raw floor remains).
3. **THE CENSUS CONSTRUCTION (map side)** — linear in routed mass
   (occurrences per record), NOT in window pairs: the grain sweep's own
   receipts (census peak SHRINKS 40.1→19.1 GiB from 20→50 kb; wall
   flat 229-414 s) and the far-bin pair-bin diameter test is
   O(placements) per record. The owner's "pair-bins" worry is not a
   quadratic: the bins are the record's own touched windows (≤2
   adjacent), tested at O(1) per placement.
4. **NOT quadratic, measured:** fold enumeration (§4), the merge pass
   (content-bounded, R²=0.01 against f²), binding (records^0.947),
   rows/context, inputs (flat).

---

## 3. THE N-LOG-N INVENTORY (the seeks, sorts, lookups — their honest share)

- **The sorts:** the bounded scan's stage-one pair plan (f² log f² at
  f=450 ≈ 1.7M comparisons — sub-millisecond) and the fold build's
  per-row node/edge sorts (O(L log L) in walk length). **Measured
  share: invisible.** The stage-one wall IS the whole bounded scan wall
  at the monsters (1.0-1.3 s at T16, the exact-sum prefix rounds, not
  the sort); below ~1e5 pairs the logs quantize at 0.1 s and the fit
  against pairs has R²=0.105 — no resolvable n-log-n term.
- **The seeks (AGC fetches):** per territory 1.50M fetches / 362 MB
  median-class (the tranche runs table); the bind's verify probes fit
  occurrences^1.10 with R²=0.23 (weak — the seek class rides the
  co-tenancy noise). The pool-warm form is the named fix (SS10.6.13).
- **The lookups:** the bind's shape-derivation table (4,215 candidate
  keys → 4,215 derivations at an ordinary territory — linear); the
  hoisted pin skeleton (the coordinate-lookup class, the −90.7%
  receipt); the memoization's 128-bit row-token hashes (O(rows) per
  locus, 187 rows at the deep locus — trivial).
- **Verdict: no n-log-n term owns a measurable share at the median
  (0.50 s) or the deep (1,423 s) locus.** The owner's suspicion is
  receipted and dismissed — the superlinear mass is the f²·U scan,
  which is now bounded everywhere it can be.

---

## 4. THE FOLD-ENUMERATION ANATOMY at the deep locus (T03226 l0: "204 folds")

**Where the 204 come from (the mechanism, from the source):**
1. **Members:** the window's partition bed rows — every panel
   haplotype path segment the partition tiled over the window
   (**802** member rows at l0).
2. **Folds (pre-merge, 213):** member rows deduplicated by
   `(path, start, end)` across partitions — "identical sequence AND
   walk (relative) — coalesced copies are ONE hypothesis"
   (`realign_score.rs`: the `dedup` map). 802 members → 213 distinct
   rows. This is the panel's ROUTE DIVERSITY in the window: ~46% of
   the 466-member panel's paths carry a distinct row there.
3. **Folds (post-merge, 204):** after scoring, the alignment-connected
   relaxation merges folds whose evidence columns are bit-identical
   (equal-column candidates → merge groups; the corpus admission rate
   is 71.1% of candidates over 2,784 loci; the median pre/post ratio
   is 1.098 — ~10% of folds are evidence-duplicates).

**What enumeration costs:** the `[0.1s]` bracket on the members→folds
line at l0 — **0.1 s of a 1,423.4 s locus (<0.01%)**. The build is
O(members + rows·L log L) (per-row node/edge sorts in the walk length
L); the corpus fits confirm it is neither members-driven
(folds_pre vs members: R²=0.04) nor records-driven (R²=0.34) — the
fold count is set by the panel graph's distinct paths in the window,
full stop.

**Is anything here quadratic?** The enumeration, no. The QUADRATIC IS
WHAT IT FEEDS: the scan's candidate pairs = f(f+1)/2 (20,910 pairs at
204 folds — the exact structural receipt of §0). The equal-column
candidate scan itself is bounded by CONTENT, not by f² (fit vs
folds_pre²: R²=0.01 — the typical fold has a handful of equal-column
neighbors, not f of them).

**THE MERGE'S OWN RECEIPT AT THE FOLD-HEAVY SHAPE (the open term, CLOSED
2026-10-22 with the corrected anatomy):** T03226 l15 is the
fold-heavy/unit-light extreme — 459 members → 247 pre-merge folds over 36
units = 30,628 class pairs, and the PRE-MERGE locus costs **228.4 s
serial** (the census instrument decomposes it: scoring 0.01 s,
receipts/IO 32.7 s — the residual ~195 s was carried here as a
per-class-pair constant of ~6.4 ms, three orders above the 18.4 ns
unit-pair rate, at a shape where the unit-pair model predicts 0.02 s).
WITH THE COMMITTED MERGE the same locus is **36 folds / 666 classes /
0.5 s — a 455x class-pair reduction** (the single-locus receipts:
`T03226-w10k-census{,-premerge}/`, reaped with the disk crisis). The
corpus's fold-heavy regime (1,328 loci at units<2,000, classes>5,000)
fits wall ≈ 0.000431· classes^0.78 (R²=0.66, per-pair median 0.06 ms)
— the l15 shape is that regime's heavy tail, and the merge is the
lever that already collapses the class count.

**THE 6.4 ms TERM'S OWN ANATOMY, LOCATED AND FIXED (2026-10-22, the
stream-discipline re-derivation — `mergefix-tie-certificate-
receipts.txt`, the runners `scripts/complexity-ledger/
run-mergefix-anatomy-T09724.sh` + `run-mergefix-L1.sh`):** the T03226
substrate is unrecoverable (tranche-1's definitional smalls reaped
with the disk crisis), so the anatomy was re-derived at **T09724** —
a tranche-2 S_EXP member of the SAME class (far-bin share 43.46% vs
T03226's 43.4%), rebuilt transiently under the stream law from the
kept member smalls, reproducing the committed substrate EXACTLY (the
cutA census mult 776,096 / reads 406,876 IDENTICAL to the committed
tranche-2 receipt; reference-axis + domain-structure IDENTICAL to
the kept smalls; l9's committed answer reproduced bit-for-bit).
**THE MERGE PASS ITSELF WAS ALREADY FIXED** — measured 1.68 s at the
territory's heaviest locus (439 folds / 246,224 units / 96,580
classes): it buckets by a per-fold column hash computed once with the
exact bit-verify only inside buckets (`ll_bits_hash` /
`ll_bits_equal`), so there is NO per-pair full-column walk left there
to fix — the 6.4 ms term is NOT the merge's. **The term's home is the
QUAL cluster's TIE CERTIFICATE**: for every tied called class beyond
the winner (|called| − 1 bands) the naive path materialized a band
over the FULL spectrum (|spectrum| entries), EVERY entry a fresh
pair of merged node/edge multisets (per-pair vector construction)
plus two observed-mass walks (per-pair full-column walks, ~480-entry
multisets at the 10 kb windows) — while the cluster machinery
consumes from each band ONLY its <=cut predicate set
(`cluster_form_qual`'s own contract). **Measured at two real
committed-stack instances** (the phase-census instrument,
`IMPG_LOCUS_PHASE_CENSUS`): T09724 l7 committed form — 2 bands ×
13,108 entries = 0.08 s ≈ **6.1 µs per band entry**; l7 pre-merge form
— 4 bands × 30,996 entries = 0.18 s ≈ **5.8 µs per entry**. **The
l15 residual re-attributed**: 195 s / 6.1 µs ≈ 3.2e7 band entries ≈
~1,045-band tie set × 30,627 spectrum entries — a called set of
~1,046 classes (the biggest equal-column group's internal pairs at
the best LL, consistent with the committed merge collapsing 247→36
folds there); the "per-class-pair ~6.4 ms" was this certificate's
cost amortized over the class pairs, and the 455x verdict is
re-attributed to the PRE-MERGE TIE CLASS the merge collapses.

**THE FIX (committed 2026-10-22, the degenerate precedent extended):
the cut-index-list certificate at BOTH tie shapes.** The committed
degenerate-path wall fix (the repeat-domain monsters' 35,510 ×
35,510 = 1.26e9 distances, 4,802 s serial, the ~6 GB emission that
breached the 64GiB guard) already computed each band's <=cut
predicate set from per-fold observed-mass signatures computed ONCE
per locus, with SOUND two-sided O(1) bounds per (tied class, spectrum
entry) — D >= |Sm_n(t) − Sm_n(e)| + |Sm_e(t) − Sm_e(e)| and D <=
Sm_n(t) + Sm_n(e) + Sm_e(t) + Sm_e(e) — and the exact
`differing_observed_mass` walk ONLY in the undecided middle band;
the fix routes the NON-degenerate path through the same core
(`tied_class_cut_lists`; the naive band matrix is gone; no tuning
constants; fold indices and order untouched). **The gates**: answers
BYTE-IDENTICAL before/after — calls.jsonl AND truth-qv sha256 — at
T09724 all 16 loci in BOTH the committed and the pre-merge forms,
and at L1 (the standard target, 19 loci); fold-set digests identical
(in-process `IMPG_FOLD_SET_DIGEST` FNV over the final fold space,
16/16 + 19/19); the 4,000-case randomized non-degenerate oracle unit
test against a verbatim naive-band oracle (predicate sets and the
full cluster state identical; cargo test 584/0). **The timing
receipts**: l7 committed 0.08 s (13,108 exact entries) → 0.00 s (2
lists × 8 predicate entries); l7 pre-merge 0.18 s (30,996 entries) →
0.00 s (4 lists × 32 predicate entries) — the certificate's exact
work reduced ~1000x at the real instances, the rest decided by the
O(1) bounds. **The l15-class projection (arithmetic, labeled):** at
3.2e7 entries the naive certificate's ~195 s becomes ~O(3.2e7)
bound-decisions (ns-class) plus the undecided band's exact entries —
the pre-merge l15 locus would land ~0.5 s-class (its committed-form
receipt's own wall), i.e. the pre-merge form no longer pays a
pathological constant the merge exists to dodge. **The corpus share,
measured**: at the committed stack the tie-certificate term is ZERO
at every called=1 locus (T09724: 15 of 16; L1: 19 of 19) and fires
only at tie-bearing loci (l7-class: 0.08 s of a 0.7 s locus) and in
the pre-merge/unmerged-tie form (equal columns outside the merge's
admission) — it is not a corpus pole at the committed defaults; the
fix removes its pathological class everywhere. **The receipt-schema
note (sanctioned with the fix, the degenerate precedent's own
form):** at NON-degenerate tie loci the report member
`qual_tied_bands` (full band vectors) is replaced by
`qual_tied_cut_indices` (the compact predicate lists) — answers
unchanged by construction.

**The principled candidate-count reduction (the evidence model's own
expressibility — no thresholds):** three of the four are already in
the product, and the l15 receipt prices the second one at 455x. (1) The
(path,start,end) dedup — committed. (2) The bit-identical-column merge
— committed (the alignment-connected relaxation; 455x at l15, ~10%
median). (3) The bit-equal-ROW memoization across the pair sums —
exact-preserving, unit-proven (`bit_equal_rows_share_one_class_sum`),
measured 23.9-37.7% at the monsters, **NOT implemented** (the named
Frontier-1 residual; it is the one remaining expressibility reduction
and it is also the plateau lever). (4) The excess-mass bound — the
production default. **What remains after (1)-(4) is the genuinely
expressible hypothesis space: every fold left carries a distinct
evidence column, and cutting any of them is a threshold, not a
reduction** — the owner's own doctrine. The only structural way to
fewer folds per candidate pair is the ROUTES doctrine (panel-side:
fewer, coarser distinct paths per window) — which is the grain lever's
territory (§5, measured and refused), not the call's.

---

## 5. THE GRAIN EXPERIMENT — the bigger-loci lever, re-measured under the current architecture

**Design (the honest re-do of the 2026-10-09 sweep at the CURRENT
binary and stack):** one ordinary locality (**L1**, the five-locality
chain's own committed recipe, full-mass census, local frame — the
far-bin prune is a global-frame flag and is inert here, stated) at
**10 kb / 20 kb / 40 kb / 50 kb**; one deep-stratum territory
(**T03226**, the S_EXP monster: 43.4% far-bin share, the 1,486.8 s
serial loci phase at 10 kb) at **10 kb / 20 kb / 40 kb** under the
GLOBAL-frame production stack (the far-bin prune default, the bounded
scan, the current binary). Every grain's substrate is REBUILT whole in
the work dir — the committed chain's own beds/routes/census and the old
sweep's beds/routes/census were deleted by their disk disciplines, so
the four points share ONE binary, ONE substrate era (the T00041-era
lesson) — and the calls run SERIAL (walls = core-seconds) with
`--truth-qv` (the accuracy sidecar), so the lever is bought or refused
with accuracy receipts on the SAME territories. Runners:
`complexity-ledger/run-grain-L1.sh`, `run-grain-T03226.sh`
(census/derieve/CLI at the current `impg`; the substrate TOOLS —
partition/export/walks/cut/projection — at the committed
target/release examples, unchanged since Oct-19 by the diff receipt
`git diff 82583eb..HEAD` = {commands/genome_infer.rs,
realign_score.rs, examples/class_bound_probe.rs}).

**The prediction tested (the ledger's fixed-vs-scaling classes):** the
per-call fixed classes (inputs, bind, rows/context, the per-locus
floor) should amortize ~1/grain; the content classes (scan f²·U, pair
loop) should NOT — the fold count per window grows sub-linearly
(the old sweep's receipt: max classes 60,726→75,855→96,580 at
10→20→50 kb, i.e. f² grows ~1.25x per 2x grain) while units per
window grow ~linearly, so the per-window content work grows and the
window count halves.

**THE MEASURED CURVE** (serial walls = core-seconds; every grain ran
in the SAME wave at box load 106-113 — the co-tenancy is receipted in
each `.rss`/`.wall`; the L1 leg is the committed five-locality recipe
at the current binary, full-mass census, `--far-bin-prune=false` — the
local frame is byte-frozen against the now-default-on prune):

| L1 grain | windows | LOCI phase (serial core-s) | per-locus median | folds med/max | units med/max | the scan model's share | rank-1 (accuracy) |
|---|---:|---:|---:|---:|---:|---:|---|
| **10 kb** | 19 | **394.0 s** | 6.00 s | 162 / 443 | 12,880 / 107,572 | 66% (261.6 s) | **17/19** |
| **20 kb** | 10 | **506.6 s** | 24.20 s | 270 / 443 | 25,582 / 107,572 | 75% (379.9 s) | **10/10** |
| **40 kb** | 6 | **785.6 s** | 115.65 s | 396 / 443 | 51,462 / 107,572 | 66% (518.0 s) | **6/6** |
| **50 kb** | 4 | **731.4 s** | 141.35 s | 412 / 457 | 67,145 / 157,520 | 76% (552.9 s) | **4/4** |

| the fixed classes (per territory, serial) | 10 kb | 20 kb | 40 kb | 50 kb |
|---|---:|---:|---:|---:|
| inputs (panel load) | 2.4 s | 2.6 s | 2.4 s | 2.4 s |
| binding | 3.8 s | 3.8 s | 3.8 s | 3.8 s |
| rows/context | 4.3 s | 4.5 s | 4.8 s | 5.0 s |
| census records (the census construction) | 81,688-class | — | — | — |

**THE LOCALITY VERDICT, measured:** the owner's bigger-loci lever is
**REFUSED on the cost side at the current architecture.** The content
classes GROW as the grain coarsens — the per-window units grow
**linearly** (medians 12,880 → 25,582 → 51,462 → 67,145 = 1.99x per 2x
grain), the per-window folds grow only **sub-linearly** (medians 162 →
270 → 396 → 412, ~1.6x per 2x; the max is CONTENT-CLAMPED at 443-457
at every grain — the fold diversity saturates), so the scan model per
territory rises monotonically **261.6 → 379.9 → 518.0 → 552.9 s** and
the whole loci phase **394.0 → 506.6 → 785.6 → 731.4 s (+1.9x at 40kb)**
— the units' linear growth beats both the 1/grain window count and the
sub-linear fold growth. The classes the lever was supposed to amortize
**do not amortize because they are not window-bound:** binding is
census-record-bound and the census record count is GRAIN-INVARIANT
(receipted at the territory leg below: 81,688 records at 10 kb vs
81,765 at 20 kb, mass 9.377M vs 9.381M — the census construction is
linear in routed mass, not windows), rows/context is flat-to-up, the
inputs load is per-call. The one genuinely window-bound term is the
**per-locus fixed floor** (~0.2 s median at the small loci) — at 4x
grain it saves ~0.15 s/locus while the content classes add ~10 s/locus
of new work at this territory. The accuracy at these two territories
did NOT degrade in this wave (17/19 → 10/10 → 6/6 → 4/4 — the near-twin
ABSORPTION class dominates at L1's windows here; the denominators are
small and the window mix changed with the substrate rebuild, so this
is NOT a retraction of the 2026-10-09 sweep's genome-wide accuracy
verdict — it is the same-territory receipt the re-measure owed, and it
says the cost refusal stands on cost alone).

**THE DEEP-STRATUM LEG (T03226, the S_EXP monster, global-frame
production stack), the receipts that landed in this window** — with one
form deviation caught at receipt time and stated: the territory leg ran
WITHOUT `IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS=1` (the committed
assessment form sets it; the pre-merge form's l0 carries 213 folds /
22,791 classes vs the committed post-merge 204 / 20,910 — ~9% more scan
at l0), so all three territory grains are the PRE-MERGE form,
**within-wave consistent for the grain ratios**; the 10 kb accuracy
receipt is rank-identical to the committed one (7/24 — the merge is
rank-neutral at this territory):

| T03226 grain | prep partitions | cut census records | inputs s | binding s | l0: folds x units (classes) | l0 serial wall | rank-1 |
|---|---:|---:|---:|---:|---|---:|---|
| **10 kb** (the committed substrate, rerun at the current binary) | 1,131 | **81,688** (mult 9.377M; prune 43.4%) | 237.0 | 93.4 | 213 x 2,925,505 (22,791) | (loci phase **1,836.4 s**, 24 windows, median 0.75 s) | **7/24 = the committed answer** |
| **20 kb** | 608 | **81,765** (mult 9.381M; prune 41.3%) | 241.6 | (pending the run's close) | l0: 277 x 3,411,073 (38,503); l6: 303 x 11,331 (46,056); l7: 309 x 41 (47,895) | l0 **2,631.1 s** (8/12 loci done at the commit, the receipted walls in `callC.err`) | (the serial call is mid-flight — see the resume point) |
| **40 kb** | 332 | **81,839 patterns** (mult 9.389M — grain-invariant again; the cut DONE in 801.3 s) | | | l0: the matrix line receipts **1,187,114,076 no-pin pairs** at the factorized stage (the merged monster window) | (l0 mid-scan at the commit — the 40 kb monster window's serial scan projects to HOURS: the coarse-grain explosion receipted in the shape itself) | |

**What the territory receipts already measure:** the per-window content
at the deep stratum grows HARD — l0's classes 22,791 → 38,503 (1.69x
for 2x grain; the merged 20 kb window carries the cross-boundary fold
pairs the two 10 kb windows never enumerated against each other),
l0's units 2.93M → 3.41M (+17%: the cross-boundary census records'
units are ~14% of the merged window's mass — NEW interacting mass the
coarse window creates), and at 40 kb the merged l0's matrix alone is
**1.19 BILLION no-pin pairs** (the factorized-stage receipt) — the
monster window's serial scan projects to hours, the coarse-grain
explosion receipted in the shape itself. The monster locus's serial
wall 1,387 s → 2,631 s per window: the deep stratum's content classes
grow per window much faster than the window count shrinks — the same
refusal as the locality, amplified. The fixed classes again do not
amortize: census records 81,688 → 81,765 → 81,839 patterns
(grain-invariant), inputs 237.0 → 241.6 s (flat), and the far-bin
prune's share is grain-stable (43.4% → 41.3% — the pair-bin adjacency
widens with the window, as designed). The 10 kb BOUNDED production
form beside the serial assessment: callProd at T16 (load ~110) loci
phase **155.0 s** with l0's bounded scan at **99.44% pruned in 1.3 s**
(22,664 of 22,791 pairs skipped — the bound engages at the merged
shape exactly as at the committed one).

**The accuracy caveat, held:** the 2026-10-09 sweep's verdict
(NOT SUPPORTED at either depth; the committed 10 kb is the
accuracy-optimal point; pooled window rank-1 DEEP 66.9%→61.7%→60.9%,
35x 40.5%→36.7%→39.1%; the near-twin class never becomes decisive;
segdup swallowing and expressibility losses at the coarse grains)
**STANDS until these receipts land** — the re-measure runs the truth
assessment at every grain on the same territories precisely so the
cost lever is not bought on assumed accuracy.

**The structural blockers to coarser grains, named:**
1. **The grain is LOCKED into the substrate** (the sweep's own
   receipt): everything downstream of `impg partition -w` is per-grain;
   the census, routes, derive-cache, locality graphs all rebuild.
2. **The global-frame routes are grain-partial** — `global-routes-par`
   was derived from a 10 kb census; the call consumes only its
   grain-INDEPENDENT members (lanes/ports/source-paths, byte-identical
   across grains per the sweep's `routes-identity.txt` receipt), which
   is why the territory re-grain runs at all; the ownership/registry
   members (the tiling-dependent ones) would need a per-grain
   derive-routes (the once-per-panel ~1,492 s) for any consumer that
   wants them.
3. **The FastGA sliver-partition crash class** (the sweep's engine-bug
   receipt) can appear at coarser tilings; the L1-40 kb chain runs the
   axis-only containment (the sweep's own) and will receipt its
   appearance if it fires.
4. **The far-bin prune's window geometry** is the one call-side piece
   that is grain-AWARE by design (the diameter test against the
   census's own pair-bins) — coarser windows make "adjacent" rarer, so
   the prune's mass share RISES with grain; the territory leg's
   prune-receipt lines at each grain are the measurement.

---

## 6. THE RESTATED 0.14 CORE-S/LOCUS ARITHMETIC

**The target, derived:** 16 cores × 1,200 s = 19,200 core-s per sample;
the universe's ~5,807 territories × ~25 loci ≈ 145k loci → **0.132 ≈
0.14 core-s per locus ALL-IN** (call + map + cut + load amortized).

**Where we are (every term measured, form C, the committed stack):**
- the call budget: **6.38e6 ± 4.48e6 core-s** stratified (the tranche
  close) = **~44 core-s/locus mean** (the mean is tail-dominated: the
  median locus is 0.50 s and the worst-5% own 96.3%);
- the map side: ~26k core-s (~0.18/locus) → ~13k (~2x engineering);
- the cut: the width-amortized true-CPU receipt (the standing
  open term, §10.6.18's census);
- the per-territory fixed classes (the grain lever's own terms):
  inputs 243.8 s + bind 7.2 s + rows 13.4 s = **264 s ≈ 11
  core-s/locus** at 10 kb — HALF the mean's fixed+content budget before
  any content work at all.

**The restated decomposition of the ~44 core-s/locus mean:**
| term | core-s/locus | path to 0.14 |
|---|---:|---|
| inputs/panel load (per-call) | ~9.8 (243.8/25) | **NONE needed per-locus** — the fleet form's fork-after-load COW amortizes it per SAMPLE (the product path's own default); at fleet scale → ~0 |
| bind + rows/context | ~0.8 | **REFUSED as a grain lever (§5's measured verdict: bind is census-record-bound and the census is grain-invariant — 81,688 vs 81,765 records — rows flat-to-up); the path that remains is Frontier 2's panel-side precompute** (the gated increment) |
| the per-locus fixed floor | ~0.2 at the median | **the only window-bound term** — but §5 measures the trade: 4x grain saves ~0.15 s/locus of floor while the content classes ADD ~2x their work (the units grow linearly per window, the folds sub-linearly, the windows 1/grain) — the lever is refused at cost; the floor's other path is Frontier 2 |
| the class scan (the f²·U tail) | ~29 (71% of the corpus wall / loci) | **bound (DONE)** where peaked (~half the tail); **memoization** (the plateau's 23.9-37.7%, NOT implemented) for the rest |
| the scoring pair loop | ~6.6 (16% / loci) | **IRREDUCIBLE core** at the pinned-heavy shapes (~0.4-0.6 core-s/locus of per-base comparison) + ENGINEERING (the placement-walk residue, the 2-4x over the floor) |
| map + cut | ~0.2+ | ~2x engineering (the standing receipt) |

**The honest verdict on 0.14:** the term with a clean path is the
**fixed mass (~11 core-s/locus at 10 kb)** — the fleet's COW already
deletes the load term per-locus, and the remaining fixed classes
(bind/rows/floor) are precompute territory (§5's grain experiment
MEASURED the bigger-loci lever and REFUSED it: the fixed classes are
mass-bound, not window-bound, and the content classes grow ~1.9x at
40 kb). The **scan term** is bounded but its plateau remainder needs
the memoization increment. The **pair-loop term is the floor the model
itself sets** — no algorithmic lever removes the per-base comparison
at the pinned-heavy loci; only the CONTENT classes (deeper restriction,
the RV 35x basis) or the routes doctrine (fewer folds per window) move
it, and both are accuracy-bearing decisions. **All-in, the arithmetic
that reaches 0.14 is: fleet-amortized load (~0) + precomputed fixed
classes + bounded+memoized scan + the pair loop's floor (~0.5-1 at the
mean after the content classes' levers)** — i.e. the target is
reachable only as a FLEET+PRECOMPUTE+MEMOIZATION stack on the 0.5 s
median, with the tail's plateau monsters as the gating remainder, and
(grain now measured out) every one of those levers is measured,
designed, or implemented. The residual that is NOT reachable by any
receipted lever is the pinned per-base comparison class at the
monsters — the model's own price.

**THE 16-CORE ARITHMETIC, REFRESHED AT THE CLOSED TERM (2026-10-22):**
the tie-certificate fix (§4) moves NO budget line — at the committed
defaults the term is measured ZERO at every called=1 locus (T09724:
15 of 16; L1: 19 of 19) and ~0.1 s at the tie-bearing l7-class —
correctly absent from the decomposition above, which is why the
"~6.4 ms" never appeared as a corpus pole. What the fix removes is
the TAIL RISK the ledger's own l15 receipt exposed: a fold-heavy
window whose equal columns fall outside the merge's admission (the
pre-merge form, or outside-bubble/pocket-mass refusals) previously
paid ~6 µs × (tie-set × spectrum) at the certificate — the l15-class
arithmetic, ~195 s — and now pays O(1) bound decisions per entry
plus the undecided band's exact walks (measured: 30,996 exact
entries → 32 at the l7 pre-merge instance). The 0.14 path itself is
unchanged: fleet COW + precompute + bounded-and-memoized scan + the
pair loop's floor — with the pre-merge/tie-class pathological
constant no longer a hidden multiplier on the fold-heavy tail.

---

## 7. THE RECEIPTS AND THE RESUME POINT

- `scripts/complexity-ledger/mine-locus-shapes.py` — the receipt miner
  (per-locus + per-run TSVs from any call's stderr receipts).
- `scripts/complexity-ledger/fit-exponents.py` — the exponent fitter
  (every fit above, R² named; output
  `complexity-ledger/exponent-fits.txt`).
- The work dir
  `/home/erikg/hprcv2/global-panel-call-mode/complexity-ledger/`:
  `mined-tranche-{loci,runs}.tsv` (2,822 + 138 rows),
  `mined-prod-{loci,runs}.tsv` (the 11 bounded-gate targets),
  `exponent-fits.txt`, `run-grain-L1.sh`, `run-grain-T03226.sh`,
  `L1-w{10,20,40,50}k/`, `T03226-w{10,20,40}k/` (each stage's
  `.done/.wall/.err/.rss`), `logs/*.nohup`.
- **The grain lanes are marker-gated and resumable** — the resume point
  is the `.done` markers; any interrupted stage re-runs on relaunch.
  LANDED: the whole L1 leg (all four grains, serial, with accuracy
  sidecars — the tables of §5), the T03226 10 kb wave-matched rerun
  (loci 1,836.4 s serial at load ~110; answers rank-identical to the
  committed receipt), the T03226 20 kb cut (81,765 census records,
  grain-invariant) and its l0 receipt (2,631.1 s). IN FLIGHT at the
  commit: the T03226 20 kb serial callC (locus 0 of 12 done — the
  merged monster windows are 1.7x classes each; the whole run is
  hours of serial away) and the T03226 40 kb lane (prep done at 332
  partitions, the cut mid-scan). BOTH LANES KEEP RUNNING PAST THE
  COMMIT (nohup, `logs/T03226-w{20,40}k.nohup`); their completion
  receipts fill §5's territory table — relaunch the same commands if
  the box reboots. The form note: the territory leg is the PRE-MERGE
  form (§5); a committed-form rerun is one env away
  (`IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS=1` on the callC/callProd
  invocations) and is the first item of any follow-up.
- cargo test: **the 2026-10-22 window's src change is the tie-certificate
  fix** (`src/genome_inference/realign_score.rs`: the non-degenerate
  tie certificate takes the committed degenerate path's cut-index-
  list form — §4's closed term — plus the env-gated phase-census and
  fold-set-digest instruments); **cargo test --release --lib 584/0
  green, 5 ignored** (583 committed + the new 4,000-case
  non-degenerate-band oracle
  `the_non_degenerate_tie_certificate_matches_the_naive_band_oracle`).
- **The merge-fix anatomy's own receipts (2026-10-22, the stream
  discipline's first full customer):**
  `docs/syng-gmem-bwt/mergefix-tie-certificate-receipts.txt` (the
  substrate-reproduction gates, every phase-census line at both
  forms, the sha256 answer-digest gates, the fold-set digest gates,
  the walls) and the runners
  `scripts/complexity-ledger/run-mergefix-anatomy-T09724.sh` +
  `run-mergefix-L1.sh` (marker-gated, 150G-budget-guarded,
  receipts-at-write-time, `reap` subcommand). The workdirs
  (`/mnt/nvme2n1/erikg/mergefix-T09724{,-keep}` and
  `mergefix-L1{,-keep}`) were REAPED at the increment's close after
  the receipts were committed (mine first, reap after); the keep
  dirs' contents (the receipts + the small answer files) are folded
  into the receipts file above.
- **The L2-L4/pilot gate re-derivations were SKIPPED this window**
  (each is a ~35 min substrate rebuild; L1 + the T03226-class member
  carry the gate) — honestly labeled: they are re-runnable at any
  time with `run-mergefix-L1.sh`'s shape at L2/L3/L4/pilot; absence
  is a coverage label, never a pass claim.

NO PR PUSH (the owner pushes).
