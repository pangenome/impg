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

**THE FIXED CLASSES ARE REAL AND GRAIN-AMORTIZABLE — the owner's
bigger-loci lever has a confirmed target:** the per-call inputs phase
(panel load) fits windows^-0.035 with R² = 0.016 — **FLAT** (median
243.8 s per territory regardless of window count), binding median
7.2 s, rows/context median 13.4 s, and a per-LOCUS floor of
**~0.20 s median** (measured at the 194 loci whose matrix is < 2e4
pairs — fold build + qual + receipts + the scan's floor). At the
median locus (wall 0.50 s) this floor is 40% of the locus; at the
p95 locus (28.2 s) it is <1%.

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

**The principled candidate-count reduction (the evidence model's own
expressibility — no thresholds):** three of the four are already in
the product. (1) The (path,start,end) dedup — committed. (2) The
bit-identical-column merge — committed (the alignment-connected
relaxation). (3) The bit-equal-ROW memoization across the pair sums —
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
territory (§5), not the call's.

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

**THE MEASURED CURVE** (serial walls, the wave's co-tenancy stated in
each `.wall` receipt; the box carried load 106-113 from the tranche
fleet throughout — every grain ran in the SAME wave):

*(the walls land in the work dir as the lanes complete; the table
below is filled from the `.wall`/phase-summary receipts at the
commit — see the resume point for any cell still in flight)*

| grain | L1: substrate chain | census | CLI (serial, core-s) | per-window folds (max) | rank-1 / loci (the accuracy bar) |
|---|---|---:|---:|---:|---|
| 10 kb | (s0-s6 walls in `L1-w10k/`) | | | | |
| 20 kb | | | | | |
| 40 kb | | | | | |
| 50 kb | | | | | |

| grain | T03226: prep+cutA | callC (serial, core-s) | loci | per-window folds (max) | rank-1 / loci |
|---|---|---:|---:|---:|---|
| 10 kb | (the committed substrate; the tranche receipt 2,820.7 s) | | | | |
| 20 kb | | | | | |
| 40 kb | | | | | |

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
| bind + rows/context | ~0.8 | **grain** (halve per window merge; the fit says bind ∝ census-records^0.95, and records merge with windows) + **Frontier 2 precompute** (the panel-side per-locus precompute's own gated increment) |
| the per-locus fixed floor | ~0.2 at the median | **grain** (amortizes ~1/grain) |
| the class scan (the f²·U tail) | ~29 (71% of the corpus wall / loci) | **bound (DONE)** where peaked (~half the tail); **memoization** (the plateau's 23.9-37.7%, NOT implemented) for the rest |
| the scoring pair loop | ~6.6 (16% / loci) | **IRREDUCIBLE core** at the pinned-heavy shapes (~0.4-0.6 core-s/locus of per-base comparison) + ENGINEERING (the placement-walk residue, the 2-4x over the floor) |
| map + cut | ~0.2+ | ~2x engineering (the standing receipt) |

**The honest verdict on 0.14:** the term with a clean path is the
**fixed mass (~11 core-s/locus at 10 kb)** — the fleet's COW already
deletes the load term per-locus, and the grain lever attacks the bind/
rows/floor classes at ~1/grain (§5's curve measures exactly how much).
The **scan term** is bounded but its plateau remainder needs the
memoization increment. The **pair-loop term is the floor the model
itself sets** — no algorithmic lever removes the per-base comparison
at the pinned-heavy loci; only the CONTENT classes (fewer units:
deeper restriction, the RV 35x basis) or the routes doctrine (fewer
folds per window) move it, and both are accuracy-bearing decisions.
**All-in, the arithmetic that reaches 0.14 is: fleet-amortized load (~0)
+ grained fixed classes (~11 → ~3-5 at 20-40 kb IF the accuracy holds)
+ bounded+memoized scan + the pair loop's floor (~0.5-1 at the mean
after the content classes' levers)** — i.e. the target is reachable
only as a FLEET+GRAIN+MEMOIZATION stack on the 0.5 s median, with the
tail's plateau monsters as the gating remainder, and every one of
those three levers is now measured, designed, or implemented. The
residual that is NOT reachable by any receipted lever is the pinned
per-base comparison class at the monsters — the model's own price.

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
  The T03226-w20k cut completed its scan (90,378 touching patterns,
  5,646,738 kept reads — the receipts in `T03226-w20k/cutA/`) and the
  serial callC follows.
- cargo test: NO src changes in this window (instruments + analysis
  only; the two scripts are new files) — the committed 583/0 green
  stands.

NO PR PUSH (the owner pushes).
