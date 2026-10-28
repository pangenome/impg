# Unit-level redundancy of the T00048 `ll_matrix`: how much of the pair loop's inner term is droppable, and whether the drop is exact

**Subject.** `/tmp/probe-fast/T00048/instrument-receipt.jsonl.ingredients.jsonl`
(795,081,652 B) — the ingredients receipt of one `impg infer call-region`
run (`manifest.json`: `anchor-realign-diplotype-calls-v1-farbin-prune`,
24 loci of the T00048 region, elapsed 382.2 s, shared panel load 198.7 s).
Sibling evidence used only for grounding: `instrument-receipt.jsonl`
(the per-locus receipt: `scoring`, `class_log_likelihoods`, `phase_census`).

**Answer in one paragraph.** At the deep locus (locus 8: 429 folds ×
45,423 units = 19,486,467 cells) the unit dimension is dominated by one
exact redundancy: **40,777 of the 45,423 units (89.8%) have a
bit-constant ll column** — every fold holds the same value, and that value
is the elsewhere constant `E = -14.901063413333642`. Collapsing the
whole matrix to its **842 distinct bit-identical columns** (841
non-constant vectors + 1 constant bucket) shrinks the pair loop's inner
term **53.95×** (4,189,590,405 → 77,661,870 mix-term evaluations, over the
92,235 unordered class pairs), and **60.2×** composed with the
already-implemented bit-equal-*row* memo. Most of that already exists in
the non-exact `sparse-skip` path (it walks the pinned union, ~1,892
terms = 24×), so the genuinely new lever is the last **2.25×**. The
reduction is **exact in real arithmetic — no approximation, the same addends — but it is NOT
bit-identical**: it removes addends from the middle of a fixed
left-to-right IEEE-754 accumulation and merges duplicate products, so
round-to-nearest sees a different association. Measured against the
emitted committed sums, **65,520 of the 65,535 kept class pairs change**,
by up to **1,362 ulp (1.59e-7)**. It therefore cannot produce
`calls.jsonl`-identical output from the exact/default path; it is a
`sparse-skip`-class, must-be-labeled mode.

---

## 1. The schema, as found (and one correction to the brief)

The brief described "~10.3M JSONL rows of (fold, unit) cells". **The file
has 24 rows — one per locus.** Each row carries the locus's *whole*
`folds × units` matrix, so the (fold, unit) cells number
**32,283,672**, not 10.3M, and the "deep locus" is a single 438.6 MB
line (55% of the file), not a bulk of rows.

One row (one locus), keys in emission order:

| key | shape | meaning |
|---|---|---|
| `folds` | array[F] | the locus's haplotype classes (POA folds): `length`, `sequence`, `members[{path_name,start,end}]`, `nodes`, `edges`, `walk` |
| `ll_matrix` | f64[F][U] | **the subject**: row = fold, column = unit; entry = that fold's log-likelihood term at that unit. Written in the dense form: an unpinned entry is the no-pin constant `E` |
| `locus` | int | locus index (= line index = `partition` in this run) |
| `partition` | int | partition index |
| `units` | array[U] | one entry per **column**: `count` (read multiplicity), `mirror`, `read` (sequence), `record`, `variant`, `w_lo` |

Verified by construction against `src/genome_inference/realign_score.rs`
(the streaming emitter, `write_ll_matrix_json` at :5763, the ingredients
spill at :11283-11330) and by re-deriving the receipt's own counters from
the file (§7). Note what is *not* in the ingredients row: **`E`**. It
lives in the sibling receipt's `scoring.elsewhere_log_prob`
(= `-14.901063413333642` at every locus of this run; read length 148,
uniform qualities). The unit's own `count` is exactly the `f64` the scan
uses — `unit_counts: Vec<f64> = units.iter().map(|u| u.count as f64)`
(:9722).

**What the units are for** (`ScoreMatrix::class_pair_sum`, :5652): the
class scan's inner term is a fixed left-to-right walk over units,

```
total = 0.0
for u in 0..U ascending:            # unit order is the accumulation order
    total += count[u] * mix_logsumexp(ll[first][u], ll[second][u])
```

executed for every unordered class pair (`F*(F+1)/2`, hom pairs
included). "Redundant units" therefore means: columns that let us *drop
addends* or *merge addends* without changing the real-valued sum.

## 2. The two reductions, defined exactly as measured

Grouping is by the **exact f64 bit pattern of the whole column**
(`np.void` view compare — the same predicate the repo's `ll_bits_equal`
uses; nothing is rounded or hashed-loose).

1. **Constant-vector units.** Column bits all equal across the locus's
   folds ⇒ the unit adds `count[u] * mix(v,v)` to *every* pair, and
   `mix(v,v) = v` **bit for bit** for every value in the file
   (`scripts/unit-redundancy/mix_self_check.py`: 24,392/24,392 distinct
   values, 0 violations — `exp(0)=1`, `ln(1)=+0.0` are exact, so
   `m + 0.0 = m`). One precomputed scalar `S = Σ_u count[u]*v` replaces
   the whole class.
2. **Identical-vector units.** Units whose full columns are bit-identical
   are one term with the summed count mass: `mass_g * mix(rep_a, rep_b)`,
   `mass_g = Σ_{u∈g} count[u]` (integer-exact in f64).

Per locus I report `N` = distinct columns, `M` = number of distinct
*constant* vectors (each distinct constant **value** is its own vector,
so `M ≤ #constant units`), and the reduced inner length
`N − M + 1` (non-constant vectors + one scalar bucket). Hence

```
X = U / (N − M + 1)
```

At every locus of this file `M = 1` and that single constant value **is
`E`**, i.e. the constant-vector units are exactly the fully-unpinned
units. That matters for placement (§8): reduction 1 is *already*
implemented, but only in the non-exact `sparse-skip` mode's
`constant_mass` term.

## 3. Per-locus table, all 24 loci, sorted by cell count

`cells = F*U`; `%U` = fraction of units with a constant vector;
`const count-mass` = `Σ count[u]` over those units (exact integer);
`X` = the unit-term shrink; `distinct rows` = bit-equal *fold rows*
(the existing row memo) and `X w/ row memo` = the composed shrink of the
whole `pairs × units` term.

| locus | F (folds) | U (units) | cells F*U | N (distinct vecs) | M (const vec-groups) | const units | % U | const count-mass | non-const units | non-const vecs | X = U/(N-M+1) | distinct rows | X w/ row memo |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 8 | 429 | 45423 | 19486467 | 842 | 1 | 40777 | 89.8 | 43059 | 4646 | 841 | 53.95 | 406 | 60.22 |
| 15 | 249 | 13612 | 3389388 | 251 | 1 | 12974 | 95.3 | 13501 | 638 | 250 | 54.23 | 241 | 57.88 |
| 10 | 190 | 14174 | 2693060 | 136 | 1 | 13845 | 97.7 | 22494 | 329 | 135 | 104.22 | 121 | 256.21 |
| 16 | 345 | 6268 | 2162460 | 268 | 1 | 5698 | 90.9 | 5951 | 570 | 267 | 23.39 | 314 | 28.23 |
| 14 | 143 | 13784 | 1971112 | 265 | 1 | 12858 | 93.3 | 13416 | 926 | 264 | 52.02 | 137 | 56.65 |
| 5 | 353 | 3058 | 1079474 | 178 | 1 | 2655 | 86.8 | 2774 | 403 | 177 | 17.18 | 321 | 20.77 |
| 20 | 241 | 826 | 199066 | 348 | 1 | 157 | 19.0 | 162 | 669 | 347 | 2.37 | 220 | 2.85 |
| 19 | 343 | 495 | 169785 | 181 | 1 | 147 | 29.7 | 150 | 348 | 180 | 2.73 | 328 | 2.99 |
| 4 | 109 | 1386 | 151074 | 241 | 1 | 754 | 54.4 | 808 | 632 | 240 | 5.75 | 91 | 8.24 |
| 23 | 143 | 910 | 130130 | 230 | 1 | 155 | 17.0 | 158 | 755 | 229 | 3.96 | 134 | 4.50 |
| 22 | 134 | 964 | 129176 | 211 | 1 | 512 | 53.1 | 531 | 452 | 210 | 4.57 | 124 | 5.33 |
| 7 | 113 | 1124 | 127012 | 228 | 1 | 655 | 58.3 | 687 | 469 | 227 | 4.93 | 98 | 6.55 |
| 1 | 155 | 716 | 110980 | 245 | 1 | 139 | 19.4 | 142 | 577 | 244 | 2.92 | 140 | 3.58 |
| 3 | 143 | 527 | 75361 | 202 | 1 | 32 | 6.1 | 33 | 495 | 201 | 2.61 | 133 | 3.01 |
| 13 | 154 | 487 | 74998 | 201 | 1 | 14 | 2.9 | 16 | 473 | 200 | 2.42 | 138 | 3.02 |
| 12 | 166 | 396 | 65736 | 158 | 1 | 38 | 9.6 | 40 | 358 | 157 | 2.51 | 146 | 3.24 |
| 11 | 99 | 621 | 61479 | 162 | 1 | 191 | 30.8 | 197 | 430 | 161 | 3.83 | 85 | 5.19 |
| 0 | 93 | 634 | 58962 | 261 | 1 | 23 | 3.6 | 23 | 611 | 260 | 2.43 | 88 | 2.71 |
| 6 | 119 | 452 | 53788 | 174 | 1 | 40 | 8.8 | 40 | 412 | 173 | 2.60 | 105 | 3.33 |
| 18 | 124 | 401 | 49724 | 202 | 1 | 37 | 9.2 | 37 | 364 | 201 | 1.99 | 115 | 2.31 |
| 21 | 101 | 380 | 38380 | 123 | 1 | 11 | 2.9 | 11 | 369 | 122 | 3.09 | 76 | 5.44 |
| 2 | 124 | 35 | 4340 | 1 | 1 | 35 | 100.0 | 37 | 0 | 0 | 35.00 | 1 | 271250.00 |
| 17 | 40 | 43 | 1720 | 2 | 1 | 42 | 97.7 | 44 | 1 | 1 | 21.50 | 3 | 2938.33 |
| 9 | 36 | 0 | 0 | 0 | 0 | 0 | 0.0 | 0 | 0 | 0 | — | 36 | — |

Totals over the file: 32,283,672 cells, 106,716 units, 5,704,298,994
full `pairs × units` terms → 160,828,194 with the unit reduction (35.47×)
→ 139,035,890 with the row memo added (41.03×).

Reading of the table:

* Redundancy is **bimodal**, not uniform. Six heavy loci (8, 15, 10, 16,
  14, 5) are 87–98% constant-vector units and give 17–104×; the thin loci
  (0, 3, 13, 21 …) are 3–9% constant and give 2.0–3.1×, where the
  identical-vector merge is the only thing worth having.
* Locus 9 is **degenerate**: 36 folds, an empty unit set, zero cells
  (the receipt marks `phase_census.degenerate` handling elsewhere).
* Locus 2 is a **fully constant matrix** (35 units, one distinct column,
  one distinct row): the reduced form has one term.
* The row memo (`distinct rows` column) is weak exactly where the unit
  dimension is strong: 429 → 406 rows at locus 8 is 1.12×, i.e. 24× less
  than the unit reduction. The two are orthogonal and compose (§6).

## 4. The deep locus in full (locus 8, line 8, 438,584,485 B)

```
folds F = 429, units U = 45,423, cells = 19,486,467
distinct columns N = 842
constant vector-groups M = 1 (units = 40,777, = 89.77% of U,
                              count-mass = 43,059, value = E exactly,
                              -14.901063413333642, bit-identical to
                              scoring.elsewhere_log_prob)
non-constant: 4,646 units in 841 distinct vectors
```

Unit-count texture: 43,223 units have `count` 1, 1,931 have 2, 170 have 3,
52 have 4, …, max 16 (count-mass 48,101). Inside the constant class:
38,893 units at count 1, 1,655 at 2, …, max 16.

Multiplicity distribution of the 842 distinct vectors (units sharing one
bit-identical column):

```
  1×: 456 vectors      2×: 154     3×: 62     4×: 41     5×: 33
  6×: 15    7×: 7      8×: 5       9×: 7      10×: 7     11×: 4
 12×: 5    13×: 2     14×: 4      15×: 3     17×: 6     18×: 2
 19×: 1    23×: 3     24×: 1      25×: 1     26×: 2     27×: 1
 28×: 1    29×: 1     30×: 1      31×: 1     34×: 1     37×: 2
 39×: 1    48×: 1     52×: 1      66×: 1     73×: 1     74×: 1
 79×: 1   156×: 1    165×: 1     251×: 1    272×: 1    969×: 1
40777×: 1   (the constant class)
```

Largest non-constant twin groups (units, count-mass, first unit index):
969/1,046/4,327 · 272/275/26,711 · 251/258/5,755 · 165/170/39,340 ·
156/160/26,698 · 79/81/40,461. 456 of the 841 non-constant vectors are
singletons; the other 385 cover 4,190 units.

Density of the matrix: only **604,977 of 19,486,467 cells (3.10%)** differ
from `E`; 89.77% of all cells sit in constant columns; the bit-distinct
column basis is 842 × 429 = **361,218 cells (1.85% of the matrix)**.

## 5. **The headline arithmetic (deep locus)**

```
U = 45,423 units,  F = 429 folds,  pairs (incl. hom) = 429*430/2 = 92,235

full unit term        : 92,235 * 45,423 = 4,189,590,405  (count*mix terms)
reduced unit term     : N - M + 1 = 842 - 1 + 1 = 842 terms per pair
                      : 92,235 *      842 =    77,661,870  terms

X = 4,189,590,405 / 77,661,870 = 45,423 / 842 = 53.9466  ->  53.95x
```

Factored by cause:

```
constant-vector drop alone : 45,423 -> 4,646      =  9.78x
identical-vector merge beyond it : 4,646 -> 841   =  5.524x   (9.78*5.524 = 54.0)
scalar added back per pair : 841 + 1 = 842        (one term, S)

scalar S = sum_u count[u]*E over the 40,777 constant units
         = -641,624.8895152748 (accumulated in unit order)
         NOTE: the single multiply 43,059 * E = -641,624.8895147332,
         i.e. 5.4e-7 apart from the accumulated scalar. Even the scalar
         has an association, and it is visible.
```

**Against what the tree already does (not just against the dense walk).**
The existing non-exact `ClassScanFastIndex` walks the union of the two
rows' *pinned* (≠ `E`) units, which at this locus is already short —
measured over 3,000 random pairs (`union_probe.py`): mean **1,891.6**,
median 1,656, p95 3,157, max 4,276, i.e. **24.0×** off the committed
dense walk. So:

* reduction 1 alone (45,423 → 4,646) is **dominated** by the existing
  sparse-skip (→ ~1,892): never do only that;
* the genuinely new part is the **identical-vector merge**, which takes
  the inner term to 842 terms — a further **1,892 / 842 = 2.25×** on top
  of the fast path, and 53.95× against the committed walk that the exact
  mode and all the bit-identity gates use.

Composed with the row memo already in the tree (429 folds → 406 bit-equal
rows, 82,621 canonical pairs — the receipt's own `phase_census.scan_memo`
reports `distinct_rows: 406, canonical_pairs: 82621`):

```
full                      : 4,189,590,405
row memo only             : 82,621 * 45,423 = 3,752,835,483   (1.116x)
row memo + unit reduction : 82,621 *    842 =    69,566,882   (60.22x)
```

For scale, the same shrink at the other heavy loci: locus 10 = 104.2×,
locus 15 = 54.2×, locus 14 = 52.0×, locus 16 = 23.4×, locus 5 = 17.2×.

## 6. Exactness: **not bit-identical** — and here is the measurement

Two different questions, two different answers:

* **Is it exact in real arithmetic?** Yes, by construction. The reduced
  form sums the *same* real-valued addends
  (`count[u] * mix(ll[f][u], ll[g][u])`, resp. `mass_g * mix(rep,rep)`
  with `mass_g = Σ count`, exact in f64) and only one exactness lemma is
  needed beyond that: `mix(v,v) = v` bit for bit, verified over all
  24,392 distinct values in the file. No term is approximated, dropped,
  or thresholded; no underflow/overflow shortcut is taken.
* **Is it bit-identical to the committed emission (⇒ same `calls.jsonl`
  bytes)?** **No.** The committed sum is a fixed left-to-right
  accumulation; the reduction is a different association of the same
  addends, and IEEE-754 addition is not associative. Concretely at locus
  8 the constant units are **interleaved** through the whole unit order
  (first at 0, last at 45,422, first non-constant at 1,653, spread over
  1,020 runs), so every later partial sum changes; and the merge replaces
  `count1*x + count2*x` by `(count1+count2)*x`, which is not the same
  rounding.

### 6.1 Deep locus, anchored on the emitted values

Ground truth = the receipt's `class_log_likelihoods` for locus 8, which
is *the committed unit-order accumulation itself* at every kept
(bounded-scan) pair — 65,535 of 92,235 pairs are kept, the rest are
`null`.

| arm | pairs | bit-mismatched vs emitted | max ulp | max abs delta |
|---|--:|--:|--:|--:|
| committed replica (gate) | 400 (random sample) | **0** | **0** | 0 |
| V1: constant drop only | 4,000 | 4,000 (100%) | 1,383 | 1.610e-7 |
| V2: identical-vector merge only | 16,384 | 16,384 (100%) | 5,828 | 6.785e-7 |
| **V3: both (the 53.95× form)** | **65,535 (all kept)** | **65,520 (99.98%)** | **1,362** | **1.586e-7** |

Examples (flat pair index, emitted, reduced, ulp):

```
[0, (0,0)]      -712235.0775261772  -712235.0775260957   700
[1, (0,1)]      -712278.3635208694  -712278.3635207687   865
[2, (1,1)]      -712905.5676338089  -712905.5676337201   763
```

The 400-pair committed **replica gate is 400/400 bit-identical to the
emitted values** (including the MAP pair), so the baseline is not a
model of the arithmetic — it *is* the emitted arithmetic.

Decisions are untouched at this locus: the argmax pair is
`(369,423)` (flat 90,045) under both, its value moving by 7.5e-8; the
committed top1–top2 gap is 2.3778 = **2.04e10 ulp**, so the 1,362-ulp
perturbation has a 1.5e7× margin. The locus has no exact class ties
(`qual_tied_list_count: 0`), so nothing decision-relevant is within the
band.

### 6.2 The rest of the file, exhaustively (no sampling)

Every unordered class pair of every locus whose `pairs × units` term fits
the offline budget (17 loci, 206,048 pairs — loci 5, 8, 10, 14, 15, 16
were not censused exhaustively; 9 has no units):

| locus | F | U | N | const units | pairs compared | bit-mismatched | max ulp | max abs delta | argmax same |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|:-:|
| 0 | 93 | 634 | 261 | 23 | 4371 | 4364 | 95 | 8.640e-11 | True |
| 1 | 155 | 716 | 245 | 139 | 12090 | 11966 | 63 | 1.146e-10 | True |
| 2 | 124 | 35 | 1 | 35 | 7750 | **0** | 0 | 0.000e+00 | True |
| 3 | 143 | 527 | 202 | 32 | 10296 | 10167 | 60 | 5.457e-11 | True |
| 4 | 109 | 1386 | 241 | 754 | 5995 | 5822 | 127 | 4.620e-10 | True |
| 6 | 119 | 452 | 174 | 40 | 7140 | 7010 | 42 | 3.820e-11 | True |
| 7 | 113 | 1124 | 228 | 655 | 6441 | 6430 | 120 | 2.183e-10 | True |
| 11 | 99 | 621 | 162 | 191 | 4950 | 4922 | 69 | 6.276e-11 | True |
| 12 | 166 | 396 | 158 | 38 | 13861 | 13298 | 39 | 3.547e-11 | True |
| 13 | 154 | 487 | 201 | 14 | 11935 | 11789 | 54 | 4.911e-11 | True |
| 17 | 40 | 43 | 2 | 42 | 820 | **6** | 1 | 1.137e-13 | True |
| 18 | 124 | 401 | 202 | 37 | 7750 | 7481 | 33 | 3.001e-11 | True |
| 19 | 343 | 495 | 181 | 147 | 58996 | 58618 | 91 | 8.276e-11 | True |
| 20 | 241 | 826 | 348 | 157 | 29161 | 28860 | 68 | 1.237e-10 | True |
| 21 | 101 | 380 | 123 | 11 | 5151 | 4712 | 31 | 2.819e-11 | True |
| 22 | 134 | 964 | 211 | 512 | 9045 | 8970 | 89 | 1.619e-10 | True |
| 23 | 143 | 910 | 230 | 155 | 10296 | 10173 | 83 | 1.510e-10 | True |

**194,588 of 206,048 pairs (94.44%) change bits; the argmax pair is
unchanged at every censused locus; the worst deviation anywhere is
127 ulp.** The two loci that do *not* (uniformly) change are the ones
whose constant units form a **prefix** of the unit order:

* locus 2 — 35/35 units constant, i.e. `S` *is* the accumulation:
  **0 of 7,750 pairs differ** (bit-identical, by identity of the op
  sequence, not by luck);
* locus 17 — 42/43 constant but the single non-constant unit sits at
  index 11 with constants *after* it: 6 of 820 pairs differ (the
  trailing-constant adds are exactly what breaks it).

That is the precise statement of when reduction 1 stays bit-exact: **the
dropped addends must form a prefix** (drop them and seed the accumulator
with `S`, and the rest of the walk is unchanged). At T00048 the only such
locus is the degenerate 35-unit one.

### 6.3 One honest caveat on the baseline of §6.2

§6.2's baseline is my committed-order replica, not the emitted value
directly (the emitted `class_log_likelihoods` were compared to the
replica as a gate in the same pass): **187,742 of 206,048 pairs (91.1%)
reproduce the emitted bits exactly; the remaining 18,306 differ by
exactly ±1 ulp** and are confined to 8 loci (0, 1, 4, 7, 11, 20, 22, 23);
the other 9 censused loci replicate 100%. I could not attribute that
residual: no unit-order permutation, no chunked/pairwise/tree summation,
no single-term ±1-ulp shift, and no `fma` contraction of
`0.5*exp(a-m) + 0.5*exp(b-m)` (both products are exact, so contraction
cannot change it) reproduces it; `libm` is shared dynamically
(`ldd` shows `libm.so.6`, so Rust's `exp`/`ln` and Python's are the same
glibc), and `mix(v,v)=v` holds over every value in the file. It does not
move any verdict here — 1 ulp against deviations of 31–1,362 ulp — and at
the deep locus, where the headline is measured, the replica is exact
400/400. Pinning it down is a one-shot job with the existing
instrumentation: the receipt already carries per-`(fold,unit)`
`pair_term`s (`realign_score.rs:11237`), so a single dumped pair would
identify the differing term.

## 7. Why the file reading is trustworthy (gates run)

* `ll_matrix` extraction (byte-slice → flat f64, reshape `F×U`, counts
  from `units[].count`) re-derives the receipt's own counters **exactly**:
  locus 8 `unit_count 45423`, `class_pairs 92235`; `distinct_rows 406`
  matches `phase_census.scan_memo.distinct_rows` for the same locus;
  locus 2's single distinct row, locus 9's empty matrix.
* `E` recovered from the file's own constant class is **bit-identical**
  to `scoring.elsewhere_log_prob` at every locus (`M = 1` and its value
  equals `E` at all 24 loci).
* The committed loop replica matches the emitted `class_log_likelihoods`
  bit for bit (§6.1, 400/400; §6.3, 91.1% over 206,048 pairs with the
  ±1-ulp residual above).
* `mix(v,v) == v` over all 24,392 distinct ll values
  (`mix_self_check.py`, 0 violations) — the lemma the constant drop
  rests on.

## 8. What this buys, and what it must be labelled

* The unit dimension is by far the bigger lever at the heavy loci (54× vs
  1.12× for the row memo at locus 8; 104× vs 2.46× at locus 10). The two
  are orthogonal and compose to 41× over the whole file, 60× at locus 8.
* **Most of the deep-locus win already exists in the tree, in the
  non-exact mode only**: `ClassScanFastIndex` (`realign_score.rs:5887`)
  skips entries equal to `E` per `(fold,unit)` and adds back one
  `elsewhere * constant_mass` term per pair (`constant_mass` = the pair's
  non-union count mass, ≥ the 43,059 mass of the all-`E` columns). Its
  inner walk is ~1,892 terms at this locus (§5) = 24.0×, and it is already
  documented as an accumulation-order change ("the fast mode's sums are
  measured and named ... never hidden"). What is genuinely new is the
  *identical-vector merge*: 1,892 → 842 terms, a further **2.25×**; and
  53.95× against the committed dense walk.
* Because it changes bits, it belongs behind a named toggle with the model
  string suffixed, exactly like `IMPG_CLASS_SCAN_ORDER=sparse-skip`
  (`ClassScanOrder` at :5818) — never in the default/exact path. The
  existing gates would reject it by construction:
  `memoized_scan_is_bit_identical_to_direct_on_randomized_rows` runs
  **one-ulp near-twin negative controls** and demands bit-for-bit
  equality per pair (`docs/syng-gmem-bwt/memoization-receipts.txt:60-75`).
  That is the correct behaviour: the row memo is bit-exact (it *shares*
  values, never reorders); the unit reductions are not.
* If a bit-identical version of reduction 1 is ever wanted, the prefix
  case (§6.2, locus 2) is the only free one, and it is not the shape the
  data takes. A compensated/double-double accumulation would not help
  either: it changes the emitted bits too.

## 9. Residual risks / open items

1. **Not byte-identical** (the whole point). Any use of these reductions
   in a shipped path is a labelled accumulation-order change; the
   fixture-scale identity gates must be re-baselined, not waived.
2. **Near-ties.** Deviation is 1–1,362 ulp. At locus 8 the decision margin
   is 2.04e10 ulp (safe by 1.5e7×) and there are no exact ties. The class
   of case that *would* break is a fixture whose top-1/top-2 (or tie
   certificate) sits inside a few thousand ulp — the repo already
   constructs one-ulp near-twins deliberately.
3. **The ±1-ulp replica residual at 8 of 17 censused loci** (§6.3) is
   unexplained. It weakens nothing measured here, but it means my replica
   is not a certified model of the scan's accumulation everywhere; the
   deep-locus headline does not depend on it (400/400 against the emitted
   values).
4. **The bounded scan's bounds are a separate account.** `pair_bounds`
   uses per-row aggregates (`E*C + A_f + A_g`, with `A` built from
   `units[].count` and per-row spread, :6100-6200). The unit reduction
   leaves those intact, but a real implementation must check that the
   merged counts feed the bound identically (the mass per unit changes
   identity, not the row totals) — not measured here.
5. Wall-clock: no claim is made. The only timing in evidence is the
   receipt's own `phase_census.classes_seconds = 3.802` for locus 8 from
   the original pipeline run (not re-measured by me), against 4.19e9
   nominal terms.

## 10. Reproduction

Scripts (committed, pure python + numpy, no pipeline run, no network, one
input file plus its sibling receipt):

```
scripts/unit-redundancy/analyze.py         per-locus table       (~18 s)
scripts/unit-redundancy/exactness.py       deep-locus arms       (~110 s)
scripts/unit-redundancy/census.py          exhaustive 17 loci    (~28 s, 24 procs)
scripts/unit-redundancy/mix_self_check.py  the mix(v,v)=v lemma  (~25 s)
```

```
ING=/tmp/probe-fast/T00048/instrument-receipt.jsonl.ingredients.jsonl
REC=/tmp/probe-fast/T00048/instrument-receipt.jsonl
OUT=docs/syng-gmem-bwt/unit-redundancy
python3 scripts/unit-redundancy/analyze.py       $ING $OUT/loci.json
python3 scripts/unit-redundancy/exactness.py     $ING $REC $OUT/exactness.json
python3 scripts/unit-redundancy/census.py 60000000 $ING $REC $OUT/census.json
python3 scripts/unit-redundancy/union_probe.py   8 3000   # the 1,892-term check
python3 scripts/unit-redundancy/mix_self_check.py $ING
```

The JSON receipts under `docs/syng-gmem-bwt/unit-redundancy/` are the
outputs of exactly those commands.

```
scripts/unit-redundancy/union_probe.py     sparse-skip inner length (~5 s)
```

Walls above are wall-clock on `erikg`'s box with **other work running**
(load average 15–18 on 256 cores; a `scatter_sort` and another `impg` were
resident during the runs), single-threaded except `census.py`'s 24-process
pool (2m29 user / 28 s wall). The measurements themselves are exact
counters (bit comparisons, integer counts), not timings.
