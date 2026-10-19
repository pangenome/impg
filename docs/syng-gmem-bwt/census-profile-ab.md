# THE CENSUS-PROFILE A/B — the measured decision input (Ruling 2, UNDECIDED)

**THE PLAIN-LANGUAGE ANSWER FIRST.** The two profiles are: **(A) the
committed-slice read set** — the cut restricted to the committed 1x
kept-reads (the evidence mass the original local-frame receipts
banked), and **(B) the global-native full mass** — every routed read
the global scan keeps (1.002x the pilot's mass, but 2.7x/3.2x/8.4x/
15.2x the committed slice's census mass at L1/L2/L3/L4). **On
accuracy the two are a wash at every measured locality**: the
truth-pair rank-1 rates are IDENTICAL at L1, L3, L4 and the pilot,
and B buys exactly ONE truth locus at L2 (rank-1 24% → 28%, locus
17: B calls the truth pair and A does not, official QV 50 vs 42);
the median official QV is within 1.5 of the other profile
everywhere (A is the higher median at L1 and L4, equal at L2, L3,
pilot), and the ≥QV30 counts differ by at most ONE locus in either
direction. Where the two profiles call different fold material
(7-9 loci per L-locality, 0 at the pilot), most differences are
sequence-equivalent copies (equal QV); of the QV-visible flips, A
is closer to truth at 6 loci and B at 3. **The cost is NOT a wash**:
profile B's loci core runs **1.77x (L1), 2.49x (L2), 5.02x (L3),
7.68x (L4)** profile A's, and **1.008x at the pilot** (the two
profiles coincide there — the pilot's committed slice already
carries essentially the whole routed mass); at the whole-call
ex-load level (loci + bind + rows + the small classes) the
multiplier is **1.57x / 1.75x / 3.07x / 3.67x / 1.008x**. The
§10.6.14 autopsy's "~1.6x call-core multiplier at L1" is confirmed
at L1 (1.57x whole-call, 1.77x loci-phase) and is *worse* further
from pattern saturation — L3/L4 pay 3-4x whole-call because their
full-mass censuses are 8.4x/15.2x the committed slice. **This
document does not decide** — the owner does; every number below is
a receipt (the runs: `rulings-ab-{A,B}-<L>`, the analysis:
`rulings-ab-analysis.txt`, the cuts: `global-cut-<L>-shadow1x`).

## The machinery (the same assessment as every banked accuracy receipt)

Both profiles of every locality ran the committed global-frame call
recipe (serial, exact class scan, the executed-Ruling-1 production
path) with the OFFICIAL QV machinery: `--truth-qv <axis> <truth2>`,
`IMPG_TRUTH_QV_ALLWAVE=1`, the likegt max-qv convention (per-truth-
hap best alignment, identity = matches/query_span, QV =
−10·log10(1−identity) capped at 50), the diplotype QV the mean of
the two per-hap QVs. Rank-1 is the truth pair's class-LL rank (the
committed `truth_rank` convention). Ten runs, exit 0 all, co-tenant
×10 (the per-phase receipts are the honest numbers; the walls
301-1006s).

## THE TABLE (per locality, per profile)

| locality | profile | census mult | loci phase | loci ratio B/A | rank-1 rate | QV median | ≥QV30 | at-cap (QV50) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| L1 | A (committed slice) | 141,114 | 75.0s | — | 4/25 (16%) | 42.6 | 24/25 | 6/25 |
| L1 | B (global full mass) | 384,056 | 132.5s | **1.77x** | 4/25 (16%) | 41.1 | 22/25 | 6/25 |
| L2 | A | 48,641 | 16.9s | — | 6/25 (24%) | 43.5 | 24/25 | 6/25 |
| L2 | B | 156,976 | 42.0s | **2.49x** | **7/25 (28%)** | 43.5 | 23/25 | 7/25 |
| L3 | A | 67,359 | 24.2s | — | 2/24 (8%) | 40.0 | 24/24 | 3/24 |
| L3 | B | 567,243 | 121.4s | **5.02x** | 2/24 (8%) | 40.0 | 23/24 | 3/24 |
| L4 | A | 18,121 | 13.6s | — | 3/25 (12%) | 45.0 | 25/25 | 5/25 |
| L4 | B | 275,405 | 104.5s | **7.68x** | 3/25 (12%) | 43.5 | 25/25 | 4/25 |
| pilot | A | 512,929 | 630.9s | — | 4/22 (18%) | 36.3 | 17/22 | 4/22 |
| pilot | B | 513,836 | 635.8s | **1.008x** | 4/22 (18%) | 36.3 | 17/22 | 4/22 |

(The whole-call ex-load ratios — loci + bind + rows + the small
classes, the mass-invariant IO classes excluded per the §10.6.14
convention: L1 1.57x, L2 1.75x, L3 3.07x, L4 3.67x, pilot 1.008x.
The pilot's two profiles are the SAME mass to within 0.2% — the
committed pilot slice already routes essentially every kept read,
so the pilot row is the A/B's own control: identical answers,
identical QVs, no differing locus, 1.008x.)

## THE NAMED LOCI WHERE THE PROFILES DIFFER

Fold-material differences (the called diplotype differs; the
per-profile verdict is that profile's own call vs its own truth
pair, and the official QV of its called material):

**L1 (8 differing loci; 4 QV-visible):**
- locus 4: A=[56,86] B=[12,59] — neither matches truth; A higher QV
  (40.0 vs 38.5)
- locus 11: A=[6,144] B=[48,121] — neither; A higher (30.4 vs 25.9)
- locus 13: A=[0,13] B=[108,108] — neither; **A higher (45.0 vs
  26.9 — the widest QV swing of the table)**
- locus 18: A=[2,94] B=[5,98] — neither; A higher (45.0 vs 40.0)
- loci 2, 5, 12, 22: differing fold material, EQUAL QV
  (sequence-equivalent copies)

**L2 (7 differing; 3 QV-visible):**
- locus 17: A=[0,0] B=[0,6] — **B RIGHT, A WRONG (B calls the truth
  pair: QV 50 vs A's 42; the one clean B win of the whole table)**
- locus 5: A=[9,133] B=[9,39] — neither matches truth; B higher
  (B's call scores 45.0, A's 43.5)
- locus 6: A=[83,101] B=[83,83] — neither; A higher (34.8 vs 28.2)
- loci 3, 16, 19, 23: equal QV

**L3 (9 differing; 3 QV-visible):**
- locus 0: A=[12,150] B=[12,206] — neither; B higher (43.5 vs 39.0)
- locus 1: A=[1,16] B=[16,42] — neither; A higher (43.5 vs 34.3)
- locus 17: A=[7,7] B=[33,49] — neither; **A higher (40.0 vs 24.4)**
- loci 8, 9, 13, 16, 19, 23: equal QV

**L4 (7 differing; 2 QV-visible):**
- locus 6: A=[0,64] B=[21,64] — neither; A higher (45.0 vs 37.0)
- locus 23: A=[0,13] B=[113,120] — neither; **A higher (50.0 vs
  38.5)**
- loci 0, 10, 14, 15, 19: equal QV

**pilot: NO differing loci (the control row).**

Score of the QV-visible flips: **A closer to truth at 9 loci**
(L1: 4, 11, 13, 18; L2: 6; L3: 1, 17; L4: 6, 23), **B closer at 3**
(L2: 17 — the one clean B win, B's call IS the truth pair — and
L2: 5; L3: 0); every other fold-material difference is a
sequence-equivalent copy (the equal-QV class — the fold material
flips but the called sequence is the same class either way).

## What this does and does not say (the owner's decision, untouched)

- The superset mass (B) does NOT buy rank-1 accuracy at 4 of 5
  localities (identical rank-1 rates at L1/L3/L4/pilot) and buys
  exactly one truth locus at L2 — while the committed slice (A) is
  closer to truth at more of the QV-visible flips than B is.
- The superset's real measured price is the §10.6.14 multiplier,
  confirmed at L1 and WORSE away from saturation: 1.57x-3.67x the
  whole call ex-load (1.77x-7.68x the loci core) at the L-localities.
- The pilot is the degenerate row: the two profiles coincide (0.2%
  mass apart) — identical answers, 1.008x.
- NOT measured here: the cut side's own cost difference (the
  restricted scan saves ~8% of the cut wall per §10.5.4 — the scan
  decodes the same shards either way), the marathon-class loci (the
  reaped substrates), and the fleet-level memory/IO effects (both
  profiles' calls held the same working set class here).

**The decision is the owner's. This table is the input.**
