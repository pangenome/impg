# THE GENOME-WIDE STRATIFIED SAMPLE — tranche 1 (the owner's methodology ruling, executed)

**THE OWNER'S RULING (the motivation, verbatim):** *"this is dangerous
to be drawing lots of broad conclusions off of six random loci."* The
five-locality evidence (`census-profile-ab.md` — the wash; the one-locus
L2-17 cost; the 1.77-7.68x cost ratios) and the §10.6.17 stacked-rulings
form (the 1.38-30.4 per-locus budget spread) were measured at FIVE
localities. This document converts those conclusions into SAMPLED
GENOME-WIDE evidence: (a) the census-profile wash confirmed or refuted
at scale — the conserved-elsewhere dictionary class measured as a
RATE, not an anecdote; (b) the per-locus CPU-work variance across
strata — the genome wall as a BOUNDED estimate, not an L1 point value;
(c) the far-bin prune's value across strata.

**THE TRANCHE'S OWN SCOPE, STATED FIRST (the honest limits, no genome
claim the tranche cannot support):**
- The substrate is the RUNTIME VERIFICATION's export-profile class
  (the c-rider locality-maps-export, K=63, syng-native) — the runtime
  receipt, NOT the pggb-class accuracy substrate of the five-locality
  A/B. The wash verdict below is measured at the EXPORT-PROFILE maps;
  the five-locality A/B was measured at the committed pggb-era maps.
  Both are the same call machinery and the same census profiles; the
  maps class is the one substrate axis the tranche does not vary.
- The marathon roster (16 members), the poison census (5), and the 16
  cut.exit terminals are OUT of the draw (named below); the wall
  estimate carries the marathon tail as the committed bound, never as
  a tranche measurement.
- Every wall below is a CONTENDED wall (the co-tenant state is
  receipted in the work dir's `logs/co-tenant.log`; another worker's
  gate fleet has been on the box since before the tranche opened).
  The CPU-work receipts are the SERIAL (RAYON=1) phase timers — the
  autopsy's own method — which are core-seconds, not walls.

## 1. THE STRATA (from the committed ledgers only — no new run feeds the selection)

| ledger | what it names |
|---|---|
| the split ledger | 5,844 ACCEPTED territories (the product's universe); 9,891 windows; 3,252 transitive-expansion rejections (the expansion gate); 795 ambiguous-axis-bridge |
| the runtime receipts | ~2,078 committed cut-stats (reads_kept / census_mult); the per-territory call exit codes; split.info.json (per-window axis/truth-row counts) |
| the grid's segdup record | territories.bed column 5 (the per-territory segdup fraction; the committed rollup's own 0.02 class boundary) |
| the marathon roster | 16 named deep-lane memory-tail members (datasheet-runtime.md) |
| the poison census | the 5-member syng segdup-tail prep-terminal class |
| the cut.exit ledger | 16 terminals (the chr2 AGC-offset class 12 + the rest) |

**THE RUNNABLE UNIVERSE: 5,807** = 5,844 − 16 marathon − 5 poison − 16
cut.exit. **THE STRATA (disjoint, priority-ordered):**

| stratum | the rule (committed-ledger terms only) | size in U | drawn |
|---|---|---:|---:|
| S_EXP EXPANSION-GATE-ADJACENT | a bed4 genomic neighbor (±2 windows, same chrom) is a transitive-expansion rejection | **4,413** | 8 |
| S_SEG SEGID-ADJACENT | the grid's own segdup class (fraction > 0.02) OR a marathon member within ±2 windows | 217 | 10 |
| S_REP REPEAT-HEAVY | the top decile of census_mult/reads_kept among the committed cut receipts (boundary computed: ≥ 1.629) | 34 | 8 |
| S_BIND BIND-AMBIGUITY | the ledger-named 7479-class: the RV campaign's exit-1 call failures | 79 | 8 |
| S_ORD ORDINARY | everything else | 1,064 | 22 |

**THE LEDGER'S OWN SURPRISE, stated:** expansion-gate adjacency is the
NORM, not the exception — 76% of the runnable universe sits within two
windows of an expansion-gate rejection (52% immediately adjacent). The
"ordinary" bulk is the minority class (18%). The allocation therefore
OVER-samples the rare strata deliberately (min 8 per stratum) and the
genome-wide estimates below WEIGHT BACK by stratum size — the honest
stratified estimator with bounds, never a raw tranche mean.

## 2. THE DRAW (the seed and the rules, reproducible)

Seed **20261017**, Python `random.Random(SEED)`, each stratum's pool
sorted by TID, `random.sample(pool, n)` — a pure function of the
committed ledgers + the seed (the draw script is
`scripts/genome-tranche/genome-tranche-draw.py`; the tranche manifest is
`genome-tranche1/tranche.tsv` in the work dir). 56 territories drawn.
One in-flight failure class found by the tranche itself (named in §3).

## 3. THE TRUTH PAIR AT GENOME SCALE — the truth2 bridge

The truth pair is (the HG002#1 axis row — `prep/<T>/axis-path.txt`,
ALREADY the committed axis-bridge output for the paternal haplotype —
and the HG002#2 row of the refined bed). **1,270 of 5,844 accepted
territories carry MORE THAN ONE HG002#2 refined-bed row** (the
cross-chromosomal segdup-hit class the split's own axis-row rule
receipts); the five-locality A/B's first-row convention would name a
WRONG truth at exactly the segdup territories the tranche most needs.
The fix is the bridge's own machinery applied to the maternal haplotype
(`tranche-truth2-bridge.py`): each territory's window (the CHM13#0
bridge row) aligned asm5 to HG002#2#chr1-22 (extracted from the AGC,
the `hg002-axis-target.fa` recipe verbatim, mappy 2.31), the truth2 row
= the refined bed's max-overlap same-contig row — the split's own
axis-row rule. **All 56 tranche members resolved "ok"** (the
no-bridge / ambiguous / no-row fail-closed classes: zero at this
sample size).

## 4. THE MACHINERY (the fleet pipeline's first genome-scale use)

Per territory: the prep MIRROR+REGEN (the campaign's own chain, run in
the tranche's own work dir — the campaign's ledgers never touched);
the reads fetch (the committed ±50kb rule); the RV local-derive cut
RE-RUN (the restriction basis — the committed kept-read set — with a
free determinism gate against the campaign's banked cut-stats receipts);
the profile-A restricted cut (the single `cut --frame global`, the only
committed restricted path); the profile-B global-native cuts + the
production calls via `fork-cut-call-batch-global` (the shared panel
load, width-32 children, the 4-concurrent T16 in-child calls); and the
official truth assessment (`call-diplotypes --truth-qv`, the likegt
max-qv cap-50 convention, RAYON=1 serial — the same assessment as every
banked receipt), both profiles.

(RESULTS TABLES FOLLOW — filled by `tranche-analyze.py` at the stage
close; see §5-§7.)

## 5. THE WASH VERDICT AT SCALE

(to be filled)

## 6. THE CPU-WORK VARIANCE AND THE BOUNDED GENOME WALL

(to be filled)

## 7. THE FAR-BIN PRUNE'S VALUE ACROSS STRATA

(to be filled)

## 8. WHAT TRANCHE 2 WOULD ADD

(to be filled)
