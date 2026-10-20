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

## 4. THE MACHINERY (the fleet pipeline's first genome-scale use) — THREE census/profile forms

**THE OWNER'S MID-TRANCHE STEERING (executed): the far-bin prune is the
product default** (commit 2a3dd46, `--far-bin-prune` — the pure
window-geometry rule; the five-locality gates: zero truth losses at
L1-L4, the pilot WINS FOUR truth loci, L1 loci −53.6%). The tranche
therefore runs THREE forms per territory, not two:

| form | the census | the call | the role |
|---|---|---|---|
| **A** | restricted (the committed kept-read set) | the committed path | the OLD default |
| **B** | global-native (the full routed mass) | the committed path | the record form |
| **C** | restricted (form A's own cut, reused — the prune is a CALL-side census-view flag) | `--far-bin-prune` | **THE NEW PRODUCT DEFAULT** |

All three assessed with the official QV machinery. The tranche verdicts
now answer the owner's actual product questions: **does the prune's
four-truth-win hold at sampled genome scale** (C-vs-A: wins, and the
gate bar of zero truth losses), and **does the wash conclusion hold
under the new default** (B-vs-C beside B-vs-A). The far-bin mass-share
analysis (§7) IS the prune's genome-scale value receipt, tied to form
C: the machinery's own pruned-records/mass log lines at every
territory, cross-checked against the census-geometry share.

Per territory: the prep MIRROR+REGEN (the campaign's own chain, run in
the tranche's own work dir — the campaign's ledgers never touched);
the reads fetch (the committed ±50kb rule); the RV local-derive cut
RE-RUN (the restriction basis — the committed kept-read set); the
form-A restricted cut (the single `cut --frame global`, the only
committed restricted path); the form-B global-native cuts + the
production calls via `fork-cut-call-batch-global` (the shared panel
load, width-32 children, the 4-concurrent T16 in-child calls); and the
official truth assessment (`call-diplotypes --truth-qv`, the likegt
max-qv cap-50 convention, RAYON=1 serial — the same assessment as every
banked receipt), all three forms.

### 4.1 THE T00041 DETERMINISM-GATE INSPECTION (the era class, receipted)

The stage-K free determinism gate (the regenerated RV cut vs the
campaign's banked cut-stats) flagged T00041. **NOT a nondeterminism —
a deterministic, DATED behavior boundary in the campaign's own
receipts:** every tranche member whose campaign substrate was built
Oct 16-17 reproduces EXACTLY (13/13 — partition region counts equal,
cut-stats identical); every member built in the EARLY Oct-15 prep
waves differs (17/17 — a systematic ~2x region-count change: T00041
54→25, T00531 340→165, T01091 327→164, T01579 1189→549, same locus.fa
bp totals). The boundary: the campaign's binary changed mid-day Oct 15
with the committed Rust source CONSTANT across it (no src commit
between 231d21a and c78c450) — the **vendor/syng submodule pointer
drifted 054acfc→68ac197 during Oct 15** (recorded uncommitted-by-design
in the campaign's own handoff; the current checkout is 68ac197), and
the partition step's region semantics changed with it. The campaign's
committed cut receipts are therefore a MIX of two substrate eras; the
mirror+regen reproduces the LATER era byte-exactly and deterministically
cannot reproduce the early era. **The tranche's own receipts are
internally consistent — ONE binary produced every cut and every call of
all three forms — which is what the A/B/C comparisons require.** The
class is named (the OCT-15-EARLY-SUBSTRATE ERA class); no redraw, no
bug (full inspection: `rvcut/determinism-gate-inspection.txt` in the
work dir).

## 5. THE WASH VERDICT AT SCALE — PARTIAL (the checkpoint receipt; n stated per stratum; S_EXP and S_BIND have ZERO receipts at this checkpoint — their cuts/calls are still landing)

**THE MID-TRANCHE BINARY BOUNDARY, receipted first (it forced a receipt
sweep):** the product-defaults session (the co-tenant worker on this
box) landed the owner's defaults-flip mid-tranche — commits 82583eb,
e03f7d5, and 431d55a (§10.6.19: the far-bin prune DEFAULT-ON in the
call path) — and rebuilt the binary twice (sample_projection 21:09Z,
impg 22:10Z). The first A/B truth receipts after 22:10Z were therefore
silently PRUNED (the model label caught it: `...-farbin-prune`); they
were CLEARED and re-made with explicit flags — **on the flipped
binary every form states its own census view: A/B =
`--far-bin-prune=false` (the pre-ruling committed unpruned view), C =
the prune ON** — and the stage scripts now pin the forms explicitly
(run-truthcalls.sh / run-truth-bounded.sh). The one cross-binary delta
receipted (the fleet's in-child call vs a later CLI call at T00531: 2
of 25 loci differ in record_count/unit_count + LL/qual, ZERO
fold-material flips) is the binary-boundary class, not a
nondeterminism.

**The partial table (serial RAYON=1, official QV machinery, cap 50;
every rate carries its denominator; the assessment substrate is the
export-profile maps, per the scope header):**

| stratum | form | territories | loci | truth-expressible | rank-1 | QV median | ≥QV30 | at-cap |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| S_ORD | A restricted | 2 | 49 | 49 | 18/49 (37%) | 45.0 | 47/49 | 19/49 |
| S_ORD | B global-native | 3 | 46 | 46 | 10/46 (22%) | 45.0 | 44/46 | 14/46 |
| S_ORD | C restricted+prune | 1 | 24 | 24 | 8/24 (33%) | 43.1 | 24/24 | 9/24 |
| S_SEG | A | 2 | 43 | 43 | 14/43 (33%) | 43.1 | 43/43 | ~9 |
| S_SEG | B | 3 | 46 | 46 | 10/46 (22%) | 45.0 | 44/46 | 14/46 |
| S_SEG | C | 1 | 24 | 24 | 8/24 (33%) | 43.1 | 24/24 | 9/24 |
| S_REP | A | 4 | 49 | 49 | 12/49 (24%) | 43.5 | 46/49 | 12/49 |
| S_REP | B | 4 | 49 | 49 | 12/49 (24%) | 45.0 | 45/49 | 12/49 |
| S_REP | C | 5 | 62 | 62 | 14/62 (23%) | 43.1 | 58/62 | 14/62 |
| S_EXP | all | 0 | — | — | — | — | — | — |
| S_BIND | all | 0 | — | — | — | — | — | — |

**THE PARTIAL WASH READS AS A WASH** (the five-locality conclusion
holding so far at sampled scale): the rank-1 rates sit within 1-11
percentage points of each other across forms in every landed stratum,
with no form systematically ahead (n is small; per-territory rates
vary more than per-form rates do).

**THE CONSERVED-ELSEWHERE RATE: 0 losses in 141 loci assessed**
(B-right-A-wrong 0; B-right-C-wrong 0). The L2-17 class was 1-in-~121
at the five localities; a 0-in-141 tranche observation is consistent
with a rare class and bounds nothing yet — the full tranche (~1,300
loci) and tranche 2 are the denominators that will measure the RATE.
The QV-visible differing loci so far are 7, mixed in direction
(e.g. T01069 locus 21: A qv 35.8 vs B 0.2 — B's call is terrible
there; T01651 locus 22: B 45.0 vs A 39.8).

**THE CENSUS-MASS RATIO B/A AT GENOME SCALE: median 1.05-1.21x**
(S_ORD 1.20x, S_SEG 1.21x, S_REP 1.05x; n=28 territories with both
cuts) — NOT the five-locality 2.7-15x. The reason is structural: the
genome-scale restriction basis (the RV campaign's kept-read sets at
35x) already carries nearly the whole routed mass, where the
five-locality committed slices were 1x-subsampled. **The wash's COST
side at genome scale is a ~1.2x census mass, not 3-15x.**

## 5b. THE PRUNE'S OWN VERDICT AT SCALE (the owner's product question)

**ZERO PRUNE TRUTH WINS, ZERO LOSSES in 98 paired C-vs-A loci — the
gate bar (zero truth losses) HOLDS at this n.** 4 QV-visible changes
(3 higher under C, 1 higher under A). The five-locality four-truth-win
(pilot) has no analog yet in the tranche's landed subset; the
high-share members (where wins live) are still being called. The
prune's wall effect at the landed low-share pairs: n=2, median C/A
1.11x — no measurable drop yet; the high-share pairs are the
outstanding receipt.

## 6. THE CPU-WORK VARIANCE AND THE BOUNDED GENOME WALL — PARTIAL (n per row; the S_EXP/S_BIND strata have NO receipts at this checkpoint)

**Per-territory ex-load call CPU-s (serial RAYON=1 truth calls; the
autopsy's own method; the panel-load inputs phase excluded):**

| stratum | n | median | mean | p90 | max |
|---|---:|---:|---:|---:|---:|
| S_ORD (form A) | 2 | 151.6 | 151.6 | 254.7 | 254.7 |
| S_SEG (form A) | 1 | 64.9 | 64.9 | 64.9 | 64.9 |
| S_REP (form A) | 3 | 379.4 | 410.9 | 675.2 | 675.2 |
| S_ORD (form C) | 1 | 59.8 | 59.8 | 59.8 | 59.8 |
| S_SEG (form C) | 1 | 64.0 | 64.0 | 64.0 | 64.0 |

**Per-locus core-s (the loci phase / loci): S_REP median 14.55 vs S_ORD
3.04 vs S_SEG 1.46 — a ~10x cross-stratum spread** in the call core,
exactly the variance class the L1-class point extrapolation missed.
The per-locus wall distribution at form C: S_ORD median 0.60s / p90
2.60s / max 36.1s / worst-5% share 62.1%; S_SEG median 0.90s / p90
3.10s / max 9.9s — the worst-5% concentration reproduces the
five-locality finding at ordinary genome territories.

**THE GENOME WALL AT THIS CHECKPOINT CANNOT BE BOUNDED YET — STATED
PLAINLY:** the stratified estimate needs every stratum's receipts,
and S_EXP (4,413 of 5,807 territories, 76% of the universe) plus
S_BIND (79) have zero receipts so far. The interim arithmetic from the
landed strata (S_ORD 152s × 1,064 + S_SEG 65s × 217 + S_REP 411s × 34
≈ 0.19e6 core-s, CI ±0.22e6 at n=6) is a LOWER-BOUND FRAGMENT, not an
estimate; if S_EXP behaves like S_ORD (the likeliest prior —
expansion-gate adjacency is the norm, and the landed S_ORD members
are themselves mostly expansion-adjacent by the ledger's own
geometry), the full-universe ex-load call total projects to
~0.6-0.8e6 core-s ≈ **~0.7-0.9h at 256 cores** — a PROJECTION, labeled
as such, to be replaced by the receipted stratified total at the
tranche close. The marathon roster (16 members, the RV-frame p90
2,912s-deepest 11,084s serial lanes), the 5 poison, and the 16
cut.exit terminals are EXCLUDED from the universe and carried as the
committed bound beside the estimate.

**THE CUT SIDE (the contended walls, co-tenant state receipted in
logs/co-tenant.log):** cutA (the single restricted cut, per-territory
panel load paid) median 724-926s by stratum; the fleet form (the fork
driver, shared panel load, width-32 children, T16) cut all 55 in ONE
wave: chunk-0's 32-row cut phase 3,908s + chunk-1's 23 rows, whole
driver wall 10,718s incl. the production calls, peak tree RSS 264GB —
the width-32 amortization holds at genome scale under ~2x co-tenancy
(the W32 receipted per-row cut share was 34.9s solo; the tranche's
contended class is receipted, not comparable wall-for-wall). The
production calls: 5 of 8 landed (T06027's failed on the cut-collision
error below; T01091/T02854's chunk-1 slots did not complete before the
child exited — named, re-runnable).

**THE RV-CUT RESTRICTION-BASIS STAGE (stage K), receipted:** 55/55
complete. The adaptive width worked exactly as the campaign designed:
4 W=8 batches guard-killed at 69.4-70.2GB (the census working set is
not predictable pre-batch), every retry at width 4 exited 0. Two stage
bugs found and fixed mid-run, both receipted: the empty-name wipe
class was NOT reproduced (the campaign's fix is in the copied chain),
and a width-class table accumulation bug in the tranche's own batch
composer (stale members re-entering slices — fixed with a per-round
reset; the mis-composed batches were guard-killed and re-cut cleanly,
zero strikes).

## 7. THE FAR-BIN PRUNE'S VALUE ACROSS STRATA — the mass-share receipt (form C's own basis)

**THE TRANSCHE AGGREGATE: 23.85% of the restricted census mass is
far-bin scattered** (10,946,991 of 45,904,103 mult over the 28
territories with cutA receipts at the checkpoint) — the mass the new
default removes. By stratum (the census-geometry share, cross-checked
AGAINST FORM C'S OWN MACHINERY RECEIPTS — the pruned records/mass log
lines agree with the census geometry in every landed case):

| stratum | n | share median | p90 | max | zero-share |
|---|---:|---:|---:|---:|---:|
| S_ORD | 17 | 8.58% | 83.61% | 86.18% | 3/17 |
| S_SEG | 5 | 3.21% | 23.00% | 23.00% | 2/5 |
| S_REP | 6 | 18.75% | 80.97% | 80.97% | 0/6 |

The five-locality shares (L1 31.1% / L3 21.7% / L4 0.02% / L2 0.0%)
are NOT outliers: at genome scale the far-bin class is COMMON and
HEAVY-TAILED — ordinary territories carry up to 86% scattered mass,
and the median is ~9%. **The prune's genome-scale value is therefore
not a niche lever but a broad census-mass cut whose wall effect
should be measured at the high-share members — those pairs are the
tranche's outstanding receipt (the landed C-vs-A wall pairs so far
are low-share members: n=2, median C/A 1.11x, no drop visible — the
84%-share members' calls are still in flight).**

## 8. WHAT TRANCHE 2 WOULD ADD + THE EXACT RESUME POINT

**TRANCHE 2 (the ~150-territory extension) would add, in value order:**
1. **S_EXP and S_BIND at depth** — the checkpoint's biggest hole: S_EXP
   is 76% of the runnable universe with zero receipts; its 8 drawn
   members are mid-flight, but a rate estimate for the universe's
   dominant stratum wants 20-30 members. S_BIND (the 7479-class, 79
   members) measures whether the GLOBAL-FRAME call mode resolves or
   reproduces the RV campaign's bind-ambiguity fail-closed rate —
   zero receipts so far.
2. **The conserved-elsewhere RATE at decision grade**: 0-in-141 so far
   bounds nothing; the L2-17-class rate needs the full ~1,300-loci
   tranche + tranche 2's ~3,000 loci to land a per-stratum rate with
   a useful CI.
3. **The marathon members themselves** (excluded from tranche 1 by
   design): the global-frame call's behavior on the 16-member memory
   tail is the wall estimate's one unsampled term.
4. **The high-share far-bin pairs**: the prune's WALL value at the
   80%+ -share members (the five-locality L1 measured -53.6%; the
   genome tranche's landed pairs are low-share so far).

**THE EXACT RESUME POINT (everything marker-gated; the master and the
bounded passes keep running past this checkpoint):**
- `run-master.sh` (pid 2080501) is mid-stage-A (the form-A restricted
  cuts, 35/55 done at the checkpoint); next it runs TA (all-form-A
  truth calls) then TB — both over the FIXED run-truthcalls.sh (every
  form pins its census view explicitly after the defaults-flip).
- The stage-TC driver (form C, the new default) is SIGSTOPped at pid
  3422015 — `kill -CONT 3422015` resumes it (or kill it and relaunch
  `bash run-stageTC.sh`; it is marker-gated and idempotent).
- The bounded A/B passes (`run-truth-bounded.sh`, fixed) are draining
  their 31-member lists; their callX.done markers make the master's
  TA/TB skip exactly what landed.
- When all three forms' markers are complete: re-run
  `python3 scripts/genome-tranche/tranche-analyze.py <WT> <out>` and
  replace this checkpoint's tables with the full-tranche tables.
- The work dir: `/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1/`
  (tranche.tsv, truth2.tsv, tranche-failed.tsv, cutA/cutB/callA/callB/
  callC/, rvcut/ incl. determinism-gate-inspection.txt, fleetB/,
  logs/ incl. co-tenant.log, analysis-partial.txt).
- NOT touched by the tranche, by design: the campaign's ledgers
  under wgs35x-fleet/runtime-verification/, the product-defaults
  worker's files, vendor/syng.

## 9. TRANCHE 2 — THE DRAW, THE SESSION'S MACHINERY FINDINGS, AND THE GROWING RECEIPTS (the 2026-10-20 increment; every n stated; the stages still running by markers)

### 9.1 THE TRANCHE-2 DRAW (the seed and the rules, reproducible)

Seed **20261020**, the tranche-1 discipline extended (the draw script
`scripts/genome-tranche/genome-tranche2-draw.py`; the strata machinery
tranche-1 VERBATIM — the same committed ledgers, the same disjoint
priority-ordered strata, the same runnable-universe exclusions):

| block | the rule | n |
|---|---|---:|
| **S_EXP** | a fresh `random.sample` of 22 from the S_EXP pool (4,413 in U) MINUS tranche-1's 8 drawn members, pool sorted by TID, seed 20261020 | 22 |
| **S_BIND** | **THE CLASS CENSUS, not a sample**: every remaining member of the ledger-named 7479-class (79 in U; tranche 1 owns 8, mid-flight there) | 71 |
| **S_MAR** | the 16-member marathon roster census (outside the runnable universe; NEVER part of a stratified estimate — tranche 2 runs them to receipt the walls or the honest guard-kill class) | 16 |

109 rows, **row order = stage priority** (S_EXP first — the universe's
76% stratum with zero tranche-1 receipts at the checkpoint; then the
S_BIND census; then the marathon tail). One in-flight deterministic
failure found and named immediately: **T00762 extract-length-mismatch**
(refined-bed row end 1bp past the AGC contig — the SAME class as
tranche-1's T09025; `tranche-failed.tsv`, all stages skip it).
**The truth2 bridge resolved all 109 members "ok"** on first pass
(the bridge script's `#tid` header-key bug found and fixed — the
first master pass died in 10 seconds with every prep failed, the
missing `rv-extract.py` copy the second cause; both fixed, committed).

### 9.2 THE SESSION'S MACHINERY FINDINGS (all receipted in the work dirs' logs)

- **THE EMPTY-BASIS CLASS (8 tranche-1 members, repaired):** the
  pre-fix rvcut batch-composition bug re-composed batches containing
  members whose `rvcut.done` already existed; `one_batch_w`'s `rm -rf`
  at composition then destroyed their `cut-reads.fasta` (the
  restriction basis) when the batch guard-killed — **T02596 T04149
  T09925 T00041 T00130 T07448 T04006 T02594** left with stale done
  markers and empty basis dirs, cutA failing Os NotFound. The repair
  (`run-recovery.sh`, committed): stale markers cleared, the as-is
  rvcut re-cut **8/8 exit 0 (batch wall 599s)**, the as-is cutA re-run
  (all 8 landed; **cutA now 55/55 of the 55 runnable** — T09025
  excluded), then a bounded A pass over the recovered members.
- **THE CONCURRENT DOUBLE-WRITE CLASS (found, repaired, ruled):** the
  checkpoint's design assumption — "the bounded passes' callX.done
  markers make the master's TA/TB skip exactly what landed" — BREAKS
  when the master's TA starts while a bounded pass still drains the
  same list (TA entered at 01:10Z Oct 20 while three bounded/stale
  passes still held lanes). **T07036 callA is the one true double-write**
  (master TA + a stale lane ran it simultaneously 01:11-01:19Z; both
  exit 0; the output cleared and redone clean). **T08776 callA's
  os-error-2 failure was the stale lane's `rm -rf` of the master's
  out-dir mid-run** (redo clean). THE RULE the session leaves behind:
  never two same-form call drivers; the SIGSTOPped stage-TC driver and
  all stale bounded passes are dead; ONE stageTC (C) + the master's
  TA/TB (A/B) own the lanes; the relaunch loops never relaunch while a
  bounded pass holds a member mid-call, and cap at 6/8 relaunches
  (deterministic fail-closed members must not loop forever).
- **A NEW TERRITORY-LEVEL FAIL-CLOSED CALL CLASS, named:** "the
  territory carries overlapping records of one path" (exit 1,
  deterministic, multiple retries): **T06027 fails in ALL THREE forms**
  (A, B, C), **T05603 in A and B** (C pending) — 2 of the 55 runnable
  tranche-1 members; the class is census/territory-level, not
  form-level. Their loci are NOT assessable by any form.
- **THE DISK STATE:** `/` fell to **15G free** mid-session (the
  co-tenant's usage + both tranches' growth; receipted in
  `logs/co-tenant.log`), then the co-tenant freed ~1T back to ~988G.
  **Tranche 2's work dir lives on `/mnt/nvme2n1` (erikg-owned local
  NVMe, 2.3T free)** — isolated from co-tenant swings on `/`;
  tranche-1's artifacts untouched (the move-to-nvme2 playbook is the
  named contingency if `/` runs out again).

### 9.3 THE RECEIPTS AT THIS INCREMENT (tranche 1's growing n; the tables supersede the checkpoint's partials; snapshot ~04:15Z)

**THE WASH (rank-1 of truth-expressible; the official QV machinery):**

| stratum | form A | form B | form C |
|---|---|---|---|
| S_ORD | 98/377 (26%), 20 terr | 32/106 (30%), 6 terr | 92/338 (27%), 17 terr |
| S_SEG | 37/157 (24%), 9 | 15/70 (21%), 4 | 36/146 (25%), 9 |
| S_REP | 24/110 (22%), 8 | 12/49 (24%), 5 | 20/86 (23%), 6 |
| S_EXP | 20/80 (25%), 7 | 0 landed | 0 landed |
| S_BIND | 12/24 (50%), 1 | 0 landed | 0 landed |

Still a wash within strata (1-9 points across forms; n small per
cell). **The first S_BIND resolution receipts: T02594 and T02018
callA exit 0 with full truth assessments** — the global-frame call
RESOLVES two members of the RV campaign's exit-1 bind-ambiguity
ledger (the class census is tranche 2's 71 + tranche 1's 8).

**THE CONSERVED-ELSEWHERE RATE — THE CLASS APPEARED AT GENOME SCALE:
1 loss in 749 loci assessed (B-right-A-wrong; T01654 locus 3, S_ORD:
B calls the truth pair [0,10] at QV 50.0, A and C both call [10,77]
at QV 45.0)** — the same locus is the 1-in-749 B-right-C-wrong loss,
so the loss holds under BOTH defaults. Rate ~0.13% (1/749), against
the five-locality 1-in-~121: consistent with a RARE class; per-stratum
0 elsewhere (S_SEG 158, S_REP 110, S_EXP 80, S_BIND 24 loci, all
zero). The denominator is still growing (B is the laggard form —
10/53 territories landed; the full A-vs-B pairing comes with TB).

**THE PRUNE'S OWN VERDICT — THE GATE BAR IS VIOLATED AT THIS n:**
across **571 paired C-vs-A loci: 10 PRUNE TRUTH WINS, 2 LOSSES**
(the checkpoint's n=98 and the n=310 interim both held zero losses).
**BOTH losses are at T03901 — the heaviest far-bin member (86.2%
share, the same territory that carries 4 of the wins):** locus 8
(A calls [0,8]=truth at QV 50 rank-1; C calls [8,8] at QV 32.1 rank
35) and locus 12 (A [0,17]=truth QV 50 rank-1; C [86,106] QV 2.5
rank 3,430). The wins at the same territory are equally dramatic
(loci 4, 9: C calls the truth pair at QV 50 rank-1 where A sits at
rank 7 and rank 18,694). **At the heavy-share members the prune is a
mixed high-stakes lever: it fixes and breaks loci at the same
territory, and the loss class lives exactly where the wins live.**
QV-visible changes: 24 loci (17 higher under C, 7 higher under A).
The other named wins: T04399 l12, T09631 l4+l6, T09593 l4+l9+l11,
T02854 (S_REP) l14, T03901 l4+l9.

**THE PRUNE'S WALL VALUE AT THE HEAVY TAIL (the outstanding receipt
from the checkpoint, now landed in part):**

| pair | far-bin share | form A loci-phase | form C loci-phase | C/A |
|---|---:|---:|---:|---:|
| T03901 | 86.2% | 729.0s | 33.1s | **0.045x (−95.5%)** |
| T04420 | 83.5% | 41.9s | 38.2s | 0.91x |

The share alone does not predict the wall effect (T03901's scattered
mass burdens the loci phase 22x; T04420's does not) — but at
far-bin share > 1% the paired median is **C/A 0.78x (n=16)** vs
**1.04x at share ~0 (n=8)**: the prune's wall value at genome scale
is real, heavy-tailed, and territory-dependent. Per-locus form-C
walls: S_ORD median 0.60s, max 5,985.7s, worst-5% share 96.9% (n=341)
— the concentration class deepens at genome scale.

**THE FAR-BIN MASS SHARE — ALL FIVE STRATA NOW RECEIPTED (n=55
censuses; the census-geometry share, form C's own machinery receipts
agreeing in 24/24 cross-checked cases):**

| stratum | n | median | max | zero-share |
|---|---:|---:|---:|---:|
| S_ORD | 21 | 8.58% | 86.18% | 4/21 |
| S_SEG | 10 | 1.61% | 23.94% | 5/10 |
| S_REP | 8 | 37.88% | 80.97% | 0/8 |
| **S_EXP** | 8 | **24.79%** | 72.84% | 2/8 |
| **S_BIND** | 8 | **32.31%** | 62.13% | 1/8 |

**TRANSCHE AGGREGATE: 34.24% of the restricted census mass is
far-bin scattered** (31.7M of 92.6M mult, n=55) — UP from 23.85% at
the checkpoint's n=28, because **the universe's dominant stratum
S_EXP is far-bin HEAVY (median 24.8%, 3x the S_ORD median)**, and
S_BIND heavier still. The prune's census-mass value is therefore
concentrated exactly in the strata tranche 1 had not yet receipted.

**THE CENSUS-MASS RATIO B/A at genome scale, all five strata
(n=55):** medians S_ORD 1.14x, S_SEG 1.21x, S_REP 1.03x, S_EXP 1.08x,
S_BIND 1.13x — the ~1.1x genome-scale cost of the record form holds
in every stratum (the five-locality 2.7-15x was the 1x-subsampled
slice artifact). Cut walls (contended): cutA median 724-926s, cutB
2,371-5,135s by stratum.

**THE BOUNDED GENOME WALL — the honest state at this increment:**
form A (the old default) now has receipts in every stratum: the
stratified full-universe estimate is **3.58e6 core-s ± 3.84e6 (95%
CI, 107% relative) ≈ 3.9h ± 4.2h at 256 cores** — S_EXP itself
n=4 (mean 566s, median 164s, max 1,870s), S_BIND n=1. **Form C's
estimate at this snapshot covers ONLY the landed strata (S_ORD n=16,
S_SEG n=8, S_REP n=4; U=1,315 of 5,807): 0.59e6 core-s ± 0.40e6 over
that fragment** — S_EXP and S_BIND have ZERO form-C receipts still
(the stageTC driver works the tranche in order; S_EXP/S_BIND C calls
are mid-queue). THE WALL REMAINS INCOMPLETELY BOUNDED: the C-form
full-universe total needs the S_EXP C receipts (the pairs so far say
C ≤ A at high share; the honest C total is a projection until they
land). The marathon roster remains carried as the committed bound
(tranche 2's S_MAR block will replace it with measured walls).

### 9.4 THE EXACT RESUME POINT (everything marker-gated; the drivers resumable)

- TRANCHE 1 (`/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1/`):
  the master is mid-TA (form A; 38/53 callA.done at the snapshot);
  next it runs TB (form B — the conserved-elsewhere denominator's
  laggard), then exits; `run-recovery.sh` (running) relaunches the
  master (capped 6) until TA/TB are complete, then a final
  `python3 scripts/genome-tranche/tranche-analyze.py <WT> <out>`
  replaces §9.3's tables. NB: the TB done-marker can never fire while
  the 2 overlapping-records members fail-closed — the honest TB
  receipt is callB.done = 53 with T06027/T05603 named.
  stageTC (form C, LANES=4) owns C exclusively and loops to
  completion.
- TRANCHE 2 (`/mnt/nvme2n1/erikg/hprcv2/global-panel-call-mode/genome-tranche2/`):
  substrate 108/109 (T00762 named), reads 109/109, truth2 109 ok,
  rvcut in flight (~19/108 at the snapshot, W=8 batches under
  co-tenancy with the adaptive width guard active), then F (the
  fleet cutB wave) → A (cutA pool) → TA → TB per the master;
  `run-stageTC.sh` (LANES=4) and `run-t2-monitor.sh` (the master
  relaunch loop) are armed and running. The marathon members' walls
  or guard-kill class receipt from their call walls/RSS logs when
  TA/TB/TC reach the S_MAR block.
- NOT touched, by design: the campaign's ledgers under
  wgs35x-fleet/runtime-verification/, the product-defaults worker's
  files, vendor/syng. NO PR PUSH (the owner pushes).
