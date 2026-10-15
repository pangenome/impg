# THE RED-TEAM AUDIT (2026-10-09) — the diplotype genotyping headline numbers, attacked

Mandate: hostile, independent audit of the campaign's headline claims on branch
`work/genome-mem-bwt-pipeline` — the yeast fleet's 880/1,132 = 77.7% truth-rank-1
perfect calls; the human HG002 claims (~78% perfect on ordinary windows at deep
coverage, 49/121 corrected at 35x, the pilot 4/22); the head-to-head claim that
likegt/COSIGT scores 0/22 where impg scores 4 PERFECT. The auditor did not build
any of the audited code and treated the campaign's own receipts as claims, not
evidence. Four attacks, each verdict stated with its own receipts.

---

## ATTACK 1 — THE LIKEGT HARNESS: **ONE REAL HARNESS BUG FOUND (the campaign's
"fair rerun" receipt of record is INVALID), but the head-to-head verdict
SURVIVES, restated on sound legs — and it now extends to the ordinary genome,
where likegt scores 1/99.**

### 1a. The harness vs likegt's own recipe — VERIFIED FAITHFUL
- The campaign ran likegt's `geno` command with its own committed conventions,
  verified in likegt's own source: `minimap2 -ax sr --secondary=no` is
  HARD-CODED in `geno.rs` (line 502) and in `hold2out.rs` (lines 1034/1191) —
  likegt's own realignment convention for both its genotyping and its
  validation paths. `gfainject` is invoked by likegt WITHOUT `--alt-hits`
  (`geno.rs: inject_bam_to_graph`) — also its own convention.
- The graph inputs were the committed pggb/FastGA locality graphs impg itself
  scores on (the same substrate both systems) — not a graph built to sabotage
  likegt. The truth rows were verified present in every graph handed to likegt
  (the truth pair is found in every full enumeration).
- The invocation (`--top 0`, full pair enumeration) is likegt's own.

### 1b. THE BUG: the committed "fair likegt rerun" genotyped on ~10–58 reads
The fair-feed protocol subsets the region BAM with
`samtools view -N fair-p$p.matching.reads.txt`. The matching-reads files carry
`samtools fastq` MATE SUFFIXES (`.../1`, `.../2` — the fetch fastqs are flat
samtools-fastq dumps), while the region BAM's QNAMEs carry NO suffix. The
`-N` name filter matched only the suffix-less residue:

- The campaign's own receipts (this audit's measurement, pilot):
  `fair-p0.bam` holds **46 records**; likegt's own log: "reads: 46 extracted,
  45 realigned". Across the pilot's 20 graphs: **10–58 reads extracted each**
  — while `fair-head-to-head.txt` reports "matched_reads 8,328–41,852"
  (those are the minimap2 match counts, not what reached likegt).
- Therefore EVERY number in the committed `fair-head-to-head.txt` — the truth
  ranks 1,770–88,142, "the truth's standing WORSE at most windows / worsened
  at 18 of 22", the top-call cosines 0.13–0.58, graph QV 0.59–3.74, "top calls
  homozygous at 6/20" — is an artifact of a starved feed and is INVALID.
  The owner's "fair-fed per the owner's correction" claim was not delivered.

### 1c. The corrected rerun (this audit, receipts `audit-v2-*`)
With suffix-stripped names the fair feed delivers (likegt's own logs):
12k–45k realigned reads per pilot graph, 23k–29k per ordinary window.

| comparison | system | truth top-1 | notes |
|---|---|---|---|
| pilot, ALL reads (campaign's committed all-reads run — SOUND: 317,698 extracted) | likegt geno | **0/22** | truth ranks 3,039–130,107 |
| pilot, corrected fair feed (this audit) | likegt geno | **0/22** | truth ranks 2,987–130,827; top cosines 0.36–0.95; at 5 windows the truth pair is the near-runner-up (Δcos ≤ 0.005) |
| ordinary windows L1–L4, corrected fair feed (this audit, 99 windows) | likegt geno | **1/99 = 1.0%** | L1 win19 only; truth ranks 1–90,656 |
| ordinary windows L1–L4 (committed impg numbers) | impg | **77/99 = 77.8%** | all 77 PERFECT |

### 1d. The COSIGT paper recipe (alt hits) — the arm likegt's geno path omits
The pre-steer standalone COSIGT run at the pilot (the campaign's own receipt,
verified: `bwa-mem2 mem -h 10000 -c 500`, `gfainject --alt-hits 10000`,
`gafpack --len-scale --weight-queries`, cosigt over the paths table) also
failed at the hotbed (top call NA18948#1+HG04204#2, cosine 0.876). This audit
additionally ran the full paper recipe at two ORDINARY windows (L1 p0/p1):
- L1 p0: top call **HG002#1 + NA21093#2** (cosine 0.9399) — the sample's real
  haplotype 1 is recovered, paired with a rival; the truth pair ranks
  **237 of 108,811** at Δcos 4e-4 behind a wall of near-identical pairs.
- L1 p1: top call HG00253#2 + HG02027#2 (no HG002 path at all); truth pair
  rank **2,389**.

### 1e. What COSIGT's own paper-validation regime says it should achieve
likegt's own validation receipts (its README/VALIDATION_REPORT): **hold-0-out
89.2% (58/65) at the HLA-F locus** — a graph of DIVERGED haplotypes with the
sample's haplotypes present; and hold-2-out FAILURES in its own results.tsv
(wrong near-twin calls for HG00096, NA20847, HG00171 when the sample's
haplotypes are held out). Nothing in the COSIGT paper regime validates the
466-near-identical-path whole-genome-panel window; this audit's measurements
above are that regime's first numbers: **the one-hit-per-read cosine pipeline
(and even the alt-hits paper recipe) cannot call the sample's diplotype at
pangenome-panel ordinary windows — it loses to walls of near-identical pairs
by ~4e-4 cosine.**

### Attack 1 verdict
WRONG: the fair-rerun receipt of record (all its specific numbers). 
SURVIVES: "likegt 0/22 at the pilot vs impg 4 PERFECT" — restated on the
sound all-reads run and this audit's corrected fair feed, and now measured at
the ordinary genome too (likegt 1/99 vs impg 77/99). The harness followed
likegt's own conventions; the failure mechanism (coverage diffusion across
near-identical panel paths) is real and internal to the graph/panel structure.

---

## ATTACK 2 — THE DENOMINATOR: **all headline rates recomputed from the
per-locus artifacts; the expressible/raw gap is 10.4 points at yeast and
~0 at human; the exclusions are honestly named but the yeast headline uses
the flattering frame; the human "ordinary genome" number carries three
sampling-frame caveats.**

Recomputed independently from `calls.jsonl.truth-qv.jsonl` artifacts (not the
aggregate scripts):

| claim | expressible frame | raw frame (all loci) | excluded |
|---|---|---|---|
| yeast fleet (17 components) | **880/1,132 = 77.7%** | **880/1,307 = 67.3%** | 175 (13.4%) |
| human pilot (1q21.2 hotbed, exact rule) | 4/18 = 22.2% | 4/22 = 18.2% | 4 segdup windows |
| human pilot (current code, alignfold) | 4/22 = 18.2% | 4/22 = 18.2% | 0 |
| human ordinary L1–L4 (deep fetches) | 77/99 = 77.8% | 77/99 = 77.8% | 0 |
| human all five (committed depths) | 81/121 = 66.9% | 81/121 = 66.9% | 0 |
| human 35x corrected cuts | 49/121 = 40.5% | 49/121 = 40.5% | 0 |

What "expressible" excludes and whether it buries failures:
- Yeast: the 175 excluded loci are the frame-gap/contig-end/subtelomeric
  class where the truth-second fold cannot be derived (SK1's rows absent from
  the window universe, multi-window partitions, contig ends). At every one of
  them the product still EMITS a call — necessarily wrong. A consumer-visible
  raw failure rate would count them: **raw 32.7% of yeast loci get a wrong or
  unverifiable call.** The exclusions do not hide isolated embarrassments;
  they are a named substrate-limitation class — but the headline "77.7%"
  should always be quoted beside "67.3% raw".
- Human: the exact-rule pilot excluded exactly the four segdup-block windows
  (the hardest windows in the locus) — the expressible rate 22.2% flattered
  by 4 points there; the committed quotes used the raw 4/22, which is honest.
  Under the current code (alignfold) all 22 windows are expressible and the
  truth still loses at the segdup block (ranks 114–118,168), so the exclusion
  no longer exists in the current product.
- The 35x and deep-fetch human frames have NO exclusions (121/121, 99/99).

Sampling-frame caveats the "78% on the ordinary genome" headline must carry:
1. **Depth**: the 77.8% is measured at 274–320x per-homolog local fetch
   depths (~8x the 35x operating point). At 35x the same 99 ordinary windows
   score **48/99 = 48.5%** (corrected cuts; the campaign's own
   depth-response verdict, confirmed here by recount).
2. **Selection**: the eligibility rule excludes ±1Mb around every GRCh38 gap
   and acen band (86.6% of autosome bp covered by the grid), and 2 of the 8
   first-pick windows were structurally REJECTED (repeat-hub refine
   explosions, 122k/167k genome-wide rows) and deterministically replaced.
   The "ordinary genome" is thus conditioned on being locality-clean; the
   rejected class's size genome-wide is unmeasured.
3. **Sample size**: 99 windows at 4 localities of a 9,963-territory grid
   (~1%); no whole-genome human receipt exists yet (the chr1 fleet is
   mid-flight at audit time).

### Attack 2 verdict
The arithmetic of every headline rate REPRODUCES exactly from the per-locus
artifacts. WRONG: nothing numeric. MISLEADING-BY-FRAME: the yeast 77.7%
headline (raw 67.3%), and the "78% on the ordinary human genome" whenever
quoted without the depth, selection and sample-size caveats above.

---

## ATTACK 3 — THE CIRCULAR GATES: **the per-record likelihood math checks out
against its documented model; the fold-merge safety property is COMPUTED, not
assumed; the --truth-qv mode cannot make a wrong call score perfect (verified
by independent aligners and by read-level allele balance) — but the
evaluation frame is structurally truth-anchored, the "PERFECT" metric is
class-identity rather than sequence accuracy, and the 264M-check sweep is a
mirror, not an oracle. One genuine crash bug found in the derive-cache
builder.**

### 3a. The per-record likelihood (realign_score.rs, hand-derived)
- Placement score: `m·a + c·b` (matches/mismatches at the reads' own
  quality-derived per-base log-rates), factorized via pinned anchor backbones
  plus skipped-base votes, with an in-process direct per-base recomputation
  asserted equal on every scored placement (the exactness receipt).
- Per-unit per-fold LL: uniform placement prior `-ln(2·(img_len-149))` +
  logsumexp over anchored placements; NEG_INF when the fold shares no anchor
  or the image clips the read.
- Generated-elsewhere marginalization: `LL(read|fold) = logsumexp(local, E)`
  with `E = 150·A - ln(2·(G-149))` — **G is the LOCUS PANEL's mean per-strain
  path length** (measured from the receipts: pilot diploid universe 532,223bp
  → E = -14.83; yeast 28.4Mb → E = -17.87), not the real 6.2Gb human genome
  (which would give E ≈ -22.6). This is a stated convention
  ("the instrument's genome universe") but its effect on the near-twin
  separations that decide the human misses is unmeasured; a smaller universe
  makes E less negative, shrinking the unplaced-differential that the
  campaign's own autopsy named as the failure mechanism.
- Class (diplotype) LL: `Σ_u count_u · mix_logsumexp(LL_A[u], LL_B[u])` — a
  proper balanced diploid mixture; exhaustive over all pairs; argmax with
  deterministic tie order. The math is coherent and matches the documented
  model.
- The extent normalization (territory images) is symmetric in form — every
  fold is clipped to what it shares with the axis row — BUT see 3d.

### 3b. The fold machinery — the safety property is computed
`alignment_connected_merge` (realign_score.rs:1781): candidate pairs are
admitted ONLY if (i) their per-unit LL columns are BIT-IDENTICAL (verified
exactly after hash bucketing — the hash is never trusted), (ii) both sit in
the axis alignment component, (iii) they share a node, AND (iv) the variant
pockets (symmetric-difference node/edge multisets) carry ZERO observed census
mass — the pocket masses are COMPUTED from the observed-mass maps, with
refusal counters emitted (window 11's 3,526 refusals on pocket mass are in
the campaign's own receipts). The operative gate is the LL-column equality,
so the merge cannot fold genuinely-distinguishable material even if the
census attribution under-observes a pocket. VERDICT: computed, not assumed.

### 3c. The --truth-qv test mode — one locus traced end-to-end
Pilot locus 0, from `calls.jsonl` + `instrument-receipt.jsonl`:
- Called folds: fold 0 (HG002#1 row [0,10057), coalesced with CHM13 + 9
  identical rows) and fold 9 (HG002#2 row [3155,13195)).
- Truth pair: fold 0 (the axis IS the truth-first by construction) + fold 9
  (the exactly-one-fold rule). The four biWFA pairs are called-vs-truth FOLD
  pairs — at a rank-1 locus they are self-alignments (fold 0 vs fold 0), so
  `edits == 0` and QV 60 BY CONSTRUCTION of a correct class pick.
- **"PERFECT" is therefore class-identity, not an independent sequence
  accuracy claim**: the emitted sequences are panel-row spells, and the
  truth is the sample's own panel rows, so a correct call is byte-identical
  by construction (perfect == rank-1 at every locus in every receipt —
  the campaign's own tables show the counts equal everywhere).
- Could a bug make a WRONG call score perfect? Only a rival fold with a
  byte-identical spelled sequence could — folds are keyed by
  (sequence, relative walk), so that is the same class. Independently
  verified (Attack 4): the four pilot perfects are byte-identical to the
  truth rows on direct extraction, and the wrong calls are NOT perfect
  (their receipt errors are if anything OVERstated: locus 12's receipt
  error 0.263 is biWFA end-to-end extent mismatch; the aligned-block
  divergence measured by minimap2-asm5 is 0.014 — the campaign's own
  autopsy distinction, confirmed).
- The QV cap: QV = 60 at error ≤ 1e-6; the metric conflates extent mismatch
  with true divergence at wrong-call loci (measured above).

### 3d. THE EVALUATION FRAME IS TRUTH-ANCHORED (the one structural circularity)
The test mode REQUIRES truth haplotype 1 == the axis lane (`ensure` at
realign_score.rs:4963 — the CLI cannot even run under a neutral axis), and
every window's territory image is the interval shared with THE AXIS ROW.
At every yeast and human receipt of record the axis IS the truth's first
haplotype (S288C, HG002#1): the truth-first fold's image is its full extent
while every rival (including truth hap 2) is clipped to its truth-shared
interval. The docs' honest-check statement "the territory is the locus's own
coordinate, not the truth's" is TRUE of the model and MISLEADING of the
evaluation: in the evaluation the locus's coordinate IS the truth's
coordinate. No receipt exists under a neutral axis; the direction of the
resulting bias is unmeasured. This does not touch the truth-free product
path, but every accuracy number in the record is measured in a frame that
favors the truth pair's first homolog.

### 3e. The checkers: mirrors, not oracles — assumption inventory
`check-realign-scoring.py` (3,064 lines, 264,718,566 checks across the sweep)
re-derives the instrument's arithmetic from the same GFAs and the SAME census
receipt, using the SAME conventions (canonical schemes, argmax
correspondence, image rules). It proves consistency and byte-identity
between runs — it cannot detect a shared-model error. Specific assumption
points: (i) the observed-mass maps are consumed, never re-measured from the
reads; (ii) the truth-second correspondence is re-derived with the same
argmax rule the instrument uses; (iii) phase-8N's zero-regression gates
compare against the campaign's own prior receipts (relative progress, not
absolute correctness); (iv) the identity gates prove determinism, not
validity. The sweep is real and impressive as a determinism/consistency
guard — it is not independent verification, which is what Attacks 1 and 4
of this audit supply.

### 3f. A genuine crash bug (found by this audit's synthetic loop)
`examples/partition_realign_likelihood.rs:697` panics
(`index out of bounds: len 10030, index 10030`) on this audit's synthetic
wgsim reads (1% error, indels) during derive-cache construction — an
unguarded segment index in the frame-run builder. The committed real-read
inputs complete (with the already-named post-write hull-validation error),
but the builder is not robust to this input class. The error-free synthetic
loop (below) completes. Audit-only finding; no fix attempted.

### Attack 3 verdict
No bug found that inflates the headline numbers through the scoring or the
QV machinery. The fold safety property is genuinely computed. The two real
findings are structural: (1) the truth-anchored evaluation frame (3d) —
unmeasured bias direction, every headline measured under it; (2) the E
universe convention (3a) — unmeasured effect on exactly the near-twin
separations that constitute the human misses.

---

## ATTACK 4 — THE EXTERNAL ORACLE: **the pilot's PERFECT calls survive
independent sequence-level and read-level verification; the wrong calls are
genuinely wrong; the truth rows are the sample's real haplotypes at read
level; the synthetic closed loop at an ordinary human locus closes at 23/25
(error-free reads) and is blocked at the errored-read variant by the 3f
crash. The yeast fleet is itself a closed loop — on ERROR-FREE simulated
reads, a caveat the record under-states.**

### 4a. Called-vs-truth by independent aligners (pilot)
Using only `locus.fa`, `calls.jsonl` and the receipt fold identities, with
minimap2-asm5 (independent of biWFA):
- The 4 PERFECT loci (0,1,2,4): called sequences byte-identical to the truth
  rows (direct extraction) — CONFIRMED.
- The near-twin misses (3,5,6,7): receipt errors 5e-5–1e-4; independent
  aligned-block divergences 1e-4–2.6e-3 — same class, real near-twins.
- The rival-dominant misses (12,13,16,17,21): genuinely wrong (never
  byte-equal), with receipt errors OVERstated vs aligned-block divergence
  (locus 12: 0.263 vs 0.014) — the biWFA end-to-end metric charges extent
  mismatch as error.

### 4b. Read-level allele support (the strongest independent check)
Pilot locus 0 (called == truth, PERFECT): the two called homologs differ by
a 1bp deletion and a 6bp T-insert. Aligning the real 317,698 reads directly
(minimap2-sr): coverage splits ~9,660 vs 9,679 reads between the two
homologs (balanced diploid); the 1bp indel shows a ~50/50 allele split in
reads spanning the site (4 of 8 sampled alignments carry `+a`); the 6T
insert is carried by 34/119 spanning reads vs h2 (homopolymer compression
limits sr alignment). **The reads — not the instrument — support the called
diplotype.** The PERFECT calls are real sample identities, not artifacts.

### 4c. The synthetic closed loop at an ordinary human locus (L1)
Truth by construction: HG002#1 (the axis lane) + HG03927#1's panel row;
202,000 wgsim reads (148bp, ~120x/homolog) through the committed census →
derive-cache → CLI chain with `--truth-qv` on the constructed pair:
- **ERROR-FREE reads (the yeast pattern): the loop CLOSES — 23/25 windows
call the constructed diplotype, every one PERFECT (QV 60); the two misses
  are genuine: locus 2 at truth rank 2 (QV 43, the near-twin class) and
  locus 22 at truth rank 69 (QV 7.5 — the same window that is the real
  sample's one rival-dominant miss at L1, rank 958/QV 7.5: a
  haplotype-independent hard window).** Receipts: `L1/synth-loop/ef-*`
  (note: the derive-cache builder again complains after the write — this
  time `multiplicities exceed the reads` — but the cache is consumable
  and the CLI completes exit 0; same builder-validation class as the
  named pre-existing condition).
- Errored reads (wgsim -e 0.01, indel rate 0.001): BLOCKED at derive-cache
  construction by the 3f panic — the loop could not be closed for the
  realistic-error variant this audit. (No fix attempted per the audit
  mandate.)

### 4d. The yeast fleet is a closed loop — on idealized reads
`prepare-balanced-diploid.py` (verified): the yeast sample is 2,439,103
reads simulated as EXACT ERROR-FREE substrings of the S288C and SK1 paths
(150bp, 15x/15x, uniform 'I' qualities, no substitutions, no indels). The
yeast 880/1,132 = 77.7% is therefore a best-case read-model number — the
product's real-error human performance is measured by the pilot and generic
fleets, and the near-coincidence of 77.7% (yeast, error-free) with 77.8%
(human, real error at 8x the depth) is a coincidence of two different
regimes, not evidence of a shared underlying accuracy. The record's
"30x = 140x = 77.7%" depth-independence claim was already refuted by the
campaign's own 35x verdict; the error-model dimension was never varied at
yeast scale.

### Attack 4 verdict
The calls survive independent verification at every level tested:
sequence (byte + independent aligner), read (allele balance), and closed
loop (23/25 at an ordinary human locus with constructed truth — the misses
are real and haplotype-independent at one of the two windows). The
material caveat discovered is 4d: the yeast headline is an error-free-read
result, and the errored human closed loop is currently blocked by a real
crash bug rather than completed.

---

## THE HONEST RESTATED NUMBERS (every one recomputed or independently
measured by this audit)

| headline | campaign claim | audit restatement |
|---|---|---|
| yeast fleet | 880/1,132 = 77.7% perfect | **CONFIRMED at the expressible frame; raw 880/1,307 = 67.3%; on ERROR-FREE simulated reads** |
| human ordinary windows | ~78% perfect | **CONFIRMED 77/99, all PERFECT — but at 274–320x fetch depths; at 35x the same windows are 48/99 = 48.5%; 4 localities of a 9,963-window grid, locality-clean selection** |
| human segdup hotbed | 18.2% (4/22, 4 PERFECT) | **CONFIRMED raw and expressible under current code** |
| human 35x | 49/121 corrected | **CONFIRMED (40.5%), raw = expressible** |
| likegt head-to-head | 0/22 vs impg 4 PERFECT | **SURVIVES, restated: likegt-as-committed 0/22 (all-reads + corrected fair feed); paper-recipe COSIGT fails at the hotbed AND at ordinary windows (truth ranks 237 and 2,389 of ~10^5 pairs); the committed FAIR-RERUN receipt is invalid (10–58 reads) and its numbers are withdrawn** |
| depth-independence | (implied by ~78% headline) | **stand refuted at human scale (77.8% at deep coverage vs 48.5% at 35x on the same windows)** |

## THE AUDIT'S OWN RECEIPTS
- `hprcv2/generic-sampling/audit-likegt-{L1,L2,L3,L4}/` — the ordinary-window
  likegt runs (v1 broken, v2 corrected: `out-v2-p*`, `v2-p*.bamrecords`),
  the paper-recipe COSIGT arm (`paper-p{0,1}/out`), and the drivers
  (`audit-run-likegt-ordinary.sh`, `audit-v2-likegt.sh`,
  `audit-run-paper-cosigt3.sh`).
- `hprcv2/hg002-locus-pilot/cosigt-run/likegt/audit-v2-fair-pilot.sh` +
  `out-v2-fair-p*` — the corrected pilot fair feed; the committed
  `fair-*` receipts stand as the bug's evidence (10–58 `bamrecords`).
- `hprcv2/hg002-locus-pilot/audit-oracle/` — the independent-alignment
  oracle inputs/outputs; the read-level allele-balance measurements are
  reproducible from `audit-oracle/h1h2.fa` + the committed fetch fastq.
- `hprcv2/generic-sampling/L1/synth-loop/` — the closed loop (diplotype.fa,
  wgsim reads, ef-census, ef-cli; the errored variant's `cache.err` panic).

## RESIDUAL RISKS / OPEN QUESTIONS FOR THE OWNER
1. The truth-anchored evaluation frame (Attack 3d): no receipt exists under a
   neutral axis; a one-component yeast rerun with a non-truth axis would
   measure the bias direction.
2. The E-branch universe convention (Attack 3a): unmeasured effect on the
   near-twin separations; a real-genome-G variant is a one-line experiment.
3. The errored-read human closed loop is blocked by the derive-cache panic
   (3f) — the crash should be fixed and the loop completed before the
   human numbers are treated as closed-loop-verified at realistic error.
4. The rejected-locality class's genome-wide frequency (2 of 8 first picks
   exploded) is unmeasured; the "ordinary genome" claim needs it.
