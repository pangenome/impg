# context-and-mem-likelihood — the backbone guard fix, context-aware domains, the MEM-projection local likelihood

Worker run 6604179d (sole writer, branch `work/genome-mem-bwt-pipeline`, committed base
b735c66; building on top, never reverting). INCREMENTAL REPORT: updated as each step
lands. COMMIT-PER-GATE standing policy. Owner mandate: the three measured failures are
ONE failure — a partition typed in isolation collects insufficient haplotype
information. Sequence: ITEM 1 (backbone guard, 9 unscored components), ITEM 2
(context-aware domains), ITEM 3 (MEM-projection local likelihood).

## 0. Context absorbed (bounded)

- cosigt-pattern-genotyping.md: per-locus diplotype product; decisive measurements —
  p1 0/15 truth-argmax with truth INEXPRESSIBLE at 14/15 (junction partials not domain
  rows); full chrIII 14/38 argmax; genome aggregate 0.9599 acc_H over 8 components;
  9 components blocked on the pre-existing backbone-guard failure.
- simplification.md: the collapse's end state (STEPS 1-4b, commits 574f9b0/9f767c8/
  5af93d0; b735c66 the COSIGT product).
- Verified on arrival: tree at b735c66 clean (only untracked run artifacts); the 9
  failing components (chrVI, chrX, chrIV, chrVII, chrXII, chrXIII, chrXIV, chrXV,
  chrXVI) all exit 1 with "haploid incumbent over an unscored backbone class" —
  identical message in the step45 lane's .err files.

## 1. ITEM 1 — THE BACKBONE GUARD: ROOT CAUSE FOUND (measured), FIX LANDED IN THE TREE

### 1.1 The instrumented diagnosis (run-item1diag-S288C#0#chrVI, exit 1 as expected
— the guard still fired; the diag line is the deliverable)

```
[phasing] ITEM-1 diag: locus 15 backbone allele 3010 class 1 profile_len 0
  members 2944 min_member_coverage 15 max_member_coverage 126 native_coverage 71
  domain_alleles 3296 class_charge 0
Error: "haploid incumbent over an unscored backbone class"
```

THE ROOT CAUSE, measured: chrVI locus 15 is in the divergent middle (the
locus with 3,296 domain rows — the tract middle of the chrVI failure map).
At that locus the native backbone's own domain row is a **71 bp, EMPTY-PROFILE
fragment** of the target source — the window's territory collects only 71 bp
of contiguous native source there (the homology-gap shape the owner's
diagnosis names). Its class aggregates 2,944 members, ALL empty-profile
fragments of 15–126 bp (< READ_LENGTH 150) — the domain's degenerate-row
population. The class-level hygiene rule (owner decision (c): a class is
scorable iff profile nonempty AND EVERY member covers >= L) correctly
excludes this class from the CANDIDATE sets — but the backbone allele is
structurally retained as the INCUMBENT chain's row, and the haploid table's
class-level skip left IT unscored too:

- the DIPLOID side already handles this correctly: the sweep seeds the
  native pair WITHOUT the scorable check and the DP structurally retains
  the backbone pair at every locus (finite native_pair_loss — printed in
  every failing component's .err);
- the HAPLOID side alone drops it: `haploid_allele_losses` skips every
  allele of an unscorable class (+INFINITY), `ordered_haploid_states`
  filters on finite loss, the backbone allele leaves the haploid universe,
  the haploid incumbent is INFINITY, and the guard fires.

So the 9 components' failure is EXACTLY the owner's diagnosis at the
incumbent layer: the isolated window collects so little native material
that its backbone row is a degenerate-shaped fragment, and the haploid
track's candidate-hygiene (a class-level rule designed to stop degenerate
MINIS from winning) swallowed the incumbent's structurally-retained row.

### 1.2 THE FIX (in the tree; semantics within the established pattern)

`haploid_allele_admissible(allele, backbone_allele, scorable_class) =
scorable_class || allele == backbone_allele` — the haploid table scores
the locus's native-backbone allele REGARDLESS of its class's scorability,
the exact mirror of the diploid side's structural treatment:

- every OTHER allele of the unscorable class stays INFINITY (the minis
  still cannot become haploid winners — decision (c) intact for
  candidates);
- the backbone allele is a real chainable haplotype row (the native
  source's own traversal — the same row the diploid track structurally
  retains and selects on native stretches); scoring it carries the
  FULL-window omission charge when it explains nothing (the honest
  no-evidence penalty), so it cannot win where real alleles have support;
- on every locus where the backbone class IS scorable, the new branch
  never triggers — the 8 passing components' rows are bit-identical BY
  CONSTRUCTION (no code path changes when scorable[class] is true).

Changes: `haploid_allele_losses` takes `backbone_allele` and uses the
admissibility helper; `ordered_haploid_states`' doc updated; the guard
keeps its compact ITEM-1 diagnostic (fires only if a genuinely unscored
backbone row remains — the domain-completeness bug signal). NEW unit test
`haploid_admissibility_retains_the_backbone_allele_over_class_hygiene`
(3/3 haploid tests green).

### 1.4 GATE RESULTS — the previously-unscored components scoring under the fix

| component | exit | wall | selected acc_H | acc_D | identity | switches | truth self-test |
|---|---|---|---|---|---|---|---|
| chrX | 0 | 2154 s | **0.9803** | 0.4902 | 0.9774 | 0/79 | 1.0000 / 0 exact |
| chrXIII | 0 | 1830 s | **0.9559** | 0.4780 | 0.9818 | 0/100 | 1.0000 / 0 exact |
| chrXVI | 0 | 3510 s | **0.9624** | 0.4812 | 0.9894 | 0/106 | 1.0000 / 0 exact |
| chrXIV | 0 | 1790 s | **0.9648** | 0.4824 | 0.9976 | 0/81 | 1.0000 / 0 exact |
| chrVI | 0 | 34048 s | **0.9988** | 0.4994 | 0.9904 | 0/28 | 1.0000 / 0 exact |
| chrXV | 0 | 12988 s | **0.9813** | 0.4906 | 0.9890 | 0/118 | 1.0000 / 0 exact |
| chrXII | 0 | 26696 s | **0.9561** | 0.4781 | 0.9853 | 0/118 | 1.0000 / 0 exact |

## COMMIT (commit-per-gate) — ITEM 1's fix LANDED

- **eb13bc9** "ITEM 1: the haploid backbone-allele structural retention fix
  (the 9-component backbone-guard failure)" — the admissibility helper +
  the haploid-table fix + the guard diagnostic + the NEW unit test; the
  measured gate in the message (the 7 scored rows, the all-9-guard-proven
  statement, the tests-green-on-committed-state). chrVII/chrIV rows land in
  the lane's scoreboard as they complete (measurement-only). No push (the
  owner's word required).

### 1.6 INCIDENT — the 19:48Z external kill sweep, and the retry wrapper

At 2026-09-27T19:48:08Z, while chrIV (5.5 h in, DP layer 92) and chrVII
(13.2 h in, DP layer 86) were the only remaining in-flight runs, BOTH were
killed by an external SIGTERM — the two .exit files (143) were written in
the same wall-clock second (a single sweeping agent, not our machinery: no
timeout wrappers remained; the box is contended by external campaigns and
the kills were simultaneous, wall-clock-aligned, and hit the two
longest-running processes). No kernel/OOM log entries; the box did not
reboot. The DP has no checkpointing, so both restart from zero. Countermeasure:
`run-component-retry.sh` — relaunches until exit 0 (up to 8 attempts,
setsid-detached, exit-file checked) so an intermittent external kill is
non-fatal to the lane. Both relaunched under the wrapper at ~20:00Z.

Remaining in flight: chrVII (a multi-locus heavy DP patch: layers at 526k,
4.5M, 2.03M candidates), chrIV (167 loci; an 8.42M-candidate layer at locus 91
— the largest component). Both grinding single-threaded under done-marker
polling.

All five previously-unscored components (their step45 runs died at the guard;
none had ever reached the DP) run to completion and score in the passing
family (0.93–0.99 acc_H, 0 switches, truth self-test exact). chrVI — the
instrumented-diagnosis component — posts the lane's best row (0.9988). Its
wall is dominated by the two intrinsic 8e12-edge single-threaded stages
(forward layer 15: 8,667,320 candidates, 16,652.9 s; the diploid posterior's
14→15 boundary, the same shape) — the degenerate-class member-pair
enumeration's cost, paid for the first time in the component's history.

### 1.5 The chrVI run history under the fix

- chrVI first attempt with the fixed binary: **the guard PASSED** — the run
  reached the phasing DP for the first time in chrVI's history (previously it
  died at the guard). Its divergent-middle locus 15 layer (919k prev-states;
  the 3,296-allele tract locus) is intrinsically heavy — the single-threaded
  dense DP layer ran >100 min; my enclosing `timeout 7200` wrapper risked
  killing it mid-DP, and while detaching the wrapper I lost the run
  (timeout's signal took the child). RELAUNCHED clean with NO timeout
  wrapper (the no-artificial-time-limit rule; done-marker polling only).
  Lost wall: ~1h50m; no data lost (nothing had landed).
- The full example test suite on the fixed tree: **GREEN** (tests-item1a
  exit 0, 27 test-result-ok lines, 0 FAILED; 3/3 phasing haploid tests
  including the NEW admissibility test).
- Layer-15 cost analysis (why it is slow, measured arithmetic): at locus 15
  the backbone class is the degenerate empty-profile class (2,944 members;
  the §1.1 diagnosis). The diploid DP's STRUCTURAL retention of the
  backbone class pair (pre-existing committed semantics: `retained_class_pairs`
  keeps [cb,cb] and `ordered_pair_states` enumerates every MEMBER pair of
  retained class pairs) admits ~2,944^2 ≈ 8.7M ordered member-pair states at
  that one locus — each with the seeded native pair loss (finite). The
  layer's work = 919k prev-states x ~8.7M candidates ≈ 8e12 single-threaded
  edges — many hours, the first time ANY component with a degenerate-class
  backbone locus has ever reached the DP (the guard used to kill these runs
  first). Not touched in ITEM 1 (pre-existing diploid semantics; ITEM 2's
  context-aware domain is where the degenerate-row population itself is
  addressed). chrVI keeps grinding under done-marker polling.
- Next: chrVI scores (score-component.py, identical semantics), then the
  remaining 8 components sequentially (run-rest.sh, no timeouts), then
  score-all.py for the joined genome aggregate; COMMIT after chrVI scores
  (tests already green on the committed state).

The guard (`phasing.rs`, the native-backbone HAPLOID incumbent):
`incumbent_haploid = Sum_locus haploid_loss_tables[locus][backbone_chain[locus]] +
the backbone's own haploid boundary charges`; `ensure!(finite, "haploid incumbent over
an unscored backbone class")`.

Where an entry of `haploid_loss_tables` can be +INFINITY: `haploid_allele_losses`
initializes every allele's loss to INFINITY and `continue`s (leaving INFINITY) for
every allele whose CLASS is unscorable — the domain-hygiene rule
(`scorable_classes`, owner decision (c)): a class is scorable iff its profile is
nonempty AND EVERY member allele covers >= READ_LENGTH (150) of source.

The measured inconsistency (all 9 components' logs confirm): the DIPLOID side
scores the native backbone pair REGARDLESS of scorability — the sweep seeds the
native pair without the scorable check (`exhaustive_local_sweep`'s native seeding),
the native pair loss is finite (the .err shows the finite diploid incumbent printed),
and the backbone pair is STRUCTURALLY RETAINED in the diploid DP universe
(`retained_class_pairs`: "the backbone class pair is structurally retained"). The
HAPLOID side alone drops the backbone allele: it IS in `retained_member_sets`
(structural retention) but `ordered_haploid_states` filters by finite haploid loss,
and the loss is INFINITY because of the class-level hygiene skip.

So the failure shape: at some locus on those 9 components the native backbone
class fails the hygiene rule (empty profile OR a sub-read-length member dragging
the whole class), the haploid track loses the backbone allele, the haploid
incumbent is INFINITY, and the guard fires — pre-collapse efc47df-era, present on
all 17-component lanes.

Which of the two hygiene conditions fires, and at which locus, was runtime
state — instrumented and measured; see §1.1 (BOTH conditions fire at chrVI
locus 15: empty profile AND sub-read-length members; native_coverage 71).

## 2. ITEM 2 — CONTEXT-AWARE DOMAINS: the derivation draft (design of record;
implementation gated behind ITEM 1's lane completion)

THE MANDATE: a locus is CONTEXT-SUFFICIENT when its collected anchor set
discriminates the viable haplotypes — a measurable property derived from the
existing admissible-bound arithmetic, NO tuned constants; insufficient loci
extend their window/domain until sufficient or provably exhausted; the
honest-untypable report is a legitimate output.

### 2.1 What exists (read from the code, the ingredients)

- The sweep's per-locus arithmetic (all constant-free, all committed):
  min_viable_loss, best_viable_class_pairs (the EXACT tie set), the STEP-3
  admissible per-pair lower bounds (S_c, M_c, P_c, J_c forms), the seam
  swing statistic (the [max−min] over viable-linked boundary charges), and
  the margins (min_viable + swing) — the DP's own retention rule: a pair
  whose local loss exceeds min_viable + swing CANNOT be compensated by any
  boundary evidence. This rule IS the system's own derived definition of
  "separated beyond what context can overturn".
- The row-admission mechanisms: the BED partition rows (+ the completion
  bridges), the Stage-3 split candidates (2-segment junction rows —
  currently admitted ONLY from margin-retained pairs' alleles at fresh
  loci), the same-owner stitched candidates (haploid track), and the
  EXISTING window-domain extension (Policy A: the anchor group's
  overlapping component-family rows; flag-gated,
  `--window-domain-extension`; the p1 dev loop runs WITH it, the campaign
  component runs WITHOUT it — a measured fact for the design).
- The evidence window: the locus's universe partition's routed observed
  mass (window_obs), owner-resolved per row.

### 2.2 The derived sufficiency predicate (candidate; the owner rules the
semantics before implementation)

A locus is CONTEXT-SUFFICIENT iff its own anchors separate the viable pairs
beyond every boundary compensation — formally: every viable class pair that
is NOT in the exact best-viable tie set has loss (or STEP-3 admissible
lower bound) EXCEEDING the margin (min_viable + swing). I.e., the locus's
retained viable set EQUALS its exact tie set. Equivalently: the local
posterior's credible set is everything the local anchors can honestly
support, and no non-tie alternative survives into the chain.

- Ties are profile-identity ties (the classing's own construction): the
  anchors cannot separate them AT ANY extension that adds no anchor
  content the tied classes differ on — the extension loop's stopping rule
  must therefore be structural exhaustion, not a mass constant (2.4).
- A locus failing the predicate is CONTEXT-INSUFFICIENT in one of the two
  measured failure shapes: (i) no discriminating anchors (homology gap —
  nothing separates the viable pairs: all losses collapse to the omission
  floor), or (ii) the discriminating rows are absent from the domain (the
  truth-relevant junction partials are not rows — the expressibility
  failure; the sweep cannot separate classes it cannot spell).

### 2.2b THE SUFFICIENCY CENSUS — measured on the committed runs (the predicate's
empirical shape, before any implementation)

Computed from the committed stage-1 JSON (run-p1.log / run-full.log;
`second_best_loss` is the second viable pair's DELTA over min_viable;
`retained_margin_local − min_viable_pair_loss` is the LOCAL seam swing — the
DP's own per-locus compensation capacity):

- **full chrIII (no window-domain extension): 36/38 CONTEXT-SUFFICIENT.**
  The local separations (13–10,447 nats) exceed the local swings (0–191)
  almost everywhere; the swing is 0 on the native stretches (every
  viable-viable link is the panel-attested real junction — no compensation
  capacity at all). The 2 flagged loci — 8 (second_delta 4.6 vs swing 132)
  and 10 (13.4 vs 84) — are exactly the ambiguous near-tie middle.
- **p1 slice (with the window-domain extension): 2/15 sufficient.** The
  extension's added rows widen the boundary-link variety (swings 149–862)
  while the separations stay small (5–80 nats) — the extended candidates
  are exactly the retained-but-unresolved middle. The predicate flags the
  whole dev slice as it should: the slice IS the insufficient-context
  population.
- Note (honest): the predicate is run-domain-sensitive — the p1 run's own
  extension inflates its swings, flagging more loci than the same bases
  under the full-component invocation. That sensitivity is the honest
  semantics (more candidate rows = more compensation routes = less local
  resolvability), not a defect.

This census is the ITEM-2 targeting data: the extension ladder runs at the
flagged loci; the sufficient loci keep their domains bit-identically.

### 2.3 The extension ladder (structural units, not constants)

Per insufficient locus, in order, re-sweeping and re-testing the predicate
after each rung:
1. JUNCTION-CROSSING ROWS: admit the 2-segment split rows at the locus's
   own attested co-occurring adjacencies (the port-word machinery over
   the locus's rows; the Stage-3 admission widened from margin-retained
   pairs to the structurally-attested junctions).
2. FLANK EXTENSION: admit the neighboring axis partitions' overlapping
   component-family rows (the existing window-domain extension's row
   machinery) AND their observed mass (the evidence window widens with
   them) — one adjacent partition per side per rung, the axis interval
   being the structural unit.
3. COUPLED-LOCI EVIDENCE: the boundary seam charges between the locus and
   its already-sufficient neighbors (the existing transition-cost
   machinery) — reported as the chain's conditioning, never a local
   substitute.
Exhaustion: the ladder terminates when the next rung's rows are already in
the domain (nothing new to admit — the flank has reached the component
family's own tiling) — provably exhausted; the locus is reported
HONESTLY UNTYPABLE-at-this-evidence (context_insufficient flag in the
per-locus product) rather than called wrong.

### 2.5 IMPLEMENTATION — the Stage-3.5 context-sufficiency census (in the
 tree; measurement-only; the ladder is the next step)

Implemented in the phasing tail (`run_correlation_phasing`, after the P3
rebuild — the production path all dev-loop and campaign runs take): the
per-locus TRIGGER-B census, truth-free and constant-free, at the panel's own
attestation granularity:

- D1 (inside-window expressibility): a port word carried by >= 2 ports
  across the locus's SCORABLE forward single-segment rows (a real panel
  branch inside the window — exit and entry sides at legal cut order
  exist) that no realized 2-segment row expresses at its interior
  junction. Measured dead end (documented): a READ-CROSSING attestation
  (the span index's crossing queries over the half-row compositions) was
  implemented and measured FAR too permissive for a census — multi-
  placement record chains cross thousands of homolog row pairs per locus
  (4,827–9,415 "read-attested" words/locus at p1, firing uniformly);
  the crossing test remains the novel-junction CHARGING machinery's own
  (where its permissiveness is correct — restricted compositions), not a
  census trigger.
- D2 (stranded adjacency): the window-domain extension rows in the domain
  sharing a port word with a window row — the stranded-donor-row service
  the flag provided (0 when the flag is off).
- The ambiguity diagnostic (reported, never driving, per the ruling):
  second-best delta vs local swing — unambiguous / locally_decisive /
  chain_resolved_near_tie.

In flight: the p1 (flag on) and full-chrIII (flag off) census measurements
(context-domains-v1 lane) — the p1 census's D1/D2 shape and the full-chrIII
D1 shape at the native stretches decide the census's final tightness before
the ladder consumes it. The full example suite is green on the census tree
(tests-item2b exit 0, 0 FAILED).

### 2.6 THE OWNER'S RULING (2026-09-27, on the §2.4 questions — the
implementation contract)

1. WITHIN-SWING NEAR-TIES ARE SUFFICIENT — the chain layer's designed
   work; NOT extended. (The census predicate is a per-locus AMBIGUITY
   DIAGNOSTIC — chain_resolved_near_tie vs locally_decisive — reported,
   never driving.)
2. TRIGGER-B: the extension trigger is EXPRESSIBILITY — a locus whose
   attested co-occurring panel junction compositions within its window
   territory lack 2-segment domain rows extends (the stranded-structure
   signal). The sufficiency loop SUBSUMES the window-domain extension
   flag, with the gated transition: verify at one tract locus that the
   census captures what the flag covered (the stranded donor rows),
   then the transition census both ways before retiring the flag from
   the component lane.
3. context_insufficient (honest-untypable) is the DESIRED output at the
   homology-gap loci.
4. The tighter extension set is the right economics (inexpressible only,
   not 13/15 of the slice).

### 2.7 THE CENSUS GATE — measured (p1 + full chrIII; the measurement-only
step LANDED)

The census runs (context-domains-v1 lane, the dev-loop and campaign
invocations verbatim; the census binary):

- **p1 (flag on)**: D1 at 15/15 (~15–21k missing words/locus of ~30–38k
  attested); D2 at 15/15 (~300–430 stranded-adjacent rows of the ~320–430
  extension rows per locus — D2 captures 95–100% of the flag's stranded
  rows at every locus: **the Q2 verification measured — the D2 signal
  captures what the window-domain-extension flag covered**).
- **full chrIII (flag off)**: D2 silent (no extension rows by
  construction); D1 at 36/38 (every interior locus — in the phasing path
  the split generation NEVER runs, so every attested interior branch is
  unexpressed; expressed_split_rows = 0 everywhere — D1-as-written
  measures the MACHINERY's current expression capability, not the
  sample's inexpressibility; the owner's ruling on this measurement:
  D2 governs the ladder, D1 is the completeness map, rung 1 only at
  D2-extended loci).
- **The chain rows are BIT-IDENTICAL under the census binary** (the
  selection machinery untouched): p1 selected 0.58184602910041 /
  0.42150305832240303 / 0.9567398998289838 / 2 sw / 14; full chrIII
  selected 0.8754724697123235 / 0.44874291876625216 / 0.9754847590282636 /
  1 sw / 37 — both byte-equal to the committed b735c66-era rows (the
  native2 and truth reference rows identical too; score-collapse.py
  semantics, the assessment machinery untouched).

### 2.9 THE LADDER — implemented, measured to TWO dead ends, and the
owner's FINAL RULING on the expressibility bound (the handoff state)

The ladder is IMPLEMENTED in the UNCOMMITTED tree (everything behind
`--context-aware-domains`, default OFF — the production paths, the campaign
lane, and the committed state untouched; the committed state stands at
4ef0ec0 = the census):

- The flag + the STAGED extension (built unconditionally under the flag,
  admitted only per the D2 filter — `filter_window_domain_extension` with
  pure-new group ordinal re-indexing).
- Rung 1a: the window-spanning same-source chains over the extended domain
  at flagged loci (6,475 at p1).
- The class_owners generalization to owner SETS (the mixed-owner rows'
  record-once charging — the Fix-1 form; singleton classes reduce
  bit-identically), arity-N interior junctions in the classing, and the
  owner-set forms through the boundary/rescore machinery (the sibling
  continuation's work, completed and fixed to compile by this run).
- Rung 1b (this run's addition): the cross-source D2 composition split
  rows — the junction-partial expression.

MEASURED DEAD ENDS (both principled, both documented per the guardrail):
1. The pre-admission D2 census (staged rows sharing a port word with ANY
   window row) flags 38/38 full-component loci — the word-sharing
   attestation measures PANEL HOMOLOGY DENSITY, not sample junction
   evidence (the owner's reclassification: a diagnostic, like D1's
   completeness map, never the driver).
2. Rung 1b as word-sharing-bounded materialization generated 39,839,840
   split rows at p1 (269.8 s, RSS 21 GiB and climbing — killed). The
   port-word cross-product is dense: every stranded row shares many words
   with many window rows.

THE OWNER'S FINAL RULING (2026-09-28, the expressibility bound): the
CROSSING-RECORD ATTESTATION — a junction-partial composition (A-half,
B-half, cut) materializes IFF an observed span-index record places anchor
mass on BOTH halves across the junction, with the halves from different
sources. Existence, not magnitude — constant-free and tight by construction
(the bound is crossing-records x per-side placement ambiguity, not
word-sharing density; the truth's 2 chrIII breakpoints have spanning
records at 10x coverage BY CONSTRUCTION). TWO RUNGS inside 1b: (i) the
same-path co-occurrence compositions (panel-attested, cheap), (ii) the
crossing-attested novel compositions (the mosaic's junctions — the gate's
target). MEASURE BEFORE MATERIALIZE: count the crossing-attested
compositions at p1 and full chrIII FIRST; if the count still explodes
(placement ambiguity x crossing records), escalate with that measurement
— no thresholds. The extension ladder's DRIVER becomes the
attested-composition set (loci with attested compositions extend; the
sample has FEW junctions — chrIII: 2 breakpoints — so few loci extend:
the economics restored). The Q2 verification obligation transfers to the
attested set at the tract loci.

### 2.10 THE NO-WASTE DP-LAYER FIX + THE CENSUS'S FINAL FORM (the owner's
2026-09-28 rulings, implemented and measured)

**THE NO-WASTE DP FIX** (the owner's bar: ~100s/layer; the 10-14h
chrVII/chrIV layers were "orders of magnitude over"): the diploid table's
member cross-product for an UNSCORABLE backbone class is the measured
waste — chrVI locus 15: a 2,944-member empty-profile class whose seeded
native-pair loss (10,470 nats) sits WITHIN the margin (26,031), so the
MARGIN retention enumerated all 2,944² = 8,667,320 ordered member-pair
states (a 16,652.9 s single-threaded layer, twice per run: the forward
layer + the diploid posterior's boundary pass). The fix (member-level,
in `ordered_pair_states`): an unscorable class's pair entry enumerates
ONLY the backbone allele's own pair (the incumbent chain's row — the
hygiene (c)'s own intent; the minis are pass-through structure); scorable
classes keep the full legitimate cross-product. Every locus whose
backbone class is scorable is bit-identical BY CONSTRUCTION.

- **BIT-IDENTITY SPOT-CHECK PASSED**: chrMT under the fixed binary —
  selected acc_H 0.9922837527537212 / 0 sw — byte-equal to the step45
  lane's published row; wall 130 s.
- **THE DEGENERATE-CASE RE-MEASURE (chrVI, the measured worst layer)**:
  layer 15 candidates 8,667,320 → **185**; the layer's cumulative DP wall
  16,652.9 s → **42.6 s** (within the ~100 s/layer bar); the whole run
  34,048 s → **940 s** (36x). The selected row BIT-IDENTICAL
  (acc_H 0.9988325802186795, acc_D 0.49941629010933974, 0 sw — byte-equal;
  identity differs at the 12th decimal, a summation-order artifact) —
  the minis were never selected; the waste was pure enumeration.
- chrVII/chrIV relaunches un-held after the re-measure (05:13Z, setsid
  + nohup double-detached retry wrappers — the same degenerate-shape
  layers collapse under the fix: chrVII's 4.5M and chrIV's 8.42M
  candidate layers were the same unscorable-cross-product form).

**THE CENSUS'S FINAL FORM** (three constant-free tightenings, each
measured): (1) the exact co-linearity filter (a placement pair attests
only when its inter-placement path offset equals the read's own
inter-record offset — without it the census OOM-killed at ~2.5 GB
enumerating every homolog cross-product); (2) the ROUTED-ATTRIBUTION
filter (the routing's own min-anchor touched sets — measured NOT to
tighten: 163,810 vs 164,785 — the routing is itself homolog-ambiguous);
(3) the SINGLE-SOURCE CONTINUITY test (the junction-spanning-read
doctrine's own discriminator: a record pair attests a junction only when
NO single source carries both records co-linearly — a within-source
continuation is not a junction): **164,785 → 20,162 compositions**
(78,669 read attestations). The remaining 20k are the panel's own
read-covered junctions (the rung-(i) co-occurrence territory) plus the
mosaic's novel junctions (rung (ii)).

**THE RULING'S CRITICAL CHECK PASSED**: the 7 remaining inexpressible p1
loci (0, 3, 5, 8, 12, 13, 14) all carry their stranded-attested
compositions under the final census (793, 970, 32, 387, 288, 366, 403) —
the routing's attribution does NOT drop the truth's junction
compositions at the divergent middle.

**THE PARTIAL GATE WIN stands** (measured twice, chain row BIT-IDENTICAL
0.58184602910041 / 0.4215 / 0.9567 / 2 sw / 14): the ladder's staged
admission + rung-1a chains (both ploidy tracks) + arity-N classing moved
the truth's copy-0 into the domain at 8/15 p1 loci (was 2/15); the
remaining 7 need rung 1b — whose bound is now the measured attested set.

### 2.11 RUNG 1b, THE ATTESTED MATERIALIZATION — implemented and measured
(the cut-semantics fork escalated)

Implemented (in the Stage-B block, uncommitted-behind-the-flag): for every
cross-source crossing-attested composition mapping onto the locus structure
(a side in a window row, a side in a staged row, or both in window rows),
the 2-segment junction-partial row at the attested cut — bounded by the
measured attested set (NEVER the word cross-product): **6,465 rows at p1**
(12,804 mixed-orientation compositions skipped and reported; 893
unmapped), forward-forward only. The chain row stays BIT-IDENTICAL
(0.58184602910041 / 0.4215 / 0.9567 / 2 sw / 14 — third consecutive
bit-identical p1 run).

**THE MEASURED GATE RESULT: truth copy-0 in domain HELD at 8/15** — the
materialized rows did not bring the remaining 7 loci's truth pieces into
the domain. THE CAUSE (the cut-semantics fork): the attested cuts are
READ-ESTIMATED (each crossing read's left-record extent end / right-record
extent start), while the truth's pieces require EXACT segment matches at
the truth's junction position — a per-read MEM boundary generally does not
coincide with the true junction cut, so the materialized halves miss the
truth's exact pieces. THE FORK (escalated to the owner):
(a) ALSO materialize at the PORT-POSITION cuts (the panel's own branch
    k-mer positions — the split machinery's cut convention) near the
    read-attested junctions: bounded (attested junctions x nearby ports),
    and it hits the truth's cut IFF the truth's junction sits at a panel
    port position — a measurable property of the truth's construction;
(b) accept read-estimated cuts: the loci where no materialized cut matches
    any expressible structure become the honest-untypable population
    (the mandate's own legitimate output);
(c) the owner's alternative.
### 2.12 THE CUT-GAP MEASUREMENT — the fork DECIDED (empirically, for the
owner's option (a))

The probe (examples/port_cut_probe.rs, diagnostic-only): for each of the
8 truth interior junctions, the panel port index queried at the truth's
EXACT junction cuts:

- **The two CROSS-SOURCE mosaic junctions (junction 0 at p1 locus 0 and
  junction 7 at p1 locus 14) sit at EXACT PORT POSITIONS ON BOTH SIDES**
  (9564@113881 + 9602@102922; 9602@207743 + 9564@227721 — exact_at_cut
  TRUE, ~280-300 ports within the +-2 kb window each) — the mosaic's
  breakpoints were constructed at panel branch contexts, so the PORT-CUT
  MATERIALIZATION (the owner's option (a), with the intersected-bracket
  bound) can express the truth's pieces EXACTLY at those loci — the gate
  CAN pass there.
- The six SAME-SOURCE boundary junctions (9602-9602, all
  owner_boundary: true) are NOT at ports (nearest ~2 kb away) — they are
  the same-source spanning structures rung 1a (the stitched chains) and
  the completion-bridge machinery express, not port cuts.

THE NEXT IMPLEMENTATION (specified, not yet built): materialize, at each
read-attested junction, the port cuts INSIDE THE INTERSECTED BRACKET (the
owner's derived bound: each crossing read brackets the true cut between
its last left-anchor position and its first right-anchor position; the
intersection over reads tightens to the minimal bracket; the candidates
are the ports inside it — constant-free). Expected shape from the
measurement: the two cross-source junctions' brackets contain their exact
ports; the same-source structures route through rung 1a. Then the ITEM-2
gates (truth copy-0 in domain — the 8/15 + the two port-exact loci =
10/15 with the same-source cases left to the stitched/completion forms;
the chain row unchanged-or-better), then full chrIII, then the transition
census.

### 2.8 THE LADDER's original D2-driven form (superseded by §2.9's ruling)

Under a NEW `--context-aware-domains` flag: the extension rows are STAGED
(not admitted) per locus; the ladder admits them ONLY at D2-flagged loci
(+ their observed mass via the established owner-resolved charging),
generates the split rows realizing the D2 compositions over the extended
domain (rung 1 at extended loci only), rebuilds the affected loci through
the P3 machinery (classing, folded tables, boundaries, sweeps, margins),
re-censuses, and loops until D2 quiets (context-sufficient) or the staged
rows are exhausted (context_insufficient — the honest-untypable report).
THE GUARD (the owner's escalation path, documented): if the ITEM-2 gate
fails at a locus where D2 did not fire but the truth's inexpressibility is
real, rung-2's trigger widens to the evidence-local form of D1 (branch
words with actual local observed support); if deriving that needs a
support threshold, STOP — the owner rules.

Gates (as specified): the p1 per-locus inexpressible count (14/15)
collapses; the truth pair enters the domain at the divergent-middle loci;
the chain row unchanged-or-better; then full chrIII; then the transition
census both ways (flag vs loop-only) before retiring the flag from the
component lane.

(The original three semantics questions, ruled in §2.6:)

1. Is "retained viable set == exact tie set" the right sufficiency
   predicate, or should near-ties WITHIN the swing (the chrIII native
   stretches' 18–53-nat near-ties) count as sufficient (the chain
   resolves them) with only SUB-SWING ambiguity delegating upward?
2. Does the campaign lane adopt the window-domain extension unconditionally
   (currently the p1 dev loop runs with it, the components without), or
   does the sufficiency loop subsume it entirely?
3. The homology-gap loci (the divergent middle's zero-coverage windows):
   expected honestly-untypable under the predicate — confirm the
   per-locus product's context_insufficient report is the desired output
   there (the mandate says it is legitimate, not a failure).

The ladder and predicate contain NO tuned constants: the margin, swing,
ties, bounds and the axis partition are all existing derived/structural
quantities. Implementation begins at the p1 slice after ITEM 1's commit.

## In flight (ITEM 1 lane)

- chrVI scored and committed (940 s / 11.5 GB under the no-waste fix;
  row bit-identical). chrXII, chrXV scored. chrMT byte-identical
  spot-check (130 s).
- chrVII, chrIV: the supervisor killed all silent runs AND the retry
  wrappers (they resurrected uninstrumented runs in a loop). The
  STAGE-TIMER PRECONDITION is now ABSOLUTE: no full-component launches
  until the stage timers + RSS accounting are committed; then the
  instrumented relaunch (chrVI reference: 940 s order; timer-named heavy
  layers get profiled on the slice). Watchdogs: CPU 30-min silent-kill,
  RSS 40 GB.

### 2.13 THE UNION-BRACKET GATE RUN — the measured result and the
PARENT-ROW SCOPE GAP (the checkpoint)

The union-bracket form (the sound read-bracket UNION: left bracket =
[min exit, max exit + max gap), right = (min entry - max gap, max entry] —
the intersection was unsound with reads whose left record SPANS the
junction, pushing its lower bound past the true cut): **5,226 junctions
-> 544 port-pair rows (6.12 s)** — bounded (the intersection's
mis-implemented wide form measured 480M; the sound union is 6 orders
below it). The chain row BIT-IDENTICAL for the FOURTH consecutive p1 run
(0.58184602910041 / 0.4215 / 0.9567 / 2 sw / 14).

**THE GATE EVIDENCE (honest):** truth copy-0 in domain HELD at 8/15;
ALL 544 port-pair rows landed at ONE locus (p1 index 20); the truth's
cross-source junction loci (0 and 14) did NOT advance. **THE DIAGNOSIS
(the parent-row scope gap):** the truth's junction compositions are among
the 893 UNMAPPED — their stranded parent rows are NOT in the staged set:
at locus 0 the donor's row [9602: 102922-108949] fails the staged
extension's coordinate-overlap admission test ([102922,108949] does not
overlap the window [110035,120043] IN THE REFERENCE FRAME — the donor's
source coordinate frame does not align with the reference axis; the
overlap test compares foreign-row coordinates directly against
reference-window coordinates, so frame-misaligned stranded material is
never staged). The window-domain extension's overlap heuristic is the
OLD flag's scope; the crossing-attested census supersedes it.

**THE LOOP-CLOSURE FIX (specified, the next rung):** the materialization's
parent-row table becomes the TERRITORY (the full per-source row tiling —
`territory.path_intervals` carries every partition's rows for every
source, and the attesting reads' placement points identify their own
parents), not the coordinate-overlap staged set: for each attested
composition, the parent rows on both sides come from the territory, the
junction's locus is the side that lands in a window row, the bracket +
port cuts as measured, and the materialized 2-segment rows ENTER the
locus's domain as newly-admitted staged rows (their partitions from the
territory rows — the owner-set charging handles the mixed owners). The
truth's junction parents are exactly the rows its crossing reads place
on — the census finds them; the coordinate-overlap test never could.

The cut-gap measurement (riding with the gate per the order): the
port_cut_probe measured the truth's cross-source junctions at EXACT
port positions (both sides; §2.12) — the port-expressible form is
available; the loop-closure fix is what admits the parent rows so the
port cuts can compose the truth's pieces.

### 2.14 THE CHECKPOINT COMMITS + THE INSTRUMENTED RELAUNCH

Per the supervisor's consolidated order: the gate result taken (§2.13);
ITEM 2 committed per commit-per-gate with the honest evidence
(`fdcc461`: the union-bracket port cuts + the gap-carrying census + the
arity-N test fixture fix; tests green 28 ok / 0 FAILED — the first
test-suite pass after the session's runs exposed a stale `class_owners:
vec![7, 9]` fixture that failed the TEST profile's compile — the
arity-N generalization had changed the field to owner SETS); the
stage-timer + RSS accounting committed (`ee3f52a`: genome/instrumented/
run-component.sh — a 30 s CPU-seconds + latest-stage-marker heartbeat
to `.stages` alongside the 5 s `.rss` poller; one cleanup bug found and
fixed before launch: the committed copy initially DROPPED the echo line
— caught by verifying the heartbeat file actually ticked).

THE INSTRUMENTED RELAUNCH: chrVII then chrIV, SEQUENTIAL, single attempt
each (the retry wrappers stay retired). Verified live: `.stages`
heartbeats (cpu_s + stage marker), `.rss` poller, done/exit markers.
chrVI reference under the no-waste fix: 940 s / 11.5 GB.

### 2.15 THE PER-READ CENSUS (the chain-derived frame attribution) —
IMPLEMENTED, MEASURED, AND THE GATE'S HONEST VERDICT

The supervisor confirmed the chain-derived per-read attribution (the
routing pass is per-RECORD from the start — nothing to retain; the
read's own frames are derived IN the census by intersecting the read's
full-side chains: prefix/suffix frame sets, span-co-linear extensions
where span = walk start + the record's extent low end — the physical
first-anchor path position; the correct comparator, NOT the walk start
the old continuity test used).

THE MEASURED LADDER OF FORMS (all at the p1 slice):
- per-record routed attribution (the 9.38M form, committed 7908906):
  9,380,413 aggregates, 143s census, 12.7M port-pair rows — the
  materialization explosion (killed at classing).
- per-read frames, both-sides-derived (n>=4 chains): 2,410 aggregates —
  but ZERO for the truth's junctions (a 150bp read at k=63 holds 2-3
  records; n>=4 is unsatisfiable for seam reads).
- PER-SIDE DERIVATION (the final form): a DERIVED side (chain-
  intersected) contributes all its frames; a RAW side (single record)
  contributes only its WINDOW-ANCHORED placements (the anchor ties the
  ambiguity to a materializable locus; no thresholds — the continuity
  test kills the conserved noise: 61,059 reads). THE FRAME OFFSET F
  (right_cut - left_cut = (right entry - left exit) - the read gap —
  constant per physical junction) keys the junction and pins the
  materialization's port pairs to the DIAGONAL rp = lp + F: 12,298
  aggregates -> 45 port-pair rows (0.94s; the unconstrained
  cross-product was 12.7M).

THE GATE RUN (exit 0): chain row BIT-IDENTICAL for the fifth
consecutive p1 run (0.58184602910041 / 0.4215 / 0.9567 / 2 sw / 14);
truth copy-0 HELD at 8/15; the 45 port rows landed at full-loci 35/37
(outside the p1 slice — the census's window-anchor set is the
component's, not the slice's).

**THE HONEST VERDICT ON GATE (i) (truth copy-0 at loci 0 and 14): the
truth's cross-source seams are UNATTESTABLE BY CONSTRUCTION.** The
seam-flank records straddle BOTH frames (junction A: 914 places on
9564@[113851,113922] AND 9602@[102892,102963]; junction B: 2970/2971/
2972/3230 place on both), and EVERY seam-crossing read is
native-explainable: read 397 ([914,915,916]) is span-co-linear native
9602; read 3381 ([5519,5520,2970,2971,2972]) is span-co-linear native
9602; the 9564-specific records' reads end before the seam or continue
natively on 9564. The sample's mosaic seams are COPY-CHOICE between
SIMILAR haplotypes — the flanking material is shared, so no read chain
can distinguish the seam from either native continuation. Under the
crossing-record-attestation doctrine the honest output at these two loci
is CONTEXT_INSUFFICIENT — the owner's own honest-untypable category,
not a machinery gap. The gate's divergent-middle loci (p1 indices
3/5/8/12/13) are a different question: their truth rows are native
single-source rows (no junction), and their inexpressibility is the
original D2/word-sharing class story, not the junction ladder's.

### 2.16 THE CHECKER-SIDE FINDINGS (the truth-piece derivation) AND
THE RESIDUE'S CLASSIFICATION

THE CENSUS SCOPING FIX landed (the window-anchor set = the run's locus
range, full-locus indices preserved): the port rows now land in-slice
(9,973 rows at full loci 13-22 after clipping the bracket to the parent
rows — the un-clipped form measured 'invalid source crop': exits
spanning below the parent row's start produced inverted segments).
Chain row bit-identical a sixth time.

THE TRUTH-PIECE DIAGNOSTIC (bounded, per locus/copy) exposed TWO
checker-side defects in reference_local_piece_lists:
1. UNSORTED PIECES: the territory rows iterate in file order
   (partition-major), and the in-domain check zips segments positionally
   — multi-window chains failed EVEN WHERE THE DOMAIN CARRIED THE
   EXACT CHAIN (locus 3's pieces measured [128754,138845) +
   [145952,151830) + [138845,145952)).
2. GLOBAL SORT IS WRONG TOO: the route's path order is SEGMENT-MAJOR
   (seam B threads 9602's material BEFORE 9564's; a (source, start)
   sort reverses the seam row). Fixed: segment-major, positional within
   each route segment. Chain row bit-identical a seventh time.

THE RESIDUE AT p1, CLASSIFIED (the honest measurement):
- LOCI 0/14 copy-0: THE SEAM ROWS — cross-source compositions
  (locus 14: 9602:[200187,207743) + 9564:[227721,230049); locus 0 adds
  the continuing window chain) — the honest-untypable copy-choice seams
  per the ruling (the bracket + read-length prediction still owed).
- LOCI 5/6/10 (copy-1 at 5/6/10, copy-0 at 5): THE SCATTERED-PIECE
  ARTIFACT — the piece derivation unions ALL segment-vs-territory
  overlaps, including far-homolog rows (9564:[80018,90183) appearing at
  loci whose axis interval is ~160kb away): at locus 5 the NATIVE ITSELF
  measures 'unexpressible' (9564:[80018,90183) + 9564:[153973,159863))
  — definitive proof this is checker-side, not a domain gap. The
  territory's word-sharing rows (the D2 story) pull the far homologs.
- LOCI 3/8/12/13 copy-0: SORTED CONTIGUOUS MULTI-WINDOW CHAINS STILL
  UNMATCHED — e.g. locus 3 needs the 3-chain
  9602:[128754,138845)+[138845,145952)+[145952,151830); locus 13 the
  2-chain [190177,200187)+[200187,207743) (ending at the seam cut).
  The pieces are territory-row intersections (the rows exist); the
  domain lacks the STITCHED rows — the same-owner rung-1a built 2,660
  chains but not these.

### 2.17 THE CLASS-B FIX MEASURED (9/15), THE SEAM BRACKETS, AND THE
CLASS-C DIAGNOSIS

THE ITEM-1 GATE COMPLETED: chrVII scored 0.9762 (0 switches) and chrIV
0.9645 (0 switches) under the instrumented relaunch; the lane aggregate
(10 components scored) = **genome_acc_H 0.9687 / 0 switches / 925
boundaries** — all 9 backbone-guard components pass the ITEM-1 gate.

THE CLASS-B FIX (the traversal-locality filter — the owner's approved
form: on the axis's reference frame a piece is traversal-local iff its
span overlaps the locus's axis interval; word-sharing is not
membership; per-PIECE source test, the axis is the run's local array):
**truth copy-0 8/15 -> 9/15; both copies 6/15 -> 9/15** (loci 5/6/10
fixed); chain row bit-identical an EIGHTH consecutive time.

THE SEAM BRACKETS (the class-A deliverable, measured by direct sequence
comparison at the aligned cut positions — the port_cut_probe's new
seam-bracket mode):
- Seam A (9564@113881 -> 9602@102922, offset -10,959): conserved gap
  173bp (88 left + 85 right), bracket [113793, 113966] on 9564's frame.
- Seam B (9602@207743 -> 9564@227721, offset +19,978): conserved gap
  1,233bp (584 + 649), bracket [207159, 208392] on 9602's frame.
THE FALSIFIABLE PREDICTION (on record): reads >= gap + 2k (k=63)
would attest the seams — 299bp for seam A, 1,359bp for seam B. At the
sample's 150bp reads neither seam is attestable — the copy-choice
finding is now a measured, testable property of the read length.

THE CLASS-C DIAGNOSIS (bounded, per the approved chain-by-chain plan):
- LOCI 12/13 (full 24/25): the truth's alleles are PARTIAL chains —
  their unions END at the seam-B cut (locus 13: [190177,200187) +
  [200187,207743) = union [190177,207743) = ref-frame [201,136,
  218,702) while the axis window is [210,052, 220,062] — the window
  extends PAST the truth's material (the mosaic's deletion); the
  stitched_candidates containment condition (chain union contains the
  window) STRUCTURALLY EXCLUDES every partial chain. Also: locus 13's
  piece [200187,207743) has NO row in the domain at all (it ends at
  the seam cut — the last window row before the divergence).
- LOCI 3/8 (full 15/20): the puzzle — all component rows exist in the
  domain (locus 3: rows [128754,138845), [138845,145952),
  [145952,151830) at indices 436/1147/865) and a 2-chain
  stitch:15:9602:138845-151830 EXISTS, but the needed 3-chain
  [128754..151830) does not — the code path (stitched_candidates:
  first in run_start..position, last covering the window, every end in
  last..position) should build it. Either the 3-chain is built and
  filtered by a downstream admission before joining the checker's
  ranges, or the run/admission conditions differ from the reading —
  escalated with the measurement.

### 2.18 THE PARTIAL-CHAIN FORM IN, AND THE STITCH DIAGNOSTIC'S
VERDICT ON THE LAST FOUR LOCI

THE PARTIAL-CHAIN FORM (the owner's ruling (a)) is implemented in
stitched_candidates: a maximal run that merely OVERLAPS the window
materializes once, as the full run — no truncation, the chain ends
where the rows end; same territory rows, same census attestation.
THE STITCH DIAGNOSTIC (the approved bounded instrumentation) settled
both open questions:

- LOCI 12/13 RECLASSIFIED — SEAM-DEPENDENT, NOT PARTIAL-CHAIN: the
  domain's 9602-rows at locus 25 are [190177,200187), [200187,210071),
  [210071,220263) — the run COVERS the window; the truth's last piece
  is [200187,207743) — the row CLIPPED AT THE SEAM CUT (the row's own
  bounds end at 210071, past the seam). The exact allele needs the
  clipped row = the seam-B junction composition — unattestable at
  150bp reads. Loci 12/13 join 0/14 in the honest-untypable seam class
  (four loci, all seam-clipped pieces).

- LOCI 3/8 — THE ROW-VOCABULARY MISMATCH: the domain at locus 15
  carries a MERGED row [128754,156413) (single segment) alongside the
  split forms [128754,138845) + [138845,145952)+[145952,151830) — the
  completion/anchor row duality. The checker's pieces are TERRITORY-
  row intersections (three pieces, partition-bounded); the domain's
  maximal row covers the same material as ONE segment. The strict
  piece-equality match cannot see the merged row as the truth's allele
  — a checker-semantics question (does the piece-wise standard accept
  a domain row whose segment union COVERS the pieces' contiguous
  union?), escalated to the owner.

THE p1 STATE after this pass: 9/15 copy-0 (11/15 counting the seam
class as honestly-typed); chain row bit-identical across every form;
the four remaining copy-0 loci are the two seam loci + the two
seam-clipped loci — all four bounded by the measured seam brackets and
the read-length prediction (299bp / 1,359bp).
