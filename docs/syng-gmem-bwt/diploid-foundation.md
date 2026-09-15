# Diploid foundation: representation, counts and identifiability

Status: **ACCEPTED for the restricted finite fixed-two-copy model, 2026-09-14**, commit `d35b00d`. Independent review and parent reproduction passed: 712 portable tests, four focused tests, all existing configured gates, and independent reconstruction of all24 paired profiles/optima. See the [frozen experiment and results](diploid-finite-experiment.md). A/A retains two copies from one panel source; heterozygous dosage and identifiable near phase are recovered; exactly indistinguishable far-phase alternatives remain paired ties. **Automatic and general diploid inference remain unvalidated.** Next is automatic joint diploid mosaic generation, before further genome-scale development under the [structured optimization plan](structured-poisson-optimization.md). The design contract below records the accepted finite scope.

## Separate questions

1. Can a hypothesis represent two homolog copies, including two identical copies from a **single panel allele/source**, without collapsing them as duplicate descriptions or requiring a second panel alias?
2. Do their exact integer occurrence counts combine correctly under one globally coupled objective?
3. Which paired genomes are distinguishable, and which phase alternatives have exactly equal count vectors?
4. After these foundations: can automatic joint diploid search recover supported hypotheses at practical scale?

Finite success on the first three is not an automatic or genome-scale diploid caller.

## Fixed-diploid mathematical contract

For two complete homolog hypotheses and equal **per-copy** depth d:

    s_f = d * (e_f(G_1) + e_f(G_2))
    R = sum_f [s_f - C_f log(1 + s_f / beta)]

Nominal depth10 per homolog means20 relative to a haploid-length reference; realized short-fixture coverage is lower because of terminal effects (7.75 per diploid base in the accepted cases). The convention must be explicit; neither changing the depth nor adding separately evaluated haploid losses creates a diploid model. Background is added once per feature, not once per homolog. Shared signals sum inside the same logarithm.

Use checked integer per-length count addition before exposure conversion. Preserve forward/RC occurrence semantics, realized zeros, unsupported positives, exact junction/terminal replay and physical copy multiplicity through caches/deduplication. Retain depth10/background0.1/L150 controls with the declared per-copy interpretation; no fitted parameters or statistical calibration claim.

## Representation design before edits

The current public assignment uses one sample–haplotype endpoint inventory and global unit canonical source-span capacity. It is not already a diploid validator. The A/A concatenation control fails with `incomplete topology assignment`; separate valid-sized controls demonstrate within-copy capacity rejection. The new wrapper validates each whole copy separately, rather than relaxing the old validator.

Propose explicit fixed homolog/copy slots for a small known-topology experiment. Resource instances must be declared once for the whole hypothesis, never minted per segment or to escape a conflict. Specify within-copy and cross-molecule capacity and the paired validation rule. Homolog-label exchange alone is not a different phase solution.

An experimental wrapper can pair individually valid haploid assignments and sum public complete-route profiles. If used, it must explicitly define its new copy-aware finite domain—not claim the old public assignment validator accepts diploids unchanged. A per-homolog source-capacity assumption is an experimental restriction, not certified homology, inferred ploidy or a general CNV model.

Prefer addition-only experimental code and unchanged public observation machinery. Any required production/backend API or capacity change needs separate approval. Preserve the old validator's repeated-source rejection as evidence of the existing limitation, not an error to swallow. Do not fork the entire search implementation.

## Freeze finite fixtures and budgets before outcomes

- **A/A:** one panel source for A, two physical copies, exactly doubled count vector. Cache reuse cannot halve dosage. The complete DNA multiset retains multiplicity two.
- **A/a:** pooled anonymous reads, with no homolog labels or placements supplied to inference. Test joint counts and dosage, not two independently optimized haploid calls.
- **Near variants:** verify that actual within-read MEM features distinguish paired phase before requiring a phased solution. Do not assume nominal read length guarantees that the operator retained phase information.
- **Far variants:** beyond the maximum read length, construct paired AB+ab and Ab+aB with exactly equal per-length count vectors. Preserve correlated paired alternatives; independent allele sets must not manufacture unsupported Cartesian combinations.
- **Copies and orientation:** distinguish repeated descriptions from physical copies, preserve RC semantics and shared-source constraints, and show that duplicate panel aliases are not necessary evidence for A/A.

An exhaustive small paired domain is appropriate for this foundational gate. Clearly label supplied finite alternatives and domain enumeration; do not claim scalable diploid candidate discovery. Define the recipes, domain, resource limits and discrete tie tolerance before running outcomes. Preserve failed attempts; no outcome-based recipe or cap tuning.

No artificial linkage, donor prior, masking, read placement or changed count policy may be used to force phase. Freeze inferred hypotheses and score-based selection before truth assessment. Truth is not an inference input.

## Evaluation and exit

Independently check exact counts, the joint objective, the finite discrete optimum and all correlated ties. Report allele dosage, copy-weighted truth **and** query coverage, and phase only where identifiable. Use molecule multisets or bijective copy matching: one emitted A must not cover both truth copies and claim100% diploid reconstruction.

Unphased or phase-block alternatives are correct outcomes when observations cannot distinguish long-range phase. Do not invent a globally phased assembly. MEM overlap means the count loss alone provides neither calibrated confidence nor independent-Poisson statistical guarantees.

Exit: reviewed and parent-reproduced finite evidence either validates the declared fixed-diploid representation/operator/ambiguity model or identifies a concrete mismatch. Automatic diploid generation, genome-scale recovery, unknown ploidy, allele imbalance, real stock identity, noisy reads, CNV and novel topology remain separate gates.

Parent acceptance is recorded in `candidate-copy-repair/diploid-foundation-parent-v1/acceptance.json`. Worker and fresh parent feature/pair/selection ledgers are byte-identical; copy-aware assessments agree. All625 pre-existing regular files and the installed executable were unchanged during implementation/reproduction. This document's acceptance update is subsequent parent documentation only.

## Next: automatic joint diploid proposals (design first)

Construct a paired hypothesis containing a non-native mosaic from panel/count inputs, rather than supplying all candidate haplotypes. Score changes against the complete paired baseline, not two independent haploid objectives. Preserve the partner's signal, exact per-copy capacity, copy multiplicity, immutable baseline/slot identities, exploratory and compound opportunities, and correlated count-equivalent phase alternatives. Reuse existing geometry where practical without copying the entire search implementation or silently changing B1/B2 behavior. Freeze the smallest synthetic construction and ambiguity experiments, budgets and failure criteria before outcomes; finite oracle checks and copy-aware assessment remain mandatory. This is the next Gate-B design task, not authorization to claim genome-scale recovery.

The known yeast mosaic remains haploid-oriented development data. A later diploid biological claim requires a new frozen challenge. Production installation, assembly emission and general diploid-support claims remain unauthorized.
