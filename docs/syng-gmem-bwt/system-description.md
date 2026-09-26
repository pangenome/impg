# The whole-genome MEM-routed genotyping system — description of record

Date: 2026-09-25. Branch `work/genome-mem-bwt-pipeline`. This document describes the
system as it stands after the chrIII campaign: what the base data is, how evidence is
projected, what the inference machinery is, how it is evaluated, what has been measured,
and what remains open. It supersedes no prior document; it sits above them (the
chronological design record lives in `partition-local-diplotype-inference.md` and the
per-arc reports under `genome/` in the session artifacts).

## 1. Goal

Infer a sample genome from raw MEM evidence against a pangenome — without closest-reference
selection, without aligning reads — recovering all sequence, dosage (copy number) and phase
supported by the evidence, while retaining genuine ambiguity as sets rather than resolving
it by lottery. The production shape is **partition-by-partition genotyping**: the pangenome
is cut into partitions, each partition's windows are genotyped locally, and a chain layer
joins the local calls. Evaluation is genotyping-native: per-partition accuracy and switch
errors, kept separate.

## 2. Base data (the inputs everything is recomputed from — zero derived caching)

- **Panel**: `yeast235.syng` (235 strains, k63 syncmer scheme, pos64 coordinates) + GBWT.
- **Sample**: `sample.membwt` (the sample's MEM-BWT index) + the reads (assessment only).
- **Partitions**: the panel's catalog groups as BEDs (`partition-pos64-w10k-d1k-to-completion`:
  window 10 kb, step 1 kb, haplotype selection mode, 19,421 groups over 3.34 Gb of panel).
- **Axis**: `reference-axis-v1.json` — each reference path's coordinate ranges mapped to the
  ONE group the reference traversal serves (the anchor group; the candidate domain
  nevertheless consults all territory-overlapping groups — see §4.2).
- **Panel routes**: the router prior (`genome-panel-routes-rebuild-v1/routes`) —
  source-level graph structure (families/lanes/ports), partition-size independent.
- **Strand-symmetric router** (2026-09-25): the territory index/steps carry BOTH frames'
  qualifying k-mers. The k63 syncmer scheme is not rc-symmetric by itself; before this fix,
  reads matching rc-stored panel copies (inverted regions) were unplaceable and their
  shares missing. Measured yield: 3.02% of all records gained placements (strictly
  additive); 65.3% of the chrIII slice's records.

## 3. The projection (read evidence → pangenome)

The routing IS a probabilistic projection of the read set onto the pangenome:

1. **Routing**: each sample record's MEM anchors are routed through the pangenome's own MEM
   structure; the record's multiplicity is divided EQUALLY among the partitions its
   placements touch (the equal-share projection). Each record is attributed once per charge
   site (the record-once rule — a record touching k charged owners contributes m_r/t_r once
   per site, not k times).
2. **Observed maps**: per partition, the routed shares form the observed feature-support
   map the local scoring consumes.
3. **Genome-wide cross-support backgrounds** (the repetitiveness normalization): the
   per-feature background is `beta_f = base + Sum_r m_r*(1 - 1/t_r_genome)` over the
   record's genome-wide touched set — a feature present in many panel contexts must earn
   against the support the sample's other copies explain. This is the likelihood-native
   form of the panplexity/complexity idea (popular features are discounted, never
   hard-filtered; the owner's tiered design — filtering the most-repetitive tier out of
   placement entirely and measuring it as dosage — is the approved destination, not yet
   implemented).
4. **Coverage/omission (instance-level)**: the haploid track charges every candidate over
   its window's observed profile. A feature's in-window observed mass is covered iff the
   chain spells an anchor of that feature overlapping one of the RECORD'S own placement
   spans (record-level coverage — one read = one observation). The same DNA piece serving
   multiple overlapping windows is ONE instance (charged once); a repeat at distinct loci
   is distinct instances (each needs its own spelling). Uncovered observed mass costs
   `beta - C*ln(beta) + ln Gamma(C+1)` per feature (the Poisson-negative-log-likelihood
   form). This closes the "ignoring is free" loophole without size thresholds or caps.
5. **Junction pricing**: a stitch between two rows is priced by the reads that physically
   span it (spanning-read evidence); co-occurring transitions (the panel's own adjacency)
   carry the panel's linkage; unsupported novel junctions are cheap-but-starved (they earn
   nothing) — the falsification census property (chimeric junctions starve) held in every
   era's census.

## 4. The inference machinery

### 4.1 Ploidy tracks
Two tracks compete; the better track wins the selection and the ploidy statement is itself
an output. The **haploid track** (one spelled chain; the second homolog empty — a ploidy
statement, not a degenerate allele) and the **diploid track** (two real slots; merged
counts before exposure; the frozen foundation conventions, panel-(b) backgrounds,
bit-identity-tested). On the haploid chrIII test sample the haploid track wins decisively
in every run.

### 4.2 Candidate domains
Per axis window: the anchor group's full BED + every component-family row from any other
group whose interval overlaps the window's coordinates (candidate-domain completeness).
The stranded-partition defect this fixed: donor pieces partitioned into groups the axis
never consulted (the chrIII tract gaps were SK1 pieces in partitions serving zero axis
windows). Same-source stitched spanning alleles (chains of adjacent same-source rows whose
union spans the window) are admitted haploid-side. Zero-length/sub-read-length mini-alleles
are excluded as scorable candidates (domain hygiene); an empty SECOND slot is a ploidy
statement, not an allele.

### 4.3 Search
- **Stage 1**: the exhaustive local sweep — every viable allele (pair/single) per locus,
  no beam search. Margin-retained candidate sets.
- **Chain DP**: exact dynamic programming over the per-window choices with admissible
  bounds, carrying the per-window surrogate objective. k-best finalist enumeration within
  the measured optimality gap.
- **The doctrine "the DP proposes, the rescore disposes"** at chain level: the MODEL OF
  RECORD (chain-level once-per-instance omission/exoneration) re-ranks the DP's finalist
  chains exactly (bit-anchored against the decompose). The DP's surrogate is a search
  heuristic; all reported scores are exact model-of-record numbers.

### 4.4 The evaluator identity
The DP runs the SAME arithmetic as the external evaluator (the decompose machinery) —
verified at bit-identity (the campaign's 1e-7-equivalence standard, restored after every
model change). No layer may instantiate its own divergent arithmetic.

## 5. Evaluation (the genotyping scoreboard)

Two metrics, never conflated:
1. **Per-partition genotype accuracy** (`acc_H`/`acc_D`): per window/locus, the chain's
   spelled haplotypes vs the truth assignment, scored by aggregate sequence agreement
   under best-case partner assignment; primary aggregate = the fraction of the component's
   true bases called correctly.
2. **Switch errors**: sequence right, homolog assignment flipped between adjacent windows —
   a count/rate, never folded into accuracy.

Truth is read assessment-side only (private-truth/, never in production paths). The
truth-mosaic reconstruction is validated bit-exact against the recipe. Machinery
self-test: the truth chain scores 1.0000/0 switches through the whole pipeline.

## 6. Measured state (chrIII, the blind mosaic: S288C background + SK1 tract
[102922,207743) on the single hybrid molecule; haploid sample at 10x)

### 6.1 The scoreboard (the w10k base, the correct metrics)
| chain | windows | acc_H | switches | donor-block recovery |
|---|---|---|---|---|
| truth (reference) | 38 | 1.0000 | 0/37 | 1.000 |
| native2 (all-native baseline) | 38 | 0.8998 | 0/37 | 0.682 |
| **v5 diploid-DP pair** | 15 | **0.9166** | 3/14 | 0.86/0.76 |
| **chain-v5 pair** | 15 | **0.8577** | 1/14 | 0.85/0.85 |
| v3-era best haploid | 38 | 0.8397 | 3/37 | 0.689 |
| model-of-record-era selected chains | 15 | 0.21–0.49 | 0–1/14 | 0.15–0.46 |

### 6.2 The findings
- **The genotyping problem is largely solved by the machinery's best chains** —
  per-partition accuracy high, switch count modest (max 3/14 = 21%, never the dominant
  failure mode). Coherence behaves as a phasing layer on an essentially-correct local
  genotyper.
- **The objective's selections anti-correlate with genotyping truth**: the model's chosen
  chains are the worst genotypers (44.7%, below the 68.2% native baseline on the donor
  block); the truth-like chain is the objective's worst. Three rounds of honest correction
  (coverage pricing, instance exoneration, strand symmetry) each WIDENED the M1-vs-pooled
  scorer opposition (318k → 597k → 798k model-prefers-patchwork; 416k → 402k → 664k
  pooled-prefers-truth). The objective cannot SEE per-partition allele correctness — it
  rewards explaining observed mass, which any sufficiently-similar allele does.
- **Failure sites** (concentrated, not general): (a) the divergent middle (axis ~154k–200k)
  where S288C<->SK1 homology breaks down (multi-kb anchor gaps) — chains spell holes or
  call native; (b) locus 23, a consistent wrong-allele miss; (c) window-spell placement
  drift — the domain extension lets a window's spell wander 1–2 loci off its own interval
  (nothing pins a spell to its window).
- **Partition size**: w2k tested (partition size the only variable) — decisively worse
  (0.584 vs 0.9166 on the same coordinates; the divergent middle collapses to 0.000–0.025).
  The 10kb window's spanning context was LOAD-BEARING: it bridges the divergent middle's
  homology gaps, which 2kb windows fall inside (their syng chains fail the anchor filters;
  catch-all universes poison the DP; sliver windows are a new failure class). One genuine
  gain (locus 23: 0.430 -> 0.564) does not change the verdict. 10kb stands; w5k untested.

### 6.3 The scorer disagreement (documented, now decoupled from bugs)
The M1 per-window factorized objective and the pooled read-level molecule scorer have
structurally different optima (measured three ways, above). The pooled credits global
spelled mass (count-arbitrage question open: 144 = 21x7 crediting); M1 charges
per-instance coverage. This is an ARCHITECTURAL question — the objective's structure vs
the read-level oracle — not a bug chain: the strand asymmetry (fixed), the union-sum
double-count (fixed), the unit errors (fixed) are all closed; what remains is the
objective's shape.

## 7. Design laws (established by measurement, in force)

1. **Constants discipline**: no tuning constants in any path; data-derived quantities
   (multiplicity, switch densities, IQR-style thresholds) are derived and reported. Wanting
   a magic number = STOP and escalate.
2. **Zero caching**: base data only (panel, sample index, axis/BEDs, routes); everything
   derived is recomputed per run.
3. **The sufficiency principle**: score only what can discriminate; documented collapse
   cases (the three measured collapses) as the guide.
4. **The block model**: blocks are the graph's own segments; a read-to-block index is the
   missing base type (fixes the API gap, evidence compression, junction observability).
5. **Measurement first**: every change is gated on per-feature arithmetic BEFORE any rerun;
   the ladder runs slice-first; full-chrIII only per validated stage; falsification
   censuses (chimeric junctions starve) must hold.
6. **Ties are sets**: count-equivalent/near-margin candidates retained as sets, never
   resolved by lottery.
7. **Un-conflated evaluation**: genotype accuracy and switch errors are separate metrics;
   sequence-reconstruction scores (the pooled external) are diagnostics, not gates.
8. **The DP proposes, the rescore disposes**: search objectives are heuristics; the model
   of record ranks finalists exactly.

## 8. Open items (the owner's queue)

1. **The selection criterion** (the biggest): the objective cannot see per-partition allele
   correctness. Candidates: (a) read-level (pooled-style base-level agreement between each
   window's spell and the reads landing there) folded into the local objective; (b) the
   pooled-form scorer as the selection criterion with M1 as evidence collector; (c) the
   owner's molecule-coherence/switch-pricing design (sequence-continuity-classified joins;
   unsupported discontinuous switches priced by the panel's own switch rate) — documented,
   unimplemented.
2. **Window-spell placement drift**: pin each window's spell to its own interval (the
   domain extension's wandering).
3. **The divergent middle**: the homology-gap region — a homology-structure limitation;
   smaller windows make it worse (measured).
4. **The tiered projection** (approved destination): highly-repetitive MEMs filtered from
   the placement channel, measured as dosage (copy number); panplexity-style auto-threshold
   on the measured multiplicity distribution.
5. **The genome-wide re-baseline**: the genome-universe index is symmetric in code; the
   campaign re-run is pending before any genome-scale campaign.
6. **The diploid validation sample**: the diploid machinery has only ever been exercised on
   a haploid sample; a genuinely diploid test sample is an owner-level test-plan decision.
7. **Durability**: the entire implementation remains uncommitted on this disk.
8. **The pooled external's count-arbitrage audit** (144 = 21x7 crediting) — queued.
