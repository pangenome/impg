# Research: Minimal positional substring cover for impg whole-genome inference

**Scope:** written research and mathematical adaptation only. No source edits, installs, inference runs, panel rebuilds or publication were performed. Paper results, inspected implementation behavior and proposals are distinguished below.

## Executive summary

The requested paper is **“Minimal positional substring cover is a haplotype threading alternative to Li and Stephens model,” Ahsan Sanaullah, Degui Zhi and Shaojie Zhang, Genome Research 33(7):1007–1014, July2023, DOI [10.1101/gr.277673.123](https://doi.org/10.1101/gr.277673.123)**. It covers an **already-known haplotype at aligned marker positions** with a minimum number of exact panel segments. impg instead infers latent source/genotype configurations from weighted, overlapping, unplaced MEM-substring counts; its problem is not solved by the paper's linear-time query theorem.

The useful transfer is **explicit genotype-compatible donor-block parsimony, with alternative optimal solutions retained**. The existing impg DP already minimizes local composite count loss plus source-continuation breaks: calling this a “soft cover” changes nothing mathematically. The practical recommendation is a bounded comparison of **hard genotype-compatible minimum-block threading**, a **loss-budget intermediate**, and the **unchanged soft DP**, first against exhaustive tiny recombinant/diploid fixtures and then the two frozen whole-genome controls. Fewer donor blocks is not biological truth or a posterior probability.

## Evidence and method

Research angles: (1) primary definitions, algorithms and software; (2) correspondence of inputs/indexes; (3) count-constrained, diploid and nonlocal formulations; (4) bounded validation against current code.

The child runtime had no registered web-search, browser, fetch or source-check tools. The supervisor acquired originals with parent tools; this report directly inspected the complete saved text of the [published PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538481/), [WABI2022 primary paper](https://doi.org/10.4230/LIPIcs.WABI.2022.4), and [CPM2024 follow-up](https://doi.org/10.4230/LIPIcs.CPM.2024.12). These are primary evidence, not search-result summaries. No child `source_check` validation was possible; important claims were checked against the originals, with limitations below.

Read-only local sources were `src/genome_inference/genotype.rs`, `threading.rs`; `docs/syng-gmem-bwt/quantitative-threading.md`, `whole-genome-haploid-results.md`; and the separate planning worktree's `partition-inference-contract.md` and `genotyping.typ`. Current-code claims use the first worktree; older planning documents supply contracts, not current implementation status. The results document identifies implementation9e1cc30. Parent's final steering records branch `work/genome-mem-bwt-pipeline` at696346e, draft PR243 stacked on242 and unmerged; this does not change the inspected implementation.

## 1. What the actual paper establishes

### 1.1 Identity, chronology, input and output

**Claim:** the published identity above matches the approximate requested title. The article identifies a January6,2023 bioRxiv precursor titled **“Minimal Positional Substring Cover: A Haplotype Threading Alternative to Li & Stephens Model.”** It records receipt January6 and acceptance June6,2023. Its foundational algorithm is **“Haplotype Threading Using the Positional Burrows-Wheeler Transform,” WABI2022**, DOI10.4230/LIPIcs.WABI.2022.4. The WABI title page renders Zhang's given name “Shaoije”; the2023 article uses “Shaojie.” **Sources:** [published article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538481/), [WABI](https://doi.org/10.4230/LIPIcs.WABI.2022.4). **Support: direct evidence. Confidence: high.** The precise preprint DOI and online-first day were not independently verified; July2023 is the verified journal issue date.

Input is a known string z of N aligned marker alleles and a panel X of M equal-length haplotype strings, with a prebuilt PBWT and required access structures. Binary strings motivate the method; an arbitrary-alphabet variant is described. The query is not an unordered genotype, read collection, count vector or unknown mosaic. Output is a minimum-cardinality set of exact positional panel substrings covering z, or failure when no cover exists. Cover intervals and all donor-labelled copying paths are distinct objects.

### 1.2 Definitions: minimum versus merely irredundant

A nonempty positional substring `(i,j,s)` uses **inclusive** indices and matches z iff `s[i..j]=z[i..j]` at the **same marker coordinates**. A matching substring elsewhere in a donor is not a positional match at i..j. Positional-substring equality requires equal endpoints and substring content, not identity of the full source strings.

A positional substring cover C covers every query marker; each member is positionally present in z and at least one panel string. **Overlap is allowed.** Its size is |C|; its length is the sum of substring lengths, counting overlap. “Minimal” in MPSC means **globally minimum size**, not merely inclusion-minimal. For example, if a full exact donor exists, two adjacent half-length segments form an irredundant cover but not an MPSC.

A locally maximal match cannot extend against its particular donor. A **set-maximal match** is not strictly contained in a longer positional match against any panel donor. It is not simply an arbitrary MEM emitted by a read mapper.

**Sources:** [2023 Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538481/), [WABI §§2–3](https://doi.org/10.4230/LIPIcs.WABI.2022.4). **Support: direct evidence; irredundant example is researcher illustration. Confidence: high.**

### 1.3 Objectives, algorithms and guarantees

| Result | Scope and qualification |
|---|---|
| Ordinary MPSC | Minimum number of exact positional segments. It exists iff every query allele occurs at its own position in at least one donor (WABI Claim1). |
| O(N) query | Given the prebuilt PBWT/access structures, virtual insertion obtains longest positional-match information. Backwards greedy covering repeatedly selects the longest match ending immediately before the covered suffix. An exchange/modularity lemma proves optimality. Panel construction is excluded. |
| Leftmost/rightmost MPSC | Across minimum covers, leftmost minimizes each ranked segment's start; rightmost maximizes its end. Both are obtainable in O(N). These are bounds on optimal covers, not posterior breakpoint intervals. |
| Set-maximal-only MPSC | An optimum made entirely of set-maximal matches exists and is obtainable in O(N). This depends on the exact positional-cover problem, not arbitrary source-link constraints or evidence weights. |
| Length-maximal MPSC | **First minimize segment count**, then maximize summed lengths among minima. The2023 paper gives O(N) using nested successor neighborhoods between required-region layers. It maximizes overlap, not distinct covered bases, which already total N. |
| h-MPSC | Each selected substring must occur in **at least h panel strings at that position**. This is **not** a cap of h distinct donors, sample dosage or evidence count. WABI gives O(N+h|C|);2023 improves this to O(N) using forward matching intervals with support>=h. |
| Cover properties | Minimum covers have at most two covering segments at any marker, strictly ordered distinct starts/ends, and N<=summed length<=2N. These are interval-cover properties, not ploidy/repeat-copy bounds. |

**Sources:** [WABI §§3–6 and proofs](https://doi.org/10.4230/LIPIcs.WABI.2022.4); [2023 Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538481/). **Support: direct evidence. Confidence: high within the stated assumptions.**

**Algorithm detail inspected:** WABI Algorithm1 virtually inserts z, computing its PBWT position and two neighboring divergence arrays. Divergence computation is amortized O(N), not necessarily constant work per marker. Backwards covering fails if a longest available match makes no progress. Longest-starting-match information yields a rightmost cover without constructing a reverse panel. Extending selected matches yields set-maximal-only covers. WABI §6 describes oversized-query, insufficient-support and PBWT boundary cases.

The2023 forward h-MPSC follows a PBWT interval of panel strings matching the current query segment, ending that segment before support falls below h. The article calls the algorithm's rightmost characterization a belief; **that characterization is not promoted here to a proved theorem**.

The arbitrary-alphabet full-PBWT formulation stores an O(MN|Σ|) transition array; binary-panel storage is O(MN). Constant-in-M query time does not make preprocessing or index storage independent of M, or guarantee constant-in-M output of every physical donor occurrence.

### 1.4 Solution-space representation and enumeration qualifications

The2023 paper derives **required regions**: nonempty marker intervals shared by every segment of the same rank across optima. It groups set-maximal intervals by their required regions. Successor sets between consecutive layers are nested, allowing an implicit linear-size encoding despite potentially quadratic explicit edges. This supports length-maximal optimization and counting/exploration of optimal cover alternatives.

The article states O(N) counting and O(N+S_c) enumeration, with S_c the number of output covers. **Researcher qualification:** explicitly writing every segment list costs at least the total number of emitted segment records; a cover count alone does not specify that output size. Very large solution counts also raise integer bit-complexity costs beyond unit-cost arithmetic. The article's “one-to-one” language around maximal-interval graphs and general covers must not be interpreted as an established bijection with **all donor-labelled LS paths**: different nonmaximal endpoints can extend into the same maximal interval, and several donors can support one interval. Supplement-level representation details were unavailable.

**Source:** [2023 MPSC graph section](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538481/). **Support: direct evidence for reported algorithms; output/bijection qualifications are interpretation. Confidence: high for source wording, medium for reconstruction without supplementary code.**

### 1.5 Relationship to Li–Stephens: probability, prior and exact cover

1. **Probabilistic LS:** a copying HMM models donor switches and mismatches; specified transition/emission probabilities support uncertainty inference. A Viterbi optimum is not a full posterior.
2. **Optimization interpretation:** with a fixed known query, forbidden mismatches, uniform donor treatment and uniform switch cost, minimum exact donor runs corresponds to the hard-match/minimum-switch objective. General distance-dependent transitions, mismatch penalties and allele-frequency terms are not preserved by cardinality-only MPSC.
3. **Exact cover:** overlapping donor matches can be trimmed into nonoverlapping exact threadings; multiple switch positions within an overlap can correspond to one cover. This is ambiguity of exact donor matching, not proof of an ancestral recombination event. MPSC assigns no posterior probability to its alternatives.

The paper motivates exact matching using large, low-error panels and long matches. It explicitly states ordinary MPSC cannot tolerate mismatches. Its empirical smoothing alters input alleles before matching; it does not establish a noisy-count likelihood.

**Source:** [2023 introduction/Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538481/). **Support: direct evidence plus restricted-model interpretation. Confidence: high.**

### 1.6 What the empirical results test

The paper studies UK Biobank chromosome21 (9,793 microarray sites; overall panel974,818 haplotypes). Its imputation benchmark selects1,000 British-only haplotypes, removes90% of sites, and uses859,022 remaining panel haplotypes. Its simple method votes with donors supporting adjacent cover segments and imputes only on unanimous votes; the comparator is Beagle5.4. It reports comparable accuracy on covered sites and improvements for some variants, **while explicitly noting abstention and that the covered sites are easier for Beagle too**. P-smoother flips approximately1.4% of panel alleles in the experiment.

This supports investigating cover alternatives for a particular marker-imputation task, not diploid MEM-count inference, sequence assembly, repeat dosage, impg accuracy or universal superiority over LS. Table1's numeric cells were embedded in an image not present in the acquired text; no uninspected accuracy values are quoted. Supplementary runtime curves were unavailable.

**Source:** [2023 Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538481/). **Support: direct evidence. Confidence: high for design and qualified qualitative results.**

### 1.7 Follow-up and implementation access

**“Solving the Minimal Positional Substring Cover Problem in Sublinear Space,” Paola Bonizzoni, Christina Boucher, Davide Cozzi, Travis Gagie and Yuri Pirola, CPM2024**, DOI [10.4230/LIPIcs.CPM.2024.12](https://doi.org/10.4230/LIPIcs.CPM.2024.12), uses run-length PBWT matching statistics and k-set-maximal matches (k is the original h-support parameter). It proves O(r) space for k-SMEMs and O(r+|Q|) for length-maximal MPSC with the corresponding match objects. r is total PBWT runs; large savings depend on r being small relative to panel size. Strict sublinearity is not guaranteed for every possible panel.

The follow-up **explicitly leaves the general k-SMEM query-time bound open**, conjecturing O(w log r). Do not combine its compressed space with the original full-PBWT O(N) query bound without proof.

The2023 software statement points to [genome.ucf.edu/MPSC](https://genome.ucf.edu/MPSC), `Supplemental_Code.zip`, and `Supplemental_Methods.pdf` (proofs and Algorithms1,6,8,9, among others). Parent attachment fetches returned **HTML challenges, not PDF/ZIP payloads**; a legacy OA API returned404 and official software-site fetching failed. These invalid attachments were not inspected as supplements. Thus supplementary algorithms/code were **located but not source-audited or executed**. The fully inspected WABI paper supplies foundational pseudocode/proofs, not proof that the2023 supplement was read.

CPM2024 identifies [muPBWT/k-smem](https://github.com/dlcgold/muPBWT/tree/k-smem) and [k-smem-live](https://github.com/dlcgold/muPBWT/tree/k-smem-live), plus software archive `swh:1:dir:d3467768a54423c8294abfc44f87f18705b3ed02`. Repository contents were not acquired here. Its authors report modifying baseline allocation and restarting after infeasible columns for benchmarking; these are not necessarily released-original semantics.

**Sources:** [CPM2024 §§3–4](https://doi.org/10.4230/LIPIcs.CPM.2024.12), [2023 software availability](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538481/). **Support: direct primary-paper evidence; failures are parent-reported acquisition evidence. Confidence: high for stated scope; implementation verification unavailable.**

## 2. Mapping the paper to impg

| Paper object | impg object | Transfer and limitation |
|---|---|---|
| Known phased query z | Latent genotype/bundle/mosaic inferred from C(f) | Cover validated fixed alleles, or explicitly retain candidate uncertainty/count loss. Positive substring counts do not define a query haplotype. |
| Common aligned markers | Irregular ownership groups/reference intervals and variable-length donor occurrences | A certified inference axis can support a chain, but BED groups are not homologous loci; donor coordinates need not equal reference coordinates. |
| Positional exact segment | Assembly-supported occurrence run retaining full path/coordinates/strand | Continuity/provenance transfer; a monotone run with gaps does not prove exact query sequence across gaps. |
| PBWT on aligned haplotypes | Panel syng GBWT on assembly walks | Both retrieve panel support; GBWT does not impose PBWT's shared-column geometry or imply its greedy/nested-neighborhood proof. |
| Matching interval/panel support | Sample weighted FM/MEM-BWT occurrence count C(f) | C(f) is not panel donor support h, positional placement, phase or read linkage. |
| >=h panel support | Multiple assemblies, physical copies, duplicate records | Donor frequency is not sample dosage or independent evidence; frequent-template restrictions bias against rare genuine configurations. |
| Exact overlap against z | Source/reference overlap or a stored boundary-context MEM | Coordinate overlap alone neither proves sequence compatibility nor creates joint observed evidence. |
| Alternative optimal covers | Optimal source paths/sequence-equivalence classes | Transfer the reporting principle while preserving member-specific continuation and optimal edges. |

**Support: interpretation grounded in primary definitions and local code/contracts. Confidence: high.** PBWT, GBWT and FM are not interchangeable biological models. No durable read identities or read-by-haplotype matrix are required for the proposals.

## 3. Current impg objective — inspected implementation

A bundle h in ownership group i contains **all** source occurrences sharing PanSN `sample#haplotype`; complete contig spelling is retained. Distinct `(source,start,end)` feature locations determine physical multiplicity; duplicate/overlapping BED representations do not add copies.

For globally owned local features,

\[
L_i(h)=\sum_{f\in F_i}[\mu_{hf}-C(f)\log\mu_{hf}],\quad
\mu_{hf}=d\,e(s_f)m_{hf}+\beta,\qquad
 e(s)=\frac{\sum_L n_L(L-s+1)_+}{\sum_L n_LL}.
\]

The omitted log-factorial is candidate-independent. Zero counts incur signal cost; magnitude distinguishes multiplicity. This is **working/composite Poisson loss**, not a calibrated posterior: overlap, MEM selection, feature dependence and background remain uncalibrated.

Each callable singly placed axis group supplies physical occurrence states s with emission `bundle.score=L_i(b(s))`. **`thread()` iterates all bundles, not `best_bundles`.** Local best sets serve independent call reporting/evaluation, not hard threading constraints. A state selects one placement while carrying a score for the **full multicopy bundle**; its other physical copies do not disappear or become separate mutually exclusive haploid genotypes.

Let k(t,s)=0 exactly when the full source path agrees, orientation agrees and is known, and both endpoints strictly increase (+) or decrease (−); otherwise k=1. The solver computes

\[
\min_{s_{1:n}}\sum_iL_i(b(s_i))+\lambda\sum_{i=2}^nk(s_{i-1},s_i).
\]

Exact forward/backward min-sum takes O(sum_i |S_{i-1}||S_i|) time and O(sum_i |S_i|) score memory, excluding stored state/output data. All optimal marginal states within absolute1e-8 are retained. They are **not freely combinable**; retain optimal edges for complete-path alternatives.

Chromosomes and no-call/poor-fit/repeated-axis groups reset the DP. Unknown orientations compete but never establish continuation; selected unknown or poor-fit states remain unresolved. Ties break emitted blocks, so optimization blocks differ from reported resolved blocks. Gaps/overlaps are retained, not filled with sequence. There is no FASTA or exact biological breakpoint claim.

**Sources:** local `genotype.rs::{call,exposure}`, `threading.rs::{thread,solve,continuation}` and quantitative documentation. **Support: direct evidence; complexity is inspection-based derivation. Confidence: high.**

## 4. Concrete formulations — researcher proposals, not paper results

### 4.1 Hard genotype-constrained cover

Declare admissible sets A_i before threading:

- **Bundle-compatible:** accepted bundle identities, e.g. all tied local minima passing diagnostics. Preserves current call semantics, not certified allele identity.
- **Allele-compatible:** after homology/spanning validation, all source states spelling the called allele. Different donors may qualify; token/spacing equality alone cannot certify intervening nucleotides.
- **Candidate-compatible:** a retained uncertainty set, not a fixed call.

A source-supported run is a chain of admissible physical states with continuation edges, not a MEM, a same-sample label or inferred gap sequence. Let B(s)=1+sum k per nonempty support/chromosome segment. Solve

\[
\min B(s)\quad\mathrm{subject\ to}\quad s_i\in A_i.
\]

Keep all block-minimum optima, or declare count loss as a secondary tie-breaker. Empty A_i means abstention/reset or infeasibility, not forced imputation/deletion. Allele-constrained cover covers **called alleles**; bundle/candidate-constrained cover covers **admissible source-supported placements**. On our gapped/overlapping axis, this is minimum-block threading, not automatically literal positional exact-substring cover.

The existing recurrence solves this finite chain after inadmissible states are removed and emissions set to zero, or with lexicographic `(breaks,loss)` costs. No PBWT theorem is needed. This **genuinely changes policy**: local constraints cannot be overturned to save blocks.

### 4.2 Count-loss budgets

Define delta_i(s)=L_i(b(s))-min_tL_i(b(t)) on a fixed scorable universe.

- **Per-group tolerance:** A_i={s:delta_i(s)<=delta_i^max}, preserving fit/abstention rules; minimize blocks within these sets.
- **Global budget:** minimize B(s) subject to sum_i delta_i(s)<=Delta, optionally with per-group caps. A global budget alone can sacrifice one strongly supported interval to save blocks.

Soft scalarization and a budget are not identical: weighted sums may miss unsupported points on a discrete loss/block Pareto frontier. For the local chain, DP over `(i,state,integer switch count)` finds minimum count loss for each switch count; choose the smallest count satisfying Delta. This is exact with floating losses, without loss discretization; naive cost gains a factor up to n. These tolerances are **not credible/confidence regions**.

### 4.3 Soft/evidence-weighted cover

\[
\min_s\sum_iw_iL_i(b(s_i))+\lambda B(s).
\]

With w_i=1 and unchanged states, continuity and resets, this is **current DP plus one constant per segment**. It can choose locally losing bundles. Writing `exp(-loss-lambda*breaks)` does not calibrate biological probabilities.

Distance-dependent transitions, donor-specific penalties and tract constraints change the objective. Constant boundary penalties are sensitive to partition subdivision. The paper's length-maximal **secondary** overlap reward does not justify charging overlapping sample evidence twice.

A run DAG is an alternative representation: edges carry legal state runs, provenance and costs. It is equivalent only if every original path is represented and losses are counted once. Pruning to maximal runs may lose weighted/constrained optima unless closure is proved; the paper's set-maximal theorem does not prove that closure for impg.

## 5. Counterexamples required in validation

These are constructed examples, not measured paper/project results.

### A. Parsimony overrides independent allele calls

Three sites have continuous donors X=`000`, Y=`111`:

| Site | X loss | Y loss |
|---|---:|---:|
| 1 | 0 | 4 |
| 2 | 3 | 0 |
| 3 | 0 | 4 |

Local winners X–Y–X spell010 with loss0 and two switches. At lambda2, score4 loses to X–X–X, score3; Y–Y–Y costs8. Hard010 constraints forbid000. A known-query MPSC covers010 exactly rather than changing its middle character. Current soft DP allows this kind of tradeoff.

Source constraints differ from allele constraints: with X=`000`, Y=`010`, both donors spell0 at the flanks and the called010 sequence admits a one-block Y cover. Freezing source calls X–Y–X forbids it. Bundle losses can differ despite the same focal allele because wider source/copy evidence differs.

### B. Equivalent sequence is not resolved ancestry

Identical sequence/multiplicity profiles in donors X and Y give two continuous one-block optima. Report sequence agreement separately from unresolved source identity; never declare the first-sorted donor certain. Equal emission vectors with different outgoing links cannot be merged into a state that invents continuity. Each row can have optimal marginal{X,Y}, while mixed X–Y paths are suboptimal under positive switch cost.

### C. Repeat double counting

One feature occurs at two dispersed locations, contributing expected10 per copy. Global C(f)=20 is equally explained by multiplicities(1,1), (2,0), (0,2) without discriminating context. The correct mean is10(m_1+m_2)+beta in **one factor**. Assigning observation20 independently to both loci favors two copies at each, total4, or at least double-weights evidence. One span listed twice in BED contributes one genomic occurrence; two distinct spans contribute two.

### D. `AB|ab` versus `Ab|aB`

Both have local genotypes{A,a} and{B,b}. A unique AB-spanning feature occurs once in `AB|ab`, zero times in `Ab|aB`. Its prediction depends on **joint phased boundary states**, not local dosage alone. RC-orbit canonicalization preserves this contrast if the orbit is unique. Missing context, zero exposure or additional repeats may destroy identifiability. This uses linkage within a stored substring, not fabricated linkage between MEMs/mates.

## 6. Diploid, copies, shared factors and structural complications

### Diploid dosage and source-bundle copies

K certified local haplotypes give K(K+1)/2 unordered pairs, including homozygotes. Let n_ih be chromosomal dosage with sum_hn_ih=2; physical within-haplotype multiplicity a_ihf is separate:

\[
\mu_{if}=d e_f\sum_hn_{ih}a_{ihf}+\beta_f.
\]

**Do not add two haploid scores:** logarithm of summed rates differs from sum of logarithms; background contributes once per feature. Freely refitting depth cannot identify proportional dosage profiles.

Current bundles do not certify homologous diploid loci. Summing two bundle vectors is a testable artificial mixture, not automatically a biological diploid genotype. Selecting one placement per bundle does not eliminate other copies it predicts. States must describe which chromosomal walks carry which physical occurrences, including tandem/dispersed repeats.

For a certified finite chain with local emissions, ordered states z_i=(s_i^(1),s_i^(2)) give O(K²) states and naive O(nK⁴) transitions, where K now denotes retained single-homolog **occurrence states**, not just donors. Sum continuation-break costs across homologs while preserving consistent pairing. Global homolog-label swap is symmetry, not different biology. This restricted product DP is polynomial by its own derivation, **not by transferring MPSC's known-query theorem**.

### One global factor per shared observation

For canonical feature f, enumerate predicted physical occurrences once:

\[
\mu_f(z)=\beta_f+d\sum_{o\in\mathcal O_f}e_{fo}(z_{N(o)})m_{fo}(z_{N(o)}),\qquad
\mathcal L(z)=\sum_f[\mu_f(z)-C(f)\log\mu_f(z)].
\]

N(o) includes all states determining the occurrence. Boundary terms require phased tuples; dispersed-copy contributions sum **inside one logarithm**. Genomic homolog copies multiply occurrences; duplicated BED descriptions do not. Different overlapping features remain correlated, so the score is composite.

The current caller explicitly excludes these factors. Admitting them invalidates unary emissions unless scopes really reduce to unary/local-pairwise factors. Adjacent-state factors can become DP edge costs; bounded contexts can augment states. Nonlocal shared factors require joint treatment or an explicit approximation. Unique ownership alone does not establish locality.

Complexity depends on induced width and retained hypotheses. Generic elimination can be exponential in width; no polynomial or NP-hardness assertion is made for an unspecified whole biological extension. An **exact tiny oracle** enumerates every legal phased configuration, realized occurrence and single global factor. An alternative finite ILP uses one-hot states/transitions, linearized joint occurrence indicators and one-hot total multiplicities selecting precomputed Poisson costs. Fix exposure per occurrence or enumerate its finite context-dependent values. Neither is a whole-panel scalability guarantee.

### RC, spacing, SV, missing/unplaced sequence and resets

Oriented-node/internal-distance tokens preserve anchor order and spacing, not nucleotide identity inside unanchored gaps. RC counts aggregate feature orientation; palindromes count once, repeated positions count repeatedly, separators forbid cross-MEM matches. Never concatenate overlapping MEMs into unobserved contexts. Physical donor orientation remains explicit after feature canonicalization.

Variable lengths, inversions, insertions and copy arrangements require coordinates and validated junctions. Same `sample#hap` on a different contig is not same-path continuation. Monotonicity/gap reporting does not certify a recombinant splice; sequence-cover claims require spelled overlap/junction validation. Unknown strand stays unresolved, and internal partition switches require refinement rather than exact-edge claims.

Missing features may mean no exposure, error, missing panel content or low coverage—not deletion. Retain accessory/unplaced calls and explicit unknown sequence separately from donor-covered sequence. Keep chromosome resets and repeated-axis exclusions in the initial comparison. Bridging unsupported regions is an additional imputation assumption, not a free increase in coverage or tract length.

## 7. Data structure and objective fitting the existing compact sample index

**Proposal:** retain the sample MEM-BWT unchanged; build a bounded derived **occurrence-state/evidence-factor view**, not a durable read-placement table or second expanded catalog.

- **Axis row:** group/reference coordinates, status, candidate bundle IDs, full physical occurrence handles, orientation/evidence, losses and admissibility flags.
- **Continuation relation:** exact source path, both endpoints and strand; retain gaps/overlaps and member-specific links. Optional run handles compress chains without erasing provenance.
- **Factor registry:** canonical tokens/stable ID, one C(f), span/exposure, sparse unique physical incidences, scope and ownership/exclusion reasons. Boundary scopes key phased tuples; dispersed factors remain shared objects.
- **Optimal solution graph:** scores plus optimal edges/predecessors, not freely selectable marginal states. Report sequence-equivalent and source-equivalent alternatives separately.

Existing compact calls already provide losses and provenance for the first hard/soft local-chain comparison; no new sample queries are needed. Omitted feature maps are **unavailable, not zero**; repeat/boundary fixtures derive them from original indexes or explicit constructions. C(f) supplies no molecule or placement identity.

**Equivalent:** unchanged soft loss plus switch cost on the same DAG. **New:** hard admissibility, budget/lexicographic objectives, phased/copy-aware states, correctly shared/joint factors and junction validation. Run compression is an optional engineering optimization after equivalence tests, not the recommendation. The41GB catalog is a resource limitation, not the report's storage-redesign deliverable.

## 8. Bounded experiment and recommendation

### Stage A: exact small fixtures

Freeze2–4 donors and4–8 certified chunks with explicit source coordinates and short chromosomes; cap the largest diploid case for exhaustive feasibility. Include examples A–D plus one/two known recombinant switches, one internal-chunk switch, a short true donor tract, reverse monotonicity, unknown/palindromic orientation, gaps/overlaps, tandem/dispersed repeats, duplicated BED rows, missing evidence and absent panel allele.

First test explicit weighted MEM multisets, then a tiny read-to-index simulation. Use150bp/10x error-free baseline and a predeclared small perturbation set (e.g.5x, 1% independent substitutions, one missing-evidence interval), three fixed seeds. These are stress fixtures, not calibrated biological models.

Compare (1) independent count calls; (2) unchanged count+DP; (3) hard best-bundle and validated-allele covers, clearly distinguished; (4) one fixed per-group loss tolerance; (5) exact joint diploid/shared-factor oracle on tiny cases. Use a declared small penalty grid such as0,1,10; inspect the complete tiny loss/break frontier rather than tune each truth.

**Correctness gates:** constrained solutions obey A_i; scores and complete alternatives match exhaustive enumeration; copy multiplicities ignore duplicate representations; each count enters once; AB phase has correct expectation; no chromosome/orientation/sequence-fabrication violations. Do not benchmark an unvalidated objective.

### Stage B: frozen whole-genome controls

After fixture gates, reuse S288C/SK1 candidate/count/axis artifacts and unchanged count parameters. Freeze settings and outputs before truth access; avoid rebuilding the full panel or materializing another huge feature catalog. This report does **not** launch these experiments. Whole-genome recombinant/diploid validation remains a separately authorized gate once the tiny formulation is correct.

Metrics:
- hard genotype/allele compatibility, changes from local winners and count-loss increase;
- sequence concordance separately from donor/source-interval agreement;
- abstention/coverage on sample-source and reference-axis denominators separately;
- donor block count/lengths on common assessable support, optimization versus emitted blocks separately;
- recovered/spurious switches, diploid phase-switch error and breakpoint interval/distance error;
- noisy/missing-evidence robustness, dosage error and equivalent-source ambiguity;
- wall/CPU time, RSS, state/edge/factor counts, sample queries and output size, separating loading from inference.

**Adoption rule:** prefer constrained/budgeted cover only if an explicit compatibility guarantee or predeclared error/coverage tradeoff is useful. Fewer blocks alone is not success. If soft cover reproduces current DP, document equivalence and keep DP. Maintain conservative nonlocal/boundary exclusion until exact fixtures justify admitting those factors.

### Current baseline with correct denominators

Both documented controls are error-free in-panel150bp/10x over all17 assembly paths, with the same S288C axis and fixed settings. True bundle among best:97.27% S288C /98.13% SK1 **sample-source bp**. Thread coverage87.73%/85.70% **reference-axis bp**. Source agreement99.81%/99.46% on **resolved, truth-assessable reference bp**. The bundle measure allows ties. These are not allele accuracy, hold-out validation or biological phase accuracy.

There are22 repeated-axis groups (111 intervals/553,177 reference bp), deliberately unthreaded;1,066,487 of11,896,777 catalog features are excluded. The ~41GB JSON catalog dominates end-to-end resources: do not attribute loading to the cover objective. **Source:** local `whole-genome-haploid-results.md`. **Support: direct documented evidence, not rerun here. Confidence: high for recorded results; biological generalization unestablished.**

## Contradictions and qualifications

1. **h-support:** primary definitions require>=h supporting strings per segment, not<=h total donors. Any donor-cap search summary is incorrect. [2023 h-MPSC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538481/).
2. **Compressed benchmark wording:** CPM2024 abstract/conclusion claim at least two orders of magnitude memory reduction; Table1 is more nuanced. Chr22 indexing79.01GB/1.77GB≈44.6x; k=1 querying82.31GB/0.98GB≈84x, not universally>=100x. Its body says “almost two orders”; full PBWT can query faster for larger k. These are published experiments, not impg predictions. [CPM2024 §4/Table1](https://doi.org/10.4230/LIPIcs.CPM.2024.12).
3. **Enumeration/identity:** interval-cover counting is not automatically explicit enumeration of every donor-labelled path. Output size and mapping semantics require the qualifications above. [2023 graph section](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538481/).
4. **Pseudocode audit:** WABI h-MPSC pseudocode prints a threshold involving k where surrounding semantics require support h;2023 required-region endpoint prose/equation also warrants an off-by-one check. Definitions/proofs, not blindly transcribed lines, should govern any implementation.
5. **Document chronology:** older local plans describe a proposed count pipeline; current code/newer results establish a restricted haploid implementation. Proposed diploid/shared-factor extensions remain unimplemented.

## Missing evidence and residual risks

- Actual2023 supplementary PDF/code and repository source inspection are unavailable. No implementation execution, licensing audit or supplemental benchmark verification is claimed.
- Preprint DOI, online-first day, imputation table numeric cells and supplementary runtime curves were not independently inspected.
- Certified homology/spanning candidates are prerequisites for allele-constrained diploid cover; complete physical source coverage is insufficient.
- Hard constraints preserve erroneous local calls; soft parsimony can erase short true tracts. Budgets are operational compromises, not calibrated uncertainty.
- Token/spacing matches, donor identity and nucleotide allele equality require separate assessment.
- Count mixtures may never resolve equivalent repeat placements/phase. Optimization cannot restore discarded read linkage.
- Shared-factor tractability depends on actual scopes/width; no whole-genome linear-time extension is proved here.
- Exposure/background/depth calibration, missing panel alleles/SVs and extraction RC recall may dominate errors. Frequent-template support and shortest covers favor represented ancestry. In-panel controls do not test this; strict hold-out requires rebuilding relevant training-only index artifacts, not just removing a final candidate.

## Sources

**Kept primary sources**
- Sanaullah, Zhi, Zhang (2023), [Genome Research/PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538481/), [DOI](https://doi.org/10.1101/gr.277673.123) — exact identity, definitions, algorithms, empirical design and software pointers; complete article text inspected.
- Sanaullah, Zhi, Zhang (2022), [WABI](https://doi.org/10.4230/LIPIcs.WABI.2022.4), [acquired primary PDF](https://drops.dagstuhl.de/storage/00lipics/lipics-vol242-wabi2022/LIPIcs.WABI.2022.4/LIPIcs.WABI.2022.4.pdf) — greedy proof, virtual insertion, pseudocode and boundaries; complete14-page text inspected.
- Bonizzoni et al. (2024), [CPM](https://doi.org/10.4230/LIPIcs.CPM.2024.12), [acquired primary PDF](https://drops.dagstuhl.de/storage/00lipics/lipics-vol296-cpm2024/LIPIcs.CPM.2024.12/LIPIcs.CPM.2024.12.pdf) — compressed-space results, open time bound and benchmark qualifications; complete16-page text inspected.
- Software pointers [MPSC](https://genome.ucf.edu/MPSC), [k-smem](https://github.com/dlcgold/muPBWT/tree/k-smem), [k-smem-live](https://github.com/dlcgold/muPBWT/tree/k-smem-live) — identified by the primary papers, **not source-audited here**.

**Kept local sources:** the six files listed under Evidence and method — current semantics, frozen controls and extension contracts. Acquired originals are under `/home/erikg/impg/target/experiments/mpsc-research-primary/{article.txt,wabi2022.pdf.txt,cpm2024.pdf.txt}`. The inspected parent manifest records URLs, payload hashes and PDF content types; its initial missing-`pdftotext` errors do not invalidate the subsequently supplied complete PDF text extractions.

**Rejected/deprioritized:** HTML-challenge attachment bodies mislabelled as PDF/ZIP; search summaries as evidence for definitions/performance; unrelated generic substring-cover literature; uninspected software contents; stale planning-status sentences as evidence against current code. No unverified speedup/theorem is used to recommend replacing DP.

## Next steps

If implementation is separately authorized, perform the finite hard/budget/soft comparison above. Supplement/code recovery is useful before a literal algorithm port, but is not a prerequisite for testing the independently specified finite-state objective. The present recommendation is an explicit **genotype-compatibility policy experiment**, not a new whole-genome caller or a storage redesign.
