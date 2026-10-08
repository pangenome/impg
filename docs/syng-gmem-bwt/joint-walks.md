# Joint physical-walk evidence (first runnable slice)

The subsequent [automatic panel-route implementation](panel-routes.md) builds all
source lanes/topology families and explores successive donor switches, using this
unchanged finite engine as an independent replay oracle.

This **finite explicit-layout** compiler/evaluator/bounded exhaustive solver uses
shared observations jointly. It is not genome-wide candidate completion, a noisy
read model, a calibrated posterior, or permission to emit an incumbent as sequence.
There are no truth options, source priors, local-positive pruning, frequency masks,
threshold changes, or unary Genotypes/DP substitutions on this route.

## Commands and artifacts

```
impg genome-infer compile-joint-walks --panel PREFIX --layout layout.json \
  --sources panel.fa --read-lengths 150,250 --out-dir compiled
# Optional --registry-catalog legacy/catalog.json imports ALL definitions/IDs,
# never old conditional profiles. Large location values are streamed/ignored.
impg genome-infer evaluate-joint-walks --panel PREFIX \
  --compiled compiled/joint-profiles.json --sample sample.membwt \
  --haploid-depth 10 --background 0.1 --assignment assignment.json --out-dir evaluated
impg genome-infer solve-joint-walks --panel PREFIX \
  --compiled compiled/joint-profiles.json --sample sample.membwt \
  --haploid-depth 10 --background 0.1 --max-assignments 100000 \
  --max-optima 1000 --tie-epsilon 1e-9 --out-dir solved
```

Output directories must be new. Existing output handling writes authoritative
running/failed/succeeded manifests and retains partial work. The compiler writes
`joint-profiles.json`, a version/model/checksum envelope containing the full layout,
source file fingerprints, actual source-bound compiler identity, exact sparse
per-walk/per-length orientation totals and admitted start counts. No v2 exposures
are reused. FNV-1a checksums detect accidental corruption, not hostile forgery.
Optional legacy catalog provenance explicitly records its old checksum header as
unverified; the full imported input has a new raw-content fingerprint.

Evaluations/searches write `factors.json` and `result.json`. Each factor contains
canonical tokens, optional original IDs, observed count (zero included), sparse
physical terms, supporting-slot scope, constant status and background-only status.
**No factors are omitted**, including constants and unsupported sample positives.
Result provenance binds compiled/sample checksums and parameters. Evaluation rates
and losses are ordered exactly like the factor ledger. Search returns correlated
complete choice vectors, never independently combinable marginal alternatives.
Legacy threading and reconstruction reject these dedicated artifacts. No joint
sequence emission is implemented or authorized, even for a unique incumbent.

## Explicit layout contract

`layout.json` is a strict JSON object:

```json
{
  "version": 1,
  "model": "explicit-physical-linear-layout-v1",
  "panel": {"checksum_algorithm":"fnv1a64-content-v1", "sidecars": []},
  "sources": [{"id":0,"name":"A#0#chr1","length":1000}],
  "instances": [{"id":"copy-A","source":0,"start":0,"end":1000}],
  "slots": [{"id":"molecule-1","alternatives":[{
    "id":"walk-A", "topology":"linear",
    "left_endpoint":"asserted-molecule-terminus",
    "right_endpoint":"asserted-molecule-terminus",
    "pieces":[{"instance":"copy-A","start":0,"end":1000,"strand":"+"}],
    "adjacencies":[]
  }]}]
}
```

Replace the illustrative empty `panel.sidecars` with the **exact PanelIdentity**
from a freshly built sample/catalog; source IDs/names/lengths must match the panel
namespace and have unique indexed FASTA/AGC bindings. Source coordinates are
zero-based half-open and source-forward; `strand` is `+` or `-`. The source list
uses dense IDs starting at zero. Each instance permanently binds a source span.
Different IDs assert different physical copies even at identical source coordinates;
repeated descriptions/views must not mint new IDs. This is an explicit hypothesis,
not inferred dosage or certified homology.

Pieces concatenate in listed traversal order. For N pieces, exactly N-1 ordered
adjacencies `{"from":i,"to":i+1,"kind":"abut"}` must be declared. There are no
implicit gaps/continuations, routing from BED, circular/branching topology, unknown
ends, or silent clipping at source boundaries. Source ends may be extended by the
next explicit piece. Ends assert molecule termini in this hypothesis, not biological
certification. Reads longer than a molecule have exactly zero admitted starts.

One choice is required per slot. Within one walk, disjoint pieces of the same
instance are permitted (including harmless cuts); overlapping descriptions are
rejected rather than double-counted. Across selected slots, use of the same instance
is mutually exclusive, even for disjoint subranges: splitting one physical resource
among multiple molecules is not supported by this slice. Distinct resources add
exposure. Duplicate alternatives are retained as alternative descriptions, never
simultaneously summed. Physical-resource conflicts connect solver components even
without shared features.

Assignments bind the exact compiled envelope checksum:

```json
{"version":1,"model":"joint-walk-poisson-v1","compiled_checksum":"COPY_FROM_COMPILED_ENVELOPE","choices":[0]}
```

Choice indexes follow frozen `layout.slots` and each slot's `alternatives` order.
Invalid/incomplete assignments and physical conflicts fail evaluation explicitly.

## Exact exposure and solver guarantees

The provider spells start-owned cores of at most 65,536 starts with read-only right
halos, across arbitrarily many pieces. Read lengths are 1..1,048,576 bp (at most 32
compiled lengths). Both raw native orientations are extracted independently. Anchor
enter/leave events partition **every** admitted start exactly; the unchanged native
MEM selection, coordinate/content deduplication and strict subwalk pruning run once
per event range. Every surviving adjacent-node occurrence counts, including repeated
and overlapping MEM count-two events. Integer totals are checked for overflow.
The two input orientations have identical counts because reversing swaps the two
native collection passes; their totals are **averaged**, not doubled.

For histogram H(L), physical exposure is
`e = sum_L H(L)*(q_forward(L)+q_reverse(L))/2 / sum_L L*H(L)`.
For each feature, `mu = background + depth * sum_selected e`, and the loss is
`mu - observed*ln(mu)`. There is exactly one logarithm per feature. Feature closure
is all imported registry definitions, all candidate-generated triples (novel junction
features even when observed zero), and every sample-positive triple from complete
count-only BWT enumeration. Partial enumeration fails, never silently narrows scope.
Unsupported sample positives remain explicit background-only residuals.

Search enumerates the complete Cartesian product lexicographically without
pruning; infeasible assignments also consume budget. No branch bounds or unary DP
are used. Resource constraints and shared factors determine reported components;
they do not yet reduce the global enumeration. Each result separates:

- an optional incumbent (none if no feasible assignment yet);
- full search exhaustion versus budget exhaustion;
- global objective certification **only within the supplied finite universe and
  f64 objective arithmetic**; exhaustive infeasibility is explicit;
- correlated epsilon-optimum support completeness versus `max-optima` exhaustion;
- complete native context replay versus **false genome-wide candidate completeness**.

An incomplete search has no claimed lower bound. A fully exhausted feasible search
reports its minimum as the lower bound. Alternative-storage exhaustion conservatively
withholds support completeness even if later improvements could make discarded ties
irrelevant. A succeeded manifest means the command completed, not that its bounded
search certified an optimum: inspect the result status and certification fields.

## Compatibility and mandatory follow-on

No sample format, dictionary or count policy changes. The additive `observed_pairs`
API and new crate-internal token validator visibility change the existing
source-hashed v2 compiler identity. Fresh v2 compilation/legacy routes remain tested;
frozen v2 profiles are not revalidated, spoofed or reused by the joint route.

Mandatory next seam: a **truth-blind whole-panel layout/candidate builder** producing
this explicit source/instance/piece/adjacency/endpoint contract while preserving all
original candidates and physical copies, and representing repeated/off-axis routing,
cut, orientation and topology uncertainty. It must not collapse to whole-template
selection. Then replace chromosome-sized alternative enumeration with dependency-
aware event/factorized layout variables and scalable coupled search, preserving exact
start ownership, full feature closure and physical-resource constraints. Sparse
profiles/features/layout alternatives are currently held in memory; only sequence
replay is bounded by cores/halos. No chromosome-scale performance/completeness claim
is made. The finite-layout implementation passed independent read-only review
with no issues and a parent release gate of **604 tests across 11 suites**, with
zero failures or ignored tests and matching source hashes before/after execution.
The worker also supplied 191 targeted passing tests, including independent native
read/BWT and exhaustive-assignment oracles. These gates validate this implementation
slice, not genome-wide inference. Whole-panel candidate generation, scalable search
and fresh truth-blind biological controls (both QV and full-truth coverage) remain
required before claiming the coverage loss is repaired.
