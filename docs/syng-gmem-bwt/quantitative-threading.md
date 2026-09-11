# Experimental haploid partition counts and source-path threading

`genome-infer genotype` is a quantitative extension, not the old presence-only
`call`/`run` model. Those bootstrap commands are unchanged. Catalog acceptance
remains false: source accounting is not validated homology, copy structure or
biological genotyping. There is no FASTA, calibrated posterior, diploid search,
read linkage, global strain ranking or nonlocal/boundary-factor inference.

Measured full-genome S288C and fixed-reference SK1 controls are recorded in
[whole-genome haploid results](whole-genome-haploid-results.md), with explicit
coverage/accuracy denominators and remaining acceptance gates.

## Reuse existing artifacts (no catalog rebuild)

```sh
export CFLAGS=-I/home/erikg/.cache/impg/native-build-deps/root/usr/include
export LIBRARY_PATH=/home/erikg/.cache/impg/native-build-deps/root/usr/lib/x86_64-linux-gnu
export LD_LIBRARY_PATH="$LIBRARY_PATH"
export CARGO_TARGET_DIR=/home/erikg/impg/target
export CARGO_BUILD_JOBS=4 RAYON_NUM_THREADS=4
cargo build --offline --release --bin impg -j4

BASE=/home/erikg/yeast/genome-mem-bwt-bootstrap-20260911T055203Z
IMPG=/home/erikg/impg/target/release/impg
# AXIS must be independently declared; OUT must not exist and its parent must exist.
AXIS=/absolute/path/to/reference-axis.json
OUT=/absolute/path/to/new-quantitative-thread-run
/usr/bin/time -v "$IMPG" genome-infer genotype \
  --panel /home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng \
  --sample "$BASE/sample/sample.membwt" \
  --catalog "$BASE/catalog/catalog.json" \
  --haploid-depth 10 --background 0.1 --max-mean-deviance 10 \
  --axis "$AXIS" --switch-penalty 10 \
  --allow-unvalidated-catalog --out-dir "$OUT"
```

Full-panel execution is parent-only. Optional `--truth truth.json` uses the
existing version-1 truth schema, opened only AFTER `calls.json` and `threads.json`
are frozen. Omit `--axis` for independent quantitative partition calls alone.
Ploidy other than 1 is rejected. No new panel loading/extraction is needed beyond
content fingerprinting when reusing the catalog and sample. JSON reads remain
buffered. New-directory reservation, synced atomic writes, and authoritative
running/succeeded/failed manifests are retained; handled failures remove both
calls and threads as well as evaluation. Model names now match the actual stage.

## Axis schema, version 1

```json
{
  "version": 1,
  "coordinate_system": "independently-declared-S288C-reference-v1",
  "intervals": [
    {
      "component": "S288C#0#chrI",
      "start": 0,
      "end": 10000,
      "group": "EXACT_CATALOG_GROUP_ID",
      "reference_occurrence": 123,
      "reference_strand": "+",
      "orientations": {"456": "-"}
    }
  ]
}
```

This is a schema example, not a real group/occurrence assertion. Unknown fields
are rejected. `reference_occurrence` is the catalog's **global numeric occurrence
ID**, not a group-local index. `orientations` is optional (default `{}`); it maps
catalog occurrence IDs in that group to externally established source traversal
strands relative to increasing axis coordinates. Do not populate it from truth.
A reference self-override must agree with `reference_strand`.

Intervals are supplied in order. Each nonempty full component identifier appears
in one contiguous run, with strictly increasing start AND end, and a single exact
reference source path throughout. Coordinates are zero-based half-open and at
most i64::MAX. Gaps and nonnested overlaps are allowed and retained explicitly.
Axis coordinates are an explicit reference coordinate system, not assumed equal
to donor coordinates. Reference occurrence must belong to the declared group.
No chromosome suffix parsing, suffix equivalence or candidate selection uses the
axis. Parent provisional source records translate by retaining
`component/start/end/group/reference_occurrence`, adding `reference_strand:"+"`,
and dropping `reference_path`/other provisional fields. All original source path
names and occurrence IDs survive in `calls.json`.

Any group appearing more than once **anywhere** on the axis is entirely excluded
from threading, with `unresolved-repeated-axis-group` at each appearance. Its
independent partition call remains. This avoids double-counted emissions or
inconsistent bundle genotypes without pretending to solve nonlocal inference.
Groups absent from the axis remain in calls and `unscaffolded_groups`.

## Quantitative score

All occurrences in a group with the same PanSN first two fields (`sample#haplotype`)
form ONE source bundle. The remaining contig spelling is never truncated. Names
without three nonempty PanSN fields fail explicitly (no guessed identity map).
Both `A#0#copy1` and `A#0#copy2` contribute to the same haploid bundle. Duplicate
BED IDs and overlapping representations remain provenance, but feature
multiplicity counts distinct `(source,start,end)` physical feature locations,
not rows. Multiple physical copies on one or several paths retain multiplicity.

Only catalog features with `owning_group=Some(group)` become factors. Each is
queried/scored once, including count zero. Excluded shared/nonlocal/boundary
features remain in the unchanged catalog and sample index, with reason totals
in calls; their sample counts are not re-queried or materialized by this mode.

For the full physical two-anchor span `b=end-start=gap+k` and aggregate read-length
histogram `(L,n_L)`,

```
e(b) = sum_L n_L * max(L-b+1,0) / sum_L n_L * L
lambda_hf = haploid_depth * e(b) * physical_multiplicity_hf + background
score_h = sum_f [lambda_hf - count_f * ln(lambda_hf)]
```

The omitted log-factorial is candidate independent. Implementation computes one
background-only baseline per group and sparse nonzero-multiplicity corrections:

```
delta_hf = d*e*m - count_f * ln1p(d*e*m/background)
score_h = background_only_score + sum_nonzero delta_hf
```

There is no dense candidate-by-all-factors scan, cosine, presence-only scoring,
heuristic candidate pruning, depth inferred from truth, or deletion from silence.
Both positive magnitude and absent factors affect scores. No positive counts or
no local factors gives no-call, with no best-bundle genotype assignment.

Default outputs retain all bundle scores, physical source-occurrence provenance,
aggregate modeled/positive feature and physical-location counts, expected/observed
totals, best-bundle ties, score gaps and mean Poisson deviance against the saturated
per-factor model. They record verified sample/catalog **payload** checksums, real
panel sidecar fingerprints, count policy, parameters and read-length histogram.
Sample checksum is FNV-1a64 over its bincode payload; catalog checksum is FNV-1a64
over its compact JSON payload, exactly the existing artifact checksum contracts.

Per-feature factors and bundle multiplicity maps are deliberately omitted by
default to avoid duplicating a potentially multi-GB reusable catalog in calls.
`--include-feature-details` is an optional, potentially huge debug output that
adds sparse multiplicities and factor IDs/counts/spans/exposures. Computation and
scores are identical either way. `feature_details_included:false` means omitted
maps are **unavailable, not zero**; deserialization defaults missing derived maps
to empty and consumers must honor that explicit flag. Original sample/catalog
artifacts are unchanged and suffice to reproduce the detail. Absolute tie
tolerance is 1e-8.
`poor-fit` means best mean deviance exceeds `--max-mean-deviance` (default 10);
this is a declared descriptive cutoff, **not a calibrated statistical test**.
Depth, positive background and cutoff must be finite and positive. Invalid or
nonfinite scores fail rather than producing JSON null confidence values.

Read-span exposure is an idealized working approximation. Counts are selected,
overlapping MEM-substring occurrences, not independent reads/molecules; feature
correlation, error background, MEM selection and orientation recall need later
calibration. Scores and score gaps are not posterior probabilities.

## Source-path DP and outputs

Each uniquely placed callable group supplies all physical source-occurrence
states, each linking its bundle and its complete source interval. Duplicate BED
rows with identical physical state/orientation merge only as state provenance.
Same-haplotype multicopy bundles are not split into mutually exclusive genotypes:
a thread state selects a source placement, while its emission is the full bundle
score accounting for ALL copies. Count calls exist independently of the axis;
regularized thread selection may differ from a local minimum, and records which
bundle was selected.

Orientation comes from reference identity, explicit axis input, both catalog
strands relative to the reference, or consistently agreeing shared signed two-node
contexts. Palindromic-only, unsupported and conflicting evidence stays unknown;
there is no majority vote. Explicit supplied orientation is authoritative.
Unknown states stay in DP candidates but cannot establish any continuation;
selected unknown states are unresolved and break emitted blocks.

For a chromosome/run of supported unique-axis groups:

```
D_i(s) = emission_i(s) + min_t [D_(i-1)(t) + transition(t,s)]
transition(t,s) = 0 iff exact source path, same known strand, and both source
                  endpoints strictly increase (+) or strictly decrease (-)
                = switch_penalty otherwise
```

DP resets at component boundaries and no-call/poor-fit/repeated-axis intervals.
It uses row-minimum-normalized emissions internally for numeric stability,
restores the constants in segment scores, and performs exact forward/backward
min-sum inference without pruning. All optimal states are reported. They are
**marginal alternatives, not freely combinable paths**; ambiguous intervals are
unresolved and break blocks instead of reporting an arbitrary certain phase.
A smoothed choice whose individual mean deviance exceeds the cutoff also remains
unresolved. Transitions cost a finite nonnegative `--switch-penalty` (default 10).

`threads.json` includes ordered axis intervals/states, optimal alternatives,
segment scores/ties, uniquely resolved chromosome blocks, explicit source and
reference gaps/overlaps, source-path switches and orientation/nonmonotonic block
breaks. Switch bounds span the two flanking partition intervals, not exact
partition edges. There is no cross-chromosome switch, backward-jump continuation,
fixed-N filler, source sequence assembly, or inversion call from smoothing.
Coverage is axis union bp and resolved union bp, not sum of candidate intervals.

Downstream evaluation separately reports best-single-bundle truth containment,
source-interval compatibility, local statuses and coarse source-path switch
recovery for adjacent single-occurrence truth groups. Repeated-axis groups are
excluded from switch assessment; unresolved/accessory results stay visible.
These are experimental source-panel diagnostics, not validated locus biology.

## Checks

```sh
# Same <=4-job native environment above; serialize native tests.
cargo test --offline --release -j4 --lib genome_inference -- --test-threads=1
cargo test --offline --release -j4 --lib sample_mem_bwt -- --test-threads=1
cargo test --offline --release -j4 --test test_genome_inference -- --test-threads=1
cargo test --offline --release -j4 --test test_genome_inference_quantitative -- --test-threads=1 --nocapture
```

The quantitative executable fixture builds a tiny panel and reads, then invokes
actual `build-sample`, `build-catalog`, and `genotype` commands. It exercises a
known A→B switch, reverse decreasing donor coordinates with preserved overlap,
three reference chromosomes, tied paths, no reads, a repeated-axis group,
unscaffolded same-haplotype copies and duplicate/overlapping BED representations.
Changed or absent truth leaves frozen calls/threads byte-identical. Failures
cover unsupported ploidy, nonpositive/nonfinite parameters, bad axis, missing
truth after freeze, and existing directories. Unit tests compare sparse versus
dense Poisson equations, show absolute counts change copy-state winners, audit
physical multiplicity despite duplicate rows, reject majority orientation votes,
and compare DP scores and every optimal marginal to exhaustive path enumeration.

Final validation passed **433 library tests, 85 CLI unit tests, the existing
bootstrap integration fixture and the new quantitative integration fixture**,
serialized under the declared four-job environment. A final quantitative fixture
run used 484 reads/242,000 bases, 242 distinct MEMs, nine groups, 1,666 retained
factors and 28 excluded features. Compact calls were 14,513 bytes; quantitative
call+thread+evaluation took 0.0084 s after fixture inputs existed (not a full-panel
performance claim). Five of nine axis intervals resolved (15,000/24,000 axis bp),
with four explicitly unresolved intervals and one of one known coarse source-path
switches recovered. Reverse continuation retained a -500 bp source overlap.

Scaling remains the existing in-memory JSON catalog/checksum load (large temporary
serialization allocations), sparse bundle-factor storage, and exact quadratic
adjacent-state DP. No full-panel measurement or serialization redesign is claimed.
No dependencies, cache formats, native policy or extraction/count policy changed.
