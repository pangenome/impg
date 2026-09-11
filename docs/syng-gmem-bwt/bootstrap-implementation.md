# Experimental sample MEM-BWT bootstrap

This is **not completed whole-genome genotype/mosaic inference**. The approved
first executable stage is `impg genome-infer`: reusable sample counts, an
occurrence-preserving ownership catalog, and **presence-compatibility source
occurrence diagnostics**. It does not certify haploid alleles, biological
homology, ploidy, dosage, posterior probabilities, phasing or assembled sequence.
It does not reuse the old per-read cosine scorer or stitcher.

| Implemented in this stage | Required later, not implemented |
|---|---|
| Reads -> weighted collection FM/BWT; arbitrary contextual count queries | Measured whole-panel scalability / compressed-rank optimization |
| BED3 ownership import or explicit JSON groups/scaffold | Validated homologous cores and spanning candidates; source-sequence validation |
| Every input interval/copy/name/orientation retained; coordinate neighbor links | Joint locus/copy configurations, chromosome-safe mosaics/sequence emission |
| Full stored-panel scan for global two-node feature ownership | Quantitative exposure/count factors, nonlocal repeats and joint boundary factors |
| Explicit unique-compatible occurrence, ambiguity, no-call and conflict | Calibrated haploid genotyping, mixtures, dosage and recombination |
| Separate downstream source-occurrence truth compatibility evaluation | Biological genotype accuracy; strict rebuilt hold-out; noisy real-read validation |

The completed yeast BEDs remain `catalog_accepted=false`. Calls require
`--allow-unvalidated-catalog`, including explicit JSON groups; there is no
backdoor that marks a catalog accepted. Ploidy other than 1 is rejected. This
acknowledgement does not make individual BED intervals mutually exclusive alleles.
Two nonidentical copies of the **same haploid assembly** may both contribute
features to one group: neither individual interval then explains all positives.
The resulting conflict means inadequate locus/evidence modeling, **not evidence
of a mixed sample**. A singleton compatible interval is not a certified genotype.

## Commands and safe artifacts

```
impg genome-infer build-sample --panel PANEL.syng --reads reads.fq.gz --out-dir NEW_SAMPLE
impg genome-infer build-catalog --panel PANEL.syng --bed-dir COMPLETED_BEDS --out-dir NEW_CATALOG
impg genome-infer call --panel PANEL.syng --sample NEW_SAMPLE/sample.membwt \
  --catalog NEW_CATALOG/catalog.json --allow-unvalidated-catalog --out-dir NEW_CALLS
impg genome-infer run --panel PANEL.syng --reads reads.fq.gz --bed-dir COMPLETED_BEDS \
  --allow-unvalidated-catalog --out-dir NEW_RUN
```

`--groups groups.json` replaces `--bed-dir`. `--reads` accepts one or more
FASTA/strict four-line FASTQ files; compression is detected by existing niffler.
Read IDs and quality strings are transient, not written. IUPAC sequence symbols
are accepted; the existing matcher only seeds fully ACGT windows. Empty files,
empty records, invalid symbols, truncated FASTQ and invalid quality lengths fail.
FASTA can be multiline. No quality filtering or frequency/sketch policy changes.

Every output directory must be **new**, and its parent must already exist.
`create_dir` reserves it atomically; even existing empty directories and symlinks
are refused. JSON/binary files publish via a synced `.incomplete` temporary file.
`manifest.json` states `running`, `succeeded` or `failed`; only `succeeded` is an
accepted stage result. On handled failures calls/evaluation are removed; partial
index/catalog objects may remain for diagnosis under FAILED status. Abrupt
termination leaves running/missing status or incomplete files, never accepted
results. Do not reuse artifacts from a failed stage as accepted evidence. No
input files, indexes or existing output directories are overwritten.

Outputs: `sample.membwt`, `catalog.json`, `calls.json`, optionally
`evaluation.json`, and `manifest.json` (only applicable files per subcommand).
The manifest records argv, elapsed time, format/model and measured byte sizes;
`run` also records sample/catalog build-and-save and call/evaluation durations.
Sample serialization is deterministic with respect to read order; calls are
independent of the truth file and contain no calibrated confidence field.

## Count contract, version 1

A MEM is an alternating vector `(node, gap, node, ..., node)` with oriented
nonzero signed i32 nodes (excluding i32::MIN), positive u32 distances between
anchor starts, and no initial gap. Tokens are u64:

* 0: unique terminal; 1: separator between independent MEM records.
* node: `2 + 2 * zigzag_i32(node)` (even).
* internal distance: `2 * gap + 1` (odd, at least 3).

Queries must start/end on nodes. An interior subwalk therefore does not inherit
the preceding anchor's gap. Separators cannot appear in valid queries. The
complete two-anchor genomic span includes the final anchor's full syncmer
length, not only the start positions. Graph walk/spacing exactness does **not**
prove nucleotide-exact sequence across unrepresented interior gaps.

Construction uses `matched_syncmers_in_sequence` and `gbwt_mems_for_walk`. The
existing matcher already tries both nucleotide orientations and selects the
higher-anchor trace (input wins ties). This stage invokes that unchanged API on
both input orientations; it does not union arbitrary anchors or change sketch
selection. Each returned MEM is normalized into input-forward coordinates.
Identical **coordinate AND node** vectors are deduplicated within a read. A
record is removed only if it is an exact contiguous coordinate-and-node subwalk
of a longer record; envelope containment alone is insufficient. Overlapping,
noncontained or distinct-content records survive separately. The underlying
GBWT API prunes contained MEMs within its single trace. Neither construction
step joins separate MEMs or retains read/mate linkage.

Remaining records are RC-canonicalized (reverse node order, negate nodes,
reverse gaps), then identical strings aggregate multiplicity. Define oriented
`C_o(f) = sum_t w_t * occ(f,M_t)` on the **stored** strings. The sample API reports
RC-orbit counts `C(f)=C_o(f)+C_o(rc(f))`, except palindromes counted once. RC
canonicalization preserves this orbit count, not the original input-oriented
count. Repeated substring positions count repeatedly. Separate overlapping MEMs
can multiply evidence from a read: counts are **MEM-substring occurrences, not
reads or molecules**. No claims about effective independence/exposure follow.

The transient construction tape has one copy of each distinct MEM plus a
separator, with its multiplicity on each suffix position. Prefix-doubling builds
a suffix array. Only BWT tokens and suffix-order cumulative weights survive;
SA, tape, record-to-read mapping, read names and placements do not. Sparse rank
position lists are rebuilt from the BWT when loading. They locate **BWT rows,
not original reads/positions**. Backward-search weighted intervals give exact
substring counts; this is not a table of every possible substring.

`sample.membwt`: magic `IMPGMEM1`, little-endian u64 payload length and checksum,
then bincode-2 standard encoding of version, panel identity, exact count-policy
string, aggregate stats/read-length histogram, BWT and cumulative weights.
Checksum/length/version/policy/structural checks reject damaged or incompatible
objects. `catalog.json` is a versioned JSON envelope containing a checksum and
typed payload. Artifact checksums and panel fingerprints are explicitly stable
**FNV-1a64**, for accidental corruption/identity detection, not cryptographic
integrity against adversarial modification. Do not treat these as SHA-256.

Panel identity fingerprints the actual contents/lengths of `.1gbwt`, `.1khash`,
`names`, `meta`, `spos`, `pstep`, resolving primary/legacy sidecars as SyngIndex
does. It includes the native dictionary and graph, not merely k/s/seed. Metadata
v2 is mandatory. Different panels with the same syncmer parameters are rejected.
Resampling changes the identity deliberately; rebuild corresponding artifacts.
No AGC linkage or source sequence validation is implied by these checksums.

## Ownership catalog and restricted global factors

BED import accepts exactly three tab-separated fields per nonempty line, one
file/group, sorted by filename, retaining **every row**, including duplicates.
Names are exact full panel names; no PanSN splitting or chromosome assumptions.
Bounds are zero-based half-open. Empty groups, unknown paths, invalid bounds,
duplicate group IDs, malformed JSON and bad orientations are errors. Bounds,
fully contained stored anchor counts and availability of two-anchor contexts are
validated; homology, spanning, source spelling and copy structure are not.
Unsupported intervals remain in the catalog with explicit status.

Versioned explicit input (optional strand and scaffold fields can be null):

```json
{"version":1,"groups":[{"id":"group1",
 "scaffold":{"component":"verified-chrI","start":0,"end":10000},
 "occurrences":[
  {"path":"A#0#chrI","start":100,"end":10100,"strand":"+"},
  {"path":"B#0#other_contig","start":45000,"end":55200,"strand":"-"}
 ]}]}
```

Scaffold coordinates are an **explicit caller-supplied evaluation ordering**, not
numeric slicing of alternative genomes or homology evidence. Without a scaffold,
groups are sorted by ID and have no invented chromosome coordinate. With one,
diagnostics sort by component/start/end; there are no transitions or joins.
Source-neighbor links are derived independently on each exact source path.
Neighbor coordinate buckets preserve all duplicate members' links. Forward
traversal increases and reverse traversal decreases coordinates; gaps/overlaps
are explicit signed distances. Unknown/mixed orientations and nonmonotonic
nested intervals are not certified as oriented continuation. Links are neighbor
diagnostics, not exhaustive recombination compatibility or assembly instructions;
there is no trimming, fixed-N gap filling or FASTA output.

The fixed feature universe is all adjacent two-anchor contexts wholly contained
in at least one candidate interval. Candidate multiplicities preserve tandem
repeats. A second scan over **every full stored panel path**, including paths not
in the catalog, retains all occurrences of these features and their containing
source occurrence IDs. Features with any outside-catalog or boundary-crossing
occurrence are excluded from this **baseline caller**, as are features contained
in multiple groups. Reasons and global counts remain available; neither sample
index data nor catalog occurrences are masked/dropped. Stored-path RC sketch
recall remains imperfect: "all stored occurrences" is not "all biological or
nucleotide occurrences". Unrepresented off-panel repeats also remain a risk.

A retained feature has exactly one global owning group and one evidence factor.
For each group the positive retained features (`C(f)>0`) constrain compatible
source occurrences to those containing **all** such features. The magnitude,
absence, exposure, depth and error probability are ignored. No positive feature
means no call; no compatible interval means conflict; multiple intervals mean
ambiguity. One erroneous positive can create a conflict (tested). Shared/global
counts are not duplicated as local full-strength emissions. Exclusion is a
conservative bootstrap restriction, not a new repeat-filter policy for the sample
index/future quantitative pipeline. Boundary features are explicitly excluded;
the `AB|ab` vs `Ab|aB` counterexample is not modeled by local dosage terms.

## Truth: post-call source occurrence diagnostics only

`--truth truth.json` is opened only after frozen calls are written. Example:

```json
{"version":1,"coordinate_system":"explicit-fixture-reference",
 "groups":[{"group":"group1","occurrences":[
  {"path":"A#0#chrI","start":100,"end":10100,"strand":"+"}
 ]}]}
```

Truth intervals must use known source paths and valid coordinates. Missing
candidate occurrences can be evaluated but unknown/off-panel source names are
currently rejected; strict hold-out evaluation is deferred. Orientation null is
compared exactly, not assumed forward. Metrics separate truth-set containment in
the catalog, truth-set containment in compatible occurrences, exact unique
source-occurrence recovery, ambiguity, no-call and conflict. They are not
biological genotype accuracy. Duplicate source rows remain separate alternatives.

Base-pair metrics require a nonempty coordinate-system label and explicit
scaffold for **every** evaluated group; otherwise they are null. Their denominator
is the component-wise coordinate union of evaluated scaffold intervals, never
the sum of panel alternatives; compatible bp uses the same union convention.
Cross-component totals use checked accumulation: overflow is an evaluation error,
with failed status and removal of calls/evaluation rather than wrapped bp metrics.
This is evaluated-scaffold coverage, not whole-sample coverage unless its
completeness was independently established. `genotype_callable_bp` is always
null in this unvalidated bootstrap, so unsupported/excluded bases are never
reported as inferred.

## Reproducible tiny runner and evidence

Use the existing native build environment (no installation required):

```sh
export CFLAGS=-I/home/erikg/.cache/impg/native-build-deps/root/usr/include
export LIBRARY_PATH=/home/erikg/.cache/impg/native-build-deps/root/usr/lib/x86_64-linux-gnu
export LD_LIBRARY_PATH="$LIBRARY_PATH"
export CARGO_TARGET_DIR=/home/erikg/impg/target
export CARGO_BUILD_JOBS=4 RAYON_NUM_THREADS=4
cargo build --offline --release --bin impg --bin gfaffix -j4
cargo test --offline --release -j4 --lib sample_mem_bwt -- --test-threads=1
cargo test --offline --release -j4 --lib genome_inference -- --test-threads=1
cargo test --offline --release -j4 --test test_genome_inference -- --test-threads=1 --nocapture
```

The integration test is the one-command fixture runner: it constructs a tiny
native panel, deterministic error-free reads, explicit multiple chromosomes,
unique signal, a deliberate identical-profile tie, no reads at a locus, two
nonidentical haploid copies, tandem repeats, off-catalog/shared/boundary contexts,
and invokes **all four actual CLI stages**. It verifies repeated-row retention,
reverse links, gaps/overlaps, standalone artifact reuse, read-order/orientation
invariance, truth independence, single-positive error sensitivity, index and
catalog corruption/incompatibility, malformed reads/bounds, acknowledgements,
ploidy rejection and failed/unsafe output behavior. Pure count tests exhaustively
compare weighted oriented/orbit counts with literal substring counting, including
spacing, overlaps, boundaries and palindromes. No real panel is opened by tests.

A measured fixture run (14 paths / 42 kb panel; 312 x 500 bp reads) wrote a
24,048-byte sample index: 156 distinct MEMs, 312 weighted records, 4,951 collection
symbols. The catalog was 604,653 bytes, with 14 source occurrences in ten groups
and 780 retained global factors. One measured `run` took 0.0398 s, including
0.0226 s sample build/save, 0.00677 s catalog build/save and 0.00154 s calls/evaluation.
The complete integration test (all success/failure/reuse scenarios, not just index
construction) took 3.08 s with 16,072 KiB maximum RSS under `/usr/bin/time -v`.
These are fixture measurements, not promised benchmark speed or sample-only RSS.
There is no assertion of whole-panel compactness from the BWT name or this fixture.

Independent review accepted this restricted bootstrap with one demonstrated
scaffold-total overflow finding, subsequently fixed by the parent. The parent
reran 428 library, 85 CLI unit and 40 integration tests successfully, including
cross-component overflow and actual CLI failed-output cleanup, plus existing
empty-output refusal. This does not establish biological catalog acceptance or
whole-panel performance. Additional catalog/extraction oracles and publication
fault-injection coverage remain review recommendations.

Complexity: sample suffix construction O(N log²N) time/O(N) temporary storage,
where N is the distinct-record tape size. Sparse rank uses O(N) positions and
O(alphabet) vectors, not a dense alphabet-by-text table; queries cost
O(pattern tokens * log N). Catalog storage is O(features + stored feature
occurrences + source intervals); whole-panel adjacent-anchor scans and JSON
serialization can be large. All occurrences and features currently reside in
memory; duplicate coordinate buckets can have many cross-product neighbor links.
This is a correct small/bootstrap implementation, **not a demonstrated 3.3 Gb
panel memory budget**. Parent must measure RSS/time before scale acceptance.

## Parent-only whole-genome error-free diagnostic invocation

After independent review, use a known in-panel haploid's **complete chromosome
set** and deterministic error-free reads prepared by the parent, not the old
four-window pilot or truth-seeded candidates. The input reads must include all
chromosomes intended for evaluation; completed BED discovery is not homology
certification. Do not run this full-panel command from the implementation lane:

```sh
export CARGO_BUILD_JOBS=4 RAYON_NUM_THREADS=4
export LIBRARY_PATH=/home/erikg/.cache/impg/native-build-deps/root/usr/lib/x86_64-linux-gnu
export LD_LIBRARY_PATH="$LIBRARY_PATH"
IMPG=/path/to/parent-accepted/impg
READS=/path/to/parent-fixed-error-free-whole-haploid.fastq.gz
OUT=/path/to/new-whole-genome-bootstrap
/usr/bin/time -v "$IMPG" genome-infer run \
  --panel /home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng \
  --bed-dir /home/erikg/yeast/partition-pos64-w10k-d1k-to-completion/results \
  --reads "$READS" --allow-unvalidated-catalog --out-dir "$OUT"
```

The `results` directory is recorded in the completed-run `report.json` argv.
No AGC or graph rendering is required.
To obtain ordered/evaluable coordinates, independently prepare the versioned
group scaffold, substitute `--groups`, and add a downstream `--truth` manifest.
Do not infer a scaffold from contig name suffixes or sum alternative intervals.
For subsequent samples reuse `catalog.json`; for rescoring reuse `sample.membwt`.
The controlled whole-genome run is a bootstrap diagnostic, not acceptance of the
requested genome genotype/mosaic pipeline. Real SK1/Y12 noisy-read claims,
mixtures and graph repairs are expressly outside this stage.

Next interfaces: validated chunk/core/flank candidates with separate spelling,
homology/spanning and copy-configuration evidence; quantitative globally owned
feature exposures/background with calibrated count factors; joint boundary and
nonlocal states; explicit chromosome coordinate/continuation contracts for
mosaics. Those must build on the same occurrence identities and weighted count
primitive rather than converting diagnostic singleton rows into genotypes.
