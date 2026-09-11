# Partition-major observation replay (experimental v2)

Production commands now implement the observation repair described in
[observation-registry-repair.md](observation-registry-repair.md). Native small
fixtures, independent review and the 590-test parent release gate pass. Real-panel
resource/sequence controls are still required. This is **error-free source replay**, not a calibrated noisy-read
model, posterior, biological absence call or complete assembly method.

## Commands

```sh
impg genome-infer build-partition-observations \
  --panel PREFIX --catalog catalog.json --sources panel.agc \
  --read-lengths 150 --out-dir observations

# Optional bounded FEATURE scope; every partition and source is still scanned.
# Repeat --feature-group to select multiple original owning groups.
impg genome-infer build-partition-observations \
  --panel PREFIX --catalog catalog.json --sources panel.agc \
  --read-lengths 150 --feature-group partition986 --out-dir observations-986

impg genome-infer genotype-partitions \
  --panel PREFIX --observations observations --sample sample.membwt \
  --haploid-depth 10 --background 0.1 --max-mean-deviance 10 \
  --allow-unvalidated-catalog --out-dir calls \
  --axis axis.json --switch-penalty 10

# Reuse frozen full-scope calls without the legacy catalog or sample index.
impg genome-infer thread-partitions \
  --panel PREFIX --observations observations --calls calls/calls.json \
  --axis axis.json --switch-penalty 10 --out-dir rethread

impg genome-infer reconstruct \
  --calls rethread/calls.json --threads rethread/threads.json \
  --panel-names PREFIX.names --sources panel.agc --out-dir reconstruction
```

`--sources` accepts one or more indexed FASTA/AGC inputs with exact full panel
source names, one unambiguous binding per path and matching lengths. FASTA index
creation uses the existing native FASTA machinery. Reconstruction must use the
same source-file contents as compilation. `--read-lengths` accepts comma-separated
or space-separated values. Every sample histogram length must be compiled.

`thread-partitions` writes both `calls.json` and `threads.json`; reconstruct using
that matched pair, especially after changing the axis. The derived calls preserve
all original non-provenance bytes (scores, statuses, rankings, rates and sample/
registry links), and record a checked raw-content link to the original calls.
Only axis-specific orientation and derivation provenance changes. Calls made
without `--axis` do no orientation work; late rethreading derives the requested
provider from the existing checked incidence stream.

No truth option exists on the compiler/new caller. Depth, background and deviance
cutoff remain independent inputs. The original `build-sample`, `genotype`, `call`
and reconstruction of legacy v1 artifacts retain their original policy/model.

Partial scopes retain **all originally owned features** of requested groups,
including zeros. Unknown/empty selections fail. Partial outputs cannot be threaded
or reconstructed. Unselected groups without retained factors are labelled
`no-call-outside-assessed-feature-scope`, not absence. An explicitly selected
group whose factors were all excluded remains assessed and reports
`no-call-no-local-features`. Omit `--axis` when scoring
a partial scope. The factor ledger still reports newly reached other-group,
boundary and context-dependent support found anywhere in the panel.

## Model and artifacts

New calls: `partition-observation-poisson-v2`; threads:
`partition-observation-source-min-sum-v2`. New calls carry typed `observations`
provenance and leave `catalog_payload_checksum` null. Legacy catalog checksum
headers are recorded as **unverified headers**; the streamed raw file has a
separate newly computed FNV-1a content checksum. FNV checks detect accidental
corruption/identity changes, not adversarial tampering.

The compiler writes:

- `inputs.json`: immutable slim original source/group/occurrence IDs, input content
  fingerprints, supported lengths, scope and source-bound compiler identity;
- `definitions.jsonl`: compact feature definitions with original IDs/owners;
- `ownership.json`: normalized fixed physical ownership, requiring complete
  source coverage and rejecting conflicting owners;
- `shards/`: atomically completed partition/core shards and checkpoint checksums;
- `incidences.jsonl`: deterministic sorted physical records, keyed by feature,
  source-forward start/end. Each has signed orientation evidence, full-span owner
  (or boundary), context-nonlocal/source-terminal flags, and per-length/per-input-orientation integer
  surviving-substring totals. Raw views/duplicate BED descriptions never add copies;
- `metadata.json` and authoritative succeeded `manifest.json`: checked content
  linkage, model/capability, lengths, scope and execution statistics.

The consumer streams definitions and incidences, queries each selected sample
feature exactly once (also excluded/zero features), and writes `factor-ledger.jsonl`
with observed counts, original/recomputed owners, exclusion reasons, physical
location counts and zero-profile counts. Calls retain summaries, not misleading
single scalar exposures or huge per-occurrence multiplicity maps. Strict signed
orientation masks are derived only when an axis is requested. One reference is
assessed per non-repeated axis group, using direct shared-feature OR semantics;
there is no all-member-pairs matrix, majority vote or transitive inference.
`orientation_references` records every assessed reference, including zero-support
references. An uncomputed or different-axis provider is rejected, not interpreted
as assessed zero/unknown. Repeated-axis groups remain excluded from DP.

The derivation uses a (group, source) max-end interval index for original-description
containment, then one byte of accumulated support per original occurrence. At most
one reference/member mask is retained per requested group member. Sample scoring
has no occurrence-containment scan or orientation aggregation.

For each physical copy, integer totals `q[L,orientation]` are orientation-averaged,
weighted by sample read-length counts, and divided by histogram total bases.
Selected bundle copy rates are **summed before the single feature logarithm**:
`mu = background + haploid_depth * sum(copy_rates)`.
Overlapping surviving MEMs may contribute twice at one physical coordinate; a
read spanning distinct copies assigns each substring only to its own synthetic
coordinate. Actual read coordinates and names remain transient and are not stored.

A factor is locally usable only when its complete physical support has one owner,
no crossing, and the **untruncated support envelope of every spanning read length**,
including zero-count contexts and lengths longer than the original source, remains
inside that fixed owning piece. Negative lower bounds are retained when checking
validity. Unvalidated source-terminal extension is reported as
`source-terminal-context` as well as `context-nonlocal`, because placing a terminal
donor piece at a novel internal junction could introduce additional spanning reads.
Conditional integer profiles still clip valid starts to original source ends;
there is no out-of-source fetching, padding or candidate expansion. Long-span features with no supported starts remain represented with zero signal
(background prediction). Cross-partition/shared/context-nonlocal factors are
excluded once, not silently assigned zero signal or solved by expanding candidates.
Conditional source profiles remain available, but no complete-template ranking
command or general nonlocal/junction solver is included.

New-model source eligibility requires its own positive signal support and score
strictly better than background by the existing epsilon. Equality within epsilon
is `ambiguous-background-equivalent`; strict background preference is
`no-call-background-preferred`. Both have empty best sets and explicitly block DP,
resetting segments. `score_gap` remains the best-versus-second source score gap
across **all** bundles, not just eligible states. DP retains all eligible alternatives, not only local
minima, and cannot reintroduce an ineligible source. Explicit axis/BED orientation
precedence, strict signed-context conflicts, repeated-axis exclusion, chromosome
resets and existing conservative reconstruction semantics remain in force.

## Exact aggregation and bounds

The existing orientation-specific native extractor supplies two separate traces.
For every synthetic spanning-read range, exact anchor enter/leave events partition
starts into runs with the same **both-view** node/relative-coordinate context up to
uniform translation. Within a run richer-view selection/ties, node/gap GBWT MEM
queries and coordinate/content subwalk pruning are invariant; the coordinate-tagged
integer result is multiplied by the run length. There is no sampling, target-only
cache, synthetic FASTA/BWT construction in production, or per-start native replay.

Each partition-major core is capped at 1 MiB. Its DNA and both raw views are
extracted **once**, with left/right read-only halos of
`max(largest_supported_read_length, K - 1)`, clipped to the actual source ends.
Discovery and endpoint responsibility include only anchor starts inside the core;
halo anchors never add ownership. Main-shard profiles reuse these raw views.
Before constructing events, both streams are binary-sliced and rebased to each
hit's bounded spanning-read envelope; neither cloning nor event construction scans
the entire core per incidence. `total()` additionally bounds event-anchor visits
by the specific read length's envelope. Per-event MEM selection/counting and the
untruncated local-validity rules are unchanged.

Endpoint hits without a live core use bounded fallback fetches, only when a
supported spanning-read start exists. First/last raw anchors of nonempty cores
are retained separately for each actual source/view;
source-contiguous endpoint stitching recovers adjacent pairs across long empty
cores without loading chromosome DNA or mixing views. Source DNA is checked against
every stored dictionary anchor encountered. Exposure contexts are at most twice
the largest supported read length. The explicit limits are 32 distinct lengths,
1,048,576 bp per length, native i32 source coordinates, 64 KiB serialized rows and
32-way merge fan-in. Arithmetic for accumulated substring contributions is checked.
Exceeding a bound fails rather than truncating or approximating.

Memory retains the loaded native panel, compact fixed-width feature lookup, slim
ownership/occurrence metadata, endpoint records and a bounded core/profile context;
not the enormous legacy location or occurrence-multiplicity maps. Axis-specific
orientation storage is linear in original occurrences, not same-group pair counts. Source-count/core
metadata and checkpoints scale with panel size. Empty core shards are not written.
Statistics report core/endpoint counts, exact profile events, checkpoint-shard bytes
and merged incidence bytes. `core_context_fetches` and `core_context_bp_fetched`
record main extraction work. `profile_execution` records reused profiles,
`fallback_context_fetches`/`fallback_context_bp_fetched`, no-start profiles, and
maximum per-hit context bases/anchors. Each fetched context supplies two raw-view
calls. Remaining work includes exact per-hit/per-length event/MEM replay and one
bounded fetch per supported endpoint incidence; heavily fragmented ownership can
increase endpoint fallback work. Full-panel runtime, peak memory and scratch capacity
remain to be measured; this is not an all-panel performance claim.

**Fresh-only, no resume:** every output directory must be new, including empty
ones. There is no `--resume`. Failed/running outputs preserve diagnostic shards but
cannot be consumed or reused as completed output. Checked metadata/model/compiler,
panel, scope, length or content mismatches fail closed. The source-bound compiler
identity intentionally invalidates artifacts after relevant observer changes.

## Native tests

`cargo test --release --test test_partition_observations -- --test-threads=1`
exercises the real CLI compiler → sample → calls → rethread → reconstruction. The
fixture recovers 87 omitted physical incidences, verifies forward/RC prediction
ties, duplicate-description conservation, boundary/context exclusions, partial
scope and integrity failures, and spells the exact 2,400 bp source. It also derives
a different axis reference from initially unthreaded calls, verifies unchanged
non-provenance data and original-file linkage, and reconstructs the exact source
fragments under that changed axis.

`cargo test --release --lib genome_inference -- --test-threads=1` includes optimized
versus exhaustive native replay (512 feature/location/length totals, 12,426 raw crop
comparisons, 111 public singleton indexes, 28 multi-copy reads), a separate native
overlapping-MEM double-count fixture, long zero-exposure endpoint recovery across
an empty core, histogram/copy-weighted consumer checks, terminal and longer-than-source validity, background eligibility/
reset checks, and legacy calling/threading/reconstruction regressions.

The raw-window oracle additionally covers K−1/K/K+1 (62/63/64 bp), Ns, and
anchors ending exactly at the read end: 70 nonempty exact-K windows and 324
end-aligned anchors. Native selection hashes K−1 bases, but exported anchor
validity and profile inclusion use the full K-base guard `p + K <= read_end`.
The shorter hash-selection footprint is not used for inclusion events.

Additional tests exercise a 2,048-member group with only 2,047 retained direct
reference/member masks (not a quadratic matrix), zero-support assessed references,
repeated-axis exclusion, indexed containment versus exhaustive descriptions,
byte-preserving derived calls, BG-equivalent reconstruction validation, reverse-index
epsilon ties (while rejecting duplicate best indices), 17 incompatible frozen-call
header/provenance cases, and an explicitly selected partial group whose factors
were all excluded. The near-tie and rethreading defects were reproduced before
fixing them; the focused fixes passed independent review and the full release gate.

Core-reuse tests compare 422 physical profiles against fresh bounded native
extraction: 418 reuse 20 core contexts, four have no supported starts, and all
profiles/events match exactly. Their largest per-hit envelope is 376 bp/28 anchors.
A separate compiler fixture verifies a positive cross-owner endpoint profile
against its shared-core counterpart (one 183 bp fallback fetch). The CLI fixture
reuses five core extractions (6,880 bp total) for 294 profiles, with no per-hit
fallback extraction; its 6,166 exact events and 87 recovered incidences are unchanged.
These are small-fixture execution counters, not full-panel performance claims.
