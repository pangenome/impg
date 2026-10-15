# Bounded partition scheduling diagnostics

Status: diagnostic groundwork, **not a scheduling repair**. No covered-window
skip or query-context change is implemented. Real-panel performance and discovery
acceptance remain pending. The supervisor approved this instrumentation-only
result after the directional counterexample below.

## Contracts and reproducible work

`commands/partition.rs` selects sample/haplotype groups by missing bases, but
windows the full length of each selected path that still has missing regions.
Fully claimed paths are no longer in `missing_regions`; they are not newly added
merely because another path in their group remains missing. Already queued
windows are dispatched even if intervening queries claim their entire source.
`total` also windows a selected whole path; `longest` windows a missing interval.
Tail windows are merged as before.

The syng wrapper ignores traversal masks/depth and queries the original context.
The filtered wrapper chains raw anchored intervals and applies unchanged anchor
and query-span-fraction thresholds. Partition then retains the legitimate seed,
merges, extends boundaries, subtracts ownership masks, and updates missing regions.
These intervals are half-open, despite their `coitrees::Interval` container.
Seed retention is not fabricated homology.

Backend contracts differ: single PAF `Impg` DFS/BFS subtract visited/masked ranges
from the input before expanding it. `MultiImpg` also initializes visited ranges
from masks and inserts the seed through `SortedRanges::insert`, queuing only
returned unvisited fragments that meet its transitive-length rule. A fully
covered seed queues no subindex expansion. Neither backend is changed, and syng
equivalence must not be generalized to them.

The actual scheduler is exercised by `partition_scheduling_tests.rs`, using a
counting fake backend whose directional/context-dependent responses are explicit:

| Fixture | Calls before / after instrumentation | Zero-core dispatches | Observation |
|---|---:|---:|---|
| A covers queued B; B returns nothing | 2 / 2 | 1 | B produces no new source |
| Selected B has 10 of 1000 bp missing (sample/haplotype/total) | 11 / 11 | 9 | Nine zero-output calls; final 100 bp context discovers the target's remaining 5 bp |
| A covers B; covered B discovers C+D | 2 / 2 | 1 | B emits 200 new bp and groups C+D |
| Duplicate A seed and two distinct B copies, one reverse-oriented response | 5 / 5 | 3 | Exact A tuple revisited; both B occurrences and intervening residual survive |
| Covered B returns a backend error | 2 / 2 | 1 | Error propagates, not converted into a skip/singleton |

The extra direct residual-only backend probe in the context test is outside the
11 scheduler calls and returns no homolog. Removing B's discovery in the
asymmetric fixture leaves the exact source union unchanged but splits C and D
into separate singleton partitions. This is why zero core is not sufficient for
safe suppression. A separate minimum-size test verifies that a 5 bp boundary
residual is absorbed, not omitted, at minimum missing size 6 (one call).

Real syng integration fixtures check masked DFS/BFS against unmasked wrapper
results and compare byte-identical BEDs with profiling off/on for both raw and
filtered queries. The two-path 6 kb CLI fixture makes three calls in both modes,
unchanged by profiling. Existing native locks and serial execution are retained.
These fixtures do not resolve the separately documented stored-path RC sketch
asymmetry or demonstrate biological catalog/genotype acceptance.

## Opt-in controls and fields

Set `IMPG_PARTITION_PROFILE_MAX_QUERIES=N` (1..10000). An absent variable disables
profiling. Optional `IMPG_PARTITION_PROFILE_SKIP_QUERIES=S` (default 0) selects a
contiguous sample after S dispatched queries. Invalid values warn and disable
profiling. These are **diagnostic sampling bounds, not execution limits**: all
queries still run. Do not enable partition `--debug` or `SYNG_EMIT_PROFILE` for
this bounded measurement; those existing dumps have independent/unbounded output.

At most N coordinate-tuple keys and 3N+1 INFO rows are retained/emitted by this
profiler. There is no per-base, per-read, per-hit or cached-result storage.
Disabled mode does not scan missing intervals or time queries; wrapper overhead
is a thread-local enabled check and conditional metric bookkeeping. Sampled
missing-core accounting binary-searches existing sorted ranges.

Rows join by `dispatched` (one-based backend dispatch ordinal):

- `partition_profile_dispatch`: logged **before** lookup, including sequence ID,
  name, half-open start/end, `query_bp`, current `missing_core_bp`, cumulative
  `queued`, `dispatched`, and `visit_in_sample`. This identifies an unfinished
  expensive query if an external timeout interrupts the run.
- `partition_profile_query`: completed call's `query_us`, `postprocess_us`,
  `mask_us` (a subset of postprocess), `output_us`, `backend_overlaps` before seed
  retention/merging, and `emitted_new_source_bp` after ownership masking and final
  merge. The emitted count includes all source paths and boundary absorption,
  not just the seed. `error=true` means backend failure, not no homology; its
  overlap/emission counts are unavailable placeholders (zero), not valid results.
- `partition_profile_syng`: successful syng `raw_hits` (returned grouped lookup
  intervals, **not native occurrence-table visits**), `raw_anchors` (summed anchor
  vector lengths), `chains` before wrapper filters, `passing_chains` after both
  filters but before sequence-ID conversion, `lookup_us`, `chaining_us` (chaining,
  threshold filtering and conversion). `None` marks fields unavailable in raw
  mode. No syng row is emitted for other backends or failed lookups.
- `partition_profile_summary`: total queued/dispatched and sampled count, unique
  query tuples **within the sample**, skip and limit. Emitted on normal return or
  Rust error unwinding, not guaranteed on SIGTERM/SIGKILL. Counts after the sample
  are available only in this summary.

`query_us` includes the synchronous wrapper call and dispatch-log overhead.
`postprocess_us` includes seed retention, merge, boundary extension, ownership
mask/update and final merge. `output_us` is the remaining loop time (length
accounting, collection or per-partition writes, normal progress logging); it does
not measure final combined-file/GFA conversion after the scheduler. Timings are
wall-clock microseconds, include small profiling overhead, and are not additive
across the nested stages. This does not separately count native high-copy table
lookups, but raw anchors and lookup/chaining times expose their cost downstream.

Classify completed successful rows separately:

1. `missing_core_bp=0`, `emitted_new_source_bp>0`: covered-source productive work;
   do not count this as safely removable.
2. `emitted_new_source_bp=0`: zero-output query in this run, not proof it could be
   omitted without preserving error/discovery behavior.
3. `visit_in_sample>1`: exact repeated `(seq_id,start,end)` dispatch. It does not
   mean a shared node, repeat copy, or merely covered coordinate. With skip > 0,
   visit 1 means first **in sample**, not necessarily first in the whole run.
4. Compare total query/lookup/chaining/postprocess/mask/output times in those
   categories. High anchor counts with costly first-time productive queries
   would support genuine matching/chaining cost, rather than simple rescheduling,
   as the next target for investigation.

## Parent-only bounded real-data replay

The worker did **not** execute the following. This reuses the recorded 1200 s
manifest's algorithm parameters, but uses a shorter 120 s diagnostic timeout and
the instrumented working-tree binary. It does not overwrite the frozen index,
binary, or previous results. A timeout is expected to mean incomplete discovery,
not success. Parent owns approval/execution and output manifest/QC.

```bash
cd /home/erikg/impg
export CFLAGS=-I/home/erikg/.cache/impg/native-build-deps/root/usr/include
export LIBRARY_PATH=/home/erikg/.cache/impg/native-build-deps/root/usr/lib/x86_64-linux-gnu
export LD_LIBRARY_PATH="$LIBRARY_PATH"
cargo build --offline --release --bin impg
out=$(mktemp -d /home/erikg/yeast/partition-scheduling-profile-XXXXXX)
sha256sum target/release/impg > "$out/binary.sha256"
# Do not inherit the test suite's RAYON_NUM_THREADS=2: original real run uses 16.
# No --debug, frequency masking, threshold changes, or larger execution budget.
env -u SYNG_PROFILE -u SYNG_EMIT_PROFILE \
  RAYON_NUM_THREADS=16 IMPG_PARTITION_PROFILE_MAX_QUERIES=256 \
  IMPG_PARTITION_PROFILE_SKIP_QUERIES=0 \
  timeout --signal=TERM --kill-after=10s 120s target/release/impg partition \
  -a /home/erikg/yeast/syng-k63-s8-seed7-acgt-only/yeast235.syng \
  -w 10000 -d 1000 \
  --starting-sequences-file /home/erikg/yeast/partition-pilot-w10k-d1k/starting-sequences.txt \
  --selection-mode haplotype --syng-min-chain-anchors 5 \
  --syng-min-chain-fraction 0.5 --no-rehome-singletons \
  --sequence-files /home/erikg/yeast/yeast235.agc -o bed --separate-files \
  --output-folder "$out/results" -t 16 > "$out/stdout.log" 2> "$out/run.log"
printf '%s\n' "$?" > "$out/exit-code.txt"
grep 'partition_profile_' "$out/run.log" > "$out/profile.log"
```

For a later sample, change only the diagnostic skip/count after determining
backend dispatch ordinals. Partition numbers are not dispatch ordinals. Skipping
profiling does not skip processing or jump the scheduler to a late genome; a
short replay may never reach the desired sample. This contiguous sample can miss
repeat visits before/after it, so do not extrapolate exact-repeat frequency to
the entire run without qualifying the sampling boundary.

## First parent real-data sample

The first256 dispatched queries were profiled during a120s unchanged-parameter
replay on the ambiguity-safe index. All256 completed successfully and contributed
new source bases; no exact repeated tuples or zero-output queries occurred in
this sample. Three had zero missing core, yet together emitted167,318 additional
source bp. This is real evidence against treating covered source as redundant
work, not just the fake-backend counterexample.

Query time totaled80.382s: lookup76.200s (94.8%) and chaining/filtering4.180s.
Postprocessing totaled0.092s, including0.077s masking; output accounting/writing
was0.135s. Nested timings must not be added to their enclosing stages. The sample
contained19,782,419 raw anchors and emitted548,152,210 source bp.

The slowest sampled query was `S288C#0#chrIV:880000-890000`: zero missing core,
357,700 raw anchors,2.270s lookup,0.182s chaining, and40,931 newly emitted bp.
`S288C#0#chrIV:990000-1000000` was also productive despite zero missing core:
341,349 raw anchors,2.074s lookup,0.169s chaining, and6,096 new bp.

Artifacts: `~/yeast/partition-scheduling-profile-4wj5_t11/` contains the frozen
binary/source patch, command/index provenance, logs, profile-summary.json and
partial-output QC. The whole120s replay timed out incomplete after350 written
partitions and748,308,974 emitted bp; bounds and coordinate-union checks passed.
The256-query sample is only the early portion, not a measurement of the later
ANM#4 phase. It does not settle late-stage scheduler behavior.

**Next measured target:** isolate occurrence lookup/coordinate recovery and
anchor grouping on these slow queries. Preserve normalized outputs, all copies,
coordinates and thresholds. No covered-core skip or exact-query reuse is
justified by this sample. The general reuse idea below remains conditional.

## Proposed next repair, not implemented or authorized

Measure repeat rates first. If exact query tuples are materially repeated,
investigate **syng-only exact-query reuse** with immutable index and unchanged
query/filter/context parameters. It requires proof that later error behavior and
all mask/boundary/minimum-size effects are equivalent under monotone ownership.
In particular, a prior result's hull is not an adequate substitute for its
occurrence coordinates, and postprocessing may depend on the current missing
regions. No capability API, cache, visited-node suppression, or coverage-driven
policy change is included here. Otherwise retain explicit discovery obligations
and seek a separate policy decision; do not silently shrink residual queries or
infer no homology from skipped work.
