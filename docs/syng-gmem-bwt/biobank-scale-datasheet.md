# The biobank-scale cost datasheet (post routes-lever)

> **SESSION-OPEN INCREMENT — the projection bandwidth campaign (opened
> before any measurement; this section is the plan of record, the
> measurements will land as a later increment).** The owner's follow-up
> ruling on the projection pole, verbatim: **"projection 25 minutes seems
> a bit long... shouldn't even be running 192 cores"** — with the owner's
> two observations to explain and exploit: (a) the committed 35x projection
> is DRAM-bandwidth-bound, not compute-bound (568,680 reads/s with only
> ~27-37 of the 192 configured cores busy — the receipt's own statement),
> and (b) a 12-lane shape would hold a PRIVATE ~46GB in-RAM copy of the
> panel syng index per lane (~550GB of duplicated RAM, all lanes
> contending for the same DRAM bandwidth — the committed run's one shared
> process is the counter-shape). The campaign's hypotheses, stated before
> measurement: **H1 (locality batching)** — a cheap pre-pass that assigns
> each read to a panel REGION (the grain: the input stream's own region
> grouping, chromosome-grain, or finer from the index's own node space),
> then bucket-by-bucket projection so concurrent threads walk the same
> index pages (hot cache) — worth several-fold on a bandwidth-bound
> process; **H2 (one shared index)** — a single process with ONE read-only
> in-RAM index (all threads share the pages) beats N private per-lane
> copies on both RAM and effective bandwidth; **H3 (the saturation
> curve)** — reads/s vs thread count at 1,2,4,8,16,32,64,128 threads on
> ONE index copy measures the bandwidth ceiling, confirms the 16-core
> observation, and quantifies the headroom locality batching must unlock.
> The measured object comes FIRST (the campaign's rule): a bounded slice
> of the 35x BAM (5-10M reads), one shared-index process, the thread-count
> sweep, cores-busy and (if cheaply available) the DRAM counters. The
> identity gate is absolute: the optimized path's occurrence records must
> be identical to the committed projection for the measured slice (same
> records, same occurrences; the committed shard format's within-batch
> packet order is the send-in-par scheduling's own nondeterminism, so the
> standard is record/multiset identity, stated); a dropped or duplicated
> multi-mapping read is a FAILURE, not an optimization. No thresholds, no
> tuning constants (the thread counts are measurement points); the
> committed census/shards are the immutable baseline; no implementation of
> the streaming restructure without the owner's gate.

The owner's ruling, verbatim: "biobank scale - the diploid genome for every
genome in a biobank relative to the HPRC; no true realignment of reads, a
projection onto a pangenome + local haplotype-aware genotyping; SOME MINOR
LOSS OF ACCURACY IS ACCEPTABLE (must be measured, not assumed); the goal is
minutes-class whole-genome inference per sample on this box. 'Getting to one
to two hours per sample would be great. I would pull the next lever. Let's
try to get that to work... It should be effectively faster than doing a
realignment, much faster... ideally on a box like this a few minutes to do a
single sample.'"

This section is the measured cost datasheet after the routes lever
(commit f3cbb90; receipts at `/home/erikg/hprcv2/routes-lever/L1/`).
Every number below is a committed receipt or a wall measured in the lever
session; arithmetic is shown; the two arithmetic-scaled cells are labeled.

## 1. The panel/sample cost split (the survey)

Per-territory step of the re-anchor-stage chain (L1, the 250kb territory,
k=63, the committed 300x fetch chain), split by WHAT IT TOUCHES:

| step | wall (measured) | peak RSS | inputs | side |
|---|---|---|---|---|
| local syng build | 3s | small | panel sequences + tiling | PANEL |
| catalog | 9-10s | 2.5GB | beds + panel syng | PANEL |
| routes (native build) | 50-61s | 36MB | panel + catalog + sources + read-lengths | PANEL (the builder's own contract: "Sample-independent builder") |
| **routes (the lever: derive-routes)** | **4.8s** | **26MB** | panel + sources + read-lengths (+ the census as the routing-bind cross-check) | PANEL (precompute; census-independent, proven) |
| re-anchor | 3.8s | 164-228MB | the committed census + the two syngs + territory | SAMPLE |
| derive-cache | 11-12s | 1.3GB | the sample FASTQ + the local panel | SAMPLE |
| CLI call | 30-40s | 3.4GB | re-anchored census + FASTQ + cache + maps | SAMPLE |

The survey's decisive finding: **the calling machinery consumes from a
routes directory ONLY `graph.json`'s lane inventory (panel path
id/name/length) plus the source paths** (`realign_score.rs:3365-3374`,
`diplotype_chain.rs:460-469`, the derive builders likewise). The
record-to-locus routing the call consumes is the **census's own re-anchored
`partitions` field** (`occ.partitions` at `realign_score.rs:3525/3911/4693/
4805`). The 9,465,618 native ports and every per-read-length start profile
the 50-61s build derives are the route-search machinery's artifacts
(`SearchPanelRoutes`/`EvaluatePanelRoutes`), never read by the calling path.

Artifact masses per territory (L1, measured): the native routes build emits
**2.6GB** (ports.bin 723MB + per-source ports 724MB + native profiles 1.2GB);
the catalog is **1.3GB**; the lookup routes artifact is **532KB**; the local
syng is **3.2MB**; the locality maps are **1.1MB**. At ~12,400 territories:
the built-routes class alone would be ~32TB and the catalog ~16TB - neither
is a durable-asset candidate; the compressed panel-side precompute is
**~4.8MB/territory (~60GB genome-wide), the maps-class**.

## 2. The routes lever (implemented, gated)

`impg genome-infer derive-routes` (the sibling of `reanchor-census`): emits
the call-consumable routes inventory BY LOOKUP - every lane from the target
panel's own name map (BY NAME, the re-anchor principle), the per-lane source
FNV recomputed over the same bytes (one linear read per lane), honest empty
sealed native artifacts (zero ports, zero profile runs, stated in
`generation_rule` and in the receipt; the search machinery fails closed on
the artifact by design), and the re-anchored census cross-checked (every
occurrence's path id must resolve in the panel name map; L1 receipt: 10,787
records / 4,133,359 occurrences / 466 paths / 25 partitions).

The gates, all passed (receipts at `/home/erikg/hprcv2/routes-lever/L1/`):

- **Unit**: the new test proves the call-consumed fields of the lookup
  graph identical to the built graph's on the shared fixture, the honest
  zeros, the fail-closed, and the absent-path bind refused
  (cargo test --release --lib: 541 passed, 0 failed).
- **Cache identity**: the L1 derive caches built under both routes
  artifacts are BYTE-IDENTICAL (112MB, cmp clean).
- **Answer-equivalence at L1**: the compressed path's `calls.jsonl` AND
  `calls.jsonl.truth-qv.jsonl` BYTE-IDENTICAL to the committed receipt of
  record (`k-sweep/L1-k63/cli`) under BOTH the built-routes control (same
  binary, same census, the native build's own artifact - the isolation
  leg) and the lever chain (lookup routes + lookup cache). **No accuracy
  delta to quantify - the gate is byte-identity, the strongest form.**
- **Census-independence**: a 100-record census head derives a `graph.json`
  BYTE-IDENTICAL to the full census's - the artifact is a pure function of
  (panel, sources, read-lengths); the census is the routing-bind
  verification only. Hence the artifact is a panel-side precompute: build
  once per territory, reuse for every sample forever.

The measured lever: **routes 50-61s -> 4.8s end-to-end** (2.7s the panel
core: syng load + the 117.0Mbp source FNV + seals; 2.1s the census gz
cross-check), RSS 26MB - or **~0 per sample** as the panel-side precompute.
The derive-cache does NOT fold from the census (the census is record-level:
record/multiplicity/anchors/occurrences - no read identity, no read
sequences); it remains the sample-side per-territory MEM derivation, now
served by the lookup artifact (cache byte-identity proven).

## 3. The post-lever per-territory datasheet

L1-class territory (250kb, 466 HPRC lanes, 25 windows):

| pole | 300x fetch chain (measured) | 35x operating point |
|---|---|---|
| re-anchor | 3.8s (0.44s lookup core, 2.5s gz-6 of the 370MB stream) | ~1s (arithmetic: the 35x census mass is ~6x smaller - cut-L1 kept 108,230 reads vs 649,439 at 300x) |
| derive-cache | 11-12s | ~2s (arithmetic, same ~6x read mass) |
| call | 30-40s (phases below) | **20s MEASURED** (the committed cut-L1 CLI receipt; the pilot-class heavy tail measured 85s at 425,507 kept reads; L2-L4 16s) |
| routes | 0 (panel-side precompute; or 4.8s per sample if derived in-path) | same |
| **sample-side total** | **~45-55s** | **~23s (20 measured + ~3 arithmetic)** |

The call's measured phase table at 300x (the `[score] phase summary` line of
the lever run): inputs 2.6s (census gz parse + panel syng + maps), FASTQ
quality scan 3.9s, binding 4.2s (12.4M verify probes, lane loads 111.6MB),
variants 0.3s, rows/context assembly 4.7s (562.6MB of lane fetches), loci
(folds+units+scoring+classes+qual+exactness+receipts) 18.9s.

## 4. The panel-side one-time precompute (per panel territory, once forever)

syng build 3s + derive-routes 4.8s ≈ **8s/territory serial**; at
~12,400 territories ≈ 27.6h serial, parallelizable 256-wide to **~7-10
minutes one-time for the whole HPRC territory set**. Durable assets kept
(the owner's disk ruling): syng 3.2MB + routes-lookup 532KB + maps 1.1MB ≈
**~4.8MB/territory**. The catalog (1.3GB) and the native routes build
(2.6GB) are AVOIDABLE in the compressed path entirely - the call never
consumes them (the gates prove it).

## 5. The whole-genome extrapolation (measured walls, arithmetic shown)

Per sample (35x operating point), post-lever:

1. **The projection, once per sample**: 24.5min measured (1,467s; 734M
   reads / 108.7Gbp; 568,680 reads/s; peak RSS 41.0GiB). The receipt's own
   statement: memory-bandwidth-bound, ~27-37 cores busy - the projection is
   NOT serial by construction, but its throughput saturates memory
   bandwidth; additional projection lanes contend for the same bandwidth.
   It leaves ~85% of the box's 256 cores free while it runs.
2. **The territory wave**: 12,400 territories x ~23s sample-side =
   285,200s ≈ 79h serial. Lane width: the per-lane peak RSS is the call's
   3.4GB (derive 1.3GB, re-anchor 0.23GB) - the box (256 cores, 1,005GiB
   RAM) fits ~250-295 RAM-lanes and 256 core-lanes; **at 256 lanes the wave
   is 79h / 256 ≈ 19min** (870GB of the 1,005GB - tight but fitting; a
   comfortable 230 lanes is ~21min).
3. **Sum**: 24.5min (projection) + 19min (wave) ≈ **44min serial-sum**;
   overlapped (the wave grinds on the ~220 cores the projection leaves
   free) the honest wall is **~25-30min per sample at 35x**.

At the 300x fetch chain the same arithmetic is ~55s x 12,400 = 189h serial
/ 256 = 44min + 24.5min ≈ ~1.1h - the "one to two hours" target met at
both depths; the distance to "a few minutes" is carried by the two
measured poles below.

## 6. The floors analysis (each remaining pole, its measured number, its lever)

- **THE PROJECTION POLE - 24.5min once per sample (measured,
  memory-bandwidth-bound)**: 568,680 reads/s over the in-RAM syng index,
  ~27-37 cores busy of 192, peak 41GiB. The honest ceiling: this is
  one full pass of the sample (108.7Gbp at 35x) through the FM index at
  memory bandwidth; it does not go faster by adding lanes. Its levers, in
  order: (a) it is ONCE per sample and overlaps the territory wave (its own
  receipt leaves ~85% of the cores free); (b) the cut-batch receipt (413s
  whole-batch at width 1) amortizes the global load+scan across a batch of
  territories - the per-territory cut cost tends to the scan share only;
  (c) the depth lever is the owner's: the receipt's own scaling (30x ->
  0.35h) - the projection wall is proportional to read mass.
- **THE CALL POLE - 20s/territory measured at 35x (16-20s typical,
  85s at the pilot-class heavy tail)**: the phase table at 300x decomposes
  it: loci scoring 18.9s (read-mass-proportional), rows/context 4.7s
  (562.6MB of panel lane fetches - panel sequence, decompressed per call),
  binding 4.2s, inputs 2.6s (the census gz parse + panel syng + maps),
  FASTQ scan 3.9s (read-mass-proportional). The fixed-cost floor at 35x is
  the ~5-7s load class. Named levers: the lane-fetch mass is panel
  sequence - a per-territory lane cache (the maps-class, MBs) removes the
  per-call decompression; the census gz parse/emission class is streamable
  (the re-anchor stage measured gz level-6 at 2.5s of its 3.8s - level-1 or
  raw jsonl cuts it to ~0.4s); the scoring mass scales with territory read
  mass and is the irreducible model work.
- **THE DERIVE-CACHE POLE - 11-12s at 300x, ~2s at 35x (arithmetic)**: a
  pure function of (FASTQ, panel); does NOT fold from the census (the
  census is record-level, no read identity/sequences). Its lever: the
  committed projection already carries per-read packets (sequence + tokens
  + walk) against the global panel - a territory-local cache derivation
  from those packets (bind-checked) would remove the local MEM
  re-derivation; NOT implemented here (the global->local token grid bind
  is the same owner-gated cross-grid class the k=31 bind measured).
- **THE RE-ANCHOR POLE - 3.8s at the 300x-class census (0.44s core)**: the
  lookup core is sub-second; the pole is the gz level-6 emission of the
  370MB stream (2.5s) - streamable (level-1/raw), the same lever as the
  call's census parse; at the 35x census mass ~1s end-to-end.
- **THE GZ/IO CLASS**: the emission/parse class (re-anchor 2.5s emit; the
  call's census parse ~2s of inputs) is streamable end to end - the one
  cheap, measured lever that touches every per-territory step.

## 7. The honest verdict

Post-lever, whole-genome per-sample at 35x on this box (256 cores, 1TB
RAM): **~44min serial-sum / ~25-30min overlapped** - inside the owner's
"one to two hours" with the answer-equivalence gate passed at
byte-identity (no accuracy delta to state). The distance to "a few
minutes" is carried by exactly two measured poles: the once-per-sample
projection (24.5min, memory-bandwidth-bound - its lever is depth and
overlap, both the owner's call) and the per-territory call (20s measured
at 35x, of which ~5-7s is the fixed load class with the named
lane-cache/streamable-census levers). The panel-side precompute is
one-time, ~7-10 minutes of box time for the whole HPRC territory set, and
~4.8MB/territory of durable assets.

Assessment-side truth only (the --truth-qv test mode); no thresholds; no
tuning constants; the committed receipts unmodified; the scoreboard
machinery untouched.

## 8. The projection bandwidth campaign (closed; the pole measured, the
##    lever committed, the next levers named)

The owner's ruling at the head of this datasheet: **"projection 25
minutes seems a bit long... shouldn't even be running 192 cores."** The
campaign measured the pole, shipped one committed lever, and answered the
question honestly. Co-tenant DRAM hammers (rpfbwt, bench-numa) ran the
whole session: **±20% noise on every rate below**; matched back-to-back
points are the primary evidence. No thresholds, no tuning constants —
every thread count is a measurement point.

**The measured object** (the campaign's own rule): a bounded slice of the
35x BAM — 19 × 409,600 mid-chromosome reads (the committed 19-fifo
feed's own stream windows) = 7,782,400 reads, exactly 9,557,434 records
in every gated run. Receipts: `/home/erikg/hprcv2/proj-bandwidth/`.

**H3 — the committed binary's saturation curve** (one shared index,
reads/s post-load, T=1..256): 71,276 / 111,972 / 175,278 / 290,502 /
500,782 / 464,039 / 549,868 / 566,308 / **204,380 → 503,609 → 592,328**
(the committed 192-thread point's first run was a cold-cache artifact;
immediate re-runs are the honest number) / 583,048. **The plateau is
~500-600k reads/s from T=16 up; the committed 192-thread configuration
is past the peak; per-busy-core efficiency collapses ~5-10x beyond
T=16.**

**H1 — the locality lever, v1 → v2 (committed as `--region-batching`,
default off, the committed path byte-for-byte).** v1 (strict per-stream
batches): +44/54/61% at T=1/2/4 but **−22% at T=192** — single-writer
burst contention, diagnosed and replaced. v2 (each 16,384-read arrival
batch sorted by input stream before the rayon dispatch — the committed
feed's own chromosome-per-stream grouping, zero pre-pass cost, packet
sends still spread across all 19 channels): s32 **709,814 reads/s** —
above the whole baseline plateau; curve 56,883 / 145,106 / 128,329 /
256,719 / 574,965 / **709,814 (T=32)** / 580,794 / 669,877 / 295,036 /
190,218.

**The identity gates — all PASS, any failure would have been a
failure.** G1: slice vs committed whole-sample shards, 19/19 streams,
packet/multiset identity with both jitter-guard legs. G2: the rebuilt
binary's default path ≡ the committed binary's output. G3: default ≡
v1-strict. g4/g5: the v2-sorted runs at T=32 and T=192 ≡ the default
control. The commutative 128-bit per-shard multiset hash: identical on
all five slice runs, 19/19 shards, records exactly 9,557,434 everywhere.

**H2 — the decisive 12-lane leg (the owner's shape, measured).** 12
processes × 16 rayon threads, each a **private 35.5GiB in-RAM index
copy** (12 × 35.5 = 427GiB resident; the 120GiB guard never came close,
minimum available 695GiB; guard math stated before launch), disjoint
record-aligned splits of the same slice (records sum exactly 9,557,434;
the union packet multiset ≡ the 19-shard multiset by the commutative
hash). **Per-lane post-load 338k-453k reads/s; aggregate ≈4.74M reads/s
— ~6.7x the best single-process point (709,814), ~8.3x the committed
receipt's 568,680.** VERDICT: **the single-process plateau is a
PER-PROCESS (single-address-space) ceiling, NOT the box's aggregate DRAM
bandwidth — the box feeds 12 lanes at 4.7M reads/s.** The committed
receipt's "DRAM-bandwidth-bound" statement was the per-process
observation; H2 corrects it at the box level (the exact per-process
mechanism — TLB/page-table/scheduler — is unmeasured, stated honestly).

**The full-sample optimized run (the real wall-clock number).** The
committed 35x runner verbatim (same BAM, same 19 region streams, same
fifos, same map-cap, same guards) + `--region-batching` at **T=32**
(the best gated point), one shared index, nohup + markers:

| | committed (T=192) | optimized (v2, T=32) | delta |
|---|---|---|---|
| projection wall (incl. load) | 1467.05s = 24.45min | **1247.5s = 20.79min** | **−15.0%** |
| index load | 175.6s | 178.4s | equal (noise) |
| post-load rate | 568,680 reads/s | **686,965 reads/s** | **+20.8%** |
| dedup wall | 118.4s | 114.5s | equal |
| end-to-end projection+dedup | 1585.4s = 26.4min | **1362.0s = 22.7min** | −14.1% |
| peak RSS | 43,013,472 kB | 43,021,228 kB | equal |

**The full-sample identity gate — PASS at every level** (the committed
projection is the immutable baseline): `occurrences.bin` (3.9GB
canonical dedup stream) **BYTE-IDENTICAL**; per-shard pattern-multiset
occ gates **19/19 IDENTICAL**; per-shard commutative packet-multiset
hashes **19/19 IDENTICAL** (734,443,930 packets, 838,221,858 records
— the 100% standard); projection-stats counters identical. Diagnostic
note, stated: index-window sampling between two full runs is NOT a
valid standard (packet order within a shard is the send-in-par
nondeterminism over arrival-batch stream interleaving — measured ±17
records per 409,600-packet window, the batch-edge jitter class,
receipts kept); the whole-file commutative multiset cancels it and is
the standard.

**The datasheet, before/after.** Thread/RAM tenant: **one shared ~41GiB
index at T=32** replaces the committed past-peak T=192 tenant at equal
RAM — and the 12-lane shape would hold **12 private copies ≈ 427-492GiB
of the box's 1,006GiB**. Per-sample end-to-end (projection + wave at
230-256 lanes, serial-sum): **~40-42min** (was ~44min; the wave's
~19-21min is now the dominant pole). The overlapped bound (streaming
restructure below, owner-gated, NOT implemented) would be
max(projection, wave) ≈ **~21min** at the s32 shape.

**How can this be faster? (the honest answer, each lever with its
measured motivation)**

1. **Done and committed**: `--region-batching` at T=32 — +21%
   post-load bandwidth, −15% whole-sample wall, byte-identical output,
   same RAM, fewer threads. The committed 192-thread configuration was
   past the peak; the same box and index do 21% more at 32 threads.
2. **The lanes shape (H2, measured on the slice)**: 12 private-index
   lanes × 16 threads aggregate **~4.74M reads/s** — the extrapolated
   full-sample projection is load ~200s + 734.4M/4.74M ≈ 155s ≈
   **~6min** (vs 20.8min measured single-process; extrapolation from
   the slice, NOT measured at full scale — stated). Cost: ~492GiB of
   private index copies, and **the box cannot co-run the wave's 870GB
   lane fleet** — either lane fleet, not both.
3. **The streaming restructure (the named NEXT lever, owner-gated, NOT
   implemented)**: region-ordered emission + per-territory routing at
   emission + a region→offset manifest would let the per-territory
   wave tail the projection stream — hiding the whole projection pole
   behind the wave (per-sample ≈ max(projection, wave) ≈ 21min).
4. **Index packing (named lever)**: the index is 35.5GiB resident
   (GBWT: 272,213,874 vertices / 683,539,246 edges) and its load is
   178s — 14% of the s32 wall and the dominant term of the lanes
   shape's per-lane economics. A packed layout shrinks the per-lane
   private copies (12 × 35.5GiB → proportionally less), cuts the load,
   and improves per-process cache/TLB locality.
5. **NUMA placement (named lever)**: per-lane rates varied 338k-453k
   (±15%) under co-tenant contention; node-affine placement of each
   lane's private copy (one socket's DRAM per lane) removes cross-socket
   DRAM traffic — the 4.74M aggregate is a floor for the un-pinned
   shape; the wave's 870GB lane fleet needs the same discipline.
6. **The wave is now the pole**: ~19-21min vs projection 20.8min
   single-process (~6min lanes-extrapolated) — the next end-to-end
   gains live in the wave's lane economics and the fixed-load levers of
   section 6.

Standing rules held throughout: the committed projection/shards were
never touched; identity failures would have been failures; no
thresholds or tuning constants; NO PR PUSH (the owner pushes); the
64GiB single-process guard and the 120GiB multi-process guard armed and
never tripped; intermediates cleaned, receipts kept at
`/home/erikg/hprcv2/proj-bandwidth/` (the gates machinery committed at
`genome/instrumented/projection-*`).
