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
