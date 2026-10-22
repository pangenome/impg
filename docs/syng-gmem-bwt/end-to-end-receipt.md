# THE END-TO-END RECEIPT — the full three-step product at the current stack, HG002 35x vs the committed panel

> **WHAT THIS IS.** The measurement step of the chain: the product
> itself exercised — not an experiment. The three-step UI verbatim
> (panel-prep once → map per sample → call-region per territory, zero
> product flags), the full accepted territory universe (5,844, the
> runtime-verification split ledger), HG002 35x against the ONE
> committed global panel (`HPRC_r2_assemblies_0.6.1.syng`, metadata v1
> — `IMPG_PANEL_LEGACY_V1_OK=1`, the committed campaign's own
> substrate acknowledgment), at the CURRENT stack: every landed lever
> in its default position (the bounded class scan, the tie-set fix,
> the bit-equal-row memoization, the precompute slices 1+2, the
> far-bin prune, the restricted census profile, the fleet form
> COW/T16/width-32).
>
> **THE STREAM DISCIPLINE IS LAW HERE** (stream-discipline.md): the
> ANSWERS are the product (calls.jsonl + the sequences emission, kept
> per territory); every receipt is computed at write time; the raw
> instrument streams are the reaped class (the tranche-2 lesson:
> ~600MB/territory of instrument receipts beside a 35KB answer —
> mined to counters, then reaped); the per-territory transients
> (substrate, precompute asset, cut census/cache/reads, the region
> fetches) are reaped in-pass at each leg's completion; the workdir
> budget is 150G per scratch root, enforced inside every runner; df
> watched by the heartbeat.
>
> **THE RUN IS MULTI-WINDOW** (the honest shape at genome scale —
> the pre-campaign estimate itself was ~6.38e6 core-s of call work):
> the runner is marker-gated and resumable; a successor RESUMES,
> never restarts. Each window checkpoint-commits the per-territory
> receipts that landed in it; the ledger states exactly what is
> measured at each n, and nothing is extrapolated beyond its receipt.
> PENDING RECEIPTS ARE MARKED **[PENDING]** until their markers land.

## 0. The run's exact shape (the commands, as the product defines them)

The campaign scripts are committed at `scripts/end-to-end/`
(commit 25ededa); the workdir is
`/home/erikg/hprcv2/end-to-end-receipt` (the receipts under
`receipts/`, the ledger under `ledger/`, the answers under
`calls/`).

- **STEP 1 — PANEL-PREP** (once per panel + once per territory):
  1. the census-producing cut (the routing bind's input; the
     committed global-frame censuses were reaped with the gates-era
     workspace): the committed global-frame cut at L1, d=1,
     global-native profile, RAYON=32, guarded at 96GiB.
  2. the whole-panel routes inventory: `impg genome-infer
     derive-routes --panel <the ONE panel> --sources <the panel AGC>
     --read-lengths 148 --census <the L1 census> --out-dir
     global-routes` — the parallel lane-hash form, once per panel.
  3. per territory: the export-profile substrate chain (the committed
     runtime-verification recipe VERBATIM on the committed substrate
     instruments — extract → syng → w10k/d1k beds → export graphs →
     panel walks → locality maps, K=63) then `impg genome-infer
     panel-prep --name <TID> --chrom … --start … --end … --sources
     <locus.fa> --component <axis path> --axis … --partition-graphs
     … --census … --read-lengths 148 --global-panel … --global-agc …
     --global-routes … --record-origin 1 --out-dir partitions/` at
     the current stack (routes REUSED loudly; the partition global
     form; the precompute asset build, its once-price receipted).
     The content origin d=1 is the uniform substrate convention (the
     delta=1 AGC-identity receipts at BOTH substrate classes: L1 466
     records, T00000 467).
- **STEP 2 — MAP** (per sample): the projection (the committed
  20.79min recipe verbatim: the local BAM, the 19 region streams,
  `-F 0x904`, the fifos, `--region-batching` at T=32, on the
  committed projection instrument), then the dedup stream, then the
  three-level identity gate vs the committed shards; per territory
  the map leg (the committed ±50kb reads fetch rule + `impg
  genome-infer map` — the local syng built into the map's own
  partition dir + the committed LOCAL-frame cut against the ONE
  global panel; the territory's `cut-reads.fasta` is the call step's
  restricted-census basis).
- **STEP 3 — CALL-REGION** (per territory, the fleet form): one
  `impg genome-infer call-region --map <chunk root> --partitions
  partitions/ --region <32-territory list> --out-dir calls/chunk-N`
  per 32-territory chunk, ZERO product flags: the parent loads the
  ONE panel and stays the load owner (fork-shared COW); the width-32
  child builds its own AGC index, batch-cuts the chunk at the
  RESTRICTED census profile (the basis auto-derived from each
  territory's map cut), and runs the 32 calls concurrently at the
  T16 rung with the far-bin prune default ON. The 64GiB guard over
  every worker (a guard kill is a named fail-closed census entry; the
  leg re-runs on resume). The chunk's reads fetch is the committed
  ±50kb rule per territory, `-M` multi-region (a superset per member
  — harmless to the census-keyed reads lookup; the map roots are the
  map legs' own outputs, symlinked — the batch-composition
  invariance receipt is the warrant).

## 1. STEP 1 receipts — panel-prep

**THE CENSUS-PRODUCING CUT** (the routing bind's input; the committed
  global-frame censuses were reaped with the gates-era workspace):
  wall **491s** at RAYON=32, guarded at 96GiB (no event), and **the
  re-derived census reproduces the committed L1 d=1 receipt at every
  counter** — reads_scanned 734,443,930 / reads_kept **304,793** /
  census_mult **384,056** / routed **3,775** / territory rows 11,649 /
  entries 2,757,852 / load 262.7s / scan 161.1s — the committed
  attempt-3/§10.1 values exactly. (Core-s not captured on this one-time
  step — the guard-pid bug of the first runner generation, §3's guard
  census; the load/scan decomposition is the machinery's own receipt.)

**THE WHOLE-PANEL ROUTES INVENTORY** (derive-routes over the ONE
  committed global panel, the parallel lane-hash form, once per
  panel; the reaped artifact re-derived):
- wall **1180.5s** (the machinery's own elapsed receipt; the committed
  parallel baseline 1505s — this box quieter), at RAYON=192.
- **core-s: 163,159** (user 161,385.7 + sys 1,773.6).
- **peak RSS: 711GB** (memavail floor 255GB during the run) — THE
  PARALLEL LANE-HASH FORM'S OWN MEMORY APPETITE AT 192 THREADS, never
  receipted before (the committed serial form's was ~45GB-class at 1
  thread; the per-thread AGC fetch buffers are the term). A finding
  for the panel-side once-price ledger: the 25-min form is
  memory-hungry, and at fleet co-tenancy it squeezed a 1TB box to
  within ~75GB of the floor. The step ran un-guarded (the first
  runner generation's guard watched the wrapper, not the worker —
  found and fixed the same window, §3's guard census).
- **THE IDENTITY RECEIPT**: lanes **38,790 — MATCHES the committed
  receipt** (the serial-vs-parallel byte comparison's own target);
  total_source_bp **1,403,221,029,691 = 1.403 Tbp — MATCHES the
  committed receipt**; families 466; the routing bind over the fresh
  census: 3,775 records / 1,448,220 occurrences, all resolved;
  panel_metadata_version=1, panel_legacy_opt_in=True (the loud
  provenance line in the run log, kept as
  `receipts/routes.runlog.txt`). The committed serial-vs-parallel
  byte-identity receipt (`serial-vs-par.routes-cmp-receipt.txt`) was
  REAPED with the gates-era workspace — NOT AVAILABLE; the
  manifest-field match above is the available identity receipt, and
  the re-derived census's exact-counter match (above) is its
  cross-check.

**THE PER-TERRITORY PREP** (the substrate chain, the committed
  export-profile recipe VERBATIM on the committed instruments; the
  first 192 legs): per-step walls at T00000 — extract 6s (1.9u+5.4s
  core-s), syng 2s (2.4), beds 5s (18.0), export 2s (45.5), walkreq
  0s, walks 1s (1.4), maps 15s (14.6) — **31s wall / ~90 core-s** at
  THREADS=4 (the campaign's own T00000 smoke: 79s serial). The
  fleet's per-step core-s receipts land in the ledger as the legs
  complete.

**THE PER-TERRITORY PANEL-PREP** (the three-step product's own STEP 1
  invocation at the current stack): **~250s wall / ~251 core-s per
  territory** (T00000: wall 251s, cpu 206.1u+45.0s), dominated by the
  GLOBAL PANEL LOAD for the precompute asset build (the ~200s class,
  §10.6.24's own receipt) — **a once-per-panel-per-TERRORITY cost the
  product currently pays per invocation: 5,844 × ~250 core-s ≈ **1.46e6
  core-s**, the panel-prep LOAD POLE — the same class the fleet form
  amortizes on the call side; a batched panel-prep is the named lever
  (a finding of this run, stated in §4)**. The precompute asset
  itself: rows.bin 90-99MB + windows.bin ~0.1MB + bind-scaffolding
  30-35MB per territory (~125-165MB), built in ~8s marginal beside
  the load; the REUSE is loud and cheap (the redo after the bed-naming
  fix re-ran panel-prep over every prepped territory at **~0-15s
  each** — the asset is content-keyed and REUSED, the partition step
  is seconds).

**THE SPARSE-COMMUNITY BED-NAMING BUG — found by this run, fixed at
  the product, gated by the run itself** (commit 0c2ae7a): the
  export-profile axis names its windows by the ORIGINAL w10k/d1k
  community ids (T00017: partition0, partition21, partition22,
  partition23) while the locality maps renumber 0..k; the panel-prep
  bed regeneration named the beds by the map index, and **17 of the
  first 22 map legs died at `partition21.bed` NotFound** (the dense
  class — T00002/T00004/T00005, and the five demo localities by the
  luck of their substrates — succeeded). THE FIX names each
  regenerated bed by the AXIS INTERVAL'S OWN GROUP, pairing maps to
  intervals BY ORDER with the component-row-extent cross-check,
  fail-closed loudly; verified across both substrate classes before
  landing (L1 dense 25; T00017 sparse 4; T00032 mixed 17; T00000
  gapped 22) and proven by the run: T00017's map leg now completes
  (369.8s, reads_kept 65,241 / census_mult 62,338 / routed 4,140).
  All prepped registries were discarded and re-run through the fixed
  product (the precompute assets REUSED, ~0-15s each) — the bug's own
  cost, receipted here.

## 2. STEP 2 receipts — map

**THE PROJECTION** (the committed 20.79min recipe verbatim: the local
  BAM, the 19 region streams, `-F 0x904`, the fifos, `--region-batching`
  at T=32, on the committed projection instrument) into a fresh
  out-dir, then the dedup stream, then THE IDENTITY GATE vs the
  committed shards — **PASS AT ALL THREE LEVELS**:
- `occurrences.bin` (the 3.9GB canonical dedup stream) **BYTE-IDENTICAL**;
- projection-stats counters **identical** (reads 734,443,930,
  records 838,221,858);
- per-shard commutative packet-multiset hashes **19/19 IDENTICAL**
  (packets + records + the 128-bit multiset sum per shard).

THE WALL: projection **1563s** (the committed baseline 1247.5s =
20.79min — this run HEAVILY CO-TENANT with the 711GB-peak routes
  derive, stated honestly) + dedup 140s (committed 114.5s) = **1703s
  end-to-end**. Core-s: sample_projection **30,811** (27,721u +
  3,090s) + the 19 samtools stream writers beside (io-bound; the
  committed recipe's own feed). Peak RSS 43,121,460kB (~41.1GB — the
  committed receipt's own class, 43.0GB). Post-load throughput ~536k
  reads/s under the same co-tenancy (the committed receipt's 686,965/s
  at T=32 is the quiet-box baseline). **The re-run's 65G shard set was
  REAPED after the gate** (the stream law: the gate is on digests;
  the committed shards — byte-identical by the gate — carry the
  product forward).

**THE PER-TERRITORY MAP LEGS** (the committed ±50kb reads fetch + the
  three-step map step: the local syng build ~3-25s + the committed
  LOCAL-frame cut against the ONE global panel, RAYON=32, the 64GiB
  worker guard live): the first legs — T00002 **430.0s**, T00004
  429.3s, T00005 416.6s (the dense class), T00017 **369.8s** (the
  fixed sparse class: reads_kept 65,241 / census_mult 62,338 /
  routed 4,140). **[PENDING]** the fleet sums as the legs complete
  (the ledger carries each leg's wall + core-s + cut-stats).

## 3. STEP 3 receipts — call-region, the whole genome

**[PENDING]** the per-territory receipts (inputs/bind/rows/loci
walls, core-s, the locus census counters, RSS peaks, the guard /
fail-closed census), at each checkpoint: legs done/total, the walls,
the running core-s sums.

## 4. The restated 16-core ledger

**[PENDING]** the per-step core-s sums (prep once-price, map, calls),
the honest 16-core projection (core-s/16, labeled projection from
measurement), the fail-closed census (which legs the guard stopped,
named), the per-stratum walls, and the comparison to the
pre-campaign estimate (the 6.38e6 core-s call estimate + the 20-min
16-core target) — what the landed levers actually bought at genome
scale, measured.

## 5. The exact resume point

**[PENDING]** the marker state, the in-flight lanes, the next
window's first actions.
