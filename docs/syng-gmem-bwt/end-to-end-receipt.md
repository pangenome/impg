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

**THE FLEET FORM IS ALIVE END-TO-END** (the first chunk, `chunk-early`
  over the 7 already-mapped chunk-0 members — the step-3 smoke, run
  ahead of the 183-chunk campaign numbering; the campaign's chunks
  proceed by markers): one `call-region` invocation, ZERO product
  flags —
- the parent loaded the ONE global panel in **180s** (33.9GB
  private, the load owner) and forked the width-32 child; **the COW
  receipt**: child-open Rss 32.28GB / Pss **16.14GB** / private **0.00GB**
  — the fork-shared panel works exactly as designed (the child pays
  only its own working set).
- the child built its own AGC index in ~0s, batch-cut the chunk at
  the RESTRICTED census profile (the basis auto-derived from each
  territory's map-leg cut-reads: e.g. T00002 11,345 reads; T00017's
  restricted census: reads_kept 60,720 / census_mult 90,101 /
  routed 4,235 — the global-frame attribution of the committed basis,
  receipted), then ran the calls concurrently at the T16 rung.
- **wall 1441s / core-s 18,123 for the 7-territory chunk** (≈2,589
  core-s/territory at this width — the parent load amortized /7; at
  the full width-32 the amortization improves ~4.5x on the load
  term, all else equal).
- **the answers**: calls.jsonl + sequences per territory (T00002 25
  loci / 37KB; T00017 4 loci / 4.3KB; … 118 loci over 7 territories,
  164KB total), model label `…-farbin-prune` — **the far-bin prune
  default ENGAGED, no flag typed**; the bounded class scan ENGAGED
  (the mined per-locus counters carry `classscan_bounded` with the
  derived underflow point −745.133 and the slack-corrected keep
  window; T00017's loci are flat-shape — kept 107,416 / skipped 0,
  the honest shape-reading, no pruning where the landscape is flat);
  the equal-column merge counters and the scan-memo census mined
  per locus into the ledger (the raw instrument streams REAPED —
  ~600MB/territory class at L1 → the mined counters + the chunk's
  gzipped log, the tranche-2 lesson applied).
- the JIT reaper reaped **~2G of transients** at the chunk's close
  (the chunk reads fetch, the scratch call root, the 7 territories'
  substrates + precompute assets + map cut dirs — the answers, the
  per-territory ledgers, the cut-stats receipts and the chunk log
  kept).

**THE GUARD / FAIL-CLOSED CENSUS, this window**: 18 GUARD events —
  map-cut workers killed over the 64GiB guard (peaks 67.3–69.7GB:
  the memory-tail cut class, the guard-as-oracle working — those
  legs re-queue and will honestly guard-kill again if the shape
  demands it; the campaign's ledger names them per territory) —
  plus the window's own infrastructural classes, all receipted in
  the doc's §1/§3 and the ledger: the double-driver EEXIST cascade
  (a zombie driver instance from the first relaunch — killed; its
  stale exit markers cleared; the bed rows rebuilt from the
  registries), and the runner-generation guard bugs (fixed in
  commits 0c2ae7a/966d0a2's runner set: the /proc ppid field-shift,
  the $()-pipe-hold that zeroed walls, the awk auto-vivification
  flood). The fail-closed census at the window-2 checkpoint: **10
  named substrate-terminals** (the 2-strike rule's own convictions
  — 8 extract off-by-one + 2 syng SIGSEGV, §3's window-2 receipts)
  plus the 26 archived stale markers (23 bed-row-loss map exits + 3
  window-2-reap map exits + the T00043 panelprep EEXIST), every one
  re-queued and cleared — not a machinery failure of the product.

**THE CAMPAIGN STATE AT THE WINDOW-2 CHECKPOINT** (the marker
  census, 2026-10-22T15:18Z, after the shepherd's fixes below):
  substrate 250/5,844, panel-prep 129/5,844, map 116/5,844, calls
  7/5,844 (chunk-early), fail-closed **10 named substrate-terminals**
  (the 2-strike poison rule's own census, below). The driver (v4 +
  the window-2 gate set) runs DETACHED at **MAP_W=16** (widened
  from 12 this window, receipted below): the measured 16-wide pace
  **32 map legs in 22 min ≈ 87/h** (vs ~64/h at 12-wide), walls
  587–691s per leg at 16-wide co-tenancy (12×T32 scan bursts on 256
  cores), memavail floor ~400G mid-wave (978G at census). The workdir
  **110G of the 150G bound** and now self-draining (the window-2
  reaps, below): the substrate's dead-weight pieces reap at
  panel-prep completion, locus.fa + the territory's reads fetch reap
  at the map leg's completion, and the chunk's JIT reaper takes the
  rest at the call — the steady state holds the in-flight window
  (~160 valve + chunk granularity) at ~0.4G/territory. **chunk-0 is
  one map leg from dispatch** (22 callable-mapped + 2 named poison +
  T00043's map in flight) — the first live exercise of the
  partial-chunk gate. **[PENDING]** the per-territory receipts as
  the 183 chunks drain.

### WINDOW 2 — the shepherd's runner-fix receipts (all found by
  the census, fixed at the runner, committed)

- **THE 23 STALE `map.exit` MARKERS — the BED-row-loss class,
  cleared and re-queued** (archived as `map.exit.stale-bedrows`): the
  13:41Z map wave dispatched against a `partitions.bed` mid-rebuild
  (the window-1 EEXIST-cascade cleanup was re-appending the lost
  rows) and failed 23 legs instantly at the region lookup — every
  row verified present now; the legs re-ran and completed at 16-wide
  (T00000 691s / T00004 587s / T00037 666s …). Without the clear,
  chunks 0–1 could never dispatch (the M predicate excludes
  map.exit forever) and ~64G of their transients could never reap.
- **THE EARLY-REAP INCIDENT — the shepherd's own mistake, receipted
  honestly**: the first cut of the budget valve reaped the substrate
  (incl. `locus.fa`) at panel-prep completion — but the MAP cut's
  own `--sources` input IS `locus.fa` (the local-frame cut
  materializes the locus.fa slices), so 7 in-flight legs died
  NotFound (T00181–187, archived `map.exit.stale-window2-reap`)
  and ~100 pending maps would have followed. REPAIRED by the prep
  chain's own re-derivation receipt: `rv-extract.py` re-derived
  `locus.fa` for all 114 damaged territories (~8s each, the
  definitional smalls survive), the partial map roots cleared, the
  legs re-queued. THE CORRECT REAP POINTS now in the runners:
  export-walks/export-graphs/partition-beds/prep-syng at
  panel-prep done; `locus.fa` + the reads fetch at MAP done (the
  map cut was locus.fa's last consumer — the call's restricted cut
  is global-frame, no locus.fa, verified in the CallRegion arm);
  locality-graphs (the precompute asset) at the chunk's JIT reaper.
  The one-time catch-up freed **37G (120G→83G)**; the workdir then
  re-grew ~17G for the re-derived locus.fa set — the price of the
  mistake, stated.
- **THE ENDGAME PARTIAL-CHUNK GATE — required for the drain, now in
  the runner set**: a permanently-failed member (a named poison:
  map.exit, substrate-terminal, or pp-terminal) previously blocked
  its whole 32-member chunk from ever calling — the campaign could
  never reach ALL TERMINAL. Now the driver's C-set gate and the
  chunk runner treat every member as TERMINAL when map.done OR a
  named poison; the chunk calls the CALLABLE SUBSET (the
  batch-composition invariance receipt is the warrant), NAMES the
  poison members in the chunk receipt, and reaps their transients
  at the chunk's close (receipts + ledger markers kept).
- **THE 2-STRIKE POISON RULES (prep + panel-prep shapes)**: a step
  that fails on a re-run is a deterministic terminal — before the
  rule, such legs retried forever (wasted cycles) and never entered
  the fail-closed census. First convictions, all verified
  deterministic across both strikes: **8 extract off-by-one
  terminals** (the rv-extract assertion `len N != bed N+1` at AGC
  getctg contig boundaries: T00013/52/57/115/175/192/212/372 — a
  REAL data finding for the substrate ledger: the refined-bed rows
  at 8 territories demand one more base than the AGC contig
  carries) and **2 syng SIGSEGV terminals** (exit 139 on the segdup
  tails: T00024/197 — the predicted syng segdup-tail prep-terminal
  class). All named in `logs/failures.log` + the ledger census.
- **T00043's EEXIST-cascade leftover repaired**: its panelprep
  completed in window 1 (the registry verified: frame=global,
  windows=21, beds present, partitions.bed row present) but the
  runner died before the marker; the re-run refused on the
  name-dup check, retrying forever. The artifact IS the product's
  own completed output — `panelprep.done` marked, the stale exit
  archived.
- **THE MAP CORE-S CAPTURE — fixed and live from the 15:18Z wave**:
  the map runner previously recorded only walls (the ledger's
  `map cut: core-s 0` gap, found by this window's census). The fix:
  `TIMEFORMAT` written directly to `map.cpu` with the command's own
  redirect INSIDE the timed subshell (in bash 5.2 the time report
  follows the timed pipeline's own redirections — the first cut of
  the fix lost the report into the reaped map.log; the quirk is
  receipted here so no successor re-finds it). The product's raw
  stdout/stderr stream is the reaped class (kept only on failure for
  the fail-closed diagnosis).
- **THE CHUNK-RECEIPT TRUNCATION FIX**: `: > "$RCP"` at the call
  step would have wiped the poison-naming lines written at the gate
  — removed (the receipt accumulates: partial-resume notes, poison
  names, then the call's own lines).

## 4. The restated 16-core ledger — MEASURED AT THE WINDOW-2
  CHECKPOINT's n (the honest partials, every n labeled; nothing
  extrapolated)

**The per-step core-s sums at the checkpoint** (the ledger miner,
  `ledger/sums-window2.txt` / `ledger/ledger-window2.tsv`):

| step | core-s measured | n (legs) | per-territory | label |
|---|---:|---:|---:|---|
| substrate chain | 34,464 | 250 | **~138** | MEASURED (extract 2,387 + syng 1,087 + beds 11,862 + export 13,178 + walkreq 18 + walks 602 + maps 5,331; THREADS=4, the campaign's own recipe; walls 13,249s) |
| panel-prep (partition + precompute asset) | 25,323 | 129 | **~196** | MEASURED — **the panel-prep LOAD POLE: the global panel load (~200s) is paid PER TERRITORY; 5,844 × 196 ≈ 1.15e6 core-s — the named lever: a batched panel-prep (the fleet-call form's own amortization class)** |
| map cut (the map leg) | **[PENDING — capture live from the 15:18Z wave]** | 116 done | walls 587–691s at 16-wide (~595s mean over 109 mapped-uncalled) | MEASURED walls; core-s lands per-leg in map.cpu from the fixed runner (§3's window-2 receipts) |
| call chunk (restricted cut + calls, COW) | 18,123 | 7 | **~2,589** | MEASURED at width-7; the width-32 amortization improves the load term ~4.5x (ARITHMETIC, labeled) |

**The honest 16-core projection from measurement** (core-s/16, the
  per-sample all-in at the measured rates, 5,844 territories):
  substrate 138 + panel-prep 196 + map (core-s pending; the wall
  595s at RAYON=32 with a ~260s serial panel load per leg bounds
  the per-leg core-s in the ~6,000–13,000 class — REPLACED by
  measurement as the waves land) + call ~1,500-2,600 — the
  all-in per-territory is **in the same ~4,000-16,000 core-s class
  as window 1's estimate, now with the map term about to sharpen**.
  **[PENDING]** the campaign's own completion receipts sharpen every
  term — especially the tail (the guard census is already receipting
  the 67-70GB memory-tail cut class). **The comparison to the
  pre-campaign estimate stands as window 1 stated it** (the
  stratified 6.38e6 core-s call-phase estimate, ~1,098 core-s per
  territory call-only at 16 cores) and sharpens at the drain, when
  the map core-s and the width-32 call chunks land.

**What the landed levers bought at genome scale, measured so far**:
- the far-bin prune + restricted census profile: ENGAGED BY DEFAULT
  in the product's own path (the model label, the restricted basis
  auto-derivation) — zero flags, as ordered;
- the bounded class scan: ENGAGED per-locus with the derived
  slack/underflow machinery in the emitted certificates;
- the COW fleet form: measured (child private 0.00GB, Pss 16.14GB) —
  the panel load amortization is real;
- the precompute slices 1+2: the assets BUILT and content-keyed
  (once-price receipted ~250 core-s/territory incl. the load);
- the map step: byte-identity reproduced at all three gate levels.

## 5. The exact resume point

- **The driver runs detached** (v4 + the window-2 gate set, `nohup
  env MAP_W=16 PP_W=8 PREP_W=12 CALL_W=2 VALVE=160 bash
  scripts/end-to-end/ee-driver.sh`, pid in the last `driver start`
  line of `logs/driver.log`): the four pools (P substrate / PP
  panel-prep / M map / C call chunks) over the full 5,844-territory
  universe in 183 width-32 chunks; **the map pool is WIDENED to 16
  lanes** (receipted in §3; the box holds: memavail floor ~400G
  mid-wave, guard 64GiB per worker unchanged — the guard is an
  oracle, not an obstacle); the valve holds prep at 160 pending-cut;
  the budget guard 150G. A successor: check `logs/driver.log` (the
  census line per iteration), run `scripts/end-to-end/ee-ledger.sh`
  for the ledger, and RESUME — never restart (the markers govern; if
  the driver died, relaunch it — it is idempotent; the runner scripts
  may be edited BETWEEN iterations only — bash re-reads a runner per
  leg, and a mid-rewrite read is the window-2 syntax-error incident,
  receipted in §3).
- **Watch**: the map pole pace (~87/h at 16 lanes — widen further
  only with a drain-and-restart, and only against the memavail
  floor); the guard census (the 67-70GB memory-tail cut class —
  count the events honestly, never raise the guard); df (the
  heartbeat is in the driver log lines; the window-2 reap points
  self-drain the workdir);
  **the first chunks now dispatch with named poison members** (the
  partial-chunk gate) — read `receipts/call-chunk-*.txt` for the
  poison names per chunk.
- **The known un-run remainder**: the campaign's chunks 0..182
  (chunk-early's 7 are marked called); the fail-closed census at
  this checkpoint names **10 substrate-terminal poison members**
  (8 extract off-by-one + 2 syng SIGSEGV, §3) — each is excluded
  and named at its chunk's call per the committed doctrine; the
  2-strike rules convert any further deterministic prep/panelprep
  failures into named terminals automatically.
