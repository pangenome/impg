# THE WHOLE-GENOME RUNTIME DATASHEET — THE CLOSE (STOPPED BY OWNER AT 37.4% TERMINAL)

> **THE OWNER'S QUESTION, VERBATIM:** "Do we want to do a whole genome run on
> HG002 to verify the run time as you've been arguing?" — **THE ANSWER IS THE
> RUN, AND THIS IS ITS RECEIPT.**
>
> **THE CLOSE, STATED AT THE TOP:** the campaign was **STOPPED BY OWNER**
> (their words: *"ok its stuck... we should stop them"*) at **2,187/5,844
> accepted territories terminal (37.4%)** — the two honest deep-lane tails
> (T02866/T04066) then landed under their own guards, final terminal census
> **2,189/5,844 (37.5%)** — a **partial-drain close, honestly
> receipted** — the close list fired by owner order, not by all-terminal
> completion. **THE VERIFICATION WAS NOT STUCK:** it was the **constrained cuts
> pole** (width-4-direct, C=1, ~29–35 territories/h ≈ **~4–5 days remaining**
> at the current machinery); the owner chose to stop rather than burn days of
> box time, **because the campaign's actual questions were ALREADY ANSWERED by
> the receipts** — every claim below is measured at genome scale. The masters
> were killed by exact pid (prep 411760, cuts 3838087, wave 2810146, plus the
> 22-day-old impg-reference zombie 1846398/1846402); the three honest tails
> (batch-127's cut, the two deep call lanes) were left to resolve under their
> own guards and **all landed honestly** (zero hand kills at the close).
>
> **THE SCOPE, STATED AT THE TOP:** this run uses the **SYNG-NATIVE EXPORT
> MAPS** (the c-rider's `locality-maps-export` at K=63) — the answers produced
> are **EXPORT-PROFILE** (the k-sweep receipt: degraded vs the pggb pair-work
> profile, 62/121 PAIR-CHANGED) — **THIS IS THE RUNTIME RECEIPT, NOT THE
> ACCURACY RECEIPT.** The accuracy-grade genome remains gated on the
> panel-maps build (the owner's standing decision). No genome-wide `--truth-qv`
> (truth assessment-side only; the runtime receipt needs no truth).

## THE FINAL CENSUS (the closing numbers; receipt: poll-once.sh 12:41:39Z + rollup-final.json)

| counter | value |
|---|---|
| split windows | **9,891** (prep/split-done.txt; the query stage's own census) |
| accepted | **5,844 (59.1%)** (prep/accepted.txt) / rejected 4,047 = 3,252 transitive-expansion + 795 ambiguous-axis-bridge (prep/rejected.tsv) |
| substrate done | **2,453** at the stop (the 12:15:02Z poll; the close reap then dropped the producer markers on the reaped never-called set per the idempotence doctrine — **2,210 markers remain**, poll-once.sh 12:57:33Z) |
| cut.done (markers) | **2,228** + 4 work-complete-unmarked (batch-127's cohort — see §2.3) |
| call terminal | **2,189 / 5,844 (37.5%; 2,187 = 37.4% at the stop order)** = callok **1,552** + fail-closed **616** + cut.exit **16** + poison **5** |
| batch census post-swap | **127/127 with markers ALL exit=0, ZERO guard kills since the 18:52:42Z width-4 swap** (cuts/batch-*/batch.exit by file-mtime cutoff epoch 1792176762); +batch-127 work-complete, wrapper markers absent |
| disk at close | df 592–616G free through the run (heartbeat min across the campaign **286G**; heartbeat stopped at the close, 3,353 lines); **the close reap then freed +128.5G apparent — df 588G → 712G avail** (prep 113.3G → 45.3G, map 60.5G → ~0 apparent); df 706–712G at the final polls (12:46–12:57Z); rvsize ~184G → ~59G apparent |

## THE CLAIMS OF RECORD vs THE MEASUREMENTS (every refutation at genome scale)

| component | claimed | measured at genome scale | verdict |
|---|---|---|---|
| the projection (once per sample) | 20.79min | **20.79min CONFIRMED** — the committed region-batching run, identity-gated at every level; the v2 T=32 post-load rate **709,814 reads/s** (the ~710k plateau; the single-process saturation curve plateaus ~500–600k reads/s from T=16 up); and the decisive 12-lane leg (12 private indexes, disjoint record-aligned splits, union multiset identity-gated): per-lane 338–453k reads/s, **AGGREGATE ~4.74M reads/s = 6.7× the best single-process point** — the box-level scaling receipt that grounds the ~3min map pole of the optimization path (§5). Receipts: docs/syng-gmem-bwt/writeup-datasheet.md; /home/erikg/hprcv2/proj-bandwidth/datasheet-section8-draft.md §H2 | **CONFIRMED — the one claim that held** |
| the per-territory prep | ~8s/territory | **79s serial smoke** (T00000: extract 7 + syng 2 + beds 48 + export 2 + walks 1 + maps 19); at close, n=2,453 substrate.done: **median chain ~64s, MEAN ~210s** — the beds segdup tail (median 22s, p90 196s, max 9,000s) dominates. The claim counted only syng+routes with the maps presumed committed. Receipts: prep/*/…*.wall (rollup-final.json `prep`); the ~7–10min 256-wide claim is refuted the same way — the mean, not the median, sets an ~8–12h prep wall at POOL_W=64 | **REFUTED, 8–60x** |
| the per-territory map-cut | ~387s serial | **the batched amortized shape**: the width-4 batch loads the global panel ONCE (242.8s for 4 territories ≈ ~61s/territory amortized), scans the shards once (158–190s shared), per-territory census 13–21s; **batch wall median 510s per 4 territories ≈ ~128s/territory calendar** — but only at BATCH_W=4/N_CONCURRENT=1 (the 64GiB guard's constraint), pacing **~29–35 territories/h**. Receipts: cuts/batch-*/batch.{wall,err}, map/*/*/cut-stats.json (2,078 territories, Σ reads_kept 2.66e9) | **REFUTED as a rate; the amortization is real but guard-constrained** |
| the call | ~23s (35x) | **103s at 1 rayon thread** (T00000, the first complete wave item) — the claimed 23s was the **whole-box** wall; 230 lanes cannot each have the whole box. The campaign aggregate (n=2,168 call-terminal, wave-logs + receipts/*/call.wall): **median 240.5s, p10 7s, p90 2,912s, max 11,085s** under 128-lane packing + co-tenancy; peak RSS median 3.6GB but **p90 67.1GB, max 74.2GB** (the memory-tail lanes). gpartition alone: median 376s, max 1,107s | **REFUTED, 4–10x at lane scale** |
| the wave | ~19–21min at 230–256 lanes (EXTRAPOLATED) | **~80min by honest 230-lane arithmetic** (9,891 × ~111s / 230) — the claim ignored per-lane thread starvation; the measured chunked wave (LANES=128, CHUNK_MULT=2) peaked at **94 call-terminals/hour** through a marathon-heavy hour and was **never the pole** — the cuts pole starved it | **REFUTED** |
| the per-sample serial-sum | ~40–42min | **the real whole-genome shape is DAYS, not 40 minutes**: ~58h of campaign calendar reached 37.4% terminal; the remaining 3,655 never-called accepted territories at the measured cuts-pole pace (~29–35/h) ≈ **~4–5 more days**. The claim's error stack: panel-presumed-committed prep + thread starvation + a wave that presumes an infinitely fast cut feed | **REFUTED — the headline verdict** |

## 1. THE PIPELINE (measured stages; every number from the markers; rollup-final.json)

### 1.1 The panel-side one-time precompute (per PANEL, reused by every sample)
- the query stage (`impg query -d 10000`, the committed stage-B recipe): 14
  shard lanes at `-t 8` (~51GB RSS each, ~630GB RAM held for the whole
  stage); MEASURED ~18–21 windows/min under co-tenancy — the genome-wide
  9,891-window stage is a **9.7h wall** (slowest shard 35,036s; Σ shard walls
  434,068s; receipts query/shard-*.wall). Two stage findings, both receipted:
  the **55–58GB/lane RSS plateau** (growth ~1–2GB/91min once plateaued), and
  the **NON-FUNCTIONAL query-lane RSS guard** (a parse error neutered every
  10s poll for the whole 8h stage — a lane over 64GiB would NOT have been
  killed; moot at the plateau; the guard must be fixed before any re-run);
  the query master **died at its own finish line** (a `[: missing ]]` parse
  error at the wait-return) — recovered by replicating its completion action
  verbatim after verifying every shard's own criteria (queries.done touched
  12:53:10Z; receipts query-stage.nohup).
- the per-territory export-profile substrate, FINAL at n=2,453 substrate.done
  (walls carry the campaign's co-tenancy): extract median 10s (max 503s) +
  syng median 4s (ordinary max 159s; the pathological poison-class max
  12,574s) + **beds median 22s, p90 196s, max 9,000s** (the partition
  -w10k -d1k beds; the segdup tail IS the prep cost) + export 3s + walks 2s
  + maps 23s. **THE HONEST DELTA vs the ~8s/territory claim:** the
  per-territory serial prep is **~64s at the median chain and ~210s at the
  mean** — 8–60x the claim. The claim counted only syng+routes with the maps
  presumed committed; the real one-time chain is
  extract+syng+beds+export+walks+maps, dominated by the beds segdup tail.
- **THE HONEST DELTA vs the ~7–10min 256-wide claim:** the beds tail does not
  parallelize away — at POOL_W=64 the accepted set is an ~8–12h wall (the
  mean, not the median, sets it).

### 1.2 The sample-side map step (the cuts from the committed shards)
- the batched cut (the amortized-load shape): 558 batch dispatches receipted
  (cuts/campaign.log): **431 pre-swap** (296 exit=0; **128 exit-137
  batch-packaging guard kills** — the width-8/16 era's doomed walls, never
  per-territory strikes; 7 exit-1 = the chr2 AGC-offset class's strikes) and
  **127 post-swap, ALL exit=0, ZERO guard kills since the swap**
  (batch-127's cohort work-complete, wrapper markers absent — §2.3). Batch
  wall median **510s per 4-territory width-4 batch**; per-territory
  cut-stats (2,078 territories): wall mean ~500s (load 279s mean per
  cut-stats, amortized to ~61s/territory at W=4; scan 158–190s shared; census
  13–21s); **Σ reads_kept 2.66e9**.
- **THE HONEST DELTA vs the ~387s/territory serial claim:** the amortized
  batch cut is ~128s/territory of calendar at W=4 — but the 64GiB guard
  forces BATCH_W=4 N_CONCURRENT=1, and the cuts pole paces the whole genome
  at **~29–35 territories/h**. The cut machinery, not the cut arithmetic, is
  the genome-scale pole.

### 1.3 The streaming call wave
- the lane-packing reality (per-call walls + RSS at the lane width, n=2,168
  call-terminal): **median 240.5s, p10 7s, p90 2,912s, max 11,085s**
  (receipts/*/call.wall); reads fetch ~0–1s (warm cache); **gpartition
  median 376s, max 1,107s** — a stage the claims never priced; peak RSS
  median 3.6GB, p90 67.1GB, max 74.2GB (receipts/*/call.peak_rss_kb,
  call.rss). **35,310 loci called across the 1,552 callok territories**
  (~22.8 loci/territory).
- **THE HONEST DELTA vs the 19–21min wave claim:** the claim's 23s/territory
  call was a whole-box rayon wall; at 1 thread/lane the call is 103s (T00000
  calibration), and 230-lane arithmetic gives **~80min, not 19–21**. The
  measured wave never became the pole: the cuts starved it.

### 1.4 THE END-TO-END WALL vs the ~40–42min serial-sum claim
- **REFUTED at genome scale: the real shape is DAYS.** The campaign ran ~58h
  (2026-10-15 ~02Z → stop 2026-10-17 ~12:15Z) to reach 2,189/5,844 terminal
  (37.4%); the remaining 3,655 never-called accepted territories at the
  measured ~29–35/h cuts-pole pace ≈ **~4–5 more days**. The one CONFIRMED
  stage is the projection (20.79min); every downstream claim carried the
  panel-presumed-committed, thread-starvation, or infinite-cut-feed error.

## 2. THE PRODUCING-PACE RECEIPTS (the campaign's own fixes, measured live)

### 2.1 The width-4-direct swap (the +54% delta receipt)
- BEFORE (width-8-doomed, 16:57:42–18:52:42Z 2026-10-16, 115min): 12 batches
  = **9 guard kills** (doomed walls 447–630s) + 3 honest 8/8; +39 cut.done =
  **20.3/h**. AFTER (width-4-direct, first 38.3 producing minutes): 5 batches
  = **5 HONEST 4/4, ZERO kills**, walls 530/457/427/444/430s (avg 457.6s);
  +20 cut.done = **31.3/h**. **DELTA: +54% producing pace, doomed walls
  eliminated** (0/5 vs 9/12 killed). Honest caveat: the window's last 16.1min
  were the 300G disk floor self-gating the master (df 299G → released at
  305G) — the calendar hour is the disk pole's, not the width's. Receipt:
  cuts/width4-watcher.log (the swap), cut.done marker mtimes (the delta).
- THE CLOSE CENSUS: **127/127 post-swap batch.exit markers ALL exit=0, zero
  guard kills since the swap** (batch-13's width-8 guard-137 predates the
  swap by 2s — campaign history); the post-swap batch RSS band ~48–61GB
  (batch-122's ~60.9GB the widest — the dense-locus cohorts named), in-band
  under the 64GiB guard.

### 2.2 The chunked wave dispatch fix (famine/burst barrier → bounded backlog burst)
- BEFORE: the unchunked wave held cut.done backlogs 5+ hours behind marathon
  lanes (famine bands **1–15 terminals/h**; the hour histogram
  01h=3/02h=44/03h=7/04h=9/05h=106/06h=14/07h=30/08h=52) then burst-cascade
  100–259/h when a snapshot broke — the famine/burst barrier. THE FIX:
  `CHUNK_MULT=2` chunked rolling dispatch in run-wave.sh (CH=2×LANES, the
  valve re-checked between chunks). AFTER (the ride-24 receipt): **94
  call-terminals/hour INDEPENDENT through a marathon-heavy hour** — the
  bounded 15-min shape receipted (76 terminals in the first 15min = the
  bounded backlog burst, then the chunk-1 tail). Receipts:
  HANDOFF ride-24; wave-logs/ (terminal mtimes).

### 2.3 The mid-dispatch valve patch's live proof + the close's one marker caveat
- THE VALVE (the chunked-dispatch prep patch) proved LIVE in
  prep-master.log: `14:56:22Z VALVE (mid-dispatch): 218 pending cut (>= 160)
  - parking the remaining 5019 of 5275 TODO` (and again 20:46:06Z) — the
  prep pole self-pauses mid-dispatch at the 160 bound, exactly as designed.
- THE CLOSE'S ONE MARKER CAVEAT, honestly receipted: **batch-127**
  (cuts/batch-1792239119-127, dispatched 12:11:59Z, the last cut in flight
  at the owner's stop) **completed its work honestly** — `[cut-batch] DONE: 4
  territories in 428s (peak RSS 52,033,076kB, 0 emission-failures)`, all four
  cut-stats.json written (batch.err + map/T04353|T04352|T04289|T04164) — but
  its wrapper subshell died with the owner-ordered master kill, so the
  wrapper-side markers (batch.exit, batch.wall, prep/<T>/cut.done, the
  campaign.log done-line) are **absent**. The cut.done marker census stays
  2,228; the work-complete census is 2,232. Nothing was fabricated.

## 3. THE FAILURE-CLASS CENSUS (counted and named at genome scale; every class terminal)

- **the 7479-class bind-ambiguity family — THE DOMINANT fail-closed class,
  482 of the 616 fail-closed calls (78%)**: `census record has no unique
  verifying key` (the pilot's census record 7479 was the namesake; the
  two-origin-shift bind, the record refusing an ambiguous key — as
  designed). Receipts: wave-logs/failures.txt; prep/*/call.{exit,err}.
- **the deep-lane memory-tail guard-137 class — 134 guard kills total (132
  exit-137 + 2 exit-1 "RSS guard exceeded"), every one the sanctioned
  64GiB guard's own kill, zero hand kills**: the marathon roster named
  (T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915/T03919/T02722/
  T04059/T03411/T04038 + the close's T04066 and T02866), with **T03641 the
  deepest at 11,084s / killed at 71,110,568 kB (71.1GB)**; the close's two
  tails landed honestly: T04066 exit 137 at 6,804s / 67,829,716 kB, T02866
  exit 137 at 7,730s / 67,575,876 kB — the fast-climb tail after a ~1.5h
  plateau. The honest distribution point beside them: **T03208 completed
  honestly at 4,615s** (the watched 37–45GB lane that never tripped).
  Receipts: prep/*/call.{rss,wall,exit.code}; wave-logs/failures.txt.
- **the syng segdup-tail poison class — 5 members, every one poison-first
  then exact-pid kill** (T00197/T01411/T02562/T00024/T03679): the syng step
  CPU-bound at ~23MB RSS with no output, 29min–3h29m burned per member vs
  the ordinary corpus's median 4s/max 159s; each poisoned in cuts/poison.tsv
  BEFORE the kill (the master's TODO filter skips them at both dispatch
  sites). Full census receipted in cuts/poison.tsv.
- **the chr2 AGC-offset class — cut.exit 16, two locus clusters, 3-strike
  terminal**: `locus record GRCh38#0#chr2:… is not identical to the AGC
  interval under either offset convention`, deterministic; chr2:20.75Mb
  (T00341/T00905/T00955/T00724/T00631/T00950/T00791/T00974) and
  chr2:109.25Mb (T01081/T01088/T01135/T01144/T01161/T01221/T01229/T01357),
  terminal after 3 batch-level failures by the designed attempt-counting
  path. Receipts: cuts/batch-1792061780-15/batch.err,
  cuts/batch-1792076665-12/batch.err, prep/*/cut.{exit,attempts}.
- **the expansion-gate fail-closed class — 3,252 split rejections** (rows ≤
  2× distinct paths AND paths ≤ 2000; measured up to 695k rows/14.9k
  paths/187Gbp for ONE window), cumulative **184,355,392 rows / 25,214,742
  paths** refused. Receipts: prep/rejected.tsv, prep/*/split.rejected.
- **the ambiguous-axis-bridge class — 795 split rejections** (>1 same-contig
  axis candidate, fail-closed by the axis row rule). Receipt: prep/rejected.tsv.
- **the axis-drift artifact and the mappy-2.31 bridge fix (the campaign's
  own correction)**: the first split rule (fixed ±250kb slack) wrongly
  rejected every territory past HG002#1's coordinate drift (measured chr1:
  +113kb at 3.3M → +3.8Mb past 150M) — the "assembly gap" ~30% class was
  ~100% drift artifact. THE FIX: the axis bridge (every window aligned to
  the HG002#1#chr1-22 target via **mappy 2.31**, the committed minimap2
  2.31 binary having been deleted from the box; 11,953 primary rows over
  9,963 windows) — acceptance **36% (624/1,736) → 61% (1,188/1,946)**,
  final genome-wide **59.1% (5,844/9,891)**, additions only, zero
  corruptions (all 392 substrate reference-axis.json components verified
  identical to the rebuilt axis-path.txt). Receipts: axis-bridges.paf,
  axis-bridge.log, rv-axis-bridge.py; the split-rebuild ledger.
- the no-evidence class: **0** (every cut produced evidence).
- the batch-packaging guard kills (width-8/16 era, never territory strikes):
  **128 exit-137 batch markers pre-swap, ZERO post-swap** (§2.1).

## 4. THE SUBSTRATE JIT REAPER + THE DISK HIGH-WATER (the streaming design's own receipts)

- **the disk high-water (wave-logs/heartbeat.log, 3,353 lines — the
  daemon receipted and stopped at the close)**: minimum df
  **286G** free across the campaign (the 300G hard floor era; the gate
  self-paused the cuts master at df 299G and released at 305G — the floor
  held, never lowered without the owner); min memavail **21GiB** (the
  14-lane query era); the one early excursion below the original 1,500G
  floor (1,373G, co-tenant) stands recorded. df at close 592–616G.
- **the substrate JIT reaper receipts** (prep-reaper.log + the early-pass
  log prep-reaper-r33.log; the ride accounting across **79 idempotent
  passes**): **~207.2G apparent freed (~204.0G class-attributed)** — the
  log-level sums: 137,058M apparent over the first 3 passes (the 129G
  first pass: the split-rejected dead extracts) + 75,317M over the 78
  cadence passes; class A ~93.6G (rejected territories' dead locus.fa +
  abandoned partial syng indexes), class B ~20.6G (cut-awaiting-call
  consumed syng index + export/walks), class C ~97.6G (built-ahead pending
  consumed intermediates). Every cut/call input and every split ledger kept
  throughout; the exclusion proof (the six-file syng sidecar rule) ran on
  every pass.
- **the ~22s on-demand regen measurement (the T03194 full-chain re-entry
  receipt)**: an ordinary territory evicted fully (beds kept) rebuilt by the
  idempotent chain (extract+syng+export+walks+maps) and RE-ENTERED the
  normal cut/call flow — the regen rate (~22–30s/territory) beats the
  ~31/h cut rate, the JIT-prep rule's foundation. (The one marker-gap found
  and fixed en route: a kept walkreq.done for a reaped walk-requests.tsv
  broke the regen — the marker goes with the output.)
- **the class-C bulge/equilibrium receipts**: the consumed-intermediate
  reaps bulge when the prep pole outpaces the cuts (the pass-68/72/76/79
  bursts, 0.7–1.5G a pass; the rides-55/56 clusters) — the reaper and the
  valve hold the built-ahead window at the lean buffer (~175–185 pending,
  bound 200), the measured steady state of the JIT-prep equilibrium.
- **THE INVENTORY FINDING (the design verdict)**: the per-territory
  materialization **should not exist**. The measured prep/ inventory
  (228G at its ride-31 census; ~111G at close) is locus.fa + a private syng
  index + export graphs/walks + locality-graphs **per territory** — a
  per-territory re-derivation of what the global panel already knows. The
  cut already runs `--global-panel`; the call still materializes. **The
  named fix is the global-panel call mode** — call against the shared
  panel index directly, with the territory's window set as the only
  per-territory input.
- **THE CLOSE-OUT REAP (the stopped-by-owner close's own reap, 12:45–12:47Z,
  receipt close-reap.log + rv-close-reap.sh)**: **freed +128.5G apparent — df
  588G → 712G avail** (CLOSE-A: 2,129 files/dirs, 69,321M regenerable
  substrate on the 3,676 never-called accepted territories + the idempotent
  residue sweep, producer markers dropped per the JIT doctrine; CLOSE-B:
  270 files, 60,596M cut-evidence bulk — cache/census/reads, the deterministic
  byte-identity-gated cut — on 90 territories, map/<T>/<T>/cut-stats.json
  KEPT as the cut receipt). **VERIFIED SURGICAL at the close** (the
  spot-checks): the input BAM and hg002-projection/full/ untouched outside
  this tree; every definitional bed/ledger intact (locus.refined.bed,
  locus.members.bed, locus.bed, split.info.json and split.rejected on ALL
  9,891 territory dirs; accepted.txt/split-done.txt/rejected.tsv);
  receipts/ complete (2,168 dirs: 1,552 calls.jsonl + sequences + 2,014
  cut-stats.json); the campaign logs, batch receipts (558 batch.exit +
  batch-127's work-complete receipt), call walls/RSS receipts, poison.tsv
  (5 members, full prose), the handoff, this datasheet and the query
  results all intact; only regenerable substrate gone (locus.fa, the syng
  slices, export-graphs/walks/maps, locality-graphs — 0 of each remain in
  prep/, verified); the 5 poison territories' syng receipts kept verbatim.
  The two leftover rider-class daemons (heartbeat ×2, the REARM df-waiter)
  were stopped at the close — logging/waiter processes only, zero work
  processes touched.

## 5. THE VERDICT (measured, not extrapolated)

- **The one confirmed claim:** the projection — 20.79min full-sample,
  identity-gated, the ~710k reads/s region-batched plateau (709,814 reads/s
  at T=32) and the 12-lane leg's ~4.74M reads/s aggregate (6.7×) — the
  measured basis of the optimization path's ~3min map term.
- **Every downstream claim refuted at genome scale** (the table above): prep
  8–60x, cut guard-constrained to ~30/h, call 4–10x at lane scale, the wave
  ~4x by honest arithmetic, and the whole-genome serial-sum **off by two
  orders of magnitude** — the real shape is **days, not 40 minutes**.
- **THE CAMPAIGN WAS STOPPED BY OWNER at 37.4% terminal — stopped, not
  stuck**: the constrained cuts pole (width-4-direct, C=1) paced ~29–35
  territories/h ≈ ~4–5 days remaining; the owner chose to stop rather than
  burn the days, the questions already answered. Final census §THE FINAL
  CENSUS above; reap receipt in close-reap.log.
- **THE OPTIMIZATION PATH — THE OWNER'S 20-MINUTE TARGET, GREEN-LIT AT THE
  CLOSE (the measured budget table; every term sized by this campaign's
  census, no term extrapolated beyond a measured receipt):**

  | stage | budget | the measured basis |
  |---|---|---|
  | the map (projection → per-territory cuts) | **~3min** | the projection is 20.79min serial at the ~710k single-process plateau, but the 12-lane leg measured **~4.74M reads/s aggregate (6.7×)** — the same lane scaling the optimization path assumes |
  | the prep | **ZERO** | the global-panel call mode (the owner's design, under implementation now) kills the per-territory materialization — the inventory finding: the substrate is a re-derivation of what the shared panel already knows |
  | the cuts + calls | **~10–15min** | IF the fixed loads die (the 242.8s/batch global-panel load amortized to a shared resident index; the gpartition median 376s and call median 240s setup terms batched away) AND the lanes scale like the projection's measured 6.7× at 12 lanes |
  | the marathon tail | **the one research risk** | the cross_lcp_merge call-core work — this campaign's census sizes it: the deep-lane memory-tail class 134 guard kills, p90 call wall 2,912s, deepest T03641 11,084s/71.1GB; if the tail does not shrink, it alone can eat the budget |

  The mechanism, in priority order (the work order the budget terms map
  to): (1) **the global-panel call mode first** (the prep-ZERO term — the
  call reads only locus.fa/axis/locality-graphs today, all regenerable
  re-derivations of the shared panel); (2) **shared-index cut lanes** (the
  242.8s global-panel load is paid once per batch today; a shared
  resident index across concurrent lanes removes the load term and
  re-opens N_CONCURRENT); (3) **batch-amortized call** (amortize the
  per-call setup as the cut's batched load already does — the
  gpartition/call setup dominates at median 376s+240s per territory);
  (4) **the cross_lcp_merge call-core work** (the marathon-tail term — the
  deep-lane memory tail and the 103s 1-thread wall both live there).

## APPENDIX: the co-tenant/noise statement
The query stage (14 lanes), the prep pool, the cut batches and the call wave
run as a pipeline — the per-stage walls carry the pipeline's own contention
(the committed record's co-tenant ±20% class). The end-to-end wall is the
pipeline's wall; a pure per-sample wave on a prepared panel would run clean
(the per-lane calibration receipt states both). This close is a
partial-drain close by owner order: every number above is measured on the
territories that ran; the never-called 3,655 are counted, named, and their
substrate reaped as regenerable (close-reap.log) — nothing about them is
extrapolated.
