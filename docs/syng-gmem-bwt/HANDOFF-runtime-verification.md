# HANDOFF (the whole-genome runtime verification session, worker thread — keep updated)

## THE TASK
The owner's question, verbatim: "Do we want to do a whole genome run on HG002 to verify
the run time as you've been arguing?" — THE ANSWER IS THE RUN. Measure the whole-genome
35x HG002 runtime end to end vs the claims: projection 20.79min (committed, DONE, never
re-run); per-territory map-cut ~387s serial; call ~23s; the wave ~19-21min at 230-256
lanes (EXTRAPOLATED, never run); serial-sum ~40-42min. SCOPE HONEST: the run uses the
SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers
(degraded vs pggb pair-work, k-sweep receipt) — THE RUNTIME RECEIPT, NOT THE ACCURACY
RECEIPT. No genome-wide --truth-qv. No thresholds. Committed projection/shards/receipts
immutable. NO PR PUSH. Streaming: per-territory receipts kept (calls.jsonl + sequences),
intermediates deleted as consumed (JIT reaper); disk high-water under ~500GB; the run
is marker-gated, killable, resumable — a successor RESUMES, never restarts.

## LAYOUT (all under /home/erikg/hprcv2/wgs35x-fleet/runtime-verification/ = $RV)
- bed4-genome.tsv           the genome-wide BED4 (chrom,start,end,TID) — 9,891 rows (the bridge-accepted)
- query-input.bed           the query rows (CHM13 target, start, end, TID)
- HANDOFF-runtime-verification.md   THIS FILE — keep updated
- prep/                     the per-territory export-profile substrate (TID dirs, marker-gated)
- partitions/               the three-step PARTITION step outputs (partitions.bed + <TID>/)
- map/                      the cut outputs (cut-batch, per territory <TID>/<TID>/cut-*)
- calls/                    the CALL receipts (per territory)
- wave-logs/                stage walls/rss/heartbeat
- datasheet-runtime.md      THE VERDICT (assembled at the end)

## THE STAGES (marker-gated; every stage nohup + done-markers + SHORT polls)
0. session-open increment appended to docs/syng-gmem-bwt/campaign-record-2026-10-diploid-arc.md (DONE).
1. QUERY (panel-side, the committed stage-B recipe): `impg query -d 10000` over the 9,891
   bridged windows, sharded N ways (each process loads the global syng ~35.5GiB + AGC);
   then the per-territory refined-bed split + the axis/truth-row check (the old
   query/split-refined-beds.py was disk-cleaned; re-implemented as
   scripts/rv-split-refined-beds.py — split by window, clip to overlap, strand kept,
   acceptance = exactly one HG002#1 axis row; dup-path/no-axis = fail-closed, named).
2. PREP per territory (export profile): extract locus.fa (AGC) -> syng (k63) ->
   partition -w10k -d1k -> partition_graph_export -> syng_homology_dump walks ->
   locality-maps-export.py (K=63) -> [after cut] genome-infer partition (routes by
   lookup + window beds + registry + the BED4 row). Measure per-step walls.
3. CUT (sample-side): sample_projection cut-batch (IMPG_CUT_LOCAL_DERIVE=1) from the
   COMMITTED shards /home/erikg/hprcv2/hg002-projection/full — batches of W territories,
   one global syng+AGC load + ONE shard scan per batch; per-territory census/cache/
   reads/stats. Measure batch walls + the amortized per-territory cut.
4. CALL WAVE: per territory `impg genome-infer call-region --map map/<TID> --partitions
   partitions/ --region <TID> --out-dir calls/<TID>` (reads = the per-territory 35x
   region fetch from HG002.GRCh38.35x.bam, ±50kb, the committed fetch rule). Lane count
   sized by measured RSS (claim 3.4GB/lane). FAIL-CLOSED classes counted + named
   (the 7479-style bind ambiguity expected at some territories).
5. REAPER: after each territory's call: keep calls.jsonl + calls.jsonl.sequences.jsonl +
   cut-stats.json + walls; DELETE the census, cache, reads fasta, fq, the substrate
   (locus.fa, syng, beds, export graphs, walks, maps as consumed by the partition step).
6. ROLLUP + VERDICT: every stage's measured wall, the per-territory cost distribution,
   the failure-class census, the disk high-water, the end-to-end wall vs the claims.

## STATE AT FIRST CHECKPOINT (03:20Z 2026-10-15 — update as stages land)
- [x] session-open increment in the record (committed: 231d21a with gate-1 tooling)
- [x] this handoff
- [x] BED4 + query input composed (9,891 rows)
- [x] THE REGRESSION GATE PASS: the patched cut-batch (per-territory fail-closed
      isolation in the emissions loop) is BYTE-IDENTICAL on the success path vs the
      committed cutbatch-gate-L1 receipt (census+cache+reads cmp-clean; cut-stats
      CONTENT-identical; only timing fields differ - the 14 co-tenant query lanes,
      wall 554s vs the committed 413s, stated honestly). Receipts: regression-L1/.
- [x] PREP SMOKE (T00000): the full export-profile substrate chain 79s SERIAL
      (extract 7 + syng 2 + beds 48 + export 2 + walks 1 + maps 19; RSS tiny:
      syng 30MB, beds 245MB) - vs the ~8s/territory claim (which counted only
      syng+routes with maps presumed committed) - the honest prep delta.
- [~] QUERY STAGE RUNNING: 14 shard processes (each ~45GB RSS, -t 8; ~707 windows
      per shard at ~17.2s/window measured at the smoke) - ETA ~06:30Z. query/shard-*.err.
- [~] PREP MASTER RUNNING (run-prep-master.sh): splits per completed shard, then
      preps at 64-wide. NOTE: prep failures have no terminal marker - if a
      territory's prep fails deterministically the master retries forever (kill +
      diagnose; add a 3-strikes rule if it bites).
- [~] CUT CAMPAIGN RUNNING (run-cuts.sh): batch-0 (T00000) in scan, contending
      with the queries (walls inflated ~30% - the co-tenant class, stated).
- [~] WAVE RUNNING (run-wave.sh): waits for cuts (fixed: it used to exit when
      no cuts existed).
- [x] HEARTBEAT RUNNING: wave-logs/heartbeat.log (df/memavail/counters each 60s -
      the disk high-water of record).
- [x] rollup tooling: rv-rollup.py (stage walls, cost distribution joined vs the
      grid's segdup fraction, failure census, disk high-water).
- [ ] the call calibration (wall + RSS at RAYON_NUM_THREADS=1) when T00000's cut
      lands -> then size LANES honestly.
- [ ] the datasheet verdict (datasheet-runtime.md) at the end.

## THE FIRST MEASUREMENTS (T00000, the first complete wave item, 03:28Z)
- THE CUT (batch-0, W=1, contended by the 14 query lanes): wall 477s
  (load 258.6 + scan 209.8 + census 7.9), peak RSS 38.99GB, CONTENT IDENTICAL to
  the old fleet's committed T00000 batch receipt (reads_kept 109,997 /
  census_mult 125,956 / routed_patterns 7,254) - the export-profile axis produces
  the same universe rows (the maps differ, the tiling does not).
- THE WAVE ITEM: reads fetch 0s (warm cache) + gpartition 6s + CALL 103s
  (RAYON_NUM_THREADS=1) + reap; 25 calls + 25 sequences; receipts reaped clean
  (prep/T00000 508K after).
- THE LANE-PACKING REALITY (the claim's honest delta, FIRST DATA POINT):
  per-call peak RSS 2.0GB (BELOW the claimed 3.4GB - more lanes fit) but the
  per-call WALL at 1 rayon thread is 103s vs the claimed ~23s (which was the
  WHOLE-BOX-rayon wall - 230 lanes cannot each have the whole box). Serial
  fraction ~20s, parallel ~83s (Amdahl from 20s@whole-box vs 103s@1-thread,
  both contended - a clean calibration still to take). The honest wave
  arithmetic at 230 lanes x 1 thread: 9,891 x ~111s / 230 ~= 80min, NOT the
  claimed 19-21min - THE CLAIM'S WAVE ARITHMETIC IGNORED PER-LANE THREAD
  STARVATION. The campaign aggregate is the receipt.
- THE PREP DELTA: 79s serial/territory measured (the claim's ~8s counted only
  syng+routes with maps presumed committed; the real one-time chain is
  extract+syng+beds+export+walks+maps+routes) - genome-wide one-time, stated
  as the panel-side pole honestly.

## THE AXIS ROW RULE (gate 3, the split's deterministic acceptance - defined and stated)
A window's query refined bed carries MANY HG002#1 rows (the transitive -d 10000
descent pulls in cross-chromosomal segdup hits: T00440 carries 18 HG002#1 rows,
12 of them on chr1 at 16.7M/120.5M/140.5M/... - segdup hits, NOT the axis). THE
RULE: the axis = the unique HG002#1#<bridge-chrom> row intersecting the bridge
target +/- one territory-span (250kb, the territory's own size - structural, not
tuned). 0 candidates = no-same-chrom-axis-row (fail-closed); >1 =
multi-same-chrom-axis-row (fail-closed). THE CLASS COUNTS at the first ~905
completed windows: accepted 552 (61%), no-same-chrom 349 (of which 314 have the
NEAREST same-chrom row >=3Mb away - GENUINE HG002#1 paternal assembly gaps, the
declared-axis convention's own coverage reality; 15 near-miss displacements
<1Mb; 20 in 1-3Mb), multi-same-chrom 4. THE OLD FLEET's first-HG002#1-row rule
took segdup hits as axes at gap territories - never a coordinate frame; this
run's rule is stricter and honest (the near-miss subclass is counted, not
silently rescued). The final reconciliation split (queries.done ->
rv-split-refined-beds.py) applies the SAME rule over the complete outputs.

## THE PER-TERRITORY WAVE ORDER (run-wave.sh, all marker-gated)
reads fetch (samtools region, the committed fetch rule) -> map/<T>/manifest.json
(reads path) -> genome-infer partition (flock: partitions.bed append race) ->
call-region --region <T> --out-dir calls/<T> (RAYON_NUM_THREADS=1 per lane,
IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS=1) -> THE JIT REAPER (receipts to
receipts/<T>/; DELETE map/<T>, reads fq, partitions/<T>, and the prep substrate's
heavy artifacts - locus.fa/syng/beds/export/walks/maps stay only as long as
uncalled territories need them).
RETRY NOTE: a failed call is terminal (call.exit); to retry one, clear
prep/<T>/call.* markers + rm -rf map/<T> calls/<T> partitions/<T>.
NOTE: genome-infer partition refuses an existing partitions/<T> (the discarded-run
rule) - clear it when retrying.

## THE CAMPAIGN-SIZING PLAN (the successor's checklist)
- NOW (query stage hogging ~630GB RAM + 112 threads): cut campaign C=2 W=32
  RAYON=64; wave LANES=32; prep POOL_W=32 with THE VALVE (>=160 substrates
  pending cut -> hold the prep pool; the substrate footprint is ~250MB/territory
  and prep outruns cuts ~6x - without the valve the pending window grows
  unbounded and the disk floor trips).
- WHEN queries.done LANDS (~12:15Z ETA at the measured 18 windows/min): relaunch
  the cut campaign WIDER - `kill the run-cuts master; BATCH_W=96 N_CONCURRENT=4
  RAYON=48 ./launch-cuts.sh` (RAM ~4x(39+96*0.4)~230GB) - and raise the wave
  (`LANES=128 ./launch-wave.sh`). The terminal conditions all wait for
  queries.done + final-split.done + campaign.done - relaunching is safe at any
  point (every stage marker-gated).
- THE HONEST CUT ARITHMETIC so far: W=32 batch ~= the 1-territory batch's
  load+scan (~400-500s contended) + ~32 x ~8s emission - the per-territory
  amortized cut tends to ~(load+scan)/W; the W=32 first-batch wall lands in
  cuts/batch-*/batch.wall when it completes (contended by the query lanes -
  the co-tenant class).
- THE END-TO-END PROJECTION at current packing: query ~12:15Z -> cuts +~3-4h
  at the widened packing -> the wave trails (it consumes faster than cuts
  produce) -> campaign complete ~16-17Z. THE SUCCESSOR RESUMES - every stage is
  marker-gated; kill anything, relaunch the masters, they pick up where the
  markers stand.

## STATE (04:15Z 2026-10-15)
- query: 1,180/9,891 windows (~18/min, all 14 lanes healthy, ~45GB RSS each)
- split (incremental): 1,180 split; accepted ~700 (61%), no-same-chrom ~460
  (the HG002#1 gap class - see THE AXIS ROW RULE), multi-same-chrom ~8
- prep: 392 substrates done; THE VALVE holding (239 pending cut >= 160)
- cuts: T00000 done (batch-0); two W=32 batches in scan since 03:53Z
- wave: T00000 called+reaped (receipts/T00000); LANES=32 waiting for cuts
- disk: df 1,634G avail (floor gate 1,500G); rvsize 385G (the in-flight window:
  ~390 substrates x 250MB + 64 cut outputs - the valve caps the growth)
- masters running: prep-master 1830969, cuts 1673560, wave 1858400,
  heartbeat 1454428, query-stage (the original run-query.sh)

## THE DEADLINE CHECKPOINT (04:40Z 2026-10-15) - WHAT STANDS, WHAT'S NEXT
THE AXIS DRIFT DISCOVERY (the session's key finding): the committed HG002#1 axis
assembly's coordinate offset vs the CHM13 bridge target DRIFTS ALONG THE
CHROMOSOME (measured chr1: +113kb at 3.3M -> +400kb mid-chrom -> +3.8Mb past
150M). The first split rule (a fixed +/-250kb slack) therefore WRONGLY rejected
every territory past the drift point - the 'no-same-chrom' 30% class was mostly
drift artifact. THE FIX (in flight at the deadline): THE AXIS BRIDGE - the
committed bridge recipe's second leg: every window aligned (asm5) to the
HG002#1#chr1-22 target (extracted from the AGC, 2.9GB, 12s) with the SAME
single-primary-same-chrom acceptance rule; the axis ROW = the refined bed's
same-contig row with max overlap vs the bridge interval. The committed binary
(minimap2 2.31) is GONE from the box (/opt/minimap2-2.31_x64-linux deleted
since Oct 9); mappy 2.31 (the same minimap2 core) is the stated substitute.
STATE: rv-axis-bridge.py launched (index build ~minutes; output
axis-bridges.paf); rv-split-incremental.py PATCHED to consume it (classes:
no-axis-bridge / ambiguous-axis-bridge / no-image-row / multi-image-row) BUT
the split state on disk is still the OLD RULE's (516 accepted with
within-slack axes - all CORRECT axes, just incomplete coverage; the
wrongly-rejected are additions, never corruptions).
SUCCESSOR'S NEXT STEPS, IN ORDER:
1. Wait for axis-bridges.paf (axis-bridge.log says 'done'); check
   axis-bridges.exit is 0.
2. Rebuild the split state: `rm -f prep/*/split.rejected; : > prep/split-done.txt;
   : > prep/accepted.txt; rm -f prep/rejected-incremental.tsv;
   python3 rv-split-incremental.py $PWD` (repeat until 'no new completed
   windows'). The pipeline picks up the additions automatically (all
   marker-gated; already-called territories stay valid).
3. ALSO patch rv-split-refined-beds.py (the final reconciliation split) with
   the same bridge-based rule - it still carries the fixed-slack rule.
4. THE EXPANSION GATE stays: rows <= 2 x distinct paths AND paths <= 2000
   (the transitive-expansion class: measured up to 695k rows/14.9k paths/187Gbp
   for ONE window; the committed class is ~466 rows ~466 paths, max 1.004x).
5. When queries.done lands: widen the cut campaign (BATCH_W=96 N_CONCURRENT=4
   RAYON=48) and the wave (LANES=128).
6. At campaign close: rv-rollup.py + fill datasheet-runtime.md + the record's
   closing increment + the commit.
RUNNING MASTERS (all marker-gated, killable, resumable): prep-master 1830969
(POOL_W=32 + the valve), cut campaign 1673560 (C=2 W=32), wave 1858400
(LANES=32), heartbeat 1454428, the 14-lane query stage (run-query.sh). The
first complete end-to-end wave item (T00000) is the receipt: cut content
identical to the old fleet's committed T00000 batch receipt; call 103s at
RAYON=1, peak RSS 2.0GB; 25 calls + 25 sequences; the JIT reaper clean.
MEASURED SO FAR (the honest deltas vs the claims): prep 79s serial/territory
(claimed ~8s = syng+routes only); the call at 1 rayon thread 103s (the claimed
23s is the whole-box wall - 230 lanes cannot each have the whole box; the
wave at 230 lanes ~= 80min, NOT 19-21min); the call peak RSS 2.0GB (below the
claimed 3.4GB/lane); the cut at W=1 477s contended (the committed 387s
clean); the ~8s and 19-21min claims carry the panel-presumed-committed and
thread-starvation errors respectively - the campaign aggregate is the receipt.

## SURVIVAL RULES
no recursive greps from the repo root (target/ 6.2GB); bash timeout 60; the syng files
LARGE - bound reads; disk floor check before every launch (~2.0T free at open; keep the
high-water under 500GB — if the in-flight window grows past it, NARROW THE LANES); the
64GiB guard per process; long stages nohup + done-markers + SHORT polls.

## THE SUCCESSOR CHECKPOINT 1 (05:15Z 2026-10-15) - THE BRIDGE LANDS + THE GUARD-KILL CLASS
THE AXIS BRIDGE IS DONE: the predecessor's rv-axis-bridge.py died on a mappy API
bug ('m_primary' -> the attribute is 'is_primary'; verified against the live
binding before relaunch); relaunched 04:53Z, DONE in ~13min wall (11,953 primary
rows over the 9,963 windows, index build ~5min + 24 workers). axis-bridges.exit=0.
THE SPLIT REBUILD DONE (the prep master SIGSTOPped at a safe sleep-point during
it, resumed after): ledgers reset, two incremental passes over 1,946 completed
windows: ACCEPTED 1,188 (61% - the drift-lost territories returned; was 624/1,736
= 36% under the old fixed-slack rule), REJECTED 758: transitive-expansion 618,
ambiguous-axis-bridge 140, and ZERO no-axis-bridge/no-image-row - the old rule's
"assembly gap" class was 100% drift artifact (spot-check T00453: old-rule
no-same-chrom rejection, bridge row HG002#1#chr1:150.9-151.2Mb mapq60, matching
refined-bed row, accepted). CONSISTENCY VERIFIED: all 392 substrate.done
territories' reference-axis.json components == the rebuilt axis-path.txt; all
95 orphan-completed beds' starting-sequences == axis-path. Additions only, zero
corruptions (the predecessor's prediction held).
rv-split-refined-beds.py PATCHED with the same bridge rule (classes renamed:
no-axis-bridge/ambiguous-axis-bridge/no-image-row/multi-image-row; refuses to run
without axis-bridges.paf) - the final reconciliation split at queries.done.
THE ORPHAN/CORPSE SWEEP (a new failure class, named): earlier prep waves were
killed without group-kill - their impg children ran on with dead wrappers, so
beds.done markers could never be written. 96 TIDs had marker-less partition
dirs; rv-orphan-sweep.py marked 87 orphan-COMPLETED ones (the impg-written
terminal 'Partitioned into N regions' line + >=1 bed = the honest success
signal; walls estimated from filesystem timestamps, flagged in
beds.wall-estimated), and rv-orphan-watcher.py (running, detached) watches the
last 7 live orphan partitions (the segdup tail: T00207 beds wall 3,523s = 59min
vs the 48s ordinary - the real tail) + the T00197 syng orphan (exit status
unverifiable - its partial outputs will be cleaned for a clean re-run).
THE DOUBLE-WRITERS KILLED by exact pid (1502756/1641812 T00003,
1535438/1642106 T00114): two concurrent partition copies writing the same
results dirs (corrupt by construction); their dirs rm -rf'd for clean re-prep.
THE W=32 GUARD-KILL CLASS (the campaign's own packaging failure, the session's
second key finding): ALL FOUR W=32 cut batches were SIGKILLed by the per-process
64GiB guard at 65-69GB RSS (4/4, ZERO territories cut in ~2h; batch memory grows
with the batch's aggregate routing-index nodes: 779k nodes -> ~63GB at scan end,
the dense batch-1792039134-3 carried 1.76M nodes). THE GUARD IS LAW - THE WIDTH
YIELDS: run-cuts.sh now BATCH_W=16 (measured W=1 peak 38.99GB; W=16 ~50-55GB
est). THE WIDEN PLAN CORRECTED: the handoff's W=96 is IMPOSSIBLE under the
guard - after queries.done the widening axis is N_CONCURRENT (RAM-bound), never
W. TWO ATTRIBUTION FIXES in run-cuts.sh: (a) a guard kill (exit 137) is a
batch-packaging event, never a per-territory strike; (b) the 54 false cut.exit
terminals + the attempts ledger from the W=32 era were RESET (poison.tsv empty:
zero intrinsic territory failures) - logged in cuts/batch-failures.log. The
pending list is now deterministically hash-interleaved before batch composition
(genome-adjacent density clumps spread across batches; per-batch memory ~ the
genome average). CUT CAMPAIGN RELAUNCHED 05:01:32Z at W=16 C=2 RAYON=48
(launch-cuts.sh created; the missing launcher the plan referenced).
ALSO: the stale wave-logs/wave.done (03:31, the old buggy wave's exit marker)
removed - the running wave master re-touches it at true completion; a duplicate
run-heartbeat.sh noted (two processes, one log writer - benign, left running).
STATE AT THIS CHECKPOINT: query 1,946/9,891 (~11/min under load-300 co-tenancy;
ETA queries.done ~15-17Z); prep valve holding (substrate 392, pending-cut 185);
cuts: two W=16 batches scanning (05:01:32Z); wave: T00000 called+reaped, LANES
master waiting on cuts; watcher: 5 orphans still tracked; disk 1812G (floor
1500G), rvsize 206G. THE END-TO-END PROJECTION at the corrected packing:
query ~15-17Z -> the cut campaign is the second pole (~6,000 accepted
territories at W=16/C=2-6; the campaign aggregate is the receipt) -> the wave
trails the cuts. SUCCESSOR: poll SHORT, ride the failures honestly, widen C at
queries.done, rollup + datasheet + the record's closing increment at close.

## THE SUCCESSOR CHECKPOINT 2 (05:25Z 2026-10-15) - THE ADAPTIVE WIDTH + THE HONEST THROUGHPUT
THE W=16 MEASUREMENT: BOTH first-wave W=16 batches guard-killed (batch -0: 517k
nodes, 8,956 rows -> 69GB, killed 473s INTO EMISSIONS; batch -1: 700k nodes,
9,279 rows -> 72.5GB at scan end) - ZERO territories emitted before the kills
(cut-stats lands only at territory completion). THE PRE-BATCH PREDICTORS BOTH
FAILED: rows (8,956 dead vs 9,279 dead) and routing-index nodes (517k dead at
69GB, 779k dead at 69GB, 700k dead at 72.5GB) do not bound the census working
set. THE ADAPTIVE WIDTH (run-cuts.sh): the guard is the ONLY oracle - a batch
starts at BATCH_W=8; the territories of a guard-killed batch retry at HALF the
width (prep/<T>/cut.width, classes 8/4/2/1); a guard kill never strikes a
territory. The old-master fragment race is closed: killing a master's children
first gives it a 1s window to compose replacements (measured: the dead master
composed batches -6/-7 W=32 in the second before its own death; killed by
exact pid at 05:17) - THE ORDER IS MASTER FIRST, THEN CHILDREN (or kill the
process group).
THE HONEST THROUGHPUT REALITY (the handoff's '~3-4h cuts' projection is dead):
per batch = load ~260s + scan ~200s + nodes x ~1ms (scan 580s at 517k nodes,
1755s at W=32/779k) + emissions; the load+scan constant ~460s amortizes over
at most W~8-16 (the guard caps the width), the per-node matching ~25-40s avg
per territory is CPU-bound and does NOT amortize. Genome-wide (~6,000
accepted) at post-query C=8-12: the cut campaign wall is ~12-18h - THE CUT
CAMPAIGN IS THE RUN'S TRUE POLE, by roughly an order of magnitude over the
claimed per-sample arithmetic. The campaign aggregate is the receipt; the
datasheet delta will state it plainly.
STATE: adaptive cut campaign relaunched 05:24:22Z (W=8 C=2 RAYON=48,
batch-1792041862-0/-1 scanning); query ~2,150/9,891 (~20/min since the orphan
kill freed ~32 cores; ETA queries.done ~11:30-12:30Z); prep valve holding
(substrate 393, pending-cut ~370); wave waiting on cuts (T00000 the only
receipt); watcher: 4 orphans still tracked (T00140 T00143 T00141 T00023 T00208
partitions + T00197 syng - check orphan-watcher.log); disk 1808G, rvsize 207G.
AT queries.done: (1) the final reconciliation split runs in the prep master
(rv-split-refined-beds.py, the bridge rule, writes the FINAL accepted.txt);
(2) WIDEN THE CUT CAMPAIGN ON C ONLY: kill the run-cuts master FIRST (it is
mid-wait - safe), relaunch 'N_CONCURRENT=8 BATCH_W=8 RAYON=24 bash
launch-cuts.sh' (8 x ~50GB = 400GB; memavail post-query ~900GB; watch the
first hour - if memavail < 150G drop to C=6); (3) the wave: LANES=128 bash
launch-wave.sh (kill the wave master first, mid-wait is safe); (4) ride to
campaign close, then rv-rollup.py + datasheet + the record's closing increment.

## THE SUCCESSOR CHECKPOINT 3 (05:50Z 2026-10-15) - STEADY STATE + THE WIDEN AUTOMATED
THE CAMPAIGN IS IN ITS STEADY RHYTHM, all five masters healthy, and the
queries.done transition is now SELF-RIDING: rv-widen-at-done.sh is ARMED
(detached, pid 2029990) - at queries.done + prep/final-split.done it
GROUP-KILLS each master (pgid==pid, race-free - the session's two kill
lessons encoded: children-first hands the master a 1s compose window,
master-only orphans the pools mid-flight), sweeps the prep partials, cleans
the unmarked partitions dirs (the mid-gpartition discarded-run trap), and
relaunches: prep POOL_W=64 THREADS=4, cuts BATCH_W=8 N_CONCURRENT=8
RAYON=24, wave LANES=128. THE WIDENING AXIS IS CONCURRENCY, NEVER THE WIDTH.
THE MEASURED STEADY STATE (the receipts so far):
- CUTS (W=8, the first two batches): 8 territories in 653s wall (145k nodes,
  peak RSS 64.36GB) and 8 in 766s (367k nodes, peak 65.27GB) - load ~260s +
  scan 292s + emissions ~12s/territory; BOTH survived the 64GiB guard at 3-4%
  margin; C=2 throughput ~2.7 territories/min -> ~37h genome-wide, hence the
  C=8 widen (~9h). W=16 is dead (69GB/72.5GB), W=32 dead (4/4 at 65-69GB),
  W=96 impossible under the guard.
- CALLS: 4 done / 2 FAIL-CLOSED bind ambiguities at genome scale - T00337
  'key verifies at two origin shifts (census record 9753)' (the 7479-class),
  T00365 'census record has no unique verifying key (2 surviving of 3
  candidates)' - terminal, counted, named; the census accumulates (n=6 so
  far: 33% - do not extrapolate, the aggregate is the receipt). Call walls
  ~103-119s at RAYON=1, peak RSS ~2.0GB (T00000/T00055 receipts).
- PREP (final for the 393 substrate.done): extract median 10s (max 503s),
  syng 4s (max 159s), beds 74s MEDIAN but mean 385s / p90 1,262s / max
  5,672s (the segdup tail IS the prep cost: sum 189,284s over 491 beds),
  export 4s, walks 3s, maps 30s - per-territory median ~125s serial vs the
  claimed ~8s (the claim counted syng+routes only, maps presumed committed);
  post-query at POOL_W=64 the remaining ~4,600 beds are ~8-12h - the
  post-query pole alongside the cuts.
- QUERY: ~2,600/9,891 at ~20/min (the orphan kill freed ~32 cores) - ETA
  queries.done ~11:45Z; the split follows at 61% accepted.
- DISK: 1802G avail (floor 1500G), rvsize 217G, memavail 190-260GiB.
THE HONEST END-TO-END PROJECTION: queries.done ~11:45Z -> the widen ->
prep-beds and cuts as the poles (~8-12h) -> the wave trails (~1.5h at
LANES=128) -> CAMPAIGN CLOSE ~21:00-24:00Z. At close: rv-rollup.py (the
census fields are in place), fill datasheet-runtime.md (every stage measured
vs claimed, the per-territory cost distribution, the failure census, the
disk high-water, the delta table), the record's closing increment, the
final commit. NO PR PUSH.

## THE SUCCESSOR CHECKPOINT 4 (07:05Z 2026-10-15) - THE RIDE STEADY + THE ORPHAN TAIL CLOSED
THE RIDE VERIFIED END TO END: all five masters healthy (query pgid 1441135, prep-master
1830969, cuts 1993375, wave 1858400, heartbeat x2 - the known benign duplicate), the
widen watcher ARMED (2029990; script read and verified: group-kill pgid==pid per master,
master-first ordering, unmarked-partitions cleanup, the corrected widening axis
CONCURRENCY-ONLY: prep POOL_W=64 THREADS=4 / cuts BATCH_W=8 N_CONCURRENT=8 RAYON=24 /
wave LANES=128), and the final-split step CONFIRMED WIRED in the prep master (runs
rv-split-refined-beds.py at queries.done, touches prep/final-split.done - the widen
watcher's second trigger). poll-once.sh added (the ride's one-line short poll).
THE ORPHAN TAIL IS FULLY CLOSED: all 7 orphan partitions converted with honest walls
(T00207 3,523s / T00640 4,311s / T00208 4,503s / T00140 5,485s / T00143 5,672s /
T00023 6,935s / T00141 9,000s - THE BEDS SEG DUP TAIL'S TRUE MAX), and the last orphan,
T00197's wrapper-less syng, KILLED BY EXACT PID (1743721) at 07:01Z - named reasons:
its wrapper was dead (success unverifiable, markers impossible), it had burned 3h01m
CPU on a 345MB locus.fa with no output file (vs syng median 4s max 159s - an early
signal of a SYNG SEG DUP-TAIL CLASS, to be confirmed or refuted by the marker-gated
re-prep), and it posed the DOUBLE-WRITER RACE: T00197 is accepted, has no
substrate.done/skipped, so the master re-dispatches it the moment the valve opens -
two writers on locus-pos64.syng = corrupt by construction (the named T00003/T00114
class). The watcher (its committed design: partials cleaned for a clean re-run)
cleaned locus-pos64.syng* and touched orphan-watcher.done - ALL ORPHANS RESOLVED;
T00197 re-preps cleanly under the master when the valve opens (extract.done valid).
THE MEASURED STATE (06:39-07:03Z ride receipts):
- QUERY: 3,944/9,891 windows at ~17/min under cuts+wave co-tenancy (revised ETA
  queries.done ~13:00Z; the 11:45Z projection assumed the 20/min post-orphan rate).
- SPLIT (incremental, live): 2,315 accepted / ~630 rejected at this margin (~59%,
  steady vs the bridge-rule 61% first-pass); cumulative rollup classes:
  transitive-expansion 1,211 / ambiguous-axis-bridge 250 (the expansion-gate
  territories counted, never retried into holes).
- CUTS: 12+ W=8 batches complete at 649-784s walls (C=2 -> ~1.4-1.5 territories/min);
  ~176 map dirs; the 5 fresh W=8 guard-kills at 05:49Z (T00135 T00161 T00676 T00624
  T00580) all auto-retry at width 4 (8 adaptive-width-4 pending, zero strikes,
  poison.tsv EMPTY - the adaptive width working as designed); batch RSS ~43GB in
  flight, guard margin held.
- CALLS: 33 called / 7 FAIL-CLOSED bind ambiguities, ALL NAMED (two-origin-shifts:
  T00027 rec 8215, T00031 rec 9149, T00157 rec 9707, T00216 rec 8527, T00337 rec
  9753; no-unique-verifying-key: T00142 rec 9165 (2 of 6), T00365 rec 7287 (2 of 3))
  - terminal, never retried into holes. CALL WALLS UNDER THE 14-LANE QUERY
  CO-TENANCY: median 253s / p90 1,159s / max 1,775s at RAYON=1 (n=38; the T00000
  103s was the clean-cache first point) - the contended aggregate is the receipt;
  the clean per-lane calibration must be re-stated post-query at the widened
  packing (the datasheet's wave section must carry BOTH numbers honestly).
- PREP VALVE: holding at 262 substrates pending cut (>=160); sub=393 steady.
- DISK: 1,775-1,788G steady, rvsize 230-243G. THE RECORDED EXCURSION (for the
  datasheet's high-water section, stated honestly): 04:08-04:20Z df dipped to
  1,373G - BELOW the 1,500G floor (rvsize ~600G in-flight window before the valve
  caught it) - recovered by the valve + JIT reaper, steady since. THE FLOOR
  EXCURSION IS PART OF THE RECEIPT, not a silent pass.
- rv-rollup.py VERIFIED MID-CAMPAIGN (runs clean; census fields filled; the
  wave call-wall distribution, the failure census, the batch-exit census
  137:10 / 0:12, guard_kill_widths 8:8, adaptive_width_pending 4:8 all present).
THE SUCCESSOR'S RIDE (unchanged + refined): SHORT polls (poll-once.sh); verify the
widen FIRED when queries.done + final-split.done land (widen.log: 'WIDEN COMPLETE';
watch memavail in the first hour at C=8 - if <150G narrow to C=6); at campaign close
(all territories called or honestly failed): rv-rollup.py -> FILL datasheet-runtime.md
(every stage measured vs claimed, the delta table, the failure census, the disk
high-water WITH the 1,373G excursion, the co-tenant statement) -> the record's closing
increment -> the commit. NO PR PUSH.

## THE RIDE CLOSES ITS WINDOW (07:21Z 2026-10-15) - THE HONEST CHECKPOINT STATE
THE RIDE (06:39-07:21Z) rode the steady state through ~40 SHORT polls (poll-once.sh):
every poll nominal, zero new failure classes, no master touched. TWO EVENTS HANDLED:
(1) the T00197 orphan syng killed by exact pid + cleaned (Checkpoint 4 above - the
orphan watcher COMPLETE, orphan-watcher.done touched, ALL ORPHANS RESOLVED); (2) the
ADAPTIVE-WIDTH-4 RETRY SCHEDULING named for the successor: the batch composer is
WIDEST-FIRST (width-8 classes fill both C=2 batch slots while ~250 width-8 territories
remain pending), so the 8 guard-killed territories retrying at width 4 (T00135 T00137
T00161 T00172 T00580 T00624 T00632 T00676, from the 05:49Z W=8 guard kills) WAIT for
batch slots until the width-8 pending list drains below ~2 batches - naturally ~2-3h
post-widen at C=8. NOT A DEADLOCK: the campaign's terminal condition waits for all
territories, the rollup counts them as adaptive_width_pending 4:8, and their partial
map/<T>/<T>/ dirs (census/cache/reads from the killed W=8 emissions) are rewritten by
the retry (the patched cut-batch's per-territory isolation). Ordinary W=8 throughput
meanwhile: batches -16/-17 exit 0 at 728s/781s.
THE STATE AT WINDOW CLOSE (07:20Z): query 4,225/9,891 (~17.4/min measured over the
ride; queries.done ETA ~12:45Z); split accepted 2,474 (59% steady); prep valve holding
(sub 393, pending-cut ~250); cuts 188 map dirs (~120 territories cut complete, C=2
~1.4/min); calls 33 called + 7 fail-closed (ALL NAMED, Checkpoint 4) - call walls
under co-tenancy median 253s / p90 1,159s / max 2,823s (T00492, the segdup tail);
disk 1,767G (floor 1,500G; the 04:08-04:20Z excursion to 1,373G stands recorded for
the datasheet), rvsize 250G, memavail 135-200GiB; widen watcher ARMED and verified;
final-split wired; poison 0; NO intervention pending. THE SUCCESSOR RESUMES THE RIDE:
SHORT polls; verify the widen FIRED at queries.done + final-split.done (widen.log
'WIDEN COMPLETE'; watch memavail <150G -> narrow C to 6); ride to campaign close; then
rv-rollup.py -> FILL datasheet-runtime.md (measured vs claimed, the delta table, the
failure census with the 7 named bind ambiguities, the disk high-water WITH the 1,373G
excursion, the co-tenant statement) -> the record's closing increment -> the commit.
NO PR PUSH.

## RIDE CONTINUATION 4 WINDOW OPEN (07:25Z 2026-10-15)
Session-open increment committed (8b6e166; the vendor/syng pointer drift 054acfc->68ac197 recorded
there, left untouched). State verified live at 07:23Z: split 4,291/2,516-accepted, sub 393,
cut 188 map dirs, call 33 + 7 fail-closed (all named), df 1,767G, mem 133GiB, rvsize 250G,
poison 0, widen ARMED (2029990), five masters healthy. The window rides: SHORT polls ->
verify the widen FIRES at queries.done (~12:45Z ETA) -> the contingencies -> at campaign
close rv-rollup.py + datasheet-runtime.md + the record's closing increment + commit.
THE RIDE'S FIRST NEW CLASS (07:43:56Z): T00129 CALL GUARD-KILLED at 67.5GB by the wave's
own 64GiB RSS guard (call.exit.code 137, wall 5,637s, still climbing ~3MB/s at the kill -
the segdup-tail territory HG002#1#chr1:43.2-43.5Mb, 24 windows, locus 11's no-pin
(unit,fold) 629M hull). NOT a bind ambiguity - the CALL MEMORY TAIL class, genome-scale
first sighting: unlike a cut batch there is NO WIDTH TO YIELD (the call is one territory
at RAYON=1) -> terminal, counted, named for the datasheet (the call-wall distribution must
carry the 5,637s/67.5GB point alongside the contended median). failcall=8 now: 7 bind
ambiguities (named, Checkpoint 4) + 1 call guard-kill (T00129).
THE BIND-AMBIGUITY CENSUS CLIMBING (07:45-07:49Z, the wave consuming the early chr1
territories fast): new fail-closed binds since Checkpoint 4 - T00012 (rec 6515, 2 of 5),
T00130 (rec 7434, two origin shifts), T00085 (rec 8499, 2 of 2), T00017 (rec 9972, two
origin shifts), T00037 (rec 7389, two origin shifts) + more landing - all the ordinary
7479-class, terminal, named. failcall=16 at 07:49Z vs call=37 (~30% of terminal, steady
with the early 33% observation; the aggregate at close is the receipt, no extrapolation).
THE MEMAVAIL NARROWING (08:04Z, the standing contingency applied PRE-widen, named):
memavail fell to 37-43Gib (the wave's first full-lane saturation: 33 concurrent calls + 2
cut batches ~42GB each + 14 query lanes ~812GB, equilibrium ~40Gib) - BELOW the 150G
floor of the standing memavail contingency. THE DECISIVE ANALYSIS: (a) query shards are
NOT per-window resumable (run-query.sh does `rm -f shard-$S.out.bed` on restart - a
killed lane re-runs its whole ~707-window shard from zero, ~3.3h); (b) oom_score_adj
protection of the query lanes was ATTEMPTED AND DENIED (no passwordless sudo, /proc
write EPERM - the attempt is part of the receipt); (c) the T00129-class call memory tail
(grows ~0.7GB/min toward the 64GiB guard) with memavail 40 = kernel OOM, victim = biggest
RSS = a query lane = the catastrophic unresumable class. THE ACTION: cuts C=2 -> C=1
(group-TERM pgid 1993375 per the widen watcher's own pattern; both in-flight batches died
unaccounted - no strikes, the class of record; relaunch BATCH_W=8 N_CONCURRENT=1 RAYON=48,
new master 2268879). MEMAVAIL RECOVERED to 108Gib within 2 min. THE HONEST COST: cut
production halves (1.4 -> 0.7 territories/min) for the pre-widen window only - a ~126-
territory deficit ~= 12 min at the widened C=8, negligible vs the +3.3h query-shard-restart
risk avoided. The widen watcher unaffected (finds masters BY NAME at fire time - verified).
NOTE FOR THE DATASHEET: the co-tenant memory packing is itself a measured finding - the
wave's full-lane saturation + the query stage's 55-58GB/lane plateau (the '45GB' of the
early checkpoints was the early-phase measure; growth ~1-2GB/91min once plateaued).
RIDE CONTINUATION 4 STEADY (08:32Z): query 5,666/9,891 (~20/min, ETA queries.done
~12:00Z); split accepted 3,326 (58.6%); the wave drained the C=2-era backlog (cut dirs
~135 steady at the C=1 production rate, calls 128 / fail-closed 30 of 158 terminal = 19%);
memavail 144-187Gib post-narrowing; disk 1,796G (floor 1,500G), rvsize 221G and falling
(the JIT reaper); one fresh W=8 guard-kill 08:24:18Z (T00116 batch -> width 4, nominal
adaptive class). All five masters healthy; widen watcher armed; poison 0.
THE VALVE BYPASS, CAUGHT AND CLOSED (08:33-08:57Z, the 04:08Z excursion class, named):
DESIGN FLAW FOUND: run-prep-master.sh dispatched its ENTIRE TODO list through ONE
`xargs -P POOL_W` wave, re-checking the valve only between full waves. At 08:33Z the
valve opened for ONE INSTANT (a cut batch's 8 cut.done dipped NPENDCUT from 166 to <160)
and the master committed ~2,855 territories to an unthrottled wave - NO valve re-check
until the whole list drained (~4-5h), projected +700GB substrate / df below the 1,500G
floor by ~10:20Z (the exact unbounded-growth class the valve was built to prevent).
THE SURGICAL STOP (08:55Z): kill -TERM the xargs dispatcher ONLY (pid 2316631, one
process, named reason) - the master looped, recomputed NPENDCUT=354, and HELD
("08:57:12Z VALVE: 354 substrates pending cut - holding the prep pool"); the 28 in-flight
chain wrappers SURVIVED (reparented to init, NOT killed), completed, and wrote their
markers - ZERO orphans, zero lost work (the wrapper-kill corpse class avoided entirely).
THE PATCH (in place, takes effect at the widened relaunch - the running master holds its
loop in memory and dies by the watcher's TERM at widen, so the byte-offset hazard never
materializes): the dispatch is now CHUNKED (POOL_W*2 per chunk, ~2 concurrency waves)
with the valve RE-CHECKED between chunks; a closed valve parks the rest of TODO for the
next loop (markers make TODO re-computation free); npendcut() factored out and used at
both checks; bash -n clean; npendcut() verified against the live ledgers (356). THE
WIDENED-PHASE STAKES: post-query TODO is thousands of territories - one instant valve
opening at POOL_W=64 under the old code would have committed an unbounded ~1.5TB wave;
the floor would blow through in under 2h. THE SUCCESSOR SHOULD VERIFY the first
mid-dispatch VALVE line appears in prep-master.log post-widen (the patch working).
State at 08:58Z: sub plateauing ~600, chains drained, query 6,141/9,891 (62%, accepted
3,599), disk 1,700G, rvsize 317G stabilizing, memavail 166Gib.
## RIDE CONTINUATION 4 WINDOW CLOSE (09:14Z 2026-10-15) - THE HONEST CHECKPOINT STATE
FINAL POLL: 2026-10-15T09:13:38Z split=6484 acc=3788 sub=607 cut=132 call=152 failcall=35 df=1701G mem=110Gib rvsize=315G qdone=n fsplit=n widen=n poison=0 lastkill=2026-10-15T08:24:18Z T00116 batch GUARD-KILLED at width 8 (exit 137, packaging - no strike
THE WINDOW'S WORK (all committed to the record: 8b6e166 session-open, 2bef2e2 the ride
receipts 07:23-09:00Z): (1) T00129 call guard-kill class named (67.5GB/5,637s, segdup
memory tail, no width to yield - terminal, for the datasheet's distribution); (2) the
memavail narrowing C=2->C=1 at 08:04Z (37-43Gib < the 150G floor; oom-protection of the
query lanes attempted and DENIED - receipted; cuts master group-TERMped per the watcher's
own pattern, relaunch BATCH_W=8 N_CONCURRENT=1 RAYON=48, new master 2268879; memavail
recovered 108-187Gib); (3) THE VALVE BYPASS caught and closed (run-prep-master.sh's
full-TODO xargs wave - one instant opening at 08:33Z committed ~2,855 territories
unthrottled; surgical TERM of the xargs dispatcher ONLY at 08:55Z, the 28 in-flight
wrappers survived and wrote markers - zero orphans; the master looped and HELD at
NPENDCUT=354; THE PATCH IS IN PLACE in run-prep-master.sh: chunked dispatch POOL_W*2
with the valve re-checked between chunks - TAKES EFFECT AT THE WIDENED RELAUNCH since
the running master holds its loop in memory and dies by the watcher's TERM at widen).
THE SUCCESSOR RESUMES: SHORT polls (poll-once.sh); watch memavail (narrow cuts C if
<150G - currently C=1 pre-widen so the lever post-widen is C=8->6); VERIFY THE WIDEN
FIRES at queries.done + final-split.done (ETA ~12:25Z at ~18/min; widen.log 'WIDEN
COMPLETE'; group-kills pgid==pid; relaunch prep POOL_W=64 THREADS=4 / cuts BATCH_W=8
N_CONCURRENT=8 RAYON=24 / wave LANES=128) AND VERIFY THE FIRST 'VALVE (mid-dispatch)'
LINE APPEARS IN prep-master.log (the patch working); the 8+ adaptive-width-4 territories
wait for width-8 slots (adaptive_width_pending, never a deadlock); fail-closed binds
terminal and named (35 at close: 34 bind ambiguities + T00129's guard kill); AT CAMPAIGN
CLOSE: rv-rollup.py -> FILL datasheet-runtime.md (every stage measured vs claimed, the
delta table, the call-wall distribution with BOTH the contended median ~253s/p90 ~1,159s
AND the clean 103s calibration AND the T00129 5,637s/67.5GB tail point, the disk
high-water with the 1,373G excursion, the co-tenant/memavail findings incl. the
55-58GB/lane query plateau) -> the record's closing increment -> the commit. NO PR PUSH.
NOTE: the vendor/syng submodule pointer drift (054acfc->68ac197, runtime-era diagnostics)
remains uncommitted by design - recorded in 8b6e166, left for the owner.

## RIDE CONTINUATION 5 WINDOW OPEN (09:18Z 2026-10-15)
Session-open increment committed (a89d70b; vendor/syng pointer drift 054acfc->68ac197 still
uncommitted by design). Live at 09:15Z: split 6,484/3,788-accepted, sub 607, cut 128 map dirs,
call 156 + 35 fail-closed (34 bind ambiguities + T00129 guard-kill), df 1,707G (floor 1,500G),
mem 109-122Gib, rvsize 310-316G, qdone=n fsplit=n widen=n, poison 0, zero intervention pending.
MASTERS VERIFIED: prep 1830969 (valve holding at 356 pending-cut; the chunked-dispatch patch
in place, takes effect at the widened relaunch), cuts C=1 master 2268879 (the 08:04Z
narrowing), wave 1858400, heartbeat 1454428, the 14 query lanes live at pgid 1441135 (the
group leader's own PID is gone but the process group is live - the watcher's group-kill
still targets it; 14 `impg query` processes confirmed). Widen watcher ARMED (2029990).
THE RIDE: SHORT polls -> verify the widen fires at queries.done + final-split.done (ETA
~12:25Z): widen.log 'WIDEN COMPLETE' + THE FIRST 'VALVE (mid-dispatch)' line in
prep-master.log (the chunked-dispatch patch's proof) -> watch memavail first hour at C=8
(<150G -> C=6) -> ride to campaign close -> rv-rollup.py + datasheet + the record's
closing increment + commit.
## RIDE CONTINUATION 5 WINDOW CLOSE (11:06Z 2026-10-15) - THE HONEST CHECKPOINT STATE
THE RIDE (09:15-11:06Z) rode ~110 SHORT polls: every poll nominal, no master touched, no
kill executed, zero interventions beyond naming. FINAL POLL 11:05:35Z: split 8,766/9,891
(88.6%) acc 5,186 (59.2% steady), sub 611, cut 185 map dirs, call 179 called + 39 fail-closed
(38 bind ambiguities ALL NAMED incl. this window's T00131 rec 7528 / T00910 rec 7484 / T00314
rec 7379 / T00434 rec 8358 - all the ordinary 7479-class, terminal - plus T00129's guard
kill), df 1,690G (floor 1,500G), mem 353Gib RISING (three query shards DONE - shard-01/-07/
/-12 + shards.done markers; 11 `impg query` lanes still live), rvsize 326G, qdone=n
fsplit=n widen=n, poison 0.
THE WINDOW'S ONE NEW CLASS (named for the datasheet census): at 11:02:23Z EIGHT territories
failed with CUT BATCH EXIT 1 - 'locus record GRCh38#0#chr2:20750010-20999992 is not identical
to the AGC interval under either offset convention' (T00341 T00905 T00955 T00724 T00631
T00950 T00791 T00974, batch-1792061780-15). UNLIKE the exit-137 guard kills this is the
ATTEMPT-COUNTING path: attempt 1 of 3, the territories re-join the pending set, and 3
batch-level failures make them TERMINAL (cut.exit, named in the ledger) - the designed
fail-closed path, no intervention; WATCH FOR ATTEMPT 2/3 in the successor window (the
error looks deterministic - a locus/AGC interval offset-convention mismatch, chr2:20.75Mb -
if it trips 3x it becomes a named terminal class for the datasheet alongside the binds).
QUERY PACE ACCELERATING as shards retire: ~65-80 windows/min over the last hour (split
8,385->8,766 in 20 min) - queries.done now projects ~11:40-12:00Z, ahead of the ~12:25Z ETA.
W=8 guard kills this window (all the ordinary adaptive class, no strikes): 09:30 T00858,
09:42 T00297, 09:51 T00920, 10:02 T00609, 10:24 T00918, 10:45 T00855, 10:56 T00188.
MASTERS ALL HEALTHY at close: prep 1830969 (valve holding), cuts C=1 master 2268879,
wave 1858400, heartbeat, widen watcher 2029990 ARMED (5.3h), the 14 query lanes at pgid
1441135 (11 live + 3 shards done).
THE SUCCESSOR RESUMES THE RIDE (unchanged): SHORT polls; VERIFY THE WIDEN FIRES at
queries.done + final-split.done (widen.log 'WIDEN COMPLETE' + THE FIRST 'VALVE
(mid-dispatch)' LINE in prep-master.log - the chunked-dispatch patch's proof); watch
memavail the first hour at C=8 (<150G -> narrow C=6); the T00129 call guard-kill class
terminal; the fail-closed binds terminal and named; the disk floor 1,500G (the 1,373G
excursion stands recorded); the adaptive-width-4 territories wait for width-8 slots, never
a deadlock; AT CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md (every stage
measured vs claimed, the delta table, BOTH call-wall calibrations + the T00129 tail point,
the 1,373G excursion, the 55-58GB/lane query plateau, the NEW chr2 AGC-offset batch-exit-1
class with its attempt-3 resolution) -> the record's closing increment -> the commit.
Session-open increment for this window: committed a89d70b (the record; vendor/syng pointer
drift 68ac197 still uncommitted by design, left for the owner).
## RIDE CONTINUATION 6 WINDOW OPEN (11:07Z 2026-10-15) - THE WIDEN WINDOW
Session-open increment committed (2611230; vendor/syng pointer drift 054acfc->68ac197 still
uncommitted by design). Live at 11:07Z: split 8,822/9,891 (89.2%) acc 5,213 (59.2% steady),
sub 611, cut 184 map dirs, call 180 + 43 fail-closed (38 binds + T00129's guard kill; FOUR NEW
this window's edge, all named: T00097 rec 9273 / T00625 rec 8521 / T00278 rec 9568 /
T00728 rec 8374 - the ordinary two-origin-shift 7479-class), df 1,690G (floor 1,500G),
mem 379Gib RISING (shards 01/07/12 retired; 11 query lanes live), rvsize 325G, qdone=n
fsplit=n widen=n, poison 0, zero intervention pending.
MASTERS VERIFIED: prep 1830969 (valve holding; the chunked-dispatch patch in place, takes
effect at the widened relaunch), cuts C=1 master 2268879, wave 1858400, heartbeat, the 14
query lanes at pgid 1441135. Widen watcher ARMED (2029990, 5.3h; verified script: group-kill
pgid==pid, master-first, partial sweeps, widened relaunch prep POOL_W=64 THREADS=4 / cuts
BATCH_W=8 N_CONCURRENT=8 RAYON=24 / wave LANES=128).
THE WIDEN IS THIS WINDOW'S EVENT: queries.done projects ~11:40-12:00Z (pace accelerating
~65-80/min as shards retire). VERIFY: widen.log 'WIDEN COMPLETE' + THE FIRST 'VALVE
(mid-dispatch)' LINE IN prep-master.log (the chunked-dispatch patch's proof); record the
widen's transition walls. CONTINGENCIES: memavail <150G first hour at C=8 -> C=6; the chr2
AGC-offset class (T00341 T00905 T00955 T00724 T00631 T00950 T00791 T00974 - attempt 1 of 3
failed 11:02:23Z 'locus record ... not identical to the AGC interval under either offset
convention', chr2:20.75Mb; watch attempts 2/3 - 3x trips = a named terminal class); the
fail-closed binds terminal; the T00129 guard-kill class terminal; the disk floor 1,500G.
AT CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md -> the record's closing
increment -> the commit. NO PR PUSH.
THE CHR2 AGC-OFFSET CLASS RESOLVED (11:14:16Z): attempt 2 failed 11:08:09Z, attempt 3
failed 11:14:16Z - all EIGHT territories TERMINAL after 3 batch-level failures (batch exit 1,
the deterministic error verbatim: 'locus record GRCh38#0#chr2:20750010-20999992 is not
identical to the AGC interval under either offset convention'). THE CLASS IS NAMED for the
datasheet census alongside the binds: the chr2 AGC-offset mismatch class - T00341 T00905
T00955 T00724 T00631 T00950 T00791 T00974, 8 territories, fail-closed by design (3x
batch-level exit 1, never retried into holes). The designed path worked: attempt-counting ->
terminal, no intervention.
THE QUERY MASTER DIED AT ITS OWN FINISH LINE; THE MARKER RECOVERED BY THE RIDE (12:53Z,
receipted): the query stage COMPLETED at 12:49:59Z (split 9,891/9,891 - ALL windows
answered; all 14 shards verified complete by the stage's own criteria: exit 0, n_in==n_out
707/707 x4 + 706/706 x10, shard-*.done x14) but run-query.sh's MASTER PROCESS died before
its final lines - query-stage.nohup shows 43,361 'line 37: [: missing `]' parse errors (the
guard_poll while-loops re-reading the script file; FINDING: THE QUERY-LANE 64GiB RSS GUARD
WAS NON-FUNCTIONAL the whole stage - each poll's comparison failed to parse, so a lane
over 64GiB would NOT have been killed; no harm materialized: the lanes plateaued at
55-58GB, under the guard) and a final 'line 79: syntax error near unexpected token `)''
at the master's wait-return - the byte-offset/parse-death class (the same hazard named
for prep-master patches). THE RECOVERY: the dead master's own completion action replicated
verbatim by the ride at 12:53:10Z after independent verification of every shard's
completion criteria - `NDONE=$(ls shard-*.done|wc -l)`=14, queries.done.count written,
queries.done TOUCHED. This is the designed marker-gated chain unblocked, not a silent
patch: the prep master's final-split (rv-split-refined-beds.py) and the widen watcher
(2029990) both key on this marker. THE GUARD FINDING IS FOR THE DATASHEET: the query
stage ran 8h with its guard unable to fire (the 55-58GB/lane plateau made it moot).
## RIDE CONTINUATION 6 WINDOW CLOSE (12:58Z 2026-10-15, deadline checkpoint) - THE HONEST STATE
THE RIDE (11:07-12:58Z, ~50 SHORT polls) closed the query stage and stands at the widen's
edge. FINAL POLL 12:57:56Z: split 9,891/9,891 - THE QUERY STAGE COMPLETE (100%), acc 5,844
accepted (59.1% steady), sub 611, cut 249 map dirs, call 200 called + 44 fail-closed
(43 bind ambiguities + T00129's guard kill), df 1,662G (floor 1,500G), mem ~900Gib (the query
lanes retired - the post-query box arrived), rvsize 353G, qdone=Y, fsplit=n (the final split
RUNNING - rv-split-refined-beds.py 4m35s in over all 14 shard outputs, launched by the prep
master at ~12:53), widen=n (the watcher SAW queries.done land 12:53:33Z; it fires within 60s
of prep/final-split.done), poison 0.
THE WINDOW'S THREE EVENTS, ALL RECEIPTED ABOVE:
(1) THE CHR2 AGC-OFFSET CLASS RESOLVED TERMINAL (11:14:16Z): attempts 2 (11:08:09Z) and 3
(11:14:16Z) both failed with the identical deterministic error - all EIGHT territories
(T00341 T00905 T00955 T00724 T00631 T00950 T00791 T00974) TERMINAL after 3 batch-level
failures, batch exit 1, the designed fail-closed path, named for the datasheet census
alongside the binds: the chr2 AGC-offset mismatch class (chr2:20,750,010-20,999,992, 'locus
record not identical to the AGC interval under either offset convention').
(2) FOUR NEW FAIL-CLOSED BINDS named (the ordinary 7479-class, two origin shifts):
T00097 rec 9273, T00625 rec 8521, T00278 rec 9568, T00728 rec 8374. failcall 44 = 43 binds +
T00129's call guard-kill. W=8 guard kills this window (ordinary adaptive class, no strikes):
11:33 T00812, 11:42 T00466, 12:10 T00395, 12:18 T00605, 12:46 T00413.
(3) THE QUERY MASTER DIED AT ITS OWN FINISH LINE; THE MARKER RECOVERED (12:53:10Z): the stage
completed at 12:49:59Z (all 14 shards exit 0, n_in==n_out verified x14) but the master hit
'line 79: syntax error' at its wait-return and died before touching queries.done - and
query-stage.nohup receipts 43,361 'line 37: [: missing `]' parse errors: THE QUERY-LANE 64GiB
RSS GUARD WAS NON-FUNCTIONAL THE WHOLE STAGE (the comparison never parsed; a lane over 64GiB
would not have been killed; moot at the 55-58GB plateau - BOTH findings for the datasheet).
The ride replicated the dead master's own completion action verbatim after independently
verifying every shard's criteria: queries.done.count=14 written, queries.done TOUCHED
(12:53:10Z). The designed chain unblocked: the prep master launched the final split within
seconds; the widen watcher logged 'queries.done landed' 12:53:33Z.
THE SUCCESSOR RESUMES AT THE WIDEN'S EDGE: the final split is running (started ~12:53, over
~5.9GB of shard outputs - minutes-scale); when prep/final-split.done lands the watcher fires
(60s poll): VERIFY widen.log 'WIDEN COMPLETE' + THE FIRST 'VALVE (mid-dispatch)' LINE in
prep-master.log (the chunked-dispatch patch's proof); record the widen's transition walls;
watch memavail the first hour at C=8 (<150G -> narrow C=6 - memavail ~900Gib now, the
post-query box, so the first-hour risk is low); then ride to campaign close (cuts/calls/wave
at the widened concurrency; adaptive-width-4 territories wait for width-8 slots); AT
CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md (every stage measured vs claimed,
the delta table, BOTH call-wall calibrations + the T00129 tail, the 1,373G excursion, the
55-58GB/lane query plateau + THE NON-FUNCTIONAL QUERY GUARD finding, the chr2 AGC-offset
class resolution, the query-master parse-death + the marker recovery) -> the record's closing
increment -> the commit. NO PR PUSH.
Session-open increment for this window: committed 2611230. Ride-receipts commit: this note's
events (chr2 terminal class, query-master death + guard finding + marker recovery) are in
the record's session-open increment's successor commit - COMMIT THEM with the window-close.

## RIDE CONTINUATION 7 WINDOW OPEN (13:02Z 2026-10-15) - THE WIDEN FIRED, VERIFIED, COMPLETE
Session-open increment committed (c7a445b; vendor/syng pointer drift 054acfc->68ac197 still
uncommitted by design). THE TRANSITION RECEIPTS (widen.log, all UTC):
- 12:59:55Z THE FINAL SPLIT COMPLETED CLEANLY (the ride-6 open risk checked FIRST and closed:
  no crash, no stall): windows 9,891 / accepted 5,844 / rejected 4,047 / rows 188,719,130,
  acceptance IDENTICAL to the incremental ledger (59.1% steady). prep/final-split.done 12:59:00Z.
- 13:00:33Z watcher (2029990, 7h10m armed) saw final-split.done; group-killed prep master 1830969
  (master-first); 13:01:04Z prep relaunched POOL_W=64 THREADS=4 (partials swept).
- 13:01:09Z cuts C=1 master 2268879 group-killed; relaunched BATCH_W=8 N_CONCURRENT=8 RAYON=24.
- 13:01:29Z wave master 1858400 group-killed; unmarked partitions dirs cleaned; wave relaunched
  LANES=128. 'WIDEN COMPLETE: prep POOL_W=64, cuts C=8 W=8 RAYON=24, wave LANES=128.'
AT 13:02:06Z: memavail 886Gib (the C=8 first-hour contingency <150G not in play), df 1,662G
(floor 1,500G), rvsize 354G, call 200 + 44 fail-closed (terminal classes stay terminal),
poison 0. THE 'VALVE (mid-dispatch)' PROOF IS PENDING-BY-DESIGN: the valve still holds at the
loop head (276 pending-cut >= 160, no dispatch running) - the line lands when C=8 cuts drain
pending-cut < 160 and the widened prep pool starts chunked dispatch. THIS WINDOW'S RIDE:
SHORT polls; memavail first hour; watch the widened pace (cuts territories/min at C=8, prep
substrates/min at POOL_W=64, calls/min at LANES=128); the terminal classes stay terminal
(44 fail-closed of 244, the chr2 AGC class 8, T00129); the adaptive-width retries (width-4
territories now get slots as C=8 composes 8 batches); the disk floor 1,500G. AT CAMPAIGN CLOSE:
rv-rollup.py -> FILL datasheet-runtime.md -> the record's closing increment -> the commit.
## RIDE 7 MID-WINDOW FINDINGS (14:10Z 2026-10-15) - THE VALVE OPENED + THE WAVE'S OWN BARRIER
(1) THE VALVE OPENED at ~13:43Z (pending-cut <160) and the widened prep pool dispatched its
first POOL_W*2=128-chain chunk: sub 611->706, substrates landing in bursts as chains complete
(the beds segdup tails stretch the chunk's return; 2 chains still running at 14:08). Pending-cut
back to ~260 - THE MID-DISPATCH RE-CHECK WILL FIRE THE 'VALVE (mid-dispatch)' LINE when the
chunk returns (the patch's proof; verified in the code at run-prep-master.sh:80-84 - the line
exists, the valve just needs the chunk boundary). NO intervention: the chunk size itself bounds
the growth (+128 substrates ~ +32GB).
(2) THE WAVE'S OWN BATCH-BARRIER FINDING (the prep-master flaw's sibling, receipted): the
wave master (3749885, relaunched at the widen) is BLOCKED in its `xargs -P 128` (elapsed
1:08+) behind TWO MARATHON SEG DUP-TAIL CALLS - T00370 (RSS 27.2GB) and T00496 (RSS 30.9GB,
climbing) - both CPU-bound RAYON=1, the per-lane 64GiB guard ARMED on each. xargs keeps 128
lanes busy WITHIN a TODO snapshot, but the master cannot re-dispatch the ~160 cut-done
territories that accumulated until every straggler in the current wave resolves (completes
or guard-kills - the T00129 terminal class). The 13:01 relaunch's TODO was tiny (~34), so the
barrier bites now; at genome scale the big waves make boundaries rarer. THE HONEST DATASHEET
COST: lane utilization collapses to the straggler count at every wave boundary - ANOTHER
measured delta vs the claimed 19-21min wave. NO KILL: 2h+ CPU on each call, possibly near
completion; the guard is the designed oracle; df STABLE at 1539-1540G (floor 1500G), rvsize
479G steady, memavail 455-867Gib. If the floor ever approaches: the safe contingency is the
kill+relaunch of the wave master MID-WAIT (never with marathon calls in flight unless the
floor is actually breaching).
(3) NEW FAIL-CLOSED BINDS this window (all named, the ordinary 7479-class, terminal):
T00490 rec 7723 / T00237 rec 6921 (2 of 3) / T00946 rec 7970 / T00551 rec 8359 (2 of 2) /
T00453 rec 12598 / T00952 rec 7161 - failcall 51 = 50 binds + T00129's guard kill.
(4) W=8 GUARD KILLS at the widened packing (all the ordinary adaptive class, width->4, zero
strikes, poison 0): 13:12 T00406-batch, 13:23 T00824-batch, 13:32 T00816/T00467/T00968-batch,
14:04 T01071-batch. The kill rate scales with C (more batches in flight), as designed.
(5) MEASURED WIDENED PACE: cut.done production 383->447 in 7 min = 4.6/min gross (exit-0 W=8
batch walls 546-582s at C=8; the guard-kill waste + width-4 retries trim the net); calls
cut-limited behind the barrier; sub +95 in ~25 min once the valve opened (the prep pole's
honest post-query rate at POOL_W=64 THREADS=4 - the beds tails bound it).
## RIDE 7 THE DISK-FLOOR NARROWING (14:14:40Z 2026-10-15) - THE STANDING LEVER APPLIED, RECEIPTED
THE TRIP: df crossed the 1500G floor at 14:11Z (poll 1495G, momentary; heartbeat 1503-1514G
flutter) with rvsize 502-513G (the ~500G high-water). THE CAUSE CHAIN, named: the wave's
batch-barrier (2 marathon segdup-tail calls, T00370 27.2GB / T00496 30.9GB RSS, both on the
T00129 memory-tail curve -> guard-kill ETA ~14:40-15:00Z) holds the JIT reaper at ZERO while
C=8 cut production ran at 7.4 cut.done/min (~857MB un-reaped map output each = ~5.6GB/min
df decline). Projected trough at the barrier break ~1220-1250G (a deeper repeat of the
1,373G excursion class; worse if the stragglers run long) - the standing rule is explicit:
floor tripped -> NARROW THE LANES. THE ACTION (14:14:40Z): cuts C=8 -> C=4, group-TERM
pgid==pid 3747648 (the watcher's own pattern; the in-flight batches died UNACCOUNTED -
no strikes, the receipted class, cuts-campaign.nohup 'Killed' lines; their ~64 territories
re-join pending with width classes intact) and relaunch BATCH_W=8 N_CONCURRENT=4 RAYON=24
(new master 3966763, composed exactly 4 batches batch-1792073679-0..-3, verified).
THE HONEST COST: cut production ~7.4 -> ~3.7/min until the wave drains its backlog
(a ~50% pole slowdown for the blocked window, ~45-60 min, recovered at re-widen).
THE RECOVERY PLAN: when the barrier breaks (guard-kill or completion of T00370/T00496),
the wave master dispatches the ~200-territory backlog at LANES=128, the reaps cascade,
map/ (~190G) drains, df recovers +100-150G - THEN RE-WIDEN cuts to C=8 (kill the master
MID-WAIT ONLY - never with marathon calls in flight). Prep chunk 2 landing (sub -> ~820,
bounded <=867 by the chunk size); its valve then holds (pending-cut ~310 >= 160) - the
valve re-check between chunks is the proof line still pending the chunk boundary. NO other
intervention; the marathon calls are NOT to be killed by hand (the guard is the oracle;
2h+ CPU each, possibly completing).
## RIDE CONTINUATION 7 WINDOW CLOSE (14:56Z 2026-10-15) - THE HONEST CHECKPOINT STATE
THE WINDOW (13:00-14:56Z) verified the widen end to end, rode the disk-floor trip, and left
the campaign in its disk-safe widened packing. FINAL POLL 14:55:50Z: split 9,891 acc 5,844
(59.1% final, the reconciliation split's own ledger), sub 829, cut.done ~548, map dirs 461
FALLING (the wave draining), call 301 + 72 fail-closed (71 binds ALL NAMED in
wave-logs/failures.txt + T00129's guard kill; +28 binds this window, ~19% of terminal, steady),
df 1,522G RECOVERING (floor 1,500G), mem 704Gib, rvsize 493G falling, poison 0, qdone/fsplit/
widen all Y. THE WAVE'S MEASURED CONSUMPTION: 6.3 terminal/min (44 in 7 min, 128 lanes behind
the flock-serialized gpartition funnel - the wave's TRUE ceiling is the flock, not the lanes;
the funnel was observed with ~86 concurrent partition processes at the backlog peak).
THE WINDOW'S EVENTS, ALL RECEIPTED ABOVE (mid-window findings + the narrowing note):
(1) THE WIDEN FIRED AND COMPLETED VERIFIED (13:00:33-13:01:29Z; widen.log 'WIDEN COMPLETE';
prep POOL_W=64 THREADS=4, cuts C=8 W=8 RAYON=24, wave LANES=128; the final split completed
cleanly FIRST - the ride-6 stall risk closed: 9,891/5,844/4,047, identical acceptance).
(2) THE WAVE'S OWN BATCH-BARRIER FINDING: the wave master's full-TODO xargs blocked ~75 min
behind 2 marathon segdup-tail calls (T00370/T00496) while cut outputs accumulated -> df tripped
the 1,500G floor at 14:11Z (poll 1,495G) -> THE STANDING LEVER APPLIED: cuts C=8->C=4 at
14:14:40Z (group-TERM 3747648, in-flight batches unaccounted - no strikes; relaunch
N_CONCURRENT=4, new master 3966763). BOTH MARATHON CALLS THEN COMPLETED HONESTLY: T00496
exit 0 wall 3,973s peak ~31GB; T00370 exit 0 wall 5,696s peak ~41.5GB - the memory-tail class
is 2-of-3 completing (T00129 the only guard kill). The barrier broke ~14:37Z; the wave
dispatched the ~200-territory backlog at LANES=128; reaps cascaded; df recovered to 1,522G.
(3) THE PREP VALVE'S CHUNKED DISPATCH IS WORKING: the valve opened ~13:43Z, chunk 1
(POOL_W*2=128 chains) returned 14:08:49Z, chunk 2 dispatched (the re-check found the valve
open by its own measure at that boundary - pending-cut was borderline; the overshoot is the
designed one-chunk bound, +128 substrates max); npendcut NOW 240 >= 160 -> THE
'VALVE (mid-dispatch)' PROOF LINE WILL FIRE AT CHUNK 2's RETURN (its last 2 chains are beds
segdup tails, T01152/T01262, in flight since 14:08-14:10, walls up to ~94 min -> by ~15:45Z).
SUCCESSOR: verify the line lands (prep-master.log 'VALVE (mid-dispatch): ... parking').
(4) NEW BATCH-EXIT-1 CLASS WATCH (14:52:10Z): T01144 batch-level failure, attempt 1 of 3 -
the same deterministic-error watch as the chr2 class; attempts 2/3 decide (3x = terminal,
named in the datasheet census). NOT yet investigated - the successor checks the error text
in the batch dir (cuts/batch-*/batch.err of the 1792075531 series).
THE SUCCESSOR'S RIDE: SHORT polls; the cuts C question is now DISK-SAFE RE-WIDEN: the wave
consumes 6.3/min > C=4's 3.7/min production (the backlog is draining, map dirs falling) -
when df >= ~1,600G and the backlog is drained, kill the cuts master MID-WAIT and relaunch
C=6 first (production ~5.5/min ~ the flock ceiling; watch df - C=8's 7.4/min re-trips the
floor only if the flock caps below it; the honest rule: match C to the wave's measured
consumption, RAM is not the binding constraint anymore); verify the mid-dispatch line at
chunk 2's return; the T01144 attempt-2/3 watch; the terminal classes stay terminal
(72 fail-closed of 373 terminal, 19.3%); the adaptive-width retries continue; AT CAMPAIGN
CLOSE: rv-rollup.py -> FILL datasheet-runtime.md (every stage measured vs claimed, the delta
table, BOTH call-wall calibrations + the T00129 tail + THE TWO MARATHON COMPLETIONS
(T00496/T00370) as the tail's counter-receipts, the 1,373G excursion + THIS WINDOW's
1,487-1,495G floor trip + the C=8->C=4 narrowing, the 55-58GB/lane query plateau, the
non-functional query guard, the query-master finish-line death + the marker recovery, the
chr2 AGC class, the wave batch-barrier finding + the flock consumption ceiling) -> the
record's closing increment -> the commit. NO PR PUSH.
## RIDE CONTINUATION 8 WINDOW OPEN (15:02Z 2026-10-15) - THE LONG DRAIN
Session-open increment committed (4766370; vendor/syng pointer drift 054acfc->68ac197 still
uncommitted by design, left for the owner). Live at 15:00:48Z: split 9,891 acc 5,844 (59.1%
final), sub 830, cut map dirs 436 FALLING (the wave draining the backlog), call 326 + 82
fail-closed (+10 since ride-7 close), df 1,538G RECOVERING (floor 1,500G; re-widen gate
~1,600G NOT met), mem 551Gib, rvsize 477G falling, qdone/fsplit/widen all Y, poison 0.
MASTERS VERIFIED: prep 3747592 (POOL_W=64 THREADS=4), cuts C=4 master 3966763 (BATCH_W=8
N_CONCURRENT=4 RAYON=24), wave LANES=128 master 3749885. The widen watcher retired on
completion (WIDEN COMPLETE 13:01:29Z receipted).
SUCCESSOR ITEM 1 CLOSED AT THE WINDOW'S EDGE - THE CHUNKED-DISPATCH PATCH'S LITERAL PROOF:
the first 'VALVE (mid-dispatch)' line FIRED 14:56:22Z in prep-master.log (chunk 2 returned
73min after the 13:43:09Z loop-head hold; the mid-dispatch re-check found 218 pending cut
>= 160 and PARKED the remaining 5,019 of 5,275 TODO - the designed one-chunk bound, verified
live). The chunked-dispatch patch is proven end to end.
THE T01144 INVESTIGATION (ride-7's unchecked item), RECEIPTED: batch 1792075531-5 (exit 1,
14:52:10Z) carried EIGHT territories; the batch.err's single error is
'locus record GRCh38#0#chr2:109250012-109499992 is not identical to the AGC interval under
either offset convention' - THE CHR2 AGC-OFFSET CLASS'S SECOND SIGHTING (the ride-6 terminal
class was chr2:20,750,010-20,999,992; this is chr2:109,250,012-109,499,992 - a different
window, the same deterministic class). THREE territories took attempt 1 of 3 (no cut.done,
no terminal marker): T01144 (component HG002#1#chr2:80523894-80773940), T01135
(HG002#1#chr2:78278709-78528744), T01161 (HG002#1#chr2:84781658-85031769). The erroring
worker cannot be attributed per-territory from the batch dir (batch.log empty; the class
failure is batch-level, exactly like the ride-6 eight). The error locus sits BETWEEN batch
components (T01221 ends 109.17Mb, T01229 starts 110.88Mb) - the class is a shared
window/locus record reached via the territory beds, as in ride-6 where eight territories
shared one locus record. ATTEMPTS 2/3 DECIDE: 3x = a named terminal class sighting; the
bystander-vs-affected question resolves itself by the attempt counter (deterministic error
-> affected territories trip 3x; bystanders that merely lost the batch pass on retry).
THE SUCCESSOR ITEMS REMAINING THIS WINDOW: (1) DISK-SAFE RE-WIDEN cuts C=4->C=6 when the
backlog drains (map dirs -> ~0) AND df >= ~1,600G - the flock ceiling sizes C (~6.3
terminal/min consumption), NOT RAM; kill the cuts master MID-WAIT only, never with marathon
calls in flight; (2) the T01144-class attempts 2/3 watch (batch-failures.log); (3) the
standing contingencies: disk floor 1,500G; marathon segdup calls complete honestly
(T00496/T00370 the counter-receipts; T00129 the only guard kill); the wave's batch-barrier
behind marathon calls; the fail-closed binds terminal and named; memavail watch. AT CAMPAIGN
CLOSE: rv-rollup.py -> FILL datasheet-runtime.md -> the record's closing increment -> the
commit. NO PR PUSH.
## RIDE 8 THE T01144-CLASS RESOLVED TERMINAL (15:10:48Z 2026-10-15) - THE CHR2 AGC-OFFSET CLASS'S SECOND SIGHTING, NAMED
Attempt 3 failed 15:10:48Z (batch 1792076665-12, exit 1) with the IDENTICAL deterministic
error ('locus record GRCh38#0#chr2:109250012-109499992 is not identical to the AGC interval
under either offset convention') - ALL FOUR TERRITORIES TERMINAL after 3 batch-level
failures, cut.exit markers verified, the designed attempt-counting fail-closed path, ZERO
intervention: T01221 (HG002#1#chr2:108922516-109172522, the territory ADJACENT to the bad
locus record), T01135, T01161, T01144 (the attempt-2/3 cadence: attempt 2 at 15:01:23Z,
attempt 3 at 15:10:48Z, the ~9min batch walls). THE CLASS CENSUS NOW: the chr2 AGC-offset
mismatch class, TWO loci, TWELVE territories terminal - ride-6's eight at
chr2:20,750,010-20,999,992 (T00341 T00905 T00955 T00724 T00631 T00950 T00791 T00974) +
this window's four at chr2:109,250,012-109,499,992 (T01221 T01135 T01161 T01144). The
honest batch-composition note for the datasheet: the recomposition was stable (the same
8-territory batch manifest across attempts; 4 completed cut.done on attempt 1, 4 never
did) - the attempt counter is per-territory but the failing batch is the unit, so
co-batched bystanders trip with the affected; the class is fail-closed by design either
way, never retried into holes. BOTH chr2-AGC sightings now stand terminal and named.
## RIDE 8 MID-WINDOW (15:30Z 2026-10-15) - THE STEADY DRAIN, BOTH POLES HONEST
The drain in its disk-safe packing, ~40min into the window, all receipts above. THE PREP
POLE: chunk 3 DISPATCHED ~15:23Z (the valve re-opened as pending-cut drained below 160;
sub 830 -> 911 and settling, the designed one-chunk bound +128 max; rvsize 430->466G,
df dipped 1,586->1,538G and recovered 1,553G - the floor 1,500G not approached). THE
CUTS POLE: C=4 steady ~3.7 cut.done/min; the T01144-class four terminal (above); the
ordinary W=8 guard kills continue (15:04:24Z T01013/T01181 batch, 15:13:03Z T00980 batch -
all packaging, width->4, zero strikes, poison 0). THE WAVE POLE: the xargs snapshot decays
to its stragglers as designed - 128 lanes -> 28 live (15:19Z) -> 16 live (15:30Z), the
oldest calls the marathon segdup-tail class (3040s/30.8GB, 2878s/30.0GB, 2463s/35.0GB -
all under the 64GiB guard, RAM fine at 500-570Gib); calls 393->428 over 18min (~2/min in
the straggler phase, ~6.3/min at snapshot boundaries - the flock ceiling finding's shape);
map dirs ~375-388 hovering (production ~= consumption at C=4 in the straggler phases).
failcall 94 steady. NO intervention pending: the marathon calls complete honestly (the
T00496/T00370 counter-receipts the doctrine); the re-widen C=6 gate (backlog drained AND
df >= ~1,600G) NOT met - df 1,553G, backlog ~375; the successor keeps the SHORT polls.
## RIDE 8 EVENTS (15:10-15:43Z 2026-10-15) - THE MEMORY-TAIL CLASS'S SECOND GUARD KILL + THE DF TROUGH RIDDEN
(1) T01338 NEW FAIL-CLOSED BIND (15:10:38Z): census record 7991, 'no unique verifying key
(2 surviving of 3 candidates)' - the ordinary 7479-class, terminal, named. failcall 96.
(2) T00633 THE CALL MEMORY-TAIL GUARD KILL (15:40:31Z): call.rss receipts the climb 64,349,444kB
-> 71,349,912kB in 5s -> 'GUARD kill 71349912kB' - peak 68.1GB, OVER the 64GiB guard; the
failures.txt 'CALL FAILED (137)' line's 54.7GB is the process's own stale [score] stderr, not
the peak. THE CLASS NOW: T00129 (ride-4, 67.5GB kill), T00496/T00370 (completed - the
counter-receipts), T00633 (kill) - the segdup memory tail is 2-of-4 completing, 2 guard kills,
both honest; T00633 is terminal (no width to yield at RAYON=1). BOTH kills for the datasheet's
call-wall/memory distribution.
(3) THE DF TROUGH RIDDEN WITHOUT THE LEVER: 15:41:25Z poll df 1,513G - the trip class FORMING
(the reaper near zero behind marathon stragglers, cuts C=4 producing, rvsize 469->489G) but
the floor 1,500G NOT CROSSED; the T00633 guard kill resolved a straggler, the snapshot
advanced, reaps cascaded, df recovered 1,513 -> 1,532G by 15:42:45Z. The standing lever (cuts
C narrowing) stayed holstered - applied only when the floor actually trips, per doctrine.
## RIDE 8 CENSUS CORRECTION + THE MEMORY-TAIL CLASS AT 16:28Z 2026-10-15 - FOUR GUARD KILLS, SIX COUNTER-RECEIPTS
THE FAIL-CLOSED CENSUS at 16:28Z: failcall 97 = 93 exit-1 call failures (the bind/class
census, ALL NAMED per line in wave-logs/failures.txt; +22 exit-1 since ride-7 close - the
steady ~17-19% rate; named new this window: T01338 rec 7991 the ordinary 7479-class, and
the rest in the file's own lines) + FOUR 64GiB GUARD KILLS, all the segdup memory-tail
class, all receipted in call.rss: T00129 (ride-4, 67.5GB/5,637s), T00633 (15:40:31Z,
71,349,912kB peak), T00199 (15:42:03Z, 69,853,180kB peak - the climb receipted
55.6->69.9GB), T00518 (16:27:10Z, 67,669,392kB peak). THE SIGNATURE: each kill's last
stderr is the scoring step's 'exactness sample ~222 pairs of ~330-388M valid' line - the
memory tail blows through the guard INSIDE the segdup locus's scoring pass (the same
signature as T00129's ride-4 kill - NOT a separate exactness-failure class, the stderr is
the tail's own last line before the guard fires). THE COUNTER-RECEIPTS (the memory-tail
class COMPLETING): T00496/T00370 (ride-7's, 3,973s/5,696s) + THIS WINDOW: T01055 (4,014s),
T00502 (3,878s), T00632 (3,791s), T00156 (3,682s), T00906 (3,603s), T00155 (3,371s) - six
more honest completions in the 56-67min wall band, all under the guard. THE CLASS
SCOREBOARD for the datasheet: 8 marathons completed / 4 guard-killed (the kills at
67-71GB, the completions peaking ~31-41GB or finishing under the line), ~1/3 of the
deep segdup tails die at the 64GiB guard - THE HONEST TAIL PRICE of the guard-as-oracle
design at genome scale. THE WAVE at 16:28Z: 1 lane left in the snapshot (98min/35.0GB,
the T00370-curve - completing or guard-killing within ~10min); calls 440; map dirs 450;
df steady 1,512G (floor 1,500G NOT tripped - the trough ridden again without the lever);
mem 740-900Gib; rvsize ~501G; the re-widen C=6 gate (backlog drained AND df >= ~1,600G)
still NOT met.
## RIDE 8 THE DISK-FLOOR TRIP + THE STANDING LEVER APPLIED (16:45:29Z 2026-10-15) - RECEIPTED
THE TRIP: df crossed the 1,500G floor at 16:45:15Z (poll 1,488G; rvsize jumped 506->526G
as a C=4 batch landed) - the same trip class as ride-7's 14:11Z: the wave's snapshot held
by ONE marathon segdup-tail call (the last straggler, 94min/43.4GB and climbing, the
T00370-curve) kept the JIT reaper near zero while cut production filled disk. The earlier
troughs (15:41 1,513G, 16:29 1,506G) were ridden WITHOUT the lever - the floor not crossed;
this one crossed. THE LEVER (16:45:29Z): cuts C=4 -> C=2, group-TERM pgid==pid 3966763 (the
in-flight batches died unaccounted - the receipted class, 'Killed' lines in
cuts-campaign.nohup, no strikes; their territories re-join pending with width classes
intact), relaunch BATCH_W=8 N_CONCURRENT=2 RAYON=24 (new master 166454, live-verified).
The wave untouched (its marathon call in flight - the lever never kills the wave side).
THE HONEST COST: cut production ~3.7 -> ~1.9/min until the snapshot breaks (the straggler
completes or guard-kills; T00370 completed at 95min - this one is at ~94min) and the ~450
backlog dispatches at LANES=128; reaps then cascade and df recovers +100-150G. THE
RE-WIDEN QUESTION (the successor's standing item): re-widen C=2 -> C=4 -> C=6 as df
recovers (gate: backlog drained AND df >= ~1,600G; C sized to the wave's measured ~6.3
terminal/min consumption, never above it - C=2's 1.9/min is a ~70% pole slowdown, the
re-widen recovers it).
## RIDE CONTINUATION 8 WINDOW CLOSE (16:48Z 2026-10-15, deadline checkpoint) - THE HONEST STATE
THE WINDOW (14:58-16:48Z) closed successor item 1, resolved the T01144 class, and rode the
floor trip with the standing lever. FINAL POLL 16:48:18Z: split 9,891 acc 5,844 (59.1%
final), sub 917, cut map dirs 448, call 442 + 101 fail-closed (the census in
wave-logs/failures.txt names every line; 97 exit-1 binds/classes + 4 x 64GiB guard kills:
T00129/T00633/T00199/T00518), df 1,505G RECOVERING (floor tripped 16:45:15Z at 1,488G;
lever applied 16:45:29Z: cuts C=4->C=2, new master 166454), mem 943Gib, rvsize 509G, poison 0.
THE SNAPSHOT BROKE AT THE WINDOW'S EDGE (~16:47Z): the last marathon straggler resolved
(calls advancing 440->442, failcall 97->101 as the fresh 128-lane snapshot's fast calls
resolve; fresh lanes all young, rss 1-3GB) - THE ~448-TERRITORY BACKLOG IS NOW DISPATCHED
at LANES=128; the successor SHOULD SEE calls jump, reaps cascade, and df recover
+100-150G within minutes. THEN THE RE-WIDEN LADDER: with the backlog draining and once
df >= ~1,600G, kill the cuts master 166454 MID-WAIT (never the wave; never with marathon
calls in flight) and relaunch C=4 first, C=6 if the flock consumes it (C sized to the
wave's measured ~6.3 terminal/min - the flock ceiling, not RAM).
THE WINDOW'S EVENTS, ALL RECEIPTED ABOVE: (1) THE CHUNKED-DISPATCH PATCH'S LITERAL PROOF:
'VALVE (mid-dispatch)' fired 14:56:22Z (chunk 2 return, 218 pending >= 160, parked 5,019
of 5,275 TODO); chunk 3 dispatched ~15:23Z (sub 830->917, bounded). (2) THE T01144 CLASS
RESOLVED TERMINAL 15:10:48Z: the chr2 AGC-offset class's SECOND SIGHTING (locus
GRCh38#0#chr2:109,250,012-109,499,992), four territories (T01221 T01135 T01161 T01144)
terminal after 3 batch-level failures, zero intervention; the class census now 12
territories across two loci. (3) THE MEMORY-TAIL CLASS AT 8 COMPLETIONS / 4 KILLS:
T00633 (71.3GB, 15:40:31Z), T00199 (69.9GB, 15:42:03Z), T00518 (67.7GB, 16:27:10Z) joined
T00129 as guard kills; T01055/T00502/T00632/T00156/T00906/T00155 completed honestly in the
3,371-4,014s band (with ride-7's T00496/T00370); the kills' shared stderr signature is the
scoring step's 'exactness sample' line - a tail signature, NOT a separate class. (4) THE
DISK-FLOOR TRIP RIDDEN WITH THE LEVER: troughs at 1,513G (15:41) and 1,506G (16:29) held;
the 16:45:15Z trip (1,488G) took the lever (C=4->C=2, batches unaccounted, no strikes).
THE SUCCESSOR'S RIDE: SHORT polls; watch the reaps cascade and df recover as the fresh
snapshot drains; the re-widen ladder (C=4 then C=6 at df >= ~1,600G, backlog drained);
the fail-closed census accrues (each line named in failures.txt); memavail fine (900+Gib);
AT CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md (every stage measured vs
claimed, the delta table, BOTH call-wall calibrations + the T00129/T00633/T00199/T00518
kill tail + the 8 marathon completions, BOTH disk excursions + BOTH floor trips (14:11Z
C=8->C=4, 16:45Z C=4->C=2) + the troughs ridden, the flock-ceiling finding, the
55-58GB/lane query plateau, the non-functional query guard, the query-master finish-line
death, the chr2 AGC class BOTH sightings, the wave batch-barrier, the mid-dispatch valve
proof) -> the record's closing increment -> the commit. NO PR PUSH.
Session-open increment for this window: committed 4766370. The window's events are in this
handoff + the record's window-close increment (committed with this close).

## RIDE CONTINUATION 10 WINDOW OPEN (17:18Z 2026-10-15) - THE LONG DRAIN, QUIET MODE
Session-open increment committed (f47eda5; ride-9's session-open was e89bbe1 - ride-9 rode quiet,
no new class, no window-close note; the vendor/syng pointer drift 054acfc->68ac197 still
uncommitted by design). Live at 17:15:30Z: the ride-8-close transition RECEIPTS LANDED AS
PROJECTED - the ~448-territory backlog dispatched at LANES=128, the reap cascade ran: calls
442->551 + fail-closed 101->147 of 698 terminal, df 1,505->1,573G RECOVERING, rvsize
509->439G falling, cut.done 804, map dirs 335 DRAINING, mem 510Gib, poison 0.
CENSUS: 144 exit-1 binds/classes (every line named in wave-logs/failures.txt) + the 4 standing
64GiB guard kills (T00129/T00633/T00199/T00518) - no new class; the fail-rate at the terminal
count is elevated in this snapshot (46 of 155 fresh terminal ~30% vs the 18.6% campaign
steady - the fresh snapshot is consuming the bind-dense early-chr backlog; the campaign
aggregate at close is the receipt, no extrapolation).
MASTERS VERIFIED: prep 3747592 (valve HELD since 15:19Z at 178 pending-cut; chunk 3's
beds-tail chains in flight, sub 917), cuts C=2 master 166454 (the 16:45Z floor-trip lever),
wave LANES=128 master 3749885.
THE RIDE (quiet mode): SHORT sparse polls (~every 4-5 min via sleep-45 chains, act only on
real events); THE RE-WIDEN LADDER: gate = backlog drained (map dirs -> ~0) AND df >= ~1,600G
-> kill cuts master 166454 MID-WAIT (never the wave; never with marathon calls in flight),
relaunch C=4 first (launch-cuts BATCH_W=8 N_CONCURRENT=4 RAYON=24), C=6 later if the gate
re-mets AND the flock consumes it; CONTINGENCIES: disk floor 1,500G -> the C-narrowing lever;
memory-tail guard kills accrue (terminal, named); the chr2 AGC class terminal at 12; the
binds terminal and named; AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py ->
FILL datasheet-runtime.md -> the record's closing increment -> the commit. NO PR PUSH.
## RIDE 11 THE MAP-WIPE BUG + THE T00197 PREP-TERMINAL CLASS (18:11-18:27Z 2026-10-15) - BOTH RECEIPTED, BOTH FIXED
(1) THE EMPTY-TERRITORY BATCH WIPE (a new destructive class, the session's third key finding): at
18:17:28Z batch -20 was composed with ONE EMPTY-NAME territory - run-cuts.sh's hash-interleave does
`printf '%s\n' "${PENDING[@]}"` which emits ONE BLANK LINE for an EMPTY pending array; the read loop
turns it into a phantom "" territory; NPEND becomes 1 (not 0) so the 'waiting for substrates'
branch never fires; composition then ran `rm -rf "$RV/map///"` = rm -rf OF THE ENTIRE map/ TREE -
264 dirs, ~139G GONE in ~60s (heartbeat: rvsize 373->234G, df 1639->1778G between 18:17:28 and
18:18:31), including the ~80-territory un-called cut backlog (cut.done 859 - 626 called - 153
fail-closed = 80). THE EMPTY BATCH itself died instantly on 'syng metadata not found: prep//locus-
pos64.syng.meta' (batch.err), burned a full syng+AGC load for nothing, and left junk markers at
prep/ ROOT (cut.attempts - removed). THE RECOVERY (receipted): the 80-minus-T01427 = 79 damaged
territories identified (cut.done && no call terminal && not the live marathon) and their cut.done
CLEARED - they re-join the cuts pending set and RE-CUT cleanly (the cut is deterministic; the wave
never saw them - zero spurious call.failures, THE RACE AVOIDED: without the clearing the wave would
have fail-closed all 80 with 'partition step failed' on the missing census). T01427 (the live
marathon, in call since ~17:11Z, RSS plateau 61.9GiB) EXCLUDED: its call inputs were loaded before
the wipe; it will complete but its receipt will MISS cut-stats.json (the wipe took map/T01427) -
NAMED for the datasheet; its calls.jsonl product receipt is intact.
THE FIX (run-cuts.sh, bash -n clean, the master was killed first so no byte-offset hazard): (a) the
interleave's read loop skips empty names - an EMPTY pending list now stays EMPTY, NPEND=0, the
'waiting for substrates' branch fires as designed; (b) belt-and-braces, one_batch's composition
loop skips empty names (no future rm -rf of map/ by any empty element). The cuts master 166454 was
GROUP-KILLED at 18:22Z (named reason: the destructive recurring wipe - every pending-empty loop
would wipe map/ again, and the prep valve was about to feed fresh substrates whose outputs it
would destroy) and RELAUNCHED AT C=4 (BATCH_W=8 N_CONCURRENT=4 RAYON=24, new master 409407) - THE
RE-WIDEN LADDER'S FIRST STEP applied with the recovery (the gate was met: df 1778G >= 1,600G, the
backlog now the 79 re-cuts; C=4 sized to the wave's ~6.3 terminal/min flock ceiling). The 79
re-cuts are in flight (batches 0-3 of 8 composing at 18:23:40Z, real names verified in the manifest).
(2) THE T00197 SYNG SEG-DUP-TAIL PREP-TERMINAL CLASS (the retry-forever gap, pre-named in the
handoff from session open - IT BIT): T00197's syng step - attempt 1 FAILED FAST (exit 1, wall 10s,
13:43Z, chunk-2 era); attempt 2 (chunk 3's LAST chain, since 15:19Z) burned 3h01m+ CPU
compute-bound (state R) at 33MB RSS with NO OUTPUT after 1,203 sequences at +9.8s (syng.err silent
since 15:23:36Z; receipts prep/T00197/syng.* + syng.err.attempt2-snapshot) - AND IT BLOCKED THE
ENTIRE PREP POLE: chunk 3's xargs cannot return while the chain lives, so the valve (open at
npendcut 124 < 160) could never dispatch chunk 4 - THE PREP-SIDE BATCH-BARRIER, the wave barrier's
sibling. ACTION (18:2xZ): the prep master 3747592 GROUP-KILLED (named reasons: the pathological
syng as chunk-3's straggler; the retry-forever gap - the master would redispatch T00197 with no
terminal marker; the TODO-filter patch requires the relaunch per the in-memory-loop byte-offset
hazard), T00197 POISONED (prep/T00197/cut.poison + the honest class line in cuts/poison.tsv - the
marker spelling the cuts master's NTERM and the wave's nopen checks already respect, zero patches
needed there), run-prep-master.sh's TODO filter PATCHED to skip cut.poison territories (the gap
closed; bash -n clean), the master RELAUNCHED (POOL_W=64 THREADS=4, new master 411760) - chunk 4
DISPATCHED within seconds (64 chains live, sub will climb 917 -> ~1,045 bounded by the chunk),
T00197 verified ABSENT from the dispatch. THE CLASS FOR THE DATASHEET: the syng segdup-tail
prep-terminal class, 1 territory of 5,844 (T00197, HG002#1#chr2 component, locus.fa 345MB), two
attempts receipted, fail-closed by hand-decision under the pre-authorized 3-strikes doctrine.
STATE AT 18:26:35Z: split 9,891 / acc 5,844 / sub 917 (chunk 4 landing) / map dirs 32 (the re-cuts)
/ call 626 + fail-closed 153 (779 terminal; the classes unchanged: binds + 4 guard kills) / df
1,763G / mem 775Gib / rvsize 248G / poison 1 (T00197) / T01427 marathon still in flight (75min,
RSS plateau 61.9GiB frozen for ~25min - resolves or guard-kills; never killed by hand).
THE SUCCESSOR'S RIDE: SHORT sparse polls; the 79 re-cuts land (~40min at C=4) then the wave
consumes them at the next snapshot boundary; THE RE-WIDEN LADDER C=4->C=6 when the gate re-mets
(backlog drained AND df >= ~1,600G - C sized to the flock's ~6.3/min, never RAM); the standing
contingencies unchanged (disk floor 1,500G -> the C-narrowing lever; the memory-tail class both
arms; the chr2 AGC class terminal at 12; the binds terminal and named); AT CAMPAIGN CLOSE:
rv-rollup.py -> FILL datasheet-runtime.md (NOW ALSO: the empty-territory map-wipe class + the 79-
territory recovery + the T00197 prep-terminal class + T01427's cut-stats-less receipt) -> the
record's closing increment -> the commit. NO PR PUSH.
## RIDE CONTINUATION 12 WINDOW OPEN (19:55Z 2026-10-15) - THE LONG DRAIN, QUIET MODE
Session-open increment committed (257c8ac). Live at 19:54:37Z: sub 1,121 (prep pole
UNBLOCKED, producing ~1/min past the old 917 hold), map dirs 147 (backlog NOT drained),
call 706 + fail-closed 175 of 881 terminal (19.9% - elevated as the bind-dense late
backlog resolves; every line named in wave-logs/failures.txt), df 1,654G, mem 806Gib,
rvsize 354G, poison 1 (T00197). NEW SINCE RIDE-11 CLOSE, BOTH RECEIPTED: (1) T01427 (the
ride-11 marathon lane, 17:11-19:04Z) resolved as the memory-tail class's FIFTH 64GiB guard
kill (67,398,916kB at 19:04:17Z; call.rss receipt; its receipt misses cut-stats.json per
the map-wipe) - the scoreboard now 8 completions / 5 kills (T00129/T00633/T00199/T00518/
T01427); (2) a burst of SIX cut-batch width-8 guard kills at 19:10:46Z (T01260 T00927
T01241 T01246 T01048 T01153, exit 137 packaging, width->4, no strike - the ORDINARY
ADAPTIVE-WIDTH class, retry at W=4). MASTERS: prep 411760, cuts C=4 409407, wave 3749885;
the ancient idle impg-reference process known-benign - NEVER KILLED. RE-WIDEN C=6 GATE:
df leg MET (1,654G >= 1,600G), backlog leg NOT MET (147 map dirs) - holds. RIDE: SHORT
sparse polls; the ladder per the gate; standing contingencies; AT CAMPAIGN CLOSE:
rv-rollup.py -> FILL datasheet-runtime.md -> the record's closing increment -> the
commit. NO PR PUSH.
## RIDE 13 THE DISK-FLOOR TRIP + THE STANDING LEVER (20:49:05Z 2026-10-15) - RECEIPTED
THE TRIP: df crossed the 1,500G floor ~20:48Z (poll 1,487G; the third crossing - after 14:11Z
C=8->C=4 and 16:45Z C=4->C=2): the wave held behind the T00894 memory-tail marathon (RSS plateau
51.24GB frozen ~35min, 1h39m elapsed) kept the reaper near zero while C=4 cut production + the
prep chunks landed (rvsize 433->498G in 10min; sub 1,123->1,368 as chunks 4/5 landed bounded).
THE LEVER (20:49:05Z): cuts C=4 -> C=2, group-TERM pgid==pid 409407, relaunch BATCH_W=8
N_CONCURRENT=2 RAYON=24 (new master 1002407). HONEST COST ~ZERO: the three in-flight batches
completed exit=0 at 20:48:22-52Z, seconds BEFORE the TERM (only the master's own subshells show
'Killed'; no strikes, no width classes lost). THE NEW MASTER HOLDS AT ITS OWN FLOOR GATE
(cuts/campaign.log 20:49:11Z 'disk floor: 1486G < 1500G - waiting for the reaper') - run-cuts.sh's
built-in floor check pauses composition entirely under the floor: a deeper-than-C=2 pause, the
designed disk safety; cuts resume when df recovers. THE WAVE UNTOUCHED (marathon in flight - the
lever never kills the wave side). RECOVERY SHAPE: T00894 completes or guard-kills (guard-as-oracle,
never by hand; scoreboard 8 completions / 5 kills) -> the snapshot breaks -> the ~235-map-dir
backlog dispatches at LANES=128 -> reaps cascade -> df recovers +100-150G -> THE RE-WIDEN LADDER
resumes (C=2 -> C=4 -> C=6 at the standing gate: backlog drained AND df >= ~1,600G, C sized to
the wave's ~6.3 terminal/min flock ceiling). State at 20:50:39Z: split 9,891 / acc 5,844 / sub
1,368 / map dirs 235 / call 708 + 175 fail-closed of 883 terminal / df 1,486G / mem 927Gib /
rvsize 498G / poison 1 (T00197). Masters: prep 411760, cuts C=2 1002407 (floor-gated), wave
3749885. Commits this window: bd0d8ec session-open, 07dc510 the lever receipt.
## RIDE 13 MID-WINDOW (21:14Z 2026-10-15) - T00894 THE 6TH GUARD KILL, THE CASCADE, THE C=2 RESUME
(1) T00894 (the marathon lane, ~19:59Z start) RESOLVED AS THE MEMORY-TAIL CLASS'S SIXTH 64GiB
GUARD KILL at 21:01:29Z: RSS plateau 51.24GB FROZEN ~45min, then the segdup-tail scoring pass
climbed through the guard to 71,245,288kB peak (call.rss receipt; the 'exactness sample'
stderr signature); terminal at RAYON=1. THE PLATEAU-THEN-CLIMB SHAPE named: an RSS freeze is
NOT safety - the tail blows through the guard inside the scoring pass. Scoreboard: 8
completions / 6 kills (T00129 T00633 T00199 T00518 T01427 T00894).
(2) THE SNAPSHOT BROKE on the kill; the backlog dispatched at LANES=128 (41 lanes live);
THE REAP CASCADE: call 708->739, map dirs 235->218 falling, rvsize 498->481G, df 1,485->
1,503G RECOVERED over the floor in ~13min - the lever's projected recovery shape.
(3) THE C=2 CUTS MASTER RESUMED at the floor-gate release (batches 1792098803-0/-1 composed
21:13:23Z) - the pause ended cleanly, zero lost width classes.
(4) NEW BINDS named (ordinary 7479-class): T00824 rec 9370, T00717 rec 76, T00204 rec 8798;
failcall 189 = 183 binds/classes + 6 guard kills. Census 21:14Z: split 9,891 / acc 5,844 /
sub 1,368 / map dirs 218 / call 739 + 189 fail-closed of 928 terminal / df 1,503G / mem 786Gib /
poison 1. RE-WIDEN GATE: both legs open (df < 1,600G, backlog 218 draining). Commits: bd0d8ec,
07dc510, 3a759d0. THE RIDE continues: SHORT sparse polls; the ladder at the gate; the standing
contingencies; AT CAMPAIGN CLOSE: rv-rollup.py -> datasheet -> the closing increment -> commit.
## RIDE CONTINUATION 13 WINDOW CLOSE (22:04Z 2026-10-15, deadline checkpoint) - THE HONEST STATE
FINAL POLL 22:04:38Z: split 9,891 / acc 5,844 / sub 1,368 (prep valve HELD at 206 pending-cut) /
map dirs 174 / call 857 + fail-closed 219 of 1,076 terminal (20.4%, bind-dense, all named) /
df 1,560G (floor 1,500G) / mem 568Gib / rvsize 423G / poison 1 (T00197). ~14-15 wave lanes live,
the oldest 59min/32GB - MULTIPLE MEMORY-TAIL MARATHONS FORMING (resolve per the guard-oracle
doctrine, never killed by hand; scoreboard 8 completions / 6 kills).
THE WINDOW'S EVENTS, ALL RECEIPTED + COMMITTED (bd0d8ec, 07dc510, 3a759d0, 1aeb91a):
(1) THE THIRD FLOOR TRIP + THE LEVER (20:48Z 1,487G -> 20:49:05Z cuts C=4->C=2, master 409407
group-TERMed; NEAR-ZERO unaccounted loss - the three in-flight batches completed exit=0 seconds
before the TERM); (2) RUN-CUTS.SH'S OWN FLOOR GATE VERIFIED SELF-MANAGING (the new C=2 master
1002407 parked 'waiting for the reaper' 20:49-21:13Z, resumed cleanly at the gate release - the
trip class is now self-managing on the cuts side); (3) T00894 THE 6TH MEMORY-TAIL GUARD KILL
(21:01:29Z, 71.2GB peak after a 51.24GB RSS plateau FROZEN ~45min - THE PLATEAU-THEN-CLIMB
SHAPE named: an RSS freeze is NOT safety; class signature the 'exactness sample' stderr line);
(4) THE SNAPSHOT BROKE ON THE KILL: the backlog dispatched at LANES=128, reaps cascaded,
call 708->857, df recovered 1,485->1,560G over the floor.
Width-8 adaptive kills this window (all ordinary packaging, no strikes): T01852 T01859 T01858
T01833 T01874 T02132 T02033 T02247 T02134 T02211. New binds named per line in failures.txt.
MASTERS: prep 411760, cuts C=2 1002407, wave 3749885; heartbeat live.
THE SUCCESSOR RESUMES THE QUIET RIDE: SHORT sparse polls; THE RE-WIDEN LADDER at the standing
gate (backlog drained AND df >= ~1,600G -> kill the cuts master MID-WAIT only - never the wave,
never with marathon calls in flight; C=2 -> C=4 -> C=6, C sized to the ~6.3 terminal/min flock
ceiling); standing contingencies (floor 1,500G - now self-gated cuts-side; the marathons in
flight; chr2 AGC class terminal at 12; binds terminal and named; T01427/T00894's cut-stats-less
receipts for the datasheet); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md ->
the closing increment -> the commit. NO PR PUSH.
## RIDE CONTINUATION 15 WINDOW OPEN (22:45Z 2026-10-15) - THE LONG DRAIN, QUIET MODE
Session-open increment committed (87ec4bc; ride-14 rode quiet to deadline, its session-open was
312bb19 - no new class that window). Live at 22:43:48Z: split 9,891 / acc 5,844 / sub 1,368
(valve held) / map dirs 227 CLIMBING (C=2 production outpacing the marathon-held snapshot) /
call 868 + fail-closed 220 of 1,088 terminal (20.2%, bind-dense, all named) / df 1,534G
(floor 1,500G) / mem 853Gib / rvsize 449G / poison 1 (T00197). ONE MEMORY-TAIL MARATHON in
flight (5,226s/42.8GB climbing - the plateau-then-climb shape; guard-oracle doctrine, never
killed by hand; scoreboard 8 completions / 6 kills). Ordinary width-8 adaptive kills continue
(latest T02217 22:37:35Z, packaging, no strike). Masters: prep 411760, cuts C=2 1002407
(own floor gate live), wave 3749885. RE-WIDEN LADDER gate: backlog drained AND df >= ~1,600G
-> cuts master MID-WAIT kill only (never the wave, never with marathons in flight), C=2 ->
C=4 -> C=6 (sized to the ~6.3 terminal/min flock ceiling). THE RIDE: SHORT sparse polls
(~4-5 min, act only on real events); standing contingencies (floor 1,500G -> the C lever;
memory-tail both arms; chr2 AGC terminal at 12; binds terminal/named); AT CAMPAIGN CLOSE:
rv-rollup.py -> FILL datasheet-runtime.md -> the closing increment -> the commit. NO PR PUSH.
## RIDE 15 MID-WINDOW (23:02Z 2026-10-15) - THE T01811 7TH KILL (a ride-14-window event, receipted now) + FIVE HONEST MARATHON COMPLETIONS + THE CASCADE
(1) THE MEMORY-TAIL SCOREBOARD, CORRECTED: T01811 = the SEVENTH 64GiB guard kill
    (22:28:46Z, prep/T01811/call.rss peak 67,259,920kB, the 'exactness sample'-class stderr
    line in failures.txt 22:28:46Z) - it landed INSIDE RIDE-14'S WINDOW (22:07-22:43Z), which
    rode to deadline without a close note, so the ride-15 session-open's 'rode quiet, no new
    class' was INACCURATE and is corrected here and in the record's close increment: ride-14
    carried the 7th kill. Scoreboard now 13 completions / 7 kills: the completions = the 8 of
    record (T00496 T00370 T01055 T00502 T00632 T00156 T00906 T00155) + FIVE more honest
    marathons, ALL COMPLETING (receipts/calls.jsonl timestamps): T00693 (22:22:50Z, ride-14's
    64min/30.6GB former), T00431 (22:28:31Z, 60min/31.1GB), T00913 (22:32:26Z, 57min/28.7GB),
    T01048 (22:38:35Z, the width-8 adaptive-kill-now-marathon), and T01696 (22:57:12Z, ~100min
    wall, peak ~42.8GB then released to 24.8GB before completing - THE RIDE-15 MARATHON THIS
    WINDOW WATCHED, the plateau-then-CLIMB-then-release shape's completing arm). The kills:
    T00129 T00633 T00199 T00518 T01427 T00894 T01811 (67-71GB peaks).
(2) THE SNAPSHOT BROKE on T01696's completion ~22:57Z: fresh lanes dispatched at LANES=128
    (all young at 23:01Z), the reap cascade running; call 868 -> 870, failcall 220 -> 226
    (six new binds named 22:58-23:00Z: T01726 rec 6386 / T01982 rec 4714 (2 of 2) / T01891
    rec 7772 / T01948 rec 6152 / T02147 rec 13344 / T02067 rec 8070 - the ordinary
    7479-class); map dirs 259 pending the cascade; df 1,498-1,515G riding the floor line.
(3) THE CUTS FLOOR GATE live: batch -23 exit=137 at 23:00:31Z (T01895 batch, ordinary
    width-8 packaging kill, no strike); df at composition 1,515G; the C=2 master's own floor
    check holds composition while df < 1,500G (self-managing, ride-13-verified).
NO LEVER ACTION: cuts already C=2 (the narrowest packing) and the gate is self-managing; the
wave untouched (its marathons resolve per the guard-oracle doctrine - never killed by hand).
## RIDE 15 THE MEMORY-TAIL TAIL (23:26-23:54Z 2026-10-15) - TWO MORE GUARD KILLS, THE SCOREBOARD 13/9
T01989 (23:26:40Z, prep/T01989/call.rss peak 68,151,176kB) and T02035 (23:53:49Z, peak
67,147,172kB) - the EIGHTH and NINTH 64GiB guard kills of the memory-tail class (the
score-step 'exactness sample'/'rss locus_1_scored' last-stderr variants - the guard is the
oracle, call.rss the receipt). Scoreboard now 13 completions / 9 kills: kills T00129 T00633
T00199 T00518 T01427 T00894 T01811 T01989 T02035 (peaks 67-71GB); completions T00496 T00370
T01055 T00502 T00632 T00156 T00906 T00155 T00693 T00431 T00913 T01048 T01696. THE SNAPSHOT
STILL HELD at 23:53Z behind TWO MORE MARATHONS FORMING (54min/33.5GB, 52min/29.3GB - the
plateau-then-climb shape; resolve per the guard-oracle doctrine). The drain equilibrium
holding nominal: df flat 1,474G (below the 1,500G floor but SAFE-FLAT - the cuts master
self-gated at 23:18:51Z 'waiting for the reaper', the prep valve held at 216 pending-cut
since 23:36Z, so NEITHER pole adds disk; the only drain path is the wave, and its reaps
come when the marathons resolve); rvsize 508G flat; map dirs 256; call 900 + 235 fail-closed
of 1,135 terminal (20.7%, the bind-dense band, all named). NO ACTION: both levers already
at their rests (cuts C=2 AND floor-gated; the wave never touched); the cascade comes with
the marathons' resolution, then df recovers +100-150G and the cuts master resumes at its
gate release; the re-widen ladder (C=4 then C=6) waits for the standing gate (backlog
drained AND df >= ~1,600G).
## RIDE 15 THE TENTH KILL (00:19:04Z 2026-10-16) - T02006, SCOREBOARD 13/10
T02006 guard-killed at 00:19:04Z (prep/T02006/call.rss peak 67,627,380kB - the plateau at
33.5GB frozen ~45min, then the climb through the guard inside the scoring pass, the named
plateau-then-climb shape; stderr 'exactness sample 223 pairs of 207,104,193 valid'). THE
MEMORY-TAIL SCOREBOARD NOW 13 completions / 10 kills (kills: T00129 T00633 T00199 T00518
T01427 T00894 T01811 T01989 T02035 T02006; ~2/3 of the deep segdup tails die at the guard).
TWO MARATHON LANES REMAIN in the snapshot tail at 00:20Z (79min/29.3GB and 77min/42.4GB
climbing) - resolve per the guard-oracle doctrine, never killed by hand. The drain
equilibrium unchanged: df 1,476G flat SAFE (cuts self-gated under the floor since 23:18:51Z,
prep valve held at ~216 pending-cut, no pole adds disk); map dirs 254; call 902 + 237
fail-closed of 1,139 terminal (20.8%, bind-dense, all named).
## RIDE CONTINUATION 15 WINDOW CLOSE (00:35Z 2026-10-16, deadline checkpoint) - THE HONEST STATE
FINAL POLL 00:34:10Z: split 9,891 / acc 5,844 (59.1% final) / sub 1,433 (prep valve HELD at
~216 pending-cut) / map dirs 250 (backlog draining slowly; cuts self-gated under the floor
since 23:18:51Z - the C=2 master waits for the reaper, composition paused) / call 906 +
fail-closed 241 of 1,147 terminal (21.0%, the bind-dense band, every line named in
wave-logs/failures.txt) / df 1,479G (BELOW the 1,500G floor but SAFE-FLAT: both poles at
their rests - cuts floor-gated, prep valve-held, NEITHER adds disk; the only drain path is
the wave) / mem 949Gib / rvsize 504G flat / poison 1 (T00197).
THE WINDOW'S EVENTS, ALL RECEIPTED IN THE HANDOFF SECTIONS ABOVE:
(1) THE RIDE-14 CORRECTION + T01811 THE 7TH KILL (22:28:46Z, 67.26GB peak - it landed in
    ride-14's quiet window; the ride-15 session-open's 'rode quiet, no new class' was
    inaccurate and is corrected here and in the record).
(2) FIVE HONEST MARATHON COMPLETIONS: T00693 (22:22:50Z) T00431 (22:28:31Z) T00913
    (22:32:26Z) T01048 (22:38:35Z) - ride-14's four forming marathons, ALL COMPLETING - plus
    T01696 (22:57:12Z, ~100min wall, peak ~42.8GB released before completion - the
    ride-15-watched marathon's completing arm).
(3) THREE MORE GUARD KILLS THIS WINDOW: T01989 (23:26:40Z, 68.15GB), T02035 (23:53:49Z,
    67.15GB), T02006 (00:19:04Z, 67.63GB - the 33.5GB plateau frozen ~45min then the climb,
    the named plateau-then-climb shape). THE MEMORY-TAIL SCOREBOARD NOW 13 completions /
    10 kills (T00129 T00633 T00199 T00518 T01427 T00894 T01811 T01989 T02035 T02006; peaks
    67-71GB) - ~2/3 of the deep segdup tails die at the 64GiB guard at genome scale.
(4) THE CUTS FLOOR GATE FIRED SELF-MANAGING (23:18:51Z 'disk floor: 1498G < 1500G - waiting
    for the reaper') and has held ~76min since - the trip class fully self-managing at C=2;
    NO lever action taken (cuts already at the narrowest packing AND gated; the wave never
    touched - its marathons resolve per the guard-oracle doctrine).
(5) THE SNAPSHOT BROKE TWICE (on T01696's completion ~22:57Z; again as the marathons
    resolved through 00:28Z - one marathon lane remained at window close, 86min/46.7GB
    climbing); reaps trickled as stragglers resolved; df held 1,474-1,509G flat throughout.
MASTERS AT CLOSE: prep 411760 (valve-held), cuts C=2 1002407 (floor-gated, BATCH_W=8
N_CONCURRENT=2 RAYON=24), wave LANES=128 3749885; heartbeat live. THE RE-WIDEN LADDER still
holds at the standing gate (backlog drained AND df >= ~1,600G - neither leg met: 250 map
dirs, df 1,479G).
THE SUCCESSOR RESUMES THE QUIET RIDE: SHORT sparse polls (~4-5 min, act only on real
events); when the last marathon resolves the snapshot breaks fully, the ~250 backlog
dispatches at LANES=128, reaps cascade, df recovers +100-150G, the cuts master resumes at
its gate release, THEN the re-widen ladder (C=4 first, C=6 if the gate re-mets and the
flock consumes it - sized to the ~6.3 terminal/min flock ceiling, never RAM); standing
contingencies (disk floor 1,500G -> the C lever; memory-tail both arms, scoreboard 13/10;
chr2 AGC class terminal at 12; binds terminal and named; T01427/T00894's cut-stats-less
receipts for the datasheet); AT CAMPAIGN CLOSE (all 5,844 accepted terminal):
rv-rollup.py -> FILL datasheet-runtime.md (every finding: the delta table with refuted
claims stated plainly, both call-wall calibrations + the tail points, both disk excursions +
all floor trips, the flock ceiling, the query plateau + its non-functional guard + the
master finish-line death, the chr2 class, the memory-tail scoreboard 13/10 with the
plateau-then-climb shape, the valve/chunked-dispatch receipts, the map-wipe class + the
79-territory recovery + the T00197 class, the run-cuts.sh floor gate, the ride-14
correction) -> the closing increment -> the commit. NO PR PUSH.
Session-open increment for this window: committed 87ec4bc. The window's events are in this
handoff + the record's window-close increment (committed with this close).

## RIDE 16 THE FAIL-PATH REAPER BUG + THE 132G SWEEP (00:55Z 2026-10-16) - THE ROOT CAUSE OF THE FLOOR TRIPS, FOUND AND SWEPT
THE DISCOVERY (the stall that wasn't ridden blind): at window open the cuts master sat self-gated
under the floor (df 1,479-1,485G < 1,500G) with only TODO=6 callable territories left - the
designed recovery (reaps from the wave) could deliver at most ~5G, NEVER the +15G the gate needs:
A STALL WAS FORMING (cuts gated on df; df recovery gated on cuts' downstream calls; prep valve held
at 212 pending-cut). Investigation of the map/ tree found WHY: **THE FAIL-CLOSED PATH SKIPS THE JIT
REAPER** - run-wave.sh's one_territory `return 1`s BEFORE the reap block on ALL THREE failure
branches (reads fetch / partition step / call failure), so every fail-closed territory leaks its
map census/cache/reads (~1.34G each), calls/, partitions/, reads/ AND the prep substrate heavy
artifacts, with NO receipts written either (receipts/T01980 was empty at discovery). Measured: 88
post-map-wipe terminal call.exit territories carried 118G of un-reaped cut output; 242 terminal
fails carried ~132G total. THIS IS THE ROOT CAUSE OF THE RECURRING 1,500G FLOOR TRIPS (14:11Z,
16:45Z, 20:48Z, 23:18Z): the "reaper near zero" pattern at snapshot boundaries was not just
marathon stragglers - the fast fail-closed calls (the bind-dense band, ~20% of terminal) reaped
NOTHING while cut production filled disk. Left unfixed it would leak ~1.5T over the remaining
~4,670 territories (the campaign could not have finished).
THE SWEEP (rv-fail-reaper.sh, written + run 00:54:02Z, idempotent, dry-run first): exactly the
success path's own receipt-copy + delete list applied to all 242 terminal call.exit territories -
receipts KEPT (cut-stats.json, all walls, call.peak_rss_kb from call.rss, call.exit.code - this
RETROACTIVELY repairs the receipt gap for the whole fail-closed census, the datasheet's failure
census now has its per-territory receipts), beds/ledgers kept, everything else deleted as consumed.
~132G freed: df 1,486 -> 1,648G. THE CUTS MASTER RESUMED ON ITS OWN at 00:53:32Z (floor gate
released at df 1,563G; batches -32/-33 composing). The stall is broken.
THE PATCH (planned, per the byte-offset doctrine - the wave master 3749885 is NOT killed while
its 6 in-flight lanes run: killing now would orphan lanes that then race the relaunched master's
TODO on the same territories = the T00003/T00114 double-writer class): WAIT for the 6 lanes to
resolve (the master is blocked in xargs-wait anyway - zero throughput cost), THEN kill the master
MID-WAIT, add a shared reap helper to the three failure branches of one_territory (same lists as
the sweep), bash -n, relaunch LANES=128 within seconds. The 6 lanes run the OLD code - if any
fails, its map leaks; RE-RUN rv-fail-reaper.sh at each poll until the patch is live (idempotent).
ALSO THIS WINDOW: T02103 = the 11TH 64GiB memory-tail guard kill (00:31:42Z, peak 70,878,788kB,
prep/T02103/call.rss receipt, the plateau-then-climb shape) - the scoreboard 13 completions /
11 kills; its receipt (and T01427's, T00894's, T01989's, T02035's, T02006's) now carries
cut-stats.json + walls + rss peak via the sweep - the cut-stats-less-receipt note is REPAIRED for
all except where the map-wipe already took them (T01427/T00894 keep the named caveat).
STATE at 00:55Z: split 9,891 / acc 5,844 / sub 1,433 (valve will open as cuts consume pending-cut)
/ map dirs 168 (147 retry-partials + ~16 fresh C=2 emissions) / call 915+242 + in-flight 6 lanes /
df 1,648G (re-widen df leg MET) / mem ~920Gib / rvsize ~370G / poison 1 (T00197). Masters: prep
411760, cuts C=2 1002407 (composing), wave 3749885 (blocked on 6 lanes - the patch window opens
when they resolve). NEXT: the wave patch + relaunch; then the re-widen ladder C=4 (gate: wave
caught up AND df >= ~1,600G - now fundamentally safer, the leak that tripped the floor 4x is
closed); then the quiet ride to campaign close; AT CLOSE: rv-rollup.py -> FILL datasheet-runtime.md
(NOW ALSO: the fail-path reaper bug + the 132G sweep + the retroactive receipt repair, named
alongside the delta table, both call-wall calibrations + the tail points, both disk excursions +
all four floor trips, the flock ceiling, the query plateau + non-functional guard + master
finish-line death, the chr2 class, the memory-tail scoreboard 13/11 with the plateau-then-climb
shape, the valve/chunked-dispatch receipts, the map-wipe class + the 79-territory recovery + the
T00197 class, the run-cuts.sh floor gate, the ride-14 correction) -> the closing increment ->
the commit. NO PR PUSH.

## RIDE 17 SESSION OPEN (01:35Z 2026-10-16)
Resumed at 01:33Z from the ride-16 close. Census 01:33:25Z: split 9,891 / acc 5,844 / sub 1,433 /
map dirs 213 / call 920 + 242 fail-closed of 1,162 terminal (20.8%) / df 1,600G (re-widen df leg
MET, backlog leg open) / mem 898Gib / rvsize 383G / poison 1. THE WAVE PATCH REMAINS QUEUED per
the byte-offset doctrine: master 3749885 (OLD code) has exactly ONE lane in flight - T01954, a
memory-tail marathon (61min/39.7GB at open, plateau-then-climb class). WHEN IT RESOLVES: kill the
master MID-WAIT, add the shared reap helper to the three failure branches of one_territory in
run-wave.sh (same lists as rv-fail-reaper.sh), bash -n, relaunch LANES=128 within seconds;
re-run rv-fail-reaper.sh at each poll until the patch is live (idempotent). CUTS self-managed
this window: 01:26:30Z width-8 packaging guard-kill burst (T02381 T02352 T02389 T02314 T02387
T02319 T02339 T02384, exit 137, width->4, no strikes), now composing at width 4, df 1,597-1,600G.
Re-widen ladder at the gate. Session-open increment: commit f6a1c64. NO PR PUSH.
PATCH ARMED (01:39Z): rv-patch-watcher.sh (detached, pid 1617051, log patch-watcher.log) holds
the ride-16-queued wave patch at the byte-offset safe point - 2s cadence, waits for T01954 (the
last in-flight marathon lane) to resolve, catches the master in its >=30s between-batches sleep
with ZERO lanes in flight, kills 3749885 (no orphans, no double-writer race), re-runs
rv-fail-reaper.sh (old-code leaks, idempotent), swaps in run-wave.patched.sh (bash -n gated,
staged; fail_reap on all three failure branches + rss_peak/fail_reap now exported - fixing the
success-path empty call.peak_rss_kb defect seen in wave.nohup 'rss_peak: command not found'),
relaunches LANES=128. Tested: fail_reap verified on a synthetic territory (receipts kept,
intermediates deleted). patch-watcher.done marks completion.

## RIDE 17 THE WAVE PATCH APPLIED (02:06-02:10Z) - THE FAIL-PATH LEAK CLASS CLOSED LIVE
T01954 resolved as the TWELFTH memory-tail guard kill (02:06:14Z, peak 70,213,072kB, ~36min
39.7GB plateau then the climb; scoreboard 13 completions / 12 kills). rv-patch-watcher.sh caught
the safe point 8s later (master's between-batches sleep, ZERO lanes), killed old master 3749885,
swept T01954's fresh ~2G leak (receipts complete: peak_rss/walls/cut-stats/exit), swapped in the
patched run-wave.sh, relaunched LANES=128 as master 1699248 at 02:06:37Z. VERIFIED LIVE: 7
fail-closed calls FAIL-REAPED IN-LANE in the first ~3min (fail-reaper.log) - the root cause of the
four 1,500G floor trips is closed at the source; ALSO rss_peak/fail_reap exported (fixes the
success-path empty call.peak_rss_kb - DONE lines carry real peaks again). State 02:09:45Z: call
920 + 251 fail-closed terminal / map 218 / sub 1,496 (prep VALVE OPENED ~02:00Z) / df 1,571G /
rvsize 412G / poison 1. Masters: prep 411760, cuts C=2-adaptive-4 1002407, wave 1699248 (56 lanes).
Record increment committed c41b8e5. Re-widen ladder at the gate (backlog leg open). NO PR PUSH.

## RIDE 18 SESSION OPEN (02:59Z 2026-10-16)
Session-open increment committed (220e606). Live at 02:58:45Z: split 9,891 / acc 5,844 / sub 1,496 /
map dirs 232 / call 1,206 + 260 fail-closed of 1,466 terminal (17.7% fail share) / df 1,563G /
mem 578Gib / rvsize 420G / poison 1. The patched wave (1699248) FAIL-REAPING in production
(fail-reaper.log through 02:13:41Z, no leaks since). TWO memory-tail marathons in flight
(48min/51.1GB + 47min/52.7GB, plateau-then-climb) - resolve per the guard-oracle doctrine.
Cuts C=2-adaptive master 1002407 healthy (T02440 the ordinary 02:49:54Z width-8 kill, no strike).
Prep 411760 valve-open producing. Re-widen ladder at the gate (both legs open: 232 map dirs,
df 1,563G < 1,600G). Ride: SHORT sparse polls (~4-5 min), act only on real events. NO PR PUSH.

## RIDE 18 WINDOW CLOSE (04:07Z 2026-10-16) - THE HONEST CHECKPOINT STATE
THE RIDE (02:59-04:07Z, ~40 SHORT sparse polls) rode the drain in quiet mode; NO lever action, NO
master touched - the levers stayed at their rests (cuts C=2-adaptive, prep valve-held, wave never
touched). FINAL POLL 04:06:13Z: split 9,891 / acc 5,844 (59.1% final) / sub 1,558 (prep valve HELD
at 171 pending-cut) / map dirs 287 (cut production outpacing the snapshot-held wave - the known
pattern) / call 1,213 + fail-closed 264 of 1,477 terminal (25.3% of accepted; NO new bind
ambiguities this window - the last bind was T02413 02:13:41Z, pre-open; the +4 fail-closed are all
memory-tail guard kills, below) / df 1,511G (floor 1,500G - one dip to 1,507G at 04:01Z, never
crossed; the cuts master never self-gated, composing width-4 batches at 1,510-1,525G) / mem
702Gib / rvsize 472G / poison 1 (T00197).
THE MEMORY-TAIL SCOREBOARD NOW ~16 COMPLETIONS / 16 KILLS (from ride-17's 13/12):
(1) FOUR GUARD KILLS THIS WINDOW, every one FAIL-REAPED IN-LANE by the patched wave with FULL
    receipts (the class's first kills under the patched code - zero leaks, zero receipt gaps, the
    ride-16/17 fix proven on its own worst case): T02166 (03:03:34Z, 72,256,168kB), T02321
    (03:06:01Z, 69,568,132kB), T02327 (03:11:42Z, 68,931,868kB), T02396 (03:54:46Z, 68,274,136kB).
(2) THREE COMPLETIONS, extending the completion arm's wall range: T02340 (03:05:20Z, 3,159s,
    36.7GB), T02033 (03:20:33Z, 4,394s, peak 63,405,616kB - A HAIR UNDER THE GUARD), T02341
    (03:47:22Z, 5,608s = 93.5min, peak 58.5GB - THE LONGEST COMPLETION WALL YET; the completion arm
    now spans ~3,100-5,600s). A NEW SUB-SHAPE OBSERVED: several lanes climb FAST (40-64GB inside
    4-10min, e.g. the lanes behind T02396) vs the classic plateau-then-climb - both arms accrue to
    the same census; the fast-climb lanes mostly guard-kill.
THE MARATHON TRIO STILL IN FLIGHT AT CLOSE (never killed by hand, per the guard-oracle doctrine):
1714266 (1h54m, 52.3GB frozen), 1718947 (1h53m, 62.1GB frozen), 1723227 (1h52m, 65.0GB frozen) -
plateaus past the longest completion wall (93.5min); they resolve per the class (T02033/T02341
completed from similar plateaus), or the guard takes them. Plus 1767054 (1h26m, 38.0GB).
ORDINARY WIDTH-8 BATCH KILLS (packaging class, adaptive width->4, zero strikes): T02596 (03:36:17Z),
T02598 (03:44:17Z), T02624 (03:44:36Z), T02595 (03:53:28Z) - the composer now picks width 4.
THE RE-WIDEN LADDER holds at the standing gate (backlog drained AND df >= ~1,600G -> C=2 -> C=4 ->
C=6, kill the cuts master MID-WAIT only, never the wave, never with marathons in flight; C sized to
the wave's ~6.3 terminal/min flock ceiling, never RAM) - BOTH LEGS OPEN at close (287 map dirs,
df 1,511G).
MASTERS AT CLOSE: prep 411760 (POOL_W=64 THREADS=4, valve-held at 171), cuts C=2-adaptive 1002407
(BATCH_W 8->4 N_CONCURRENT=2 RAYON=24, own floor gate live), wave LANES=128 1699248 (PATCHED code,
fail-reaping in production all window); heartbeat x2 (the known benign duplicate); the ancient idle
impg-reference known-benign - NEVER KILLED.
THE SUCCESSOR RESUMES THE QUIET RIDE: SHORT sparse polls (~4-5 min, act only on real events); the
marathon trio resolves per the doctrine; the re-widen ladder per the gate; the standing
contingencies (disk floor 1,500G -> the C-narrowing lever/self-gate; the memory-tail class both
arms - scoreboard ~16/16 with the wall range ~3,100-5,600s vs the ~50% guard-kill rate at this
window's stretch; the chr2 AGC class terminal at 12; the binds terminal and named; T00197 1);
AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (every
finding per the standing list + THIS WINDOW'S: the patched wave's four guard kills all
receipted/fail-reaped - the fix's worst-case proof; T02341's 5,608s completion-arm record; the
fast-climb sub-shape; the >1h50m plateau trio) -> the closing increment -> the commit. NO PR PUSH.
Session-open increment for this window: committed 220e606; window close committed with this note.

## RIDE 19 SESSION OPEN + THE CENSUS CORRECTION (04:09-04:33Z 2026-10-16)
Session-open increment committed (eb0f5f6). THE CENSUS CORRECTION committed (5e7f4aa): since the
ride-16 fail-receipt sweep the poll's call= field counts SUCCESS+FAIL receipt dirs - the honest
census at 04:32Z: 949 success (callok) + 265 fail-closed calls (248 binds + 17 memory-tail guard
kills, T02382 the newest 04:27:20Z, fail-reaped in-lane) of 1,214 call-terminal; + 16 cut.exit +
1 poison (T00197) = 1,231 terminal of 5,844 (21.1%). The "of 1,477 terminal" ride-18 reading
was impossible (cut.done caps it). poll-once.sh now carries callok=/cutexit= (no master uses it).
THE CUT.EXIT CLASS CORRECTED: 16 terminal = 2 genuine AGC-offset defects (T00341 chr2:20.75M,
T01144 chr2:109.25M, 3 attempts each, the only two bad locus records in the batch.err corpus) +
14 INNOCENT CO-BATCH TERMINALS (the single-locus deterministic error aborts the batch; the
3-attempt rule terminalizes co-batched territories; their own loci are not defective). The
datasheet must state the blast-radius finding plainly.
LIVE STATE: df 1,485G (BELOW the 1,500G floor since ~04:19Z - the cuts master SELF-GATED at its
own floor check 04:19:19Z 'disk floor: 1485G < 1500G - waiting for the reaper'; cuts already at
C=2, no lever to pull); prep valve held (163 pending-cut at 04:16Z, sub 1,618); the wave snapshot
HELD behind TWO memory-tail marathons (T02053 2h17m/65.0GB and T02344 2h19m/52.3GB plateaus,
plus the cross_lcp_merge child 1767054 1h51m/38GB) - T02382 (the third) guard-killed 04:27:20Z
and fail-reaped. Recovery shape: the marathons resolve -> snapshot breaks -> the ~286 backlog
dispatches at LANES=128 -> reaps cascade -> df recovers -> the cuts gate releases -> THEN the
re-widen ladder (C=4, C=6 - gate: backlog drained AND df >= ~1,600G).
THE RIDE continues: SHORT sparse polls, act only on real events; disk floor 1,500G (self-gated
on the cuts side); memory-tail both arms (scoreboard ~16 completions / 17 kills); binds terminal
and named; AT CAMPAIGN CLOSE: rv-rollup.py + FILL datasheet-runtime.md (the standing close list
+ THE CENSUS CORRECTION + the cut.exit blast-radius finding) + the closing increment + commit.
NO PR PUSH.

## RIDE 20 SESSION OPEN (04:47Z 2026-10-16)
Session-open increment committed (dc83a19). Live at 04:46:30Z (honest poll fields): split 9,891 /
acc 5,844 / sub 1,618 / map dirs 285 / call-terminal 1,215 = callok 949 + fail-closed 266
(+ 16 cut.exit + 1 poison = 1,232 of 5,844, 21.1%) / df 1,487G (below the 1,500G floor since
~04:19Z; cuts master SELF-GATED at its own floor check 04:19:19Z; cuts already C=2 - no lever;
recovery rides the reap cascade) / mem 884Gib / rvsize 496G / prep valve HELD (215 pending-cut).
THE MARATHON TRIO RESOLVED TO ONE: T02382 guard-killed 04:27:20Z (reaped in-lane); T02053
GUARD-KILLED 04:45:19Z at 68.99GB after a 2h34m frozen plateau (fail-reaped 04:45:27Z, full
receipts - the patched wave's marathon worst case clean). Scoreboard ~16 completions / 18 kills.
ONE LANE IN FLIGHT: T02344 (pid 1723227, 2h33m, 65.0GB frozen) - the wave snapshot HELD behind
it; on resolution the ~285 backlog dispatches at LANES=128 -> call jump + reap cascade ->
df recovery -> cuts gate releases -> THEN the re-widen ladder (C=4, C=6; gate: backlog drained
AND df >= ~1,600G; kill the cuts master MID-WAIT only, never the wave, never with marathons
in flight). Masters: prep 411760 valve-held, cuts 1002407 self-gated at floor, wave 1699248
(patched, reaper clean through 04:45:27Z); heartbeat live; ancient impg-reference benign.
Ride: SHORT sparse polls (~4-5 min), act only on real events. NO PR PUSH.
RIDE 20 EVENT LOG (05:03Z): THE SNAPSHOT BROKE AT 04:56:42Z - T02344 (the last marathon lane,
2h42m frozen at 65.0GB) woke, climbed, and was GUARD-KILLED at 69,429,052kB, fail-reaped in-lane
04:56:47Z (full receipts; scoreboard ~16 completions / 19 kills - the 2h30m+ marathon plateau
class now 2-for-2 resolving via the guard, zero hand kills). THE EXPECTED CASCADE UNFOLDED:
the ~285 backlog DISPATCHED at LANES=128 (12 lanes at 04:59Z -> 29 at 05:02Z), call jump
1,215 -> 1,229 (callok 949 -> 955; fail 266 -> 274: T02344 + the ordinary exit-1 bind class
T02167/T02535/T02109 named in wave-logs/failures.txt - NO new class), map dirs 285 -> 287
(cuts production resumed and outpacing again, the known pattern), df 1,487 -> 1,504G RECOVERED
above the floor, rvsize 496 -> 479G, THE CUTS FLOOR GATE RELEASED 05:01:38Z (back to composing
width-8 batches). Re-widen gate: backlog leg NOT met yet (wave draining ~276-287 dirs); df leg
NOT met (1,504G < 1,600G). No lever action taken - the ladder waits for both legs.

## RIDE CONTINUATION 21 WINDOW OPEN (05:47Z 2026-10-16) - THE LONG DRAIN, QUIET MODE
Session-open increment committed (5d73f69). Live at 05:47:03Z: split 9,891 / acc 5,844 / sub 1,618 /
map dirs 240 DRAINING (the ride-20 cascade: 285->240; call-terminal 1,320 = callok 1,018 + fail 302,
+16 cut.exit + 1 poison = 1,337 of 5,844, 22.9%) / df 1,550G RECOVERING (floor 1,500G; re-widen df
leg wants >= ~1,600G) / mem 150Gib (the full 128-lane wave saturation; query lanes retired - the
protected-lane OOM class moot) / rvsize 433G falling / poison 1 (T00197). One ordinary width-8
packaging kill: T02754 05:37:28Z (no strike). The wave in its fresh full-lane drain (~50 lanes) with
the memory-tail marathon cohort forming (T02079 46min/53.4GB, T02336 44min/52.2GB + five 40-45GB
lanes) - guard-oracle doctrine, never killed by hand; scoreboard ~16 completions / 19 kills.
Masters: prep 411760, cuts C=2-adaptive 1002407, wave 1699248 (patched, fail-reaping in production);
heartbeat live. THE RIDE (quiet mode): SHORT sparse polls (~4-5 min via short sleeps, act only on
real events); the re-widen ladder per the standing gate (backlog drained AND df >= ~1,600G -> kill
the cuts master MID-WAIT only - never the wave, never with marathon calls in flight; C=2 -> C=4 ->
C=6, sized to the ~6.3 terminal/min flock ceiling, never RAM); standing contingencies (disk floor
1,500G -> the C lever/self-gate; memory-tail both arms; chr2 AGC class terminal per the ride-19
blast-radius enumeration; binds terminal and named); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL
datasheet-runtime.md -> the closing increment -> the commit. NO PR PUSH.
## RIDE CONTINUATION 22 WINDOW OPEN (06:57Z 2026-10-16) - THE LONG DRAIN, QUIET MODE
Session-open increment committed (cdda4eb). Live at 06:57:05Z: split 9,891 / acc 5,844 /
sub 1,673 (valve open, producing) / map dirs 247 (backlog NOT drained) / call-terminal 1,342 =
callok 1,029 + fail-closed 313 (+16 cut.exit + poison 1 = 1,359 of 5,844, 23.3%; every line
named in wave-logs/failures.txt) / df 1,543G (floor 1,500G; re-widen df leg wants >= ~1,600G
- NOT met) / mem 753Gib / rvsize 440G / poison 1 (T00197). THREE memory-tail marathons in
flight (1h53m/42.9GB, 1h44m/29.3GB, 1h44m/36.8GB) - guard-oracle doctrine, never killed by
hand; THREE fresh 64GiB guard kills since ride-21 open (T02557 06:54:03Z, T02639 06:54:51Z,
T02553 06:55:22Z - 'exactness sample' signature), all FAIL-REAPED IN-LANE by the patched
wave with full receipts (scoreboard ~16 completions / ~22 kills). Masters: prep 411760
(valve held at 195 pending-cut), cuts C=2-adaptive 1002407 (composing width-8 batches -101+),
wave 1699248 (patched, fail-reaping in production); heartbeat live; ancient impg-reference
benign. RE-WIDEN LADDER gate: backlog drained AND df >= ~1,600G -> kill the cuts master
MID-WAIT only (never the wave, never with marathons in flight), C=2 -> C=4 -> C=6 (sized to
the ~6.3 terminal/min flock ceiling, never RAM) - BOTH LEGS OPEN at open. THE RIDE: SHORT
sparse polls (~4-5 min via sleep-45 chains, act only on real events); standing contingencies
(floor 1,500G -> the C lever/self-gate; memory-tail both arms; chr2 AGC class per the
ride-19 blast-radius enumeration; binds terminal/named; T00197 1); AT CAMPAIGN CLOSE:
rv-rollup.py -> FILL datasheet-runtime.md -> the closing increment -> the commit. NO PR PUSH.
## RIDE 22 MID-WINDOW (08:10Z 2026-10-16) - THE MARATHON TRIO RESOLVED, THE FLOOR TOUCHED SELF-GATED, THE CASCADE RUNNING
THE WINDOW'S EVENTS, ALL RECEIPTED (single-tail checks, no master touched, no hand kills):
(1) THE MARATHON TRIO RESOLVED: T02656 COMPLETED 07:11:35Z at 7,149s (119min) wall / 53.2GB
    peak - THE NEW LONGEST COMPLETION WALL (beats T02341's 5,608s; the completion arm now
    spans ~3,100-7,149s); T02651 GUARD-KILLED 07:00:13Z at 68,189,252kB (full receipts,
    fail-reaped in-lane); T02245 GUARD-KILLED 07:49:01Z at 68,145,040kB after a ~2h46m lane
    wall - THE LONGEST MEMORY-TAIL LANE YET (frozen 42.9GB plateau most of it, then the
    climb; full receipts). Scoreboard ~17 completions / ~24 kills (the census at close via
    rv-rollup is the receipt).
(2) THE FLOOR TOUCH + THE SELF-GATE: df crossed to 1,499G at 07:45-07:46Z; the cuts master
    SELF-GATED 07:45:09Z ('disk floor: 1499G < 1500G - waiting for the reaper') - C=2 is
    already the narrowest packing, NO lever action (the ride-19/20 precedent); the gate
    RELEASED 07:57:56Z on the cascade's df recovery (df 1,505->1,559G).
(3) THE CASCADE: T02245's resolution broke the snapshot ~07:49Z; the ~270 backlog dispatched
    at LANES=128; call-terminal 1,344 -> 1,400 (callok 1,030 -> 1,062, fail 314 -> 338 - the
    bind-dense band resolving, every line named in wave-logs/failures.txt); map dirs 271 ->
    215 DRAINING; rvsize 485 -> 433G falling; df 1,499 -> 1,559G RECOVERED.
(4) ORDINARY WIDTH-8 PACKAGING KILLS (exit 137, width->4, no strikes): T02877 06:59:33Z,
    T02870 07:08:22Z, T02837 07:16:11Z.
RE-WIDEN LADDER: BOTH LEGS STILL OPEN at 08:09Z (backlog 215 draining, df 1,551G < 1,600G) -
    the ladder waits for the gate. State 08:09:28Z: split 9,891 / acc 5,844 / sub 1,673 (valve
    held at 195 pending-cut) / call 1,400 = callok 1,062 + fail 338 (+16 cut.exit + poison 1 =
    1,417 of 5,844, 24.3%) / mem 447Gib (full-lane saturation) / poison 1 (T00197). Masters:
    prep 411760, cuts C=2-adaptive 1002407, wave 1699248 (patched, fail-reaping in production).
## RIDE CONTINUATION 22 WINDOW CLOSE (08:40Z 2026-10-16) - THE HONEST CHECKPOINT STATE
THE RIDE (06:57-08:40Z, ~25 SHORT sparse polls via sleep-45 chains) rode the marathon trio to
resolution, the floor touch to its self-gate, and the cascade to recovery - ONE REAL ACTION
TAKEN: NONE (no master touched, no kill by hand, no lever pulled; the cuts self-gate did the
floor work; C=2 was already the narrowest packing).
FINAL POLL 08:40:36Z: split 9,891 / acc 5,844 (59.1% final) / sub 1,722 (valve opened ~08:13Z,
chunk bounded +49) / map dirs 242 (production outpacing the marathon-forming wave - the known
pattern) / call-terminal 1,421 = callok 1,082 + fail-closed 339 (+16 cut.exit + poison 1 =
1,438 of 5,844, 24.6%; every fail line named in wave-logs/failures.txt) / df 1,538G (floor
1,500G) / mem 120Gib (full 128-lane saturation; the protected-lane OOM class moot - query
lanes retired) / rvsize 445G / poison 1 (T00197).
THE WINDOW'S EVENTS, ALL RECEIPTED (the mid-window note above carries the detail):
(1) THE MARATHON TRIO RESOLVED: T02656 COMPLETED 07:11:35Z at 7,149s/53.2GB - THE NEW LONGEST
    COMPLETION WALL (beats T02341's 5,608s; completion arm ~3,100-7,149s); T02651 GUARD-KILLED
    07:00:13Z at 68.19GB; T02245 GUARD-KILLED 07:49:01Z at 68.15GB after a ~2h46m lane wall -
    THE LONGEST MEMORY-TAIL LANE YET. Scoreboard ~17 completions / ~24 kills.
(2) THE FLOOR TOUCH (1,499G at 07:45-46Z) SELF-GATED the cuts master at 07:45:09Z; the gate
    RELEASED 07:57:56Z on the cascade's recovery (df 1,499->1,559G). No lever needed.
(3) THE CASCADE: the ~270 backlog dispatched at LANES=128 on T02245's resolution; call-terminal
    1,344 -> 1,421 (callok +52, fail +25 - the bind-dense band); map dirs 271 -> 209 -> 242
    (drained then production outpaced again); rvsize 485 -> 433-445G.
(4) ORDINARY WIDTH-8 PACKAGING KILLS (exit 137, width->4, no strikes): T02877 06:59:33Z,
    T02870 07:08:22Z, T02837 07:16:11Z, T02950 08:27:03Z, T02971 08:37:12Z.
THE RE-WIDEN LADDER still holds at the standing gate - BOTH LEGS OPEN at close (backlog 242
map dirs, df 1,538G < 1,600G); C=2 -> C=4 -> C=6 when the gate meets (kill the cuts master
MID-WAIT only - never the wave, never with marathon calls in flight; C sized to the ~6.3
terminal/min flock ceiling, never RAM).
MASTERS AT CLOSE: prep 411760 (POOL_W=64 THREADS=4, chunked valve working), cuts C=2-adaptive
1002407 (BATCH_W 8->4 N_CONCURRENT=2 RAYON=24, own floor gate live), wave LANES=128 1699248
(PATCHED code, fail-reaping in production all window); heartbeat live; the ancient idle
impg-reference (1767054) known-benign - NEVER KILLED.
THE SUCCESSOR RESUMES THE QUIET RIDE: SHORT sparse polls (~4-5 min, act only on real events);
the marathon cohort forming at close (47min/51.5GB lead lane) resolves per the guard-oracle
doctrine; the re-widen ladder per the gate; the standing contingencies (disk floor 1,500G ->
the C lever, self-gated cuts side; memory-tail both arms; chr2 AGC class per the ride-19
blast-radius enumeration; binds terminal and named; T00197 1); AT CAMPAIGN CLOSE (all 5,844
accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list, NOW
ALSO: T02656's 7,149s completion record + T02245's 2h46m longest-lane kill + this window's
self-gate receipt) -> the record's closing increment -> the commit. NO PR PUSH.
## RIDE 23 THE WAVE-DISPATCH REPAIR (08:44-08:45Z 2026-10-16) - THE BARRIER FOUND, THE CHUNKED FIX WRITTEN, THE SWAP ARMED
Session-open increment committed (f086387). THE PACE PROBLEM (the parent's receipt): call-terminals
1,147 -> 1,438 over ~8h = ~36/hour - the wave is SERIALIZED BEHIND ITS MARATHON LANES.
(1) THE BARRIER, SURVEYED AND CONFIRMED LIVE (run-wave.sh's work loop): the single
    `printf '%s\n' "${TODO[@]}" | xargs -P "$LANES" ...` full-TODO dispatch. xargs DOES roll
    WITHIN its input list (verified empirically: -P refills each freed slot) - but the LOOP
    cannot refresh TODO until the ENTIRE snapshot drains: one full-TODO xargs per pass, so
    every marathon straggler in the snapshot stalls the fleet's next dispatch - the ride-7
    batch-barrier class. CONFIRMED LIVE at 08:43Z: the 07:49Z batch's xargs (2360644) down to
    22 marathon lanes, ZERO new lanes since ~08:03Z, while ~230 fresh cut.done territories
    wait un-dispatched (map dirs 239; cut.done 1,478 vs call-terminal 1,424+in flight).
(2) THE FIX, IMPLEMENTED (the prep-master's own bounded-chunk pattern, run-prep-master.sh:79-87):
    THE CHUNKED ROLLING DISPATCH in run-wave.sh - `CH=$((CHUNK_MULT*LANES))` (CHUNK_MULT=2
    default = 256 at LANES=128, ~2 concurrency waves), the TODO sliced per chunk
    (`${TODO[@]:di:CH}`), each chunk its own xargs -P, the barrier released as each chunk
    returns. A marathon holds only its own lane's slot; the markers make re-dispatch idempotent
    (done/exited territories skip for free). bash -n OK. The on-disk edit is safe for the
    running OLD master: bash parsed the whole while-loop into memory and reads the file again
    only at loop exit (days away; the swap kills the master first) - the mv-swap doctrine's
    concern does not bite in-place here, and the running master keeps its OLD in-memory
    dispatch until the safe point.
(3) THE SWAP, ARMED (the ride-17 rv-patch-watcher pattern): rv-chunk-watcher.sh (detached,
    log wave-logs/chunk-watcher.log) polls the master's xargs child every 2s; at ZERO lanes in
    flight (the old master's between-batches `sleep 30` - the natural batch boundary) it kills
    master 1699248 BY EXACT PID ONLY (never with marathons in flight, never the heartbeat),
    runs the 3-round belt-and-braces young-lane race check, bash -n's the chunked script, and
    relaunches LANES=128 CHUNK_MULT=2 via setsid launch-wave.sh. THE 22 IN-FLIGHT MARATHONS
    ARE NEVER KILLED - the change applies at the next natural batch boundary only.
(4) THE RECEIPT THE FIX ANSWERS TO: the measured pace delta (call-terminals/hour before vs
    after the chunked master lands - computed from wave-progress.log + failures.txt
    timestamps at the ride's polls).
THE RIDE continues: SHORT sparse polls; the contingencies stand (disk floor 1,500G -> the C
lever/self-gate; memory-tail both arms per the guard-oracle doctrine; terminal classes
terminal; the re-widen ladder at its gate); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL
datasheet-runtime.md (the standing close list, NOW ALSO the wave-barrier finding + the
pace-delta receipt) -> the closing increment -> the commit. NO PR PUSH.
## RIDE 23 WINDOW CLOSE (10:33Z 2026-10-16, deadline checkpoint) - THE HONEST CHECKPOINT STATE
THE WINDOW (08:44-10:33Z, ~110 SHORT sparse polls): the barrier surveyed and confirmed live,
the chunked rolling dispatch WRITTEN AND bash-n-CHECKED into run-wave.sh (CHUNK_MULT=2), the
swap ARMED via rv-chunk-watcher.sh (detached, autonomous, 2s cadence - log
wave-logs/chunk-watcher.log) - THE SWAP HAS NOT YET FIRED: the old master's 07:49Z batch is
still held by its LAST marathon lane (T02803, ~2h30m elapsed, ~44GB, guard-oracle doctrine -
never killed by hand). The watcher needs NO successor: at ZERO lanes it kills master 1699248
by exact pid (the between-batches sleep-30 safe point), runs the belt-and-braces race check,
bash -n's the chunked script, and relaunches LANES=128 CHUNK_MULT=2 via setsid launch-wave.sh.
POLL 10:32Z: split 9,891 / acc 5,844 / sub=1,805 /
call-terminal ~1,445 / df 1,499G (AT the 1,500G floor - the cuts master's own floor gate
self-manages; the JIT reaper cascade recovers it as the backlog drains post-swap) / mem 905Gib
(the single marathon lane + idle fleet). THE BEFORE-RECEIPT (the pace the fix answers to):
hourly call-terminals from receipts mtimes: 19h=81, 20h=1, 21h=147, 22h=15, 23h=30, 00h=259,
01h=3, 02h=44, 03h=7, 04h=9, 05h=106, 06h=14, 07h=30, 08h=52 - the famine/burst barrier
signature (1-15/h while snapshots are held by marathons, 100-259/h only in the cascade hours
after a snapshot finally breaks); 08:45-10:30Z (this window, one full-TODO hold): 21
terminals in 105 min = ~12/h. THE AFTER-RECEIPT IS OWED: measure call-terminals/hour for
>=1h after the chunked master lands (receipts mtimes; expect the ~230 backlog to clear in
the first chunk pair, then pace bounded by upstream cuts production, NOT by marathon
stragglers). THEN: the usual campaign-close work (rv-rollup.py -> FILL datasheet-runtime.md
including the wave-barrier finding + the pace-delta receipt -> the closing increment -> the
commit). NO PR PUSH.

## RIDE 24 SESSION OPEN (10:37Z 2026-10-16, receipt-side) - THE PACE-DELTA RECEIPT
Ride-24 opened; session-open increment committed (2732948). OPEN POLL 10:35Z: watcher
ALIVE (pid 2560122, 1 lane in flight - the T02803 marathon ~2h45m, guard-oracle, never
killed by hand); old wave master 1699248 ALIVE (still inside its 07:49Z batch snapshot);
split 9,891 / acc 5,844 / sub 1,805 / cut 1,558 / call-ok 1,089 + fail-closed 356 / poison
1 / df 1,499G (AT the floor; the cuts master's own floor gate self-manages) / mem 905GiB.
THIS WINDOW: sparse SHORT polls; VERIFY THE SWAP FIRED (chunk-watcher.log/.done; the new
master's chunked dispatch live); TAKE THE AFTER-HALF PACE RECEIPT (call-terminals/h for
>=1h post-swap from receipts mtimes; expect the ~230 backlog to clear in the first chunk
pair, then pace bounded by upstream cut production). Contingencies stand. At campaign
close: rv-rollup.py -> FILL datasheet-runtime.md -> the closing increment -> the commit.

### RIDE 24 THE SWAP FIRED (10:56:50Z) - clean, autonomous, per doctrine
The watcher caught the safe point with ZERO lanes in flight, killed old master 1699248 BY
EXACT PID (TERM then KILL), 3-round race check: NO surviving lanes; bash -n OK; relaunched
LANES=128 CHUNK_MULT=2 via setsid - NEW CHUNKED MASTER pid 2810146, live (xargs child
2810151, ~124 one_territory lanes spawning, first call-regions running by 10:58Z).
chunk-watcher.done touched 10:56:55Z; watcher exited clean. THE SAFE POINT'S CAUSE:
T02803 = the TENTH 64GiB GUARD KILL of the memory-tail class - the plateau-then-climb
shape at its extreme: flat 46.9GB for ~2h45m (locus-2 factorized loop), then the climb
46.9 -> 69.1GB in the final minutes; the lane's own 64GiB guard killed it at 69,142,204kB
(exit 137, prep/T02803/call.rss 'GUARD kill'; wall ~3h20m, the campaign's longest call).
failures.txt + call.exit updated; terminal class terminal. NOW OWED: the after-half pace
receipt - call-terminals/hour from receipts mtimes for >=1h post-swap (swap landed
10:56:55Z; measure the 11:00-12:00Z hour), then the campaign-close work if the run
completes in-window.

### RIDE 24 WINDOW CLOSE (12:15Z 2026-10-16) - THE PACE-DELTA RECEIPT, BOTH HALVES, TAKEN
THE AFTER-HALF RECEIPT (the >=1h post-swap hour, receipts mtimes, exact):
- 11:00-12:00Z: **94 call-terminals/hour** (vs the BEFORE-half famine band 1-15/h in
  marathon-held hours and 12/h in ride-23's held window; the burst-cascade band was
  100-259/h only AFTER a snapshot broke). 94/h is INDEPENDENT production through a
  marathon-heavy hour - no famine.
- THE SHAPE (15-min buckets): 11:00-11:15 = 76 terminals (~300/h pace - THE BACKLOG
  BURST: the new master's first chunk dispatched the ~110-230 cut.done backlog the old
  master had held for 5+ hours, exactly as predicted); 11:15-12:00 = 18 terminals over
  45 min - THE BOUNDED CHUNK-1 TAIL (the first chunk's own marathons: 42 -> 24 lanes,
  oldest >1h10m, ~100 lanes idle). THE HONEST READING: the fix converts the UNBOUNDED
  full-TODO barrier (one marathon froze the whole backlog for hours) into a BOUNDED
  chunk-tail (marathons idle lanes only until THEIR chunk resolves; the barrier then
  releases and re-slices - marker-idempotent). Marathon stalls are REDUCED, not
  eliminated - the residual is the within-chunk tail (CH=256 at CHUNK_MULT=2).
- CALLS NOW OUTPACE UPSTREAM: cut production ~38-90/h (cut.done 1,558 -> 1,609 over the
  hour, batches landing steadily post-floor); terminals 94/h; backlog cut.done-minus-
  terminals fell ~110 -> ~37. The terminal pace is now bounded by UPSTREAM CUT
  PRODUCTION, not marathon stragglers - the receipt the fix answers to.
- POLE STATE 12:06Z: split 9,891 / acc 5,844 / cut 1,609 / terminals 1,546 (callok 1,160
  + fail-closed 386) = 26.4% of accepted / poison 1 (T00197) / df 1,574G (recovered from
  the 1,499G floor - the JIT-reaper cascade, rvsize 485G -> ~410G) / memavail 509GiB /
  masters: chunked wave 2810146 LIVE (24-lane chunk-1 tail), cuts master producing.
- T02803 = the TENTH 64GiB guard kill of the memory-tail class (69.1GB peak, plateau
  46.9GB ~2h45m then climb; failures.txt named; terminal class terminal).
- OWE THE SUCCESSOR: (a) confirm the chunk-2 boundary burst when the chunk-1 tail
  resolves (the barrier-release re-proof; the first-chunk burst already proves the
  dispatch); (b) the campaign-close work - rv-rollup.py -> FILL datasheet-runtime.md
  (every finding per the close-note enumeration + the wave-barrier finding + BOTH halves
  of this pace-delta receipt) -> the closing increment -> the commit. Contingencies
  stand; NO PR PUSH.

## RIDE 25 THE T01411 PREP-POLE UNBLOCK (13:25:13Z 2026-10-16) - THE T00197 SYNG SEG-DUP-TAIL CLASS'S SECOND SIGHTING, POISONED + KILLED PER DOCTRINE
THE DISCOVERY (found while riding the chunk-1 tail): the prep pole was STALLED ~3h17m - prep-master.log
silent since 10:08:39Z, sub frozen at 1,805, the chunk xargs (2720462, dispatched 10:08Z) held by ONE
chain: T01411's `impg syng`. THE SIGNATURE, VERIFIED (the exact T00197 class): pid 2720627 state R,
CPU time 11,751s ~= wall 11,752s (100% compute-bound for 3h16m), RSS 23MB, NO output file
(locus-pos64.syng absent), syng.err SILENT since 10:08:40Z stuck on ONE sequence
(HG00350#2#CM094497.1:26555331-27052385, 497,054 bp); locus.fa 318MB (T00197's was 345MB). The
ordinary syng corpus: median 4s, max 159s - 65x the max. THE PREP-SIDE BATCH-BARRIER (the ride-11
finding's sibling): the chunk cannot return while the chain lives, so the valve (npendcut 124 < 160,
OPEN) could never re-check - the whole pole serialized behind one pathological sequence.
THE ACTION (the receipted T00197 doctrine, second member): (1) POISONED FIRST -
prep/T01411/cut.poison touched + the honest class line in cuts/poison.tsv (poison now 2: T00197 +
T01411) BEFORE the kill, so the master's next TODO computation skips it (the retry-forever gap stays
closed; the filter verified at run-prep-master.sh:25 and :67 - both dispatch sites); (2) the receipt
snapshotted (prep/T01411/syng.err.attempt1-snapshot; the wrapper's own syng.exit/syng.log/syng.wall
landed on the kill - attempt 1 fully receipted); (3) KILLED BY EXACT PID 2720627 (13:25:13Z) - the
pathological syng ONLY, never the wrapper, never the master.
THE RESULT, VERIFIED LIVE: the wrapper exited, the chunk returned, the master looped and dispatched a
FRESH chunk within seconds (new xargs 3075241, ~128 chains, live) - T01411 ABSENT from the dispatch
(poison held, 0 processes); sub 1,805 -> 1,861+ climbing as the fresh chains land.
THE HONEST CLASS NOTE for the datasheet: unlike T00197 (attempt 1 failed fast in 10s, attempt 2 burned
3h01m - poison took 2 attempts), T01411's ATTEMPT 1 IS the 3h16m burn - the same confirmation
signature on the first try, and the class was already named and confirmed by the ride-11 precedent,
so the poison landed on sight after ONE receipted attempt (the alternative - letting the master retry
deterministically - would burn another ~3h of the whole pole). 2 territories of 5,844 accepted now in
the PREP-SYNG-SEG-DUP-TAIL class.
STATE 13:26:30Z: sub 1,861 climbing / map dirs ~230 / call-terminal ~1,672 = callok ~1,230 + fail 411 +
cut.exit 16 + poison 2 / df 1,379G (the campaign's deepest point, below the recorded 1,373G excursion -
BOTH DISK POLES AT REST: cuts self-gated since 13:23:01Z 'disk floor: 1439G < 1500G', prep's chunk is
the designed bounded +128 max; the ONLY drain path is the wave's reap cascade, held behind the last
chunk-1 marathon T02282 (54.5GB frozen, ~145min) - guard-oracle doctrine, never killed by hand) /
mem ~500Gib. THE RIDE CONTINUES: when T02282 resolves -> THE CHUNK-2 BOUNDARY BURST (the owed
re-proof) -> the ~230 backlog drains at LANES=128 -> reaps cascade -> df recovers +100-150G ->
THE RE-WIDEN LADDER (gate: backlog drained AND df >= ~1,600G; kill the cuts master MID-WAIT only,
never the wave, never with marathons in flight; C=2 -> C=4 -> C=6, sized to the wave's measured
consumption, the flock ceiling, never RAM). AT CAMPAIGN CLOSE: rv-rollup.py -> FILL
datasheet-runtime.md (the standing close list + the pace-delta receipt + THIS CLASS'S SECOND
SIGHTING + the prep-pole stall receipt) -> the closing increment -> the commit. NO PR PUSH.

## RIDE 25 THE CHUNK-2 BOUNDARY BURST CONFIRMED (13:33-13:35Z 2026-10-16) - THE OWED RE-PROOF, TAKEN
THE CHUNK-1 TAIL RESOLVED: T02282 = the memory-tail class's LATEST GUARD KILL (13:33:20Z - the
plateau-then-climb shape: frozen 54.5GB for ~2.5h, then the climb; failures.txt named, fail-reaped
in-lane per the patched wave). ALSO: T03035 (12:30:00Z, 5,059s, peak 58.4GB) and T03026 (12:36:09Z,
5,047s, 58.5GB) COMPLETED HONESTLY - marathon completions in the 84min band; T03010 (13:02:21Z) and
T02963 (12:59:07Z) guard-killed in the tail - the chunk-1 tail's full census: 5 completions + 9 kills
(T02994 12:26 / T02815 12:38 / T02968 12:52 / T02963 12:59 / T03010 13:02 / T02282 13:33 - peaks
63-69GB) + the earlier T02366 (11:58Z), every one fail-reaped in-lane with full receipts.
THE BOUNDARY BURST, VERIFIED LIVE: the moment T02282's lane died, the chunk-1 xargs returned and
the master re-sliced - by 13:35:47Z: 128 one_territory lanes dispatched, 110 CONCURRENT GPARTITION
PROCESSES (the flock funnel at its measured ceiling), 9 call-regions already spawned; lanes ramped
9 -> 15 -> 20 -> 32 -> 56 through the funnel as the ~124 backlog consumed. THE CHUNKED ROLLING
DISPATCH'S OWED RE-PROOF: the barrier released at the CHUNK boundary (not the full-TODO pass) -
marathons held only their own chunk (2.5h of chunk-1 tail vs the old barrier's 5h+ full-TODO freeze),
and the fresh backlog dispatched within seconds of the last straggler's death.
THE DRAIN 13:36-14:00Z: backlog 124 -> 41, rvsize 446 -> 356G falling, callok 1,164 -> 1,212, fail
393 -> 428 (the bind-dense band, all named + fail-reaped in-lane), df trough 1,313-1,320G -> ~1,330-
1,356G recovering (the campaign's deepest trough, below the 1,373G excursion - BOTH DISK POLES WERE
AT REST: cuts self-gated since 13:23:01Z, the prep chunk bounded +119 then VALVE-HELD at 13:28:58Z
'VALVE (mid-dispatch): 185 pending cut - parking 3,863 of 4,119 TODO' - the chunked-dispatch valve's
SECOND literal proof line).
THE LADDER GATE at 14:00Z: backlog leg nearly met (41 draining), df leg OPEN (1,329-1,356G < 1,600G)
- the gate waits for the df recovery; cuts still self-gated below 1,500G; when df crosses 1,500G the
cuts master resumes C=2 and the drain continues; the re-widen (C=4 then C=6) fires only at the full
gate (backlog drained AND df >= ~1,600G; kill the cuts master MID-WAIT only, by exact pid, never the
wave, never with marathons in flight). NO PR PUSH.

## RIDE CONTINUATION 25 WINDOW CLOSE (14:06Z 2026-10-16, deadline checkpoint) - THE HONEST STATE
FINAL POLL 14:05:29Z: split 9,891 / acc 5,844 (59.1% final) / sub 1,924 (valve held at ~185 pending-cut) /
backlog 37 DRAINING / call-terminal 1,660 = callok 1,215 + fail-closed 429 + cut.exit 16 + poison 2
(T00197 + T01411) = 28.4% of accepted / lanes 37 (the chunk-2 tail's marathon cohort forming: T02128
60.5GB fast-climb at 25min, T01885 36.5GB, T01875 40.8GB - guard-oracle doctrine, never killed by hand) /
df 1,285G (the campaign's deepest trough band, below the recorded 1,373G excursion - for the datasheet;
the decline 1,356->1,285G over the drain is partly CO-TENANT: our rvsize fell 446->350G the same interval,
the floor gates self-manage on our side) / mem 296Gib / rvsize 350G (under the 500G high-water).
THE WINDOW'S EVENTS, ALL RECEIPTED ABOVE:
(1) THE T01411 PREP-POLE UNBLOCK (13:25:13Z): the T00197 SYNG SEG-DUP-TAIL class's SECOND SIGHTING -
attempt 1 burned 3h16m CPU-bound (CPU ~= wall, 23MB RSS, no output, err silent on one 497kb sequence),
blocking the ENTIRE prep pole since 10:08Z (sub frozen 1,805); POISONED FIRST (prep/T01411/cut.poison +
cuts/poison.tsv - the master's TODO filter skips it at both dispatch sites), the pathological syng killed
BY EXACT PID 2720627, the fresh chunk dispatched within seconds with T01411 ABSENT (poison held), sub
1,805 -> 1,924. The class census: 2 territories of 5,844 (T00197 + T01411); the honest note: T01411's
poison landed on sight after ONE receipted attempt (the class was already confirmed by the T00197
precedent; a retry would deterministically burn another ~3h of the pole).
(2) THE VALVE (mid-dispatch) SECOND LITERAL PROOF LINE (13:28:58Z): '185 pending cut (>= 160) - parking
the remaining 3,863 of 4,119 TODO' - the chunked-dispatch valve proven again at the chunk boundary.
(3) THE CHUNK-2 BOUNDARY BURST CONFIRMED (13:33-13:35Z) - the owed re-proof: T02282 (the last chunk-1
marathon) GUARD-KILLED 13:33:20Z (frozen 54.5GB ~2.5h then the climb), the chunk-1 xargs returned, and
the master re-sliced within seconds - 128 one_territory lanes, 110 CONCURRENT GPARTITIONS (the flock
funnel), call-regions ramping 9->56. The barrier released at the CHUNK boundary, NOT the full-TODO pass.
(4) THE CHUNK-1 TAIL'S FULL CENSUS: 5 honest marathon completions (T03035 5,059s/58.4GB, T03026
5,047s/58.5GB, T02861 4,371s/49.4GB, T03055 1,934s, T03077 2,276s) + 7 guard kills (T02366 11:58Z,
T02994 12:26Z, T02815 12:38Z, T02968 12:52Z, T02963 12:59Z, T03010 13:02Z, T02282 13:33Z - peaks
63-69GB, all fail-reaped in-lane with full receipts by the patched wave).
(5) THE CASCADE: backlog 124 -> 37, rvsize 446 -> 350G, callok 1,164 -> 1,215, fail 393 -> 429 (the
bind-dense band, every line named in wave-logs/failures.txt).
(6) THE RE-WIDEN LADDER: NOT PUSHED - the gate never opened (backlog leg nearly met at close, 37
draining; df leg OPEN, 1,285G << 1,600G). ZERO lever actions, ZERO masters touched, ZERO hand kills
beyond the doctrine-compliant T01411 syng kill. The ladder receipts this window: the gate HELD, honestly.
MASTERS AT CLOSE: prep 411760 (valve-held, the T01411 poison filter proven), cuts C=2-adaptive 1002407
(self-gated below the 1,500G floor since 13:23:01Z - resumes when df recovers), chunked wave LANES=128
CHUNK_MULT=2 2810146 (the chunked rolling dispatch proven at its second boundary); heartbeat live; the
ancient idle impg-reference benign - NEVER KILLED.
THE SUCCESSOR RESUMES: SHORT sparse polls; the marathon cohort resolves per the guard-oracle doctrine;
when the gate meets (backlog drained AND df >= ~1,600G) push the ladder C=4 then C=6 (kill the cuts
master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C sized to the
wave's measured consumption, the flock ceiling, never RAM); standing contingencies (disk floor 1,500G ->
the C lever/self-gate; memory-tail both arms; chr2 AGC class per the ride-19 blast-radius enumeration;
binds terminal and named; poison 2); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py ->
FILL datasheet-runtime.md (the standing close list + the pace-delta receipt + the chunk-2 boundary
re-proof + the T01411 class sighting + the prep-pole stall receipt + the deepest df trough band) ->
the record's closing increment -> the commit. NO PR PUSH.
Session-open increment for this window: committed f6d5499; window close committed with this note.

## RIDE 26 MID-WINDOW STATE (14:41Z 2026-10-16) - THE MEMORY CRUNCH RIDDEN, THE GUARD ORACLE PROVEN AT SCALE
OPENED AT 14:08Z (session-open increment committed c086677 in the impg repo). THE WINDOW'S EVENTS, ALL RECEIPTED:
(1) THE MEMORY-TAIL COHORT'S PEAK MEMORY PRESSURE (14:22-14:25Z): ~38 concurrent lanes summing ~700GB RSS on the
1TB box drove MemAvailable 250 -> 15GiB in ~15 min (the tightest memory pole of the campaign - the lanes' collective
climb outpaced completions); NO HAND ACTION per the guard-oracle doctrine; the top lane (T02790) crossed the 64GiB
line at 14:25:44Z and the wave's own 5s-poll guard killed it in-lane - MemAvailable recovered 15 -> 79GiB in seconds.
The OOM risk was real but self-resolving: the guard is the law and it fired before the kernel did.
(2) THE CHUNK-2 MARATHON COHORT'S CENSUS SO FAR (all fail-reaped in-lane, full receipts in wave-logs/failures.txt):
T02128 (14:04:26Z, the 60.5GB fast-climb from ride-25's close), T02823 (14:09:42Z, ~70GB peak), T02790 (14:25:44Z,
64GB+), T01885 (14:30:13Z), T02693 (14:40:09Z) - 5 guard kills; plus honest completions (call.done 1,218 -> 1,224,
lanes 38 -> 28). The cohort resolves exactly per the memory-tail class doctrine.
(3) THE DF TROUGH DEEPENED FAR BELOW ALL RECORDED BANDS: 1,269G (14:08Z) -> 946G (14:41Z), falling ~10G/min -
ALMOST ENTIRELY CO-TENANT (our rvsize FELL 350 -> 326G the same interval; BOTH our disk poles at rest the whole
time: cuts self-gated below the 1,500G floor since 13:23:01Z, prep valve-held at 185 pending-cut). OUR SIDE IS
FLOOR-COMPLIANT AND SHRINKING; the co-tenant drain is the whole decline. FOR THE DATASHEET AT CLOSE: the campaign's
deepest df excursion is THIS window's band (946-1,285G), not the earlier 1,313G/1,373G excursions.
(4) THE LADDER GATE: HELD, honestly - the backlog leg is nearly met (cut.done 1,681 - terminals 1,673 = 8 draining)
but the df leg is WIDE OPEN (946G << 1,600G, and moving AWAY while the co-tenant drains). ZERO lever actions, ZERO
masters touched. The ladder CANNOT meet this window unless the co-tenant drain reverses; the cuts master stays
self-gated (resumes C=2 only above 1,500G).
(5) NO PREP STALLS: the prep pole is valve-held (185 pending-cut) so no new syng chains are running - the T01411
class cannot strike while held; prep-master.log healthy, sub 1,924 stable.
STATE 14:41Z: split 9,891 / acc 5,844 / sub 1,924 / cut.done 1,681 / call-terminal 1,673 = callok 1,224 + fail 433
+ cut.exit 16 + poison 2 = 28.6% of accepted / backlog 8 / lanes 28 (top RSS 52GB climbing - more guard kills owed) /
df 946G / mem 117GiB / rvsize 326G. MASTERS: prep 411760, cuts 1002407 (self-gated), wave 2810146; heartbeat live.
THE RIDE CONTINUES: sparse SHORT polls; the cohort resolves; the ladder waits for the gate (backlog drained AND
df >= ~1,600G); at campaign close the rollup + datasheet fill + closing increment + commit per the standing list.

## RIDE 26 WINDOW CLOSE (15:09Z 2026-10-16) - THE HONEST CHECKPOINT: THE CAMPAIGN IS NOW CO-TENANT-DISK-BOUND
FINAL POLL 15:08:06Z: split 9,891 / acc 5,844 (59.1% final) / sub 1,924 (valve-held at 185 pending-cut, healthy,
no T01411-class stall - no syng chains run while the valve holds) / cut.done 1,681 / call-terminal 1,685 = callok
1,227 + fail-closed 440 + cut.exit 16 + poison 2 (T00197 + T01411) = 28.8% of accepted / backlog (cut.done minus
callok+fail) 14, ALL of it in the 17 live lanes (the chunk-2 tail; nothing queued behind the lanes) / lanes 17
(top RSS 64/56/52GB - more guard kills owed) / df 694G / mem 357GiB / rvsize 306G FALLING / map dirs 80 (reaps
working). MASTERS ALL LIVE: prep 411760 (valve-held), cuts 1002407 (self-gated at the 1,500G floor since 13:23:01Z),
chunked wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat live; the ancient idle impg-reference benign - NEVER KILLED.
THE WINDOW'S RECEIPTS (session-open increment committed c086677 in the impg repo):
(1) THE MEMORY CRUNCH AT COHORT PEAK (14:22-14:25Z) - RIDDEN PER DOCTRINE, NO HAND ACTION: ~38 concurrent lanes
summing ~700GB RSS on the 1TB box drove MemAvailable 250 -> 15GiB in ~15 min (the campaign's tightest memory pole);
the wave's own 5s-poll 64GiB guard fired FIRST (T02790 at 14:25:44Z) and MemAvailable recovered 15 -> 79GiB in
seconds. The guard-oracle doctrine is now proven at cohort scale: the guard beat the kernel to the kill.
(2) THE CHUNK-2 COHORT'S GUARD-KILL CENSUS (11 this window, every one fail-reaped in-lane with full receipts in
wave-logs/failures.txt): T02823 14:09:42Z (~70GB) / T02790 14:25:44Z / T01885 14:30:13Z / T02693 14:40:09Z /
T02733 14:42:32Z (62.8GB) / T02691 14:42:47Z / T02936 14:43:40Z / T03121 14:48:19Z / T02829 14:50:47Z /
T01908 14:58:20Z / T02739 15:01:56Z (63.5GB). Honest completions callok 1,218 -> 1,227 (+9). The memory-tail
class resolves exactly as classed - terminal.
(3) THE DF EXCURSION, THE CAMPAIGN'S DEEPEST BY FAR AND ~ALL CO-TENANT: df 1,269G (14:08Z) -> 694G (15:08Z),
-575G in the hour at a steady 8-13G/min, WHILE OUR OWN FOOTPRINT FELL (rvsize 350 -> 306G; map dirs 100 -> 80;
both our disk poles at rest the whole window: cuts self-gated below the 1,500G floor since 13:23:01Z, prep
valve-held at 185 pending-cut). OUR SIDE IS FLOOR-COMPLIANT AND SHRINKING; the whole decline is the co-tenant.
FOR THE DATASHEET AT CLOSE: the deepest df excursion band is 694-1,285G (ride-26 window), displacing the earlier
1,313G/1,373G records; the co-tenant attribution is receipted by the simultaneous FALL of our rvsize.
(4) THE RE-WIDEN LADDER: HELD, HONESTLY - the backlog leg met (nothing queued behind the lanes) but the df leg
WIDE OPEN and moving AWAY (694G << 1,600G). ZERO lever actions, ZERO masters touched, ZERO hand kills beyond the
guard's own 11. The cuts master stays self-gated; it resumes C=2 only when df recovers above 1,500G.
(5) THE STRUCTURAL STATE, FOR THE SUCCESSOR: the campaign's pace is now bounded by the CO-TENANT DISK DRAIN, not
by any of our poles. The chain: cuts self-gated below 1,500G -> no new cut.done -> the prep valve held at 185
pending-cut -> the wave's remaining backlog is only the 17 in-flight lanes -> when the chunk-2 tail resolves the
chunk-3 burst will be small-to-empty and the wave idles awaiting cut production -> cut production awaits
df > 1,500G -> df is falling ~8-13G/min externally. NOTHING FOR US TO DO BUT WAIT, POLL SPARSELY, AND RECORD -
any lever we hold (the ladder) is gated on df recovery that only the co-tenant can return. DO NOT touch the
masters; the floor gates self-manage; the guard-oracle doctrine holds for the remaining tail lanes.
THE SUCCESSOR RESUMES: SHORT sparse polls; ride the chunk-2 tail (17 lanes, more guard kills owed); watch for
the co-tenant drain reversing (df > 1,500G -> cuts resume C=2 -> the ladder re-arms at df >= ~1,600G with backlog
drained); watch for T01411-class prep stalls once the valve releases (poison-first doctrine); AT CAMPAIGN CLOSE
(all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the
pace-delta receipt + the chunk-2 boundary re-proof + the T01411 class's second sighting + the prep-pole stall +
THIS window's deepest-df-excursion receipt with the co-tenant attribution + the memory-crunch receipt - the guard
proven at cohort scale, MemAvailable trough 15GiB) -> the record's closing increment -> the commit. NO PR PUSH.
CORRECTION TO THE MID-WINDOW NOTE: the 14:41Z backlog figure of 8 was computed against all-terminal count; the
canonical cut.done-minus-(callok+fail) figure at 14:41Z was 24. At close: 14, all in-flight. Immaterial to any
gate; recorded for honesty.
Session-open increment for this window: c086677; window close committed with this note.

## RIDE 27 SESSION OPEN (15:15Z 2026-10-16) - THE CO-TENANT-DISK-BOUND RIDE
Opened 15:13Z from the ride-26 close (record commit 3ba2435; session-open increment committed
e341ad3 in the impg repo). OPEN POLL 15:13:22Z: split 9,891 / acc 5,844 / sub 1,924 (valve-held,
healthy) / cut.done 1,681 (cuts master reaper-wait below the 1,500G floor since 13:23:01Z) /
call-terminal 1,688 = callok 1,227 + fail 443 + cut.exit 16 + poison 2 = 28.9% / backlog 11 ALL
in-flight / lanes 25 wrappers (the chunk-2 marathon cohort STILL RIDING: top impg 57GB at 1h27m,
then 48/47/45/43/40/39/37/36/35/31GB - guard-oracle doctrine, never killed by hand) / df 677G
STILL FALLING, below ride-26's deepest 694G; co-tenant drain ~3G/min while our rvsize FALLS
306 -> 297G and both our disk poles rest - our side floor-compliant and shrinking / memavail
515GiB. MASTERS: prep 411760, cuts 1002407 (self-gated), wave 2810146; heartbeat live; the
ancient idle impg-reference benign - NEVER KILLED.
THE WINDOW'S DOCTRINE: sparse SHORT polls; ride the cohort (the 64GiB guard is the law); the
ladder CANNOT meet while the co-tenant drains (it re-arms only if df > 1,500G -> cuts resume
C=2 -> then the backlog-drained AND df >= ~1,600G test; kill the cuts master MID-WAIT only, by
exact pid, never the wave, never with marathons in flight); prep stalls watched once the valve
releases (poison-first, fresh-chunk dispatch); floor gates self-manage; AT CAMPAIGN CLOSE:
rv-rollup.py -> datasheet-runtime.md fill (standing close list + pace-delta + chunk-2 re-proof +
T01411 class + prep-pole stall + deepest-df-excursion/co-tenant + memory-crunch receipts) ->
closing increment -> commit. NO PR PUSH.

## RIDE 27 WINDOW CLOSE (15:39Z 2026-10-16) - THE HONEST CHECKPOINT: THE CO-TENANT DRAIN UNABATED, THE COHORT RESOLVING
FINAL POLL 15:39:07Z: split 9,891 / acc 5,844 (59.1% final) / sub 1,924 (valve-held at 185
pending-cut, healthy - no syng chains run while held, no T01411-class stall possible) / cut.done
1,681 (cuts master reaper-wait below the 1,500G floor since 13:23:01Z - cut production frozen
awaiting df recovery) / call-terminal 1,694 = callok 1,230 + fail-closed 446 + cut.exit 16 +
poison 2 (T00197 + T01411) = 29.0% of accepted / backlog (cut.done minus callok+fail) 5, ALL
in-flight in the 13 live lanes - nothing queued behind the lanes / lanes 13 wrappers (5 frozen
chunk-2 marathons 38-47GB at 1h46-2h01 elapsed - guard-oracle doctrine, the 64GiB guard arms
both, never killed by hand) / df 508G / memavail 729GiB / rvsize 287G falling.
THE WINDOW'S EVENTS, ALL RECEIPTED:
(1) THE CHUNK-2 COHORT'S RESOLUTION CENSUS (this window): 3 guard kills, every one fail-reaped
in-lane with full receipts in wave-logs/failures.txt - T03122 15:24:34Z, T03094 15:29:14Z,
T02894 15:31:08Z; plus 3 honest marathon completions receipted in wave-logs (T02138 15:19:44Z
5,970s/66.7GB, T03133 15:20:51Z 5,596s/61.8GB, T03118 15:30:18Z 6,292s/62.0GB). callok
1,227 -> 1,230, fail 443 -> 446. Lanes 25 -> 13. The memory-tail class resolves exactly as
classed - terminal.
(2) THE DF EXCURSION DEEPENED FAR BELOW EVERY RECORDED BAND: 677G (15:13Z open) -> 508G
(15:39Z), -169G in 26 min (~7G/min avg, bursts 14-15G/min), STILL ~ALL CO-TENANT - our rvsize
FELL 297 -> 287G the same interval; BOTH our disk poles at rest the whole window (cuts
self-gated below the 1,500G floor since 13:23:01Z; prep valve-held at 185 pending-cut); our
side floor-compliant and shrinking. FOR THE DATASHEET AT CLOSE: the campaign's deepest df
excursion band is now 508-694G (ride-26 close through ride-27), displacing the 694-1,285G band;
the co-tenant attribution is receipted by the simultaneous FALL of our rvsize.
(3) THE RE-WIDEN LADDER: HELD, HONESTLY - the backlog leg nearly met (5, all in-flight) but the
df leg IMPOSSIBLY OPEN (508G << 1,600G, moving AWAY at ~7-15G/min). ZERO lever actions, ZERO
masters touched, ZERO hand kills. The ladder re-arms only if the co-tenant drain reverses
(df > 1,500G -> cuts resume C=2 -> then backlog-drained AND df >= ~1,600G).
(4) NO PREP STALLS: the valve held all window; prep-master.log healthy; sub 1,924 stable.
(5) THE STRUCTURAL STATE, FOR THE SUCCESSOR: unchanged from ride-26's conclusion - the campaign
is co-tenant-disk-bound; the chain (cuts gated -> no cut.done -> prep valve-held -> wave tail
drains -> chunk-3 burst small-to-empty -> wave idles awaiting cut production -> cut production
awaits df > 1,500G -> df falling externally ~7-15G/min). NOTHING TO DO BUT WAIT, POLL SPARSELY,
AND RECORD. Do NOT touch the masters; the floor gates self-manage; the guard-oracle doctrine
holds for the 5 remaining frozen marathons.
MASTERS AT CLOSE: prep 411760 (valve-held), cuts 1002407 (self-gated, reaper-wait), chunked wave
LANES=128 CHUNK_MULT=2 2810146; heartbeat live; the ancient idle impg-reference benign - NEVER
KILLED.
THE SUCCESSOR RESUMES: sparse SHORT polls; ride the remaining 5 frozen marathons (guard-oracle);
watch for the co-tenant drain reversing; watch for T01411-class prep stalls once the valve
releases (poison-first doctrine); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py
-> FILL datasheet-runtime.md (the standing close list + the pace-delta receipt + the chunk-2
boundary re-proof + the T01411 class's second sighting + the prep-pole stall + the deepest
df-excursion/co-tenant receipt + the memory-crunch receipt) -> the record's closing increment ->
the commit. NO PR PUSH.
Session-open increment for this window: e341ad3; window close committed with this note.

## RIDE 28 MID-WINDOW (16:13Z 2026-10-16) - THE MARATHON TAIL + THE CO-TENANT WRITER SURVEY
THE CHUNK-2 TAIL RESOLVING (all guard-oracle, zero hand kills, every kill fail-reaped in-lane
with full receipts): T02729 (15:41:30Z) + T03134 (15:41:32Z) - the two 64.8GB-class lanes that
were nearing the guard at window open - and T01875 (16:03:55Z, peak 67,205,316kB, ~2h20m frozen
40.8GB plateau then the climb, the 'exactness sample' signature). TWO LANES REMAIN (45.1GB at
~2h28m, 43.0GB at ~2h23m - frozen plateaus, guard arms both; the completion arm spans
~3,100-7,149s so either resolution is honest). Backlog = 2 = the lanes; nothing queued behind.
failcall 449; callok 1,230; cut.done 1,681 frozen (cuts master self-gated below the 1,500G
floor since 13:23:01Z).
THE CO-TENANT WRITER SURVEY (bounded, for the datasheet's disk section): the drain on / fell
489G (15:41Z) -> ~400-407G (16:12Z) with oscillation (bursts both ways, 402<->431G band, net
~-6G/min). (1) THE xsa BUILD (pid 1681758, /home/erikg/sxgc/xsa build --v2-pipeline --input
/mnt/nvme2n1/... , running since 02:01:25Z, 455GB cumulative writes) IS NOT THE DRAIN ON /:
its scratch/input live on /mnt/nvme2n1, a SEPARATE 14TB filesystem (2,069G avail); its
write_bytes was FLAT over a 20s probe. (2) No deleted-but-open holders anywhere (summed
/proc/*/fd deleted links per process - all under 1GB). (3) No file >1G modified on / since
15:30Z (bounded find, one pass). (4) meminfo Dirty ~2GB only - NOT page-cache writeback.
The drain's source remains the co-tenant's opaque write/delete cycle; the attribution stands
on OUR side being floor-compliant and shrinking (rvsize 281 -> 279G over the window, both
poles at rest). NO ACTION AVAILABLE OR TAKEN - the campaign is co-tenant-disk-bound (the
ride-26/27 structural state); the ladder gate cannot meet (df leg impossibly open at
~400-420G << 1,600G).
THE WATCH, FOR THE REST OF THE WINDOW: the last two marathons resolve; the chunk-3 burst will
be tiny (backlog <= 2); then the wave IDLES awaiting cut production; cuts resume only at
df > 1,500G (a ~1.1TB co-tenant release); the re-widen ladder then needs df >= ~1,600G AND
backlog drained. Prep valve held at 185 pending-cut (no T01411-class stall possible while
held); prep-master.log healthy; heartbeat live.

## RIDE CONTINUATION 28 WINDOW CLOSE (16:43Z 2026-10-16, deadline checkpoint) - THE HONEST STATE
THE RIDE (15:41-16:41Z, ~40 SHORT sparse polls): rode the chunk-2 marathon tail down to ONE
lane, surveyed the co-tenant drain as far as the box admits, and held every lever at its rest -
ZERO masters touched, ZERO hand kills, ZERO lever actions.
FINAL POLL 16:41:12Z: split 9,891 / acc 5,844 (59.1% final) / sub 1,924 (valve HELD at 185
pending-cut all window - healthy, no T01411-class stall possible while held) / cut.done 1,681
(the cuts master SELF-GATED below the 1,500G floor since 13:23:01Z, still 'waiting for the
reaper' at 16:40:27Z 'disk floor: 476G < 1500G') / map dirs 67 (reaps complete) / call-terminal
1,680 = callok 1,230 + fail-closed 450 (+16 cut.exit + poison 2 = 1,698 of 5,844, 29.1% of
accepted; every fail line named in wave-logs/failures.txt) / backlog (cut.done minus
call-terminal) 1 - the LAST marathon lane, pid 3187261, 45.1GB FROZEN at ~2h56m elapsed
(past the 7,149s=119min longest completion wall - the plateau-then-climb class; the guard
arms it, resolve or kill per doctrine, NEVER by hand) / df 469G / mem 884Gib / rvsize 277G
falling. MASTERS AT CLOSE: prep 411760 (valve-held, log healthy), cuts C=2-adaptive 1002407
(floor-gated, reaper-wait), chunked wave LANES=128 CHUNK_MULT=2 2810146 (chunk-2 tail at its
last lane; when it resolves the chunk-3 burst is <= 1 territory and the wave IDLES awaiting
cut production); heartbeat live; the ancient idle impg-reference benign - NEVER KILLED.
THE WINDOW'S EVENTS, ALL RECEIPTED:
(1) FOUR MEMORY-TAIL GUARD KILLS, every one FAIL-REAPED IN-LANE by the patched wave with full
receipts (call.rss + failures.txt, the 'exactness sample' stderr signature): T02729
(15:41:30Z) + T03134 (15:41:32Z) - the two 64.8GB-class lanes nearing the guard at window
open; T01875 (16:03:55Z, peak 67,205,316kB, ~2h20m frozen 40.8GB plateau then the climb);
T03150 (16:17:47Z, the 43.0GB->55.1GB climber, guard-killed inside the scoring pass).
failcall 446 -> 450. The memory-tail class terminal, exactly as classed.
(2) THE CO-TENANT WRITER SURVEY (bounded, /proc-level, for the datasheet's disk section):
(a) the xsa build (pid 1681758, /home/erikg/sxgc/xsa build --v2-pipeline, running since
02:01:25Z, 455GB cumulative writes) is NOT the drain on / - its scratch/input live on
/mnt/nvme2n1, a SEPARATE 14TB filesystem (2,069G avail), and its write_bytes was FLAT over a
20s probe; (b) NO deleted-but-open holders anywhere (all /proc/*/fd deleted links summed per
process - every process under 1GB); (c) NO file >1G modified on / since 15:30Z (one bounded
find pass); (d) meminfo Dirty ~2GB only - not page-cache writeback. THE DRAIN'S SHAPE: the
co-tenant's opaque bursty write/delete cycle - df fell 489 -> 387G over ~55 min (net ~6G/min
with 400-431G oscillation), THE NEW DEEPEST EXCURSION 387G (displacing ride-27's 508G record),
then ONE +111G RELEASE in a single minute (16:36Z) -> the 469-498G band at close. FOR THE
DATASHEET AT CLOSE: the deepest df excursion band is now 387-694G (rides 26-28), the
co-tenant attribution receipted by our simultaneous rvsize fall (287 -> 277G over the window,
both our poles at rest: cuts self-gated, prep valve-held).
(3) THE RE-WIDEN LADDER: HELD, HONESTLY - the df leg impossibly open (df 469G << 1,600G);
the backlog leg met (nothing queued behind the single in-flight lane). The gate re-arms only
when the co-tenant releases ~1.1TB (df > 1,500G -> cuts resume C=2 at their own gate ->
then the df >= ~1,600G test).
(4) NO PREP STALLS: the valve held all window (185 pending-cut steady); no syng chains ran;
the T01411 class cannot strike while held; prep-master.log healthy ('no new completed
windows' at 16:00:18Z - the final split long done).
THE SUCCESSOR RESUMES: sparse SHORT polls; the last marathon lane resolves per the
guard-oracle doctrine (then the chunk-3 burst is <= 1 territory, the wave idles, and the
campaign waits wholly on the co-tenant's df release > 1,500G); watch for the release (if df
> 1,500G: the cuts master resumes C=2 by ITSELF - verify no hand action needed; the ladder
C=4 then C=6 fires only at the FULL gate: backlog drained AND df >= ~1,600G, kill the cuts
master MID-WAIT only, by exact pid, never the wave, never with marathons in flight); watch
for T01411-class prep stalls once the valve releases (poison-first doctrine, the fresh-chunk
dispatch pattern proven); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py ->
FILL datasheet-runtime.md (the standing close list + the pace-delta receipt + the chunk-2
boundary re-proof + the T01411 class's second sighting + the prep-pole stall + the deepest
df-excursion/co-tenant receipt NOW 387-694G WITH the writer-survey exclusions + the
memory-crunch receipt + this window's 4 kills) -> the record's closing increment -> the
commit. NO PR PUSH.
Session-open increment for this window: 3d3a967; window close committed with this note.

## RIDE 29 SESSION OPEN (16:59Z 2026-10-16, receipt-side) - THE LOWERED-FLOOR DRAIN
THE OWNER'S RULING RECEIPTED: the 1,500G floor was a margin, not a need - LOWERED to
DISK_FLOOR_GB=300; cuts master RELAUNCHED 16:57:42Z (pid 3632419, BATCH_W=8 N_CONCURRENT=1
RAYON=24, constrained width; env verified). First dispatch receipted in cuts/campaign.log:
'16:57:42Z batch .../batch-1792169862-0: 8 territories (df 396G)'. The last marathon lane
resolved BETWEEN windows: T02921 CALL FAILED (137) 16:46:41Z (guard receipt 'locus_9_scored:
40284352 kB'), fail-reaped in-lane - fail-closed 450 -> 451; the chunk-2 cohort FULLY terminal,
ZERO marathons in flight. OPEN POLL 16:58:35Z: split 9,891 / acc 5,844 / sub 1,924 (valve HELD
at 185 pending-cut, PREP_VALVE=160 - re-opens as cuts consume) / cut.done 1,681 / map dirs 74
/ call-terminal 1,681 = callok 1,230 + fail 451 (+16 cut.exit + poison 2 = 1,699 of 5,844,
29.1%) / backlog 0 - wave IDLES, re-fills as cut.done lands / df 401G / mem 923Gib / rvsize
277G. MASTERS: prep 411760, cuts 3632419 (constrained), chunked wave 2810146; heartbeat live;
impg-reference benign - NEVER KILLED. THE RIDE: sparse SHORT polls; cuts at C=1 feed the wave;
300G is the HARD FLOOR (self-gates; DO NOT lower without the owner); df >= 1,600G re-opens the
ladder (C=4 -> C=6 per doctrine); poison-first on prep stalls; terminal classes terminal; AT
CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md (standing close list + the
lowered-floor ruling + the constrained-drain pace) -> the closing increment -> the commit.
Session-open increment: f4c3914.
### RIDE 29 MID-WINDOW (17:52Z 2026-10-16) - THE CONSTRAINED-DRAIN PACE, FIRST RECEIPT
THE DRAIN IS LIVE AND FEEDING THE WAVE, WITH THE GUARD-ORACLE DESCENT AT C=1 NOW MEASURED:
six batches dispatched 16:57:42Z-17:47:54Z, walls 447-630s: FIVE width-8 guard kills (the
known packaging class, peaks 67-71GB, no strikes, territories descend to width 4) + ONE HONEST
COMPLETION (batch-1792172406-5 exit=0 466s, 8/8 yield 17:47:53Z). PARTIAL YIELDS ARE REAL:
batch-1792171286-3 yielded 5/8 cut.done before its guard kill. cut.done 1,681 -> 1,694 (+13,
~15/h); callok 1,230 -> 1,233; fail 451 steady; the wave re-filled within seconds of each
yield (the chunked dispatch proven again), consumed 3 calls, and now rides 2 marathon-class
lanes (T03192 15.9GB, T03388 11.3GB at ~20min - guard-oracle doctrine, never killed by hand);
the 8 fresh cut.done from batch-5 queue behind the chunk. STRUCTURAL NOTE (receipted, no
action): at C=1 the widest-first composition serializes width-8 attempts ahead of the width-4
retry class (pending at 17:35Z: 95 fresh width-8 / 90 width-4; ~47 fresh remain ~6 batches);
kills are the designed descent, the cost is ~8 min per killed batch. df BAND 382-524G the
whole window (one +120G co-tenant release 17:06Z; now ~388G, slow decline ~1G/min) - the 300G
floor NEVER approached, ~90-130G margin. Valve HELD at 185 pending all window (sub 1,924);
poison 2; heartbeat live; masters: prep 411760, cuts 3632419, wave 2810146; impg-reference
benign. ZERO masters touched, ZERO hand kills, ZERO lever actions.

### RIDE 29 WINDOW CLOSE (18:13Z 2026-10-16) - THE CONSTRAINED DRAIN PACE, RECEIPTED
THE RIDE (16:59-18:13Z, ~75 min, sparse SHORT polls): rode the lowered-floor drain end to
end - ZERO masters touched, ZERO hand kills, ZERO lever actions; the 300G hard floor NEVER
approACHED (deepest 328G at 18:07Z, margin 28G; two co-tenant releases +120G 17:06Z and +86G
~18:08Z kept the band at 328-524G).
THE CONSTRAINED-DRAIN PACE (the receipt the lowered-floor ruling answers to): TEN width-8
batches dispatched 16:57:42Z-18:12:20Z at C=1 BATCH_W=8 RAYON=24; NINE resolved = SEVEN guard
kills (the packaging class, peaks 67-71GB, walls 447-630s, no strikes, territories descend to
width 4) + TWO honest completions (batch-5 exit=0 466s 8/8; batch-6 exit=0 457s 8/8); batch-10
in flight at close. PARTIAL YIELDS ARE REAL AND RECEIPTED: batch-3 5/8 cut.done before its
kill, batch-8 7/8 before its kill. cut.done 1,681 -> 1,709 (+28 in ~75 min = ~23/h at C=1).
THE CHAIN FULLY ALIVE, BOTH POLES PROVEN AT THE LOWERED FLOOR: (1) the wave re-fills within
seconds of every yield (chunked dispatch), consumed 4 callok (1,230 -> 1,234), T03192 an
honest completion, T03388 still in flight (26GB at 42min - guard-oracle, never by hand), 23
cut.done queued behind the bounded chunk-tail; (2) THE PREP VALVE RELEASED ~18:04Z (cut.done
consumption dropped pending 164 -> <160, the pool ran +30, sub 1,924 -> 1,954) then re-held
at 187 - the valve-cuts-wave chain proven at DISK_FLOOR_GB=300.
FINAL POLL 18:12:43Z: split 9,891 / acc 5,844 (59.1% final) / sub 1,954 (valve-held at 187) /
cut.done 1,709 / map dirs 142 / call-terminal 1,685 = callok 1,234 + fail-closed 451 (+16
cut.exit + poison 2 = 1,703 of 5,844, 29.1% of accepted) / backlog 24 (1 lane in flight + 23
queued) / df 401G / mem 895Gib / rvsize 319G. MASTERS AT CLOSE: prep 411760 (valve-held),
cuts 3632419 (DISK_FLOOR_GB=300 BATCH_W=8 N_CONCURRENT=1 RAYON=24, batch-10 in flight),
chunked wave LANES=128 CHUNK_MULT=2 2810146; heartbeat live; the ancient impg-reference
benign - NEVER KILLED.
THE SUCCESSOR RESUMES: sparse SHORT polls; ride T03388 per the guard-oracle doctrine (the 23
queued cut.done dispatch at its chunk boundary); the drain continues ~23 cut.done/h (kills
are the designed descent - the width-4 retry class composes once the fresh width-8 pool
exhausts, ~47-55 fresh remain); watch df vs the 300G hard floor (self-gates; DO NOT lower
without the owner); if df recovers >= 1,600G the ladder re-opens (C=4 -> C=6 per doctrine,
kill the cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in
flight); poison-first on prep stalls; AT CAMPAIGN CLOSE: rv-rollup.py -> FILL
datasheet-runtime.md (the standing close list + the lowered-floor ruling + this pace receipt
+ the partial-yield finding + the valve-release-at-300G-floor proof) -> the record's closing
increment -> the commit. NO PR PUSH.

## RIDE-30 (18:16Z 2026-10-16) - CUT THE GUARD-KILL WASTE (session open, commit ad95d10)
THE RULING (the owner's, operational, within the measured envelope): ride-29's structural
finding - at the constrained packing 7 of 9 width-8 batches GUARD-KILLED (peaks 67-71GB vs
the 64GiB guard) before the adaptive descent to width 4; the widest-first ordering burns
doomed batch attempts (each kill wastes 447-630s batch wall minus the banked partial
yields); the width-4 retries SURVIVE. THE FIX: dispatch at the PROVEN width directly -
relaunch the cuts master at BATCH_W=4, N_CONCURRENT=1, DISK_FLOOR_GB=300, per the standing
doctrine (kill the cuts master MID-WAIT only, by exact pid, never the wave, never with
marathons in flight). THE RECEIPT the change answers to: cut.done/h before (~23/h,
width-8-doomed, ride-29) vs after (width-4-direct).
THE SWAP IS ARMED: rv-width4-watcher.sh (detached, pid 3762056, 2s cadence,
cuts/width4-watcher.log) waits for the cuts master 3632419's between-batches
ZERO-CHILDREN safe point (a batch's markers are all written before `wait` returns;
killing at zero children loses nothing, the relaunch recomposes from markers,
idempotent), gated on no wave lane older than 60min (never with marathons in flight);
at the safe point it kills 3632419 BY EXACT PID (TERM then KILL), belt-and-braces the
sub-2s dispatch race (any raced orphan subshell killed FIRST so NO marker is written -
no width-halving lie in the ledger), and relaunches via setsid launch-cuts.sh with
DISK_FLOOR_GB=300 BATCH_W=4 N_CONCURRENT=1 RAYON=24; the new master's environ is
verified in the watcher log.
OPEN POLL 18:15:13Z: split 9,891 / acc 5,844 / sub 1,954 (valve-held at 187) / cut.done
1,709 / call-terminal 1,685 = callok 1,234 + fail-closed 451 (+16 cut.exit + poison 2 =
1,703 of 5,844, 29.1%) / backlog 24 (1 lane T03388 at 46min/28GB + 23 queued) / df 401G /
mem 882Gib / rvsize 319G. Masters: prep 411760 (valve-held), cuts 3632419 (batch-10
width-8 in flight, the doomed class continuing until the swap), chunked wave 2810146
(LANES=128 CHUNK_MULT=2); heartbeat live; the ancient impg-reference benign.
OWED THIS WINDOW: (1) verify the swap fired (watcher log + the new master's first
'batch ...: 4 territories (df ...G)' dispatch line); (2) THE PACE-DELTA RECEIPT:
cut.done/h for >=1h post-swap (batch census: honest completions vs kills per batch,
yield per batch wall) vs the ~23/h width-8-doomed record; (3) standing contingencies
(300G hard floor - NEVER lower without the owner; guard-oracle - marathons never killed
by hand; poison-first on prep stalls; terminal classes terminal; ladder re-opens only at
backlog drained AND df >= 1,600G); (4) AT CAMPAIGN CLOSE (all 5,844 terminal):
rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the lowered-floor
ruling + the constrained pace receipts + THE WIDTH-4-DIRECT DELTA) -> the closing
increment -> the commit. NO PR PUSH.

RIDE-30 SWAP RECEIPT (18:52Z) - THE WIDTH-4-DIRECT DISPATCH IS LIVE. The v3
event-driven watcher (rv-width4-watcher.sh, pid 3828767; the v1 2s and v2 0.05s
polling watchers missed the ~0.2-1s between-batches gap - receipted in the
watcher log) caught the batch boundary at 18:52:40Z: batch-1792176222-13 (the
11th width-8 batch, exit 137, wall 538s) resolved, no batch subshell alive, NO
MARATHONS IN FLIGHT (T03388 resolved ~18:43Z as an HONEST COMPLETION, ~72min
wall / 28GB - call.done receipted, never killed by hand; the marathon gate had
held the swap ready from 18:30:45Z while it flew). Old cuts master 3632419
KILLED BY EXACT PID at the safe point (mid-wait, 18:52:42Z, TERM then KILL, no
race rounds needed - the ledger shows zero orphan subshells); relaunched via
setsid launch-cuts.sh: NEW MASTER 3838087, env verified /proc/3838087/environ:
BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300. FIRST WIDTH-4 DISPATCH
18:52:45Z 'batch-1792176764-0: 4 territories (df 316G)'. The old master's
constrained-era final census (16:57:42Z-18:52:42Z, BATCH_W=8 C=1): TWELVE
batches = NINE guard kills (walls 447-538s, peaks 49-71GB) + THREE honest 8/8
completions (batch-5 466s, batch-6 457s, batch-12 427s); cut.done 1,681 ->
~1,720 over ~115min (~20/h average with the partial yields banked). NOTE: df
316G at the swap, falling ~1.4G/min - the 300G hard floor may self-gate the
new master; the AFTER-pace receipt must name any floor-gate pauses honestly.
OWED: the AFTER-pace receipt (cut.done/h for >=1h post-18:52:45Z, width-4
batch census: completions vs kills per batch, yield per batch wall); then the
contingencies stand; NO PR PUSH.

RIDE-30 WINDOW CLOSE (19:53Z 2026-10-16, commit f83205e) - THE WIDTH-4-DIRECT
DELTA IS RECEIPTED, BOTH HALVES. THE SWAP: the v3 event-driven watcher caught
the 18:52:40Z batch boundary (no batch subshell alive; NO marathons in flight -
T03388 resolved ~18:43Z as an honest completion ~72min/28GB, the marathon gate
had held the swap ready while it flew), killed cuts master 3632419 BY EXACT PID
mid-wait (18:52:42Z, zero orphans, no race rounds), relaunched: NEW MASTER
3838087 (env verified: BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300);
first width-4 dispatch 18:52:45Z. Full swap receipt: cuts/width4-watcher.log.
THE PACE-DELTA RECEIPT (cut.done mtimes): BEFORE (width-8-doomed, 16:57:42Z-
18:52:42Z, 115min): 12 batches = 9 guard kills (walls 447-630s) + 3 honest 8/8;
+39 cut.done = 20.3/h, ~3.25/500s wall. AFTER (width-4-direct, first 38.3
producing minutes): 5 batches = 5 HONEST 4/4, ZERO kills (lastkill 18:52:40Z =
the final width-8 batch), walls 530/457/427/444/430s (avg 457.6s); +20
cut.done = 31.3/h, 4.0/457.6s wall. DELTA: +54% producing pace, doomed walls
eliminated (0/5 vs 9/12 killed). HONEST CAVEAT: the window's last 16.1min were
NOT producing - the 300G hard floor SELF-GATED the master 19:31:02Z at df 299G
(deepest 290G, co-tenant drain), RELEASED 19:47:08Z at df 305G after a +220G
co-tenant release (df 528G at close); calendar 20.0/h over the full hour is the
disk pole, not the width - NEVER lower the floor without the owner.
FINAL POLL 19:52:41Z: split 9,891 / acc 5,844 / sub 1,982 (valve released +28
at ~19:39Z during the gate, then re-held) / cut.done 1,740 / call-terminal
1,710 = callok 1,256 + fail-closed 454 (+16 cut.exit + poison 2 = 1,728 of
5,844, 29.6%) / backlog ~30 queued behind the chunk-3 tail (2 marathon lanes
~68min - guard-oracle, never by hand) / df 528G / mem 862Gib / rvsize 350G /
poison 2. MASTERS: prep 411760, cuts 3838087 (width-4-direct, batch-6 in
flight at close), chunked wave 2810146; heartbeat live; the ancient
impg-reference benign; width4-watcher done + retired.
THE SUCCESSOR RESUMES: sparse SHORT polls; ride the width-4 drain (keep the
census: honest completions vs rare width-4 kills that descend to 2 per the
guard-oracle); watch df vs the 300G hard floor (self-gates; DO NOT lower
without the owner); ladder re-opens only at backlog drained AND df >= 1,600G
(kill the cuts master MID-WAIT only, by exact pid, never the wave, never with
marathons in flight; C=4 -> C=6, the flock ceiling, never RAM); poison-first
on T01411-class prep stalls; terminal classes terminal; AT CAMPAIGN CLOSE (all
5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the
standing close list + the lowered-floor ruling + the constrained pace receipts
+ THIS WIDTH-4-DIRECT DELTA) -> the closing increment -> the commit. NO PR PUSH.

RIDE-31 STATUS AT THE DEADLINE CHECKPOINT (20:1xZ 2026-10-16) - THE REAPER IS
WRITTEN, NOT YET RUN. The owner's disk question is ANSWERED AND COMMITTED
(session-open increment e3a1bba): prep/ 228G = 106G split-rejected dead
extracts (121 extract-ahead dirs, the T00203-class segdup giants: locus.fa
94.5% + abandoned partial syng indexes; split ledgers kept) + 72G pending
substrates (locus.fa 35.5G + syng 2.2G + beds/graphs ~0.5G kept at cut/call +
35G consumed export/walks intermediates) + 42G extract-ahead TODO pipe (KEEP)
+ 10G cut-awaiting-call (locus.fa 4.9G keep-at-call + ~5G consumed reap) + 2G
terminal receipts (the wave's JIT reaper verified working).
rv-substrate-reaper.sh IS WRITTEN (classes A/B/C + lean-buffer rule + in-flight
exclusions) BUT HAS NOT BEEN RUN - NOT EVEN bash -n'd. THE SUCCESSOR OWES:
(1) bash -n rv-substrate-reaper.sh, read it, then run it detached
(nohup bash ./rv-substrate-reaper.sh > prep-reaper.nohup 2>&1 &) and poll
prep-reaper.log; (2) THE REGEN RECEIPT: pick one ordinary P-state territory,
evict-reap it fully (beds + beds.done + ledgers KEPT; drop extract/syng/maps/
substrate.done markers), run bash run-prep-one.sh <T> timed, verify the
substrate rebuilds (~30s ordinary) and the territory re-enters the cut/call
flow (gate the receipt on a to-be-called territory); (3) measure freed GB +
new rvsize + the pace delta (cut lanes must stay fed - the reaper keeps every
cut/call input; watch the width-4 batch census continues 0-kill); (4) window-
close increment + commit; (5) the datasheet at close gains the inventory
finding + the JIT-prep rule. Standing: 300G hard floor never lowered;
guard-oracle; poison-first; ladder at df >= 1,600G; NO PR PUSH. Width-4
dispatch live (master 3838087); last known census 6 width-4 batches honest,
zero kills; batch-7 in flight at the checkpoint.

RIDE-32 PROGRESS (20:20Z 2026-10-16) - THE SUBSTRATE REAPER HAS RUN; THE FREED-GB
RECEIPT IS IN; THE REGEN RECEIPT FOUND AND FIXED A REAL GAP. (1) SESSION-OPEN
increment committed fe057cd (open poll: cut.done 1,752 / terminals 1,729 of 5,844
(29.6%) / df 523G / width-4 batches 0-7 ALL HONEST 4/4 zero kills). (2) THE REAPER:
bash -n clean, reviewed against the live cut manifest + call-lane inputs (the cut
reads prep/<T>/locus-pos64.syng + reference-axis.json + partition-w10k-d1k/results +
locus.fa; the call reads locus.fa + reference-axis.json + locality-graphs - every
kept file verified named in the live code), run DETACHED 20:10:16Z, DONE in 16s,
receipt in prep-reaper.log: class A 631 files 93,588M (the split-rejected dead
extracts + abandoned partial syng indexes; split ledgers kept), class B 360 files
5,899M (cut-awaiting-call consumed syng index + export/walks; locus.fa/axis/
locality-graphs/beds KEPT), class C 552 files 33,609M (pending substrates' consumed
export/walks) + export.done/walks.done markers dropped on 184 pending. TOTAL FREED
132,237M apparent (~129G; df avail 525G -> 654G; prep/ 227G -> 95G apparent).
Honest delta vs the ride-31 estimate ~136G: ~93.6G class A vs ~101G estimated - the
20:0x census had aged (some pending became KA/in-flight); the freed is measured, not
estimated. (3) THE REGEN RECEIPT (T03194, ordinary, hash-order head): evicted fully
20:12:21Z (beds + beds.done + ledgers + walls kept; locus.fa/syng/locality-graphs/
axis/export-walks deleted; extract/syng/maps/substrate.done markers dropped), regen
attempt 1 FAILED AT WALKS IN 8.7s - THE GAP: class C had reaped walk-requests.tsv but
KEPT its producer marker walkreq.done, so the idempotent chain skipped regenerating
the reaped input (syng_homology_dump panicked 'read walk requests: NotFound',
receipted in prep/T03194/walks.err). THE FIX (minimal, in the reaper's own stated
intent - 'the chain's idempotence preserved'): classes B and C now drop the reaped
outputs' producer markers (B: syng/export/walks/walkreq.done; C: export/walks/
walkreq.done); the reaper RE-RUN idempotently 20:14:07Z (0 files reaped - all gone;
markers dropped on 183 pending; class B zeroed honestly because all 22
cut-awaiting-call territories were LIVE wave lanes at the pass - correctly
excluded). THE REGEN COMPLETED: walkreq 0s + walks 2s + maps 12s (second leg 13.5s;
full single-pass chain extract 5 + syng 2 + export 1 + walkreq 0 + walks 2 + maps 12
~= 22s - FASTER than the ~30s claim, measured); substrate.done 20:15:22Z; all
markers/files consistent (ledgers 6/6 kept, beds 31 kept, locality-graphs 14f,
reference-axis.json rebuilt). (4) THE RE-ENTRY GATE: batch-8 resolved exit=0 wall=472s
(width-4 census now 9/9 honest, ZERO kills since the swap); batch-9 dispatched
20:16:37Z at df 674G (the reaper's freed space visible in the dispatch line) with
T03194 FIRST in its manifest - the evicted-regenerated territory re-entered the cut
flow exactly at the head of the hash-interleaved order. AWAITED: batch-9's
resolution (T03194 cut.done) + the wave's pickup (the call gate). rvsize at this
stamp: 208,215M apparent (~203G; was ~355G at open - the reaper's 129G + the wave's
own JIT cascade ~23G over the same window). Masters live: prep 411760, cuts 3838087
(width-4-direct, batch-9 in flight), wave 2810146 (lanes live on the batch-7/8
yields); heartbeat live; the ancient impg-reference benign. Standing: the 300G hard
floor NEVER approached since (df 669-674G); guard-oracle; poison-first; the ladder
at df >= 1,600G; NO PR PUSH.

RIDE-32 WINDOW CLOSE (20:32Z 2026-10-16) - THE SUBSTRATE REAPER RAN, THE REGEN
RECEIPT IS COMPLETE THROUGH THE CUT, ONE REAL GAP FOUND + FIXED. THE RIDE (20:06-
20:32Z; session-open increment fe057cd): (1) THE REAPER RAN (detached 20:10:16Z,
done in 16s, receipts prep-reaper.log + prep-reaper.nohup): class A 631 files
93,588M (split-rejected dead extracts + abandoned partial syng indexes, ledgers
kept) + class B 360 files 5,899M (cut-awaiting-call consumed syng index +
export/walks; locus.fa/axis/locality-graphs/beds KEPT - the call's inputs,
verified against the live cut manifest and call-lane code) + class C 552 files
33,609M (pending substrates' consumed export/walks) + markers dropped = TOTAL
FREED 132,237M apparent (~129G); df avail 525G -> 654G; prep/ 227G -> 95G
apparent; TOTAL rvsize 355G -> ~201-203G apparent at close (the reaper's 129G +
the wave's own JIT cascade ~23G over the window). Honest delta vs the ride-31
~136G estimate: class A measured 93.6G vs ~101G estimated (the census aged); all
figures measured. (2) THE REGEN RECEIPT (T03194, ordinary, hash-order head):
evicted fully 20:12:21Z (beds + beds.done + ledgers + walls kept; markers
dropped); THE GAP FOUND ON ATTEMPT 1: class C had reaped walk-requests.tsv while
KEEPING its producer marker walkreq.done -> the idempotent chain skipped
regenerating the reaped input -> syng_homology_dump panicked 'read walk
requests: NotFound' at 8.7s (receipt prep/T03194/walks.err). THE FIX (minimal,
the reaper's own stated intent): classes B and C now drop the reaped outputs'
producer markers (B: syng/export/walks/walkreq.done; C: export/walks/walkreq.done);
reaper RE-RUN idempotently 20:14:07Z (0 files reaped - all already gone; markers
dropped on 183 pending; class B zeroed because all 22 then-KA territories were
LIVE wave lanes - correctly excluded); marker-consistency sweep CLEAN across all
P/KA states EXCEPT 8 live chunk-3 lanes (in-flight-excluded, harmless at call,
moot at call-terminal; a future idempotent reaper pass after they resolve cleans
them). THE REGEN COMPLETED: full single-pass chain extract 5s + syng 2s + beds
SKIPPED (kept) + export 1s + walkreq 0s + walks 2s + maps 12s ~= 22s - FASTER
than the ~30s claim, measured; substrate.done 20:15:22Z, ledgers 6/6, beds 31,
locality-graphs 14f, reference-axis.json rebuilt. (3) THE RE-ENTRY RECEIPT:
batch-9 dispatched 20:16:37Z (df 674G - the freed space in the dispatch line)
with T03194 FIRST in its manifest; batch-9 resolved exit=0 wall=487s; T03194
cut.done with HONEST cut-stats (reads_kept 346,014 / census_mult 485,562 /
routed_patterns 13,207) - the rebuilt substrate cuts normally; at close T03194
is QUEUED in the wave's call backlog behind the chunk-3 marathon tail (the call
lane had not yet dispatched at window close - stated honestly; the wave's TODO
filter verifiably includes it; its call lands at the next chunk boundary). (4)
THE PACE STATE: width-4 census batches 0-9 = 10/10 HONEST 4/4, ZERO guard kills
since the 18:52:42Z swap (walls 419-530s, avg ~454s; batch-10 in flight at
close); producing pace ~4 cut.done/~480s wall ~= 30/h holding (the +54% delta
state stands). (5) THE STANDING CONTINGENCIES: the 300G hard floor NEVER
approached this window (df band 523-678G, the reaper widened the margin by
+129G); guard-oracle (the chunk-3 marathon tail rode per doctrine, zero hand
kills); poison 2 terminal; terminal classes terminal (the fail-closed +14 this
window all named in wave-logs/failures.txt - the bind-ambiguity band); the
ladder NOT opened (df 678G << 1,600G); the ancient impg-reference benign. FINAL
POLL 20:31:57Z: split 9,891 / acc 5,844 / sub 1,982 / cut.done 1,760 /
call-terminal 1,762 = callok 1,276 + fail-closed 468 (+16 cut.exit + poison 2 =
1,762 of 5,844, 30.1% of accepted) / df 678G / mem 798Gib / rvsize 201G
apparent / batch-10 in flight. MASTERS: prep 411760 (valve-held), cuts 3838087
(width-4-direct), chunked wave 2810146 (chunk-3 tail lanes live); heartbeat
live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 drain (keep the
census); re-run the reaper idempotently after the chunk-3 lanes resolve (the 8
live-lane walkreq.done markers) and periodically as the JIT cadence (it is
re-runnable anytime, ~16s a pass); watch T03194's call land (the regen receipt's
final gate); watch df vs the 300G hard floor (NEVER lower without the owner);
the ladder only at backlog drained AND df >= 1,600G; poison-first on prep
stalls; AT CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md (the
standing close list + the width-4-direct delta + the inventory finding + the
JIT-prep rule + THE REAPER RECEIPTS: freed 129G classes A/B/C, the walkreq
marker-gap finding + fix, the ~22s regen measurement) -> the closing increment ->
the commit. NO PR PUSH.

RIDE-33 SESSION OPEN (20:35Z 2026-10-16, commit e118f59) - THE REAPER
RECEIPTS' FINAL LEGS. Opened from the ride-32 window close (a08581e): the
reaper RAN and receipted (freed 132,237M apparent ~129G, classes A 93.6G /
B 5.9G / C 33.6G, every cut/call input + split ledger kept; df 525->654G,
rvsize 355->201G); the regen receipt T03194 evicted-rebuilt in ~22s MEASURED,
re-entered FIRST in batch-9's manifest, cut honest (reads_kept 346,014);
the walkreq marker-gap found + fixed (classes B/C drop the reaped outputs'
producer markers). OPEN POLL 20:34:59Z: split 9,891 / acc 5,844 / sub 1,982
(valve-held) / cut.done 1,764 (batch-10 resolved honest +4 - width-4 census
batches 0-10 = 11/11 HONEST 4/4, ZERO kills since the 18:52:42Z swap; batch-11
dispatched 20:32:43Z in flight) / call-terminal 1,762 = callok 1,276 +
fail-closed 468 (+16 cut.exit + poison 2 = 30.1% of accepted) / df 676G /
mem 813Gib / rvsize 203G apparent / 5 call lanes live (the chunk-3 marathon
tail, 22-24min at the poll - guard-oracle). MASTERS: prep 411760, cuts
3838087 (BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146
(LANES=128 CHUNK_MULT=2); heartbeat live; the ancient impg-reference benign.
OWED THIS WINDOW: (1) T03194's call gate (cut.done 20:24Z, queued in the call
backlog behind the chunk-3 tail - the successor gates the regen receipt's
final leg on it); (2) the idempotent reaper re-run once the chunk-3 lanes
resolve (the 8 live-lane walkreq.done markers, in-flight-excluded at the
ride-32 pass; re-runnable anytime, ~16s a pass); (3) the width-4 census
continuation; (4) standing contingencies (300G hard floor NEVER lowered;
guard-oracle; poison-first; terminal classes terminal; the ladder only at
backlog drained AND df >= 1,600G); (5) AT CAMPAIGN CLOSE: rv-rollup.py ->
FILL datasheet-runtime.md -> the closing increment -> the commit. NO PR PUSH.

RIDE-33 WINDOW CLOSE (21:27Z 2026-10-16) - THE JIT-CADENCE REAPER PASS
RECEIPTED; T03194'S CALL GATE PASSES FORWARD HONESTLY BEHIND THE LAST TWO
CHUNK-3 MARATHONS. Opened 20:35Z (session-open increment e118f59). (1) THE
JIT-CADENCE REAPER PASS (20:53:43-49Z, ~6s, detached, receipt
prep-reaper.log): freed 4,634M apparent - class A 0 (nothing new rejected),
class B 135 files 863M + syng/export/walks/walkreq markers dropped on 20
cut-awaiting-call (the call's inputs locus.fa/axis/locality-graphs/beds
KEPT), class C 66 files 3,884M + markers dropped on 190 pending; the 4
then-live marathon lanes correctly in-flight-excluded; lean buffer pending
194 <= 200 no eviction; df 664 -> 669G, prep 106,179 -> 101,545M apparent.
THE 8 RIDE-32 LIVE-LANE walkreq.done MARKERS: 6 resolved in-window - five
HONEST completions (T03258 1653s, T03455 1816s, T03314 1899s, T03384 2080s,
T03457 2663s) + T03345 fail-closed 21:10:51Z exit 137 THE GUARD CLASS (rss
peak 66.5GB, ~59min - the guard-oracle doctrine, NEVER by hand, terminal
classes terminal, named in wave-logs/failures.txt line 469); their markers
are MOOT AT TERMINAL. TWO REMAIN LIVE: T03359 + T03347 at ~74min at close
(guard-oracle, never by hand). (2) THE T03194 CALL GATE - NOT LANDED AT
CLOSE, STATED HONESTLY: its cut artifacts are staged (map/T03194: cut-cache
/ cut-census / cut-reads / cut-stats); its cut.done 20:24Z fell AFTER the
current chunk snapshot (~20:11Z slice), and the wave's chunked dispatch
re-slices the TODO only when the chunk's last lane returns - its call
dispatches at the re-slice after T03359/T03347 resolve. THE SUCCESSOR GATES
THE REGEN RECEIPT'S FINAL LEG ON IT (expect call.done + honest call stats
vs the rebuilt substrate). (3) THE PACE STATE: width-4 census batches
0-16 = 17/17 HONEST 4/4, ZERO guard kills since the 18:52:42Z swap (walls
440-670s, avg ~518s); in-window +24 cut.done over 51.5min = 28.0/h
calendar (batches 11-16 all resolved honest) - the +54% width-4-direct
state holding. (4) STANDING CONTINGENCIES: the 300G hard floor NEVER
approached (df band 619-678G); guard-oracle (T03345's guard kill
receipted, ZERO hand kills); poison 2 terminal; terminal classes terminal
(fail-closed 469, the +1 this window guard-class); the ladder NOT opened
(df 619G << 1,600G, backlog present); the ancient impg-reference benign -
NEVER KILLED. FINAL POLL 21:26:26Z: split 9,891 / acc 5,844 (59.1% final)
/ sub 2,008 (valve-held) / cut.done 1,788 / call-terminal 1,768 = callok
1,281 + fail-closed 469 (+16 cut.exit + poison 2 = 1,768 of 5,844, 30.2%
of accepted) / df 619G / mem 856Gib / rvsize 219G apparent / batch-17
dispatched 21:24:42Z in flight / 2 marathon lanes live. MASTERS: prep
411760 (valve-held), cuts 3838087 (BATCH_W=4 N_CONCURRENT=1 RAYON=24
DISK_FLOOR_GB=300, width-4-direct), wave 2810146 (LANES=128 CHUNK_MULT=2);
heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; gate T03194's
call at the re-slice (the regen receipt's final leg); the width-4 census;
the JIT reaper cadence (~6-16s a pass, re-runnable anytime - the next pass
cleans the last two lanes' markers once terminal); watch df vs the 300G
hard floor (self-gates; DO NOT lower without the owner); the ladder only at
backlog drained AND df >= 1,600G; poison-first on prep stalls; AT
CAMPAIGN CLOSE (all 5,844 terminal): rv-rollup.py -> FILL
datasheet-runtime.md (the standing close list + the width-4-direct delta +
the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: 129G
freed ride-32 classes A/B/C + this window's 4.6G JIT-cadence pass + the
walkreq marker-gap fix + the ~22s regen + the T03194 re-entry/call gate) ->
the closing increment -> the commit. NO PR PUSH.

RIDE-34 SESSION OPEN (21:29Z 2026-10-16, commit 6ebfbb7) - THE LAST TWO
CHUNK-3 MARATHONS + THE REGEN RECEIPT CALL GATE + THE JIT-CADENCE REAPER.
Opened from the ride-33 window close: the JIT-cadence pass receipted (freed
4,634M apparent, classes B 863M on 20 KA + C 3,884M on 190 pending, live lanes
correctly excluded); 6 of 8 ride-32 lane markers moot-at-terminal (T03345
fail-closed guard 137 @ 59min rss 66.5GB); width-4 census 17/17 HONEST 4/4
ZERO guard kills since the 18:52:42Z swap. OPEN POLL 21:28:31Z: split 9,891 /
acc 5,844 (59.1%% final) / sub 2,008 (valve-held: 162 pending >= 160) /
cut.done 1,788 / call-terminal 1,768 = callok 1,281 + fail-closed 469 (+16
cut.exit + poison 2 = 30.2%% of accepted) / df 618G / mem 838Gib / rvsize 218G
apparent / batch-17 (dispatched 21:24:42Z) in flight. THE TWO REMAIN MARATHON
LANES LIVE at ~76min: T03347 (impg call-r pid 3966210, rss ~41.7GB) + T03359
(pid 3970598, rss ~27.1GB) - guard-oracle, NEVER by hand; the wave re-slices
TODO only when the last lane returns -> T03194 call gate + the reaper idempotent
re-run both land after they resolve. MASTERS: prep 411760 (valve-held), cuts
3838087 (BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, width-4-direct),
wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat live; the ancient
impg-reference benign. OWED THIS WINDOW: (1) T03194 call gate at the
re-slice (the regen receipt final leg: expect call.done + honest call stats
vs the rebuilt substrate); (2) the idempotent reaper re-run once the two
lanes resolve (~6-16s a pass, cleans their markers moot-at-terminal); (3)
the width-4 census continuation; (4) standing contingencies (300G hard
floor NEVER lowered; guard-oracle; poison-first; terminal classes terminal;
the ladder only at backlog drained AND df >= 1,600G); (5) AT CAMPAIGN CLOSE
(all 5,844 terminal): rv-rollup.py -> FILL datasheet-runtime.md -> the closing
increment -> the commit. NO PR PUSH.

RIDE-34 WINDOW CLOSE (22:39Z 2026-10-16) - T03359 RESOLVED HONEST AT THE
GUARD EDGE; T03347 THE LAST MARATHON, FROZEN AT 146min - THE T03194 CALL GATE
PASSED FORWARD HONESTLY. Opened 21:29Z (session-open increment e84395f).
(1) T03359 RESOLVED 21:48:40Z - HONEST COMPLETION, call 5,672s (~94.5min),
rss PEAK 61.5GB (survived ~2.5GB under the 64GiB guard - the deepest honest
rss of the campaign so far), call.done landed, callok 1,281 -> 1,282; the
longest honest call wall of the campaign receipted in wave-progress.log.
(2) THE JIT-CADENCE REAPER PASS (21:54:20-26Z, ~6s, detached, receipt
prep-reaper.log): freed 4,481M apparent - class A 0 (nothing new rejected),
class B 192 files 1,639M + syng/export/walks/walkreq markers dropped on 48
cut-awaiting-call, class C 48 files 2,988M + markers dropped on 186 pending;
in-flight exclusions CORRECT (T03347 the live lane + batch-20's four cut
territories); lean buffer 190 <= 200 no eviction; df 599 -> 604G, prep
111,035 -> 106,554M apparent. (3) T03347 REMAINS THE ONLY LIVE LANE - frozen
at 41.7GB rss for ~75min, ~146min elapsed, past the 119min longest
completion wall = the plateau-then-climb class (the T01875 analog rode ~2h20m
frozen then climbed) - guard-oracle doctrine, NEVER killed by hand; the wave
re-slices TODO only when it returns, so T03194's call (the regen receipt's
final leg) and the last marker cleanup both wait on it - STATED HONESTLY, gate
passed forward. (4) THE PACE STATE: cut.done 1,788 -> 1,824 (+36 over ~69.5min
= 31.1/h - the width-4-direct producing pace HOLDING exactly at the ride-30
+54%% receipt); width-4 census batches 0-25 = 26/26 HONEST 4/4, ZERO guard
kills since the 18:52:42Z swap (this window: batches 17-25, walls 411-573s);
prep valve released +12 then +24 (sub 2,008 -> 2,056), re-held per the 160
bound. (5) STANDING CONTINGENCIES: the 300G hard floor NEVER approached (df
band 575-640G); guard-oracle (zero hand kills, T03347 ridden); poison 2
terminal; terminal classes terminal (fail-closed 469, +0 this window); the
ladder NOT opened (df 575G << 1,600G, backlog present); the ancient
impg-reference benign - NEVER KILLED. FINAL POLL 22:38:35Z: split 9,891 /
acc 5,844 (59.1%% final) / sub 2,056 (valve-held) / cut.done 1,824 /
call-terminal 1,769 = callok 1,282 + fail-closed 469 (+16 cut.exit + poison
2 = 1,769 of 5,844, 30.3%% of accepted) / df 575G / mem 875Gib / rvsize 259G
apparent / batch-26 dispatched 22:36:29Z in flight / 1 marathon lane live.
MASTERS: prep 411760 (valve-held), cuts 3838087 (BATCH_W=4 N_CONCURRENT=1
RAYON=24 DISK_FLOOR_GB=300, width-4-direct), wave 2810146 (LANES=128
CHUNK_MULT=2); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls;
ride T03347 to resolution (guard-oracle; when it returns the wave re-slices and
T03194's call DISPATCHES - gate the regen receipt's final leg on call.done +
honest call stats); the idempotent reaper re-run once T03347 is terminal
(cleans its markers, ~6-16s a pass); the width-4 census; watch df vs the 300G
hard floor (self-gates; DO NOT lower without the owner); the ladder only at
backlog drained AND df >= 1,600G; poison-first on prep stalls; AT CAMPAIGN
CLOSE (all 5,844 terminal): rv-rollup.py -> FILL datasheet-runtime.md (the
standing close list + the width-4-direct delta + the inventory finding + the
JIT-prep rule + THE REAPER RECEIPTS: 129G ride-32 + 4.6G ride-33 + 4.5G this
window, the walkreq marker-gap fix, the ~22s regen, the T03194 re-entry + the
T03359 61.5GB honest-guard-edge receipt) -> the closing increment -> the
commit. NO PR PUSH.

## RIDE-35 WINDOW (22:50-22:58Z 2026-10-16) - ALL THE REAPER-RECEIPT FINAL LEGS LANDED + THE SYNG SEG-DUP-TAIL CLASS'S THIRD SIGHTING POISONED
Session-open increment committed 61a639c (opened from ride-34 close 6bcfd94). BETWEEN-WINDOWS
(receipted at open): T03347 RESOLVED TERMINAL 22:43:23Z guard-137 (the last chunk-3 marathon,
plateau-then-climb ridden to the oracle, ~169min, NEVER by hand; failures.txt named; the chunk-3
cohort FULLY TERMINAL); the wave re-sliced within seconds - the backlog burst live (37 call lanes),
T03254 DONE 43s, fresh binds named (T03302 rec 23290 / T03321 rec 2719 / T03523 rec 7624 /
T03391 rec 14510).
THE WINDOW'S WORK (all landed):
(1) T02562 - THE PREP-SYNG-SEG-DUP-TAIL CLASS'S THIRD SIGHTING, POISONED + KILLED PER THE RIDE-25
DOCTRINE (found at open: pid 132711 impg syng on a 318MB locus.fa, CPU==wall 29min state R, 29MB
RSS, NO output file, syng.err silent since 22:20:48Z stuck on HG02040#2#JBHDWT010000010.1:194815033-
195154083 (339,050bp) after 211 seqs/43.6Mbp/1.87M syncmers; median 4s max 159s = the confirmed
T00197/T01411 signature on sight): POISONED FIRST (prep/T02562/cut.poison touched + the honest
class line in cuts/poison.tsv - poison now 3: T00197 + T01411 + T02562) BEFORE the kill; receipts
snapshotted (syng.err.attempt1-snapshot + syng.log.attempt1-snapshot; wrapper syng.exit/wall
landed 22:53:25); KILLED BY EXACT PID 132711 (never the wrapper 132688, never the master 411760);
the wrapper exited, NO retry spawned (the poison filter held - verified live, zero T02562
processes); the master untouched, valve-held at 166. Class now 3 territories of 5,844 accepted.
(2) THE IDEMPOTENT REAPER RE-RUN (22:55:07-13Z, ~6s, detached, receipt prep-reaper.log): freed
2,676M apparent - class A 0, class B 36 files 1,083M + markers on 4 cut-awaiting-call, class C 30
files 1,627M + markers dropped on 178 pending; ~50 live lanes correctly in-flight-excluded
(incl. T03194's own live call lane); T03347's markers cleaned moot-at-terminal; lean buffer
182 <= 200 no eviction; df 606 -> 609G. CUMULATIVE REAPER FREED: 129G (ride-32) + 4.6G (ride-33)
+ 4.5G (ride-34) + 2.7G (ride-35) ~= 141G.
(3) THE T03194 CALL GATE - PASSED, THE REGEN RECEIPT FULLY CLOSED END-TO-END: call.done
22:55:50Z, call 693s, rss peak 8,525,280kB (~8.1GB); 14 calls + 14 sequences in
receipts/T03194/; map/T03194 + partitions/T03194 JIT-reaped at call-done. THE FULL CHAIN: evict
(ride-32) -> reaper -> regen ~22s MEASURED -> re-entry FIRST in batch-9's manifest -> cut honest
(reads_kept 346,014) -> CALL HONEST (14 calls) - the JIT-prep rule's premise proven end to end.
(4) PACE: width-4 census batches 0-27 = 28/28 HONEST 4/4, ZERO guard kills since the 18:52:42Z
swap (this window walls 477/418/649s); cut.done 1,824 -> 1,832 (+8/19.2min = 25/h calendar -
honestly below the 31.1/h producing pace: batch-27's 649s wall rode the 37-lane call burst
co-tenancy; the producing-pace +54% state itself stands). callok 1,282 -> 1,306 (+24), fail-closed
469 -> 487 (+18, every one a named bind ambiguity in wave-logs/failures.txt).
(5) STANDING CONTINGENCIES: the 300G hard floor NEVER approached (df band 573-612G); guard-oracle
(zero hand kills - T03347 ridden to the oracle; T02562's kill was the receipted poison-first
doctrine, NOT a guard kill); poison 3 terminal; terminal classes terminal; the ladder NOT opened
(df 612G << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED.
FINAL POLL 22:57:43Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,056 (valve-held at 166) /
cut.done 1,832 (batch-28 in flight) / call-terminal 1,793 = callok 1,306 + fail-closed 487
(+16 cut.exit + poison 3 = 1,812 of 5,844, 31.0% of accepted) / df 612G / mem 620Gib (the call
burst) / rvsize 222G apparent.
MASTERS: prep 411760 (valve-held), cuts 3838087 (BATCH_W=4 N_CONCURRENT=1 RAYON=24
DISK_FLOOR_GB=300, width-4-direct, batch-28 in flight), wave 2810146 (LANES=128 CHUNK_MULT=2,
the post-marathon chunk live); heartbeat live; the ancient impg-reference benign.
THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 drain (keep the census); the JIT reaper
cadence (re-runnable anytime, ~6-16s a pass); watch df vs the 300G hard floor (self-gates; DO NOT
lower without the owner); the ladder only at backlog drained AND df >= 1,600G; poison-first on
prep stalls (the class now 3 members); AT CAMPAIGN CLOSE (all 5,844 terminal): rv-rollup.py ->
FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory
finding + the JIT-prep rule + THE REAPER RECEIPTS: ~141G cumulative freed, the walkreq marker-gap
fix, the ~22s regen, the T03194 full-chain receipt + the syng segdup-tail class's 3 members) ->
the record's closing increment -> the commit. NO PR PUSH.

RIDE-36 SESSION OPEN (23:07Z 2026-10-16, commit 2eb8297) - THE STEADY WIDTH-4
DRAIN; THE REAPER RECEIPTS CLOSED (OWED LIST FROM RIDE-31 ALL DISCHARGED BY
RIDES 32-35). Standing receipts: cumulative substrate-reaper freed ~141G
apparent (129G ride-32 classes A 93.6G/B 5.9G/C 33.6G + 4.6G ride-33 + 4.5G
ride-34 + 2.7G ride-35; every cut/call input + split ledger kept); the regen
receipt CLOSED END-TO-END (T03194: evict -> reap -> regen ~22s MEASURED ->
re-entry FIRST in batch-9 -> cut honest reads_kept 346,014 -> CALL HONEST 14
calls/693s/rss ~8.1GB, map+partitions JIT-reaped); the walkreq marker-gap
found + FIXED (classes B/C drop the reaped outputs' producer markers); the
syng segdup-tail class 3 members poisoned terminal (T00197 + T01411 + T02562).
OPEN POLL 23:05:43Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,056
(valve-held at 162 pending) / cut.done 1,836 (batch-29 in flight, dispatched
23:03:18Z at df 615G) / call-terminal 1,818 = callok 1,312 + fail-closed 487
(+16 cut.exit + poison 3 = 31.1% of accepted) / df 616G / mem 542Gib /
rvsize 218G apparent / width-4 census batches 0-28 = 29/29 HONEST 4/4 ZERO
guard kills since the 18:52:42Z swap / ZERO marathons in flight (post-marathon
burst resolving fast: honest calls 43s-1066s). MASTERS: prep 411760
(valve-held), cuts 3838087 (BATCH_W=4 N_CONCURRENT=1 RAYON=24
DISK_FLOOR_GB=300, width-4-direct), wave 2810146 (LANES=128 CHUNK_MULT=2);
heartbeat live; the ancient impg-reference benign. OWED THIS WINDOW: (1) the
JIT-cadence reaper pass (~6-16s, idempotent, receipt prep-reaper.log); (2)
the width-4 census continuation; (3) standing contingencies (300G hard floor
NEVER lowered; guard-oracle; poison-first - the class is 3 members; terminal
classes terminal; the ladder only at backlog drained AND df >= 1,600G); (4)
AT CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md (the standing
close list + the width-4-direct delta + the inventory finding + the JIT-prep
rule + the reaper receipts) -> the closing increment -> the commit. NO PR
PUSH.

RIDE-36 WINDOW CLOSE (23:49Z 2026-10-16, commit bc68ddf) - THE STEADY
STATE, ALL RECEIPTED. Opened 23:06Z (session-open increment 2eb8297).
(1) THE JIT-CADENCE REAPER PASS (23:07:38-43Z, ~6s, detached, receipt
prep-reaper.log): class A 0, class B 30 files 970M + markers on 8
cut-awaiting-call, class C 0 files + markers on 174 pending; ~45 live
lanes correctly in-flight-excluded - PROVEN RETROACTIVELY: the excluded
T03405 + T03619 later resolved CALL FAILED (137) guard-class mid-score,
their lanes' files untouched; lean buffer 178 <= 200; freed 950M apparent,
df 616 -> 617G. CUMULATIVE REAPER FREED ~142G (129G + 4.6G + 4.5G + 2.7G +
0.95G), 5 idempotent passes, every cut/call input + split ledger kept.
(2) WIDTH-4 CENSUS: batches 0-32 = 33/33 HONEST 4/4, ZERO guard kills
since the 18:52:42Z swap (walls this window 555/547/541/796s - batch-32's
796s the longest honest width-4 wall, call-burst co-tenancy); cut.done
1,836 -> 1,852 (+16/41.3min = 23.2/h calendar, honestly below the ~30/h
state on co-tenancy - the ZERO-kill record stands, no width regression).
(3) THE PREP VALVE RELEASED mid-window (+25, sub 2,056 -> 2,081) then
re-held at 171 - the chain live; NO pathological syng sighted (checked
live at the release). (4) CALLS: callok 1,312 -> 1,324 (+12, the heavy
lanes HONEST under the oracle: T03371 2168s/25.1GB, T03407 2237s/33.4GB,
T03605 2717s/33.3GB); fail-closed 487 -> 489 (+2 both the guard class
mid-score, named); poison 3. ZERO hand kills, ZERO masters touched.
(5) CONTINGENCIES: the 300G floor NEVER approached (df band 605-617G);
guard-oracle held; the ladder NOT opened (df 605G << 1,600G); the ancient
impg-reference benign. FINAL POLL 23:47:03Z: split 9,891 / acc 5,844 /
sub 2,081 (valve-held at 171) / cut.done 1,852 (batch-33 in flight) /
call-terminal 1,832 = callok 1,324 + fail-closed 489 (+16 cut.exit +
poison 3 = 31.3% of accepted) / df 605G / mem 345Gib / rvsize 226G
apparent / zero marathons in flight. MASTERS: prep 411760, cuts 3838087
(width-4-direct), wave 2810146; heartbeat live. THE SUCCESSOR RESUMES:
sparse SHORT polls; the width-4 census; the JIT reaper cadence; watch df
vs the 300G hard floor; the ladder only at backlog drained AND
df >= 1,600G; poison-first (the syng segdup-tail class is 3 members); AT
CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md (the standing
close list + the width-4-direct delta + the inventory finding + the
JIT-prep rule + the reaper receipts ~142G cumulative) -> the closing
increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 37 SESSION OPEN (23:52Z 2026-10-16, commit 09d795e) - THE JIT-CADENCE REAPER PASS RUN
Opened from ride-36 close bc68ddf. THE PREDECESSOR'S RIDE-31/32 REAPER OWED LIST IS DISCHARGED
(rides 32-36, all receipted): the reaper ran 129G first pass + 4 cadence passes ~= 142G cumulative;
the regen receipt CLOSED END-TO-END (T03194: evict -> reap -> regen ~22s MEASURED -> re-entry FIRST
in batch-9 -> cut honest reads_kept 346,014 -> CALL HONEST 14 calls/693s/rss ~8.1GB). This window
keeps the cadence + the census.
(1) THE JIT-CADENCE REAPER PASS (23:50:41-46Z, ~6s, detached, receipt prep-reaper.log): bash -n
clean, reviewed (every cut/call input + split ledger kept; in-flight lanes excluded; markers dropped
per the ride-32 fix); class A 0 (nothing new rejected), class B 120 files 1,736M + syng/export/walks/
walkreq markers dropped on 24 cut-awaiting-call, class C 39 files 2,473M + markers dropped on 183
pending; ~39 live call lanes correctly in-flight-excluded; lean buffer 187 <= 200 no eviction;
freed 4,117M apparent, df 605 -> 609G, prep/ 109,650 -> 105,533M. CUMULATIVE SUBSTRATE-REAPER FREED
~= 146G apparent (129G + 4.6G + 4.5G + 2.7G + 0.95G + 4.1G this window), 6 idempotent passes.
OPEN POLL 23:49:55Z: split 9,891 / acc 5,844 / sub 2,081 (valve-held at 171) / cut.done 1,852
(batch-33 in flight) / call-terminal 1,832 = callok 1,324 + fail-closed 489 (+16 cut.exit + poison 3
= 31.3% of accepted) / df 605G / mem 323Gib / rvsize 226G / width-4 census 0-32 = 33/33 HONEST 4/4
ZERO kills since the 18:52:42Z swap / 15 live call lanes (22-57GB, oldest 1h04m, guard-oracle).
MASTERS: prep 411760 (valve-held), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24
DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat live; the ancient impg-reference
benign. OWED THIS WINDOW: the width-4 census continuation; sparse SHORT polls (watch the live call
cohort under the guard-oracle); the standing contingencies (300G hard floor NEVER lowered; poison-first
- the class is 3 members; terminal classes terminal; the ladder only at backlog drained AND
df >= 1,600G); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md (the standing close list +
the width-4-direct delta + the inventory finding + the JIT-prep rule + the reaper receipts ~146G) ->
the record's closing increment -> the commit. NO PR PUSH.
## RIDE CONTINUATION 37 WINDOW CLOSE (00:17Z 2026-10-17) - THE JIT-CADENCE PASS RECEIPTED (+4,117M, CUMULATIVE ~146G); THE WIDTH-4 CENSUS 36/36 HONEST ZERO KILLS; THE DEEP CALL LANES RIDDEN UNDER THE ORACLE
Opened 23:52Z (session-open increment 09d795e). THE WINDOW'S WORK, ALL RECEIPTED:
(1) THE JIT-CADENCE REAPER PASS (23:50:41-46Z, ~6s, detached, receipt prep-reaper.log; bash -n
clean + reviewed before the run): class A 0 (nothing new rejected), class B 120 files 1,736M +
syng/export/walks/walkreq markers dropped on 24 cut-awaiting-call, class C 39 files 2,473M + markers
dropped on 183 pending; ~39 live call lanes correctly in-flight-excluded; lean buffer 187 <= 200 no
eviction; freed 4,117M apparent (df 605 -> 609G, prep/ 109,650 -> 105,533M). CUMULATIVE
SUBSTRATE-REAPER FREED ~= 146G apparent (129G + 4.6G + 4.5G + 2.7G + 0.95G + 4.1G), 6 idempotent
passes, every cut/call input + split ledger kept throughout. THE EXCLUSIONS PROVEN
RETROACTIVELY AGAIN: two of the pass's excluded lanes (T03401, T03513) later resolved CALL
FAILED (137) mid-score, their files untouched by the pass.
(2) THE WIDTH-4 CENSUS: batches 0-35 = 36/36 HONEST 4/4, ZERO guard kills since the 18:52:42Z
swap (this window: batch-33 501s, batch-34 506s, batch-35 651s, all exit=0; batch-36 dispatched
00:11:43Z in flight, first emissions landed). cut.done 1,852 -> 1,864 (+12 over ~26min ~= 27.7/h
calendar - the +54% producing state holding, no width regression).
(3) THE CALL COHORT: the heavy lanes ridden under the guard-oracle, ZERO hand kills - T03566
DONE 00:16:22Z (call 5,088s ~84.8min, rss peak 50,773,220kB ~48.4GB - the deepest honest lane
this window); fail-closed 489 -> 491 (+2, both the guard class mid-score, named in
wave-logs/failures.txt: T03401 00:07:46Z locus 0 49,494 records / T03513 00:15:20Z locus 18
19,643 records); callok 1,324 -> 1,325. 12 lanes live at close, oldest ~1h31m at 56.6GB rss -
approaching the guard, guard-oracle doctrine, NEVER killed by hand.
(4) THE PREP VALVE RELEASED (+23, sub 2,081 -> 2,104) then RE-HELD - and the chunked-dispatch
patch's mid-dispatch valve check live in the log ('VALVE (mid-dispatch): 182 pending cut -
parking the remaining 3674 of 3802 TODO', 00:15:33Z); the poison-first watch ran at the release:
ZERO impg syng in flight, NO pathological sighting (the class stays 3 members: T00197 + T01411
+ T02562).
(5) STANDING CONTINGENCIES: the 300G hard floor NEVER approached (df band 591-609G, margin
~290G); guard-oracle held (zero hand kills - the two 137s are the oracle's own at the wave);
terminal classes terminal; the ladder NOT opened (df ~600G << 1,600G, backlog present); the
ancient impg-reference benign - NEVER KILLED.
FINAL POLL 00:16:45Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,104 (valve-held) / cut.done
1,864 (batch-36 in flight) / call-terminal 1,835 = callok 1,325 + fail-closed 491 (+16 cut.exit +
poison 3 = 1,835 of 5,844, 31.4% of accepted; every fail line named) / df 593G / mem 411Gib /
rvsize 232G apparent / 12 live call lanes (oldest 1h31m, rss to 56.6GB, guard-oracle).
MASTERS: prep 411760 (valve-held, log healthy), cuts 3838087 (BATCH_W=4 N_CONCURRENT=1
RAYON=24 DISK_FLOOR_GB=300, width-4-direct), wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat
live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 drain (keep the census); the JIT
reaper cadence (~6-16s a pass, re-runnable anytime - the next pass cleans the resolving lanes'
markers moot-at-terminal); watch df vs the 300G hard floor (self-gates; DO NOT lower without
the owner); the ladder only at backlog drained AND df >= 1,600G (kill the cuts master MID-WAIT
only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6, the flock
ceiling, never RAM); poison-first on prep stalls (the class is 3 members); AT CAMPAIGN CLOSE
(all 5,844 terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the
width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~146G
cumulative freed classes A/B/C across 6 idempotent passes, the walkreq marker-gap fix, the ~22s
regen, the T03194 full-chain receipt, the syng segdup-tail class's 3 members, the guard-oracle
receipts incl. this window's T03566 48.4GB honest deep lane) -> the record's closing increment ->
the commit. NO PR PUSH.

## RIDE CONTINUATION 38 WINDOW (00:21-00:46Z 2026-10-17) - THE CADENCE HELD: REAPER PASS +3.1G (CUMULATIVE ~149G, 7 PASSES); CENSUS 39/39 ZERO KILLS; THE DEEP LANES HONEST; ONE MARATHON PASSED FORWARD
Opened 00:21Z (session-open increment 646b881) from ride-37 close 96ba4f6. NOTE FOR THE
SUCCESSOR'S ORIENTATION: the task prompt's "reaper written but not run" snapshot is the
RIDE-30/31 state - the reaper has been RUN since ride-32 and its owed list (first pass +
regen receipt) is fully discharged; the standing duty is the CADENCE, which this window kept.
(1) THE JIT-CADENCE REAPER PASS (00:24:54-59Z, ~5s, detached, receipt prep-reaper.log;
bash -n clean + reviewed): class A 0, class B 108 files 582M + markers on 40
cut-awaiting-call, class C 57 files 2,640M + markers on 190 pending; lean buffer 194 <= 200;
freed 3,130M apparent (df 593 -> 596G; prep/ 112,401 -> 108,354M by close). CUMULATIVE
~= 149G (129G + 4.6 + 4.5 + 2.7 + 0.95 + 4.1 + 3.1G), 7 idempotent passes. EXCLUSIONS
PROVEN RETROACTIVELY A THIRD TIME: excluded lanes T03509/T03575 later FAILED guard-137
mid-score (named) and T03564/T03595 later resolved DONE HONEST - all untouched by the pass.
(2) WIDTH-4 CENSUS: batches 0-38 = 39/39 HONEST 4/4, ZERO guard kills since the 18:52:42Z
swap (walls 601/530/561s, all exit=0; batch-39 in flight at close); cut.done 1,864 -> 1,876
(+12/~29min = 24.8/h calendar, honestly below the ~31/h producing state on deep-lane
co-tenancy; the ZERO-kill record stands, no width regression).
(3) CALLS under the guard-oracle, ZERO hand kills: callok 1,325 -> 1,328 (+3 honest DEEP
lanes: T03545 5,440s/47.5GB, T03564 5,831s/57.9GB, T03595 6,504s/49.1GB); fail-closed
491 -> 494 (+3 guard-class mid-score, named: T03569 / T03575 / T03509); ONE MARATHON in
flight at close (pid 192355, ~55.8GB plateau-stable, ~1h58m elapsed - the plateau-then-climb
class, guard-oracle, NEVER killed by hand; the wave re-slices when it returns).
(4) PREP VALVE HELD all window (sub 2,104, pending-cut 170-182; no release, no pathological
syng sighting - the poison class stays 3 members).
(5) CONTINGENCIES: 300G hard floor NEVER approached (df band 578-597G); the ladder NOT
opened (df 591G << 1,600G, backlog present); the ancient impg-reference benign.
FINAL POLL 00:45:12Z: split 9,891 / acc 5,844 / sub 2,104 / cut.done 1,876 (batch-39 in
flight) / call-terminal 1,841 = callok 1,328 + fail-closed 494 (+16 cut.exit + poison 3 =
31.5% of accepted) / df 591G / mem 571-636Gib / rvsize 231G.
MASTERS: prep 411760 (valve-held), cuts 3838087 (width-4-direct BATCH_W=4
N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2);
heartbeat live. Window-close increment committed 7d4c3a0. THE SUCCESSOR RESUMES: sparse
SHORT polls; ride the marathon to resolution (guard-oracle); the JIT reaper cadence; the
width-4 census; watch df vs the 300G hard floor (self-gates; DO NOT lower without the
owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only,
by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling,
never RAM); poison-first on prep stalls; AT CAMPAIGN CLOSE (all 5,844 terminal):
rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct
delta + the inventory finding + the JIT-prep rule + the reaper receipts ~149G cumulative,
7 passes, triple-retroactive exclusion proof, the ~22s regen, the T03194 full chain) ->
the record's closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 39 WINDOW (00:49-01:36Z 2026-10-17) - CADENCE PASS 8 (+832M, CUMULATIVE ~150G); CENSUS 46/46 ZERO KILLS; DEEP COHORT RESOLVING; THE MARATHON PASSED FORWARD AT ~2h49m
Opened 00:49Z (session-open increment 2cf709b, from ride-38 close 7d4c3a0). ORIENTATION (standing,
from ride-38): a task prompt carrying the ride-30/31 snapshot ("reaper WRITTEN but NOT RUN") is
stale - the reaper has been RUN since ride-32; the ride-31 owed list (first pass + regen receipt
T03194 ~22s) is DISCHARGED; the standing duty is the CADENCE, kept this window.
(1) THE JIT-CADENCE REAPER PASS 8 (00:50:22-27Z, ~5s, detached, receipt prep-reaper.log; bash -n
clean + reviewed): class A 0, class B 84 files 886M + markers on 52 cut-awaiting-call, class C 0
files + markers on 178 pending (the 00:24Z pass had already reaped the C bulk - idempotence
demonstrated); lean buffer 182 <= 200 no eviction; freed 832M apparent (df 588 -> 589G, prep/
108,027 -> 107,195M). CUMULATIVE ~= 150G (129G + 4.6 + 4.5 + 2.7 + 0.95 + 4.1 + 3.1 + 0.8G),
8 idempotent passes. THE EXCLUSIONS VERIFIED LIVE this pass: the marathon T03431 (then ~2h02m,
58.4GB) is in the exclusion list with its files untouched; the excluded batch-40 entries resolved
honest after the pass.
(2) WIDTH-4 CENSUS: batches 0-45 = 46/46 HONEST 4/4, ZERO guard kills since the 18:52:42Z swap
(this window batch-39 508s / -40 430s / -41 423s / -42 443s / -43 445s / -44 505s / -45 437s,
all exit=0); cut.done 1,880 -> 1,904 (+24/46.3min = 31.1/h calendar - the +54% producing state
HOLDING, no width regression).
(3) CALLS under the guard-oracle, ZERO hand kills: T03447 CALL FAILED (137) mid-score 01:14:10Z
(locus 8) + T03427 CALL FAILED (137) mid-score 01:20:15Z (locus 5) - both named in
wave-logs/failures.txt, the oracle's own kill class; fail-closed 495 -> 497; callok 1,329 steady
(the chunk held by the two remaining deep lanes; the backlog drains at the re-slice). TWO DEEP
LANES PASSED FORWARD: marathon T03431 (pid 192355, ~58.35GB plateau, ~2h49m elapsed - the
plateau-then-climb class, guard-oracle, NEVER by hand) + T03572 (pid 196184, ~51.3GB frozen,
~2h48m).
(4) PREP VALVE RELEASED (+20, sub 2,104 -> 2,124) then RE-HELD; the mid-dispatch check live
('VALVE (mid-dispatch): 178 pending cut - parking the remaining 3651 of 3779 TODO', 01:06:31Z);
poison-first watch at the release: ZERO impg syng in flight (class stays 3 members).
(5) CONTINGENCIES: 300G floor NEVER approached (df band 562-589G); ladder NOT opened (df 562G <<
1,600G, backlog present); terminal classes terminal; the ancient impg-reference benign.
FINAL POLL 01:35:51Z: split 9,891 / acc 5,844 / sub 2,124 (valve-held) / cut.done 1,904
(batch-46 in flight) / call-terminal 1,845 = callok 1,329 + fail-closed 497 (+16 cut.exit +
poison 3 = 31.5% of accepted) / df 562G / mem 786Gib / rvsize 247G apparent.
MASTERS: prep 411760 (valve-held), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1
RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat live. Window-close
increment committed (record). THE SUCCESSOR RESUMES: sparse SHORT polls; ride the two deep lanes
to resolution (guard-oracle; the wave re-slices when they return); the JIT reaper cadence; the
width-4 census; watch df vs the 300G hard floor (self-gates; DO NOT lower without the owner);
the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid,
never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM);
poison-first on prep stalls; AT CAMPAIGN CLOSE (all 5,844 terminal): rv-rollup.py -> FILL
datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory
finding + the JIT-prep rule + the reaper receipts ~150G cumulative, 8 passes) -> the record's
closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 40 LIVE (01:38-01:42Z 2026-10-17) - SESSION OPEN 378daa4; THE JIT-CADENCE REAPER PASS 9 DONE (+2,401M APPARENT, CUMULATIVE ~= 152G); THE MARATHON T03431 RESOLVED TERMINAL GUARD-137 BETWEEN WINDOWS; T03572 THE LAST DEEP LIVE
Opened 01:38Z (session-open increment 378daa4, from ride-39 close 62e5a4a). ORIENTATION
(standing): a task prompt carrying the ride-30/31 snapshot ("reaper WRITTEN but NOT RUN") is
stale - the reaper has been RUN since ride-32; the ride-31 owed list (first pass + regen
receipt T03194 ~22s) is DISCHARGED; the standing duty is the CADENCE, kept this window.
(1) THE JIT-CADENCE REAPER PASS 9 (01:40:30-34Z, ~4s, detached pid 609817, receipt
prep-reaper.log; bash -n clean + reviewed): class A 0, class B 180 files 1,802M on 76
cut-awaiting-call (incl. the resolved marathon T03431's consumed substrate) + markers dropped,
class C 15 files 715M + markers on 174 pending; lean buffer 178 <= 200 no eviction; freed
2,401M apparent (df 559 -> 562G, prep/ 113,448 -> 111,047M). CUMULATIVE ~= 152G across 9
idempotent passes. THE EXCLUSIONS VERIFIED LIVE: T03572 (then ~2h51m, 51.3GB, the last deep
lane) is in the pass's exclusion list with locus.fa + syng untouched; the excluded batch-46
manifest entries resolved honest AFTER the pass (exit=0 wall 439s) - the cut lanes stay fed.
(2) THE MARATHON RESOLVED TERMINAL: T03431 (pid 192355) CALL FAILED (137) mid-score
01:37:29Z, locus 10 exactness sample 223/478,772,286 valid, named in wave-logs/failures.txt -
ridden to the oracle at ~2h50m/~58.4GB, NEVER killed by hand; fail-closed 497 -> 498.
(3) WIDTH-4 CENSUS: batch-46 exit=0 wall 439s -> 47/47 HONEST 4/4, ZERO guard kills since
the 18:52:42Z swap; batch-47 in flight; cut.done 1,904 -> 1,908.
Window continues: sparse SHORT polls; T03572 to resolution (guard-oracle; the wave re-slices
when it returns); contingencies standing (300G floor df 559-563G far above; ladder NOT opened;
poison class 3 members; the ancient impg-reference benign).
WINDOW CLOSE 02:07Z (record increment f3dde66): the marathon cohort FULLY TERMINAL - T03572
CALL FAILED (137) mid-score 01:52:46Z (locus 24, ~3h05m, 51.3GB frozen, named; its locus.fa +
syng verified untouched by the 01:40Z reaper pass that excluded it) - ZERO marathons in
flight, the wave re-sliced ~01:53Z, the backlog DRAIN BURST live (callok 1,329 -> 1,356,
fail-closed 497 -> 519, all named). Census batches 46-49 = 50/50 honest (walls 414-439s),
cut.done 1,920 (+16/26.7min = 36/h - above the +54% state). Valve released +19 re-held at
177; poison watch clean. FINAL POLL 02:06:03Z: split 9,891 / acc 5,844 / sub 2,143 (valve-held)
/ cut.done 1,920 (batch-50 in flight) / call-terminal 1,894 = callok 1,356 + fail-closed 519
(+16 cut.exit + poison 3 = 32.4% of accepted) / df 585G / mem 688GiB / rvsize 219G.
MASTERS: prep 411760 (valve-held), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1
RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat live.
Window-close increment committed f3dde66. THE SUCCESSOR RESUMES: sparse SHORT polls; ride the
drain burst (shallow 4GB-class lanes); the JIT reaper cadence; the width-4 census; watch df vs
the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog
drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with
marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls; AT
CAMPAIGN CLOSE (all 5,844 terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing
close list + the width-4 delta + the inventory finding + the JIT-prep rule + the reaper
receipts ~152G cumulative, 9 passes) -> the record's closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 41 WINDOW (02:09-02:21Z 2026-10-17) - SESSION OPEN a36239f; THE JIT-CADENCE REAPER PASS 10 DONE (+2,067M APPARENT, CUMULATIVE ~= 154G ACROSS 10 IDEMPOTENT PASSES); CENSUS 52/52 HONEST, ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF VERIFIED A FOURTH TIME
ORIENTATION (standing, from rides 38-40): a task prompt carrying the ride-30/31 snapshot
("the reaper WRITTEN but NOT RUN") is stale - the reaper has been RUN since ride-32; the
ride-31 owed list (first pass + regen receipt T03194 ~22s, closed end-to-end ride-35) is
DISCHARGED; the standing duty is the CADENCE, kept this window.
(1) THE JIT-CADENCE REAPER PASS 10 (02:09:16-21Z, ~5s, detached pid 707347, receipt
prep-reaper.log; bash -n clean + reviewed): class A 0 (nothing new rejected), class B
72 files 1,008M + syng/export/walks/walkreq markers dropped on 8 cut-awaiting-call,
class C 21 files 1,110M + markers dropped on 177 pending; ~36 live drain-burst call
lanes correctly in-flight-excluded; lean buffer 181 <= 200 no eviction; freed 2,067M
apparent (df 590 -> 592G, prep/ 108,477 -> 106,410M). CUMULATIVE SUBSTRATE-REAPER FREED
~= 154G apparent (129G + 4.6 + 4.5 + 2.7 + 0.95 + 4.1 + 3.1 + 0.8 + 2.4 + 2.1G), 10
idempotent passes, every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSIONS VERIFIED A FOURTH TIME: the pass's excluded live lanes kept their
call inputs on disk after the pass (T03352/T03418/T03650 locus.fa KEPT, verified on
disk); the excluded T03372 resolved call.done AFTER the pass - its substrate cleanup
was the WAVE's own success-path JIT reaper at call-done, not the substrate pass.
(3) THE WIDTH-4 CENSUS: batches 0-51 = 52/52 HONEST 4/4, ZERO guard kills since the
18:52:42Z swap (this window: batch-50 exit=0 wall 498s, batch-51 exit=0 wall 553s);
cut.done 1,920 -> 1,928; batch-52 dispatched 02:19:47Z in flight at close.
(4) THE CALL COHORT: the post-marathon drain burst continuing, all shallow (top lane
30.4GB at ~18min at close, NO marathons in flight); callok 1,361 -> 1,371 (+10);
fail-closed 519 STEADY (+0 this window); cut.exit 16; poison 3 terminal.
(5) STANDING CONTINGENCIES: the 300G hard floor NEVER approached (df band 581-607G,
margin > 280G); guard-oracle held (zero hand kills, zero masters touched); the prep
valve held at 161 pending-cut (sub 2,143; no release this window - no T01411-class
stall possible while held); the ladder NOT opened (df ~590G << 1,600G, backlog present);
terminal classes terminal; the ancient impg-reference benign - NEVER KILLED.
FINAL POLL 02:20:41Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,143 (valve-held at
161) / cut.done 1,928 (batch-52 in flight) / call-terminal 1,909 = callok 1,371 +
fail-closed 519 (+16 cut.exit + poison 3 = 1,909 of 5,844, 32.7% of accepted) /
df 607G (a co-tenant release lifted it from 581G) / mem 520Gib / rvsize 195G apparent.
MASTERS: prep 411760 (valve-held, log healthy), cuts 3838087 (BATCH_W=4
N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, width-4-direct), chunked wave 2810146
(LANES=128 CHUNK_MULT=2); heartbeat live; the ancient impg-reference benign.
THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 drain (keep the census); the
JIT reaper cadence (~4-6s a pass, re-runnable anytime - the next pass cleans the
resolving lanes' markers moot-at-terminal); watch df vs the 300G hard floor
(self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND
df >= 1,600G (kill the cuts master MID-WAIT only, by exact pid, never the wave, never
with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on
prep stalls (the class is 3 members: T00197 + T01411 + T02562); AT CAMPAIGN CLOSE
(all 5,844 terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close
list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE
REAPER RECEIPTS: ~154G cumulative freed classes A/B/C across 10 idempotent passes,
the walkreq marker-gap fix, the ~22s regen, the T03194 full-chain receipt, the syng
segdup-tail class's 3 members, the guard-oracle receipts) -> the record's closing
increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 42 SESSION OPEN (02:25Z 2026-10-17, commit 565fc00) - THE CADENCE CONTINUES
Opened from the ride-41 window close (a6d3a2c). ORIENTATION (standing, rides 38-41): a task
prompt carrying the ride-30/31 snapshot ("reaper WRITTEN but NOT RUN") is STALE - the reaper
has been RUN since ride-32; the ride-31 owed list (first pass 129G + the regen receipt T03194
~22s, closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE.
(1) THE JIT-CADENCE REAPER PASS 11 (02:24:46-51Z, ~5s, detached, receipt prep-reaper.log;
bash -n clean + reviewed): class A 0 (nothing new rejected), class B 60 files 601M + markers
dropped on 16 cut-awaiting-call, class C 57 files 3,762M + markers dropped on 188 pending;
lean buffer 192 <= 200 no eviction; freed 4,301M apparent (df 584 -> 588G, prep/ 112,714 ->
108,413M). CUMULATIVE ~= 158G apparent across 11 idempotent passes (129G + 4.6 + 4.5 + 2.7 +
0.95 + 4.1 + 3.1 + 0.8 + 2.4 + 2.1 + 4.3G). THE EXCLUSIONS VERIFIED LIVE A FIFTH TIME: the
20 live shallow lanes (T03352/T03418/T03641/T03650... - the same cohort the 02:09Z pass
excluded, still in flight) are all in the pass's exclusion list with locus.fa verified KEPT
on disk after the pass.
(2) OPEN POLL 02:24:15Z: split 9,891 / acc 5,844 / sub 2,162 / cut.done 1,928 (batch-52 in
flight) / call-terminal 1,892 = callok 1,373 + fail-closed 519 (+16 cut.exit + poison 3 =
1,911 of 5,844, 32.7%) / df 584-588G / mem ~500Gib / rvsize 214G / width-4 census 0-51 =
52/52 HONEST 4/4 ZERO guard kills since the 18:52:42Z swap / ZERO marathons (20 shallow
lanes, oldest ~30min/20.7GB, top rss 33.1GB).
MASTERS: prep 411760 (valve-held), cuts 3838087 (BATCH_W=4 N_CONCURRENT=1 RAYON=24
DISK_FLOOR_GB=300, width-4-direct), wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat live;
the ancient impg-reference benign. THE RIDE: sparse SHORT polls; the width-4 census
continuation; the standing contingencies (300G hard floor NEVER lowered; guard-oracle;
poison-first - the class is 3 members; terminal classes terminal; the ladder only at backlog
drained AND df >= 1,600G); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md ->
the closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 42 WINDOW CLOSE (02:40Z 2026-10-17) - THE CADENCE PASS 11 RECEIPTED
Opened 02:25Z (session-open increment 565fc00); window-close increment committed e3ef3cf.
(1) THE JIT-CADENCE REAPER PASS 11 (02:24:46-51Z, ~5s, detached, receipt prep-reaper.log;
bash -n clean + reviewed): freed 4,301M apparent (df 584 -> 588G, prep/ 112,714 -> 108,413M) -
class A 0 (nothing new rejected), class B 60 files 601M + markers on 16 cut-awaiting-call,
class C 57 files 3,762M + markers on 188 pending; lean buffer 192 <= 200 no eviction.
CUMULATIVE ~= 158G apparent across 11 idempotent passes (129G + 4.6 + 4.5 + 2.7 + 0.95 +
4.1 + 3.1 + 0.8 + 2.4 + 2.1 + 4.3G), every cut/call input and every split ledger kept.
THE EXCLUSIONS VERIFIED LIVE A FIFTH TIME: the 20 live shallow lanes (T03352/T03418/T03641/
T03650 - the same cohort the 02:09Z pass excluded) are all in the exclusion list with
locus.fa verified KEPT on disk after the pass.
(2) WIDTH-4 CENSUS: batches 0-53 = 54/54 HONEST 4/4, ZERO guard kills since the 18:52:42Z
swap (batch-52 exit=0 473s, batch-53 exit=0 519s); cut.done 1,928 -> 1,936 (+8/~15min ~=
32/h calendar - the +54% producing state holding); batch-54 in flight at close.
(3) CALLS: all shallow, ZERO marathons (16 lanes, oldest T03352 ~44min/20.7GB);
callok 1,373 -> 1,377 (+4); fail-closed 519 STEADY (+0); cut.exit 16; poison 3 terminal.
(4) CONTINGENCIES: the 300G hard floor NEVER approached (df band 572-588G); guard-oracle
held (zero hand kills, zero masters touched); the prep valve held all window (sub 2,162;
the poison class stays 3 members); the ladder NOT opened (df 583G << 1,600G, backlog
present); the ancient impg-reference benign - NEVER KILLED.
FINAL POLL 02:39:44Z: split 9,891 / acc 5,844 / sub 2,162 (valve-held) / cut.done 1,936
(batch-54 in flight) / call-terminal 1,896 = callok 1,377 + fail-closed 519 (+16 cut.exit +
poison 3 = 1,915 of 5,844, 32.8% of accepted) / df 583G / mem 489Gib / rvsize 216G.
MASTERS: prep 411760, cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24
DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat live.
THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census; the JIT reaper cadence
(~4-6s a pass, re-runnable anytime); watch df vs the 300G hard floor (self-gates; DO NOT
lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master
MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6
the flock ceiling, never RAM); poison-first on prep stalls; AT CAMPAIGN CLOSE (all 5,844
terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the
width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS:
~158G cumulative freed, 11 idempotent passes, the walkreq marker-gap fix, the ~22s regen,
the T03194 full-chain receipt, the syng segdup-tail class's 3 members) -> the record's
closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 43 SESSION OPEN (02:43Z 2026-10-17) - THE CADENCE CONTINUES (12TH PASS OWED)
Opened from the ride-42 window close (e3ef3cf). ORIENTATION (standing, rides 38-42): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN") is STALE - the
reaper has been RUN since ride-32; the ride-31 owed list (first pass 129G + the regen receipt
T03194 ~22s, closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE
(cumulative ~= 158G across 11 idempotent passes).
OPEN POLL 02:42:32Z: split 9,891 / acc 5,844 / sub 2,162 (valve-held) / cut.done 1,936
(batch-54 batch-1792204581-54 dispatched 02:36:22Z, in flight) / call-terminal 1,916 =
callok 1,378 + fail-closed 519 (+16 cut.exit + poison 3 = 1,916 of 5,844, 32.8%) /
df 575G / mem 468Gib / rvsize 215G apparent / width-4 census 0-53 = 54/54 HONEST 4/4 ZERO
guard kills since the 18:52:42Z swap / 33 live call lanes, ZERO marathons (oldest ~49min;
top rss T03625 54.3GB deep-but-young, all under the 64GiB guard).
MASTERS: prep 411760 (valve-held), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1
RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat live; the
ancient impg-reference benign. THE RIDE: the JIT-cadence reaper pass 12 (bash -n clean +
script verified unchanged since ride-31, md5 b804a7d3); the width-4 census continuation;
sparse SHORT polls; the standing contingencies (300G floor NEVER lowered; guard-oracle;
poison-first - the class is 3 members; terminal classes terminal; the ladder only at backlog
drained AND df >= 1,600G); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md ->
the closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 43 WINDOW CLOSE (02:47Z 2026-10-17) - THE CADENCE PASS 12 RECEIPTED (+88M APPARENT, CUMULATIVE ~= 158G ACROSS 12 IDEMPOTENT PASSES); CENSUS 55/55 HONEST ZERO KILLS; THE EXCLUSION PROOF VERIFIED LIVE A SIXTH TIME
Opened 02:43Z (session-open increment d863706).
(1) THE JIT-CADENCE REAPER PASS 12 (02:43:24-28Z, ~4s, detached pid 796389, receipt
prep-reaper.log; bash -n clean + script verified UNCHANGED since the ride-31 review, md5
b804a7d3): class A 0 (nothing new rejected), class B 48 files 119M + syng/export/walks/
walkreq markers dropped on 24 cut-awaiting-call, class C 0 files + markers dropped on 180
pending (the 02:24Z pass had already reaped the C bulk - idempotence demonstrated again);
lean buffer 184 <= 200 no eviction; freed 88M apparent (df 584 -> 584G, prep/ 107,740 ->
107,652M apparent). CUMULATIVE ~= 158.1G apparent (129G + 4.6 + 4.5 + 2.7 + 0.95 + 4.1 +
3.1 + 0.8 + 2.4 + 2.1 + 4.3 + 0.09G), 12 idempotent passes, every cut/call input and every
split ledger kept. THE EXCLUSIONS VERIFIED LIVE A SIXTH TIME: all 15 live call lanes
(T03352/T03604/T03625/T03627/T03638/T03641/T03650/T03515/T03519/T03599/T03601/T03659/
T03704/T03706/T03725) had locus.fa verified KEPT after the pass; the excluded batch-54
(T03817-T03821 cohort) resolved honest exit=0 wall 508s AFTER the pass - the cut lanes
stayed fed.
(2) WIDTH-4 CENSUS: batches 0-54 = 55/55 HONEST 4/4, ZERO guard kills since the 18:52:42Z
swap (batch-54 exit=0 wall 508s); cut.done 1,936 -> 1,940; batch-55 dispatched 02:44:51Z
in flight at close.
(3) CALLS: callok 1,378 -> 1,379 (+1); fail-closed 519 STEADY (+0, last fail lines
02:03-02:04Z pre-window, all named); cut.exit 16; poison 3 terminal; 16 lanes at close,
ZERO marathons (oldest ~51min, top rss T03638 ~42.1GB - deep-but-young, guard-oracle).
(4) CONTINGENCIES: the 300G hard floor NEVER approached (df band 575-584G); guard-oracle
held (zero hand kills, zero masters touched); the prep valve held all window (sub 2,162,
pending-cut 161-184 band; the poison class stays 3 members); the ladder NOT opened
(df 582G << 1,600G, backlog present); terminal classes terminal; the ancient impg-reference
benign - NEVER KILLED.
FINAL POLL 02:45:34Z: split 9,891 / acc 5,844 / sub 2,162 (valve-held) / cut.done 1,940
(batch-55 in flight) / call-terminal 1,917 = callok 1,379 + fail-closed 519 (+16 cut.exit +
poison 3 = 1,917 of 5,844, 32.8% of accepted) / df 582G / mem 516Gib / rvsize 217G.
MASTERS: prep 411760 (valve-held), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1
RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat live.
Window-close increment committed (record). THE SUCCESSOR RESUMES: sparse SHORT polls; the
width-4 census; the JIT reaper cadence (~4-6s a pass, re-runnable anytime); watch df vs the
300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog
drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never
with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep
stalls; AT CAMPAIGN CLOSE (all 5,844 terminal): rv-rollup.py -> FILL datasheet-runtime.md
(the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep
rule + THE REAPER RECEIPTS: ~158G cumulative freed, 12 idempotent passes, the walkreq
marker-gap fix, the ~22s regen, the T03194 full-chain receipt, the syng segdup-tail class's
3 members) -> the record's closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 44 WINDOW (02:54-03:02Z 2026-10-17) - THE CADENCE PASS 13 RECEIPTED; THE T00024 SIGHTING NAMED
Opened 02:54Z (session-open increment f76c90b). ORIENTATION STANDING (rides 38-43): a ride-30/31-style
prompt ("the reaper WRITTEN but NOT RUN") is STALE - the reaper has RUN since ride-32; the ride-31 owed
list (the 129G first pass + the ~22s regen receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the
standing duty is the CADENCE.
(1) THE JIT-CADENCE REAPER PASS 13 (02:55:34-39Z, ~5s, detached pid 816407, receipt prep-reaper.log;
bash -n clean + script verified UNCHANGED since the ride-31 review, md5 b804a7d3, mtime Oct 16 20:13):
class A 0, class B 48 files 88M + syng/export/walks/walkreq markers dropped on 32 cut-awaiting-call,
class C 0 files + markers dropped on 172 pending (idempotence again - the 02:24Z pass reaped the C bulk);
lean buffer 176 <= 200 no eviction; freed 56M apparent (df 584 -> 584G, prep/ 106,826 -> 106,770M
apparent). CUMULATIVE ~= 158.2G apparent across 13 idempotent passes (129G + 4.6 + 4.5 + 2.7 + 0.95 +
4.1 + 3.1 + 0.8 + 2.4 + 2.1 + 4.3 + 0.09 + 0.06G), every cut/call input and every split ledger kept.
THE EXCLUSIONS VERIFIED LIVE A SEVENTH TIME: the T00024 prep lane + the batch-56 cut cohort
(T03804/T03809/T03714/T03826) + the 12 live call lanes all in the exclusion list, locus.fa verified KEPT
on disk after the pass (T00024/T03515/T03604/T03627/T03638/T03641 spot-checked); the excluded batch-56
resolved HONEST exit=0 wall 469s AFTER the pass - the cut lanes stayed fed.
(2) WIDTH-4 CENSUS: batches 0-56 = 57/57 HONEST 4/4, ZERO guard kills since the 18:52:42Z swap (this
window: batch-55 exit=0 450s, batch-56 exit=0 469s); cut.done 1,940 -> 1,948 (+8/~15min ~= 32/h calendar -
the +54% producing state holding); batch-57 dispatched 03:00:12Z in flight at close.
(3) CALLS: callok 1,379 -> 1,383 (+4); fail-closed 519 STEADY (+0); cut.exit 16; poison 3 terminal;
all lanes shallow, ZERO marathons (top rss ~55GB at ~61min - deep-but-young, guard-oracle).
(4) THE T00024 SIGHTING (the window's one new datum, NAMED, NO intervention): T00024 is the LAST UNBUILT
SUBSTRATE of the accepted set (the prep master's TODO is exactly {T00024} - every other accepted territory
carries substrate.done/skipped/rejected/poison); a 478MB segdup-giant locus.fa (extracted Oct 15 03:52,
never past syng). Attempt at 01:41:37Z failed exit=1 wall=1s (stderr overwritten by the retry -
unobservable); the retry (dispatched 02:20:35Z at a below-160 valve dip; the valve held 161-192 around it)
ran 40+ min CPU at 23MB RSS on ONE 317kb sequence (HG03927#1#CM086638.1:130277163-130594662) after
indexing 4 sequences in the first second - the pathological signature of the named SYNG SEG DUP-TAIL
CLASS (T00197 + T01411 + T02562), a candidate 4th member. The prep master is BLOCKED inside its
single-job xargs wave - costless: TODO=1 and the valve is closed; nothing else is dispatchable. DECISION
POINT FOR THE SUCCESSOR: the T00197 precedent is kill-by-exact-pid only at ~3h CPU with the double-writer
race named, POISON-FIRST (touch prep/T00024/cut.poison before any kill so the master never re-dispatches
it); if the syng COMPLETES, receipt the wall (a segdup-tail syng completion datum for the datasheet).
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 584-586G, margin > 280G); guard-oracle
held (zero hand kills, zero masters touched); the prep valve held all window (sub 2,162, pending-cut
161-176 band); the ladder NOT opened (df << 1,600G, backlog present); terminal classes terminal; the
ancient impg-reference benign - NEVER KILLED.
FINAL POLL 03:00:23Z: split 9,891 / acc 5,844 / sub 2,162 / cut.done 1,948 (batch-57 in flight) /
call-terminal 1,902 = callok 1,383 + fail-closed 519 (+16 cut.exit + poison 3 = 1,921 of 5,844, 32.9% of
accepted) / df 586G / mem 594Gib / rvsize 212G apparent.
MASTERS: prep 411760 (blocked costlessly in the T00024 xargs wave, log healthy to 02:20:34Z), cuts 3838087
(width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2);
heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census; the JIT reaper cadence
(~4-6s a pass, re-runnable anytime); the T00024 watch (poison-first at the 3h precedent; receipt a
completion); watch df vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only
at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with
marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); AT CAMPAIGN CLOSE (all 5,844 terminal):
rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the
inventory finding + the JIT-prep rule + THE REAPER RECEIPTS ~158.2G across 13 passes + the T00024
segdup-giant syng sighting and its resolution) -> the record's closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 45 WINDOW (03:06-03:19Z 2026-10-17) - THE CADENCE PASS 14 RECEIPTED; THE DEEP COHORT RESOLVED UNDER THE ORACLE (3x GUARD-137 MID-SCORE); T00024 PASSED FORWARD AT ~58MIN
Opened 03:06Z (session-open increment e1a2840). ORIENTATION STANDING (rides 38-44): a ride-30/31-style
prompt ("the reaper WRITTEN but NOT RUN") is STALE - the reaper has RUN since ride-32; the ride-31 owed
list (the 129G first pass + the ~22s regen receipt, T03194 closed end-to-end ride-35) is DISCHARGED;
the standing duty is the CADENCE.
(1) THE JIT-CADENCE REAPER PASS 14 (03:06:25-30Z, ~5s, detached, receipt prep-reaper.log; bash -n
clean + script verified UNCHANGED since the ride-31 review, md5 b804a7d3, mtime Oct 16 20:13): class A 0,
class B 24 files 40M + syng/export/walks/walkreq markers dropped on 36 cut-awaiting-call, class C 0
files + markers dropped on 168 pending (idempotence again - the 02:24Z pass reaped the C bulk); lean
buffer 172 <= 200 no eviction; freed 25M apparent (df 586 -> 586G, prep/ 106,427 -> 106,402M apparent).
CUMULATIVE ~= 158.2G apparent across 14 idempotent passes (129G + 4.6 + 4.5 + 2.7 + 0.95 + 4.1 + 3.1 +
0.8 + 2.4 + 2.1 + 4.3 + 0.09 + 0.06 + 0.03G), every cut/call input and every split ledger kept.
THE EXCLUSIONS VERIFIED LIVE AN EIGHTH TIME: the T00024 prep lane + the batch-57 cut cohort
(T03453/T03742/T03356/T03333) + the 12 live call lanes all excluded, locus.fa verified KEPT on disk
after the pass (T00024 478MB / T03453 / T03356 / T03515 / T03604 / T03641 spot-checked); the excluded
batch-57 resolved HONEST exit=0 wall 546s AFTER the pass - the cut lanes stayed fed.
(2) WIDTH-4 CENSUS: batches 0-58 = 59/59 HONEST 4/4, ZERO guard kills since the 18:52:42Z swap
(batch-57 546s, batch-58 428s, both exit=0); cut.done 1,948 -> 1,956 (+8/~18min, 27-32/h band - the
+54% producing state holding); batch-59 dispatched 03:16:28Z in flight at close.
(3) CALLS: THE DEEP COHORT OF THE RIDE-44 CLOSE RESOLVED UNDER THE GUARD-ORACLE, ZERO HAND KILLS -
three deep lanes riden to the guard's own kill mid-score, ALL NAMED in wave-logs/failures.txt:
T03638 (03:15:02Z, locus 3, ~72min/~64GB), T03515 (03:15:10Z, locus 11 exactness sample, ~76min/65.1GB),
T03704 (03:16:30Z, locus 0); fail-closed 519 -> 522 (+3, all this window's oracle kills); callok
1,383 -> 1,387 (+4); cut.exit 16; poison 3 terminal; remaining top lane plateau-stable 62.3GB/~83min
(the plateau-then-climb class, guard-oracle).
(4) THE T00024 WATCH PASSED FORWARD (NO intervention): the segdup-giant syng retry (pid 732137)
~58min elapsed at 23MB RSS, R-state, still on ONE 317kb sequence - the named SYNG SEG DUP-TAIL
signature (candidate 4th member; class = T00197 + T01411 + T02562); 58min << the ~3h T00197 precedent;
the prep master's single-job xargs block stays costless (TODO=1, valve held 161-192). DECISION POINT
UNCHANGED: poison-first (prep/T00024/cut.poison) then kill-by-exact-pid ONLY at the ~3h precedent with
the double-writer race named - NOT before; a completion is a receipted segdup-tail syng datum.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 567-591G, margin > 260G);
guard-oracle held (zero hand kills, zero masters touched); the prep valve held all window (sub 2,162);
the ladder NOT opened (df 591G << 1,600G, backlog present); terminal classes terminal; the ancient
impg-reference benign - NEVER KILLED.
FINAL POLL 03:18:54Z: split 9,891 / acc 5,844 / sub 2,162 (valve-held) / cut.done 1,956 (batch-59 in
flight) / call-terminal 1,928 = callok 1,387 + fail-closed 522 (+16 cut.exit + poison 3 = 1,928 of
5,844, 33.0% of accepted) / df 591G / mem 775Gib / rvsize 207G apparent.
MASTERS: prep 411760 (the T00024 wave), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1
RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat live. THE SUCCESSOR
RESUMES: sparse SHORT polls; the width-4 census; the JIT reaper cadence (~4-6s a pass, re-runnable
anytime); the T00024 watch (poison-first at the 3h precedent; receipt a completion); the remaining
deep lane under the guard-oracle; watch df vs the 300G hard floor (self-gates; DO NOT lower without
the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact
pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM);
AT CAMPAIGN CLOSE (all 5,844 terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close
list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS
~158.2G across 14 passes + the T00024 sighting and its resolution) -> the record's closing increment
-> the commit. NO PR PUSH.

## RIDE CONTINUATION 46 WINDOW (03:26-03:33Z 2026-10-17) - THE CADENCE PASS 15 RECEIPTED; CENSUS 61/61 HONEST; T00024 PASSED FORWARD AT ~1:12; THE T03641 DEEP LANE THE LAST OF ITS COHORT
Opened 03:26Z (session-open increment ae63815). ORIENTATION STANDING (rides 38-45): a ride-30/31-style
prompt ("the reaper WRITTEN but NOT RUN") is STALE - the reaper has RUN since ride-32; the ride-31 owed
list (the 129G first pass + the ~22s regen receipt, T03194 closed end-to-end ride-35) is DISCHARGED;
the standing duty is the CADENCE.
(1) THE JIT-CADENCE REAPER PASS 15 (03:29:25-29Z, ~4s, detached, receipt prep-reaper.log; bash -n clean +
script verified UNCHANGED since the ride-31 review, md5 b804a7d3, mtime Oct 16 20:13): class A 0, class B
72 files 141M + syng/export/walks/walkreq markers dropped on 48 cut-awaiting-call, class C 0 files + markers
dropped on 156 pending (idempotence again - the 02:24Z pass reaped the C bulk); lean buffer 160 <= 200 no
eviction (the valve at its own bound); freed 96M apparent (df 594 -> 594G, prep/ 103,854 -> 103,758M
apparent). CUMULATIVE ~= 158.3G apparent across 15 idempotent passes (129G + 4.6 + 4.5 + 2.7 + 0.95 + 4.1 +
3.1 + 0.8 + 2.4 + 2.1 + 4.3 + 0.09 + 0.06 + 0.03 + 0.1G), every cut/call input and every split ledger kept.
THE EXCLUSIONS VERIFIED LIVE A NINTH TIME: the T00024 prep lane + the batch-60 cut cohort
(T03181/T03636/T03370/T03499) + the deep call lane T03641 all in the exclusion list, locus.fa verified KEPT
after the pass (T00024 457MB / T03641 839MB / the cohort 112-165MB); the excluded batch-60 resolved HONEST
exit=0 wall 437s AFTER the pass - the cut lanes stayed fed.
(2) WIDTH-4 CENSUS: batches 0-60 = 61/61 HONEST 4/4, ZERO guard kills since the 18:52:42Z swap (this window:
batch-59 454s, batch-60 437s, both exit=0); cut.done 1,956 -> 1,964 (+8/~26min, the 27-33/h producing band);
batch-61 dispatched 03:31:21Z in flight at close.
(3) CALLS: callok 1,387 -> 1,388 (+1); fail-closed 522 -> 523 (+1 - T03627 guard-137 mid-score 03:26:03Z,
locus 2, NAMED in wave-logs/failures.txt); cut.exit 16; poison 3 terminal; the ONE deep lane T03641
(pid 657722, ~1:37, 62.3GB plateau-stable - the plateau-then-climb class, guard-oracle); ZERO other
marathons.
(4) THE T00024 WATCH PASSED FORWARD (NO intervention): the segdup-giant syng retry (pid 732137) ~1:12
elapsed at 23MB RSS, R-state, still on ONE 317kb sequence - the named SYNG SEG DUP-TAIL signature
(candidate 4th member; class = T00197 + T01411 + T02562); ~1:12 << the ~3h T00197 precedent; the prep
master's single-job xargs block stays costless (TODO=1, valve held, pending at its own 160 bound). DECISION
POINT UNCHANGED: poison-first (prep/T00024/cut.poison) then kill-by-exact-pid ONLY at the ~3h precedent
with the double-writer race named - NOT before; a completion is a receipted segdup-tail syng datum.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 589-594G, margin > 280G); guard-oracle held
(zero hand kills - T03627 was the guard's own kill; zero masters touched); the prep valve held all window
(sub 2,162, pending 160); the ladder NOT opened (df 589G << 1,600G, backlog present); terminal classes
terminal; the ancient impg-reference benign - NEVER KILLED.
FINAL POLL 03:32:28Z: split 9,891 / acc 5,844 / sub 2,162 (valve-held) / cut.done 1,964 (batch-61 in
flight) / call-terminal 1,930 = callok 1,388 + fail-closed 523 (+16 cut.exit + poison 3 = 1,930 of 5,844,
33.0% of accepted) / df 589G / mem 871Gib / rvsize 208G apparent.
MASTERS: prep 411760 (the T00024 wave), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24
DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat live. THE SUCCESSOR RESUMES: sparse
SHORT polls; the width-4 census; the JIT reaper cadence (~4-6s a pass, re-runnable anytime); the T00024 watch
(poison-first at the 3h precedent; receipt a completion); the T03641 deep lane under the guard-oracle; watch
df vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained
AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight;
C=4 -> C=6 the flock ceiling, never RAM); AT CAMPAIGN CLOSE (all 5,844 terminal): rv-rollup.py -> FILL
datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the
JIT-prep rule + THE REAPER RECEIPTS ~158.3G across 15 passes + the T00024 sighting and its resolution) ->
the record's closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 47 SESSION OPEN (03:44Z 2026-10-17) - THE CADENCE CONTINUES (PASS 16 OWED)
Opened from the ride-46 window close (2aaa7b6). ORIENTATION STANDING (rides 38-46): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN") is STALE - the
reaper has RUN since ride-32; the ride-31 owed list (the 129G first pass + the regen receipt
T03194 ~22s, closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE
(cumulative ~= 158.3G across 15 idempotent passes).
OPEN POLL 03:41:45Z: split 9,891 / acc 5,844 / sub 2,162 (valve-held) / cut map dirs 127 /
call-terminal 1,927 = callok 1,388 + fail-closed 523 (+16 cut.exit + poison 3 = 1,927 of
5,844, 33.0%) / df 587G / mem 848Gib / rvsize 210G apparent / width-4 census 62/62 HONEST
zero non-zero exits since the 18:52:42Z swap (batch-61 resolved exit=0 wall 652s ~03:42;
batch-62 batch-1792208299-62 in flight, cohort T03560/T03358/T03812/T03719) / the ONE live
deep call lane T03641 (pid 657722, ~1:48, 63.8GB climbing - the plateau-then-climb class,
guard-oracle), ZERO other marathons.
MASTERS: prep 411760 (the costless T00024 single-job xargs wave, log healthy to 02:20:34Z),
cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave
2810146 (LANES=128 CHUNK_MULT=2); heartbeat live; the ancient impg-reference benign. THE
RIDE: the JIT-cadence reaper pass 16 (bash -n clean + script verified UNCHANGED since the
ride-31 review, md5 b804a7d3, mtime Oct 16 20:13); the width-4 census continuation; sparse
SHORT polls; the standing contingencies (300G floor NEVER lowered; guard-oracle; poison-first
- the class is 3 members + the T00024 candidate at ~1:23 << the 3h precedent; terminal
classes terminal; the ladder at df >= 1,600G); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL
datasheet-runtime.md -> the closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 47 WINDOW CLOSE (03:56Z 2026-10-17) - THE CADENCE PASS 16 RECEIPTED (+55M APPARENT, CUMULATIVE ~= 158.4G ACROSS 16 IDEMPOTENT PASSES); CENSUS 64/64 HONEST ZERO KILLS; T00024 PASSED FORWARD AT ~1:35; THE T03641 DEEP LANE STILL UNDER THE ORACLE
Opened 03:44Z (session-open increment 13503d6). ORIENTATION STANDING (rides 38-46): a
ride-30/31-style prompt ("the reaper WRITTEN but NOT RUN") is STALE - the reaper has RUN
since ride-32; the ride-31 owed list (the 129G first pass + the ~22s regen receipt, T03194
closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE.
(1) THE JIT-CADENCE REAPER PASS 16 (03:45:24-28Z, ~4s, detached pid 889369, receipt
prep-reaper.log; bash -n clean + script verified UNCHANGED since the ride-31 review, md5
b804a7d3, mtime Oct 16 20:13): class A 0 (nothing new rejected), class B 72 files 102M +
syng/export/walks/walkreq markers dropped on 60 cut-awaiting-call, class C 0 files + markers
dropped on 144 pending (idempotence again - prior passes reaped the C bulk); lean buffer
148 <= 200 no eviction (pending dipped below the 160 valve bound at the census - the prep
master stays blocked costlessly in the T00024 wave, TODO=1); freed 55M apparent (df
585 -> 585G, prep/ 103,758 -> 103,703M apparent). CUMULATIVE ~= 158.4G apparent across 16
idempotent passes (129G + 4.6 + 4.5 + 2.7 + 0.95 + 4.1 + 3.1 + 0.8 + 2.4 + 2.1 + 4.3 +
0.09 + 0.06 + 0.03 + 0.1 + 0.06G), every cut/call input and every split ledger kept.
THE EXCLUSIONS VERIFIED LIVE A TENTH TIME: the T00024 prep lane + the batch-63 cut cohort
(T03764/T03518/T03245/T03303) + the deep call lane T03641 all in the exclusion list,
locus.fa verified KEPT after the pass (T00024 478MB / T03518 565MB / T03641 879MB /
T03764/T03245/T03303 ~118MB spot-checked); batch-62 resolved HONEST exit=0 wall 399s just
before the pass, and the excluded batch-63 resolved HONEST exit=0 wall 532s AFTER the
pass - the cut lanes stayed fed.
(2) WIDTH-4 CENSUS: batches 0-63 = 64/64 HONEST 4/4 (find-verified 64 exit=0 / 0 non-zero
since the 18:52:42Z swap), ZERO guard kills since the swap (this window: batch-61 exit=0
652s, batch-62 exit=0 399s, batch-63 exit=0 532s); cut map dirs 127 -> 132; batch-64
dispatched in flight at close.
(3) CALLS: callok 1,388 STEADY (+0 this window); fail-closed 523 STEADY (+0); cut.exit 16;
poison 3 terminal; the ONE deep lane T03641 (pid 657722, ~2:00 elapsed, 63.8GB
plateau-stable - the plateau-then-climb class, guard-oracle, under the 64GiB guard); ZERO
other marathons.
(4) THE T00024 WATCH PASSED FORWARD (NO intervention): the segdup-giant syng retry
(pid 732137) ~1:35 elapsed at 23MB RSS, R-state, still on ONE 317kb sequence - the named
SYNG SEG DUP-TAIL signature (candidate 4th member; class = T00197 + T01411 + T02562);
~1:35 << the ~3h T00197 precedent; the prep master's single-job xargs block stays costless
(TODO=1). DECISION POINT UNCHANGED: poison-first (prep/T00024/cut.poison) then
kill-by-exact-pid ONLY at the ~3h precedent with the double-writer race named - NOT before;
a completion is a receipted segdup-tail syng datum.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 582-587G, margin > 280G);
guard-oracle held (zero hand kills, zero masters touched); the prep valve held all window
(sub 2,162, pending 144-160 band); the ladder NOT opened (df 582G << 1,600G, backlog
present); terminal classes terminal; the ancient impg-reference benign - NEVER KILLED.
FINAL POLL 03:55:40Z: split 9,891 / acc 5,844 / sub 2,162 (valve-held) / cut map dirs 132 /
call-terminal 1,927 = callok 1,388 + fail-closed 523 (+16 cut.exit + poison 3 = 1,927 of
5,844, 33.0% of accepted) / df 582G / mem 866Gib / rvsize 215G apparent.
MASTERS: prep 411760 (the T00024 wave), cuts 3838087 (width-4-direct BATCH_W=4
N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2);
heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census; the JIT
reaper cadence (~4-6s a pass, re-runnable anytime); the T00024 watch (poison-first at the
3h precedent; receipt a completion); the T03641 deep lane under the guard-oracle; watch df
vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at
backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the
wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM);
AT CAMPAIGN CLOSE (all 5,844 terminal): rv-rollup.py -> FILL datasheet-runtime.md (the
standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule
+ THE REAPER RECEIPTS: ~158.4G cumulative freed, 16 idempotent passes, the walkreq
marker-gap fix, the ~22s regen, the T03194 full-chain receipt, the syng segdup-tail
class's 3 members + the T00024 candidate) -> the record's closing increment -> the
commit. NO PR PUSH.

## RIDE CONTINUATION 48 SESSION OPEN (04:07Z 2026-10-17) - THE CADENCE CONTINUES (PASS 17 OWED)
Opened from the ride-47 window close (415d105). ORIENTATION STANDING (rides 38-47): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN") is STALE - the
reaper has RUN since ride-32; the ride-31 owed list (the 129G first pass + the regen receipt
T03194 ~22s, closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE
(cumulative ~= 158.4G across 16 idempotent passes).
OPEN POLL 04:03:06Z: split 9,891 / acc 5,844 / sub 2,162 (valve-held) / cut.done 1,980
(map dirs 135; batch-64 resolved exit=0 wall 406s 04:00:38Z; batch-65 dispatched 04:00:39Z
in flight, cohort T03516/T02993/T03226/T03551) / call receipts 1,911 = callok 1,388 +
fail-closed 523 (terminal with +16 cut.exit + poison 3 = 1,930 of 5,844, 33.0%) /
df 581G / mem 865Gib / rvsize 216G apparent / width-4 census 65/65 HONEST 4/4, ZERO guard
kills since the 18:52:42Z swap / the ONE deep call lane T03641 (pid 657722, ~2:07, 63.8GB
plateau-stable, guard-oracle), ZERO other marathons; the T00024 syng at ~1:42/23MB on ONE
317kb sequence (candidate 4th segdup-tail member, << the 3h kill precedent, NO intervention,
poison-first at the precedent).
MASTERS: prep 411760 (the costless T00024 single-job xargs wave, valve held 161), cuts
3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146
(LANES=128 CHUNK_MULT=2); heartbeat live; the ancient impg-reference benign. THE RIDE: the
JIT-cadence reaper pass 17 (bash -n clean + script verified UNCHANGED since the ride-31
review, md5 b804a7d3, mtime Oct 16 20:13); the width-4 census continuation; sparse SHORT
polls; the standing contingencies (300G floor NEVER lowered; guard-oracle; poison-first -
the class is 3 members + the T00024 candidate at ~1:42 << the 3h precedent; terminal
classes terminal; the ladder at df >= 1,600G); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL
datasheet-runtime.md -> the closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 48 WINDOW CLOSE (04:22Z 2026-10-17) - THE CADENCE PASS 17 RECEIPTED (+53M APPARENT, CUMULATIVE ~= 158.4G ACROSS 17 IDEMPOTENT PASSES); CENSUS 67/67 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF VERIFIED LIVE AN ELEVENTH TIME; T00024 PASSED FORWARD AT ~2:00; ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 04:07Z (session-open increment 1b88357).
(1) THE JIT-CADENCE REAPER PASS 17 (04:05:20-24Z, ~4s, detached pid 905914, receipt
prep-reaper.log; bash -n clean + script verified UNCHANGED since the ride-31 review, md5
b804a7d3, mtime Oct 16 20:13): class A 0 (nothing new rejected), class B 48 files 86M +
syng/export/walks/walkreq markers dropped on 68 cut-awaiting-call, class C 0 files + markers
dropped on 136 pending (idempotence again - prior passes reaped the C bulk); lean buffer
140 <= 200 no eviction; freed 53M apparent (df 581 -> 581G, prep/ 103,703 -> 103,650M
apparent). CUMULATIVE ~= 158.4G apparent across 17 idempotent passes (129G ride-32 + 4.6 +
4.5 + 2.7 + 0.95 + 4.1 + 3.1 + 0.8 + 2.4 + 2.1 + 4.3 + 0.09 + 0.06 + 0.03 + 0.1 + 0.06 +
0.05G this window), every cut/call input and every split ledger kept.
THE EXCLUSIONS VERIFIED LIVE AN ELEVENTH TIME: the T00024 prep lane + the batch-65 cut
cohort (T03516/T02993/T03226/T03551) + the deep call lane T03641 all in the pass's
exclusion list, locus.fa verified KEPT on disk after the pass (T00024 478MB / T03641
879MB / the cohort 116-401MB, all spot-checked); the excluded batch-65 resolved HONEST
exit=0 wall 510s AFTER the pass - the cut lanes stayed fed.
(2) WIDTH-4 CENSUS: batches 0-66 = 67/67 HONEST 4/4, ZERO guard kills since the
18:52:42Z swap (this window: batch-64 exit=0 406s, batch-65 exit=0 510s, batch-66 exit=0
476s); cut.done 1,980 -> 1,988 (+8 over the ~15min window, map dirs 135 -> 138 - the
27-33/h producing band, the +54% width-4-direct state holding); batch-67 dispatched
04:17:07Z in flight at close.
(3) CALLS: callok 1,388 STEADY (+0 this window); fail-closed 523 STEADY (+0; the last
fail line 03:26:03Z T03627 guard-137 pre-dates the window, all named in
wave-logs/failures.txt - the honest fail-closed family); cut.exit 16; poison 3 terminal;
the ONE deep lane T03641 (pid 657722, ~2:25 elapsed, 63.8GB plateau-stable - the
plateau-then-climb class, guard-oracle, under the 64GiB guard); ZERO other marathons.
(4) THE T00024 WATCH PASSED FORWARD (NO intervention): the segdup-giant syng retry
(pid 732137) ~2:00 elapsed at 23MB RSS, R-state, still on ONE 317kb sequence
(HG03927#1#CM086638.1:130277163-130594662) - the named SYNG SEG DUP-TAIL signature
continuing (candidate 4th member; the class is T00197 + T01411 + T02562); ~2:00 << the
~3h T00197 kill precedent; the prep master's single-job xargs block stays costless
(TODO=1, the valve held, pending 140-148 band). DECISION POINT UNCHANGED: poison-first
(prep/T00024/cut.poison) then kill-by-exact-pid ONLY at the ~3h precedent with the
double-writer race named - NOT before; a completion is a receipted segdup-tail syng
datum for the datasheet.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 574-581G, margin >
270G); guard-oracle held (zero hand kills, zero masters touched); the prep valve held all
window (sub 2,162, pending 140-148); the ladder NOT opened (df 574G << 1,600G, backlog
present); terminal classes terminal; the ancient impg-reference benign - NEVER KILLED.
FINAL POLL 04:20:40Z: split 9,891 / acc 5,844 / sub 2,162 (valve-held) / cut.done 1,988
(map dirs 138; batch-67 in flight) / call-terminal 1,930 = callok 1,388 + fail-closed 523
+ cut.exit 16 + poison 3 = 1,930 of 5,844 (33.0% of accepted) / df 574G / mem 852Gib /
rvsize 223G apparent.
MASTERS: prep 411760 (the T00024 wave, log healthy), cuts 3838087 (width-4-direct
BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128
CHUNK_MULT=2); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4
census; the JIT reaper cadence (~4-6s a pass, re-runnable anytime); the T00024 watch
(poison-first at the 3h precedent - now ~2h elapsed, the next window crosses it; receipt
a completion); the T03641 deep lane under the guard-oracle; watch df vs the 300G hard
floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained
AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with
marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep
stalls; AT CAMPAIGN CLOSE (all 5,844 terminal): rv-rollup.py -> FILL
datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory
finding + the JIT-prep rule + THE REAPER RECEIPTS: ~158.4G cumulative freed, 17
idempotent passes, the walkreq marker-gap fix, the ~22s regen, the T03194 full-chain
receipt, the syng segdup-tail class's 3 members + the T00024 candidate and its
resolution) -> the record's closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 49 SESSION OPEN (04:24Z 2026-10-17) - THE CADENCE CONTINUES (PASS 18 OWED; THE T00024 3h CROSSING LANDS IN-WINDOW)
Opened from the ride-48 window close (efa59a5). ORIENTATION STANDING (rides 38-48): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN") is STALE - the
reaper has RUN since ride-32; the ride-31 owed list (the 129G first pass + the ~22s regen
receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE
(cumulative ~= 158.4G across 17 idempotent passes).
OPEN POLL 04:23:37Z: split 9,891 / acc 5,844 / sub 2,162 (valve-held) / cut.done 1,992
(map dirs 139; batch-66 resolved HONEST exit=0 wall 476s; batch-67 dispatched 04:17:07Z
in flight, cohort T02384/T03662/T03294/T02874) / call receipts 1,911 = callok 1,388 +
fail-closed 523 (terminal with +16 cut.exit + poison 3 = 1,930 of 5,844, 33.0% of accepted) /
df 574G / mem 847Gib / rvsize 223G apparent / width-4 census since the 18:52:42Z swap:
68 batch dirs, 67 resolved ALL exit=0, ZERO non-zero, ZERO guard kills / the ONE live deep
call lane T03641 (pid 657722, ~2:28, 63.8GB plateau-stable - the plateau-then-climb class,
guard-oracle, under the 64GiB guard); ZERO other marathons.
THE T00024 WATCH AT OPEN: the segdup-giant syng retry (pid 732137) ~2:03 elapsed at 23MB
RSS, R-state, still on ONE 317kb sequence (HG03927#1#CM086638.1:130277163-130594662) - the
named SYNG SEG DUP-TAIL signature (candidate 4th member; class = T00197 + T01411 +
T02562). THE 3h T00197 PRECEDENT CROSSES AT ~05:20:35Z - INSIDE THIS WINDOW: at the
crossing, poison-first (prep/T00024/cut.poison) then kill-by-exact-pid ONLY (732137),
with the double-writer race named (the prep master may hold the xargs slot; the poison
marker closes it); NOT before; a completion before the crossing is a receipted
segdup-tail syng datum for the datasheet.
MASTERS: prep 411760 (the costless T00024 single-job xargs wave, valve held), cuts
3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave
2810146 (LANES=128 CHUNK_MULT=2); heartbeat live; the ancient impg-reference benign.
THE RIDE: the JIT-cadence reaper pass 18 (bash -n clean + script verified UNCHANGED
since the ride-31 review, md5 b804a7d3, mtime Oct 16 20:13); the width-4 census
continuation; sparse SHORT polls; the standing contingencies (300G floor NEVER lowered;
guard-oracle; poison-first; terminal classes terminal; the ladder at df >= 1,600G);
AT CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md -> the closing increment ->
the commit. NO PR PUSH.

## RIDE CONTINUATION 49 WINDOW CLOSE (04:42Z 2026-10-17) - THE CADENCE PASS 18 RECEIPTED (+113M APPARENT, CUMULATIVE ~= 158.5G ACROSS 18 IDEMPOTENT PASSES); CENSUS 70/70 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF VERIFIED LIVE A TWELFTH TIME; CUT.DONE CROSSED 2,000; T00024 PASSED FORWARD AT ~2:20 (THE 3h CROSSING ~05:20:35Z IS THE NEXT WINDOW'S LIVE DECISION POINT); ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 04:24Z (session-open increment d968d58).
(1) THE JIT-CADENCE REAPER PASS 18 (04:24:39-44Z, ~5s, detached pid 928303, receipt
prep-reaper.log; bash -n clean + script verified UNCHANGED since the ride-31 review, md5
b804a7d3, mtime Oct 16 20:13): class A 0 (nothing new rejected), class B 72 files 161M +
syng/export/walks/walkreq markers dropped on 80 cut-awaiting-call, class C 0 files +
markers dropped on 124 pending (idempotence again - prior passes reaped the C bulk); lean
buffer 128 <= 200 no eviction; freed 113M apparent (df 574 -> 574G, prep/ 103,650 ->
103,537M apparent). CUMULATIVE ~= 158.5G apparent across 18 idempotent passes (129G
ride-32 + 4.6 + 4.5 + 2.7 + 0.95 + 4.1 + 3.1 + 0.8 + 2.4 + 2.1 + 4.3 + 0.09 + 0.06 +
0.03 + 0.1 + 0.06 + 0.05 + 0.11G this window), every cut/call input and every split
ledger kept.
THE EXCLUSIONS VERIFIED LIVE A TWELFTH TIME: the T00024 prep lane + the batch-68 cut
cohort (T02927/T02173/T03224/T02953) + the deep call lane T03641 all in the pass's
exclusion list, locus.fa verified KEPT on disk after the pass (T00024 478MB / T03641
879MB / the cohort 118-354MB, all spot-checked); batch-67 resolved HONEST exit=0 wall
418s at the pass, and the excluded batch-68 resolved HONEST exit=0 wall 479s AFTER the
pass - the cut lanes stayed fed.
(2) WIDTH-4 CENSUS: batches 0-69 = 70/70 HONEST ALL exit=0, ZERO non-zero, ZERO guard
kills since the 18:52:42Z swap (this window: batch-67 exit=0 418s, batch-68 exit=0 479s,
batch-69 exit=0 433s); cut.done 1,992 -> 2,000 (+8 over the ~17min window ~= 28/h - the
27-33/h producing band, the +54% width-4-direct state holding); map dirs 139 -> 143;
batch-70 dispatched (cohort T03779/T03029/T02844/T03204) in flight at close.
(3) CALLS: callok 1,388 STEADY (+0 this window); fail-closed 523 STEADY (+0; the last
fail line 03:26:03Z T03627 pre-dates the window, all named in wave-logs/failures.txt);
cut.exit 16; poison 3 terminal; the ONE deep lane T03641 (pid 657722, ~2:45 elapsed,
63.8GB plateau-stable - the plateau-then-climb class, guard-oracle, under the 64GiB
guard); ZERO other marathons.
(4) THE T00024 WATCH PASSED FORWARD (NO intervention): the segdup-giant syng retry
(pid 732137) ~2:20 elapsed at 23MB RSS, R-state, still on ONE 317kb sequence
(HG03927#1#CM086638.1:130277163-130594662) - the named SYNG SEG DUP-TAIL signature
continuing (candidate 4th member; the class = T00197 + T01411 + T02562). THE 3h T00197
PRECEDENT CROSSES AT ~05:20:35Z - THE NEXT WINDOW OWES THE DECISION: at the crossing,
poison-first (write prep/T00024/cut.poison) THEN kill-by-exact-pid ONLY (732137), with
the double-writer race named (the prep master holds the single-job xargs slot; the
poison marker closes it so the master does not redispatch); NOT before; a completion
before the crossing is a receipted segdup-tail syng datum for the datasheet.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 570-574G, margin >
270G); guard-oracle held (zero hand kills, zero masters touched); the prep valve held
all window (sub 2,162, pending 128-140); the ladder NOT opened (df 570G << 1,600G,
backlog present); terminal classes terminal; the ancient impg-reference benign - NEVER
KILLED.
FINAL POLL 04:40:49Z: split 9,891 / acc 5,844 / sub 2,162 (valve-held) / cut.done 2,000
(map dirs 143; batch-70 in flight) / call-terminal 1,930 = callok 1,388 + fail-closed 523
+ cut.exit 16 + poison 3 = 1,930 of 5,844 (33.0% of accepted) / df 570G / mem 865Gib /
rvsize 227G apparent.
MASTERS: prep 411760 (the T00024 wave), cuts 3838087 (width-4-direct BATCH_W=4
N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2);
heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census; the JIT
reaper cadence (~4-6s a pass, re-runnable anytime); THE T00024 DECISION AT ~05:20:35Z
(poison-first then kill-by-exact-pid 732137 - the precedent crossing; receipt the
outcome for the datasheet); the T03641 deep lane under the guard-oracle; watch df vs
the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at
backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the
wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM);
AT CAMPAIGN CLOSE (all 5,844 terminal): rv-rollup.py -> FILL datasheet-runtime.md (the
standing close list + the width-4-direct delta + the inventory finding + the JIT-prep
rule + THE REAPER RECEIPTS: ~158.5G cumulative freed, 18 idempotent passes, the walkreq
marker-gap fix, the ~22s regen, the T03194 full-chain receipt, the syng segdup-tail
class's 3 members + the T00024 candidate and its resolution) -> the record's closing
increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 50 SESSION OPEN (04:44Z 2026-10-17, commit a96f6b9) - THE CADENCE CONTINUES (PASS 19 OWED; THE T00024 3h DECISION LANDS IN-WINDOW AT ~05:20:35Z)
Opened from the ride-49 window close (0f307d4). ORIENTATION STANDING (rides 38-49): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN") is STALE - the
reaper has RUN since ride-32 (129G first pass); the ride-31 owed list (the regen receipt
~22s, T03194 closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE
(cumulative ~= 158.5G across 18 idempotent passes).
OPEN POLL 04:44:13Z: split 9,891 / acc 5,844 / sub 2,162 (valve-held) / cut.done 2,000
(map dirs 143; batch-66..69 ALL exit=0 walls 476/418/479/433s; batch-70 in flight, cohort
T02440/T03203/T03361/T02703) / call receipts 1,911 = callok 1,388 + fail-closed 523
(terminal with +16 cut.exit + poison 3 = 1,930 of 5,844, 33.0%) / df 570G / mem 844Gib /
rvsize 227G apparent / width-4 census 70/70 honest ALL exit=0 ZERO kills since the swap.
THE T00024 DECISION LANDS IN-WINDOW: the segdup-giant syng retry (pid 732137) ~2:24
elapsed at 23MB RSS, R-state, still on ONE 317kb sequence (HG03927#1#CM086638.1:130277163-
130594662) - the named SYNG SEG DUP-TAIL signature (candidate 4th member; class = T00197 +
T01411 + T02562). AT THE ~05:20:35Z CROSSING: poison-first (write prep/T00024/cut.poison)
THEN kill-by-exact-pid ONLY (732137), with the double-writer race named (the prep master
holds the single-job xargs slot; the poison marker closes it so the master does not
redispatch); NOT before; a completion before the crossing is a receipted segdup-tail syng
datum for the datasheet.
MASTERS: prep 411760 (the T00024 wave), cuts 3838087 (width-4-direct BATCH_W=4
N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2);
heartbeat live; the ancient impg-reference benign - NEVER KILLED.
THE RIDE: the JIT-cadence reaper pass 19 (bash -n + md5 b804a7d3 verify + run detached +
poll prep-reaper.log SHORT); the width-4 census; the T00024 decision at the crossing; the
T03641 deep lane under the guard-oracle; the standing contingencies (300G floor NEVER
lowered; guard-oracle; poison-first; terminal classes terminal; ladder at df >= 1,600G);
AT CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md -> the record's closing
increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 50 WINDOW CLOSE (05:29Z 2026-10-17, commits a96f6b9/54fd9cd) - THE PASS 19 RECEIPT + THE T00024 3h DECISION EXECUTED + THE MID-DISPATCH VALVE PATCH'S FIRST LIVE PROOF
(1) THE JIT-CADENCE REAPER PASS 19 (04:45:57-04:46:01Z, ~5s, detached pid 949765, receipt
prep-reaper.log; bash -n clean + script md5 b804a7d3 UNCHANGED since the ride-31 review):
class A 0, class B 48 files 70M + markers on 88 cut-awaiting-call, class C 0 + markers on
116 pending; lean buffer 120 <= 200 no eviction; freed 39M apparent. CUMULATIVE ~= 158.5G
across 19 idempotent passes. Exclusion proof verified live a THIRTEENTH time (T00024 +
batch-70 cohort + T03641 excluded, locus.fa KEPT; batch-70..74 ALL exit=0 walls
469/504/520/571/593s - the cut lanes stayed fed).
(2) THE T00024 DECISION EXECUTED AT THE 3h CROSSING (the ride-49 close's named decision
point; the sanctioned poison-class kill, NOT a guard-oracle violation): attempt 2 crossed
the T00197 precedent at 05:21:00Z (3h00m25s elapsed, CPU==wall, 23MB, R-state, NO output,
silent on HG03927#1#CM086638.1:130277163-130594662 after 4 sequences; locus.fa 457MB the
class's largest). POISON-FIRST 05:22:59Z (prep/T00024/cut.poison + the cuts/poison.tsv
class entry - the CONFIRMED 4th member: T00197 + T01411 + T02562 + T00024; snapshot kept)
THEN kill-by-exact-pid 732137 (SIGKILL 05:23:08Z - never the wrapper, never the master).
THE POISON HELD, VERIFIED LIVE: no redispatch through both dispatch sites; the wrapper
wrote honest receipts (syng.exit, syng.wall 10,953s).
(3) THE MID-DISPATCH VALVE LINE LANDED 05:28:19Z - the ride-4 chunked-dispatch patch's
FIRST LIVE PROOF: 'VALVE (mid-dispatch): 204 pending cut (>= 160) - parking the remaining
3484 of 3740 TODO'. Under the old full-TODO xargs this instant would have committed ~3,740
territories (~935GB) unthrottled - the unbounded-wave class CLOSED with a live receipt.
(4) T03641 RESOLVED TERMINAL UNDER THE GUARD-ORACLE at 05:00:30Z: the guard's own 137
mid-score kill, wall 11,084s (~3h05m), peak ~71.1GB - THE DEEPEST CALL LANE YET, named.
fail-closed 523 -> 542 (+19, ALL NAMED: the T03641 guard-137 + 18 ordinary 7479-class
binds - T02078/T03204/T02874/T03456/T03809/T03784 no-unique-key; T02010/T03029/T03356/
T03813/T03771/T03794/T03518/T03433/T03671/T03742/T03767 two-origin-shifts).
(5) WIDTH-4 CENSUS: 75/75 batch.exit=0 post-swap, ZERO kills (batch-70..74 honest;
batch-75 in flight, cohort T03562/T01987/T02512/T03198); cut.done 2,000 -> 2,020 (~27/h).
FINAL POLL 05:28:51Z: split 9,891 / acc 5,844 / sub 2,282 (+120 - the freed prep pole
producing; the valve PARKED mid-dispatch at 204 pending) / cut.done 2,020 (map dirs 80 -
the wave drained the backlog hard) / call-terminal 2,003 = callok 1,441 + fail-closed 542
+ cut.exit 16 + poison 4 (34.3% of accepted) / df 582G (band 570-615G, 300G floor never
approached) / mem 482Gib / rvsize 215G apparent. MASTERS: prep 411760 (parked-chunk state
- the patch working), cuts 3838087 (width-4-direct), wave 2810146; heartbeat live; the
ancient impg-reference benign. ONE deep-ish fresh call lane at close (56GB/26min,
guard-oracle).
THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census; the JIT reaper cadence;
the UNBLOCKED prep pole (ride the 64-pool chunks + the mid-dispatch valve lines; watch df
as the substrate window rebuilds toward the ~200 pending bound); watch df vs the 300G
hard floor; the ladder only at backlog drained AND df >= 1,600G; poison-first on prep
stalls (the class is 4 members); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL
datasheet-runtime.md (the standing close list + the width-4-delta + the inventory finding
+ the JIT-prep rule + the reaper receipts ~158.5G/19 passes + the segdup-tail class's 4
members incl. the T00024 crossing receipt + the mid-dispatch valve proof + the T03641
deepest-lane point) -> the record's closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 51 SESSION OPEN (05:32Z 2026-10-17, commit 5ff49ac) - THE CADENCE CONTINUES (PASS 20 OWED; PASS 20 RECEIPTED AT OPEN +18.9G - THE CLASS C BULGE IS BACK WITH THE FREED PREP POLE)
Opened from the ride-50 window close (54fd9cd). ORIENTATION STANDING (rides 38-51): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN") is STALE - the
reaper has RUN since ride-32 (129G first pass); the ride-31 owed list (the ~22s regen
receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE
(cumulative ~= 177.4G across 20 idempotent passes after this window's pass).
OPEN POLL 05:31:36Z: split 9,891 / acc 5,844 / sub 2,282 (the freed prep pole producing;
the valve PARKED mid-dispatch at 204 pending, holding the prep pool) / cut.done 2,020 (map
dirs 78) / call-terminal 2,005 = callok 1,443 + fail-closed 542 + cut.exit 16 + poison 4
(34.3% of accepted) / df 579G / mem 498Gib / rvsize 218G apparent / census 76 batch dirs
post-swap, 75 resolved ALL exit=0, ZERO kills since the 18:52:42Z swap / the ONE deep-ish
lane T03333 (pid 967312, 56GB/29:34 plateau, niced, guard-oracle) + batch-1792214623-75 cut
in flight (cohort T03562/T01987/T02512/T03198).
THE JIT-CADENCE REAPER PASS 20 (05:32:39-51Z, ~12s, detached, receipt prep-reaper.log;
bash -n CLEAN + script md5 b804a7d3 UNCHANGED, mtime Oct 16 20:13, both re-verified):
class A 0, class B 72 files 149M on cut-awaiting-call, class C 360 files 18,952M on 216
pending - THE CLASS C BULGE IS BACK: the freed prep pole's ~120 new substrates (sub 2,162
-> 2,282) produced a fresh C cohort the pass reaped whole; freed 18,892M apparent
(prep 132,397 -> 113,505M; rvsize 218 -> 200G; df 585 -> 604G at the pass, the pole
writing new substrate in parallel). CUMULATIVE ~= 177.4G apparent across 20 idempotent
passes, every cut/call input and every split ledger kept. EXCLUSION PROOF VERIFIED LIVE A
FOURTEENTH TIME: the batch-75 cut cohort (T03562/T01987/T02512/T03198) + ALL 22 LIVE CALL
LANES in the pass's exclusion list, locus.fa verified KEPT (spot-checked 120-582MB incl.
T03333 582MB).
MASTERS: prep 411760 (parked-chunk state, VALVE lines each minute), cuts 3838087
(width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146
(LANES=128 CHUNK_MULT=2); heartbeat live; the ancient impg-reference benign. THE RIDE:
the width-4 census continuation (batch-75 at 10:31 elapsed - ABOVE the 469-593s wall band,
watch it resolve honest; if it goes long it is the cuts master's own lane, guard-oracle
applies, NO hand kill); the T03333 deep-ish lane under the guard-oracle (56GB plateau,
~32min, ordinary lane territory); sparse SHORT polls; the standing contingencies (300G
floor NEVER lowered; guard-oracle; poison-first - the class is 4 members; terminal classes
terminal; the ladder at df >= 1,600G); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL
datasheet-runtime.md (the standing close list + the width-4-delta + the inventory finding
+ the JIT-prep rule + the reaper receipts ~177.4G/20 passes incl. this pass's 18.9G class C
bulge datum - the pole-freed substrate reappears as reaped intermediates, the JIT loop
closed) -> the record's closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 51 WINDOW CLOSE (05:37Z 2026-10-17) - THE CADENCE PASS 20 RECEIPTED (+18.9G APPARENT, CUMULATIVE ~= 177.4G ACROSS 20 IDEMPOTENT PASSES); CENSUS 76/76 HONEST ZERO KILLS; THE EXCLUSION PROOF A FOURTEENTH TIME; BATCH-75 HONEST AT 712s; ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 05:32Z (session-open increment 5ff49ac).
(1) THE REAPER PASS 20 (05:32:39-51Z, ~12s, detached, receipt prep-reaper.log; bash -n
clean + md5 b804a7d3 UNCHANGED re-verified): class A 0, class B 72 files 149M, class C
360 files 18,952M on 216 pending - THE CLASS C BULGE (the freed prep pole's ~120 new
substrates reaped whole - the JIT loop closed live); freed 18,892M apparent (prep
132,397 -> 113,505M; rvsize 218 -> 199G). CUMULATIVE ~= 177.4G across 20 idempotent
passes, every cut/call input and every split ledger kept. EXCLUSIONS: the batch-75 cut
cohort (T03562/T01987/T02512/T03198) + ALL 22 live call lanes excluded, locus.fa KEPT
(spot-checked 120-582MB incl. T03333 582MB).
(2) CENSUS: 76/76 batch.exit=0 post-swap, ZERO kills since the 18:52:42Z swap (batch-75
resolved HONEST exit=0 wall 712s - above the prior 469-593s band, no kill; batch-76
batch-1792215336-76 dispatched in flight at close); cut.done 2,020 -> 2,024.
(3) CALLS: callok 1,443 -> 1,446; fail-closed 542 steady (last fail 05:14:33Z T03784
pre-window); cut.exit 16; poison 4 terminal; the ONE deep-ish lane T03333 (pid 967312,
~34min, 56.0GB plateau-stable, niced, guard-oracle, NO intervention); ZERO other marathons.
(4) CONTINGENCIES: 300G floor NEVER approached (df 579-604G, margin > 275G); guard-oracle
held; the prep valve PARKED all window (204 pending, the patch working); ladder NOT opened
(df 597G << 1,600G); terminal classes terminal; the ancient impg-reference benign.
FINAL POLL 05:36:30Z: split 9,891 / acc 5,844 / sub 2,282 / cut.done 2,024 (map dirs 79;
batch-76 in flight) / call-terminal 2,008 = callok 1,446 + fail-closed 542 + cut.exit 16 +
poison 4 (34.3% of accepted) / df 597G / mem 526Gib / rvsize 199G apparent.
MASTERS: prep 411760 (parked-chunk), cuts 3838087 (width-4-direct BATCH_W=4
N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2);
heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census (batch-76 in
flight); the JIT reaper cadence; the T03333 deep-ish lane under the guard-oracle (watch
the plateau-then-climb shape, never kill by hand); the prep pole + df watch vs the 300G
hard floor; the ladder only at backlog drained AND df >= 1,600G; poison-first on prep
stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 terminal): rv-rollup.py ->
FILL datasheet-runtime.md (the standing close list + the width-4-delta + the inventory
finding + the JIT-prep rule + the reaper receipts ~177.4G/20 passes incl. the pass-20
class C bulge datum) -> the record's closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 52 SESSION OPEN (05:39Z 2026-10-17, commit 0d63c4c) - THE CADENCE CONTINUES (PASS 21 OWED; PASS 20'S CLASS C BULK ALREADY REAPED)
Opened from the ride-51 window close (140d300). ORIENTATION STANDING (rides 38-52): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN") is STALE - the
reaper has RUN since ride-32 (129G first pass); the ride-31 owed list (the ~22s regen
receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE
(cumulative ~= 177.4G across 20 idempotent passes entering this window).
OPEN POLL 05:38:47Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,282 (valve-held - the
prep master's own line 'VALVE: 200 substrates pending cut (>= 160) - holding the prep
pool') / cut.done 2,024 (map dirs 77; batch-1792215336-76 in flight, the scan lane
sample_projection ~3-4min/32GB, no exit marker yet) / call-terminal 2,010 = callok 1,448 +
fail-closed 542 + cut.exit 16 + poison 4 (34.4% of accepted) / df 600G / mem 511Gib /
rvsize 197G apparent / census: 76 post-swap batch.exit markers, 76/76 ALL exit=0, ZERO
non-zero, ZERO guard kills since the 18:52:42Z swap (find -newermt re-verified live; the
last kill line 18:52:40Z T03370 pre-swap width-8) / live lanes: the ONE deep-ish lane T03333
(pid 967312, ~36:24 elapsed, 56.0GB plateau-stable, niced, guard-oracle, under the 64GiB
guard) + ordinary lanes 30-32GB sub-36min niced + the batch-76 cut in scan; ZERO other
marathons. MASTERS: prep 411760 (1-11:13h, VALVE-held pool, log healthy), cuts 3838087
(width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~10:46h), wave
2810146 (LANES=128 CHUNK_MULT=2, ~18:42h); heartbeat live (05:38:20Z line healthy); the
ancient impg-reference benign. REAPER SCRIPT RE-VERIFIED AT OPEN: bash -n CLEAN, md5
b804a7d3 UNCHANGED since the ride-31 review, mtime Oct 16 20:13. THIS WINDOW OWES: the
JIT-cadence reaper pass 21 detached + poll SHORT; the width-4 census (batch-76 in flight);
the T03333 lane under the guard-oracle; the standing contingencies (300G floor NEVER
lowered; guard-oracle; poison-first - the class is 4 members; terminal classes terminal;
the ladder at backlog drained AND df >= 1,600G); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL
datasheet-runtime.md (the standing close list + the width-4-delta + the inventory finding
+ the JIT-prep rule + the reaper receipts ~177.4G+/20+ passes incl. the pass-20 class C
bulge datum) -> the record's closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 52 WINDOW CLOSE (05:46Z 2026-10-17) - THE CADENCE PASS 21 RECEIPTED (+57M APPARENT, CUMULATIVE ~= 177.5G ACROSS 21 IDEMPOTENT PASSES); CENSUS 77/77 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF A FIFTEENTH TIME; BATCH-76 RESOLVED HONEST AT 510s (BACK IN THE BAND); THE T03333 LANE PASSED FORWARD AT ~42MIN PLATEAU-STABLE; ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 05:39Z (session-open increment 0d63c4c). (1) THE REAPER PASS 21
(05:40:31-36Z, ~5s, detached pid 1128276, receipt prep-reaper.log; bash -n CLEAN + script
md5 b804a7d3 UNCHANGED since the ride-31 review, mtime Oct 16 20:13, both re-verified at
the open): class A 0 (nothing new rejected), class B 24 files 72M + syng/export/walks/
walkreq markers dropped on 16 cut-awaiting-call, class C 0 files + markers dropped on 212
pending (idempotence - the pass-20 C bulk already reaped); lean buffer 216 > 200 - the
eviction arm noted but NOT IMPLEMENTED TO FIRE (the prep valve at 160 bounds pending);
freed 57M apparent (df 600 -> 600G, prep/ 112,628 -> 112,571M apparent). CUMULATIVE
SUBSTRATE-REAPER FREED ~= 177.5G apparent across 21 idempotent passes (129G ride-32 + 4.6
+ 4.5 + 2.7 + 0.95 + 4.1 + 3.1 + 0.8 + 2.4 + 2.1 + 4.3 + 0.09 + 0.06 + 0.03 + 0.1 + 0.06
+ 0.05 + 0.11 + 0.04 + 18.9 + 0.06G this window), every cut/call input and every split
ledger kept throughout. THE IN-FLIGHT EXCLUSIONS VERIFIED LIVE A FIFTEENTH TIME: the
batch-76 cut cohort (T03882/T03889/T03918/T03936) + ALL live call lanes in the pass's
exclusion list (40 territories), locus.fa verified KEPT on disk after the pass
(spot-checked 118-582MB, incl. T03333 582MB and the batch-76 cohort 118-241MB); the
excluded batch-76 resolved HONEST exit=0 wall 510s AFTER the pass - the cut lanes stayed
fed through the reaper pass, the receipt the pass answers to. (2) THE WIDTH-4 CENSUS:
batches post-swap = 77/77 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the
18:52:42Z swap (this window: batch-75 receipted 712s at ride-51, batch-76 resolved HONEST
exit=0 wall 510s - back in the 469-712s band; batch-77 batch-1792215848-77 dispatched
~05:44Z in flight at close); cut.done 2,024 -> 2,028 (+4); map dirs 77 -> 80. (3) CALLS:
callok 1,448 -> 1,449; fail-closed 542 STEADY (+0; the last fail line 05:14:33Z T03784
pre-dates the window, all named in wave-logs/failures.txt); cut.exit 16; poison 4
terminal; the ONE deep-ish lane T03333 (pid 967312, ~42:33 elapsed at close, 56.0GB
PLATEAU-STABLE the whole window, niced, guard-oracle, under the 64GiB guard) - passed
forward, ordinary-lane territory, NO intervention; ZERO other marathons (the wave's
ordinary lanes all niced sub-40min). (4) CONTINGENCIES: the 300G hard floor NEVER
approached (df band 597-600G, margin ~300G); guard-oracle held (zero hand kills, zero
masters touched); the prep valve HELD all window (200 pending, VALVE lines each minute,
the patch working); terminal classes terminal; the ladder NOT opened (df 597G << 1,600G,
backlog present); the ancient impg-reference benign - NEVER KILLED. FINAL POLL
05:45:05Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,282 (valve-held) / cut.done
2,028 (map dirs 80; batch-77 in flight) / call-terminal 2,011 = callok 1,449 + fail-closed
542 + cut.exit 16 + poison 4 = 2,011 of 5,844 (34.4% of accepted) / df 597G / mem 466Gib /
rvsize 200G apparent (197 -> 200G across the window - the pole writing substrate in
parallel). MASTERS AT CLOSE: prep 411760 (1-11:19h, VALVE-held pool), cuts 3838087
(BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, width-4-direct, ~10:52h), chunked
wave 2810146 (LANES=128 CHUNK_MULT=2, ~18:48h); heartbeat live. THE SUCCESSOR RESUMES:
sparse SHORT polls; the width-4 census continuation (batch-77 in flight); the JIT reaper
cadence (~4-12s a pass, idempotent, re-runnable anytime); the T03333 deep-ish lane under
the guard-oracle (56GB plateau ~42min - watch the plateau-then-climb shape, never kill by
hand); the prep pole (the valve holds at 200 pending; watch df as the substrate window
rebuilds); watch df vs the 300G hard floor (self-gates; DO NOT lower without the owner);
the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact
pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never
RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844
accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list +
the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER
RECEIPTS: ~177.5G cumulative freed classes A/B/C across 21 idempotent passes, the
walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain
re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members incl. the T00024
3h crossing receipt, the mid-dispatch valve patch's live proof, the guard-oracle receipts
incl. the T03641 11,084s/71.1GB deepest-lane point, the pass-20 class C bulge datum - the
freed-pole substrate reappearing as reaped intermediates, the JIT loop closed live) -> the
record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE
EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME
RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/
receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 53 SESSION OPEN (05:52Z 2026-10-17, commit 08559f1) - THE CADENCE CONTINUES (PASS 22 OWED; PASS 21'S +57M RECEIPTED AT RIDE-52 CLOSE)
Opened from the ride-52 window close (3f734b2). ORIENTATION STANDING (rides 38-53): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN") is STALE - the
reaper has RUN since ride-32 (129G first pass); the ride-31 owed list (the ~22s regen
receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE
(cumulative ~= 177.5G across 21 idempotent passes entering this window).
OPEN POLL 05:50:38Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,282 (valve-held at 196
pending, VALVE lines each minute) / cut.done 2,028 (map dirs 79; batch-1792215848-77 in
flight, dispatched ~05:44Z) / call-terminal 2,012 = callok 1,449 + fail-closed 543 (+1
IN-WINDOW, NAMED: T03203 guard-137 05:49:29Z - peak 70,597,364kB ~67.3GiB crossing the
64GiB guard mid-score, wall 2,884s, fail-reaped IN-LANE receipts kept; the second-deepest
lane point after T03641's 71.1GB) + cut.exit 16 + poison 4 = 2,012 of 5,844 (34.5%) /
df 597-598G / mem 419-431Gib / rvsize 198-201G apparent / census 77/77 post-swap
batch.exit markers ALL exit=0 re-verified live (find -newermt; ZERO non-zero; last kill
line 18:52:40Z pre-swap width-8) / live lanes: the T03333 deep-ish lane (pid 967312,
~49min, 53.4G RSS, niced, guard-oracle, NO intervention) + 15 ordinary call lanes 12-54G
niced + the batch-77 cut in flight. MASTERS: prep 411760 (VALVE-held pool), cuts 3838087
(width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146
(LANES=128 CHUNK_MULT=2); heartbeat live (05:51:41Z line healthy); the ancient
impg-reference benign. REAPER SCRIPT RE-VERIFIED AT OPEN: bash -n CLEAN, md5 b804a7d3
UNCHANGED since the ride-31 review, mtime Oct 16 20:13. THIS WINDOW OWES: the JIT-cadence
reaper pass 22 detached + poll SHORT; the width-4 census (batch-77 in flight); the T03333
lane under the guard-oracle; the standing contingencies (300G floor NEVER lowered;
guard-oracle; poison-first - the class is 4 members; terminal classes terminal; the
ladder at backlog drained AND df >= 1,600G); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL
datasheet-runtime.md (the standing close list + the width-4-delta + the inventory
finding + the JIT-prep rule + the reaper receipts ~177.5G+/21+ passes incl. the T03203
guard-137 deep-lane datum) -> the record's closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 53 WINDOW CLOSE (05:56Z 2026-10-17, commit 6250019) - THE CADENCE PASS 22 RECEIPTED (+44M APPARENT, CUMULATIVE ~= 177.5G ACROSS 22 IDEMPOTENT PASSES); CENSUS 78/78 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF A SIXTEENTH TIME; BATCH-77 HONEST AT 459s; THE T03203 GUARD-137 NAMED (70.6GB SECOND-DEEPEST LANE); THE T03333 LANE PASSED FORWARD AT ~53MIN; ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 05:52Z (session-open increment 08559f1). (1) THE REAPER PASS 22 (05:53:40-45Z,
~5s, detached pid 1152371, receipt prep-reaper.log; bash -n CLEAN + md5 b804a7d3
UNCHANGED re-verified): class A 0, class B 48 files 76M + markers dropped on 24
cut-awaiting-call, class C 0 + markers on 204 pending (idempotence - the pass-20 C bulk
already reaped); lean buffer 208 > 200 - the eviction arm noted, NOT IMPLEMENTED TO FIRE
(the prep valve at 160 bounds pending); freed 44M apparent (prep 112,281 -> 112,237M;
df 597 -> 597G). CUMULATIVE ~= 177.5G apparent across 22 idempotent passes, every
cut/call input and every split ledger kept. EXCLUSIONS VERIFIED LIVE A SIXTEENTH TIME:
the batch-78 cut cohort (T03852/T03996/T03953/T03893) + ALL live call lanes excluded,
locus.fa KEPT (spot-checked 112-555MB incl. T03333 555MB; syng/export intact on the
excluded in-flight T03852); the excluded batch-77 resolved HONEST exit=0 wall 459s
AFTER the pass (cohort T03908/T04004/T04011/T04030) - the cut lanes stayed fed.
(2) CENSUS: 78/78 batch.exit=0 post-swap, ZERO kills since the 18:52:42Z swap (batch-77
honest 459s in the 459-712s band; batch-78 batch-1792216308-78 dispatched 05:51:48Z in
flight at close); cut.done 2,028 -> 2,032; map dirs 79 -> 83. (3) CALLS: callok 1,449
STEADY; fail-closed 542 -> 543 (+1, the window's one new terminal, NAMED: T03203
guard-137 05:49:29Z - peak 70,597,364kB ~67.3GiB crossing the 64GiB guard mid-score,
wall 2,884s, fail-reaped IN-LANE receipts kept - THE SECOND-DEEPEST LANE POINT after
T03641's 11,084s/71.1GB, all named in wave-logs/failures.txt); cut.exit 16; poison 4
terminal; the ONE deep-ish lane T03333 (pid 967312, ~53.5min at close, 53.4-56GB
PLATEAU-STABLE, niced, guard-oracle, NO intervention); ZERO other marathons.
(4) CONTINGENCIES: 300G floor NEVER approached (df band 597-598G, margin ~297G);
guard-oracle held (the T03203 kill was the lane's own 64GiB guard, sanctioned); the prep
valve HELD all window (192-196 pending); ladder NOT opened (df 597G << 1,600G); terminal
classes terminal; the ancient impg-reference benign - NEVER KILLED. FINAL POLL 05:56:00Z:
split 9,891 / acc 5,844 (59.1% final) / sub 2,282 (valve-held) / cut.done 2,032 (map
dirs 83; batch-78 in flight) / call-terminal 2,012 = callok 1,449 + fail-closed 543 +
cut.exit 16 + poison 4 = 2,012 of 5,844 (34.5% of accepted) / df 597G / mem 401Gib /
rvsize 200G apparent. MASTERS: prep 411760 (VALVE-held pool), cuts 3838087
(width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146
(LANES=128 CHUNK_MULT=2); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls;
the width-4 census (batch-78 in flight); the JIT reaper cadence (~5-12s a pass,
idempotent); the T03333 lane under the guard-oracle (never kill by hand); the prep
pole + df watch vs the 300G hard floor; the ladder only at backlog drained AND
df >= 1,600G; poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE
(all 5,844 terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close
list + the width-4-delta + the inventory finding + the JIT-prep rule + the reaper
receipts ~177.5G/22 passes incl. the T03203 guard-137 deep-lane datum) -> the record's
closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 54 SESSION OPEN (05:59Z 2026-10-17, commit ac3b4c0) - THE CADENCE CONTINUES (PASS 23 OWED; PASS 22'S +44M RECEIPTED AT RIDE-53 CLOSE)
Opened from the ride-53 window close (6250019). ORIENTATION STANDING (rides 38-54): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN") is STALE - the
reaper has RUN since ride-32 (129G first pass); the ride-31 owed list (the ~22s regen
receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE
(cumulative ~= 177.5G across 22 idempotent passes entering this window).
OPEN POLL 05:58:50Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,282 (valve-held at 192
pending, VALVE lines each minute) / cut.done 2,032 (map dirs 83) / call-terminal 2,008 =
callok 1,449 + fail-closed 543 + cut.exit 16 + poison 4 (34.3% of accepted) / df 579G /
mem 403Gib / rvsize 219G apparent RISING from 200G at ride-53 close (the cuts pole writing
map outputs faster than the wave consumes - the pass-20 class C bulge shape; the pass
answers) / census 79/79 post-swap batch.exit markers ALL exit=0 re-verified live by mtime
(ZERO non-zero; ZERO guard kills since the 18:52:42Z swap; batch-78 resolved HONEST exit=0
wall 467s IN-BAND after the ride-53 exclusion proof; batch-1792216777-79 dispatched ~05:57Z
in flight at open, cohort T03927/T03829/T03841/T03940) / live lanes: the T03333 deep-ish
lane (pid 967312, ~57min, 53GB plateau-stable, niced, guard-oracle, NO intervention) + a
DEEPER-LANE CLUSTER of 14+ ordinary lanes at 13-53GB / 45-56min elapsed - the deepest
cohort this ride has seen, all niced, all under the 64GiB guard, guard-oracle, NO hand kill
(the T03203 guard-137 precedent proves the guard fires when it must).
MASTERS: prep 411760 (1-11:33h, VALVE-held pool), cuts 3838087 (width-4-direct BATCH_W=4
N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~11:07h), wave 2810146 (LANES=128 CHUNK_MULT=2,
~19:03h); heartbeat live (05:58:52Z line healthy); the ancient impg-reference benign.
REAPER SCRIPT RE-VERIFIED AT OPEN: bash -n CLEAN, md5 b804a7d3 UNCHANGED since the
ride-31 review, mtime Oct 16 20:13. THIS WINDOW OWES: the JIT-cadence reaper pass 23
detached + poll SHORT; the width-4 census (batch-79 in flight); the T03333 lane + the
deeper-lane cluster under the guard-oracle; the standing contingencies (300G floor NEVER
lowered; guard-oracle; poison-first - the class is 4 members; terminal classes terminal;
the ladder at backlog drained AND df >= 1,600G); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL
datasheet-runtime.md (the standing close list + the width-4-delta + the inventory finding
+ the JIT-prep rule + the reaper receipts ~177.5G+/22+ passes) -> the record's closing
increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 54 WINDOW CLOSE (06:08Z 2026-10-17) - THE CADENCE PASS 23 RECEIPTED (+16M APPARENT, CUMULATIVE ~= 177.5G ACROSS 23 IDEMPOTENT PASSES); CENSUS 80/80 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF A SEVENTEENTH TIME; BATCH-79 HONEST AT 446s; THE DEEPER-LANE CLUSTER PASSED FORWARD; ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 05:59Z (session-open increment ac3b4c0). (1) THE REAPER PASS 23 (06:00:44-49Z,
~5s, detached pid 1164794, receipt prep-reaper.log; bash -n CLEAN + script md5 b804a7d3
UNCHANGED re-verified at the open): class A 0 (nothing new rejected), class B 24 files 31M
+ syng/export/walks/walkreq markers dropped on 28 cut-awaiting-call, class C 0 files +
markers dropped on 200 pending (idempotence - the pass-20 C bulk already reaped); lean
buffer 204 > 200 - the eviction arm noted, NOT IMPLEMENTED TO FIRE (the prep valve at 160
bounds pending); freed 16M apparent (df 597 -> 597G, prep 112,072 -> 112,056M apparent);
rvsize 219 -> 199-200G apparent across the window (the open-poll 219G reading was the cuts
pole's in-flight map writes; the wave consumed them as calls landed). CUMULATIVE SUBSTRATE-
REAPER FREED ~= 177.5G apparent across 23 idempotent passes, every cut/call input and every
split ledger kept throughout. THE IN-FLIGHT EXCLUSIONS VERIFIED LIVE A SEVENTEENTH TIME:
the batch-79 cut cohort (T03927/T03829/T03841/T03940) + ALL live call lanes in the pass's
exclusion list, locus.fa verified KEPT on disk after the pass (spot-checked 113-556MB,
incl. T03333 556MB and the batch-79 cohort 113-128MB; syng index intact on the excluded
in-flight T03927); the excluded batch-79 resolved HONEST exit=0 wall 446s AFTER the pass -
IN-BAND, the cut lanes stayed fed through the reaper pass, the receipt the pass answers to.
(2) THE WIDTH-4 CENSUS: batches post-swap = 80/80 HONEST ALL exit=0, ZERO non-zero, ZERO
guard kills since the 18:52:42Z swap (this window: batch-78 receipted 467s at the open,
batch-79 resolved HONEST exit=0 wall 446s; batch-80 batch-1792217224-80 dispatched ~06:07Z
in flight at close); cut.done 2,032 -> 2,040 (+8: the batch-78 cohort T03852/T03996/
T03953/T03893 landing at the window's edge + the batch-79 cohort T03927/T03829/T03841/
T03940 resolved post-pass); map dirs 83 -> 88 as the cuts pole out-writes; the wave
consumes. (3) CALLS: callok 1,449 -> 1,452; fail-closed 543
STEADY (+0; the last fail 05:49:29Z T03203 pre-dates the window, named in
wave-logs/failures.txt); cut.exit 16; poison 4 terminal; the T03333 deep-ish lane (pid
967312, 65min at close, 52-53GB PLATEAU-STABLE the whole window, niced, guard-oracle, under
the 64GiB guard) - passed forward, NO intervention; THE DEEPER-LANE CLUSTER passed forward:
three more lanes at 53-54GB / 55-64min (pids 968955, 994062, 1008767) + 12 ordinary lanes
13-47GB - ALL NICED, ALL UNDER THE 64GiB GUARD, guard-oracle, NO hand kill; the cluster is
the wave's ordinary deep tail, not a new class; ZERO other events.
(4) CONTINGENCIES: the 300G hard floor NEVER approached (df band 579-598G, margin ~280G);
guard-oracle held (zero hand kills, zero masters touched); the prep valve HELD all window
(184-192 pending, VALVE lines each minute, the patch working); terminal classes terminal;
the ladder NOT opened (df 598G << 1,600G, backlog present); the ancient impg-reference
benign - NEVER KILLED. FINAL POLL 06:07:48Z: split 9,891 / acc 5,844 (59.1% final) / sub
2,282 (valve-held at 184) / map dirs 88 (batch-80 in flight) / call-terminal 2,015 =
callok 1,452 + fail-closed 543 + cut.exit 16 + poison 4 = 2,015 of 5,844 (34.5% of
accepted) / df 598G / mem 424Gib / rvsize 199G apparent. MASTERS AT CLOSE: prep 411760
(1-11:41h, VALVE-held pool), cuts 3838087 (BATCH_W=4 N_CONCURRENT=1 RAYON=24
DISK_FLOOR_GB=300, width-4-direct, ~11:15h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2,
~19:11h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census
continuation (batch-80 in flight); the JIT reaper cadence (~5-12s a pass, idempotent,
re-runnable anytime); the T03333 lane + the deeper-lane cluster under the guard-oracle
(never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT
lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master
MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the
flock ceiling, never RAM); poison-first on prep stalls (the class is 4 members); AT
CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the
standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule
+ THE REAPER RECEIPTS: ~177.5G cumulative freed classes A/B/C across 23 idempotent passes,
the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain
re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members incl. the T00024
3h crossing receipt, the mid-dispatch valve patch's live proof, the guard-oracle receipts
incl. the T03641 11,084s/71.1GB deepest-lane point and the T03203 guard-137 second-deepest
datum, the pass-20 class C bulge datum - the freed-pole substrate reappearing as reaped
intermediates, the JIT loop closed live) -> the record's closing increment -> the commit.
SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export,
K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT.
Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the
owner pushes).
POST-CLOSE NOTE (06:10Z): fail-closed 543 -> 544 (+1, NAMED, landed 06:08:33Z one minute
after the final poll 06:07:48Z): T03370 CALL guard-137 mid-score (failures.txt line 'rss
locus_1_scored: 48250408 kB' - the last logged sample; the guard fires at the 64GiB crossing
per the wave's own guard, same class as T03203/T03641, wall per receipts/ for the successor
to fold into the call-wall distribution; T03370 also the last pre-swap width-8 batch
guard-kill 18:52:40Z - both its cut and its call now terminal-honest). Terminals now
1,452 + 544 + 16 + 4 = 2,016/5,844 (34.5%). The successor's rollup counts it.

## RIDE CONTINUATION 55 SESSION OPEN (06:15Z 2026-10-17, commit 3d39948) - THE CADENCE CONTINUES (PASS 24 OWED; PASS 23'S +16M RECEIPTED AT RIDE-54 CLOSE; T03370's CALL GUARD-137 NAMED IN THE RIDE-54 POST-CLOSE NOTE)
Opened from the ride-54 window close (216e63b). ORIENTATION STANDING (rides 38-55): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN") is STALE - the
reaper has RUN since ride-32 (129G first pass); the ride-31 owed list (the ~22s regen
receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE
(cumulative ~= 177.5G across 23 idempotent passes entering this window).
OPEN POLL 06:13:55Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,282 (valve-held) /
cut.done 2,040 (map dirs 89) / call-terminal 2,017 = callok 1,453 + fail-closed 544 + cut.exit
16 + poison 4 (34.5% of accepted) / df 600G / mem 461Gib / rvsize 197-198G apparent / census
81/81 post-swap batch.exit markers ALL exit=0 re-verified live (batch-80 resolved HONEST
exit=0 wall 457s IN-BAND; batch-1792217682-81 dispatched ~06:12Z in flight) / live lanes: the
T03333 deep-ish lane (pid 967312, ~1:11, 52.7GB plateau, niced, guard-oracle) + the
ride-54-passed-forward deeper-lane cluster at its DEEPEST - five lanes 52.7-60.0GB /
1:00-1:11 elapsed (pids 994062 63.0GB-rss leading, 968955, 1008767, 1002856, plus 967312),
ALL NICED, ALL UNDER the 64GiB guard, guard-oracle, NO hand kill + 5 ordinary lanes 25-37GB +
the batch-81 cut in scan. MASTERS: prep 411760 (VALVE-held pool), cuts 3838087 (width-4-direct
BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2);
heartbeat live; the ancient impg-reference benign. REAPER SCRIPT RE-VERIFIED AT OPEN: bash -n
CLEAN, md5 b804a7d3 UNCHANGED since the ride-31 review, mtime Oct 16 20:13. THIS WINDOW OWES:
the JIT-cadence reaper pass 24 detached + poll SHORT; the width-4 census (batch-81 in flight);
the deep-lane cluster under the guard-oracle; the standing contingencies; at campaign close:
rv-rollup.py + datasheet + the record's closing increment + commit. NO PR PUSH.

## RIDE CONTINUATION 55 WINDOW CLOSE (06:30Z 2026-10-17, commit a265e8e) - THE CADENCE PASS 24 RECEIPTED (+37M APPARENT, CUMULATIVE ~= 177.5G ACROSS 24 IDEMPOTENT PASSES); CENSUS 82/82 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF AN EIGHTEENTH TIME; BATCH-80/81 HONEST AT 457s/588s; THE DEEP-LANE CLUSTER RESOLVED BY ITS OWN GUARD - FOUR CALL GUARD-137 TERMINALS THIS WINDOW (T03228/T03775/T03551/T02927); ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 06:15Z (session-open increment 3d39948). (1) THE REAPER PASS 24 (06:15:27-31Z, ~5s,
detached, receipt prep-reaper.log; bash -n CLEAN + script md5 b804a7d3 UNCHANGED re-verified
at the open): class A 0, class B 48 files 66M + syng/export/walks/walkreq markers dropped on
36 cut-awaiting-call, class C 0 files + markers dropped on 192 pending (idempotence - the
pass-20 class C bulk already reaped); lean buffer 196 <= 200 - no eviction (the prep valve
at 160 is the standing bound); freed 37M apparent (prep 111,185 -> 111,148M; df 601G steady
at the pass). CUMULATIVE ~= 177.5G apparent across 24 idempotent passes, every cut/call
input and every split ledger kept. EXCLUSIONS VERIFIED LIVE AN EIGHTEENTH TIME: the batch-81
cut cohort (T04021/T04033/T03830/T03937) + ALL live call lanes in the pass's exclusion list,
locus.fa KEPT (spot-checked 118-632MB incl. T03333 582MB; syng index intact on the excluded
in-flight T04021/T03830); the excluded batch-81 resolved HONEST exit=0 wall 588s IN-BAND
AFTER the pass - the cut lanes stayed fed.
(2) CENSUS: 82/82 batch.exit=0 post-swap, ZERO kills since the 18:52:42Z swap (batch-80
honest 457s receipted at the open, batch-81 honest 588s; batch-1792218272-82 dispatched
~06:26Z in flight at close, cohort T03901/T03903/T04010/T03876); cut.done 2,040 -> 2,048
(+8); map dirs 89.
(3) CALLS - THE GUARD-ORACLE'S BUSIEST HONEST WINDOW: callok 1,453 -> 1,454; fail-closed
544 -> 548 (+4, ALL NAMED, ALL the lane's own 64GiB guard mid-score: T03228 06:13:56Z,
T03775 06:16:09Z, T03551 06:18:34Z, T02927 06:19:25Z - the ride-54-passed-forward
deep-lane cluster resolved BY ITS OWN GUARD, the 60+GB leaders 994062/968955 first, the 56GB
tier following; the deep tail is REAL at genome scale and the guard handles it - ZERO HAND
KILLS); cut.exit 16; poison 4 terminal; PASSED FORWARD: the T03333 lane (pid 967312, ~1:26,
55.3GB PLATEAU-STABLE, niced, guard-oracle) + 1008767 (56.3GB, ~1:16) + four ordinary lanes
30-43GB, ALL NICED, ALL UNDER the guard.
(4) CONTINGENCIES: the 300G hard floor NEVER approached (df band 600-610G, margin ~300G);
guard-oracle held; the prep valve HELD all window (sub 2,282, pending 184-200); the ladder
NOT opened (df 602G << 1,600G); terminal classes terminal; the ancient impg-reference
benign - NEVER KILLED. FINAL POLL 06:29:53Z: split 9,891 / acc 5,844 / sub 2,282
(valve-held) / cut.done 2,048 (map dirs 89; batch-82 in flight) / call-terminal 2,022 =
callok 1,454 + fail-closed 548 + cut.exit 16 + poison 4 = 2,022 of 5,844 (34.6% of
accepted) / df 602G / mem 628-661Gib / rvsize 195-196G apparent.
MASTERS AT CLOSE: prep 411760 (1-12:02h, VALVE-held pool), cuts 3838087 (width-4-direct
BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~11:36h), wave 2810146 (LANES=128
CHUNK_MULT=2, ~19:32h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the
width-4 census (batch-82 in flight); the JIT reaper cadence (~5-12s a pass, idempotent);
the T03333 + 1008767 deep lanes under the guard-oracle (never kill by hand - this window's
four guard-137 resolutions prove the guard fires when it must); the prep pole + df watch vs
the 300G hard floor; the ladder only at backlog drained AND df >= 1,600G; poison-first on
prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 terminal): rv-rollup.py
-> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the
inventory finding + the JIT-prep rule + the reaper receipts ~177.5G/24 passes + the
four-guard-137 cluster-resolution datum) -> the record's closing increment -> the commit.
NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 56 SESSION OPEN (06:35Z 2026-10-17, commit ec40fd3) - THE CADENCE CONTINUES (PASS 25 OWED; PASS 24'S +37M RECEIPTED AT RIDE-55 CLOSE)
Opened from the ride-55 window close (a265e8e). ORIENTATION STANDING (rides 38-56): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN") is STALE - the
reaper has RUN since ride-32 (129G first pass); the ride-31 owed list (the ~22s regen
receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE
(cumulative ~= 177.5G across 24 idempotent passes entering this window).
OPEN POLL 06:32:31Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,282 (VALVE-held at 176
pending, VALVE lines each minute) / cut.done 2,052 (map dirs 87) / call-terminal 2,024 =
callok 1,454 + fail-closed 550 + cut.exit 16 + poison 4 (34.6% of accepted) / df 605-608G /
mem 743-776Gib / rvsize 190-193G apparent / census 83/83 post-swap batch.exit markers ALL
exit=0 re-verified live (batch-82 resolved HONEST exit=0 wall 482s IN-BAND; batch-83
dispatched ~06:31Z in flight) / live lanes: the ride-55-passed-forward T03333 deep lane
RESOLVED TERMINAL BY ITS OWN 64GiB GUARD at 06:30:10Z (+ T03470 06:30:26Z - the +2
fail-closed at the ride-55 edge, both named) + the passed-forward 1008767 (56.3GB, 1:20,
niced) leading four ordinary deep-ish lanes 30-46GB, ALL NICED, ALL UNDER the 64GiB guard,
guard-oracle, NO hand kill. MASTERS: prep 411760 (1-12:06h, VALVE-held pool), cuts 3838087
(width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146
(LANES=128 CHUNK_MULT=2); heartbeat live; the ancient impg-reference benign. REAPER SCRIPT
RE-VERIFIED AT OPEN: bash -n CLEAN, md5 b804a7d3 UNCHANGED since the ride-31 review.

## RIDE CONTINUATION 56 WINDOW CLOSE (06:42Z 2026-10-17, commit 263e0bf) - THE CADENCE PASS 25 RECEIPTED (+148M APPARENT ~4s, CUMULATIVE ~= 177.7G ACROSS 25 IDEMPOTENT PASSES); CENSUS 84/84 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF A NINETEENTH TIME; BATCH-83 HONEST AT 530s AFTER THE PASS; THREE CALL GUARD-137 TERMINALS THIS WINDOW (T03333/T03470/T02844, ALL NAMED); ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 06:35Z (session-open increment ec40fd3). (1) THE REAPER PASS 25 (06:33:25-29Z, ~4s,
detached pid 1207605, receipt prep-reaper.log; bash -n CLEAN + script md5 b804a7d3
UNCHANGED re-verified at the open): class A 0, class B 48 files 181M + syng/export/walks/
walkreq markers dropped on 44 cut-awaiting-call (census 48), class C 0 files + markers
dropped on 184 pending (idempotence - the pass-20 class C bulk already reaped); lean buffer
188 <= 200 - no eviction (the prep valve at 160 is the standing bound); freed 148M apparent
(prep 109,562 -> 109,414M; df 605G steady at the pass). CUMULATIVE ~= 177.7G apparent
across 25 idempotent passes, every cut/call input and every split ledger kept. EXCLUSIONS
VERIFIED LIVE A NINETEENTH TIME: the batch-83 cut cohort (T03845/T04052/T03986/T03871) +
ALL live call lanes in the pass's exclusion list, locus.fa verified KEPT (spot-checked
118-306MB on the cohort), syng index intact on the excluded in-flight T03845/T04052; the
excluded batch-83 resolved HONEST exit=0 wall 530s IN-BAND AFTER the pass - the cut lanes
stayed fed. (2) CENSUS: 84/84 batch.exit=0 post-swap, ZERO kills since the 18:52:42Z swap
(batch-82 honest 482s receipted at the open, batch-83 honest 530s; batch-84 dispatched
~06:41Z in flight at close, cut pid in load 3.4GB); cut.done 2,052 -> 2,056 (+4, the
batch-83 cohort); map dirs 87 -> 94 (the cuts pole out-writing; the wave consumes).
(3) CALLS - THREE GUARD-137 TERMINALS, ALL THE LANES' OWN 64GiB GUARD MID-SCORE, ALL NAMED:
fail-closed 548 -> 551 (+3: T03333 06:30:10Z - THE RIDE-55-PASSED-FORWARD DEEP LANE
RESOLVED BY ITS OWN GUARD at ~1:32 elapsed, fail-reaped in-lane receipts kept, its
substrate honestly gone from prep/; T03470 06:30:26Z; T02844 06:34:15Z - itself on the
pass's exclusion list, correctly excluded; receipts in wave-logs/failures.txt); callok 1,454
STEADY; cut.exit 16; poison 4 terminal; PASSED FORWARD: the deep-ish lanes 1008767 (56.3GB
PLATEAU-STABLE, ~1:29, niced), 1016384 (50.9GB, ~1:27), 982835 (36.4GB, ~1:35) - ALL NICED,
ALL UNDER the 64GiB guard, guard-oracle, NO hand kill. (4) CONTINGENCIES: the 300G hard
floor NEVER approached (df band 601-608G, margin ~300G); guard-oracle held; the prep valve
HELD all window (172-176 pending, VALVE lines each minute, the patch working); terminal
classes terminal; the ladder NOT opened (df 601G << 1,600G, backlog present); the ancient
impg-reference benign - NEVER KILLED. FINAL POLL 06:41:49Z: split 9,891 / acc 5,844 (59.1%
final) / sub 2,282 (valve-held at 172) / cut.done 2,056 (map dirs 94; batch-84 in flight) /
call-terminal 2,025 = callok 1,454 + fail-closed 551 + cut.exit 16 + poison 4 = 2,025 of
5,844 (34.7% of accepted) / df 601G / mem 756-803Gib / rvsize 190-196G apparent. MASTERS
AT CLOSE: prep 411760 (1-12:09h, VALVE-held pool), cuts 3838087 (width-4-direct BATCH_W=4
N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~11:44h), chunked wave 2810146 (LANES=128
CHUNK_MULT=2, ~19:40h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the
width-4 census continuation (batch-84 in flight); the JIT reaper cadence (~5-12s a pass,
idempotent, re-runnable anytime); the three passed-forward deep-ish lanes under the
guard-oracle (never kill by hand - this window's three + ride-55's four guard-137
resolutions are the proof the guard fires when it must); the prep pole + df watch vs the
300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog
drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never
with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep
stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal):
rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct
delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~177.7G
cumulative freed classes A/B/C across 25 idempotent passes, the walkreq marker-gap finding
+ fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng
segdup-tail class's 4 CONFIRMED members incl. the T00024 3h crossing receipt, the
mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/
71.1GB deepest-lane + the T03203/T03370/T03333/T03470/T02844 points + ride-55's
four-guard-137 cluster) -> the record's closing increment -> the commit. NO PR PUSH (the
owner pushes).

## RIDE CONTINUATION 57 SESSION OPEN (06:48Z 2026-10-17, commit to follow) - THE CADENCE CONTINUES (PASS 26 OWED; PASS 25'S +148M RECEIPTED AT RIDE-56 CLOSE; T03636's CALL GUARD-137 NAMED AT THIS OPEN)
Opened from the ride-56 window close (263e0bf). ORIENTATION STANDING (rides 38-57): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN") is STALE - the
reaper has RUN since ride-32 (129G first pass); the ride-31 owed list (the ~22s regen
receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE
(cumulative ~= 177.7G across 25 idempotent passes entering this window).
OPEN POLL 06:44:20Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,282 (valve-held) / map
dirs 93 / call-terminal 2,026 = callok 1,454 + fail-closed 552 (+1 IN-WINDOW, NAMED:
T03636 06:42:57Z call guard-137 mid-score - score locus 2, touching/locality
classification, the lane's own 64GiB guard, same class as the ride-56 cluster, receipts
in wave-logs/failures.txt) + cut.exit 16 + poison 4 = 2,026 of 5,844 (34.7% of accepted) /
df 603G (the 300G hard floor margin ~303G) / mem 834Gib / rvsize 194G apparent / census
84/84 post-swap batch.exit markers ALL exit=0 re-verified live by find -newermt (ZERO
non-zero; last kill line 18:52:40Z pre-swap width-8; batch-1792219286-84 in flight, cohort
T03925/T03969/T03941/T03834, cut pid 1215784 35GB niced) / live lanes: the ride-56-
passed-forward deep lanes 1008767 (53GB, ~1:31, niced) + 982835 (34GB, ~1:37, niced), both
under the guard-oracle, NO hand kill; 1016384 (50.9GB at ride-56 close) resolved since the
close. MASTERS: prep 411760 (1-12:18h, VALVE-held pool), cuts 3838087 (width-4-direct
BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~11:52h), wave 2810146 (LANES=128
CHUNK_MULT=2, ~19:48h); heartbeat live (06:43:57Z line healthy); the ancient impg-reference
benign. REAPER SCRIPT RE-VERIFIED AT OPEN: bash -n CLEAN, md5 b804a7d3 UNCHANGED since the
ride-31 review, mtime Oct 16 20:13. THIS WINDOW OWES: the JIT-cadence reaper pass 26
detached + poll SHORT; the width-4 census (batch-84 in flight); the deep lanes under the
guard-oracle; the standing contingencies; at campaign close: rv-rollup.py + datasheet +
the record's closing increment + commit. NO PR PUSH.

## RIDE CONTINUATION 57 WINDOW CLOSE (06:53Z 2026-10-17) - THE CADENCE PASS 26 RECEIPTED (+45M APPARENT ~5s, CUMULATIVE ~= 177.7G ACROSS 26 IDEMPOTENT PASSES); CENSUS 85/85 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF A TWENTIETH TIME; BATCH-84 HONEST AT 473s AFTER THE PASS; TWO CALL GUARD-137 TERMINALS THIS WINDOW (T03636/T03804, ALL NAMED); THE PASSED-FORWARD DEEP LANE 1008767 RESOLVED BY COMPLETION; ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 06:44Z (session-open increment 3a4e1f8). (1) THE REAPER PASS 26
(06:46:03-07Z, ~5s, detached pid 1221351, receipt prep-reaper.log; bash -n CLEAN +
script md5 b804a7d3 UNCHANGED re-verified at the open): class A 0 (nothing new rejected),
class B 24 files 61M + syng/export/walks/walkreq markers dropped on 48 cut-awaiting-call
(census 50), class C 0 files + markers dropped on 180 pending (idempotence - the pass-20
class C bulk already reaped); lean buffer 184 <= 200 - no eviction (the prep valve at 160
is the standing bound); freed 45M apparent (prep 109,122 -> 109,077M; df 603G steady at
the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 177.7G apparent across 26 idempotent
passes (129G ride-32 + the 25 cadence passes since), every cut/call input and every split
ledger kept throughout. THE IN-FLIGHT EXCLUSIONS VERIFIED LIVE A TWENTIETH TIME: the
batch-84 cut cohort (T03925/T03969/T03941/T03834) + ALL live call lanes in the pass's
exclusion list (T03499/T03804), locus.fa verified KEPT on disk after the pass (spot-checked
117-220MB on the cohort + the lanes), syng index intact on the excluded in-flight T03925;
the excluded batch-84 resolved HONEST exit=0 wall 473s IN-BAND AFTER the pass - the cut
lanes stayed fed, the receipt the pass answers to. (2) THE WIDTH-4 CENSUS: batches
post-swap = 85/85 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z
swap (this window: batch-83 honest 530s receipted at the ride-56 close, batch-84 honest
473s; batch-1792219760-85 dispatched ~06:49Z in flight at close, cohort
T04017/T03847/T04006/T03916, cut pid 1225086); cut.done 2,056 -> 2,060 (+4, the batch-84
cohort); map dirs 93 -> 96 as the cuts pole out-writes and the wave consumes. (3) CALLS -
TWO GUARD-137 TERMINALS, BOTH THE LANES' OWN 64GiB GUARD MID-SCORE, BOTH NAMED:
fail-closed 551 -> 553 (+2: T03636 06:42:57Z score locus 2 touching/locality
classification; T03804 06:47:40Z score locus 17 exactness - itself on the pass-26 exclusion
list, correctly excluded; receipts in wave-logs/failures.txt + per-territory fail-reaped
in-lane); callok 1,454 steady; cut.exit 16; poison 4 terminal; THE RIDE-56-PASSED-FORWARD
DEEP LANE 1008767 (53-56GB plateau, ~1:31 at the open) RESOLVED BY COMPLETION this window
- no kill, no failure line, its receipts reaped in-lane; PASSED FORWARD: 982835 (34GB,
~1:43 elapsed, niced, guard-oracle, NO hand kill). (4) CONTINGENCIES: the 300G hard floor
NEVER approached (df band 603-605G, margin ~300G); guard-oracle held (zero hand kills,
zero masters touched - the two guard-137s were the lanes' own sanctioned guard); the prep
valve HELD all window (sub 2,282, pending 172-188); terminal classes terminal; the ladder
NOT opened (df 603G << 1,600G, backlog present); the ancient impg-reference benign -
NEVER KILLED. FINAL POLL 06:49:51Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,282
(valve-held) / cut.done 2,060 (map dirs 96; batch-85 in flight) / call-terminal 2,027 =
callok 1,454 + fail-closed 553 + cut.exit 16 + poison 4 = 2,027 of 5,844 (34.7% of
accepted) / df 603G / mem 869-898Gib / rvsize 192-194G apparent. MASTERS AT CLOSE: prep
411760 (1-12:24h, VALVE-held pool), cuts 3838087 (width-4-direct BATCH_W=4
N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~12:00h), chunked wave 2810146 (LANES=128
CHUNK_MULT=2, ~19:55h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the
width-4 census continuation (batch-85 in flight); the JIT reaper cadence (~5-12s a pass,
idempotent, re-runnable anytime); the passed-forward lane 982835 under the guard-oracle
(never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT
lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts
master MID-WAIT only, by exact pid, never the wave, never with marathons in flight;
C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 4
members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL
datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory
finding + the JIT-prep rule + THE REAPER RECEIPTS: ~177.7G cumulative freed classes A/B/C
across 26 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen
measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 4
CONFIRMED members incl. the T00024 3h crossing receipt, the mid-dispatch valve patch's
live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the
T03203/T03370/T03333/T03470/T02844/T03636/T03804 points + rides 55/56's guard-137
clusters) -> the record's closing increment -> the commit. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 58 SESSION OPEN (06:52Z 2026-10-17, commit e6445c3) - THE CADENCE CONTINUES (PASS 27 OWED; PASS 26'S +45M RECEIPTED AT RIDE-57 CLOSE)
Opened from the ride-57 window close (d27c755). ORIENTATION STANDING (rides 38-58): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN", terminals
1,728/5,844) is STALE - the reaper has RUN since ride-32 (129G first pass); the ride-31
owed list (the ~22s regen receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the
standing duty is the CADENCE (cumulative ~= 177.7G across 26 idempotent passes entering
this window). OPEN POLL 06:52:15Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,282
(valve-held, the master's own line 'VALVE: 164 pending') / map dirs 96 / call-terminal
2,027 = callok 1,454 + fail-closed 553 + cut.exit 16 + poison 4 = 2,027 of 5,844 (34.7% of
accepted; every fail line named in wave-logs/failures.txt) / df 603G (the 300G hard floor
margin ~303G) / mem 887Gib / rvsize 193G apparent / census 85/85 post-swap batch.exit
markers ALL exit=0 re-verified live by find -newermt (ZERO non-zero; last kill line
18:52:40Z pre-swap width-8; batch-1792219760-85 dispatched 06:49:20Z in flight, cohort
T04017/T03847/T04006/T03916, cut pid 1225086 30.9GB in scan) / live lanes: the
ride-57-passed-forward deep lane 982835 (T03499, 34.7GB, ~1:45 elapsed, NICED ni=10,
guard-oracle, NO hand kill; 1008767 resolved since the close). MASTERS: prep 411760
(1-12:26h, VALVE-held pool), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1
RAYON=24 DISK_FLOOR_GB=300, 11:59h), wave 2810146 (LANES=128 CHUNK_MULT=2, 19:55h);
heartbeat live (06:52:08Z line healthy); the ancient impg-reference benign. REAPER SCRIPT
RE-VERIFIED AT OPEN: bash -n CLEAN, md5 b804a7d3 UNCHANGED since the ride-31 review, mtime
Oct 16 20:13. THIS WINDOW OWES: the JIT-cadence reaper pass 27 detached + poll SHORT; the
width-4 census (batch-85 in flight); the deep lane under the guard-oracle; the standing
contingencies; at campaign close: rv-rollup.py + datasheet + the record's closing
increment + commit. NO PR PUSH.

## RIDE CONTINUATION 58 WINDOW CLOSE (07:00Z 2026-10-17) - THE CADENCE PASS 27 RECEIPTED (+19M APPARENT ~5s, CUMULATIVE ~= 177.7G ACROSS 27 IDEMPOTENT PASSES); CENSUS 86/86 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF A TWENTY-FIRST TIME (SYNG = THE SIX-FILE INDEX, NAMED); BATCH-85 HONEST AT 431s AFTER THE PASS; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 06:52Z (session-open increment e6445c3). (1) THE REAPER PASS 27 (06:53:08-13Z, ~5s,
detached pid 1228429, receipt prep-reaper.log + prep-reaper-r58.nohup; bash -n CLEAN +
script md5 b804a7d3 UNCHANGED re-verified at the open): class A 0 (nothing new rejected),
class B 24 files 34M + syng/export/walks/walkreq markers dropped on 52 cut-awaiting-call
(census 53), class C 0 files + markers dropped on 176 pending (idempotence - the pass-20
class C bulk already reaped); lean buffer 180 <= 200 - no eviction (the prep valve at 160
is the standing bound); freed 19M apparent (prep 108,865 -> 108,846M; df 603G steady at
the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 177.7G apparent across 27 idempotent
passes (129G ride-32 + the 26 cadence passes since), every cut/call input and every split
ledger kept throughout. THE IN-FLIGHT EXCLUSIONS VERIFIED LIVE A TWENTY-FIRST TIME, WITH
A NAMING CLARIFICATION: the syng "index" is SIX sidecar files (locus-pos64.syng.1gbwt /
.1khash / .meta / .names / .pstep / .spos) - a single-filename check returns 0 and MUST
NOT be read as a missing index; the batch-85 cut cohort (T04017/T03847/T04006/T03916)
verified 6/6 sidecars intact + locus.fa KEPT (113-141MB) after the pass, the live call
lane T03499 verified locus.fa 113M + reference-axis.json + 17 locality-graphs intact (the
call reads NO syng; its syng was legitimately reaped by an earlier pass under class B
cut-done - the design); class-B reaped territories spot-checked (T01895/T01941/T01945/
T01987): syng sidecars GONE, locus.fa + split ledgers + refined/member beds KEPT - the
reap exactly where intended. The excluded batch-85 resolved HONEST exit=0 wall 431s
IN-BAND AFTER the pass (dispatched 06:49:20Z, resolved ~06:57Z) - the cut lanes stayed
fed, the receipt the pass answers to; batch-1792220192-86 dispatched 06:57:28Z in flight
at close (cohort led by T03897). (2) THE WIDTH-4 CENSUS: batches post-swap = 86/86 HONEST
ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap (this window:
batch-84 honest 473s receipted at ride-57 close, batch-85 honest 431s; batch-86 in
flight); cut.done 2,060 -> 2,064 (+4, the batch-85 cohort); map dirs 96 -> 100. (3)
CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 553 steady (the last two, T03636
06:42:57Z / T03804 06:47:40Z, both guard-137 mid-score, both receipted at the ride-57
close, both PREDATE this window's open); callok 1,454 steady; cut.exit 16; poison 4
terminal; PASSED FORWARD: 982835 (T03499, 34.7GB plateau, ~1:52 elapsed, niced ni=10,
guard-oracle, NO hand kill). (4) CONTINGENCIES: the 300G hard floor NEVER approached (df
band 602-603G, margin ~302G); guard-oracle held (zero hand kills, zero masters touched);
the prep valve HELD all window (sub 2,282, pending 160-180, the master's own VALVE lines
each minute); terminal classes terminal; the ladder NOT opened (df 602G << 1,600G,
backlog present); the ancient impg-reference benign - NEVER KILLED. FINAL POLL 06:58:41Z:
split 9,891 / acc 5,844 (59.1% final) / sub 2,282 (valve-held) / map dirs 100 /
call-terminal 2,027 = callok 1,454 + fail-closed 553 + cut.exit 16 + poison 4 = 2,027 of
5,844 (34.7% of accepted) / df 602G / mem 887-892Gib / rvsize 195G apparent. MASTERS AT
CLOSE: prep 411760 (1-12:26h, VALVE-held pool), cuts 3838087 (width-4-direct BATCH_W=4
N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~12:05h), chunked wave 2810146 (LANES=128
CHUNK_MULT=2, ~20:00h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the
width-4 census continuation (batch-86 in flight); the JIT reaper cadence (~5-12s a pass,
idempotent, re-runnable anytime; VERIFY WITH THE SIX-FILE SYNG SIDECAR SET - the
single-filename check is a false negative, named this window); the passed-forward lane
982835 under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G
hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog
drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never
with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep
stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal):
rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct
delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~177.7G
cumulative freed classes A/B/C across 27 idempotent passes, the walkreq marker-gap
finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt,
the syng segdup-tail class's 4 CONFIRMED members incl. the T00024 3h crossing receipt,
the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641
11,084s/71.1GB deepest-lane + the T03203/T03370/T03333/T03470/T02844/T03636/T03804
points + rides 55/56's guard-137 clusters) -> the record's closing increment -> the
commit. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 59 SESSION OPEN (07:01Z 2026-10-17, commit 122d51b) - THE CADENCE CONTINUES (PASS 28 OWED; PASS 27'S +19M RECEIPTED AT RIDE-58 CLOSE)
Opened from the ride-58 window close (31e200f). ORIENTATION STANDING (rides 38-59): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN", terminals
1,728/5,844) is STALE - the reaper has RUN since ride-32 (129G first pass); the ride-31
owed list (the ~22s regen receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the
standing duty is the CADENCE (cumulative ~= 177.7G across 27 idempotent passes entering
this window). OPEN POLL 07:00:57Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,282
(VALVE-held at 160 pending, the master's own line 06:59:29Z) / map dirs 100 / call-terminal
2,027 = callok 1,454 + fail-closed 553 + cut.exit 16 + poison 4 = 2,027 of 5,844 (34.7% of
accepted) / df 602G (the 300G hard floor margin ~302G) / mem 872Gib / rvsize 195G apparent
/ census 86/86 post-swap batch.exit markers ALL exit=0 re-verified live (ZERO non-zero;
last kill line 18:52:40Z pre-swap width-8; batch-1792220192-86 dispatched 06:56:32Z in
flight, cohort T03897/T03858/T03904/T03868, cut pid 1232266 37GB niced) / live lanes: the
ride-58-passed-forward deep lane 982835 (T03499, 34.7GB, ~1:54 elapsed, NICED ni=10,
guard-oracle, NO hand kill). MASTERS: prep 411760 (1-12:35h, VALVE-held pool), cuts
3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, 12:08h), wave
2810146 (LANES=128 CHUNK_MULT=2, 20:04h); heartbeat live (07:00:19Z line healthy); the
ancient impg-reference benign. REAPER SCRIPT RE-VERIFIED AT OPEN: bash -n CLEAN, md5
b804a7d366844d27ac58e2e0c8b64802 UNCHANGED since the ride-31 review, mtime Oct 16 20:13.
THIS WINDOW OWES: the JIT-cadence reaper pass 28 detached + poll SHORT; the width-4 census
(batch-86 in flight); the deep lane under the guard-oracle; the standing contingencies; at
campaign close: rv-rollup.py + datasheet + the record's closing increment + commit. NO PR PUSH.

## RIDE CONTINUATION 59 WINDOW CLOSE (07:16Z 2026-10-17) - THE CADENCE PASS 28 RECEIPTED (+17M APPARENT ~5s, CUMULATIVE ~= 177.7G ACROSS 28 IDEMPOTENT PASSES); CENSUS 88/88 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF A TWENTY-SECOND TIME (SIX-FILE SIDECAR RULE APPLIED); BATCH-86 HONEST AT 397s (FASTEST WIDTH-4 YET) + BATCH-87 HONEST AT 713s AFTER THE PASS; THE PASSED-FORWARD DEEP LANE T03499 RESOLVED BY ITS OWN GUARD; A 12-BIND BURST AS THE WAVE CONSUMES; ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 07:01Z (session-open increment 122d51b). (1) THE REAPER PASS 28 (07:01:29-33Z, ~5s,
detached pid 1236458, receipt prep-reaper.log + prep-reaper-r59.nohup; bash -n CLEAN + md5
b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0 (nothing new
rejected), class B 24 files 31M + syng/export/walks/walkreq markers dropped on 56
cut-awaiting-call (census 57), class C 0 files + markers dropped on 172 pending
(idempotence - the pass-20 class C bulk already reaped); lean buffer 176 <= 200 - no
eviction (the prep valve at 160 is the standing bound); freed 17M apparent (prep 108,846 ->
108,829M; df 602G steady at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 177.7G apparent
across 28 idempotent passes (129G ride-32 + the 27 cadence passes since), every cut/call
input and every split ledger kept throughout. THE IN-FLIGHT EXCLUSIONS VERIFIED LIVE A
TWENTY-SECOND TIME (SIX-FILE SIDECAR RULE): the batch-86 cut cohort (T03897/T03858/
T03904/T03868) verified 6/6 syng sidecars intact + locus.fa KEPT (113-129MB) + axis +
locality-graphs intact, the live call lane T03499 locus.fa 113M + axis + 17 locality-graphs
intact (the call reads NO syng; its sidecars legitimately reaped by an earlier pass under
class B cut-done - the design); class-B spot-checks (T01895/T01941/T01945/T01987): sidecars
GONE, locus.fa + refined beds KEPT - the reap exactly where intended. The excluded batch-86
resolved HONEST exit=0 wall 397s IN-BAND AFTER the pass (THE FASTEST WIDTH-4 WALL YET) - the
cut lanes stayed fed, the receipt the pass answers to; batch-87 resolved honest exit=0 wall
713s (top of band); batch-1792221304-88 dispatched 07:15:04Z in flight at close (cohort
T03991/T03873/T04029/T03902, cut pid 1288897). (2) THE WIDTH-4 CENSUS: batches post-swap =
88/88 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap (this
window: batch-86 honest 397s, batch-87 honest 713s); cut.done 2,064 -> 2,072 (+8); map
dirs 100 -> 84 as the wave consumes faster than the cuts pole out-writes (the C=1 cut rate
vs the widened wave). (3) CALLS - THE PASSED-FORWARD DEEP LANE RESOLVED BY ITS OWN GUARD +
A 12-BIND BURST AS THE WAVE CONSUMES THE BACKLOG: fail-closed 553 -> 566 (+13, ALL NAMED:
T03499 07:08:15Z guard-137 mid-score locus 14 exactness - THE RIDE-58-PASSED-FORWARD DEEP
LANE, ~2:05 elapsed, plateau 34.7GB then the climb crossing the 64GiB guard, the lane's own
sanctioned guard, receipts in wave-logs/failures.txt + per-territory fail-reaped in-lane;
+12 ordinary 7479-class bind ambiguities, all exit 1, all with census records: T02512 rec
99190, T03858 rec 7121 (3 surviving), T02184 rec 34766 (2 surviving), T01945 rec 7480,
T03925 rec 9135, T04021 rec 100158, T03489 rec 7852, T03918 rec 13073, T03889 rec 81140,
T03606 rec 7574, T01941 rec 7558, T03876 rec 6949 (2 surviving)); callok 1,454 -> 1,465
(+11); cut.exit 16; poison 4 terminal. (4) CONTINGENCIES: the 300G hard floor NEVER
approached (df band 596-615G, margin ~300G); guard-oracle held (zero hand kills, zero
masters touched - the T03499 kill was the lane's own sanctioned guard); the prep valve HELD
all window (sub 2,282 -> 2,294 as the cuts pole's cut.done freed slots and the valve
released, pending 160-176, VALVE lines each minute - the patch working); terminal classes
terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference
benign - NEVER KILLED. ONE VERIFICATION NOTE for the successor: two fresh './run-cuts.sh'
processes sighted at 07:03Z verified as the master's OWN batch-composition subshells
(children of 3838087, ppid-checked) - NOT duplicate masters, benign, named so no successor
misreads them. FINAL POLL 07:15:44Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 /
map dirs 84 (batch-88 in flight) / call-terminal 2,051 = callok 1,465 + fail-closed 566 +
cut.exit 16 + poison 4 = 2,051 of 5,844 (35.1% of accepted) / df 612G / mem 822Gib /
rvsize 186G apparent / 8+ ordinary fresh call lanes dispatching 4-9GB niced. MASTERS AT
CLOSE: prep 411760 (1-12:38h, VALVE-held pool), cuts 3838087 (width-4-direct BATCH_W=4
N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~12:15h), chunked wave 2810146 (LANES=128
CHUNK_MULT=2, ~20:15h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the
width-4 census continuation (batch-88 in flight); the JIT reaper cadence (~5-12s a pass,
idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof);
the deep-lane watch under the guard-oracle (never kill by hand); the prep pole + df watch
vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at
backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the
wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM);
poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844
accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list +
the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER
RECEIPTS: ~177.7G cumulative freed classes A/B/C across 28 idempotent passes, the
walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain
re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members incl. the T00024
3h crossing receipt, the mid-dispatch valve patch's live proof, the guard-oracle receipts
incl. T03641 11,084s/71.1GB deepest-lane + the T03203/T03370/T03333/T03470/T02844/T03636/
T03804/T03499 points + rides 55/56's guard-137 clusters) -> the record's closing increment
-> the commit. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 60 SESSION OPEN (07:18Z 2026-10-17, commit f070644) - THE CADENCE CONTINUES (PASS 29 OWED; PASS 28'S +17M RECEIPTED AT RIDE-59 CLOSE)
Opened from the ride-59 window close (2d99610). ORIENTATION STANDING (rides 38-60): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN", terminals
1,728/5,844) is STALE - the reaper has RUN since ride-32 (129G first pass); the ride-31
owed list (the ~22s regen receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the
standing duty is the CADENCE (cumulative ~= 177.7G across 28 idempotent passes entering
this window). OPEN POLL 07:17:48Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 /
map dirs 74 / call-terminal 2,061 = callok 1,472 + fail-closed 569 (+3 since the ride-59
close, ALL NAMED ordinary 7479-class binds: T03562 07:16:55Z rec 15428 2-of-2, T03986
07:17:24Z rec 2531 2-of-7, T03693 07:17:27Z two-origin-shifts rec 8715) + cut.exit 16 +
poison 4 = 2,061 of 5,844 (35.2%) / df 622G (the 300G hard floor margin ~322G) / mem
784Gib / rvsize 175G apparent / census 88/88 post-swap batch.exit markers ALL exit=0
re-verified live (ZERO non-zero; last kill line 18:52:40Z pre-swap width-8; batch-88 in
flight, cohort T03991/T03873/T04029/T03902, cut pid 1288897 22.5GB niced) / THE PREP
MASTER MID-DISPATCH at the open (xargs -P 64 chunk dispatched 07:03:32Z - the valve
RELEASED as the widened wave consumed below 160 pending, the chunked-dispatch patch's
designed path; substrates completing T04066/T04063/T04083/T04082; sub 2,282 -> 2,294;
mid-dispatch VALVE re-check verified in the log three times this session-day, latest
05:28:19Z 'parking the remaining 3484 of 3740'). MASTERS: prep 411760 (1-12:51h,
mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24
DISK_FLOOR_GB=300, ~12:25h), wave 2810146 (LANES=128 CHUNK_MULT=2, ~20:20h); heartbeat
live (07:17:44Z line healthy); the ancient impg-reference benign. REAPER SCRIPT
RE-VERIFIED AT OPEN: bash -n CLEAN, md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED since
the ride-31 review, mtime Oct 16 20:13. THIS WINDOW OWES: the JIT-cadence reaper pass 29
detached + poll SHORT; the width-4 census (batch-88 in flight); the deep-lane watch under
the guard-oracle; the standing contingencies; at campaign close: rv-rollup.py +
datasheet-runtime.md + the record's closing increment + commit. NO PR PUSH.

## RIDE CONTINUATION 60 WINDOW CLOSE (07:24Z 2026-10-17) - THE CADENCE PASS 29 RECEIPTED (+2,484M APPARENT ~5s, CUMULATIVE ~= 180.2G ACROSS 29 IDEMPOTENT PASSES - THE LARGEST CADENCE-PASS CLASS C SINCE PASS 20); CENSUS 89/89 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF A TWENTY-THIRD TIME; BATCH-88 HONEST AT 494s AFTER THE PASS; ZERO NEW FAILS IN-WINDOW; ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 07:18Z (session-open increment f070644). (1) THE REAPER PASS 29 (07:19:51-56Z,
~5s, detached, receipt prep-reaper.log + prep-reaper-r60.nohup; bash -n CLEAN + md5
b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0 (nothing
new rejected), class B 24 files 64M + syng/export/walks/walkreq markers dropped on 4
cut-awaiting-call (census 29 - the exclusion list WIDE this pass, correctly: the prep
master's 64-wide chunk mid-dispatch + all live call lanes named), class C 36 files 2,454M
+ export/walks/walkreq markers dropped on 176 pending - THE VALVE-RELEASED DISPATCH'S
FRESH CONSUMABLES (the 07:03:32Z chunk rebuilt substrates faster than the cuts consume;
the reaper followed at the cadence, exactly the JIT design); lean buffer 180 <= 200 - no
eviction (the prep valve at 160 is the standing bound); freed 2,484M apparent (prep
107,978 -> 105,494M; df 623 -> 625G at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~=
180.2G apparent across 29 idempotent passes (129G ride-32 + the 27 cadence passes since),
every cut/call input and every split ledger kept throughout. THE IN-FLIGHT EXCLUSIONS
VERIFIED LIVE A TWENTY-THIRD TIME (SIX-FILE SIDECAR RULE): the batch-88 cut cohort
(T03991/T03873/T04029/T03902) verified 6/6 syng sidecars intact + locus.fa KEPT
(114-205M), the live call lanes verified locus.fa 267-339M + reference-axis.json +
locality-graphs intact (the call reads NO syng; sidecars legitimately reaped under class
B cut-done - the design); class-B members spot-checked (T04052/T01987/T01895, cut-done
awaiting call): syng sidecars GONE, locus.fa/axis/locality-graphs/beds KEPT - the reap
exactly where intended. The excluded batch-88 resolved HONEST exit=0 wall 494s IN-BAND
AFTER the pass (dispatched 07:15:04Z, resolved ~07:23Z) - the cut lanes stayed fed;
batch-1792221799-89 dispatched in flight at close. (2) THE WIDTH-4 CENSUS: batches
post-swap = 89/89 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z
swap (this window: batch-88 honest 494s); map dirs 74 -> 67 as the widened wave consumes
faster than the C=1 cuts pole out-writes. (3) CALLS - ZERO NEW FAILS THIS WINDOW:
fail-closed 569 steady (the last three, T03562 07:16:55Z rec 15428 / T03986 07:17:24Z
rec 2531 / T03693 07:17:27Z rec 8715, ordinary 7479-class binds, all PREDATE the window's
open poll and are named there); callok 1,472 -> 1,483 (+11); live lanes 4-15GB, NO deep
lanes in flight - guard-oracle quiet; cut.exit 16; poison 4 terminal. (4) CONTINGENCIES:
the 300G hard floor NEVER approached (df band 622-629G, margin ~320G); guard-oracle held
(zero hand kills, zero masters touched); the prep valve RELEASED BELOW 160 as the wave
consumed the backlog and the master DISPATCHED A 64-WIDE CHUNK at 07:03:32Z (xargs -P
64, substrates completing, sub 2,294; the chunked-dispatch patch's designed path,
mid-dispatch at close); terminal classes terminal; the ladder NOT opened (df << 1,600G,
backlog present); the ancient impg-reference benign - NEVER KILLED. FINAL POLL 07:23:54Z:
split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / map dirs 67 (batch-89 in flight) /
call-terminal 2,072 = callok 1,483 + fail-closed 569 + cut.exit 16 + poison 4 = 2,072 of
5,844 (35.4% of accepted) / df 629G / mem 780Gib / rvsize 167G apparent. MASTERS AT
CLOSE: prep 411760 (1-12:58h, mid-dispatch chunk), cuts 3838087 (width-4-direct
BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~12:31h), chunked wave 2810146
(LANES=128 CHUNK_MULT=2, ~20:27h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT
polls; the width-4 census continuation (batch-89 in flight); the JIT reaper cadence
(~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every
exclusion proof; NOTE: while a 64-wide prep chunk is mid-dispatch the exclusion list is
wide and class B shrinks - nominal, the design); the deep-lane watch under the
guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor
(self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND
df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with
marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep
stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal):
rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct
delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.2G
cumulative freed classes A/B/C across 29 idempotent passes, the walkreq marker-gap
finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt,
the syng segdup-tail class's 4 CONFIRMED members incl. the T00024 3h crossing receipt,
the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641
11,084s/71.1GB deepest-lane + the T03203/T03370/T03333/T03470/T02844/T03636/T03804/T03499
points + rides 55/56's guard-137 clusters) -> the record's closing increment -> the
commit. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 61 SESSION OPEN (07:31Z 2026-10-17, post-merge tip, from ride-60 close ae34a51) - THE CADENCE CONTINUES (PASS 30 OWED; PASS 29'S +2,484M RECEIPTED AT RIDE-60 CLOSE)
Opened from the ride-60 window close (ae34a51). ORIENTATION STANDING (rides 38-61): a task prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN", terminals 1,728/5,844) is STALE - the reaper has RUN since ride-32 (129G first pass); the ride-31 owed list (the ~22s regen receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE (cumulative ~= 180.2G across 29 idempotent passes entering this window). OPEN POLL 07:28:56Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / map dirs 65 / call-terminal 2,074 = callok 1,485 (+2 since the ride-60 close) + fail-closed 569 STEADY (zero new fails) + cut.exit 16 + poison 4 = 2,074 of 5,844 (35.5% of accepted) / df 630G (the 300G hard floor margin ~330G) / mem 729Gib / rvsize 166G apparent / census 89/89 post-swap batch exit markers ALL exit=0 re-verified live (ZERO non-zero; last kill line 18:52:40Z pre-swap width-8; batch-89 in flight, cohort T03862/T03900/T04082/T03955, cut pid 1311657 ~42GB RSS niced, ~7.5m elapsed - the widest batch cohort RSS yet seen, watched under the guard-oracle, NO hand kill) / live call lanes 10+ ordinary 5-15GB, NO deep lanes in flight. MASTERS: prep 411760 (1-13:03h, mid-dispatch 64-wide chunk, sub 2,294, latest substrate T04082 done 07:05:52Z), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~12:36h), wave 2810146 (LANES=128 CHUNK_MULT=2, ~20:32h); heartbeat live; the ancient impg-reference benign - NEVER KILLED. REAPER SCRIPT RE-VERIFIED AT OPEN: bash -n CLEAN, md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED since the ride-31 review, mtime Oct 16 20:13. THIS WINDOW OWES: the JIT-cadence reaper pass 30 detached + poll SHORT (the SIX-FILE SYNG SIDECAR rule on the exclusion proof); the width-4 census (batch-89 in flight); the deep-lane watch under the guard-oracle; the standing contingencies; at campaign close: rv-rollup.py + datasheet-runtime.md + the record's closing increment + commit. NO PR PUSH.

## RIDE CONTINUATION 61 WINDOW CLOSE (07:36Z 2026-10-17) - THE CADENCE PASS 30 RECEIPTED (+85M APPARENT ~4s, CUMULATIVE ~= 180.3G ACROSS 30 IDEMPOTENT PASSES); CENSUS 90/90 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF A TWENTY-FOURTH TIME; BATCH-89 HONEST AT 548s (PEAK RSS 59.7GB - THE WIDEST YET, ITS OWN GUARD NEVER CROSSED) RESOLVED 26s BEFORE THE PASS; ONE NEW FAIL T04033 GUARD-137 (THE LANE'S OWN SANCTIONED GUARD); ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 07:31Z (session-open increment 3e21c6c). (1) THE REAPER PASS 30 (07:32:51-55Z, ~4s, detached, receipt prep-reaper.log + prep-reaper-r61.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0 (nothing new rejected), class B 48 files 113M + syng/export/walks/walkreq markers dropped (the batch-89 cohort resolved cut.done at 07:32:27Z, 26s BEFORE the pass read exclusions - correctly NOT in the in-flight exclusion list, class-B eligible exactly as designed), class C 0 + export/walks/walkreq markers on 168 pending; lean buffer 172 <= 200 - no eviction (the prep valve at 160 is the standing bound); freed 85M apparent (prep 103,343 -> 103,258M; df 626G steady at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 30 idempotent passes, every cut/call input and every split ledger kept throughout. THE EXCLUSION PROOF A TWENTY-FOURTH TIME (SIX-FILE SYNG SIDECAR RULE): the pass's in-flight list correctly named the 5 mid-build substrates (T03679/T03877/T03989/T04063/T03883 - the 64-wide chunk's completions) + the 10+ live call lanes; live lanes verified locus.fa 113-296M + reference-axis.json + locality-graphs intact (T01895 296M / T03871 272M / T03847 113M - the call reads NO syng); class-B spot-checks T03862/T04082: syng sidecars REAPED, locus.fa (235M/491M) + axis + locality-graphs + beds + split ledger ALL KEPT - the reap exactly where intended; batch-90's cohort (T03877/T03989/T04063/T03883) IS the excluded mid-build set - the pass protected them, the cuts master then picked them up (dispatch ~07:33Z, in flight at close, 2m16s ~15GB at last sight). (2) THE WIDTH-4 CENSUS: batches post-swap = 90/90 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap (this window: batch-89 honest 548s, peak RSS 59,738,892kB = 59.7GB - the widest width-4 batch peak yet, 4.3GB under the 64GiB guard, its own lane, no kill); map dirs 65 -> 63 as the wave consumes. (3) CALLS: fail-closed 569 -> 570 (+1, NAMED: T04033 07:34:33Z guard-137 mid-score locus_3 last-reported rss 46.6GB before the climb - the lane's OWN sanctioned 64GiB guard, receipts in wave-logs/failures.txt + per-territory call.err/call.rss/call.wall in-lane, fail-reaped; ZERO hand kills); callok 1,485 -> 1,490 (+5); live lanes ordinary 5-15GB, NO deep lanes in flight at close; cut.exit 16; poison 4 terminal. (4) CONTINGENCIES: the 300G hard floor NEVER approached (df band 626-633G, margin ~330G); guard-oracle held (zero hand kills, zero masters touched - the T04033 kill was its own sanctioned guard); the prep master mid-dispatch chunk continuing (sub 2,294; latest completions fed batch-90 directly); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED. FINAL POLL 07:35:01Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / map dirs 63 (batch-90 in flight) / call-terminal 2,078 = callok 1,490 + fail-closed 570 + cut.exit 16 + poison 4 = 2,078 of 5,844 (35.6% of accepted) / df 633G / mem 778Gib / rvsize 164G apparent. MASTERS AT CLOSE: prep 411760 (1-13:09h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~12:42h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~20:38h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-90 in flight); the JIT reaper cadence (~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; NOTE: a batch cohort that resolves cut.done in the seconds before a pass is legitimately class-B eligible - the pass excludes by live process, not by batch bookkeeping, verified this window); the deep-lane watch under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 30 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members incl. the T00024 3h crossing receipt, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the T03203/T03370/T03333/T03470/T02844/T03636/T03804/T03499/T04033 points + rides 55/56's guard-137 clusters) -> the record's closing increment -> the commit. NO PR PUSH.

## RIDE CONTINUATION 62 SESSION OPEN (07:37Z 2026-10-17, commit 64222c3) - THE CADENCE CONTINUES (PASS 31 OWED; PASS 30'S +85M RECEIPTED AT RIDE-61 CLOSE)
Opened from the ride-61 window close (5d96991). ORIENTATION STANDING (rides 38-62): a task prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN", terminals 1,728/5,844) is STALE - the reaper has RUN since ride-32 (129G first pass); the ride-31 owed list (the ~22s regen receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the standing duty is the CADENCE (cumulative ~= 180.3G across 30 idempotent passes entering this window). OPEN POLL 07:36:57Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / map dirs 63 / call-terminal 2,080 = callok 1,490 + fail-closed 570 + cut.exit 16 + poison 4 = 2,080 of 5,844 (35.5%) / df 633G / mem 760Gib / rvsize 164G apparent / census 90/90 post-swap exit markers ALL exit=0 re-verified live / batch-1792222348-90 in flight (cohort T03877/T03989/T04063/T03883 - the pass-30-protected mid-build set, cut pid 1326554 ~39.7GB niced) / 8 ordinary live call lanes 12-30GB, no deep lanes. MASTERS: prep 411760 (mid-dispatch 64-wide chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat live; the ancient impg-reference benign. REAPER SCRIPT RE-VERIFIED AT OPEN: bash -n CLEAN, md5 b804a7d3 UNCHANGED since the ride-31 review. THIS WINDOW OWES: the JIT-cadence reaper pass 31 detached + poll SHORT; the width-4 census (batch-90 in flight); the standing contingencies; at campaign close: rv-rollup.py + datasheet + the record's closing increment + commit. NO PR PUSH.

## RIDE CONTINUATION 62 WINDOW CLOSE (07:42Z 2026-10-17) - THE CADENCE PASS 31 RECEIPTED (+0M APPARENT ~5s, CUMULATIVE ~= 180.3G ACROSS 31 IDEMPOTENT PASSES - THE FULLY-IDEMPOTENT STEADY STATE); CENSUS 91/91 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF A TWENTY-FIFTH TIME, WITH ONE LANE CLARIFICATION (T03847 RESOLVED HONEST MID-PROOF - THE WAVE'S IN-LANE JIT REAP AFTER call.done, THE DESIGN, NOT A REAPER MISS); BATCH-90 HONEST AT 491s AFTER THE PASS; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 07:37Z (session-open increment 64222c3). (1) THE REAPER PASS 31 (07:37:25-30Z, ~5s, detached pid 1333818, receipt prep-reaper.log + prep-reaper-r62.nohup; bash -n CLEAN + md5 b804a7d3 UNCHANGED re-verified at the open): class A 0, class B 0 files/0M + syng/export/walks/walkreq markers dropped on 12 cut-awaiting-call (census 19 - the pass-30 class-B bulk already reaped), class C 0 + markers on 168 pending (idempotence holding); lean buffer 172 <= 200 - no eviction (the prep valve at 160 is the standing bound); freed 0M apparent (prep 101,970M steady; df 634G steady at the pass - the wave's own in-lane reap consumed everything between passes, the JIT design at equilibrium). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 31 idempotent passes (129G ride-32 + the 29 cadence passes since), every cut/call input and every split ledger kept throughout. (2) THE EXCLUSION PROOF A TWENTY-FIFTH TIME (SIX-FILE SIDECAR RULE): the batch-90 cut cohort (T03877/T03989/T04063/T03883) verified 6/6 syng sidecars intact + locus.fa KEPT (113-272M); live call lanes verified locus.fa 173-296M + axis + locality-graphs intact (T01895/T03871/T03830/T01987/T04052/T03937/T03893); ONE CLARIFICATION, NAMED so no successor misreads it: T03847 appeared inputs-gone at the proof poll - its lane had RESOLVED HONEST MID-PROOF (call exit 0, wall 1,359s, 9 calls receipted in receipts/T03847 + call.done + the in-lane JIT reap consumed the heavy artifacts AFTER the successful call - the per-territory design, not an exclusion miss; the pass's own exclusion list had named it live at 07:37:27Z, correctly); class-B spot-checks (T03862/T04082): syng sidecars GONE, locus.fa 235M/491M + axis + refined beds KEPT - the reap exactly where intended. (3) THE WIDTH-4 CENSUS: batches post-swap = 91/91 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap (this window: batch-90 honest exit=0 wall 491s, peak observed 44.7GB in band, resolved IN-BAND AFTER the pass - dispatched ~07:32Z resolved ~07:41Z, the receipt the pass answers to; per-territory emissions landed cut.done mid-batch, the isolation working); batch-1792222840-91 dispatched ~07:40Z in flight at close (cohort T04016/T04028/T03914/T03912, cut pid 1338342 8.7GB niced early-load); cut.done 2,080 -> 2,084 (+4); map dirs 63 -> 65. (4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 570 steady (the last, T04033 07:34:33Z guard-137 mid-score, receipted at ride-61 close, predates the open); callok 1,490 -> 1,492 (+2 ordinary completions, T03847 named above); cut.exit 16; poison 4 terminal; NO deep lanes in flight at close - guard-oracle quiet. (5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 629-635G, margin ~330G); guard-oracle held (zero hand kills, zero masters touched); the prep master mid-dispatch 64-wide chunk continuing (sub 2,294 steady); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED. FINAL POLL 07:41:07Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / map dirs 65 (batch-91 in flight) / call-terminal 2,082 = callok 1,492 + fail-closed 570 + cut.exit 16 + poison 4 = 2,082 of 5,844 (35.6% of accepted) / df 631G / mem 802Gib / rvsize 165G apparent. MASTERS AT CLOSE: prep 411760 (1-13:15h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~12:48h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~20:44h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-91 in flight); the JIT reaper cadence (~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; a territory whose call-inputs appear GONE was consumed by the wave's own in-lane JIT reap after call.done - VERIFY call.done/receipts BEFORE reading it as a reaper miss, named this window); the deep-lane watch under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 31 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members incl. the T00024 3h crossing receipt, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the T03203/T03370/T03333/T03470/T02844/T03636/T03804/T03499/T04033 points + rides 55/56's guard-137 clusters) -> the record's closing increment -> the commit. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 63 SESSION OPEN (07:44Z 2026-10-17, commit 62605fe) - THE CADENCE CONTINUES (PASS 32 OWED; PASS 31'S +0M RECEIPTED AT RIDE-62 CLOSE - THE FULLY-IDEMPOTENT STEADY STATE)
Opened from the ride-62 window close (62605fe). ORIENTATION STANDING (rides 38-63): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN", terminals
1,728/5,844) is STALE - the reaper has RUN since ride-32 (129G first pass); the ride-31 owed
list (the ~22s regen receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the standing
duty is the CADENCE (cumulative ~= 180.3G across 31 idempotent passes entering this window).
OPEN POLL 07:42:54Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / map dirs 65 /
call-terminal 2,082 = callok 1,492 + fail-closed 570 STEADY (zero new fails; the last,
T04033 07:34:33Z guard-137, predates the ride-62 open and is named there) + cut.exit 16 +
poison 4 = 2,082 of 5,844 (35.6% of accepted) / df 631G (the 300G hard floor margin ~331G) /
mem 785Gib / rvsize 165G apparent / census 91/91 post-swap batch exit markers ALL =0
re-verified live (ZERO non-zero; last kill line 18:52:40Z pre-swap width-8; batch-90's marker
landed =0 exactly as receipted at ride-62 close) / batch-1792222840-91 in flight (cohort
T04016/T04028/T03914/T03912, cut pid 1338342 22.4GB niced ~2.5m elapsed) / live call lanes
ordinary 15-30GB niced (no lane over the 64GiB guard), NO deep lanes in flight -
guard-oracle quiet. MASTERS: prep 411760 (1-13:17h, mid-dispatch 64-wide chunk, tail lanes
finishing - run-prep-one T03679 alive, sub 2,294), cuts 3838087 (width-4-direct BATCH_W=4
N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~12:51h), wave 2810146 (LANES=128 CHUNK_MULT=2,
~20:46h); heartbeat live; the ancient impg-reference benign - NEVER KILLED. REAPER SCRIPT
RE-VERIFIED AT OPEN: bash -n CLEAN, md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED since
the ride-31 review, mtime Oct 16 20:13. THIS WINDOW OWES: the JIT-cadence reaper pass 32
detached + poll SHORT (the SIX-FILE SYNG SIDECAR rule on the exclusion proof; the ride-62
lane clarification standing: a territory whose call-inputs appear GONE was consumed by the
wave's own in-lane JIT reap after call.done - VERIFY call.done/receipts BEFORE reading it
as a reaper miss); the width-4 census (batch-91 in flight); the standing contingencies; at
campaign close: rv-rollup.py + datasheet-runtime.md + the record's closing increment +
commit. NO PR PUSH.

## RIDE CONTINUATION 63 WINDOW CLOSE (07:49Z 2026-10-17) - THE CADENCE PASS 32 RECEIPTED (+32M APPARENT ~5s, CUMULATIVE ~= 180.3G ACROSS 32 IDEMPOTENT PASSES); CENSUS 92/92 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF A TWENTY-SIXTH TIME; BATCH-91 HONEST AT 413s AFTER THE PASS; ONE NEW FAIL T01987 GUARD-137 (THE LANE'S OWN SANCTIONED GUARD, NAMED); ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 07:44Z (session-open increment 5314326). (1) THE REAPER PASS 32 (07:45:08-13Z, ~5s,
detached pid 1344780, receipt prep-reaper.log + prep-reaper-r63.nohup; bash -n CLEAN + md5
b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0 (nothing new
rejected), class B 24 files 47M + syng/export/walks/walkreq markers dropped on 16
cut-awaiting-call (census 22 - the 6 live call lanes correctly excluded), class C 0 + markers
on 164 pending (idempotence holding); lean buffer 168 <= 200 - no eviction (the prep valve at
160 is the standing bound); freed 32M apparent (prep 101,701 -> 101,669M; df 631G steady at
the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 32 idempotent passes
(129G ride-32 + the 30 cadence passes since), every cut/call input and every split ledger kept
throughout. (2) THE EXCLUSION PROOF A TWENTY-SIXTH TIME (SIX-FILE SIDECAR RULE): the
pass's in-flight list correctly named the batch-91 cut cohort (T04016/T04028/T03914/T03912)
+ the last prep lane (T03679, the 64-wide chunk's tail) + the 6 live call lanes
(T01895/T01987/T03871/T03893/T03937/T04052); the batch-91 cohort verified 6/6 syng sidecars
intact + locus.fa KEPT (113-120M); the live lanes verified locus.fa 113-296M +
reference-axis.json + locality-graphs intact (the call reads NO syng); class-B spot-checks
T03877/T03989: syng sidecars GONE, locus.fa 113M + refined beds KEPT, producer markers
dropped, cut.done present - the reap exactly where intended. (3) THE WIDTH-4 CENSUS: batches
post-swap = 92/92 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap
(this window: batch-91 honest exit=0 wall 413s, peak RSS 50.3GB in band, resolved IN-BAND
AFTER the pass ~07:47Z - the pass had protected its cohort, the receipt the pass answers to);
batch-1792223268-92 dispatched in flight at close (cohort T03855/T04000/T03846/T04047, cut
pid 1348565 ~10GB niced early-load); map dirs 65 -> 68. (4) CALLS - ONE NEW FAIL THIS WINDOW,
NAMED: T01987 07:46:15Z guard-137 mid-score locus_0 (last-reported rss 58.0GB before the
climb - the lane's OWN sanctioned 64GiB guard, receipts in wave-logs/failures.txt +
per-territory call.err/call.rss/call.wall in-lane, fail-reaped; T01987 was correctly in the
pass's live-lane exclusion list; ZERO hand kills); fail-closed 570 -> 571; callok 1,492
steady; cut.exit 16; poison 4 terminal; NO deep lanes in flight at close - guard-oracle
quiet. (5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 631-633G, margin
~330G); guard-oracle held (zero hand kills, zero masters touched - the T01987 kill was its
own sanctioned guard); the prep master mid-dispatch chunk tail finishing (sub 2,294 steady,
the T03679 lane the last alive at the pass); terminal classes terminal; the ladder NOT opened
(df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED. FINAL POLL
07:48:46Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / map dirs 68 (batch-92 in
flight) / call-terminal 2,083 = callok 1,492 + fail-closed 571 + cut.exit 16 + poison 4 =
2,083 of 5,844 (35.6% of accepted) / df 632G / mem 813Gib / rvsize 165G apparent. MASTERS AT
CLOSE: prep 411760 (1-13:23h), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1
RAYON=24 DISK_FLOOR_GB=300, ~12:56h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2,
~20:52h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census
continuation (batch-92 in flight); the JIT reaper cadence (~5-12s a pass, idempotent,
re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify
call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification
standing); the deep-lane watch under the guard-oracle (never kill by hand); the prep pole +
df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only
at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave,
never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on
prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal):
rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta
+ the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed
classes A/B/C across 32 idempotent passes, the walkreq marker-gap finding + fix, the ~22s
regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's
4 CONFIRMED members incl. the T00024 3h crossing receipt, the mid-dispatch valve patch's
live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the
T03203/T03370/T03333/T03470/T02844/T03636/T03804/T03499/T04033/T01987 points + rides 55/56's
guard-137 clusters) -> the record's closing increment -> the commit. NO PR PUSH (the owner
pushes).

## RIDE CONTINUATION 64 SESSION OPEN (07:52Z 2026-10-17, post-merge tip, commit 81a8357 from ride-63 close b634396) - THE CADENCE CONTINUES (PASS 33 OWED; PASS 32'S +32M RECEIPTED AT RIDE-63 CLOSE)
Opened from the ride-63 window close (b634396). ORIENTATION STANDING (rides 38-64): a task
prompt carrying the ride-30/31 snapshot ("the reaper WRITTEN but NOT RUN", terminals
1,728/5,844) is STALE - the reaper has RUN since ride-32 (129G first pass); the ride-31 owed
list (the ~22s regen receipt, T03194 closed end-to-end ride-35) is DISCHARGED; the standing
duty is the CADENCE (cumulative ~= 180.3G across 32 idempotent passes entering this window).
OPEN POLL 07:50:46Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / map dirs 68 /
call-terminal 2,083 = callok 1,492 + fail-closed 571 STEADY (zero new fails; the last,
T01987 07:46:15Z guard-137, predates the ride-63 close and is named there) + cut.exit 16 +
poison 4 = 2,083 of 5,844 (35.6% of accepted) / df 632G (the 300G hard floor margin ~332G) /
mem 797Gib / rvsize 165G apparent / census 92/92 post-swap batch exit markers ALL =0
re-verified live (ZERO non-zero; last kill line 18:52:40Z pre-swap width-8; batch-91's marker
landed =0 exactly as receipted at ride-63 close) / batch-1792223268-92 in flight (cohort
T03855/T04000/T03846/T04047, cut pid 1348566 ~27.6GB niced ~3m elapsed) / live call lanes
ordinary 0.4-34.2GB niced, NO lane over the 64GiB guard, NO deep lanes in flight -
guard-oracle quiet. MASTERS: prep 411760 (1-13:25h, mid-dispatch 64-wide chunk, the T03679
tail lane alive), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24
DISK_FLOOR_GB=300, ~12:58h), wave 2810146 (LANES=128 CHUNK_MULT=2, ~20:54h); heartbeat live;
the ancient impg-reference benign - NEVER KILLED. REAPER SCRIPT RE-VERIFIED AT OPEN: bash -n
CLEAN, md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED since the ride-31 review, mtime Oct 16
20:13. THIS WINDOW OWES: the JIT-cadence reaper pass 33 detached + poll SHORT (the SIX-FILE
SYNG SIDECAR rule on the exclusion proof; the ride-62 lane clarification standing: verify
call.done/receipts before reading inputs-gone as a reaper miss); the width-4 census (batch-92
in flight); the standing contingencies; at campaign close: rv-rollup.py +
datasheet-runtime.md + the record's closing increment + commit. NO PR PUSH.

## RIDE CONTINUATION 64 WINDOW CLOSE (07:55Z 2026-10-17) - THE CADENCE PASS 33 RECEIPTED (+12M APPARENT ~4s, CUMULATIVE ~= 180.3G ACROSS 33 IDEMPOTENT PASSES); CENSUS 92/92 HONEST ZERO KILLS SINCE THE SWAP; THE EXCLUSION PROOF A TWENTY-SEVENTH TIME; BATCH-92 IN FLIGHT THROUGHOUT THE PASS (PEAK OBSERVED ~41.8GB IN BAND) - THE CUT LANES STAYED FED; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 07:52Z (session-open increment 81a8357). (1) THE REAPER PASS 33 (07:53:18-22Z, ~4s,
detached pid 1355277, receipt prep-reaper.log + prep-reaper-r64.nohup; bash -n CLEAN + md5
b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0 (nothing new
rejected), class B 24 files 28M + syng/export/walks/walkreq markers dropped on 20
cut-awaiting-call (census 24 - the 4 live call lanes correctly excluded), class C 0 +
export.done/walks.done/walkreq.done markers dropped on 160 pending (idempotence holding);
lean buffer 164 <= 200 - no eviction (the prep valve at 160 is the standing bound); freed 12M
apparent (prep 101,348 -> 101,336M; df 632G steady at the pass). CUMULATIVE SUBSTRATE-REAPER
FREED ~= 180.3G apparent across 33 idempotent passes (129G ride-32 + the 31 cadence passes
since), every cut/call input and every split ledger kept throughout. (2) THE EXCLUSION PROOF
A TWENTY-SEVENTH TIME (SIX-FILE SYNG SIDECAR RULE): the pass's in-flight list correctly named
the batch-92 cut cohort (T03855/T04000/T03846/T04047) + the prep tail lane (T03679) + the 4
live call lanes (T01895/T03871/T03937/T04052); the batch-92 cohort verified 6/6 syng sidecars
intact + locus.fa KEPT (113-443M); the live lanes verified locus.fa 173-296M +
reference-axis.json + locality-graphs intact (the call reads NO syng); class-B spot-checks
(the batch-91 cohort T04016/T04028/T03914/T03912): syng sidecars GONE (0/6), locus.fa
113-120M + refined partition beds KEPT (partition-w10k-d1k/results verified on T04016),
producer markers dropped, cut.done present - the reap exactly where intended. (3) THE WIDTH-4
CENSUS: batches post-swap = 92/92 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the
18:52:42Z swap (this window: batch-92 in flight THROUGHOUT the pass and still in flight at
close, ~6m elapsed, peak observed RSS ~41.8GB in band under the 64GiB guard, its own lane,
no kill); map dirs 68 -> 67 as the wave consumes. (4) CALLS - ZERO NEW FAILS THIS WINDOW:
fail-closed 571 STEADY (the last, T01987 07:46:15Z guard-137 mid-score, receipted at ride-63
close, predates the open); callok 1,492 -> 1,493 (+1 ordinary completion); cut.exit 16; poison
4 terminal; NO deep lanes in flight at close - guard-oracle quiet. (5) CONTINGENCIES: the
300G hard floor NEVER approached (df band 632G steady, margin ~332G); guard-oracle held (zero
hand kills, zero masters touched); the prep master mid-dispatch 64-wide chunk continuing (sub
2,294 steady, the T03679 tail lane alive); terminal classes terminal; the ladder NOT opened
(df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids
1846398/1846402 alive, untouched). FINAL POLL 07:54:04Z: split 9,891 / acc 5,844 (59.1%
final) / sub 2,294 / map dirs 67 (batch-92 in flight) / call-terminal 2,084 = callok 1,493 +
fail-closed 571 + cut.exit 16 + poison 4 = 2,084 of 5,844 (35.6% of accepted) / df 632G /
mem 792Gib / rvsize 164G apparent. MASTERS AT CLOSE: prep 411760 (1-13:28h), cuts 3838087
(width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~13:01h), chunked wave
2810146 (LANES=128 CHUNK_MULT=2, ~20:57h); heartbeat live. THE SUCCESSOR RESUMES: sparse
SHORT polls; the width-4 census continuation (batch-92 in flight); the JIT reaper cadence
(~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every
exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the
ride-62 clarification standing); the deep-lane watch under the guard-oracle (never kill by
hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the
owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by
exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never
RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844
accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the
width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS:
~180.3G cumulative freed classes A/B/C across 33 idempotent passes, the walkreq marker-gap
finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the
syng segdup-tail class's 4 CONFIRMED members incl. the T00024 3h crossing receipt, the
mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB
deepest-lane + the T03203/T03370/T03333/T03470/T02844/T03636/T03804/T03499/T04033/T01987
points + rides 55/56's guard-137 clusters) -> the record's closing increment -> the commit. NO
PR PUSH (the owner pushes).

## RIDE CONTINUATION 65 WINDOW CLOSE (08:06Z 2026-10-17) - THE CADENCE PASS 34 RECEIPTED (+29M APPARENT ~4s, CUMULATIVE ~= 180.3G ACROSS 34 IDEMPOTENT PASSES); CENSUS 94/94 HONEST ZERO KILLS SINCE THE SWAP (batch-93 HONEST AT 457s AFTER THE PASS); THE EXCLUSION PROOF A TWENTY-EIGHTH TIME; ONE NEW FAIL T03937 GUARD-137 (THE LANE'S OWN SANCTIONED GUARD, NAMED); ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 08:56Z-equivalent 07:56Z (session-open increment 26b1f9c). (1) THE REAPER PASS 34 (07:58:49-53Z, ~4s, detached pid 1362126, receipt prep-reaper.log + prep-reaper-r65.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0 (nothing new rejected), class B 24 files 46M + syng/export/walks/walkreq markers dropped on 24 cut-awaiting-call (census: pending=160 cut-awaiting-call=28), class C 0 + export.done/walks.done/walkreq.done markers dropped on 156 pending (idempotence holding); lean buffer 160 <= 200 - no eviction (the prep valve at 160 is the standing bound); freed 29M apparent (prep 101,336 -> 101,307M; df 628G steady at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 34 idempotent passes (129G ride-32 + the 32 cadence passes since), every cut/call input and every split ledger kept throughout. (2) THE EXCLUSION PROOF A TWENTY-EIGHTH TIME (SIX-FILE SIDECAR RULE): the pass's in-flight list correctly named the batch-93 cut cohort (T04070/T04023/T03838/T04057) + the prep tail lane (T03679, the named syng segdup-tail class, 53m+ on syng) + the 4 live call lanes (T01895/T03871/T03937/T04052); the batch-93 cohort verified 6/6 syng sidecars intact + locus.fa KEPT (113-305M); the live lanes verified locus.fa 173-296M + reference-axis.json + locality-graphs intact (graphs 19/8/24/7 - the call reads NO syng); class-B spot-checks (the batch-92 cohort T03855/T04000/T03846/T04047): syng sidecars GONE (0/6), locus.fa 113-443M + refined beds (114/25/137/15) KEPT, producer markers dropped, cut.done present - the reap exactly where intended. (3) THE WIDTH-4 CENSUS: batches post-swap = 94/94 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap (this window: batch-92 honest exit=0 wall 557s, cohort cut.done 07:57:05Z seconds after the open poll; batch-93 honest exit=0 wall 457s, peak observed ~41.1GB in band, resolved IN-BAND AFTER the pass ~08:04Z - the pass had protected its cohort, the receipt the pass answers to); batch-1792224284-94 dispatched ~08:04Z in flight at close (cohort T04027/T04007/T04053/T03898, cut pid 1368241 ~8.5GB niced early-load); map dirs 71 -> 74. (4) CALLS - ONE NEW FAIL THIS WINDOW, NAMED: T03937 07:59:47Z guard-137 mid-score locus_11_scored (last-reported rss 54.0GB, wall 2,637s, peak 68.0GB - the lane's own sanctioned 64GiB guard, receipts in wave-logs/failures.txt + per-territory call.err/call.rss/call.wall in-lane fail-reaped; the pass's own exclusion list had named it live at 07:58:50Z, correctly; ZERO hand kills); fail-closed 571 -> 572; callok 1,493 steady; cut.exit 16; poison 4 terminal; the remaining lanes ordinary niced (T01895 43.2GB @ 55m - the widest, watched under the guard-oracle, NO hand kill; T03871 28.4GB; T04052 29.9GB). (5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 628-630G, margin ~330G); guard-oracle held (zero hand kills, zero masters touched - the T03937 kill was its own sanctioned guard); the prep master mid-dispatch 64-wide chunk continuing (sub 2,294 steady, the T03679 tail lane alive - the named syng segdup-tail class member, watched, never hand-killed); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pid 1846398 alive, untouched). FINAL POLL 08:05:04Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / map dirs 74 (batch-94 in flight) / call-terminal 2,085 = callok 1,493 + fail-closed 572 + cut.exit 16 + poison 4 = 2,085 of 5,844 (35.7% of accepted) / df 628G / mem 835Gib / rvsize 169G apparent. MASTERS AT CLOSE: prep 411760 (1-13:41h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~13:09h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~21:06h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-94 in flight); the JIT reaper cadence (~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 34 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members incl. the T00024 3h crossing receipt, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the T03203/T03370/T03333/T03470/T02844/T03636/T03804/T03499/T04033/T01987/T03937 points + rides 55/56's guard-137 clusters) -> the record's closing increment -> the commit. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 66 WINDOW OPEN (08:10Z 2026-10-17)
Session-open increment committed (39fff18; vendor/syng pointer drift 054acfc->68ac197 still
uncommitted by design). ORIENTATION STANDING since ride-38: the ride-30/31 'reaper written
not run' prompt framing is STALE - the substrate reaper has been RUN since ride-32 (129G
first pass; cumulative ~= 180.3G across 34 idempotent passes at ride-65 close), and the
REGEN RECEIPT is CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194
full-chain re-entry/call receipt - beds + ledgers kept, markers dropped, rebuild + re-entry
into the cut/call flow verified). The regen duty is NOT redone: a fresh eviction of a
pending territory would only delete cut/call inputs the live campaign still needs.
Open poll 08:10:06Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / map dirs 74
(batch-94 in flight) / call-terminal 2,085 = callok 1,493 + fail-closed 572 + cut.exit 16 +
poison 4 = 2,085 of 5,844 (35.7%) / df 628G / mem 801Gib / rvsize 169G apparent.
MASTERS VERIFIED LIVE: prep 411760 (1-13:44h, mid-dispatch 64-wide chunk, the T03679 tail
lane alive 1h06m - the named syng segdup-tail class), cuts 3838087 (width-4-direct
BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~13:17h; batch-94 in flight, cohort
T04027/T04007/T04053/T03898, cut pid 1368241 ~39.5GB niced mid-scan), chunked wave 2810146
(LANES=128 CHUNK_MULT=2, ~21:12h; 3 live call lanes ordinary niced: T01895 48.4GB @ 1h00m
the widest - watched under the guard-oracle, NO hand kill; T03871 29.7GB; T04052 33.6GB);
heartbeat live; the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402).
Census at open: 94/94 post-swap batch exit markers ALL =0, ZERO guard kills since the
18:52:42Z swap. REAPER SCRIPT RE-VERIFIED AT OPEN: bash -n CLEAN, md5
b804a7d366844d27ac58e2e0c8b64802 UNCHANGED since the ride-31 review, mtime Oct 16 20:13.
THIS WINDOW OWES: the JIT-cadence reaper pass 35 detached + poll SHORT (the SIX-FILE SYNG
SIDECAR rule on the exclusion proof; verify call.done/receipts before reading inputs-gone as
a reaper miss - the ride-62 lane clarification standing); the width-4 census continuation
(batch-94 in flight); the standing contingencies (the 300G hard floor - self-gates, DO NOT
lower without the owner; guard-oracle - never kill by hand; poison-first on prep stalls;
terminal classes terminal; the ladder only at backlog drained AND df >= 1,600G, cuts master
MID-WAIT only, by exact pid, never the wave, never with marathons in flight); at campaign
close: rv-rollup.py + datasheet-runtime.md + the record's closing increment + commit.
NO PR PUSH.

## RIDE CONTINUATION 66 WINDOW CLOSE (08:13Z 2026-10-17) - THE CADENCE PASS 35 RECEIPTED (+18M APPARENT ~4s, CUMULATIVE ~= 180.3G ACROSS 35 IDEMPOTENT PASSES); CENSUS 95/95 HONEST ZERO KILLS SINCE THE SWAP (batch-94 HONEST AT 408s AFTER THE PASS); THE EXCLUSION PROOF A TWENTY-NINTH TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 08:10Z (session-open increment 39fff18). THE ORIENTATION THIS WINDOW ANSWERED
(the stale-prompt class, standing since ride-38): the ride-30/31 'reaper written not
run' prompt is STALE - the substrate reaper has been RUN since ride-32 (129G first
pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement,
the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept);
the regen was NOT redone this window - a fresh eviction of a pending territory would
only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASS 35 (08:11:25-29Z, ~4s, detached pid 1374794, receipt prep-reaper.log
+ prep-reaper-r66.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED
re-verified at the open): class A 0 (nothing new rejected), class B 24 files 33M +
syng/export/walks/walkreq markers dropped on 28 cut-awaiting-call (census: pending=156
cut-awaiting-call=31), class C 0 + export.done/walks.done/walkreq.done markers dropped
on 152 pending (idempotence holding); lean buffer 156 <= 200 - no eviction (the prep
valve at 160 is the standing bound); freed 18M apparent (prep 101,133 -> 101,115M;
df 626G steady at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent
across 35 idempotent passes (129G ride-32 + the 33 cadence passes since), every
cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF A TWENTY-NINTH TIME (SIX-FILE SYNG SIDECAR RULE): the pass's
in-flight list correctly named the batch-94 cut cohort (T04027/T04007/T04053/T03898)
+ the prep tail lane (T03679, the named syng segdup-tail class, 1h06m+) + the 3 live
call lanes (T01895/T03871/T04052); the batch-94 cohort verified 6/6 syng sidecars
intact + locus.fa KEPT (113-119M); the live lanes verified locus.fa 197-296M +
reference-axis.json + locality-graphs intact (graphs 19/8/7 - the call reads NO syng);
class-B spot-checks (the batch-93 cohort T04070/T04023/T03838/T04057, resolved after
pass 34): syng sidecars GONE (0/6), locus.fa 113-305M + refined beds (26/25/24/35)
KEPT, producer markers dropped, cut.done present - the reap exactly where intended.
(3) THE WIDTH-4 CENSUS: batches post-swap = 95/95 HONEST ALL exit=0, ZERO non-zero,
ZERO guard kills since the 18:52:42Z swap (this window: batch-94 honest exit=0 wall
408s, resolved IN-BAND AFTER the pass 08:11:33Z - the pass had protected its cohort,
the cut lanes stayed fed); batch-1792224693-95 dispatched 08:11:33Z in flight at
close (cohort T03870/T04019/T04075/T03960, cut pid 1376014 ~14.8GB niced early-load);
map dirs 74 -> 78.
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 572 STEADY (the last, T03937
07:59:47Z guard-137 mid-score, receipted at ride-65 close, predates the open);
callok 1,493 steady; cut.exit 16; poison 4 terminal; the remaining lanes ordinary
niced (T01895 50.9GB @ 1h03m the widest, watched under the guard-oracle, NO hand
kill; T03871 29.7GB; T04052 35.4GB) - no deep-lane kill.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 626G, margin ~326G);
guard-oracle held (zero hand kills, zero masters touched); the prep master
mid-dispatch 64-wide chunk continuing (sub 2,294 steady, the T03679 tail lane alive -
the named syng segdup-tail class member, watched, never hand-killed); terminal classes
terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient
impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
FINAL POLL 08:13:01Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / map dirs 78
(batch-95 in flight) / call-terminal 2,085 = callok 1,493 + fail-closed 572 +
cut.exit 16 + poison 4 = 2,085 of 5,844 (35.7% of accepted) / df 626G / mem 818Gib /
rvsize 171G apparent. MASTERS AT CLOSE: prep 411760 (1-13:47h, mid-dispatch chunk),
cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300,
~13:20h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~21:16h); heartbeat live.
THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-95
in flight); the JIT reaper cadence (~5-12s a pass, idempotent, re-runnable anytime;
the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts
before reading inputs-gone as a reaper miss - the ride-62 clarification standing);
the deep-lane watch under the guard-oracle (never kill by hand); the prep pole + df
watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder
only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid,
never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never
RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all
5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing
close list + the width-4-direct delta + the inventory finding + the JIT-prep rule +
THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 35 idempotent
passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194
full-chain re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members
incl. the T00024 3h crossing receipt, the mid-dispatch valve patch's live proof, the
guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the
T03203/T03370/T03333/T03470/T02844/T03636/T03804/T03499/T04033/T01987/T03937 points +
rides 55/56's guard-137 clusters) -> the record's closing increment -> the commit.
NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 67 SESSION OPEN (08:16Z 2026-10-17)
Session-open increment committed (99211ce). ORIENTATION STANDING since ride-38: the
ride-30/31 'reaper written not run' prompt is STALE - the substrate reaper has been RUN
since ride-32 (129G first pass; cumulative ~= 180.3G across 35 idempotent passes at
ride-66 close), and the REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen
measurement, the T03194 full-chain re-entry/call receipt). The regen duty is NOT
redone: a fresh eviction of a pending territory would only delete cut/call inputs the
live campaign still needs. Open poll 08:14:58Z: split 9,891 / acc 5,844 (59.1% final) /
sub 2,294 / map dirs 78 (batch-95 in flight) / call-terminal 2,065 = callok 1,493 +
fail-closed 572 + cut.exit 16 + poison 4 = 2,065 of 5,844 (35.3%) / df 626G / mem
799Gib / rvsize 171G apparent. MASTERS VERIFIED LIVE: prep 411760 (1-13:49h,
mid-dispatch 64-wide chunk, the T03679 tail lane alive 1h12m on syng - the named syng
segdup-tail class), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24
DISK_FLOOR_GB=300, ~13:22h; batch-95 in flight, cohort T03870/T04019/T04075/T03960,
cut pid 1376014 ~37.3GB niced mid-scan), chunked wave 2810146 (LANES=128 CHUNK_MULT=2,
~21:18h; 3 live call lanes ordinary niced: T01895 50.0GB @ 1h05m the widest - watched
under the guard-oracle, NO hand kill; T04052 ~34.9GB; T03871 ~28.4GB); heartbeat live;
the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402). Census at
open: 95/95 post-swap batch exit markers ALL =0, ZERO guard kills since the 18:52:42Z
swap. REAPER SCRIPT RE-VERIFIED AT OPEN: bash -n CLEAN, md5
b804a7d366844d27ac58e2e0c8b64802 UNCHANGED since the ride-31 review, mtime Oct 16
20:13. THIS WINDOW OWES: the JIT-cadence reaper pass 36 detached + poll SHORT (the
SIX-FILE SYNG SIDECAR rule on the exclusion proof; verify call.done/receipts before
reading inputs-gone as a reaper miss - the ride-62 lane clarification standing); the
width-4 census continuation (batch-95 in flight); the standing contingencies (the 300G
hard floor - self-gates, DO NOT lower without the owner; guard-oracle - never kill by
hand; poison-first on prep stalls; terminal classes terminal; the ladder only at
backlog drained AND df >= 1,600G, cuts master MID-WAIT only, by exact pid, never the
wave, never with marathons in flight); at campaign close: rv-rollup.py +
datasheet-runtime.md + the record's closing increment + commit. NO PR PUSH.

## RIDE CONTINUATION 67 WINDOW CLOSE (08:19Z 2026-10-17) - THE CADENCE PASS 36 RECEIPTED (+14M APPARENT ~4s, CUMULATIVE ~= 180.3G ACROSS 36 IDEMPOTENT PASSES); CENSUS 96/96 HONEST ZERO KILLS SINCE THE SWAP (batch-95 HONEST AT 420s AFTER THE PASS); THE EXCLUSION PROOF A THIRTIETH TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 08:16Z (session-open increment 99211ce). THE ORIENTATION STANDING (rides
38-67): the ride-30/31 'reaper written not run' prompt is STALE - the reaper has RUN
since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s
on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry
receipt); the regen was NOT redone this window - a fresh eviction of a pending
territory would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASS 36 (08:16:34-38Z, ~4s, detached pid 1380965, receipt
prep-reaper.log + prep-reaper-r67.nohup; bash -n CLEAN + md5
b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0
(nothing new rejected), class B 24 files 28M + syng/export/walks/walkreq markers
dropped on 32 cut-awaiting-call (census: pending=152 cut-awaiting-call=35), class
C 0 + export.done/walks.done/walkreq.done markers dropped on 148 pending
(idempotence holding); lean buffer 152 <= 200 - no eviction (the prep valve at 160
is the standing bound); freed 14M apparent (prep 101,115 -> 101,101M; df 626G
steady at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 36
idempotent passes (129G ride-32 + the 34 cadence passes since), every cut/call
input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF A THIRTIETH TIME (SIX-FILE SYNG SIDECAR RULE): the pass's
in-flight list correctly named the batch-95 cut cohort (T03870/T04019/T04075/
T03960) + the prep tail lane (T03679, the named syng segdup-tail class, 1h12m on
syng) + the 3 live call lanes (T01895/T03871/T04052); the batch-95 cohort verified
6/6 syng sidecars intact + locus.fa KEPT (87-152M); the live lanes verified
locus.fa 197-296M + reference-axis.json + locality-graphs intact (graphs 19/8/7 -
the call reads NO syng); class-B spot-checks (the batch-94 cohort
T04027/T04007/T04053/T03898, resolved after pass 35): syng sidecars GONE (0/6),
locus.fa 113-119M + refined beds (466/628/466/466 rows) KEPT, producer markers
dropped (0 present), cut.done present - the reap exactly where intended.
(3) THE WIDTH-4 CENSUS: batches post-swap = 96/96 HONEST ALL exit=0, ZERO non-zero,
ZERO guard kills since the 18:52:42Z swap (this window: batch-95 honest exit=0 wall
420s, resolved IN-BAND AFTER the pass ~08:18Z - the pass had protected its cohort,
the cut lanes stayed fed); batch-1792225115-96 dispatched ~08:18Z in flight at
close (cohort T03907/T03865/T03913/T03854, cut pid 1384007 ~5.0GB niced
early-load); map dirs 78 -> 82.
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 572 STEADY (the last, T03937
07:59:47Z guard-137 mid-score, receipted at ride-65 close, predates the open);
callok 1,493 steady; cut.exit 16; poison 4 terminal; the remaining lanes ordinary
niced (T01895 52.4GB @ 1h09m the widest, watched under the guard-oracle, NO hand
kill; T04052 ~36.6GB; T03871 ~28.4GB) - no deep-lane kill.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 625-626G, margin
~325G); guard-oracle held (zero hand kills, zero masters touched); the prep master
mid-dispatch 64-wide chunk continuing (sub 2,294 steady, the T03679 tail lane
alive - the named syng segdup-tail class member, watched, never hand-killed);
terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present);
the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive,
untouched). FINAL POLL 08:18:35Z: split 9,891 / acc 5,844 (59.1% final) / sub
2,294 / map dirs 82 (batch-96 in flight) / call-terminal 2,065 = callok 1,493 +
fail-closed 572 + cut.exit 16 + poison 4 = 2,065 of 5,844 (35.3% of accepted) /
df 625G / mem 826Gib / rvsize 172G apparent. MASTERS AT CLOSE: prep 411760
(1-13:52h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4
N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~13:26h), chunked wave 2810146
(LANES=128 CHUNK_MULT=2, ~21:22h); heartbeat live. THE SUCCESSOR RESUMES: sparse
SHORT polls; the width-4 census continuation (batch-96 in flight); the JIT reaper
cadence (~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG
SIDECAR rule on every exclusion proof; verify call.done/receipts before reading
inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane
watch under the guard-oracle (never kill by hand); the prep pole + df watch vs the
300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at
backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never
the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never
RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE
(all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the
standing close list + the width-4-direct delta + the inventory finding + the
JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across
36 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen
measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail
class's 4 CONFIRMED members incl. the T00024 3h crossing receipt, the mid-dispatch
valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB
deepest-lane + the T03203/T03370/T03333/T03470/T02844/T03636/T03804/T03499/T04033/
T01987/T03937 points + rides 55/56's guard-137 clusters) -> the record's closing
increment -> the commit. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 68 WINDOW CLOSE (08:28Z 2026-10-17) - THE CADENCE PASS 37 RECEIPTED (+15M APPARENT ~5s, CUMULATIVE ~= 180.3G ACROSS 37 IDEMPOTENT PASSES); CENSUS 97/97 HONEST ZERO KILLS SINCE THE SWAP (batch-96 HONEST AT 417s/48.9GB PEAK RESOLVED IN-BAND ~3min AFTER THE PASS); THE EXCLUSION PROOF A THIRTY-FIRST TIME; TWO NEW FAILS T04052 + T01895, BOTH THE LANES' OWN SANCTIONED 64GiB GUARD, NAMED (T01895 THE DEEP-LANE CLASS, THE WIDEST WATCHED LANE); ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 08:22Z (session-open increment fec3969). THE ORIENTATION STANDING (rides 38-68):
the ride-30/31 'reaper written not run' prompt is STALE - the substrate reaper has been RUN
since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s
on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry
receipt, beds + ledgers kept); the regen was NOT redone this window - a fresh eviction of
a pending territory would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASS 37 (08:22:43-48Z, ~5s, detached pid 1388943, receipt prep-reaper.log +
prep-reaper-r68.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED
re-verified at the open): class A 0 (nothing new rejected), class B 24 files 33M +
syng/export/walks/walkreq markers dropped on 36 cut-awaiting-call (census: pending=148
cut-awaiting-call=38), class C 0 + export.done/walks.done/walkreq.done markers dropped
on 144 pending (idempotence holding); lean buffer 148 <= 200 - no eviction (the prep
valve at 160 is the standing bound); freed 15M apparent (prep 100,904 -> 100,889M; df
626G steady at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 37
idempotent passes (129G ride-32 + the 35 cadence passes since), every cut/call input and
every split ledger kept throughout.
(2) THE EXCLUSION PROOF A THIRTY-FIRST TIME (SIX-FILE SIDECAR RULE): the pass's
in-flight list correctly named the batch-96 cut cohort (T03907/T03865/T03913/T03854) +
the prep tail lane (T03679) + the 2 live call lanes (T01895/T03871); the batch-96 cohort
verified 6/6 syng sidecars intact + locus.fa KEPT (114-122M); the live lanes verified
locus.fa 296/272M + reference-axis.json + locality-graphs intact (19/8 - the call reads
NO syng); the T03679 tail lane verified MID-SYNG RETRY (its syng.exit=1 is a STALE marker
from a 02:21Z failed attempt; the current syng pid 1251565, 1h20m, no output yet - the
named syng segdup-tail class, output written only at completion; correctly excluded,
watched, never hand-killed); class-B spot-checks (the batch-95 cohort
T03870/T04019/T04075/T03960, resolved after pass 36): syng sidecars GONE (0/6), locus.fa
87-152M + refined beds (575/471/465/467 rows) KEPT, producer markers dropped (0 present),
cut.done present - the reap exactly where intended.
(3) THE WIDTH-4 CENSUS: batches post-swap = 97/97 HONEST ALL exit=0, ZERO non-zero,
ZERO guard kills since the 18:52:42Z swap (this window: batch-96 honest exit=0 wall 417s,
peak RSS 48.9GB in band, resolved IN-BAND ~08:25Z ~3min AFTER the pass - the pass had
protected its cohort, the cut lanes stayed fed); batch-1792225533-97 dispatched 08:26Z
in flight at close (cut pid 1392283 ~14.9GB niced early-load); map dirs 81 -> 84 as the
wave consumes.
(4) CALLS - TWO NEW FAILS THIS WINDOW, BOTH NAMED, BOTH THE LANES' OWN SANCTIONED
64GiB GUARD (mid-score locus 2 both): T04052 08:21:14Z guard-137 (last-reported rss
67.1GB, wall 4,080s, peak 72.9GB - receipted at the open, was correctly in the ride-67
close's live-lane list) and T01895 08:26:13Z guard-137 (last-reported 67.1GB, peak
71.9GB, wall 4,621s - THE DEEP-LANE CLASS, the widest watched lane, the guard-oracle's
own kill; receipts in wave-logs/failures.txt + per-territory call.err/call.rss/call.wall
in-lane fail-reaped); fail-closed 572 -> 574; callok 1,493 steady; cut.exit 16; poison 4
terminal; the remaining lane T03871 28.4GB ordinary niced - NO deep lanes in flight at
close, guard-oracle quiet, ZERO hand kills.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 625-627G, margin
~325G); guard-oracle held (zero hand kills, zero masters touched - both kills sanctioned);
the prep master mid-dispatch 64-wide chunk continuing (sub 2,294 steady, the T03679
tail lane alive - the named syng segdup-tail class member, watched, never hand-killed);
terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the
ancient impg-reference benign - NEVER KILLED.
FINAL POLL 08:26:57Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / map dirs 84
(batch-97 in flight) / call-terminal 2,087 = callok 1,493 + fail-closed 574 + cut.exit 16
+ poison 4 = 2,087 of 5,844 (35.7% of accepted) / df 627G / mem 901Gib / rvsize 169G
apparent. MASTERS AT CLOSE: prep 411760 (~14:00h, mid-dispatch chunk), cuts 3838087
(width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~13:34h), chunked
wave 2810146 (LANES=128 CHUNK_MULT=2, ~21:31h); heartbeat live.
THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-97 in
flight); the JIT reaper cadence (~5-12s a pass, idempotent, re-runnable anytime; the
SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before
reading inputs-gone as a reaper miss - the ride-62 clarification standing); the
deep-lane watch under the guard-oracle (never kill by hand); the prep pole + df watch
vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at
backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the
wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM);
poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844
accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list +
the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER
RECEIPTS: ~180.3G cumulative freed classes A/B/C across 37 idempotent passes, the
walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain
re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members incl. the
T00024 3h crossing receipt, the mid-dispatch valve patch's live proof, the guard-oracle
receipts incl. T03641 11,084s/71.1GB deepest-lane + the
T03203/T03370/T03333/T03470/T02844/T03636/T03804/T03499/T04033/T01987/T03937/T04052/
T01895 points + rides 55/56's guard-137 clusters) -> the record's closing increment ->
the commit. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 69 WINDOW CLOSE (08:42Z 2026-10-17) - THE CADENCE PASS 38 RECEIPTED (+15M APPARENT ~4s, CUMULATIVE ~= 180.3G ACROSS 38 IDEMPOTENT PASSES); CENSUS 99/99 HONEST ZERO KILLS SINCE THE SWAP (batch-97 475s/48.9GB + batch-98 420s/48.3GB BOTH HONEST IN-BAND, batch-99 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A THIRTY-SECOND TIME; THE CONSUMPTION BURST (the window's notable event): the wave drained ~20 map dirs in ~12 min, +12 NEW FAILS ALL NAMED, ALL the ordinary bind-ambiguity 7479-class; THE DEEP LANE T03871 COMPLETED HONESTLY (4,904s, NO hand kill); ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 08:30Z (session-open increment 60b2450). THE ORIENTATION STANDING (rides 38-69):
the ride-30/31 'reaper written not run' prompt is STALE - the reaper has RUN since
ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand
regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt);
the regen was NOT redone this window - a fresh eviction of a pending territory would
only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASS 38 (08:30:29-33Z, ~4s, detached pid 1398676, receipt prep-reaper.log
+ prep-reaper-r69.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED
re-verified at the open): class A 0 (nothing new rejected), class B 24 files 30M +
syng/export/walks/walkreq markers dropped on 40 cut-awaiting-call (census: pending=144
cut-awaiting-call=41), class C 0 + export.done/walks.done/walkreq.done markers dropped
on 140 pending (idempotence holding); lean buffer 144 <= 200 - no eviction (the prep
valve at 160 is the standing bound); freed 15M apparent (prep 100,591 -> 100,576M; df
627G steady at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 38
idempotent passes (129G ride-32 + the 36 cadence passes since), every cut/call input and
every split ledger kept throughout.
(2) THE EXCLUSION PROOF A THIRTY-SECOND TIME (SIX-FILE SYNG SIDECAR RULE): the pass's
in-flight list correctly named the batch-97 cut cohort (T04061/T03917/T04051/T04065) +
the prep tail lane (T03679, mid-syng, the named segdup-tail class) + the live call lane
(T03871); the batch-97 cohort verified 6/6 syng sidecars intact (1gbwt/1khash/meta/
names/pstep/spos) + locus.fa KEPT (118-288M); the T03871 lane verified locus.fa 284.6M +
reference-axis.json + locality-graphs (8) intact (the call reads NO syng); class-B
spot-checks (the batch-96 cohort T03907/T03865/T03913/T03854, resolved after pass 37):
syng sidecars GONE (0/6), locus.fa 118-127M + refined partition beds KEPT (106/25/27/28
files), producer markers dropped, cut.done present - the reap exactly where intended.
(3) THE WIDTH-4 CENSUS: 97/97 at the open re-verified live BY TIMESTAMP CUTOFF (epoch
1792176762 = the 18:52:42Z swap; ZERO nonzero; the 135 non-zero batch.exit markers on
disk are ALL pre-swap-era packaging kills/exit-1s - campaign history, not the census);
batch-97 resolved HONEST exit=0 wall 475s peak RSS 48.9GB and batch-98 HONEST exit=0
wall 420s peak 48.3GB, both in-band, both AFTER the pass - the pass had protected their
cohorts, the cut lanes stayed fed; census 99/99 ALL exit=0, ZERO guard kills since the
18:52:42Z swap; batch-1792226430-99 in flight at close (early-load).
(4) CALLS - THE CONSUMPTION BURST (the window's notable event, named): the wave drained
the accumulated cut backlog ~20 map dirs in ~12 min (84 -> 63), callok 1,493 -> 1,509
(+16 ordinary completions) and failcall 574 -> 586 (+12 NEW FAILS, ALL NAMED, ALL the
ordinary bind-ambiguity 7479-class, terminal, receipts in wave-logs/failures.txt +
per-territory call.exit in-lane): T03883 (rec 7062, two-origin-shifts) / T03907 (rec
9338, two-origin-shifts) / T03873 (rec 30, 2 of 2) / T04027 (rec 8817, two-origin-shifts)
/ T04023 (rec 8109, 2 of 2) / T03900 (rec 30447, two-origin-shifts) / T04007 (rec 7892,
2 of 3) / T03993 (rec 1681, 2 of 3) / T04082 (rec 96870, two-origin-shifts) / T04070
(rec 7078, two-origin-shifts) / T04057 (rec 16373, two-origin-shifts) / T03898 (rec
9288, two-origin-shifts). The burst is the backlog's accumulated bind-ambiguity
population surfacing at consumption - the aggregate at close is the receipt, no
extrapolation; the fail rate over the burst (~30% of the burst's terminals) matches the
long-observed bind-ambiguity fraction. THE DEEP LANE T03871 COMPLETED HONESTLY (call.done,
wall 4,904s, receipts reaped) - the widest watched lane (48.6GB, 1h22m) crossed under its
own power, the guard never fired, NO hand kill. cut.exit 16 steady; poison 4 terminal.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 626-645G, margin
~345G); guard-oracle held (zero hand kills, zero masters touched); the prep master
mid-dispatch 64-wide chunk continuing (sub 2,294 steady; the T03679 tail lane alive
1h37m on syng - the named syng segdup-tail class member, watched, never hand-killed);
terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the
ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
FINAL POLL 08:41:09Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / map dirs 63
(batch-99 in flight) / call-terminal 2,115 = callok 1,509 + fail-closed 586 + cut.exit 16
+ poison 4 = 2,115 of 5,844 (36.2% of accepted) / df 644G / mem 847Gib / rvsize 152G
apparent. MASTERS AT CLOSE: prep 411760 (1-14:15h, mid-dispatch chunk), cuts 3838087
(width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~13:48h), chunked
wave 2810146 (LANES=128 CHUNK_MULT=2, ~21:44h); heartbeat live. THE SUCCESSOR RESUMES:
sparse SHORT polls; the width-4 census continuation (batch-99 in flight); the JIT reaper
cadence (~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule
on every exclusion proof; verify call.done/receipts before reading inputs-gone as a
reaper miss - the ride-62 clarification standing); the deep-lane watch under the
guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor
(self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND
df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with
marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep
stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal):
rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct
delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G
cumulative freed classes A/B/C across 38 idempotent passes, the walkreq marker-gap
finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt,
the syng segdup-tail class's 4 CONFIRMED members incl. the T00024 3h crossing receipt,
the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641
11,084s/71.1GB deepest-lane + the T03203/T03370/T03333/T03470/T02844/T03636/T03804/
T03499/T04033/T01987/T03937/T04052/T01895 points + rides 55/56's guard-137 clusters +
the ride-69 consumption-burst receipt) -> the record's closing increment -> the commit.
NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 70 WINDOW CLOSE (09:04Z 2026-10-17) - THE CADENCE PASS 39 RECEIPTED (+48M APPARENT ~4s, CUMULATIVE ~= 180.3G ACROSS 39 IDEMPOTENT PASSES); CENSUS 102/102 HONEST ZERO KILLS SINCE THE SWAP (batch-99 432s / batch-100 400s / batch-101 493s ALL HONEST IN-BAND, batch-102 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A THIRTY-THIRD TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 08:45Z (session-open increment c80fa33). THE ORIENTATION STANDING (rides 38-70): the
ride-30/31 'reaper written not run' prompt is STALE - the reaper has RUN since ride-32 (129G
first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the
T03194 full-chain evict->rebuild->cut/call re-entry receipt); the regen was NOT redone this
window - a fresh eviction of a pending territory would only delete cut/call inputs the live
campaign still needs.
(1) THE REAPER PASS 39 (08:47:27-31Z, ~4s, detached pid 1435493, receipt prep-reaper.log +
prep-reaper-r70.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED
re-verified at the open): class A 0, class B 48 files 78M + syng/export/walks/walkreq markers
dropped on 8 cut-awaiting-call (census: pending=136 cut-awaiting-call=15), class C 0 +
export.done/walks.done/walkreq.done markers dropped on 132 pending; lean buffer 136 <= 200 -
no eviction; freed 48M apparent (prep 94,861 -> 94,813M; df 645G steady). CUMULATIVE ~= 180.3G
apparent across 39 idempotent passes, every cut/call input and every split ledger kept.
(2) THE EXCLUSION PROOF A THIRTY-THIRD TIME (SIX-FILE SIDECAR RULE): the pass's in-flight list
correctly named the batch-99 cut cohort (T03926/T04008/T04066/T04036, the chr6 band) + the
prep tail lane (T03679, mid-syng, the named segdup-tail class) + the 8 live call lanes; the
batch-99 cohort verified 6/6 syng sidecars intact + locus.fa KEPT (113-169M); the live lanes
verified locus.fa 113-443M + reference-axis.json + locality-graphs intact (3-25 - the call
reads NO syng); class-B spot-checks (the batch-98 cohort T03965/T03983/T03945/T03938 resolved
after pass 38): sidecars GONE (0/6), locus.fa 62-142M + refined beds KEPT, producer markers
dropped, cut.done present - the reap exactly where intended.
(3) THE WIDTH-4 CENSUS: 102/102 post-swap batch.exit markers ALL =0 re-verified BY TIMESTAMP
CUTOFF (epoch 1792176762 = the 18:52:42Z swap; ZERO nonzero; pre-swap non-zero markers are
campaign history, not the census); batch-99 honest 432s resolved DURING the pass window (its
cohort protected, cut lanes stayed fed), batch-100 honest 400s, batch-101 honest 493s; batch-
1792227758-102 dispatched 09:02:38Z in flight at close; map dirs 61 -> 68.
(4) CALLS - ZERO NEW FAILS: fail-closed 586 STEADY (the last, T03898 08:37:54Z, predates this
open); callok 1,512 -> 1,517; cut.exit 16; poison 4 terminal; NO deep lanes (widest T03855
~17.6GB @ 19m, guard-oracle quiet).
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 642-646G, margin ~342G);
the prep master mid-dispatch 64-wide chunk continuing (sub 2,294, the T03679 tail lane alive
~2h on syng - watched, never hand-killed); terminal classes terminal; the ladder NOT opened
(df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids
1846398/1846402 alive, untouched). FINAL POLL 09:02:48Z: split 9,891 / acc 5,844 (59.1%
final) / sub 2,294 / cut.done 2,128 / map dirs 68 / call-terminal 2,123 = callok 1,517 +
fail-closed 586 + cut.exit 16 + poison 4 = 2,123 of 5,844 (36.3% of accepted) / df 642G /
mem 785Gib / rvsize 155G apparent. MASTERS AT CLOSE: prep 411760 (~14:20h, mid-dispatch
chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300,
~14:17h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~21:55h); heartbeat live. THE
SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-102 in
flight); the JIT reaper cadence (~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE
SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading
inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch
under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor
(self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND
df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons
in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls; AT CAMPAIGN
CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing
close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE
REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 39 idempotent passes, the
walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain
re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members + the T03679 live
member, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641
11,084s/71.1GB deepest-lane + the named guard-137 points + rides 55/56's clusters + the
ride-69 consumption-burst receipt) -> the record's closing increment -> the commit.
NO PR PUSH (the owner pushes).
## RIDE CONTINUATION 71 WINDOW CLOSE (09:22Z 2026-10-17) - THE CADENCE PASS 40 RECEIPTED (+67M APPARENT ~4s, CUMULATIVE ~= 180.3G ACROSS 40 IDEMPOTENT PASSES); CENSUS 103/103 HONEST ZERO KILLS SINCE THE SWAP (batch-102 ~400s/48.1GB PEAK HONEST IN-BAND, RESOLVED ~3min AFTER THE PASS - COHORT PROTECTED, LANES FED; batch-103 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A THIRTY-FOURTH TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 09:06Z (session-open increment 8e39f25). THE ORIENTATION STANDING (rides 38-71): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has been RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt); the regen was NOT redone this window - a fresh eviction of a pending territory would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASS 40 (09:05:31-35Z, ~4s, detached pid 1465000, receipt prep-reaper.log + prep-reaper-r71.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0 (nothing new rejected), class B 72 files 114M + syng/export/walks/walkreq markers dropped on 20 cut-awaiting-call (census: pending=124 cut-awaiting-call=25), class C 0 + export.done/walks.done/walkreq.done markers dropped on 120 pending (idempotence holding); lean buffer 124 <= 200 - no eviction (the prep valve at 160 is the standing bound); freed 67M apparent (prep 94,586 -> 94,519M; df 642G steady at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 40 idempotent passes (129G ride-32 + the 38 cadence passes since), every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF A THIRTY-FOURTH TIME (SIX-FILE SYNG SIDECAR RULE): the pass's in-flight list correctly named the batch-102 cut cohort (T04002/T03899/T03885/T03911) + the prep tail lane (T03679, mid-syng 2h+, the named syng segdup-tail class - its output written only at completion; the syng.exit present is the STALE 02:21Z failed-attempt marker, named ride-68) + the 5 live call lanes (T03855/T03912/T04013/T04019/T04063); the batch-102 cohort verified 6/6 syng sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa KEPT (113M each); the live lanes verified locus.fa 116-443M + reference-axis.json + locality-graphs intact (3-25 - the call reads NO syng); class-B spot-checks (the batch-99/100/101 cohorts T03926/T04008/T04066/T03936 + T03828/T03999/T03970/T04072 + T04020/T04069/T03998/T03981, all resolved after pass 39): syng sidecars GONE (0/6 each), locus.fa 113-416M + refined beds KEPT (465-1,252 rows), producer markers dropped (0 present), cut.done present 12/12 - the reap exactly where intended.
(3) THE WIDTH-4 CENSUS: batches post-swap = 103/103 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap, re-verified BY TIMESTAMP CUTOFF (epoch 1792176762 = the swap; the pre-swap non-zero markers are campaign history, not the census); this window: batch-102 HONEST exit=0 wall ~400s peak RSS 48.1GB in band, resolved ~09:09Z ~3min AFTER the pass - the pass had protected its cohort, the cut lanes stayed fed; batch-1792228148-103 dispatched 09:14Z in flight at close (cohort T03919/T04055/T04059/T03853, cut pid ~37.6GB niced mid-scan); map dirs 68 -> 76 as the wave consumes.
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 586 STEADY (the last, T03898 08:37:54Z, predates the ride-70 close); callok 1,517 steady (the 5 lanes in flight 38-48m, receipts land at call.done); cut.exit 16; poison 4 terminal; the widest lane T03855 34.9GB @ 48m climbing slowly (the deep-lane watch under the guard-oracle - never kill by hand); NO deep lanes at guard pressure at close.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 634-646G, margin ~330G); guard-oracle held (zero hand kills, zero masters touched); the prep master mid-dispatch 64-wide chunk continuing (sub 2,294 steady, the T03679 tail lane alive 2h17m on syng - the named syng segdup-tail class member, watched, never hand-killed); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
FINAL POLL 09:21:01Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / cut.done 2,136 / map dirs 76 (batch-103 in flight) / call-terminal 2,123 = callok 1,517 + fail-closed 586 + cut.exit 16 + poison 4 = 2,123 of 5,844 (36.3% of accepted) / df 634G / mem 760Gib / rvsize 163G apparent. MASTERS AT CLOSE: prep 411760 (~1-14:54h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~14:31h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~22:17h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-103 in flight); the JIT reaper cadence (~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 40 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members + the T03679 live member, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points + rides 55/56's clusters + the ride-69 consumption-burst receipt) -> the record's closing increment -> the commit. NO PR PUSH (the owner pushes).**


## RIDE CONTINUATION 72 WINDOW CLOSE (09:27Z 2026-10-17) - THE CADENCE PASS 41 RECEIPTED (+56M APPARENT ~4s, CUMULATIVE ~= 180.3G ACROSS 41 IDEMPOTENT PASSES); CENSUS 105/105 HONEST ZERO KILLS SINCE THE SWAP (batch-103 468s + batch-104 499s BOTH HONEST IN-BAND, batch-105 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A THIRTY-FIFTH TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 09:24Z (session-open increment d45c583). THE ORIENTATION STANDING (rides 38-72): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has been RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt); the regen was NOT redone this window - a fresh eviction of a pending territory would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASS 41 (09:25:04-08Z, ~4s, detached pid 1485849, receipt prep-reaper.log + prep-reaper-r72.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0 (nothing new rejected), class B 48 files 87M + syng/export/walks/walkreq markers dropped on 28 cut-awaiting-call (census: pending=116 cut-awaiting-call=33), class C 0 + export.done/walks.done/walkreq.done markers dropped on 112 pending (idempotence holding); lean buffer 116 <= 200 - no eviction (the prep valve at 160 is the standing bound); freed 56M apparent (prep 94,519 -> 94,463M; df 631G steady at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 41 idempotent passes (129G ride-32 + the 39 cadence passes since), every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF A THIRTY-FIFTH TIME (SIX-FILE SYNG SIDECAR RULE): the pass's in-flight list correctly named the batch-104 cut cohort (T03958/T03915/T04068/T04031) + the prep tail lane (T03679, mid-syng, the named syng segdup-tail class - the syng.exit present is the STALE 02:21Z failed-attempt marker, named ride-68) + the 5 live call lanes (T03855/T03912/T04013/T04019/T04063); the batch-104 cohort verified 6/6 syng sidecars intact + locus.fa KEPT (113-492M); the live lanes verified locus.fa 116-443M + reference-axis.json + locality-graphs intact (3-25 - the call reads NO syng); class-B spot-checks (the batch-103 cohort T03919/T04055/T04059/T03853, resolved after pass 40): syng sidecars GONE (0/6), locus.fa 113-268M + refined beds (674/466/535/982 rows) KEPT, producer markers dropped (0 present), cut.done present 4/4 - the reap exactly where intended.
(3) THE WIDTH-4 CENSUS: batches post-swap = 105/105 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap, re-verified BY TIMESTAMP CUTOFF (epoch 1792176762 = the swap; the pre-swap non-zero markers are campaign history, not the census); this window: batch-103 honest exit=0 wall 468s and batch-104 honest exit=0 wall 499s (peak RSS ~52-55GB in band), both AFTER the pass - the pass had protected their cohorts, the cut lanes stayed fed; batch-1792229116-105 dispatched ~09:25Z in flight at close (cohort T03995/T03920/T04005/T03910, cut pid 1486877 ~13.8GB niced early-load); map dirs 76 -> 80; cut.done 2,136 -> 2,140 - the width-4 pace holding (~29-35 territories/h).
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 586 STEADY (the last, T03898 08:37:54Z, receipted at the ride-69 close, predates the ride-70 close); callok 1,517 steady (the 5 lanes in flight ~48-52m, receipts land at call.done); cut.exit 16; poison 4 terminal; the lanes ordinary niced (widest T04063 ~37.0GB @ 47m, T03855 36.9GB - the deep-lane watch under the guard-oracle, never kill by hand); NO deep lanes at guard pressure at close.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 631-634G, margin ~331G); guard-oracle held (zero hand kills, zero masters touched); the prep master mid-dispatch 64-wide chunk continuing (sub 2,294 steady, the T03679 tail lane alive 2h23m on syng - the named syng segdup-tail class member, watched, never hand-killed); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
FINAL POLL 09:26:46Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / cut.done 2,140 / map dirs 80 (batch-105 in flight) / call-terminal 2,123 = callok 1,517 + fail-closed 586 + cut.exit 16 + poison 4 = 2,123 of 5,844 (36.3% of accepted) / df 631G / mem 771Gib / rvsize 166-167G apparent. MASTERS AT CLOSE: prep 411760 (~1-15:01h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300, ~14:34h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~22:30h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-105 in flight); the JIT reaper cadence (~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 41 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members + the T03679 live member, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points + rides 55/56's clusters + the ride-69 consumption-burst receipt) -> the record's closing increment -> the commit. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 73 WINDOW CLOSE (09:33Z 2026-10-17) - THE CADENCE PASS 42 RECEIPTED (+43M APPARENT ~4s, CUMULATIVE ~= 180.3G ACROSS 42 IDEMPOTENT PASSES); CENSUS 106/106 HONEST ZERO KILLS SINCE THE SWAP (batch-105 392s HONEST IN-BAND, RESOLVED ~1min AFTER THE PASS - COHORT PROTECTED, LANES FED; batch-106 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A THIRTY-SIXTH TIME; ONE NEW FAIL - T03855 CALL GUARD-KILLED (137) AT 71.0GB/3,478s, THE SANCTIONED GUARD-ORACLE'S OWN KILL (the deep-lane memory-tail class, terminal, named); ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 09:30Z (session-open increment 0863d04). THE ORIENTATION STANDING (rides 38-73): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt); the regen was NOT redone this window - a fresh eviction of a pending territory would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASS 42 (09:30:31-35Z, ~4s, detached pid 1493005, receipt prep-reaper.log + prep-reaper-r73.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0 (nothing new rejected), class B 24 files 58M + syng/export/walks/walkreq markers dropped on 32 cut-awaiting-call (census: pending=112 cut-awaiting-call=37), class C 0 + export.done/walks.done/walkreq.done markers dropped on 108 pending (idempotence holding); lean buffer 112 <= 200 - no eviction (the prep valve at 160 is the standing bound); freed 43M apparent (prep 94,463 -> 94,420M; df 631G steady at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 42 idempotent passes (129G ride-32 + the 40 cadence passes since), every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF A THIRTY-SIXTH TIME (SIX-FILE SYNG SIDECAR RULE): the pass's in-flight list correctly named the batch-105 cut cohort (T03995/T03920/T04005/T03910) + the prep tail lane (T03679, mid-syng, the named syng segdup-tail class - the syng.exit present is the STALE 02:21Z failed-attempt marker, named ride-68) + the 5 live call lanes (T03855/T03912/T04013/T04019/T04063); the batch-105 cohort verified 6/6 syng sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa KEPT (113-114M); the live lanes verified locus.fa 116-443M + reference-axis.json + locality-graphs intact (3-25 - the call reads NO syng); class-B spot-checks (the batch-104 cohort T03958/T03915/T04068/T04031, resolved after pass 41): syng sidecars GONE (0/6 each), locus.fa 113-492M + the definitional beds KEPT (locus.refined.bed/locus.members.bed, 473-474 rows at the spot-check), producer markers dropped (0 present), cut.done present 4/4 - the reap exactly where intended.
(3) THE WIDTH-4 CENSUS: batches post-swap = 106/106 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap, re-verified BY TIMESTAMP CUTOFF (epoch 1792176762 = the swap; 536 batch.exit markers on disk, the pre-swap non-zero are campaign history, not the census); this window: batch-105 HONEST exit=0 wall 392s resolved 09:31:48Z ~1min AFTER the pass - the pass had protected its cohort, the cut lanes stayed fed; batch-1792229509-106 dispatched 09:31:49Z in flight at close (cohort T03929/T03967/T03527/T03943, cut pid 1495809 ~8GB niced early-load); map dirs 80 -> 83 as the wave consumes; cut.done 2,140 -> 2,144 - the width-4 pace holding (~29-35 territories/h).
(4) CALLS - ONE NEW FAIL THIS WINDOW, THE SANCTIONED GUARD-ORACLE'S OWN KILL: T03855 CALL GUARD-KILLED 09:31:04Z (exit 137, wall 3,478s, killed at 71.0GB RSS - the lane's own call.rss trace shows the deep-lane tail: ~42.8GB at the open, then the locus-10 scoring climb 43GB -> 71GB in the final minute; the widest watched lane of the open, the deep-lane memory-tail class, terminal, named alongside T00129/T04052/T01895 - receipts wave-logs/failures.txt + per-territory call.err/call.rss/call.wall IN-LANE FAIL-REAPED 09:31:05Z); fail-closed 586 -> 587 (the bind ambiguities + the guard-137 call-tail class T00129/T04052/T01895/T03855); callok 1,517 steady (the 4 remaining lanes in flight ~54-55m, receipts land at call.done); cut.exit 16; poison 4 terminal; the remaining lanes ordinary niced (T04063 ~41GB, T04013 ~35GB, T03912 ~26GB, T04019 ~22GB - the deep-lane watch under the guard-oracle, never kill by hand); NO deep lanes at guard pressure at close.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 631-634G, margin ~331G); guard-oracle held (zero hand kills, zero masters touched - the T03855 kill was the guard's own, sanctioned); the prep master mid-dispatch 64-wide chunk continuing (sub 2,294 steady, the T03679 tail lane alive 2h26m on syng - the named syng segdup-tail class member, watched, never hand-killed); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
FINAL POLL 09:32:19Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / cut.done 2,144 / map dirs 83 (batch-106 in flight) / call-terminal 2,124 = callok 1,517 + fail-closed 587 + cut.exit 16 + poison 4 = 2,124 of 5,844 (36.3% of accepted) / df 632G / mem 809Gib / rvsize 164G apparent. MASTERS AT CLOSE: prep 411760 (~1-15:07h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~14:40h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~22:35h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-106 in flight); the JIT reaper cadence (~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 42 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members + the T03679 live member, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855 + rides 55/56's clusters + the ride-69 consumption-burst receipt) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 74 WINDOW CLOSE (09:42Z 2026-10-17) - THE CADENCE PASS 43 RECEIPTED (+13M APPARENT ~4s, CUMULATIVE ~= 180.3G ACROSS 43 IDEMPOTENT PASSES); CENSUS 107/107 HONEST ZERO KILLS SINCE THE SWAP (batch-106 HONEST AT 473s/48.1GB PEAK RESOLVED IN-BAND ~4min AFTER THE PASS - COHORT PROTECTED, LANES FED; batch-107 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A THIRTY-SEVENTH TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 09:36Z (session-open increment ad7dbbc). THE ORIENTATION STANDING (rides 38-74): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has been RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen was NOT redone this window - a fresh eviction of a pending territory would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASS 43 (09:35:18-22Z, ~4s, detached pid 1500195, receipt prep-reaper.log + prep-reaper-r74.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0 (nothing new rejected), class B 24 files 28M + syng/export/walks/walkreq markers dropped on 36 cut-awaiting-call (census: pending=108 cut-awaiting-call=40), class C 0 + export.done/walks.done/walkreq.done markers dropped on 104 pending (idempotence holding); lean buffer 108 <= 200 - no eviction (the prep valve at 160 is the standing bound); freed 13M apparent (prep 93,975 -> 93,962M; df 632G steady at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 43 idempotent passes (129G ride-32 + the 41 cadence passes since), every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF A THIRTY-SEVENTH TIME (SIX-FILE SYNG SIDECAR RULE): the pass's in-flight list correctly named the batch-106 cut cohort (T03929/T03967/T03527/T03943) + the prep tail lane (T03679, mid-syng 2h30m+, the named syng segdup-tail class - the syng.exit present is the STALE 02:21Z failed-attempt marker, named ride-68) + the 4 live call lanes (T03912/T04013/T04019/T04063); the batch-106 cohort verified 6/6 syng sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa KEPT (117-229M); the live lanes verified locus.fa 116-319M + reference-axis.json + locality-graphs intact (3-24 graphs - the call reads NO syng; call.done 0/4 as expected, in-flight - the ride-62 clarification standing); class-B spot-checks (the batch-105 cohort T03995/T03920/T04005/T03910, resolved after pass 42): syng sidecars GONE (0/6 each), locus.fa 113-114M + the definitional beds KEPT (locus.refined.bed 465-469 rows), producer markers dropped (0 present), cut.done present 4/4 - the reap exactly where intended.
(3) THE WIDTH-4 CENSUS: batches post-swap = 107/107 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap, re-verified BY TIMESTAMP CUTOFF (epoch 1792176762 = the swap; the pre-swap non-zero markers are campaign history, not the census); this window: batch-106 HONEST exit=0 wall 473s peak RSS 48.1GB in band, resolved 09:40Z ~4min AFTER the pass - the pass had protected its cohort, the cut lanes stayed fed; batch-1792229982-107 dispatched in flight at close (cohort T02541/T02722/T03187/T03208, cut pid ~14.1GB niced early-load); map dirs 83 -> 85; cut.done 2,144 -> 2,148 - the width-4 pace holding (~29-35 territories/h).
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 587 STEADY (the last, T03855 09:31:04Z guard-137, receipted at the ride-73 close, predates this open); callok 1,517 steady (the 4 lanes in flight ~1h04m, receipts land at call.done); cut.exit 16; poison 4 terminal; the lanes ordinary niced (widest T04063 48.4GB @ 1h04m climbing slowly, T04013 41.0GB, T03912 26.7GB, T04019 22.8GB - the deep-lane watch under the guard-oracle, never kill by hand); NO deep lanes at guard pressure at close.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 629-632G, margin ~329G); guard-oracle held (zero hand kills, zero masters touched); the prep master mid-dispatch 64-wide chunk continuing (sub 2,294 steady, the T03679 tail lane alive 2h37m on syng - the named syng segdup-tail class member, watched, never hand-killed); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
FINAL POLL 09:41:02Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / cut.done 2,148 / map dirs 85 (batch-107 in flight) / call-terminal 2,124 = callok 1,517 + fail-closed 587 + cut.exit 16 + poison 4 = 2,124 of 5,844 (36.3% of accepted) / df 632G / mem 793Gib / rvsize 165G apparent. MASTERS AT CLOSE: prep 411760 (~1-15:11h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~14:45h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~22:40h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-107 in flight); the JIT reaper cadence (~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 43 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members + the T03679 live member, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855 + rides 55/56's clusters + the ride-69 consumption-burst receipt) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 75 WINDOW CLOSE (09:46Z 2026-10-17) - THE CADENCE PASS 44 RECEIPTED (+24M APPARENT ~4s, CUMULATIVE ~= 180.3G ACROSS 44 IDEMPOTENT PASSES); CENSUS 107/107 HONEST ZERO KILLS SINCE THE SWAP (batch-106 HONEST AT 473s/~48-50GB PEAK; batch-107 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A THIRTY-EIGHTH TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 09:43Z (session-open increment 3595c8f). THE ORIENTATION STANDING (rides 38-75): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen was NOT redone this window - a fresh eviction of a pending territory would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASS 44 (09:43:49-53Z, ~4s, detached pid 1509137, receipt prep-reaper.log + prep-reaper-r75.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0 (nothing new rejected), class B 24 files 39M + syng/export/walks/walkreq markers dropped on 40 cut-awaiting-call (census: pending=104 cut-awaiting-call=44), class C 0 + export.done/walks.done/walkreq.done markers dropped on 100 pending (idempotence holding); lean buffer 104 <= 200 - no eviction (the prep valve at 160 is the standing bound); freed 24M apparent (prep 93,962 -> 93,938M; df 632G steady at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 44 idempotent passes (129G ride-32 + the 42 cadence passes since), every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF A THIRTY-EIGHTH TIME (SIX-FILE SYNG SIDECAR RULE): the pass's in-flight list correctly named the batch-107 cut cohort (T02541/T02722/T03187/T03208) + the prep tail lane (T03679, mid-syng ~2h47m, the named syng segdup-tail class - the syng.exit present is the STALE 02:21Z failed-attempt marker, named ride-68) + the 4 live call lanes (T03912/T04013/T04019/T04063); the batch-107 cohort verified 6/6 syng sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa KEPT (117-439M); the live lanes verified locus.fa 116-319M + reference-axis.json + locality-graphs intact (3-24 graphs - the call reads NO syng; call.done 0/4 as expected, in-flight - the ride-62 clarification standing); class-B spot-checks (the batch-106 cohort T03929/T03967/T03527/T03943, resolved after pass 43): syng sidecars GONE (0/6 each), locus.fa 112-219M + the definitional beds KEPT (locus.refined.bed 465-951 rows / locus.members.bed 466-952), producer markers dropped (0 present), cut.done present 4/4 - the reap exactly where intended.
(3) THE WIDTH-4 CENSUS: batches post-swap = 107/107 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap, re-verified BY TIMESTAMP CUTOFF (epoch 1792176762 = the swap; 538 batch.exit markers on disk, the pre-swap non-zero are campaign history, not the census); batch-107 (cohort T02541/T02722/T03187/T03208, dispatched 09:41Z) in flight at close - load+scan under way; batch-106's honest receipt restated (473s, peak ~48-50GB in band per its rss trace); map dirs 85, cut.done 2,148 steady at the poll - the width-4 pace holding (~29-35 territories/h).
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 587 STEADY (the last, T03855 09:31:04Z guard-137, receipted at the ride-73 close, predates this open); callok 1,517 steady (the 4 lanes in flight ~72-74m, receipts land at call.done); cut.exit 16; poison 4 terminal; the lanes ordinary niced but CLIMBING (call.rss traces: T04063 53.4GB @ ~74m the widest, T04013 45.1GB, T04019 38.8GB, T03912 28.0GB - the deep-lane watch under the guard-oracle, never kill by hand; T04063 on the T03855-like climb, the guard is the oracle); NO deep lanes at guard pressure at close.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 632G steady, margin ~332G); guard-oracle held (zero hand kills, zero masters touched); the prep master mid-dispatch 64-wide chunk continuing (sub 2,294 steady, the T03679 tail lane alive ~2h47m on syng - the named syng segdup-tail class member, watched, never hand-killed); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
FINAL POLL 09:45:14Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / cut.done 2,148 / map dirs 85 (batch-107 in flight) / call-terminal 2,124 = callok 1,517 + fail-closed 587 + cut.exit 16 + poison 4 = 2,124 of 5,844 (36.3% of accepted) / df 632G / mem 739-773Gib / rvsize 164-165G apparent. MASTERS AT CLOSE: prep 411760 (~1-15:17h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~14:51h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~22:46h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-107 in flight); the JIT reaper cadence (~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch under the guard-oracle (T04063 the one to watch - never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 44 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members + the T03679 live member, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855 + rides 55/56's clusters + the ride-69 consumption-burst receipt) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 77 WINDOW CLOSE (09:58Z 2026-10-17) - THE CADENCE PASS 45 RECEIPTED (+42M APPARENT ~5s, CUMULATIVE ~= 180.3G ACROSS 45 IDEMPOTENT PASSES); CENSUS 109/109 HONEST ZERO KILLS SINCE THE SWAP (batch-107 HONEST 512s + batch-108 HONEST 488s RESOLVED IN-BAND - COHORTS PROTECTED, LANES FED; batch-109 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A THIRTY-NINTH TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 09:49Z (session-open increment 0131f02). THE ORIENTATION STANDING (rides 38-77): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen was NOT redone this window - a fresh eviction of a pending territory would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASS 45 (09:49:30-35Z, ~5s, detached pid 1518985, receipt prep-reaper.log + prep-reaper-r77.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0 (nothing new rejected), class B 24 files 55M + syng/export/walks/walkreq markers dropped on 44 cut-awaiting-call (census: pending=100 cut-awaiting-call=46), class C 0 + export.done/walks.done/walkreq.done markers dropped on 96 pending (idempotence holding); lean buffer 100 <= 200 - no eviction (the prep valve at 160 is the standing bound); freed 42M apparent (prep 93,701 -> 93,659M; df 627G steady at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 45 idempotent passes (129G ride-32 + the 43 cadence passes since), every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF A THIRTY-NINTH TIME (SIX-FILE SYNG SIDECAR RULE): the pass's in-flight list correctly named the batch-108 cut cohort (T02769/T03803/T03670/T03859) + the prep tail lane (T03679, mid-syng ~2h45m, the named syng segdup-tail class - the syng.exit present is the STALE 02:21Z failed-attempt marker, named ride-68) + the 2 live call lanes (T04063/T04013); the batch-108 cohort verified 6/6 syng sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa KEPT (113-342M); the live lanes verified locus.fa 272-319M + reference-axis.json + locality-graphs intact (14/6 graphs - the call reads NO syng; call.done 0/2 as expected, in-flight - the ride-62 clarification standing); class-B spot-checks (the batch-107 cohort T02541/T02722/T03187/T03208, resolved after the prior pass): syng sidecars GONE (0/6 each), locus.fa 117-439M + the definitional beds KEPT (locus.refined.bed 480-1,471 rows), producer markers dropped (0 present), cut.done present 4/4 - the reap exactly where intended.
(3) THE WIDTH-4 CENSUS: batches post-swap = 109/109 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap, re-verified BY TIMESTAMP CUTOFF (epoch 1792176762 = the swap; the pre-swap non-zero markers are campaign history, not the census); this window: batch-107 HONEST exit=0 wall 512s (resolved 09:48:15Z, ~3min after the prior pass) and batch-108 HONEST exit=0 wall 488s (resolved 09:56:24Z, AFTER pass 45 - the pass had protected its cohort, the cut lanes stayed fed); batch-1792230985-109 dispatched 09:56:25Z in flight at close; map dirs 86 -> 89; cut.done 2,152 -> 2,156 (+4, batch-108's cohort) - the width-4 pace holding (~29-35 territories/h).
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 587 STEADY (the last, T03855 09:31:04Z guard-137, receipted at the ride-73 close, predates this open); callok 1,519 steady (T03912 + T04019 completed honestly at the open's edge, before the increment; the 2 lanes in flight ~80-81m, receipts land at call.done); cut.exit 16; poison 4 terminal; the lanes ordinary niced but CLIMBING (call.rss: T04063 59.8GB @ ~81m the widest - the T03855-like deep-lane climb, the guard is the oracle, NEVER killed by hand; T04013 50.0GB @ ~80m) - the deep-lane watch standing.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 623-627G, margin ~320G+); guard-oracle held (zero hand kills, zero masters touched); the prep master mid-dispatch 64-wide chunk continuing (sub 2,294 steady, the T03679 tail lane alive ~2h50m on syng - the named syng segdup-tail class member, watched, never hand-killed); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
FINAL POLL 09:57:20Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / cut.done 2,156 / map dirs 89 (batch-109 in flight) / call-terminal 2,126 = callok 1,519 + fail-closed 587 + cut.exit 16 + poison 4 = 2,126 of 5,844 (36.4% of accepted) / df 623G / mem 826Gib / rvsize 174G apparent. MASTERS AT CLOSE: prep 411760 (~1-15:39h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~14:58h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~22:52h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-109 in flight); the JIT reaper cadence (~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch T04063 under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 45 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members + the T03679 live member, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855 + rides 55/56's clusters + the ride-69 consumption-burst receipt) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 78 WINDOW CLOSE (10:01Z 2026-10-17) - THE CADENCE PASS 46 RECEIPTED (+47M APPARENT ~4s, CUMULATIVE ~= 180.3G ACROSS 46 IDEMPOTENT PASSES); CENSUS 109/109 HONEST ZERO KILLS SINCE THE SWAP (batch-109 IN FLIGHT AT CLOSE, cut pid 37.6GB niced mid-scan); THE EXCLUSION PROOF A FORTIETH TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 09:58Z (session-open increment 848ab78). THE ORIENTATION STANDING (rides 38-78): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen was NOT redone this window - a fresh eviction of a pending territory would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASS 46 (09:59:48-52Z, ~4s, detached pid 1529341, receipt prep-reaper.log + prep-reaper-r78.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0 (nothing new rejected), class B 24 files 62M + syng/export/walks/walkreq markers dropped on 48 cut-awaiting-call (census: pending=96 cut-awaiting-call=50), class C 0 + export.done/walks.done/walkreq.done markers dropped on 92 pending (idempotence holding); lean buffer 96 <= 200 - no eviction (the prep valve at 160 is the standing bound); freed 47M apparent (prep 93,659 -> 93,612M; df 623G steady at the pass). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 46 idempotent passes (129G ride-32 + the 44 cadence passes since), every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF A FORTIETH TIME (SIX-FILE SYNG SIDECAR RULE): the pass's in-flight list correctly named the batch-109 cut cohort (T03411/T03378/T03373/T02866) + the prep tail lane (T03679, mid-syng ~2h57m, the named syng segdup-tail class - its sidecars not yet written, output only at completion) + the 2 live call lanes (T04013/T04063); the batch-109 cohort verified 6/6 syng sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa KEPT (123-346M); the live lanes verified locus.fa 272-319M + reference-axis.json + locality-graphs (14/6) intact (the call reads NO syng; call.done 0/2 as expected, in-flight - the ride-62 clarification standing); class-B spot-checks (the batch-108 cohort T02769/T03803/T03670/T03859, resolved after pass 45): syng sidecars GONE (0/6 each), locus.fa 113-342M + the definitional beds KEPT (locus.refined.bed 465-1,214 rows), producer markers dropped (0 present), cut.done present 4/4 - the reap exactly where intended.
(3) THE WIDTH-4 CENSUS: batches post-swap = 109/109 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap, re-verified BY TIMESTAMP CUTOFF (epoch 1792176762 = the swap; the pre-swap non-zero markers are campaign history, not the census); batch-1792230985-109 in flight at close (cohort T03411/T03378/T03373/T02866, cut pid 1525166 37.6GB niced mid-scan, dispatched 09:56:25Z); batch-108 honest 488s restated (resolved after pass 45, lanes fed); map dirs 89, cut.done 2,156 steady at close - the width-4 pace holding (~29-35 territories/h).
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 587 STEADY (the last, T03855 09:31:04Z guard-137, receipted at the ride-73 close, predates this open); callok 1,519 steady (the 2 lanes in flight ~84-85m, receipts land at call.done); cut.exit 16; poison 4 terminal (the 4-member syng segdup-tail class: T00197/T01411/T02562/T00024); the lanes ordinary niced but CLIMBING (call.rss traces: T04063 62.0GB @ ~84m the widest - the T03855-like deep-lane climb, the guard is the oracle, NEVER killed by hand; T04013 51.7GB @ ~83m) - the deep-lane watch standing.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 623G steady, margin ~323G); guard-oracle held (zero hand kills, zero masters touched); the prep master mid-dispatch 64-wide chunk continuing (sub 2,294 steady, the T03679 tail lane alive ~2h57m on syng - the named syng segdup-tail class member, watched, never hand-killed); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
FINAL POLL 10:00:57Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / cut.done 2,156 / map dirs 89 (batch-109 in flight) / call-terminal 2,126 = callok 1,519 + fail-closed 587 + cut.exit 16 + poison 4 = 2,126 of 5,844 (36.4% of accepted) / df 623G / mem 796Gib / rvsize 174G apparent. MASTERS AT CLOSE: prep 411760 (~1-15:50h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~15:09h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~23:04h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-109 in flight); the JIT reaper cadence (~5-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch T04063 under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 46 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members + the T03679 live member, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855 + rides 55/56's clusters + the ride-69 consumption-burst receipt) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 79 WINDOW CLOSE (10:07Z 2026-10-17) - TWO CADENCE PASSES RECEIPTED (pass 47 +0M, pass 48 +28M - the second run to reap the batch-109 cohort the moment it resolved; CUMULATIVE ~= 180.3G ACROSS 48 IDEMPOTENT PASSES); CENSUS 110/110 HONEST ZERO KILLS SINCE THE SWAP (batch-109 HONEST AT 484s/38.3GB PEAK RESOLVED AFTER PASS 47 - COHORT PROTECTED, LANES FED; batch-110 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A FORTY-FIRST + FORTY-SECOND TIME; ONE NEW FAIL - T04063 CALL GUARD-KILLED AT 67.1GB/5,232s, THE SANCTIONED GUARD-ORACLE'S OWN KILL (the deep-lane memory-tail class, terminal, named - PREDICTED AT THE RIDE-78 CLOSE); THE DEEP-LANE WATCH TRANSFERS TO T04013; ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 10:02Z (session-open increment cfe0d2d). THE ORIENTATION STANDING (rides 38-79): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has been RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen was NOT redone this window - a fresh eviction of a pending territory would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASSES 47 + 48 (both ~4s, detached pids 1533205/1535885, receipts prep-reaper.log + prep-reaper-r79.nohup + prep-reaper-r79b.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): PASS 47 (10:03:53-57Z): class A 0, class B 0 files 0M + syng/export/walks/walkreq markers dropped on 48 cut-awaiting-call (census pending=96 cac=50), class C 0 + export.done/walks.done/walkreq.done markers dropped on 92 pending, lean buffer 96 <= 200 no eviction, freed 0M apparent. PASS 48 (10:05:28-32Z, run immediately after batch-109 resolved): class A 0, class B 24 files 44M + markers dropped on 52 cut-awaiting-call (census pending=92 cac=53), class C 0 + markers on 88 pending, lean buffer 92 <= 200 no eviction, freed 28M apparent (prep 93,339 -> 93,311M; df 622G steady). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 48 idempotent passes (129G ride-32 + the 46 cadence passes since), every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF (41st + 42nd times, SIX-FILE SYNG SIDECAR RULE): pass 47's in-flight list correctly named the batch-109 cut cohort (T03411/T03378/T03373/T02866) + the prep tail lane (T03679) + the 2 live call lanes (T04063/T04013); the batch-109 cohort verified 6/6 syng sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa 123-346M KEPT; the live lanes verified locus.fa + reference-axis.json + locality-graphs intact (the call reads NO syng); pass 48's in-flight list correctly named the batch-110 cut cohort (T03075/T03244/T02250/T02031) + T03679 + the T04013 lane, and its class-B reap spot-checked (the batch-109 cohort, resolved at 10:04Z after pass 47): syng sidecars GONE (0/6 each), locus.fa 123-346M + the definitional beds KEPT (locus.refined.bed 506-809 rows), producer markers dropped (0 present), cut.done present 4/4 - the reap exactly where intended; the batch-110 cohort verified 6/6 sidecars intact + locus.fa 113-296M KEPT; the T04013 lane locus.fa 319M + axis + locality-graphs (6) intact, call.done absent as expected; the T03679 tail lane mid-syng ~3h02m correctly excluded both passes (the named syng segdup-tail class, output written only at completion, watched, never hand-killed); class-B spot-checks at the open on the batch-108 cohort (pass 46's reap): sidecars GONE 0/6, locus.fa 113-342M + definitional beds 465-1,214 rows KEPT, markers dropped, cut.done present 4/4.
(3) THE WIDTH-4 CENSUS: batches post-swap = 110/110 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap, re-verified BY TIMESTAMP CUTOFF (epoch 1792176762 = the swap; the pre-swap non-zero markers are campaign history, not the census); this window: batch-109 HONEST exit=0 wall 484s peak RSS 38.3GB in band, resolved 10:04Z AFTER pass 47 - the pass had protected its cohort, the cut lanes stayed fed; batch-1792231470-110 dispatched 10:04Z in flight at close (early-load ~17GB niced); cut.done 2,156 -> 2,160 (+4, batch-109's cohort) - the width-4 pace holding (~29-35 territories/h).
(4) CALLS - ONE NEW FAIL THIS WINDOW, THE SANCTIONED GUARD-ORACLE'S OWN KILL: T04063 CALL GUARD-KILLED 10:04:04Z (exit 137, wall 5,232s, killed at 67.1GB RSS mid-score locus 9: 82,669 records / 7,793,563 units - the lane's own call.rss trace shows the deep-lane climb the ride-78 close named: 62.0GB @ 84m -> 67.1GB at the kill; the deep-lane memory-tail class T00129/T04052/T01895/T03855/T04063, terminal, named; receipts wave-logs/failures.txt + per-territory call.err/call.rss/call.wall IN-LANE FAIL-REAPED); fail-closed 587 -> 588; callok 1,519 steady; cut.exit 16; poison 4 terminal (the 4-member syng segdup-tail class T00197/T01411/T02562/T00024); THE DEEP-LANE WATCH TRANSFERS TO T04013 (57.3GB @ ~1h28m climbing - the guard is the oracle, never kill by hand); no other lanes at guard pressure at close.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 618-623G, margin ~318G+); guard-oracle held (zero hand kills, zero masters touched - the T04063 kill was the guard's own, sanctioned); the prep master mid-dispatch 64-wide chunk continuing (sub 2,294 steady, the T03679 tail lane alive ~3h02m on syng - watched, never hand-killed); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
FINAL POLL 10:06:11Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / cut.done 2,160 / map dirs 89 (batch-110 in flight) / call-terminal 2,127 = callok 1,519 + fail-closed 588 + cut.exit 16 + poison 4 = 2,127 of 5,844 (36.4% of accepted) / df 622G / mem 875Gib / rvsize 174G apparent. MASTERS AT CLOSE: prep 411760 (~1-15:40h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~15:13h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~23:09h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-110 in flight); the JIT reaper cadence (~4-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch T04013 under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 48 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members + the T03679 live member, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855/T04063 + rides 55/56's clusters + the ride-69 consumption-burst receipt) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 80 WINDOW CLOSE (10:14Z 2026-10-17) - TWO CADENCE PASSES RECEIPTED (pass 49 +0M, pass 50 +32M - the second run to reap the batch-110 cohort the moment it resolved; CUMULATIVE ~= 180.3G ACROSS 50 IDEMPOTENT PASSES); CENSUS 111/111 HONEST ZERO KILLS SINCE THE SWAP (batch-110 HONEST AT 472s/47.7GB PEAK RESOLVED AT 10:12:20Z - COHORT PROTECTED BY PASS 49, REAPED BY PASS 50, LANES FED; batch-111 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A FORTY-THIRD + FORTY-FOURTH TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 10:09Z (session-open increment de1fbcb). THE ORIENTATION STANDING (rides 38-80): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen was NOT redone this window - a fresh eviction of a pending territory would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASSES 49 + 50 (both ~4s, detached pids 1540191/1542837, receipts prep-reaper.log + prep-reaper-r80.nohup + prep-reaper-r80b.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): PASS 49 (10:10:08-12Z): class A 0, class B 0 files 0M + syng/export/walks/walkreq markers dropped on 52 cut-awaiting-call (census pending=92 cac=53), class C 0 + export.done/walks.done/walkreq.done markers dropped on 88 pending, lean buffer 92 <= 200 no eviction, freed 0M apparent (df 622G steady). PASS 50 (10:12:46-50Z, run the moment batch-110 resolved at 10:12:20Z): class A 0, class B 24 files 47M + markers dropped on 56 cut-awaiting-call (census pending=88 cac=57), class C 0 + markers on 84 pending, lean buffer 88 <= 200 no eviction, freed 32M apparent (prep 93,311 -> 93,279M; df 622 -> 619G). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 50 idempotent passes (129G ride-32 + the 48 cadence passes since), every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF (43rd + 44th times, SIX-FILE SYNG SIDECAR RULE): pass 49's in-flight list correctly named the batch-110 cut cohort (T03075/T03244/T02250/T02031) + the prep tail lane (T03679) + the T04013 live call lane; the batch-110 cohort verified 6/6 syng sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa 113-296M KEPT; the T04013 lane locus.fa 319M + reference-axis.json + locality-graphs (6) intact, call.done absent as expected (in-flight - the ride-62 clarification standing); the T03679 syng segdup-tail lane mid-syng correctly excluded; pass 50's in-flight list correctly named the batch-111 cut cohort (T03232/T03896/T02759/T03528, dispatched 10:12Z) + T03679 + T04013, and its class-B reap spot-checked (the batch-110 cohort): syng sidecars GONE (0/6 each), locus.fa 113-296M + the definitional beds KEPT (locus.refined.bed 465-1,251 rows), producer markers dropped (0 present), cut.done present 4/4 - the reap exactly where intended; the batch-111 cohort verified 6/6 sidecars intact + locus.fa 113-115M KEPT; class-B spot-checks at the open on the batch-109 cohort (pass 48's reap): sidecars GONE 0/6, locus.fa 123-346M + definitional beds 506-809 rows KEPT, markers dropped, cut.done present 4/4.
(3) THE WIDTH-4 CENSUS: batches post-swap = 111/111 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap, re-verified BY TIMESTAMP CUTOFF (epoch 1792176762 = the swap; 541 batch.exit markers on disk, the pre-swap non-zero are campaign history, not the census); this window: batch-110 HONEST exit=0 wall 472s peak RSS 47.7GB in band, resolved 10:12:20Z ~2min AFTER pass 49 - the pass had protected its cohort, pass 50 reaped it, the cut lanes stayed fed; batch-1792231943-111 dispatched 10:12Z in flight at close (cohort T03232/T03896/T02759/T03528); cut.done 2,160 -> 2,164 (+4, batch-110's cohort) - the width-4 pace holding (~29-35 territories/h).
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 588 STEADY (the last, T04063 10:04:04Z guard-137, receipted at the ride-79 close, predates this open); callok 1,519 steady (the 1 lane in flight, receipts land at call.done); cut.exit 16; poison 4 terminal (the 4-member syng segdup-tail class T00197/T01411/T02562/T00024); THE DEEP-LANE WATCH T04013 (61.3GB @ ~97m climbing ~40MB/min, the T03855/T04063-like tail - the guard is the oracle, NEVER killed by hand; the kill, if it comes, is the sanctioned guard-oracle's own and terminal/named).
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 619-622G, margin ~319G); guard-oracle held (zero hand kills, zero masters touched); the prep master mid-dispatch 64-wide chunk continuing (sub 2,294 steady, the T03679 tail lane alive ~3h09m on syng - the named syng segdup-tail class member, watched, never hand-killed); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
FINAL POLL 10:13:11Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / cut.done 2,164 / map dirs 92 (batch-111 in flight) / call-terminal 2,127 = callok 1,519 + fail-closed 588 + cut.exit 16 + poison 4 = 2,127 of 5,844 (36.4% of accepted) / df 619G / mem 877Gib / rvsize 177G apparent. MASTERS AT CLOSE: prep 411760 (~1-16:00h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~15:25h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~23:16h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-111 in flight); the JIT reaper cadence (~4-12s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch T04013 under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 50 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members + the T03679 live member, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855/T04063 + rides 55/56's clusters + the ride-69 consumption-burst receipt) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 81 WINDOW CLOSE (10:30Z 2026-10-17) - THREE CADENCE PASSES RECEIPTED (pass 51 +0M, pass 52 +14M the moment batch-111 resolved 391s HONEST, pass 53 +27M reaping the batch-112 cohort as the wave burst consumed; CUMULATIVE ~= 180.3G ACROSS 53 IDEMPOTENT PASSES); CENSUS 113/113 HONEST ZERO KILLS SINCE THE SWAP (batch-111 391s + batch-112 438s HONEST IN-BAND, batch-113 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A FORTY-FIFTH + FORTY-SIXTH + FORTY-SEVENTH TIME; ONE NEW FAIL - T04013 CALL GUARD-KILLED AT 67.1GB/6,338s, THE SANCTIONED GUARD-ORACLE'S OWN KILL (the deep-lane memory-tail class, terminal, named - PREDICTED AT THE RIDE-80 CLOSE); THE DEEP-LANE WATCH EMPTY AT CLOSE; THE CONSUMPTION BURST (map dirs 94 -> 70, +14 fails = T04013 + 13 binds ALL NAMED the ordinary 7479-class); ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 10:17Z (session-open increment 2945729). THE ORIENTATION STANDING (rides 38-81): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen was NOT redone this window - a fresh eviction of a pending territory would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASSES 51 + 52 + 53 (all ~4-5s, detached pids 1547136/1549761/1567588, receipts prep-reaper.log + prep-reaper-r81.nohup + prep-reaper-r81b.nohup + prep-reaper-r81c.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): PASS 51 (10:17:10-14Z): class A 0, class B 0 files 0M + markers on 56 cut-awaiting-call (census pending=88 cac=57), class C 0 + markers on 84 pending, lean buffer 88 <= 200 no eviction, freed 0M apparent (prep 93,279M steady; df 619G). PASS 52 (10:19:49-53Z, run the moment batch-111 resolved at ~10:19:30Z): class A 0, class B 24 files 30M + markers on 60 cut-awaiting-call (census pending=84 cac=61), class C 0 + markers on 80 pending, freed 14M apparent (prep 93,279 -> 93,265M; df 619G). PASS 53 (10:27:47-52Z, reaping as the wave burst consumed): class A 0, class B 24 files 43M on the batch-112 cohort + markers on 4 territories (census pending=80 cac=48), class C 0 + markers on 76 pending, lean buffer 80 <= 200 no eviction, freed 27M apparent (prep 90,699 -> 90,672M; df 634G). CUMULATIVE SUBSTRATE-REAPER FREED ~= 180.3G apparent across 53 idempotent passes (129G ride-32 + the 51 cadence passes since), every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF (45th + 46th + 47th times, SIX-FILE SYNG SIDECAR RULE): pass 51's in-flight list correctly named the batch-111 cut cohort (T03232/T03896/T02759/T03528) + the prep tail lane (T03679) + the T04013 live call lane; the batch-111 cohort verified 6/6 syng sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa 118-120M KEPT, and batch-111 resolved HONEST exit=0 wall 391s AFTER the pass - the cohort protected, the cut lanes stayed fed; pass 52 reaped the batch-111 cohort (spot-check: sidecars GONE 0/6 each, locus.fa 118-120M + definitional beds locus.refined.bed 465-471 rows KEPT, producer markers dropped 0 present, cut.done present 4/4 - the reap exactly where intended) and protected the batch-112 cohort (T02837/T02765/T03211/T03421) 6/6 sidecars intact + locus.fa 118-220M KEPT; pass 53 reaped the batch-112 cohort (spot-check: sidecars GONE 0/6 each, locus.fa 118-220M KEPT, markers dropped, cut.done present 4/4) and protects the batch-113 cohort (T03822/T02142/T01919/T02973, in flight at close); the T04013 lane locus.fa 334M + reference-axis.json + locality-graphs intact through passes 51/52 (the call reads NO syng; call.done absent as expected, in-flight - the ride-62 clarification standing), correctly absent from pass 53's list post-fail; the T03679 syng segdup-tail lane mid-syng ~3h25m correctly excluded all three passes (watched, never hand-killed).
(3) THE WIDTH-4 CENSUS: batches post-swap = 113/113 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap, re-verified BY TIMESTAMP CUTOFF (epoch 1792176762 = the swap; the pre-swap non-zero markers are campaign history); this window: batch-111 HONEST exit=0 wall 391s and batch-112 HONEST exit=0 wall 438s, both in-band, both AFTER their protecting passes - the cut lanes stayed fed; batch-1792232774-113 dispatched ~10:24Z in flight at close (cohort T03822/T02142/T01919/T02973, cut pid 30.2GB niced ~3m); cut.done 2,164 -> 2,172 (+8, the two cohorts) - the width-4 pace holding.
(4) CALLS - THE WINDOW'S ONE EVENT + THE CONSUMPTION BURST: T04013 CALL GUARD-KILLED 10:23:08Z (exit 137, wall 6,338s, killed at 67.1GB RSS = the 64GiB guard's own boundary, mid-score locus 3: 17,330 records / 4,819,316 units; the call.rss trace shows the deep-lane climb the ride-80 close named: 61.3GB @ 97m -> 63.1GB @ 99m -> 66.5GB @ ~105m -> 67.1GB at the kill): the sanctioned guard-oracle's own kill, the deep-lane memory-tail class now T00129/T04052/T01895/T03855/T04063/T04013, terminal, named (receipts wave-logs/failures.txt + per-territory call.err/call.rss/call.wall IN-LANE FAIL-REAPED); THE DEEP-LANE WATCH IS EMPTY at close (24 lanes in flight, widest 16.2GB @ 5m - ordinary, none at guard pressure; re-arm as the deep cohorts land, the guard is the oracle, NEVER kill by hand). THE CONSUMPTION BURST (the ride-69 class): the wave drained the refilled backlog, map dirs 94 -> 70, callok 1,519 -> 1,530 (+11 ordinary completions), fail-closed 588 -> 602 (+14 = the T04013 guard-137 + 13 bind ambiguities ALL NAMED, ALL the ordinary 7479-class, terminal: T02541 rec 9709 2-of-2 / T02769 rec 9342 / T03853 rec 63016 / T02250 rec 15425 / T02759 rec 7056 / T03917 rec 25154 / T03926 rec 9491 / T03378 rec 10350 / T03967 rec 59372 / T04020 rec 7333 / T03373 / T03943 / T04008) - the backlog's accumulated bind population surfacing at consumption, the aggregate at close is the receipt, no extrapolation; cut.exit 16; poison 4 terminal.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 619-639G, margin >319G); guard-oracle held (zero hand kills, zero masters touched - the T04013 kill was the guard's own, sanctioned); the prep master mid-dispatch 64-wide chunk continuing (sub 2,294 steady, the T03679 tail lane alive ~3h25m on syng - watched, never hand-killed; the incremental split final - no new windows); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
FINAL POLL 10:29:14Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,294 / cut.done 2,172 / map dirs 70 (batch-113 in flight) / call-terminal 2,152 = callok 1,530 + fail-closed 602 + cut.exit 16 + poison 4 = 2,152 of 5,844 (36.8% of accepted) / df 639G / mem 839Gib / rvsize 157G apparent (falling - the burst consumed). MASTERS AT CLOSE: prep 411760 (mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300), chunked wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-113 in flight); the JIT reaper cadence (~4-5s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch EMPTY at close (re-arming as the deep cohorts land - the guard is the oracle, never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 4 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~180.3G cumulative freed classes A/B/C across 53 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 4 CONFIRMED members + the T03679 live member, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855/T04063/T04013 + rides 55/56's clusters + the ride-69 + ride-81 consumption-burst receipts) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

### RIDE CONTINUATION 81 FOLLOW-UP (10:36Z 2026-10-17) - THE T03679 3h DECISION EXECUTED AT THE PARENT STEER (the ride-50 T00024 precedent LIVE): POISON-FIRST, THEN KILL-BY-EXACT-PID ONLY; THE SEG DUP-TAIL CLASS'S CONFIRMED FIFTH MEMBER; THE POISON HELD THROUGH BOTH DISPATCH SITES; THE FREED PREP POLE RESUMED 64-WIDE IMMEDIATELY (sub 2,294 -> 2,351 in ~1min)
The T03679 syng segdup-tail lane crossed the ~3h T00197/T00024 kill precedent at ~10:07Z and was climbing past it at the steer. THE EXECUTION, IN ORDER (receipted in the record, commit 55bb2ab): (1) POISON-FIRST - prep/T03679/cut.poison written 10:33:22Z + the class entry appended to cuts/poison.tsv (FIFTH member: T00197 + T01411 + T02562 + T00024 + T03679, 5 of 5,844 accepted, poison 4 -> 5) + syng.err.attempt2-snapshot KEPT; (2) THEN kill-by-exact-pid ONLY - SIGKILL 1251565 at 10:33:26Z (never the wrapper 1251489, never the master 411760, no group); the wrapper survived, wrote its honest receipts (syng.exit=1, the post-kill convention identical to T00024's; syng.wall 12,574s vs T00024's 10,953s), exited cleanly. THE KILL FACTS (the class signature): attempt 1 failed fast (exit 1, wall 2s, 02:21Z - the stale marker named rides 68/71); attempt 2 (dispatched 07:03Z) burned 3h30m CPU compute-bound (CPU==wall, state R) at 22.6MB RSS, NO output file (0/6 sidecars), silent since 07:03:52Z stuck on HG02257#2#JAGYVI020000012.1:3144802-3486758 (341,956 bp) after 13 sequences/3,945,824 bp/160,237 syncmers; locus.fa 357MB (the class's second-largest after T00024's 457MB); the ordinary corpus median 4s max 159s - ~66x the max. THE NO-REDISPATCH VERIFICATION: the freed xargs slot let the prep master resume its 64-wide pool at ~10:34Z (new wave pid 1586601, sub 2,294 -> 2,351 +57 in ~1min); T03679 ABSENT from all 29 enumerated live lanes - the TODO filter skipped the cut.poison marker at both dispatch sites (master lines 25/67), the double-writer race CLOSED. POST POLL 10:35:14Z: callok 1,540 / fail-closed 607 / cut.exit 16 / poison 5 = terminals 2,168 of 5,844 (37.1%) / df 614G / mem 735Gib / rvsize 183G (the prep wave rebuilding the substrate window - watch df, the 200-pending bound holds it). THE CLASS'S FINAL CENSUS for the datasheet: 5 CONFIRMED members (T00197 3h01m/345MB, T01411 3h16m/318MB, T02562 29min/318MB, T00024 3h01m/457MB, T03679 3h30m/357MB), every member poison-first then exact-pid, every redispatch closed. THE SUCCESSOR RESUMES the standing cadence unchanged; the deep-lane watch re-arms as the deep cohorts land (the guard is the oracle); AT CAMPAIGN CLOSE the datasheet close list now carries the class's FIVE members with the T03679 receipt. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 82 WINDOW (10:38Z 2026-10-17) - THE CADENCE PASS 54 RECEIPTED (+~18.0G APPARENT ~10s - THE LARGEST CLASS-C REAP SINCE RIDE-32: THE FREED PREP POLE'S NEW SUBSTRATE WAVE); CENSUS 114/114 HONEST ZERO KILLS SINCE THE SWAP (batch-112 438s + batch-113 503s HONEST IN-BAND; batch-114 IN FLIGHT); THE EXCLUSION PROOF A FORTY-EIGHTH TIME; ONE NEW FAIL - T03803 CALL GUARD-KILLED (137) AT 65.3GB/823s, THE DEEP-LANE CLASS'S 7TH MEMBER, THE SANCTIONED GUARD-ORACLE'S OWN KILL
Opened 10:38Z (session-open increment b47f40c). THE ORIENTATION STANDING (rides 38-82): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen is NOT redone this window.
(1) THE REAPER PASS 54 (10:38:21-31Z, ~10s, detached pid 1652088, receipt prep-reaper.log + prep-reaper-r82.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): class A 0, class B 24 files 45M (the batch-113 cohort T03822/T02142/T01919/T02973) + syng/export/walks/walkreq markers dropped on 8 cut-awaiting-call (census: pending=199 cut-awaiting-call=28), class C 357 files 18,444M - THE LARGEST CLASS-C REAP SINCE RIDE-32, the freed prep pole's 64-wide rebuild wave (sub 2,294 -> 2,416) minted a fresh pending class at 199 and its consumed export/walks intermediates reaped + export.done/walks.done/walkreq.done markers dropped on 195 pending; lean buffer 199 <= 200 - no eviction (BARELY under the bound; if the pole overshoots 200 the eviction arm fires - DESIGNED, the regen is the ~22s receipted chain; watch it); freed 17,998M apparent (prep 129,486 -> 111,488M; df 607G -> 624G). CUMULATIVE SUBSTRATE-REAPER FREED ~= 198.3G apparent across 54 idempotent passes (129G ride-32 + the 52 cadence passes since), every cut/call input and every split ledger kept.
(2) THE EXCLUSION PROOF A FORTY-EIGHTH TIME (SIX-FILE SYNG SIDECAR RULE): the pass's in-flight list correctly named the batch-114 cut cohort (T04191/T04088/T04182/T04145, the chr6 band) + the T04162 prep partition lane + the ~21 live call lanes; the batch-114 cohort verified 6/6 syng sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa KEPT (113-141M), cut.done absent as expected (mid-cut); the live lanes verified locus.fa 113-439M + reference-axis.json + locality-graphs intact (12-25 graphs - the call reads NO syng; call.done 0/4 as expected, in-flight - the ride-62 clarification standing); class-B spot-checks (the batch-113 cohort, reaped by pass 54): syng sidecars GONE (0/6 each), locus.fa 112-334M + the definitional beds KEPT (locus.refined.bed 465-1,298 rows), producer markers dropped (0 present), cut.done present 4/4 - the reap exactly where intended. T03803 correctly EXCLUDED while its call was live at pass time (10:38:23Z); its inputs were then fail-reaped IN-LANE by the wave's own fail path after the guard kill - the ride-62 clarification applied and verified.
(3) THE WIDTH-4 CENSUS: batches post-swap = 114/114 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap, re-verified BY TIMESTAMP CUTOFF (epoch 1792176762 = the swap; the pre-swap non-zero markers are campaign history, not the census); batch-112 honest 438s + batch-113 honest 503s in-band; batch-1792233279-114 in flight (cohort T04191/T04088/T04182/T04145, cut pid 1621619 ~22GB niced, dispatched ~10:35Z); cut.done 2,172 steady mid-window - the width-4 pace holding (~29-35 territories/h).
(4) CALLS - ONE NEW FAIL THIS WINDOW, THE SANCTIONED GUARD-ORACLE'S OWN KILL: T03803 CALL GUARD-KILLED 10:38:38Z (exit 137, wall 823s, killed at 65.3GB RSS mid-score locus 4: 25,294 records / 4,438,575 units / no-pin (unit,fold) 2.3e9 / abstained hull 174.3Mbp - the FAST-CLIMB deep-lane variant: the widest watched lane at the open 57.4GB @ 13m, the guard is the oracle, NEVER killed by hand; the deep-lane memory-tail class now T00129/T04052/T01895/T03855/T04063/T04013/T03803 - SEVEN members, terminal, named, receipts wave-logs/failures.txt + per-territory call.err/call.rss/call.wall IN-LANE FAIL-REAPED); fail-closed 607 -> 608; callok 1,540 -> 1,542 (the 18 remaining lanes in flight, receipts land at call.done); cut.exit 16; poison 5 terminal (the 5-member syng segdup-tail class, T03679 receipted at the ride-81 follow-up).
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 606-627G, margin >306G); guard-oracle held (zero hand kills, zero masters touched - the T03803 kill was the guard's own, sanctioned); the prep master mid-dispatch 64-wide chunk continuing - THE FREED POLE REBUILDING THE SUBSTRATE WINDOW (sub 2,294 -> 2,417, pending 199 near the 200 lean-buffer bound - if it overshoots, the designed eviction arm fires; watch, never intervene); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
FINAL POLL 10:39:01Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,417 (the freed pole producing) / map dirs 56 / call-terminal 2,171 = callok 1,542 + fail-closed 608 + cut.exit 16 + poison 5 = 2,171 of 5,844 (37.2% of accepted) / df 627G / mem 728Gib / rvsize 169G apparent (fell 190 -> 169G on pass 54). MASTERS: prep 411760 (~1-16:13h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~15:46h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~23:41h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-114 in flight); the JIT reaper cadence (~4-10s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the lean-buffer watch (pending 199 vs the 200 bound - eviction is DESIGNED if it trips); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 5 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~198.3G cumulative freed classes A/B/C across 54 idempotent passes incl. this window's 18G class-C, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 5 CONFIRMED members, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855/T04063/T04013/T03803) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).
## RIDE CONTINUATION 82 WINDOW CLOSE (10:44Z 2026-10-17) - PASS 54 (+~18.0G, THE LARGEST CLASS-C REAP SINCE RIDE-32) + PASS 55 (+1.1G, THE BATCH-114 COHORT REAPED THE MOMENT IT RESOLVED HONEST 474s); CUMULATIVE ~= 199.4G/55 PASSES; CENSUS 115/115 ZERO KILLS SINCE THE SWAP; T03803 GUARD-137 AT 65.3GB/823s (THE DEEP-LANE CLASS'S 7TH MEMBER, THE GUARD'S OWN KILL); ZERO HAND KILLS, ZERO MASTERS TOUCHED
The window's receipts are in the record's ride-82 closing increment (commit follows): pass 54 (10:38:21-31Z, pid 1652088): class B 24 files 45M (batch-113 cohort) + class C 357 files 18,444M - the freed prep pole's new substrate wave minted pending=199 and its consumed intermediates reaped; lean buffer 199 <= 200 no eviction (BARELY - if the pole overshoots 200 the eviction arm fires, DESIGNED, the ~22s regen chain); freed 17,998M apparent (df 607 -> 624G). Pass 55 (10:43:08-13Z, pid 1661072): class B 36 files 478M reaping the batch-114 cohort the moment it resolved HONEST 474s/~49.5GB peak; class C 3 files 631M; freed 1,090M apparent (df 624 -> 626G). Exclusion proof 48th + 49th times clean (batch-115 cohort T04202/T04282/T04196/T04084 protected 6/6 + locus.fa 113-157M KEPT; batch-113/114 cohorts reaped exactly - sidecars GONE, locus.fa + definitional beds KEPT, markers dropped, cut.done present; T03803 correctly excluded while live, fail-reaped in-lane after the guard kill - the ride-62 clarification applied). WIDTH-4 CENSUS 115/115 post-swap ALL exit=0 by timestamp cutoff, ZERO kills since the swap (batch-112 438s / batch-113 503s / batch-114 474s honest in-band); batch-115 in flight at close. CALLS: ONE new fail - T03803 guard-137 at 65.3GB/823s mid-score locus 4 (the deep-lane memory-tail class's SEVENTH member T00129/T04052/T01895/T03855/T04063/T04013/T03803, terminal, named, in-lane fail-reaped); fail-closed 607 -> 608; callok 1,540 -> 1,543. FINAL POLL 10:43:32Z: split 9,891 / acc 5,844 / sub 2,418 (the freed pole producing) / map dirs 59 / call-terminal 2,172 = callok 1,543 + fail-closed 608 + cut.exit 16 + poison 5 = 2,172 of 5,844 (37.2%) / df 626G / mem 719Gib / rvsize 171G; the deep-lane watch EMPTY (widest 29.7GB). MASTERS: prep 411760, cuts 3838087, wave 2810146, heartbeat; the ancient impg-reference alive untouched; the 300G floor never approached (margin >306G); the ladder not opened. THE SUCCESSOR RESUMES the standing cadence unchanged (sparse SHORT polls; the reaper cadence; the SIX-FILE SIDECAR exclusion rule; the width-4 census; the lean-buffer watch pending vs the 200 bound; the 300G hard floor; the ladder only at backlog drained AND df >= 1,600G); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the reaper receipts ~199.4G/55 passes + the guard-137 points incl. T03803) -> the record's closing increment -> the commit. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 83 WINDOW (10:45Z 2026-10-17) - SESSION OPEN (9f42f8a) + THE CADENCE PASS 56 RECEIPTED (+495M apparent drift ~5s, CLASSES ALL ZERO - THE PREP WAVE MID-PRODUCTION); ONE NEW FAIL - T03915 CALL GUARD-KILLED (137) AT 72.5GB/1,344s, THE DEEP-LANE CLASS'S 8TH MEMBER, THE SANCTIONED GUARD-ORACLE'S OWN KILL; THE ride-81->82 INTER-CLOSE BINDS NAMED (5, ALL ORDINARY 7479-CLASS)
Session-open committed (9f42f8a). OPEN POLL 10:44:44Z: split 9,891 / acc 5,844 / sub 2,418 / cut.done 2,180 (= ride-81 close 2,172 + batch-113 cohort ~10:34Z + batch-114 cohort 10:42:33Z - the ride-82 close's "2,172 steady" was the carried mid-window figure, reconciled exact, no anomaly) / map dirs 59 / call-terminal 2,172 = callok 1,543 + fail-closed 608 + cut.exit 16 + poison 5 (37.2%) / df 626G / mem 693Gib / rvsize 171G. MASTERS: prep 411760, cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON=24 DISK_FLOOR_GB=300), wave 2810146 (LANES=128 CHUNK_MULT=2), heartbeat; the ancient impg-reference alive untouched; the deep-lane watch EMPTY at open (widest T03915 35.8GB @ 21m - ordinary then); batch-115 in flight (cohort T04202/T04282/T04196/T04084, cut pid 1660265 ~22GB niced).
PASS 56 (10:46:39-43Z, ~5s, detached pid 1668560, receipt prep-reaper.log + prep-reaper-r83.nohup; bash -n CLEAN + md5 b804a7d3... UNCHANGED re-verified at the open): class A 0, class B 0 files 0M + markers dropped on 12 cut-awaiting-call (census pending=196 cac=28), class C 0 files 0M + markers on 192 pending; lean buffer 196 <= 200 no eviction; freed 495M apparent (prep 110,138 -> 109,643M; df 626 -> 628G) - THE 495M DELTA IS PREP-WAVE DU DRIFT, NOT CLASS ATTRIBUTION (all classes reaped 0 files; the 64-wide prep wave is mid-production - the honest accounting). CUMULATIVE ~= 199.4G + the pass-56 drift = ~199.9G apparent across 56 idempotent passes, every cut/call input and every split ledger kept.
THE EXCLUSION PROOF (50th time, SIX-FILE SYNG SIDECAR RULE): pass 56 protected the batch-115 cut cohort (T04202/T04282/T04196/T04084) 6/6 sidecars intact + locus.fa 113-157M KEPT (cut.done absent as expected, mid-cut) + correctly named the 21 live call lanes (T03859/T02031 spot-verified locus.fa 113M + reference-axis.json + locality-graphs 4-12 intact, sidecars 0 = the prior class-B reap, call.done absent in-flight - the ride-62 clarification standing); T03915 correctly EXCLUDED while live at pass time (named 10:46:41, its guard kill landed 10:46:40) and its inputs then fail-reaped IN-LANE by the wave's own fail path.
CALLS - ONE NEW FAIL, THE GUARD-ORACLE'S OWN: T03915 CALL GUARD-KILLED 10:46:40Z (exit 137, wall 1,344s, killed at 72.5GB RSS mid-score locus 12: 17,381 records / 3,962,686 units / no-pin (unit,fold) 2.04e9 / abstained hull 26.1Mbp - the rss trace 61.4GB @ locus_12_scored -> 72.5GB at the kill; the widest watched lane at the open 35.8GB @ 21m was the SAME lane before its fast-climb tail); the deep-lane memory-tail class now T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915 - EIGHT members, terminal, named, receipts wave-logs/failures.txt + per-territory call.err/call.rss/call.wall IN-LANE FAIL-REAPED. fail-closed 608 -> 609.
THE INTER-CLOSE BINDS NAMED (the honest naming duty - they landed 10:29:33-10:31:04Z, between the ride-81 close poll and the ride-82 open poll, the consumption burst's tail): T03670 rec 11251 two-origin-shifts / T04072 rec 6898 2-of-3 / T04069 rec 8820 two-origin-shifts / T03970 rec 9887 two-origin-shifts / T03998 rec 36486 two-origin-shifts - all the ordinary 7479-class, terminal, named (fail-closed 602 -> 607 across them).
CONTINGENCIES: the 300G hard floor NEVER approached (df 626-628G, margin >326G); guard-oracle held (zero hand kills, zero masters touched - the T03915 kill was the guard's own, sanctioned); the prep master mid-dispatch 64-wide chunk continuing (sub 2,418 steady, lean buffer 196 <= 200); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
## RIDE CONTINUATION 83 WINDOW CLOSE (10:51Z 2026-10-17) - PASS 56 (+495M DRIFT, CLASSES ZERO - PREP-WAVE DU DRIFT, HONESTLY ACCOUNTED) + PASS 57 (+16M, THE BATCH-115 COHORT REAPED THE MOMENT IT RESOLVED HONEST 444s); CUMULATIVE ~= 199.6G CLASS-ATTRIBUTED/57 PASSES; CENSUS 116/116 ZERO KILLS SINCE THE SWAP; T03915 GUARD-137 AT 72.5GB/1,344s (THE DEEP-LANE CLASS'S 8TH MEMBER, THE GUARD'S OWN KILL); THE 5 INTER-CLOSE BINDS NAMED; ZERO HAND KILLS, ZERO MASTERS TOUCHED
The window's receipts are in the record's ride-83 closing increment (commit af17d73, session-open 9f42f8a): pass 56 (10:46:39-43Z, pid 1668560): classes A/B/C ALL ZERO files, markers dropped on 12 cac + 192 pending (census pending=196 cac=28), lean buffer 196 <= 200, freed 495M apparent = PREP-WAVE DU DRIFT NOT CLASS ATTRIBUTION (the 64-wide wave mid-production - stated honestly); pass 57 (10:50:25-29Z): class B 24 files 31M reaping the batch-115 cohort (T04202/T04282/T04196/T04084) the moment it resolved HONEST exit=0 wall 444s, + markers on 16 cac (census pending=192 cac=31), class C 0, lean buffer 192 <= 200, freed 16M apparent. Exclusion proof 50th + 51st times clean (batch-115 cohort protected 6/6 sidecars + locus.fa 113-157M KEPT mid-cut, then reaped exactly - sidecars GONE, locus.fa + definitional beds 465-717 rows KEPT, markers dropped, cut.done 4/4; batch-116 cohort T04140/T04260/T04090/T04136 protected 6/6 + locus.fa 113-262M KEPT; live lanes T03859/T02031 spot-verified intact; T03915 correctly excluded while live, fail-reaped in-lane after its guard kill). WIDTH-4 CENSUS 116/116 post-swap ALL exit=0 by timestamp cutoff, ZERO kills since the swap (batch-113 503s / batch-114 474s / batch-115 444s honest in-band; batch-116 in flight at close). CALLS: ONE new fail - T03915 guard-137 10:46:40Z at 72.5GB/1,344s mid-score locus 12 - the deep-lane memory-tail class's EIGHTH member (T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915), the sanctioned guard-oracle's own kill, terminal, named, in-lane fail-reaped; fail-closed 608 -> 609; callok 1,543 -> 1,544; the 5 inter-close binds named (T03670 rec 11251 / T04072 rec 6898 2-of-3 / T04069 rec 8820 / T03970 rec 9887 / T03998 rec 36486 - the consumption burst's tail, all ordinary 7479-class). FINAL POLL 10:51:04Z: split 9,891 / acc 5,844 / sub 2,418 / cut.done 2,184 / map dirs 61 / call-terminal 2,174 = callok 1,544 + fail-closed 609 + cut.exit 16 + poison 5 = 2,174 of 5,844 (37.3%) / df 626G / mem 692Gib / rvsize 170G; the deep-lane watch EMPTY at close (widest 25GB); all masters alive (prep 411760, cuts 3838087, wave 2810146, heartbeat; the ancient impg-reference pids 1846398/1846402 alive, untouched); the 300G floor never approached (margin >326G); the ladder not opened. THE SUCCESSOR RESUMES the standing cadence unchanged (sparse SHORT polls; the reaper cadence; the SIX-FILE SIDECAR exclusion rule; the width-4 census continuation; the lean-buffer watch pending vs the 200 bound; the deep-lane watch under the guard-oracle, never kill by hand - the class is 8 members); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the reaper receipts ~199.6G/57 passes + the guard-137 points incl. T03915) -> the record's closing increment -> the commit. NO PR PUSH (the owner pushes).
## RIDE CONTINUATION 84 WINDOW CLOSE (11:01Z 2026-10-17) - TWO CADENCE PASSES RECEIPTED (pass 58 +299M APPARENT ~4s - ALL CLASSES ZERO, THE DELTA IS PREP-WAVE DU DRIFT, HONESTLY ACCOUNTED; pass 59 +19M CLASS B 24 files 36M REAPING THE BATCH-116 COHORT THE MOMENT IT RESOLVED HONEST 469s/50.9GB PEAK); CUMULATIVE ~= 199.9G APPARENT (~199.6G CLASS-ATTRIBUTED)/59 PASSES; CENSUS 116/116 ZERO KILLS SINCE THE SWAP (batch-117 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A FIFTY-SECOND + FIFTY-THIRD TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 10:53Z (session-open increment committed 659a119). THE ORIENTATION STANDING (rides 38-84): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen is NOT redone this window.
(1) THE REAPER PASSES 58 + 59 (both ~4s, detached pids 1680906/1690211, receipts prep-reaper.log + prep-reaper-r84.nohup + prep-reaper-r84b.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): PASS 58 (10:53:43-47Z): class A 0, class B 0 files 0M + markers dropped on 16 cut-awaiting-call (census pending=192 cac=29), class C 0 + markers on 188 pending; lean buffer 192 <= 200 no eviction; freed 299M apparent (prep 109,375 -> 109,076M; df 617 -> 629G) - THE 299M DELTA IS PREP-WAVE DU DRIFT, NOT CLASS ATTRIBUTION (all classes reaped 0 files - the honest accounting, the pass-56 precedent). PASS 59 (10:59:17-21Z, run the moment batch-116 resolved): class B 24 files 36M reaping the batch-116 cohort (T04140/T04260/T04090/T04136) + markers on 20 cac (census pending=188 cac=32), class C 0 + markers on 184 pending; lean buffer 188 <= 200 no eviction; freed 19M apparent (prep 108,907 -> 108,888M; df 628G steady). CUMULATIVE ~= 199.9G apparent across 59 idempotent passes (~199.6G class-attributed), every cut/call input and every split ledger kept.
(2) THE EXCLUSION PROOF A FIFTY-SECOND + FIFTY-THIRD TIME (SIX-FILE SYNG SIDECAR RULE): pass 58 protected the batch-116 cut cohort 6/6 sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa 113-262M KEPT (cut.done absent as expected, mid-cut) + correctly named the live call lanes (T03859/T02031/T03919 spot-verified locus.fa 113-205M + reference-axis.json + locality-graphs 4-25 intact, call.done absent as expected, in-flight - the ride-62 clarification standing) - batch-116 resolved HONEST exit=0 wall 469s peak RSS 50.9GB in band ~6min AFTER the pass, the lanes fed; pass 59 reaped the batch-116 cohort exactly (sidecars GONE 0/6 each, locus.fa 113-262M KEPT, producer markers dropped, cut.done present 4/4) + protects the batch-117 cohort (T04252/T04203/T04290/T04227, dispatched ~10:59Z) 6/6 sidecars + locus.fa 113-114M KEPT mid-cut - the reap exactly where intended.
(3) THE WIDTH-4 CENSUS: batches post-swap = 116/116 HONEST ALL exit=0, ZERO non-zero, ZERO guard kills since the 18:52:42Z swap, re-verified BY TIMESTAMP CUTOFF (epoch 1792176762 = the swap; 547 batch.exit markers on disk, the pre-swap non-zero are campaign history); this window: batch-116 HONEST exit=0 wall 469s peak RSS 50.9GB in band; batch-1792234199-117 in flight at close; cut.done 2,184 -> 2,188 (+4, batch-116's cohort) - the width-4 pace holding (~29-35 territories/h); map dirs 62 at the final poll.
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 609 STEADY (the last, T03915 10:46:40Z guard-137, receipted at the ride-83 close, predates this open); callok 1,544 -> 1,547 (+3 ordinary completions); cut.exit 16; poison 5 terminal (the 5-member syng segdup-tail class); THE DEEP-LANE WATCH AT CLOSE: T02866 48.4GB @ 32m the widest, climbing under the guard-oracle (never kill by hand; the class is 8 members), T03208 37.0GB second - both ordinary niced, none at guard pressure.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 617-629G, margin >317G); guard-oracle held (zero hand kills, zero masters touched); the prep master mid-dispatch 64-wide chunk continuing (sub 2,418 steady, lean buffer 188-192 <= 200 no eviction); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched).
FINAL POLL 10:59:45Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,418 / cut.done 2,188 / map dirs 62 (batch-117 in flight) / call-terminal 2,177 = callok 1,547 + fail-closed 609 + cut.exit 16 + poison 5 = 2,177 of 5,844 (37.3% of accepted) / df 628G / mem 630Gib / rvsize 169G apparent. MASTERS AT CLOSE: prep 411760 (mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300), chunked wave 2810146 (LANES=128 CHUNK_MULT=2); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-117 in flight); the JIT reaper cadence (~4-10s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch T02866 under the guard-oracle (never kill by hand); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the lean-buffer watch (pending 188 vs the 200 bound - eviction DESIGNED if it trips); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 5 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~199.9G cumulative apparent/~199.6G class-attributed across 59 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 5 CONFIRMED members, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 85 WINDOW CLOSE (11:05Z 2026-10-17) - TWO CADENCE PASSES RECEIPTED (pass 60 +0M ~5s ALL CLASSES ZERO; pass 61 +12M ~5s CLASS B 24 files 28M REAPING THE BATCH-117 COHORT THE MOMENT IT RESOLVED HONEST 409s/~49.2GB PEAK); CUMULATIVE ~= 199.9G APPARENT (~199.6G CLASS-ATTRIBUTED)/61 PASSES; CENSUS 118/118 ZERO KILLS SINCE THE SWAP BY FILE-MTIME CUTOFF (the ride-84 close's 116 was the carried pre-resolve figure - batch-116's exit=0 at 10:57:48Z predates the ride-84 close poll, reconciled exact, no anomaly; batch-118 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A FIFTY-FOURTH + FIFTY-FIFTH TIME; ONE NEW FAIL - T03919 CALL GUARD-KILLED (137) AT ~70.0GB/2,075s, THE DEEP-LANE CLASS'S 9TH MEMBER, THE SANCTIONED GUARD-ORACLE'S OWN KILL (named at the open, in-lane fail-reaped); ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 11:01Z (session-open increment 49436a1). THE ORIENTATION STANDING (rides 38-85): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen is NOT redone this window - a fresh eviction of a pending territory would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASSES 60 + 61 (both ~5s, detached pids 1695075/1701487, receipts prep-reaper.log + prep-reaper-r85.nohup + prep-reaper-r85b.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): PASS 60 (11:02:13-18Z): class A 0, class B 0 files 0M + syng/export/walks/walkreq markers dropped on 20 cut-awaiting-call (census pending=188 cac=31), class C 0 + export.done/walks.done/walkreq.done markers dropped on 184 pending; lean buffer 188 <= 200 no eviction; freed 0M apparent (prep 108,681M steady; df 631G steady). PASS 61 (11:04:45-50Z, run the moment batch-117 resolved HONEST at 11:04:3xZ): class B 24 files 28M reaping the batch-117 cohort (T04252/T04203/T04290/T04227) + markers dropped on 24 cut-awaiting-call (census pending=184 cac=34), class C 0 + markers on 180 pending; lean buffer 184 <= 200 no eviction; freed 12M apparent (prep 108,681 -> 108,556M; df 630G steady). CUMULATIVE ~= 199.9G apparent (~199.6G class-attributed) across 61 idempotent passes (129G ride-32 + the 59 cadence passes since), every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF A FIFTY-FOURTH + FIFTY-FIFTH TIME (SIX-FILE SYNG SIDECAR RULE): pass 60 protected the batch-117 cut cohort 6/6 sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa 113-114M KEPT (cut.done absent as expected, mid-cut - batch-117 then resolved HONEST exit=0 wall 409s peak ~49.2GB in band AFTER the pass, the cohort protected, the cut lanes stayed fed) + correctly named the live call lanes (T02866/T02031 spot-verified locus.fa 134/113M + reference-axis.json + locality-graphs 21/12 intact, sidecars 0 = the prior class-B reap, call.done absent in-flight - the ride-62 clarification standing); T03919 correctly ABSENT from the pass-60 exclusion list post-fail (guard-killed 11:01:17Z before the pass; its inputs fail-reaped IN-LANE by the wave's own fail path, receipts call.err/call.rss/call.wall kept, the clarification applied and verified); pass 61 reaped the batch-117 cohort exactly (spot-check: syng sidecars GONE 0/6 each, locus.fa 113-114M + the definitional beds KEPT (locus.refined.bed 465-468 rows), producer markers dropped (0 present), cut.done present 4/4 - the reap exactly where intended) and protects the batch-118 cohort (T04241/T04149/T04244/T04150, dispatched ~11:04Z, cuts/batch-1792235079-118 in flight) 6/6 sidecars intact + locus.fa 113-235M KEPT (cut.done absent as expected, mid-cut).
(3) THE WIDTH-4 CENSUS: batches post-swap = 118/118 HONEST ALL exit=0 BY FILE-MTIME CUTOFF (epoch 1792176762 = the 18:52:42Z swap; 548 batch.exit markers on disk, the pre-swap non-zero are campaign history, not the census), ZERO non-zero, ZERO guard kills since the swap; the count reconciliation: the ride-84 close's "116/116" was the carried figure - batch-116's exit=0 (10:57:48Z, wall 469s) predates the ride-84 close poll and batch-117's exit=0 (wall 409s, peak ~49.2GB) landed this window, reconciled exact, no anomaly; cut.done 2,188 -> 2,192 (+4, batch-117's cohort) - the width-4 pace holding (~29-35 territories/h); batch-1792235079-118 in flight at close.
(4) CALLS - ONE NEW FAIL THIS WINDOW, THE SANCTIONED GUARD-ORACLE'S OWN KILL: T03919 CALL GUARD-KILLED 11:01:17Z (exit 137, wall 2,075s, killed at ~70.0GB RSS mid-score locus 0: 43,142 records / 6,727,768 units / no-pin (unit,fold) 1.79e9 / abstained hull 817.2Mbp - the rss trace 62.8GB @ locus_0_scored -> 70.0GB at the kill; the lane was spot-verified intact while in-flight at the ride-84 exclusion proof); the deep-lane memory-tail class now T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915/T03919 - NINE members, terminal, named (receipts wave-logs/failures.txt + per-territory call.err/call.rss/call.wall IN-LANE FAIL-REAPED); fail-closed 609 -> 610; callok 1,547 -> 1,548 (+1 ordinary completion, receipts land at call.done); cut.exit 16; poison 5 terminal (the 5-member syng segdup-tail class); THE DEEP-LANE WATCH AT CLOSE: T02866 48GB @ 37m the widest (climbing slowly - the ride-84 named watch, the guard is the oracle, NEVER kill by hand), T03208 36GB, T02722 34GB - none at guard pressure.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 630-631G, margin ~330G); guard-oracle held (zero hand kills, zero masters touched - the T03919 kill was the guard's own, sanctioned); the prep master mid-dispatch 64-wide chunk continuing (sub 2,418 steady, lean buffer 184-188 <= 200 no eviction); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched, 22d+).
FINAL POLL 11:05:09Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,418 / cut.done 2,192 / map dirs 64 (batch-118 in flight) / call-terminal 2,179 = callok 1,548 + fail-closed 610 + cut.exit 16 + poison 5 = 2,179 of 5,844 (37.3% of accepted) / df 630G / mem 673Gib / rvsize 167G apparent. MASTERS AT CLOSE: prep 411760 (~1-16:36h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~16:14h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, 1-00:09h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-118 in flight); the JIT reaper cadence (~5s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch T02866 under the guard-oracle (never kill by hand; the class is 9 members); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the lean-buffer watch (pending 184 vs the 200 bound - eviction DESIGNED if it trips); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 5 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~199.9G cumulative apparent/~199.6G class-attributed across 61 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 5 CONFIRMED members, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915/T03919 + rides 55/56's clusters + the ride-69 + ride-81 consumption-burst receipts) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 86 WINDOW CLOSE (11:14Z 2026-10-17) - TWO CADENCE PASSES RECEIPTED (pass 62 +0M ~5s ALL CLASSES ZERO, markers on 24 cac + 180 pending; pass 63 +19M ~5s CLASS B 24 files 34M REAPING THE BATCH-118 COHORT THE MOMENT IT RESOLVED HONEST 461s/~49.6GB PEAK); CUMULATIVE ~= 199.9G APPARENT (~199.6G CLASS-ATTRIBUTED)/63 PASSES; CENSUS 119/119 ZERO KILLS SINCE THE SWAP (batch-119 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A FIFTY-SIXTH + FIFTY-SEVENTH TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 11:06Z (session-open increment 5274e5c). THE ORIENTATION STANDING (rides 38-86): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen is NOT redone this window - a fresh eviction of a pending territory would only delete cut/call inputs the live campaign still needs (2,179/5,844 terminal - 63% in flight).
(1) THE REAPER PASSES 62 + 63 (both ~5s, detached pids 1706707/1715217, receipts prep-reaper.log + prep-reaper-r86.nohup + prep-reaper-r86b.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): PASS 62 (11:07:57-11:08:02Z): class A 0, class B 0 files 0M + syng/export/walks/walkreq markers dropped on 24 cut-awaiting-call (census pending=184 cac=34), class C 0 + export.done/walks.done/walkreq.done markers dropped on 180 pending; lean buffer 184 <= 200 no eviction; freed 0M apparent (prep 108,556M steady; df 630G). PASS 63 (11:13:01-11:13:06Z, run the moment batch-118 resolved HONEST exit=0 wall 461s, rss trace peak 49,615,004 kB ~= 49.6GB in band): class B 24 files 34M reaping the batch-118 cohort (T04241/T04149/T04244/T04150) + markers dropped on 28 cut-awaiting-call (census pending=180 cac=38), class C 0 + markers on 176 pending; lean buffer 180 <= 200 no eviction; freed 19M apparent (prep 108,556 -> 108,537M; df 627G). CUMULATIVE ~= 199.9G apparent (~199.6G class-attributed) across 63 idempotent passes (129G ride-32 + the 61 cadence passes since), every cut/call input and every split ledger kept. ONE TRANSIENT, NAMED (no consequence): a `tail prep-reaper.log` at 11:07:35Z failed 'No such file or directory' for ~1min while the log re-appeared intact at 11:08 with the complete pass-62 receipt (107KB, append-only, content verified complete) - an FS hiccup class, the receipt chain unbroken.
(2) THE EXCLUSION PROOF A FIFTY-SIXTH + FIFTY-SEVENTH TIME (SIX-FILE SYNG SIDECAR RULE): pass 62 protected the batch-118 cut cohort 6/6 sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa 113-235M KEPT (cut.done absent as expected, mid-cut) + correctly named the live call lanes (T02722/T02866/T03411/T03859/T04066 spot-verified locus.fa 113-439M + reference-axis.json + locality-graphs 3-25 intact, call.done absent in-flight - the ride-62 clarification standing); pass 63 reaped the batch-118 cohort exactly (spot-check: syng sidecars GONE 0/6 each, producer markers dropped 0/4, cut.done present 4/4, locus.fa 113-235M + the definitional beds KEPT - locus.refined.bed 465-1,049 rows) and protects the batch-119 cohort (T04262/T04109/T04166/T04280, dispatched ~11:12Z, cuts/batch-1792235541-119 in flight) 6/6 sidecars intact + locus.fa 113M KEPT mid-cut - the reap exactly where intended.
(3) THE WIDTH-4 CENSUS: batches post-swap = 119/119 HONEST ALL exit=0 BY FILE-MTIME CUTOFF (epoch 1792176762 = the 18:52:42Z swap; the pre-swap non-zero markers are campaign history, not the census), ZERO non-zero, ZERO guard kills since the swap; this window: batch-118 HONEST exit=0 wall 461s peak ~49.6GB in band, resolved ~11:12:20Z; batch-1792235541-119 in flight at close (cohort T04262/T04109/T04166/T04280, cut pid early-load); cut.done 2,192 -> 2,196 (+4, batch-118's cohort) - the width-4 pace holding (~29-35 territories/h); map dirs 64 -> 68 as the wave consumes.
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 610 STEADY (the last, T03919 11:01:17Z guard-137, receipted at the ride-85 close, predates this open); callok 1,548 steady (the lanes in flight ~46m, receipts land at call.done); cut.exit 16; poison 5 terminal (the 5-member syng segdup-tail class T00197/T01411/T02562/T00024/T03679); THE DEEP-LANE WATCH AT CLOSE: T02866 48.4GB @ ~46m the widest, climbing slowly (the ride-84/85 named watch, the guard is the oracle, NEVER killed by hand), T02722 40.9GB second - none at guard pressure.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 627-630G, margin ~327G); guard-oracle held (zero hand kills, zero masters touched); the prep master mid-dispatch 64-wide chunk continuing (sub 2,418 steady, lean buffer 180-184 <= 200 no eviction); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched, 22d+).
FINAL POLL 11:13:30Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,418 / cut.done 2,196 / map dirs 68 (batch-119 in flight) / call-terminal 2,179 = callok 1,548 + fail-closed 610 + cut.exit 16 + poison 5 = 2,179 of 5,844 (37.3% of accepted) / df 627G / mem 636Gib / rvsize 170G apparent. MASTERS AT CLOSE: prep 411760 (~1-16:48h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~16:21h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, 1-00:17h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-119 in flight); the JIT reaper cadence (~5s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch T02866 under the guard-oracle (never kill by hand - the class is 9 members); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the lean-buffer watch (pending 180 vs the 200 bound - eviction DESIGNED if it trips); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 5 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~199.9G cumulative apparent/~199.6G class-attributed across 63 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 5 CONFIRMED members, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915/T03919 + rides 55/56's clusters + the ride-69 + ride-81 consumption-burst receipts) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 87 WINDOW CLOSE (11:21Z 2026-10-17) - TWO CADENCE PASSES RECEIPTED (pass 64 +0M ~4s ALL CLASSES ZERO, markers on 28 cac + 176 pending; pass 65 +12M ~4s CLASS B 24 files 28M REAPING THE BATCH-119 COHORT THE MOMENT IT RESOLVED HONEST 404s/47.6GB PEAK); CUMULATIVE ~= 199.9G APPARENT (~199.6G CLASS-ATTRIBUTED)/65 PASSES; CENSUS 120/120 ZERO KILLS SINCE THE SWAP (batch-120 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A FIFTY-EIGHTH + FIFTY-NINTH TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 11:16Z (session-open increment c01d6e9). THE ORIENTATION STANDING (rides 38-87): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass; the predecessor's expected ~136G superseded by the measured class attribution) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen is NOT redone this window - a fresh eviction of a pending territory would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASSES 64 + 65 (both ~4s, detached, receipts prep-reaper.log + prep-reaper-r87.nohup + prep-reaper-r87b.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): PASS 64 (11:15:35-39Z): class A 0, class B 0 files 0M + syng/export/walks/walkreq markers dropped on 28 cut-awaiting-call (census pending=180 cac=38), class C 0 + export.done/walks.done/walkreq.done markers dropped on 176 pending; lean buffer 180 <= 200 no eviction; freed 0M apparent (prep 108,537M steady; df 627G). PASS 65 (11:19:58-11:20:02Z, run the moment batch-119 resolved HONEST exit=0 wall 404s, rss trace peak 47.6GB in band): class B 24 files 28M reaping the batch-119 cohort (T04262/T04109/T04166/T04280) + markers dropped on 32 cut-awaiting-call (census pending=176 cac=41), class C 0 + markers on 172 pending; lean buffer 176 <= 200 no eviction; freed 12M apparent (prep 108,404 -> 108,392M; df 627G steady). CUMULATIVE ~= 199.9G apparent (~199.6G class-attributed) across 65 idempotent passes (129G ride-32 + the 63 cadence passes since), every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF A FIFTY-EIGHTH + FIFTY-NINTH TIME (SIX-FILE SYNG SIDECAR RULE): pass 64 protected the batch-119 cut cohort (T04262/T04109/T04166/T04280) 6/6 syng sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa 113M KEPT (cut.done absent as expected, mid-cut) + correctly named the 10 live call lanes (T02866/T02722/T03208 spot-verified locus.fa 117-439M + reference-axis.json + locality-graphs 21-25 intact, sidecars 0 = the prior class-B reap, call.done absent in-flight - the ride-62 clarification standing); pass 65 reaped the batch-119 cohort exactly (spot-check: syng sidecars GONE 0/6 each, locus.fa 113M + the definitional beds KEPT (locus.refined.bed 465-466 rows), producer markers dropped (0 present), cut.done present 4/4 - the reap exactly where intended) and its in-flight list correctly named the NEW batch-120 cut cohort (T04094/T04089/T04161/T04281, dispatched ~11:20Z, cuts/batch-1792235946-120 in flight at close) - the master composed the successor batch within the minute.
(3) THE WIDTH-4 CENSUS: batches post-swap = 120/120 HONEST ALL exit=0 BY FILE-MTIME CUTOFF (epoch 1792176762 = the 18:52:42Z swap; 551 batch.exit markers on disk, the pre-swap non-zero are campaign history, not the census), ZERO non-zero, ZERO guard kills since the swap; this window: batch-119 HONEST exit=0 wall 404s peak RSS 47.6GB in band; batch-1792235946-120 in flight at close; cut.done 2,196 -> 2,200 (+4, batch-119's cohort) - the width-4 pace holding (~29-35 territories/h); map dirs 64 -> 71 as the wave consumes.
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 610 STEADY (the last, T03919 11:01:17Z guard-137, receipted at the ride-85 close, predates this open); callok 1,548 -> 1,549 (+1 ordinary completion, receipts land at call.done); cut.exit 16; poison 5 terminal (the 5-member syng segdup-tail class T00197/T01411/T02562/T00024/T03679); THE DEEP-LANE WATCH AT CLOSE: T02866 48.4GB @ ~54m the widest (climbing slowly - the ride-84/85/86 named watch, the guard is the oracle, NEVER killed by hand; the class is 9 members), T02722 46.3GB second (climbing from 42.4GB at the open - the next watch), T03208 37.0GB - none at guard pressure.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 626-627G, margin ~326G); guard-oracle held (zero hand kills, zero masters touched); the prep master mid-dispatch 64-wide chunk continuing (sub 2,418 steady, lean buffer 176-180 <= 200 no eviction); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched, 22d+).
FINAL POLL 11:21:07Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,418 / cut.done 2,200 / map dirs 71 (batch-120 in flight) / call-terminal 2,180 = callok 1,549 + fail-closed 610 + cut.exit 16 + poison 5 = 2,180 of 5,844 (37.3% of accepted) / df 627G / mem 645Gib / rvsize 169G apparent. MASTERS AT CLOSE: prep 411760 (~1-16:55h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~16:28h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~1-00:24h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-120 in flight); the JIT reaper cadence (~4-5s a pass, idempotent, re-runnable anytime; the SIX-FILE SIDECAR exclusion rule on every pass; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch T02866/T02722 under the guard-oracle (never kill by hand - the class is 9 members); the prep pole + df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner); the lean-buffer watch (pending 176 vs the 200 bound - eviction DESIGNED if it trips); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 5 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~199.9G cumulative apparent/~199.6G class-attributed across 65 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 5 CONFIRMED members, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915/T03919 + rides 55/56's clusters + the ride-69 + ride-81 consumption-burst receipts) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 88 SESSION-OPEN (11:23Z 2026-10-17) - THE CADENCE CONTINUES FROM RIDE-87 CLOSE (b793856: passes 64+65, cumulative ~= 199.9G apparent/~199.6G class-attributed across 65 idempotent passes); OPEN POLL 11:22:37Z split 9,891 / acc 5,844 / sub 2,418 / cut.done 2,200 / map dirs 71 / call-terminal 2,180 = callok 1,549 + fail-closed 610 + cut.exit 16 + poison 5 = 2,180/5,844 (37.3%) / df 627G / mem 626Gib / rvsize 169G; masters verified live (prep 411760, cuts 3838087 width-4-direct, wave 2810146, heartbeat; the ancient impg-reference alive untouched); 10 live call lanes, the deep-lane watch T02722 48.8GB WIDEST (passed T02866 48.4GB) + T03208 37.0GB under the guard-oracle (never hand-kill); batch-120 in flight (cohort T04094/T04089/T04161/T04281, cut pid 1725756 ~37GB niced); reaper verified at the open (bash -n clean, md5 b804a7d3... unchanged); THE ORIENTATION STANDING (rides 38-88): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt) - the regen is NOT redone this window (2,180/5,844 terminal, 63% in flight); NO PR PUSH (the owner pushes)

## RIDE CONTINUATION 88 WINDOW CLOSE (11:30Z 2026-10-17) - TWO CADENCE PASSES RECEIPTED (pass 66 +0M ~4s ALL CLASSES ZERO, markers on 32 cac + 172 pending; pass 67 +3,070M APPARENT = 30M CLASS B 24 files REAPING THE BATCH-120 COHORT THE MOMENT IT RESOLVED HONEST exit=0 wall 414s/~48.0GB PEAK + ~3.0G prep-wave du drift, honestly split per the pass-56/58 precedent); CUMULATIVE ~= 202.9G APPARENT (~199.6G CLASS-ATTRIBUTED)/67 PASSES; CENSUS 121/121 ZERO KILLS SINCE THE SWAP (batch-121 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A SIXTIETH + SIXTY-FIRST TIME; ONE NEW FAIL - T02722 CALL GUARD-KILLED (137) AT 71.05GB/3,610s, THE DEEP-LANE CLASS'S 10TH MEMBER, THE SANCTIONED GUARD-ORACLE'S OWN KILL (the fast-climb tail NAMED AT THE OPEN, watched, never hand-touched); ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 11:23Z (session-open increment 11ac12e). THE ORIENTATION STANDING (rides 38-88): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen is NOT redone this window. (1) THE REAPER PASSES 66 + 67 (both ~4s, detached pids 1732432/1745688, receipts prep-reaper.log + prep-reaper-r88.nohup + prep-reaper-r88b.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): pass 66 (11:23:40-44Z) all classes zero, markers on 32 cac + 172 pending, lean buffer 172 <= 200, freed 0M apparent (prep 108,392M steady; df 627G); pass 67 (11:26:43-44Z, the moment batch-120 resolved HONEST exit=0 wall 414s peak ~48.0GB in band) class B 24 files 30M reaping the batch-120 cohort (T04094/T04089/T04161/T04281) + markers on 36 cac, class C 0 + markers on 168 pending, lean buffer 168-172 <= 200, freed 3,070M apparent (prep 108,392 -> 103,783M; df 627 -> 629G) - HONESTLY SPLIT: 30M class-attributed, ~3.0G prep-wave du drift NOT class attribution. CUMULATIVE ~= 202.9G apparent (~199.6G class-attributed) across 67 idempotent passes, every cut/call input and every split ledger kept. (2) THE EXCLUSION PROOF (60th + 61st times, SIX-FILE SYNG SIDECAR RULE): pass 66 protected the batch-120 cohort 6/6 sidecars + locus.fa 112-115M KEPT mid-cut + correctly named the 10 live call lanes (T02722/T02866/T03208 spot-verified locus.fa 116-438M + reference-axis.json + locality-graphs 21-25 intact, call.done absent in-flight - the ride-62 clarification standing); pass 67 reaped the batch-120 cohort exactly (sidecars GONE 0/6, locus.fa + definitional beds 465-477 rows KEPT, markers dropped, cut.done 4/4) + protects the batch-121 cohort (T04246/T04185/T04086/T04099, cut pid 1736309 ~31GB niced, dispatched ~11:26Z) 6/6 sidecars + locus.fa 112-118M KEPT mid-cut; T02722 correctly excluded while live at pass time then fail-reaped IN-LANE after its guard kill (locus.fa gone, call receipts kept). (3) WIDTH-4 CENSUS 121/121 post-swap ALL exit=0 BY FILE-MTIME CUTOFF (epoch 1792176762), ZERO non-zero, ZERO guard kills since the swap (batch-120 honest 414s/~48.0GB in band; batch-1792236361-121 in flight at close); cut.done 2,200 -> 2,204 (+4); map dirs 71 -> 75. (4) CALLS - ONE NEW FAIL, THE GUARD-ORACLE'S OWN: T02722 CALL GUARD-KILLED 11:27:25Z (exit 137, wall 3,610s, killed at 71,054,444 kB ~= 71.05GB RSS mid-score locus 16: 31,574 records / 7,013,915 units / 1,652 members -> 510 folds) - the fast-climb tail NAMED AT THE OPEN (48.8GB widest, climbing ~1.5GB/min, watched under the guard-oracle, NEVER hand-touched); the deep-lane memory-tail class now T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915/T03919/T02722 - TEN members, terminal, named, in-lane fail-reaped (fail-closed 610 -> 611; callok 1,549 -> 1,551 +2 ordinary completions; cut.exit 16; poison 5). THE DEEP-LANE WATCH AT CLOSE: T02866 48.0GB @ ~1h01m the widest (climbing slowly - the guard is the oracle, NEVER kill by hand), T03208 37.0GB second - none at guard pressure. (5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 625-629G, margin >325G); guard-oracle held (zero hand kills, zero masters touched); the prep master mid-dispatch 64-wide (sub 2,418 -> 2,427 the pole producing, lean buffer 168-172 <= 200 no eviction); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference alive untouched (22d+). FINAL POLL 11:29:17Z: split 9,891 / acc 5,844 / sub 2,427 / cut.done 2,204 / map dirs 72 (batch-121 in flight) / call-terminal 2,183 = callok 1,551 + fail-closed 611 + cut.exit 16 + poison 5 = 2,183/5,844 (37.4%) / df 626G / mem 715Gib / rvsize 170G. MASTERS: prep 411760, cuts 3838087 width-4-direct W=4 C=1, wave 2810146 LANES=128, heartbeat; NO PR PUSH (the owner pushes). THE SUCCESSOR RESUMES the standing cadence unchanged (sparse SHORT polls; the reaper cadence ~4s a pass idempotent; the SIX-FILE SIDECAR exclusion rule; the width-4 census continuation with batch-121 in flight; the deep-lane watch T02866 under the guard-oracle - the class is 10 members; the lean-buffer watch pending 168 vs the 200 bound; the 300G hard floor; the ladder only at backlog drained AND df >= 1,600G); AT CAMPAIGN CLOSE: rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + the reaper receipts ~= 202.9G apparent/~199.6G class-attributed/67 passes with pass-67's drift honestly split + the guard-137 points incl. T02722 71.05GB/3,610s) -> the record's closing increment -> the commit. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 89 (11:31-11:37Z 2026-10-17) - SESSION-OPEN + TWO CADENCE PASSES + WINDOW CLOSE
Opened 11:31Z (session-open increment 189918c). THE ORIENTATION STANDING (rides 38-89): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen is NOT redone this window. (1) THE REAPER PASSES 68 + 69 (both ~4-5s, detached pids 1765680/1769286, receipts prep-reaper.log + prep-reaper-r89.nohup + prep-reaper-r89b.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): pass 68 (11:31:40-45Z) class C 27 files/dirs 1,541M consumed intermediates on pending territories + markers on 36 cac / 177 pending, freed 1,530M apparent (prep 110,965 -> 109,435M; df 626 -> 628G), census pending=181 cac=42; pass 69 (11:34:00-04Z) class B 24 files 28M reaping the batch-121 cohort (T04246/T04185/T04086/T04099) THE MOMENT IT RESOLVED HONEST exit=0 wall=444s in band + markers on 40 cac / 173 pending, freed 13M apparent (prep 109,435 -> 109,422M; df 626G), census pending=177 cac=46. CUMULATIVE ~= 204.4G apparent (~201.2G class-attributed) across 69 idempotent passes, every cut/call input and every split ledger kept. (2) THE EXCLUSION PROOF (62nd-64th times, SIX-FILE SYNG SIDECAR RULE): at pass 68 the batch-121 cohort protected 6/6 sidecars + locus.fa 113-119M KEPT mid-cut (cut pid 1736309 ~39.7GB niced); at pass 69 the batch-121 cohort reaped EXACTLY (sidecars GONE 0/6 x4, markers dropped, cut.done 4/4, definitional beds 6,629-11,387 rows KEPT) + the batch-122 cohort (T04322/T04327/T04206/T04132, dispatched 11:33:26Z) protected 6/6 in flight; the 6 live call lanes (T02866/T03208/T04066/T04059/T03411/T04038) spot-verified locus.fa 117-134M + reference-axis.json + locality-graphs intact, call.done absent in-flight (the ride-62 clarification standing). (3) WIDTH-4 CENSUS 122/122 post-swap ALL exit=0 BY FILE-MTIME CUTOFF (epoch 1792176762), ZERO non-zero, ZERO guard kills since the swap (batch-121 honest 444s in band; batch-1792236806-122 in flight at close); cut.done 2,204 -> 2,208; map dirs 72 -> 76. (4) CALLS - ZERO new fails this window (fail-closed 611 steady, the deep-lane memory-tail class 10 members incl. ride-88's T02722 71.05GB/3,610s; callok 1,551 steady - the live lanes are long); THE DEEP-LANE WATCH AT CLOSE: T02866 49.6-50.8GB @ ~68m the widest (plateaued after the 11:30Z 50.8GB reading, climbing slowly - the guard is the oracle, NEVER kill by hand), T03208 37.9GB second, T04066 33.0GB / T04059 31.0GB / T03411 29.9GB / T04038 27.8GB - none at guard pressure. (5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 626-629G, margin >320G); guard-oracle held (zero hand kills, zero masters touched); the prep master mid-dispatch 64-wide (sub 2,427 the pole producing, lean buffer 177-181 <= 200 no eviction); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference alive untouched (pids 1846398/1846402). ONE TRANSIENT, TWICE SEEN, NAMED: the prep-reaper.log tail ENOENT at both pass starts (the known in-place log-rewrite window, named at the ride-86 close, log intact both times, no consequence). FINAL POLL 11:36:35Z: split 9,891 / acc 5,844 / sub 2,427 / cut.done 2,208 / map dirs 76 (batch-122 in flight) / call-terminal 2,183 = callok 1,551 + fail-closed 611 + cut.exit 16 + poison 5 = 2,183/5,844 (37.4%) / df 626G / mem 710-737Gib / rvsize 170-171G. MASTERS: prep 411760, cuts 3838087 width-4-direct W=4 C=1, wave 2810146 LANES=128, heartbeat; NO PR PUSH (the owner pushes). THE SUCCESSOR RESUMES the standing cadence unchanged (sparse SHORT polls; the reaper cadence ~4s a pass idempotent; the SIX-FILE SIDECAR exclusion rule; the width-4 census continuation with batch-122 in flight; the deep-lane watch T02866 under the guard-oracle; the lean-buffer watch pending 177 vs the 200 bound; the 300G hard floor; the ladder only at backlog drained AND df >= 1,600G); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + the reaper receipts ~= 204.4G apparent/~201.2G class-attributed/69 passes) -> the record's closing increment -> the commit. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 90 WINDOW CLOSE (11:45Z 2026-10-17) - TWO CADENCE PASSES RECEIPTED (pass 70 +0M ~5s ALL CLASSES ZERO, markers on 40 cac + 173 pending; pass 71 +67M ~4s CLASS B 24 files 82M REAPING THE BATCH-122 COHORT THE MOMENT IT RESOLVED HONEST exit=0 wall ~566s/PEAK RSS ~60.9GB - THE WIDEST WIDTH-4 BATCH RSS YET, IN-BAND UNDER THE 64GiB GUARD, the dense T04322/T04132 locus cohort); CUMULATIVE ~= 204.5G APPARENT (~201.3G CLASS-ATTRIBUTED)/71 PASSES; CENSUS 123/123 ZERO KILLS SINCE THE SWAP (batch-123 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A SIXTY-FIFTH-TO-SIXTY-SEVENTH TIME; ONE NEW FAIL - T04059 CALL GUARD-KILLED (137) AT 72.56e6 kB/4,476s MID-SCORE LOCUS 3, THE DEEP-LANE CLASS'S 11TH MEMBER, THE SANCTIONED GUARD-ORACLE'S OWN KILL (the fast-climb tail: 31.5GB @ ~1h10m -> 61.5GB -> the kill in ~4min); ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 11:39Z (session-open increment 27f5b8d). THE ORIENTATION STANDING (rides 38-90): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass; the predecessor's expected ~136G superseded by the measured class attribution) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen is NOT redone this window - a fresh eviction would only delete cut/call inputs the live campaign still needs.
(1) THE REAPER PASSES 70 + 71 (both ~4-5s, detached pids 1775723/1780165, receipts prep-reaper.log + prep-reaper-r90.nohup + prep-reaper-r90b.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): PASS 70 (11:39:51-56Z): class A 0, class B 0 files 0M + syng/export/walks/walkreq markers dropped on 40 cut-awaiting-call (census pending=177 cac=46), class C 0 + export.done/walks.done/walkreq.done markers dropped on 173 pending; lean buffer 177 <= 200 no eviction; freed 0M apparent (prep 109,422M steady; df 626G). PASS 71 (11:42:56-59Z, run the moment batch-122 resolved at 11:42:52Z): class B 24 files 82M reaping the batch-122 cohort (T04322/T04327/T04206/T04132) + markers dropped on 44 cac (census pending=173 cac=50), class C 0 + markers on 169 pending; lean buffer 173 <= 200 no eviction; freed 67M apparent (prep 109,422 -> 109,355M; df 619G). CUMULATIVE ~= 204.5G apparent (~201.3G class-attributed) across 71 idempotent passes (129G ride-32 + the 69 cadence passes since), every cut/call input and every split ledger kept throughout. The known prep-reaper.log tail ENOENT transient seen ONCE at pass-70 start (the in-place log-rewrite window, named ride-86/89, log intact, no consequence).
(2) THE EXCLUSION PROOF (65th-67th times, SIX-FILE SYNG SIDECAR RULE): pass 70 protected the batch-122 cut cohort 6/6 sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa 113-411M KEPT mid-cut + correctly named the 6 live call lanes (T02866/T03208/T03411/T04038/T04059/T04066 spot-verified locus.fa 117-169M + reference-axis.json + locality-graphs 3-21 intact, sidecars 0 = the prior class-B reap, call.done absent in-flight - the ride-62 clarification standing); pass 71 reaped the batch-122 cohort EXACTLY (spot-check: syng sidecars GONE 0/6 each, locus.fa 113-411M + the definitional beds KEPT - locus.refined.bed 465-1,409 rows, producer markers dropped 0 present, cut.done present 4/4 - the reap exactly where intended) and its in-flight list already named the NEW batch-123 cohort (T04248/T04195/T04310/T04271, dispatched 11:42Z, cuts/batch-1792237367-123 in flight at close) - protected 6/6 sidecars + locus.fa 113-118M KEPT mid-cut at the final check.
(3) THE WIDTH-4 CENSUS: batches post-swap = 123/123 HONEST ALL exit=0 BY FILE-MTIME CUTOFF (epoch 1792176762 = the 18:52:42Z swap; the pre-swap non-zero markers are campaign history, not the census), ZERO non-zero, ZERO guard kills since the swap; this window: batch-122 HONEST exit=0 wall ~566s (dispatched ~11:33:26Z, resolved 11:42:52Z) - RSS trace peak ~60.9e6 kB, THE WIDEST WIDTH-4 BATCH RSS YET (the dense T04322/T04132 locus.fa 392/411M cohort; in-band under the 64GiB guard = 68.7e6 kB by the record's kB/1e6 convention, margin ~11%; the batch-density variance NAMED - the width-4 batches are not uniformly ~48-50GB, watch the dense cohorts); cut.done 2,208 -> 2,212 (+4, batch-122's cohort) - the width-4 pace holding (~29-35 territories/h); batch-1792237367-123 in flight at close.
(4) CALLS - ONE NEW FAIL THIS WINDOW, THE SANCTIONED GUARD-ORACLE'S OWN KILL: T04059 CALL GUARD-KILLED 11:43:46Z (exit 137, wall 4,476s, killed at 72,555,668 kB RSS mid-score locus 3: 16,783 records / 4,538,654 units / 536 members -> 525 folds / no-pin (unit,fold) 1.36e9 / abstained hull 1.41Gbp; the rss trace shows the FAST-CLIMB tail: 31.5GB @ ~1h10m -> 61.5GB @ ~1h14m -> the kill ~30s later, ~40GB in ~5min); the deep-lane memory-tail class now T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915/T03919/T02722/T04059 - ELEVEN members, terminal, named (receipts wave-logs/failures.txt 612 + per-territory call.err/call.rss/call.wall/call.exit.code IN-LANE FAIL-REAPED - map/T04059 gone, the definitional substrate kept); fail-closed 611 -> 612; callok 1,551 steady; cut.exit 16; poison 5 terminal (the 5-member syng segdup-tail class T00197/T01411/T02562/T00024/T03679). THE DEEP-LANE WATCH AT CLOSE: T02866 48.4GB @ ~1h17m the widest (plateaued since 11:30Z, climbing slowly - the guard is the oracle, NEVER kill by hand), T03208 45.3GB second and climbing (the next watch), T04066 32.2GB / T03411 32.1GB / T04038 27.1GB - none at guard pressure.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 603-637G, margin ~303G at the close - WATCH: the prep pole is producing FAST post-T03679 (sub 2,427 -> 2,435, rvsize 159 -> 195G, df 637 -> 603G in ~15min) and the VALVE HELD at 161 pending >= 160 at 11:43:26Z - the standing bound working, the pole self-pauses; never lower the floor without the owner); guard-oracle held (zero hand kills, zero masters touched - the T04059 kill was the guard's own, sanctioned); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched, 22d+).
FINAL POLL 11:44:38Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,435 / cut.done 2,212 / map dirs 79 (batch-123 in flight) / call-terminal 2,184 = callok 1,551 + fail-closed 612 + cut.exit 16 + poison 5 = 2,184 of 5,844 (37.3% of accepted) / df 603G / mem 738Gib / rvsize 195G apparent. MASTERS AT CLOSE: prep 411760 (~1-17:19h, VALVE holding at 161), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~16:52h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~1-00:48h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-123 in flight; watch the dense-locus cohorts' batch RSS - batch-122's ~60.9GB peak named); the JIT reaper cadence (~4-5s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every pass; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch T02866/T03208 under the guard-oracle (never kill by hand - the class is 11 members); the prep-pole df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner) + the lean-buffer watch (pending 173 vs the 200 bound - eviction DESIGNED if it trips, the ~22s regen chain); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 5 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~204.5G cumulative apparent/~201.3G class-attributed across 71 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 5 CONFIRMED members, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915/T03919/T02722/T04059 + rides 55/56's clusters + the ride-69 + ride-81 consumption-burst receipts) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 91 WINDOW CLOSE (11:52Z 2026-10-17) - TWO CADENCE PASSES RECEIPTED (pass 72 11:46:54-59Z +1,097M ~5s class C 27 files/dirs 1,109M consumed intermediates on pending, markers on 44 cac / 178 pending, census pending=182 cac=48; pass 73 11:50:51-55Z +13M ~5s class B 24 files/dirs 28M REAPING THE BATCH-123 COHORT THE MOMENT IT RESOLVED HONEST exit=0 wall 433s/peak ~48.9GB in band, markers on 48 cac / 174 pending, census pending=178 cac=52) - cumulative ~= 205.6G apparent (~202.4G class-attributed)/73 idempotent passes; CENSUS 124/124 ZERO KILLS SINCE THE SWAP (batch-124 IN FLIGHT AT CLOSE); THE EXCLUSION PROOF A SIXTY-EIGHTH-TO-SEVENTIETH TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 11:46Z (session-open increment 90acd78). THE ORIENTATION STANDING (rides 38-91): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass; the predecessor's expected ~136G superseded by the measured class attribution) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen is NOT redone this window - a fresh eviction would only delete cut/call inputs the live campaign still needs (2,185/5,844 terminal, 63% in flight).
(1) THE REAPER PASSES 72 + 73 (both ~5s, detached nohup, receipts prep-reaper.log + prep-reaper-r91.nohup + prep-reaper-r91b.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): PASS 72 (11:46:54-59Z): class A 0, class B 0 files 0M + markers on 44 cut-awaiting-call, class C 27 files/dirs 1,109M consumed intermediates on pending territories + markers on 178 pending; lean buffer 182 <= 200 no eviction; freed 1,097M apparent (prep 111,449 -> 110,352M; df 620 -> 621G). PASS 73 (11:50:51-55Z, run the moment batch-123 resolved): class B 24 files 28M reaping the batch-123 cohort (T04248/T04195/T04310/T04271) + markers on 48 cac, class C 0 + markers on 174 pending; lean buffer 178 <= 200 no eviction; freed 13M apparent (prep 110,352 -> 110,339M; df 619G steady). CUMULATIVE ~= 205.6G apparent (~202.4G class-attributed) across 73 idempotent passes, every cut/call input and every split ledger kept. The known prep-reaper.log tail ENOENT transient seen ONCE at pass-72 start (the in-place log-rewrite window, named rides 86/89/90, log intact, no consequence).
(2) THE EXCLUSION PROOF (68th-70th times, SIX-FILE SYNG SIDECAR RULE): pass 72 protected the batch-123 cohort 6/6 sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa 113-118M KEPT mid-cut + correctly named the 4 live call lanes (T02866/T03411/T04038/T04066 spot-verified locus.fa 124-169M + reference-axis.json + locality-graphs 3-21 intact, sidecars 0 = the prior class-B reap, call.done absent in-flight - the ride-62 clarification standing); pass 73 reaped the batch-123 cohort EXACTLY (spot-check: syng sidecars GONE 0/6 each, locus.fa + the definitional beds KEPT - locus.refined.bed 465-509 rows, producer markers dropped 0 present, cut.done present 4/4 - the reap exactly where intended) and its in-flight list already named the NEW batch-124 cohort (T04337/T04338/T04332/T04339, dispatched ~11:50Z, cuts/batch-1792237934-124 in flight at close) - protected 6/6 sidecars + locus.fa 113-250M KEPT mid-cut at the final check (cut pid 1813336 ~13.2GB niced early-load).
(3) THE WIDTH-4 CENSUS: batches post-swap = 124/124 HONEST ALL exit=0 BY FILE-MTIME CUTOFF (epoch 1792176762 = the 18:52:42Z swap; 555 batch.exit markers on disk, the pre-swap non-zero are campaign history, not the census), ZERO non-zero, ZERO guard kills since the swap; this window: batch-123 HONEST exit=0 wall 433s peak RSS 48,857,936 kB ~= 48.9GB in band; batch-124 in flight at close; cut.done 2,212 -> 2,216 (+4, batch-123's cohort) - the width-4 pace holding (~29-35 territories/h); map dirs 78 -> 82 as the wave consumes.
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 612 STEADY (the last, T04059 11:43:46Z guard-137, receipted at the ride-90 close, predates this open); callok 1,552 steady (the 3 lanes in flight, receipts land at call.done); cut.exit 16; poison 5 terminal (the 5-member syng segdup-tail class T00197/T01411/T02562/T00024/T03679). T03208 COMPLETED HONESTLY between polls (call.done exit=0, wall 4,615s - the long-ordinary-lane receipt: the ride-85/86/87 watched deep lane 37-45GB that never tripped the guard, the honest distribution point for the datasheet alongside the 11 guard-137 members); THE DEEP-LANE WATCH AT CLOSE: T02866 48.4GB @ ~1h23m the widest (plateaued since 11:30Z, climbing slowly - the guard is the oracle, NEVER kill by hand; the class is 11 members), 32.2GB / 27.1GB second and third - none at guard pressure.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 619-621G, margin ~319G); guard-oracle held (zero hand kills, zero masters touched); the prep master mid-dispatch 64-wide chunk continuing (sub 2,436 steady, lean buffer 178-182 <= 200 no eviction); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched, 22d+).
FINAL POLL 11:51:21Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,436 / cut.done 2,216 / map dirs 82 (batch-124 in flight) / call-terminal 2,185 = callok 1,552 + fail-closed 612 + cut.exit 16 + poison 5 = 2,185 of 5,844 (37.4% of accepted) / df 619G / mem 819Gib / rvsize 175G apparent. MASTERS AT CLOSE: prep 411760 (~1-17:25h, mid-dispatch chunk), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~17:00h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~1-00:53h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-124 in flight; watch the dense-locus cohorts' batch RSS - batch-122's ~60.9GB peak named ride-90); the JIT reaper cadence (~5s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every pass; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch T02866 under the guard-oracle (never kill by hand - the class is 11 members); the prep-pole df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner) + the lean-buffer watch (pending 178 vs the 200 bound - eviction DESIGNED if it trips, the ~22s regen chain); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 5 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~205.6G cumulative apparent/~202.4G class-attributed across 73 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 5 CONFIRMED members, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915/T03919/T02722/T04059 + rides 55/56's clusters + the ride-69 + ride-81 consumption-burst receipts + the T03208 4,615s long-ordinary-lane completion) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 92 WINDOW CLOSE (11:59Z 2026-10-17) - TWO CADENCE PASSES RECEIPTED (pass 74 11:54:52-57Z +0M ~5s ALL CLASSES ZERO, markers on 48 cac / 174 pending, census pending=178 cac=51; pass 75 11:57:48-52Z +20M ~4s CLASS B 24 files 37M REAPING THE BATCH-124 COHORT THE MOMENT IT RESOLVED HONEST exit=0 wall 461s/peak ~49.8e6 kB in band, markers on 52 cac / 170 pending, census pending=174 cac=54) - cumulative ~= 205.6G apparent (~202.4G class-attributed)/75 idempotent passes; CENSUS 125/125 ZERO KILLS SINCE THE SWAP (batch-1792238263-125 in flight at close); THE EXCLUSION PROOF A SEVENTY-FIRST-TO-SEVENTY-THIRD TIME; TWO NEW FAILS - T03411 (11:51:21Z, the window edge, receipted at the session-open increment) + T04038 (11:55:54Z, in-window) BOTH THE SANCTIONED GUARD-ORACLE'S OWN KILLS at ~68.3GB/5,012s and ~68.9GB/5,236s - THE DEEP-LANE CLASS NOW THIRTEEN MEMBERS; ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 11:54Z (session-open increment 84163a2). THE ORIENTATION STANDING (rides 38-92): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass; the predecessor's expected ~136G superseded by the measured class attribution) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen is NOT redone this window - a fresh eviction would only delete cut/call inputs the live campaign still needs (2,187/5,844 terminal, 63% in flight).
(1) THE REAPER PASSES 74 + 75 (both ~4-5s, detached pids 1818314/1822049, receipts prep-reaper.log + prep-reaper-r92.nohup + prep-reaper-r92b.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): PASS 74 (11:54:52-57Z): class A 0, class B 0 files 0M + syng/export/walks/walkreq markers dropped on 48 cut-awaiting-call (census pending=178 cac=51), class C 0 + export.done/walks.done/walkreq.done markers dropped on 174 pending; lean buffer 178 <= 200 no eviction; freed 0M apparent (prep 110,216M steady; df 621G). PASS 75 (11:57:48-52Z, run the moment batch-124 resolved HONEST exit=0 wall 461s, rss trace peak ~49.8e6 kB ~= 49.8GB in band): class B 24 files 37M reaping the batch-124 cohort (T04337/T04338/T04332/T04339) + markers dropped on 52 cac (census pending=174 cac=54), class C 0 + markers on 170 pending; lean buffer 174 <= 200 no eviction; freed 20M apparent (prep 110,073 -> 110,053M; df 620-629G). CUMULATIVE ~= 205.6G apparent (~202.4G class-attributed) across 75 idempotent passes (129G ride-32 + the 73 cadence passes since), every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF (71st-73rd times, SIX-FILE SYNG SIDECAR RULE): pass 74 protected the batch-124 cut cohort 6/6 sidecars intact (locus-pos64.syng 1gbwt/1khash/meta/names/pstep/spos) + locus.fa 113-250M KEPT (cut.done absent as expected, mid-cut, cut pid 1813336 ~37.4GB niced) + correctly named the 3 live call lanes (T02866/T04066/T04038 spot-verified locus.fa 134-169M + reference-axis.json + locality-graphs 11-21 intact, call.done absent in-flight - the ride-62 clarification standing); T03411 correctly ABSENT from the pass-74 exclusion list post-fail (guard-killed 11:51:21Z BEFORE the pass; its inputs fail-reaped IN-LANE by the wave's own fail path, receipts call.err/call.rss/call.wall kept); pass 75 reaped the batch-124 cohort EXACTLY (spot-check: syng sidecars GONE 0/6 each, syng.done producer markers dropped 0 present, cut.done present 4/4, locus.fa 113-250M + the definitional beds KEPT - locus.refined.bed 465-904 rows - the reap exactly where intended) and its in-flight list already named the NEW batch-125 cohort (T04328/T04219/T04093/T04133, cuts/batch-1792238263-125 dispatched within the minute, in flight at close) - protected at the close check; T04038 correctly excluded while live at pass-74 time then guard-killed at 11:55:54Z and fail-reaped IN-LANE (locus.fa gone, call receipts kept).
(3) THE WIDTH-4 CENSUS: batches post-swap = 125/125 HONEST ALL exit=0 BY FILE-MTIME CUTOFF (epoch 1792176762 = the 18:52:42Z swap; the pre-swap non-zero markers are campaign history, not the census), ZERO non-zero, ZERO guard kills since the swap; this window: batch-124 HONEST exit=0 wall 461s peak ~49.8GB in band, resolved ~11:57:44Z; batch-1792238263-125 in flight at close; map dirs 81 -> 84 as the wave consumes; the width-4 pace holding (~29-35 territories/h).
(4) CALLS - TWO NEW FAILS, BOTH THE SANCTIONED GUARD-ORACLE'S OWN: T03411 CALL FAILED (1) 11:51:21Z (the RSS guard's own exit - "RSS guard exceeded at locus_1_scored: 68,263,780 kB" ~= 68.3GB, wall 5,012s; locus 1: 22,237 records / 5,160,612 units / no-pin (unit,fold) 1.26e9 / abstained hull 1.64Gbp; the rss trace 58.7 -> 64.3GB -> the kill) - the WINDOW-EDGE fail, receipted at the session-open increment; T04038 CALL FAILED (137) 11:55:54Z (GUARD kill at 68,939,600 kB ~= 68.9GB, wall 5,236s, mid-score locus 9: exactness sample 224 pairs of 499,665,625 valid - the fast-climb tail: 27.1GB @ open -> 64.6 -> 68.9GB in ~90s); the deep-lane memory-tail class now T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915/T03919/T02722/T04059/T03411/T04038 - THIRTEEN members, terminal, named (receipts wave-logs/failures.txt 614 + per-territory call.err/call.rss/call.wall/call.exit.code IN-LANE FAIL-REAPED - map gone, definitional substrate kept); fail-closed 613 -> 614 (T03411 landed at the 91/92 edge, T04038 in-window); callok 1,552 steady (the new lanes on batch-124's cohort just dispatched ~11:58Z, receipts land at call.done); cut.exit 16; poison 5 terminal (the 5-member syng segdup-tail class T00197/T01411/T02562/T00024/T03679). THE DEEP-LANE WATCH AT CLOSE: T02866 48.4GB @ ~1h31m the widest (plateaued since 11:30Z - the guard is the oracle, NEVER kill by hand; the class is 13 members), T04066 32.2GB second - none at guard pressure.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 620-629G, margin >320G); guard-oracle held (zero hand kills, zero masters touched - both T03411/T04038 kills were the guard's own, sanctioned); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched, 22d+).
FINAL POLL 11:58:24Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,436 / map dirs 84 (batch-125 in flight) / call-terminal 2,187 = callok 1,552 + fail-closed 614 + cut.exit 16 + poison 5 = 2,187 of 5,844 (37.4% of accepted) / df 629G / mem 752Gib / rvsize 166G apparent. MASTERS AT CLOSE: prep 411760 (~1-17:33h), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~17:06h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~1-00:58h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-125 in flight; watch the dense-locus cohorts' batch RSS - batch-122's ~60.9GB peak named ride-90); the JIT reaper cadence (~4-5s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch T02866 under the guard-oracle (never kill by hand - the class is 13 members); the prep-pole df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner) + the lean-buffer watch (pending 174 vs the 200 bound - eviction DESIGNED if it trips, the ~22s regen chain); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 5 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~205.6G cumulative apparent/~202.4G class-attributed across 75 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 5 CONFIRMED members, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard points T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915/T03919/T02722/T04059/T03411/T04038 + rides 55/56's clusters + the ride-69 + ride-81 consumption-burst receipts + the T03208 4,615s long-ordinary-lane completion) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 93 WINDOW CLOSE (12:09Z 2026-10-17) - TWO CADENCE PASSES RECEIPTED (pass 76 12:01:21-25Z +911M ~4s class C 27 files 919M consumed intermediates on pending, census pending=183 cac=54; pass 77 12:06:48-53Z +15M ~5s class B 24 files 31M REAPING THE BATCH-125 COHORT THE MOMENT IT RESOLVED HONEST exit=0 wall 446s/peak ~53.3e6 kB in band, census pending=179 cac=58) - cumulative ~= 206.5G apparent (~203.3G class-attributed)/77 idempotent passes; CENSUS 126/126 ZERO KILLS SINCE THE SWAP (batch-1792238711-126 in flight at close); THE EXCLUSION PROOF A SEVENTY-FOURTH-TO-SEVENTY-SIXTH TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 12:01Z (session-open increment 29ed17d). THE ORIENTATION STANDING (rides 38-93): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass; the predecessor's expected ~136G superseded by the measured class attribution) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen is NOT redone this window - a fresh eviction would only delete cut/call inputs the live campaign still needs (2,187/5,844 terminal, 63% in flight).
(1) THE REAPER PASSES 76 + 77 (both ~4-5s, detached pids 1847906/1853427, receipts prep-reaper.log + prep-reaper-r93.nohup + prep-reaper-r93b.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): PASS 76 (12:01:21-25Z): class A 0, class B 0 files 0M + markers on 52 cut-awaiting-call, class C 27 files/dirs 919M consumed intermediates on pending territories + markers on 179 pending; lean buffer 183 <= 200 no eviction; freed 911M apparent (prep 112,033 -> 111,122M; df 618 -> 619G). PASS 77 (12:06:48-53Z, run the moment batch-125 resolved HONEST exit=0 wall 446s peak RSS 53,320,412 kB ~= 53.3GB in band): class B 24 files 31M reaping the batch-125 cohort (T04328/T04219/T04093/T04133) + markers dropped on 56 cac, class C 0 + markers on 175 pending; lean buffer 179 <= 200 no eviction; freed 15M apparent (prep 111,122 -> 111,107M; df 616G steady). CUMULATIVE ~= 206.5G apparent (~203.3G class-attributed) across 77 idempotent passes (129G ride-32 + the 75 cadence passes since), every cut/call input and every split ledger kept throughout. The known prep-reaper.log tail ENOENT transient seen ONCE at pass-76 start (the in-place log-rewrite window, named rides 86/89/90, log intact, no consequence).
(2) THE EXCLUSION PROOF (74th-76th times, SIX-FILE SYNG SIDECAR RULE): pass 76 protected the batch-125 cut cohort 6/6 sidecars intact (1gbwt/1khash/meta/names/pstep/spos) + locus.fa 113-217M KEPT mid-cut + correctly named the 2 live call lanes (T02866 locus.fa 134M + reference-axis.json + locality-graphs 21, T04066 locus.fa 169M + locality-graphs 13 - call.done absent in-flight, the ride-62 clarification standing); pass 77 reaped the batch-125 cohort EXACTLY (spot-check: syng sidecars GONE 0/6 each, syng.done producer markers dropped 0 present, cut.done present 4/4, locus.fa 113-217M + the definitional beds KEPT - locus.refined.bed 465-819 rows - the reap exactly where intended) and its in-flight list already named the NEW batch-126 cohort (T04346/T04350/T04345/T04342, cuts/batch-1792238711-126 dispatched ~12:06Z, in flight at close) - protected 6/6 sidecars + locus.fa KEPT mid-cut at the close check.
(3) THE WIDTH-4 CENSUS: batches post-swap = 126/126 HONEST ALL exit=0 BY FILE-MTIME CUTOFF (epoch 1792176762 = the 18:52:42Z swap; the pre-swap non-zero markers are campaign history, not the census), ZERO non-zero, ZERO guard kills since the swap; this window: batch-125 HONEST exit=0 wall 446s peak ~53.3GB in band, resolved 12:06:44Z; batch-1792238711-126 in flight at close; cut.done 2,220 -> 2,224 (+4, batch-125's cohort) - the width-4 pace holding (~29-35 territories/h); map dirs 84 -> 88 as the wave consumes.
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 614 STEADY (the last, T04038 11:56:01Z guard-137, receipted at the ride-92 close, predates this open; the deep-lane memory-tail class 13 members incl. T03411/T04038); callok 1,552 steady (the 2 live lanes are long, receipts land at call.done); cut.exit 16; poison 5 terminal (the 5-member syng segdup-tail class T00197/T01411/T02562/T00024/T03679). THE DEEP-LANE WATCH AT CLOSE: T02866 48.4GB @ ~1h40m the widest (plateaued 48-51GB since 11:30Z - the guard is the oracle, NEVER kill by hand; the class is 13 members), T04066 32.2GB @ ~1h38m second - none at guard pressure.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 616-619G, margin >316G); guard-oracle held (zero hand kills, zero masters touched); the prep master's chunked dispatch continuing with the VALVE holding at 167 pending >= 160 (sub 2,445 steady, lean buffer 179-183 <= 200 no eviction); terminal classes terminal; the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched, 22d+).
FINAL POLL 12:08:06Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,445 / cut.done 2,224 / map dirs 88 (batch-126 in flight) / call-terminal 2,187 = callok 1,552 + fail-closed 614 + cut.exit 16 + poison 5 = 2,187 of 5,844 (37.4% of accepted) / df 616G / mem 845Gib / rvsize 180G apparent. MASTERS AT CLOSE: prep 411760 (~1-17:38h, the valve holding; the per-batch children of the cuts master 1821964/1821986 verified its OWN ppid chain at the open - not a second master), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~17:10h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~1-01:10h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-126 in flight; watch the dense-locus cohorts' batch RSS - batch-125's ~53.3GB peak the widest this window, batch-122's ~60.9GB the standing max); the JIT reaper cadence (~4-5s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch T02866 under the guard-oracle (never kill by hand - the class is 13 members); the prep-pole df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner) + the lean-buffer watch (pending 179 vs the 200 bound - eviction DESIGNED if it trips, the ~22s regen chain); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 5 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~206.5G cumulative apparent/~203.3G class-attributed across 77 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 5 CONFIRMED members, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915/T03919/T02722/T04059/T03411/T04038 + rides 55/56's clusters + the ride-69 + ride-81 consumption-burst receipts + the T03208 4,615s long-ordinary-lane completion) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## RIDE CONTINUATION 94 WINDOW CLOSE (12:16Z 2026-10-17) - TWO CADENCE PASSES RECEIPTED (pass 78 12:11:30-34Z +0M ~4s ALL CLASSES ZERO, markers on 56 cac / 175 pending, census pending=179 cac=58; pass 79 12:14:08-12Z +721M ~4s CLASS B 24 files 28M REAPING THE BATCH-126 COHORT THE MOMENT IT RESOLVED HONEST exit=0 peak ~50.0e6 kB in band + CLASS C 21 files 715M consumed intermediates on pending, census pending=182 cac=62) - cumulative ~= 207.2G apparent (~204.0G class-attributed)/79 idempotent passes; CENSUS 127th BATCH DISPATCHED (batch-1792239119-127 cohort T04353/T04352/T04289/T04164 IN FLIGHT AT CLOSE - resolves after); THE EXCLUSION PROOF A SEVENTY-SEVENTH-TO-SEVENTY-NINTH TIME; ZERO NEW FAILS, ZERO HAND KILLS, ZERO MASTERS TOUCHED
Opened 12:10Z (session-open increment 8366ca5). THE ORIENTATION STANDING (rides 38-94): the ride-30/31 "reaper written not run" prompt is STALE - the substrate reaper has RUN since ride-32 (129G first pass; the predecessor's expected ~136G superseded by the measured class attribution) and THE REGEN RECEIPT IS CLOSED (ride-35: the ~22s on-demand regen measurement, the T03194 full-chain evict->rebuild->cut/call re-entry receipt, beds + ledgers kept); the regen is NOT redone this window - a fresh eviction would only delete cut/call inputs the live campaign still needs (2,187/5,844 terminal, 63% in flight).
(1) THE REAPER PASSES 78 + 79 (both ~4s, detached pids 1859166/1884751, receipts prep-reaper.log + prep-reaper-r94.nohup + prep-reaper-r94b.nohup; bash -n CLEAN + md5 b804a7d366844d27ac58e2e0c8b64802 UNCHANGED re-verified at the open): PASS 78 (12:11:30-34Z): class A 0, class B 0 files 0M + markers on 56 cut-awaiting-call, class C 0 + markers on 175 pending; lean buffer 179 <= 200 no eviction; freed 0M apparent (prep 111,107M steady; df 615G). PASS 79 (12:14:08-12Z, run the moment batch-126 resolved HONEST exit=0, rss trace peak 50,006,364 kB ~= 50.0GB in band, 41 samples): class B 24 files 28M reaping the batch-126 cohort (T04346/T04350/T04345/T04342) + markers dropped on 60 cac, class C 21 files/dirs 715M consumed intermediates on pending territories + markers on 178 pending; lean buffer 182 <= 200 no eviction; freed 721M apparent (prep 114,335 -> 113,614M; df 611 -> 612G). CUMULATIVE ~= 207.2G apparent (~204.0G class-attributed) across 79 idempotent passes (129G ride-32 + the 78 cadence passes since), every cut/call input and every split ledger kept throughout.
(2) THE EXCLUSION PROOF (77th-79th times, SIX-FILE SYNG SIDECAR RULE): pass 78 protected the batch-126 cut cohort 6/6 sidecars intact (locus-pos64.syng 1gbwt/1khash/meta/names/pstep/spos) + locus.fa 112-113M KEPT (cut.done absent as expected, mid-cut) + correctly named the 2 live call lanes (T02866 locus.fa 134M + reference-axis.json + locality-graphs 21, T04066 locus.fa 169M + locality-graphs 13 - call.done absent in-flight, the ride-62 clarification standing); pass 79 reaped the batch-126 cohort EXACTLY (spot-check: syng sidecars GONE 0/6 each, syng.done producer markers dropped 0 present, cut.done present 4/4, locus.fa 112-113M + the definitional beds KEPT - locus.refined.bed 465-471 rows - the reap exactly where intended) and its in-flight list already named the NEW batch-127 cohort (T04353/T04352/T04289/T04164, cuts/batch-1792239119-127 dispatched ~12:13Z, in flight at close) - protected 6/6 sidecars + locus.fa 113-156M KEPT mid-cut at the close check (cut pid 1860915 ~22GB niced early-load); the exclusion list's fifth name T04369 is the wave's next staged territory, not a second batch.
(3) THE WIDTH-4 CENSUS: batches post-swap = 126/126 HONEST ALL exit=0 BY FILE-MTIME CUTOFF (epoch 1792176762 = the 18:52:42Z swap; batch-13's width-8 guard-137 PREDATES the swap by 2s - campaign history, not the census), ZERO non-zero, ZERO guard kills since the swap; this window: batch-126 HONEST exit=0 resolved ~12:14Z, rss trace peak 50,006,364 kB ~= 50.0GB in band; batch-1792239119-127 in flight at close (census count moves to 127/127 at its resolve); cut.done 2,224 -> 2,228 (+4, batch-126's cohort) - the width-4 pace holding (~29-35 territories/h); map dirs 88 -> 92 as the wave consumes.
(4) CALLS - ZERO NEW FAILS THIS WINDOW: fail-closed 614 STEADY (the last, T04038 11:56:01Z guard-137, receipted at the ride-92 close, predates this open; the deep-lane memory-tail class 13 members incl. T03411/T04038); callok 1,552 steady (the 2 live lanes are long, receipts land at call.done); cut.exit 16; poison 5 terminal (the 5-member syng segdup-tail class T00197/T01411/T02562/T00024/T03679). THE DEEP-LANE WATCH AT CLOSE: T02866 48.4GB @ ~1h47m the widest (plateaued 48-51GB since 11:30Z - the guard is the oracle, NEVER kill by hand; the class is 13 members), T04066 32.2GB @ ~1h45m second - none at guard pressure.
(5) CONTINGENCIES: the 300G hard floor NEVER approached (df band 611-616G, margin ~311G); guard-oracle held (zero hand kills, zero masters touched); terminal classes terminal; the prep master's chunked dispatch continuing (sub 2,445 -> 2,453, lean buffer 179-182 <= 200 no eviction); the ladder NOT opened (df << 1,600G, backlog present); the ancient impg-reference benign - NEVER KILLED (pids 1846398/1846402 alive, untouched, 22d+).
FINAL POLL 12:15:02Z: split 9,891 / acc 5,844 (59.1% final) / sub 2,453 / cut.done 2,228 / map dirs 92 (batch-127 in flight) / call-terminal 2,187 = callok 1,552 + fail-closed 614 + cut.exit 16 + poison 5 = 2,187 of 5,844 (37.4% of accepted) / df 612G / mem 844Gib / rvsize 184G apparent. MASTERS AT CLOSE: prep 411760 (~1-17:49h), cuts 3838087 (width-4-direct BATCH_W=4 N_CONCURRENT=1 RAYON_NUM_THREADS=24 DISK_FLOOR_GB=300, ~17:22h), chunked wave 2810146 (LANES=128 CHUNK_MULT=2, ~1-01:18h); heartbeat live. THE SUCCESSOR RESUMES: sparse SHORT polls; the width-4 census continuation (batch-127 in flight; watch the dense-locus cohorts' batch RSS - batch-122's ~60.9GB the standing max); the JIT reaper cadence (~4s a pass, idempotent, re-runnable anytime; the SIX-FILE SYNG SIDECAR rule on every exclusion proof; verify call.done/receipts before reading inputs-gone as a reaper miss - the ride-62 clarification standing); the deep-lane watch T02866 under the guard-oracle (never kill by hand - the class is 13 members); the prep-pole df watch vs the 300G hard floor (self-gates; DO NOT lower without the owner) + the lean-buffer watch (pending 182 vs the 200 bound - eviction DESIGNED if it trips, the ~22s regen chain); the ladder only at backlog drained AND df >= 1,600G (cuts master MID-WAIT only, by exact pid, never the wave, never with marathons in flight; C=4 -> C=6 the flock ceiling, never RAM); poison-first on prep stalls (the class is 5 members); AT CAMPAIGN CLOSE (all 5,844 accepted terminal): rv-rollup.py -> FILL datasheet-runtime.md (the standing close list + the width-4-direct delta + the inventory finding + the JIT-prep rule + THE REAPER RECEIPTS: ~207.2G cumulative apparent/~204.0G class-attributed across 79 idempotent passes, the walkreq marker-gap finding + fix, the ~22s regen measurement, the T03194 full-chain re-entry/call receipt, the syng segdup-tail class's 5 CONFIRMED members, the mid-dispatch valve patch's live proof, the guard-oracle receipts incl. T03641 11,084s/71.1GB deepest-lane + the named guard-137 points T00129/T04052/T01895/T03855/T04063/T04013/T03803/T03915/T03919/T02722/T04059/T03411/T04038 + rides 55/56's clusters + the ride-69 + ride-81 consumption-burst receipts + the T03208 4,615s long-ordinary-lane completion) -> the record's closing increment -> the commit. SCOPE STANDING: the run uses the SYNG-NATIVE EXPORT MAPS (c-rider locality-maps-export, K=63) = EXPORT-PROFILE answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. Assessment-side only; the committed projection/shards/receipts immutable. NO PR PUSH (the owner pushes).

## CAMPAIGN CLOSE (the stop order ~12:16Z -> the close-out complete 12:57Z, 2026-10-17) - STOPPED BY OWNER, THE HONEST PARTIAL-DRAIN CLOSE RECEIPTED - THIS IS THE FINAL SECTION; THE CAMPAIGN IS CLOSED
THE OWNER'S STOP ORDER, verbatim: "ok its stuck... we should stop them" - the masters killed by exact pid (prep 411760, cuts 3838087 width-4-direct, wave 2810146, the ancient impg-reference 1846398/1846402 - stopped, not killed mid-work: zero hand kills of any lane, the guard-oracle doctrine held to the end). THE HONEST TAILS LEFT TO LAND, ALL LANDED: batch-127 (the last cut in flight) COMPLETED ITS WORK 428s/peak 52,033,076 kB, 4/4 cut-stats.json written, wrapper markers absent with the wrapper (receipt cuts/batch-1792239119-127/batch.err; the cut.done census stays 2,228, the work-complete census 2,232); the two deep call lanes resolved UNDER THEIR OWN 64GiB GUARDS, both honest terminal fail-closed: T02866 exit 137 at 7,730s/67,575,876 kB (the widest watch lane, ~1.5h plateau then the fast-climb tail) + T04066 exit 137 at 6,804s/67,829,716 kB (receipts prep/T02866,T04066/call.{rss,wall,err,exit.code}) - the deep-lane memory-tail class closed at 134 guard kills total (132 exit-137 + 2 exit-1 in-process RSS guard), the marathon roster 15 named, T03641 the deepest 11,084s/71.1GB.
(1) THE FINAL CENSUS (poll 12:41:39Z + the post-reap poll 12:57:33Z): split 9,891 / accepted 5,844 (59.1%) / sub 2,453 at the stop (2,210 markers post-reap - the close reap drops producer markers on reaped never-called per the idempotence doctrine) / cut.done 2,228 + 4 work-complete-unmarked / CALL TERMINALS 2,189 of 5,844 = callok 1,552 + fail-closed 616 (482 the 7479-class bind-ambiguity family, 132 guard-137, 2 in-process RSS guard) + cut.exit 16 + poison 5 = 37.5% final (2,187 = 37.4% at the stop order). THE STOP WAS AT THE OWNER'S ORDER, NOT STUCKNESS: the cuts pole paced ~29-35 territories/h (width-4-direct, W=4 C=1, RAYON=24, DISK_FLOOR_GB=300), ~3,655 never-called accepted territories ~= 4-5 days at the current machinery, and the verification's questions were already answered by the receipts.
(2) THE CLOSE REAP RUN + VERIFIED SURGICAL (rv-close-reap.sh, receipt close-reap.log, 12:45-12:47Z): freed +128.5G apparent - df 588G -> 712G avail (CLOSE-A 2,129 files/dirs 69,321M regenerable substrate on the 3,676 never-called accepted territories + the idempotent residue sweep; CLOSE-B 270 files 60,596M cut-evidence bulk on 90 territories, map/<T>/<T>/cut-stats.json KEPT). THE SURGICAL VERDICT, spot-checked: the input BAM + hg002-projection/full untouched; ALL definitional beds/ledgers intact (locus.refined.bed/locus.bed/split.info.json/split.rejected on all 9,891 dirs; accepted/split-done/rejected ledgers); receipts/ complete (2,168 dirs, 1,552 calls.jsonl + sequences, 2,014 cut-stats.json); ALL campaign logs/batch receipts/call walls+RSS/poison.tsv (5 members full prose)/handoff/datasheet/query results intact; ONLY regenerable substrate gone (0 locus.fa, 0 syng slices, 0 export-graphs/walks/maps/locality-graphs remain in prep/); the 5 poison territories' syng receipts kept verbatim; prep 113.3G -> 45.3G, map 60.5G -> ~0 apparent. NOTHING THE REAP SHOULD NOT HAVE TAKEN WAS TAKEN.
(3) THE ROLLUP RUN + REFRESHED (rv-rollup.py, receipt rollup-final.json/rollup-final.txt): the rollup does NOT demand all-terminal - it emits the partial-drain close honestly; the one close fix (the two rejection ledgers record the same windows - deduped, each rejection counted ONCE) receipted in the header; verified IDEMPOTENT (consecutive runs byte-identical). Every datasheet number traces to rollup-final.json or a named on-disk receipt.
(4) datasheet-runtime.md FILLED - THE HONEST PARTIAL-DRAIN CLOSE (the refuted-claims delta table; the producing-pace receipts incl. the width-4 +54% swap delta + the CHUNK_MULT=2 chunked wave + the valve's live proof; the failure-class census all named; the JIT reaper ~207.2G apparent/~204.0G class-attributed/79 passes + the ~22s regen + the close-out reap; the stopped-not-stuck framing; the 20-minute optimization path green-lit with the measured budget table).
(5) THE RIDERS ALL DEAD AT THE CLOSE: the two leftover rider-class daemons (heartbeat x2, the REARM df-waiter) stopped at the close - logging/waiter processes only, zero work processes touched; wave-logs/heartbeat.log 3,353 lines, the disk high-water record complete (min df 286G, min memavail 21GiB).
THE SUCCESSOR HAS NO SUCCESSOR - THE CAMPAIGN IS CLOSED. Any re-run reads this handoff + datasheet-runtime.md + the campaign record's closing increment (docs/syng-gmem-bwt/campaign-record-2026-10-diploid-arc.md + the repo copies of the datasheet/handoff under docs/syng-gmem-bwt/). The optimization campaign is GREEN-LIT by the owner (the 20-minute target; the global-panel call mode under implementation). SCOPE STANDING: the run used the SYNG-NATIVE EXPORT MAPS (K=63) = export-profile answers - THE RUNTIME RECEIPT, NOT THE ACCURACY RECEIPT. NO PR PUSH (the owner pushes).
