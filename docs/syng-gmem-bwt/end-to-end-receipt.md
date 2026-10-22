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

**[PENDING]** the census cut, the routes inventory (wall, core-s, the
identity comparison), the per-territory prep + panel-prep walls and
core-s sums.

## 2. STEP 2 receipts — map

**[PENDING]** the projection wall/core-s/throughput vs the committed
20.79min byte-identity receipt; the identity gate; the per-territory
map leg walls/core-s.

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
