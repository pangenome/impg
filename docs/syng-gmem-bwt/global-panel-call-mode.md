# THE GLOBAL-PANEL CALL MODE — lever 1 of the 20-minute whole-genome budget

> **THE MANDATE (the owner's green light):** "we should be calling whole
> genome in 20 minutes or less." The measured budget from the closed
> runtime-verification campaign: map ~3 min (the 12-lane projection
> projection basis: 4.74M reads/s aggregate, byte-identical outputs),
> prep ZERO (the global panel), cuts + calls ~10-15 min, the marathon
> tail bounded (lever 4, the cross_lcp_merge segdup core). This
> document is lever 1's design of record: **the global-panel call
> mode** — cut and call read the ONE committed whole-panel syng
> (`/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng`) instead of the
> per-territory substrate (locus.fa + the local six-file syng sidecar
> set + the local routes), with the per-locus cut a coordinate-scoped
> operation against the global index — **conversion-not-materialization**,
> the same principle the re-anchor work proved (the 0.44s core
> transform; k-independent panel coordinates).

## 1. What is being deleted, and why it is safe

Today's per-territory chain (the runtime verification's measured shape):

| per-territory artifact | built by | wall (receipted) | fate under lever 1 |
|---|---|---:|---|
| locus.fa (the member-subset fasta) | extract (AGC getctg) | ~7s | **deleted from the cut/call path** — the AGC is fetched by coordinate instead |
| local syng (six-file sidecar set) | `impg syng` over locus.fa | ~3s + load per consumer | **deleted** — the global panel's own grid is the frame |
| node map (global→local node correspondence) | the cut's AGC-verified walk pairing | the cut's `load_seconds` pole (204-259s/territory measured) | **deleted** — no frame conversion exists because no second frame exists |
| local routes (532KB, derive-routes lookup) | `derive-routes --panel <local syng>` | ~4.8s | **deleted per territory** — ONE whole-panel routes inventory, built once, panel-side |
| window beds + locality maps (the window grid) | the panel-side substrate | panel-side, once per panel | **kept verbatim** — the global panel does not define the window partition |

The safety argument is the re-anchor principle: the syng node id is an
**opaque interned k-mer token**, and every semantic quantity the call
computes — the window grid, the member-row extents, the occurrence
intervals, the territory images, the spelled sequences — is
**coordinate- and content-carried**, not id-carried. The committed
per-territory path already proves the point in the opposite
direction: the cut's node map is a pure bijection over the touched
node sets, and the census/cache are frame-consistent within whichever
frame they were emitted. Removing the second frame removes the
conversion, nothing else.

## 2. The frame translation (the only new machinery)

The committed locality substrate is stated in the **territory-local
frame**: every member row's path name is a locus.fa record header
`path:start-end` and every coordinate is slice-relative (0-based
within the record). The record header IS the frame statement — the
same convention the cut's own provenance parse (`parse_locus_header`)
and the region-as-interface resolution
(`resolve_region_loci`: "the component's own path frame: its name ends
chr:start-end") already consume.

The global frame of the same row is exactly:

```
local:   ("HG002#1#chr1:4632437-4883593", start, end)
global:  ("HG002#1#chr1", 4632437 + start, 4632437 + end)
```

- **Forward (maps on load):** every partition-map member row and the
  component lane are translated at load by splitting the name at the
  final `:` and adding the interval offset. Fail-closed on malformed
  names (a map member whose name carries no `path:start-end` interval
  is a substrate error, loudly named).
- **Backward (product emission):** the called records are emitted in
  the partition's own reference frame — the committed
  REFERENCE-ANCHORED output shape — by the inverse map: a fold's
  global row `(path, g_start, g_end)` translates to the containing
  record `(path:s-e, g_start - s, g_end - s)`. Fail-closed if no
  record interval contains the fold row (a named structural error,
  receipted, never guessed).
- The window ids, partition ids, loci ids, and the census `partitions`
  field are frame-independent (the same axis file defines them in both
  frames).
- `--region chr:start-end` resolution is UNCHANGED: it consumes the
  component's own local frame today (the offset subtraction is
  already committed behavior) and the loci it returns are
  frame-independent.

The translation is implemented once, in
`src/genome_inference/territory_frame.rs`, and consumed by both sides
of the cut/call seam.

## 3. The interface change (the three-step UI stands)

**PARTITION** (`genome-infer partition`) gains `--global-panel
<prefix>` + `--global-routes <dir>`: when set, the local syng build
and the per-territory derive-routes are SKIPPED (prep ZERO's panel
side), the window beds are still regenerated verbatim from the
committed locality maps, and the registry (`partition.json`) records
the global pair plus `"frame": "global"`. The BED4 artifact is
unchanged.

**MAP** (`sample_projection cut`) gains `--frame global` (default
`local` = the committed path byte-for-byte): in the global frame the
cut takes NO `--panel` and NO `--sources` —
```
sample_projection cut --frame global \
  --global-panel <global syng> --agc <panel AGC> \
  --axis <axis json> --bed-directory <window beds> \
  --component <axis path> --shards-dir <projection> --out-dir <cut>
```
and does, in order: (1) load the global panel + AGC once; (2) translate
the universe rows to the global frame (§2); (3) build the territory
index against the GLOBAL panel over the translated rows — the per-row
slices for the reverse-frame extraction are fetched from the AGC by
coordinate (the same bytes the locus.fa extraction produced, by the
cut's own AGC-identity receipt: delta=1 uniform across 466-467
records, the committed offset convention); (4) derive the IDENTITY
EXPRESSIBILITY MAP — the canonical anchor walk of each record's own
content (fetched from the AGC at the record's header interval)
against the global panel's grid, every node mapping to itself; (5)
scan the shards through the committed remap flow over that map: a
record survives iff every anchor is one of the substrate's own
canonical anchors — exactly the strict local cut's evidence rule
(the node map's real role is the EXPRESSIBILITY FILTER, not the id
conversion; the read-edge qualification class is dropped by
construction either way), the survivors verbatim; (6) route every
touching pattern through the unchanged committed routing machinery
against the global-frame territory index. The emitted
census/cache/reads are in the GLOBAL frame: occurrence path ids are
global panel path ids, canonical tokens are global node ids, own
walks are the projection's own.

**The measurement that forced the expressibility rule** (the honest
course-correction, receipted): the first global-frame cut kept the
projection's records with NO expressibility filter and produced a
census 8.5x fatter than the committed one (reads_kept 891,230 vs
108,230; census_mult 976,657 vs 115,012 at L1; the fat dominated by
single-anchor high-multiplicity patterns concentrated in repeat
windows) — the read-edge qualification class: a short record's
syncmer qualification is computed on the READ's truncated windows,
so the committed slice-context rule drops almost all of them. The
call would have consumed a different evidence set, not a converted
one — the conversion principle requires the committed filter, restated
in the panel's own frame.

**CALL** (`genome-infer call-region`) gains `--global-panel <prefix>`
+ `--global-routes <dir>`: when set, every selected partition is
called with the global panel and the ONE whole-panel routes inventory
instead of the registry's local syng/routes, the locality maps are
translated on load (§2), and the product emission is translated back to
the partition's own frame. The underlying
`call-diplotypes` gains the same capability as `--global-frame` (the
engine switch), so the committed chain invocations and the
global-panel path run the SAME machinery.

**The whole-panel routes inventory** is built ONCE per panel, panel
side (the routes lever's own recipe, commit f3cbb90):
```
impg genome-infer derive-routes --panel <global syng> \
  --sources <panel AGC> --read-lengths 148 \
  --census <any global-frame cut census> --out-dir <global routes>
```
The call consumes from it ONLY the lane inventory (graph.json's lanes
+ source paths — the lookup lever's own survey). The census argument
is the routing-bind cross-check input; the artifact is a pure
panel-side function of the panel (the census-independence gate, f3cbb90).

## 4. What loads become shared, and the honest memory model

Per call process (the global-panel path):

| load | size | sharing |
|---|---:|---|
| `.1gbwt` (GBWT) | 9.5GB | **private per process** (C-struct, read by `syngBWTread`) — and **never queried by the call** (named follow-up below) |
| `.1khash` (k-mer dictionary) | 5.4GB | **private per process** (C-struct) — queried by the reverse-frame syncmer extraction |
| `.pstep` / `.spos` (path steps, checkpoints) | 8.9GB | **mmap'd, page-cache shared** — the call's walks (`walk_path_range`) touch only the called territories' member-path intervals |
| `.names` | 2.3MB | mmap'd, shared |
| the global routes graph.json | ~MBs | read once, shared page cache |
| the AGC (sequence source) | 3.3GB | shared page cache; the call fetches per-lane spans on demand (§5) |
| per-territory census/cache/reads/maps | ~1-3GB at 35x | per-territory working set, unchanged from today |

So the honest model is: **~15GB private per call process + ~12GB
shared page cache (the mmap'd sidecars + AGC, touched once for the
whole fleet)**, against today's per-lane ~0.1GB local panel + ~1-3GB
working set. The 20-minute budget's lane count is therefore
memory-bounded by the private sidecars: at the box's ~1TB, ~40-60
concurrent call lanes fit (measured in §8), versus the wave's 96-256
lanes of small local panels. Two named follow-ups sharpen it:
(a) **the gbwt-free call load** — the call never queries the GBWT
(grep-verified: the call path touches only `walk_path_range`,
`raw_matched_syncmers`, and the name map); a load variant that skips
the 9.5GB GBWT drops the private per-process fixed cost to ~5.4GB and
removes the dominant load wall; (b) **per-lane span loads of the
khash** are not possible without an index restructure (the dictionary
is a whole-structure hash), so 5.4GB private is the floor without a
new sidecar.

The fixed-load amortization: the global-panel call process calls K
territories per process (the cut-batch precedent), so the ~15GB fixed
load and the load wall amortize over K; across processes the mmap'd
sidecars and the AGC hit the page cache (the projection campaign's
measured shape: one shared read-only index, the pages faulted once).

## 5. The per-locus cut scoping, at the coordinate level

"Coordinate-scoped" means exactly: the territory index (the routing
structure the cut builds) is the global panel's own walks over the
translated row intervals — `walk_path_range(path, row.start - SLACK,
row.end + SLACK)` — plus the reverse-frame extraction over the AGC
slice fetched by coordinate. No sequence is materialized beyond the
row slices the committed cut already held in memory (locus.fa's
records vs the AGC's row slices: the same bytes, fetched on demand);
no second index is built; no node correspondence is derived. The
cut's per-territory `load_seconds` pole (the node map's AGC-verified
walk pairing, 204-259s measured at L1/T00000) disappears ENTIRELY —
there is nothing to convert.

## 6. The identity gate (the acceptance standard)

For the gate set (the localities with committed byte-identical call
receipts on disk: L1-L4 against `$PROJ/cut-L*-cli/calls.jsonl`), the
global-panel path is accepted iff the product `calls.jsonl` (and the
sequences emission) is BYTE-IDENTICAL to the committed per-territory
receipt. The run shape: fresh global-frame cut from the committed 35x
projection shards → global-panel call with the committed locality
maps, reads, and loci → `cmp` against the receipt.

The honest risk list — the named structural divergences that could
break byte-identity, each measured by the gate, never hidden:

1. **The edge-band anchor classes.** A locus.fa record's own syng
   computes syncmers with truncated context windows at the record's
   head/tail; the global panel computes the same k-mers with full
   chromosome context. Within ~w+k of a record edge the two grids can
   disagree: a record-edge syncmer that is not a global node (a DROPPED
   step globally), and a conserved-elsewhere k-mer that the global
   dictionary resolves but the local one lacks (an ADDED step
   globally). Fold node multisets (the QUAL spectrum's material
   signature) and the pin skeletons can therefore differ at the row
   edges; the likelihood's vote/extent machinery is coordinate-carried
   and unaffected. If the cluster machinery moves a QUAL at any gated
   locus, byte-identity fails LOUDLY at that locus and the delta is
   receipted (per-locus field diff), per the task's honesty rule.
2. **The census record-set margin.** The committed local cuts drop
   records carrying global-only anchors (the strict remap) or
   re-derive locally (the fleet's local-derive mode); the global-frame
   cut keeps the projection's records verbatim. The marginal reads are
   the edge-band class of (1) plus reads routed through territory
   slack. The §9 precedent (the three-step datasheet: a census that
   differed at the margin — reads_kept 108,230 → 94,242 — still
   produced byte-identical calls) says the called product is robust to
   census margins, but the precedent is not the proof; the gate is.
3. **The pair-bins convention.** The cut's census pair-bins compare
   contig suffixes; in the global frame the component's contig is the
   chromosome name (the local frame's is `chr:start-end`). The census
   pair_bins field is receipt-only (the call never reads it); the
   global cut derives the contig from the TRANSLATED component so the
   bins stay meaningful.
4. **The derive-routes routing bind.** The whole-panel routes
   inventory requires every census occurrence path id to resolve in
   the global panel's name map — a stronger bind than the local one
   (38,790 lanes vs the territory's hundreds). A global-frame census
   occurrence on a path outside the name map is impossible by
   construction (the routing resolved it against the global panel).
5. **The load wall.** The global syng load (~3 min measured in the
   projection campaign; the GBWT parse dominates) is a per-process
   fixed cost the local path does not pay. Without amortization (K
   territories per process) the global-panel call is SLOWER per
   territory; the lever's win is conditional on the batch shape, and
   the measurement (§8) states it as such.

## 7. What does NOT change

- The scoring machinery (`realign_score`) is untouched: same selection
  conventions, same fold/cluster/QUAL machinery, same receipts, same
  fail-closed binds. The engine switch (`--global-frame`) translates
  inputs at the boundaries (maps on load, rows on emission) and
  nothing between them.
- The committed CLI paths are immutable: every flag defaults to the
  committed behavior (`--frame local`, no `--global-panel` = the
  per-territory path byte-for-byte).
- The three-step UI shape, the BED4 artifact, the partition registry
  keys, the region-as-interface semantics, the truth discipline
  (test-mode only), no thresholds, no tuning constants (the SLACK and
  window margins are the committed machinery's own, restated not
  invented).

## 8. The measurement plan (duty 4; receipts land here as measured)

1. **Per-territory call wall, global vs local:** the same locality
   (L1) called both ways over the committed inputs; the phase table
   (load / binding / scoring) restated. Expected: the machinery wall
   ~equal (the same code, the same reads), the load wall moved from
   the per-territory syng build+load (small) to the per-process
   global load (large, amortizable).
2. **The lane-scaling measurement:** N concurrent call processes
   against the shared index (N=1..4 over the L1-L4 global cuts;
   co-tenant discipline stated), per-lane wall + peak RSS vs the solo
   run — the projection campaign's per-process-ceiling question asked
   of the call path.
3. **The budget restatement:** map ~3min + prep 0 + cuts (the
   global-frame cut-batch, amortized) + calls (the amortized shape)
   vs the 20-minute line, with what lever 2 (shared-index cut lanes)
   and lever 4 (the cross_lcp_merge marathon core) must still deliver.

## 9. Status

- [x] Design (this document).
- [x] Implementation: the frame module
  (`src/genome_inference/territory_frame.rs`), the `--global-frame`
  call engine switch (`call-diplotypes`), the `--global-panel`/
  `--global-routes` plumbing (`call-region`, `partition`), the
  `--frame global` cut (`sample_projection cut`), tests (the frame
  module's unit tests + the mini-panel conversion proof: the global
  territory index routes exactly what the committed local substrate
  routes, modulo the coordinate rename).
- [x] Identity gates: L1-L4 + the pilot, ALL RUN (§10.3 — the verdict
  is DELTA at every locality, never byte-identity: the divergence is
  the two named §6 risk classes, measured per-field and receipted;
  the elsewhere branch's record-context fix reproduces the committed
  E exactly at L1-L3, with L4/the pilot carrying the named
  panel-extra-record class; the §9 precedent (a 13%
  census margin with byte-identical calls) is NOT confirmed at the
  2.8x margin — the called fold material flips at 2-5 loci per
  locality. THE OWNER'S STANDING QUESTION the gates answered: the
  census-superset class is accuracy-relevant, not cosmetic).
- [x] Measurements + budget restatement (§10.3: the lane-scaling
  curve N=1/2/4/5 with the corrected memory model — peak Pss
  34-41GB private per lane, ~1.7-1.9GB shared, the ~24-29-lane ceiling —
  and the 20-minute budget restated from receipts only, with the
  named levers the line still needs).

## 10. The measured receipts (the state at this increment)

**The L1 global-frame cut** (the three-step regenerated window beds,
the committed 35x projection shards, 64 rayon threads, the box shared
with the live campaign):

| | committed local cut-L1 | global-frame cut (attempt 2) |
|---|---:|---:|
| wall | 606s | 533s (attempt 1, no filter: 674s) |
| load pole | 205s (local syng + node map + AGC verify) | 325s (the ONE global panel load; NO node map, NO locus.fa, NO AGC-identity walk pairing) |
| territory rows | 11,649 | 11,649 (the translation is exact on rows) |
| territory entries | 2,696,729 | 2,751,623 (+2.0% — the global grid's superset class at row edges) |
| expressibility anchors | 13,130 (node_map pairs) | 13,418 (the identity map; +2.2% — the record-edge class) |
| peak RSS | 38.0GB | 69.3GB (the +15GB panel; above the campaign's 64GiB per-process guard — the batch shape must carry it, or the gbwt-free load must land first) |

**The two course-corrections the measurements forced** (each receipted,
none hidden):

1. **The expressibility rule** (attempt 1 -> 2): keeping the
   projection's records with NO filter produced a census 8.5x fatter
   than the committed one (reads_kept 891,230 vs 108,230; census_mult
   976,657 vs 115,012) — dominated by single-anchor high-multiplicity
   patterns (the read-edge qualification class: a short record's
   syncmer qualification is computed on the READ's truncated windows).
   The committed node map's real role is the EXPRESSIBILITY FILTER;
   the global frame now derives the same set from the global panel's
   own grid (the identity map over the territory records' canonical
   anchors) and rides the committed remap flow.
2. **The AGC fetch convention** (attempt 2 -> 3): the frame-offset
   probe (examples/frame_offset_probe.rs, over the committed record
   HG002#1#chr1:4632437-4883593) proved the library fetch is
   0-BASED HALF-OPEN: the canonical anchor walk of fetch(s, s+len)
   aligns EXACTLY with the panel's own walk steps (82/82 anchors),
   while the shifted candidates are off-by-one (83/83). Attempt 2's
   +1-late reverse-arm steps corrupted the territory index (78% of
   CHM13 occurrence placements +1 vs the committed census; 1,848
   touching patterns dropped at routing). The committed local cut's
   "delta=1" receipt was the LOCUS.FA extraction convention (the
   substrate records are one base later than their headers), not the
   fetch convention — the global frame, reading the panel and the AGC
   directly, is freed of both.

**The identity gate: OPEN.** The attempt-3 census correspondence check
(the physical occurrence comparison, frame-independent) is the next
gate step; then the whole-panel routes inventory (derive-routes over
the global panel + AGC, once, panel-side), then the global-panel call
(recipe: run-gate-L1.sh in the work directory), then `cmp` against
$PROJ/cut-L1-cli/calls.jsonl. The named residual risks after the
convention fix: the expressibility superset (+2.2% anchors), the
territory-entry superset (+2.0%), and the record-edge anchor band —
each measurable at the census and the call.

**Attempt 3's receipts** (after the convention fix): routed 3,775 of
3,797 touching patterns (99.4% — the routing-drop class collapsed
from 1,848 to 22), wall 558.9s, reads_kept 304,793, census_mult
384,056. The physical-occurrence probe still shows a −1 delta vs the
committed local census — now UNDERSTOOD and it is NOT a bug: the
committed substrates' locus.fa content sits ONE BASE LATER than the
record headers (the delta=1 AGC-identity receipt, uniform across 933
records of both substrate classes: L1 466, T00000 467), so the
substrate frame's origin is s+1, while the global census sits at the
panel's TRUE positions. The committed call never needed the absolute
truth (a closed frame); the global-frame call DOES — the cut and the
call must translate by the SAME per-record content origin, and the
emission's back-translation cancels it exactly.

**The resume point** (the exact next steps, nothing else outstanding):
0. THE CONTENT-ORIGIN PARAMETER (the one implementation step that
   remains before the gate can run): the substrates' extraction
   convention is carried by their locus.fa, which the global frame
   deletes by design — so it must be stated as the cut/call's
   `--record-origin <d>` parameter (default 0 = the header's own
   interval; the committed substrates: d=1, per the delta=1 receipts —
   a convention restated from the receipts, not tuned). The cut
   applies it to the universe-row translation + the expressibility
   content fetch; the call applies it to the map translation (the
   frame table's origins = s + d); the emission's back-translation is
   origin-independent (it cancels). With d=1 the cut's territory, the
   call's folds, and the census occurrences all sit at the substrate's
   TRUE content, reproducing the committed grid exactly.
1. Re-run the L1 cut with --record-origin 1 (~9min) and re-run the
   occurrence-delta probe (the −1 class must collapse to exact
   matches against the committed census translated at s+1).
2. `impg genome-infer derive-routes --panel <global syng> --sources
   <panel AGC> --read-lengths 148 --census <the attempt-3 cut census>
   --out-dir <work>/global-routes` (~30min panel-side, once).
3. The L1 global-panel call via run-gate-L1.sh stage B (~5min incl.
   the panel load), then the byte-identity `cmp` (stage C).
4. On a pass: extend to L2-L4 (each needs a fresh partition step for
   its regenerated window beds, then the global cut ~9-11min each,
   then the call), then the lane-scaling measurement (N concurrent
   call processes against the shared index, N=1..4 over L1-L4), then
   the 20-minute budget restatement with levers 2 and 4's remaining
   share.
5. On a fail: receipt the per-locus delta honestly (the named
   structural candidates are all edge-band classes; the §9
   three-step precedent — a 13% census margin with byte-identical
   calls — is the prior, not the proof).

## 10.1 The content-origin increment's receipts (steps 0-1 of the
resume point: DONE; the exact state at the handoff)

**The `--record-origin <d>` parameter LANDED** (commit c899569):
the cut (`sample_projection cut --record-origin`, refused in the
local frame), the call (`call-diplotypes --record-origin`, requires
`--global-frame`), the three-step registry (partition records it,
call-region states it), the frame module (translate_row + the frame
table at `start + d`; the back-translation CANCELS d exactly — 4 new
unit tests; the mini-panel proof gains the ONE-BASE-LATER substrate
class). Tests 558/0. BOTH committed substrate conventions stand,
explicit and fail-closed: L1-L4 d=1 (delta=1 receipts), the pilot d=0
(its own delta=0 receipt — the smoke-cut-pilot's offset_variants).

**THE LEGACY-PANEL OPT-IN LANDED** (commit 56588b4, escalated and
supervisor-ratified): PanelIdentity's metadata-v2 gate (8f6a385)
refused the ONE committed global panel (metadata v1, May 13-14) at
BOTH derive-routes and call-diplotypes — the gate's rationale ("legacy
ambiguity-derived dictionaries are not acceptable for NEW samples")
does not cover this use: the panel is the committed evidence base's
OWN dictionary (the 35x shards' records carry ITS node ids; a fresh
rebuild renumbers every id and byte-identity becomes unattainable by
construction). `IMPG_PANEL_LEGACY_V1_OK=1` is the fail-closed-by-
default opt-in; every accepted use prints the loud provenance line;
derive-routes' manifest receipts panel_metadata_version +
panel_legacy_opt_in. THE RETIREMENT PATH (named follow-up): a
metadata re-stamp of the SAME dictionary (provenance upgrade, no
rebuild, no id churn) would retire the opt-in honestly.

**The d=1 cuts (the fresh window beds for L2-L4 via the fresh
partition steps — receipts in $WORK/partitions; L1 uses the committed
three-step beds)**:

| cut (frame global) | wall | routed | reads_kept | census_mult | committed local (routed / mult) |
|---|---:|---:|---:|---:|---|
| L1 --record-origin 1 | 563.4s | 3,775/3,797 | 304,793 | 384,056 | 3,514 / 115,012 |
| L2 --record-origin 1 | 542.5s | 3,542/3,546 | 141,920 | 156,976 | 3,412 / 43,051 |
| L3 --record-origin 1 | 504.2s | 3,660/3,709 | n/a | 567,243 | 3,450 / 17,887 |
| pilot --record-origin 0 | 485.0s | 12,286/12,310 | n/a | 513,836 | 12,218 / 512,323 |

(The global-frame census is a SUPERSET of the committed one — the
expressibility/record-edge class — reads_kept 2.8x at L1 but routed
patterns only +7%; the pilot, whose committed cut shares the same
smoke shards and exact-header convention, is nearly exact: +0.6%
routed, +0.3% mult. L4's cut is NOT yet run.)

**The occurrence-delta probe at L1 d=1** ($WORK/probe-L1-d1.txt):
the −1 off-by-one class is GONE (zero collapsible keys); 2,886 of
3,015 committed physical occurrence patterns match EXACTLY at the
translated content origin; the residual: 129 committed-only
patterns (mult 14,191) and 362 global-only (the record-edge superset
class) — the named structural risk #1/#2 of §6, now measured. The
census MASS divergence (committed 111,926 vs global 369,851 over the
common pattern set) is the real gate risk: the global census keeps
2.8x the reads (the superset expressibility class); whether the
calls are robust to it is exactly what the gate decides — the §9
precedent (13% margin, byte-identical) is the prior, not the proof.

**THE EXACT RESUME POINT** (nothing else in flight except the named
background jobs; all receipts in $WORK=/home/erikg/hprcv2/global-
panel-call-mode):
1. derive-routes over the d=1 L1 census → $WORK/global-routes WAS IN
   FLIGHT at the handoff (setsid, started ~15:47Z, ~30min expected,
   IMPG_PANEL_LEGACY_V1_OK=1; on completion verify graph.json +
   manifest panel_legacy_opt_in; if the process died, rerun
   $WORK/run-derive-routes.sh — idempotent, it rm -rf's first).
2. The L4 global cut: `$WORK/run-global-cut.sh L4 global-cut-L4-d1`
   (the L4 beds exist at $WORK/partitions/L4/beds; d=1).
3. THE GATES: `$WORK/run-global-call.sh L1 $WORK/global-cut-L1-d1 1
   $PROJ/cut-L1-cli/calls.jsonl` (likewise L2/L3/L4 vs cut-L*-cli);
   the PILOT: the committed-path baseline call failed once (the
   pilot's committed routes dir was reaped) — rerun run-smoke-cli-
   cut.sh with --routes /home/erikg/hprcv2/writeup-codification/
   substrate-pilot/routes (the v2 rebuild of the same substrate,
   built over the same smoke-cut-pilot census), then
   `$WORK/run-global-call.sh pilot $WORK/global-cut-pilot-d0 0
   <that-baseline>/calls.jsonl pilot-d0` AND a second cmp against
   $P/cli-product-hap1-alignfold/calls.jsonl (the on-disk committed
   receipt — expected to differ by the evidence-set margin, the L1
   precedent cut-L1.calls-diff.txt; receipt both honestly).
4. Lane scaling: N=1 (the L1 gate call is the solo baseline), then
   N=2/4 concurrent run-global-call.sh lanes (each RAYON threads=1,
   same box co-tenancy); the poller now logs Pss/RssAnon/shared per
   lane — the ~15GB private + page-cache model must hold or be
   corrected.
5. Restate the 20-min budget from the receipts; write §10 fully;
   commit per gate, explicit paths only; NO PR PUSH.
## 10.2 The gate-preparation increment's receipts (duties 1-3's
measured state; the gates' own receipts land in §10.3 as they run)

**The L4 d=1 global cut** (the last of the four committed-substrate
cuts; `run-global-cut.sh L4 global-cut-L4-d1`, the $WORK/partitions/L4
beds, exit 0): wall **554.2s**, routed **3,698/3,712** touching
patterns, reads_kept 240,407, census_mult 275,405, territory rows
11,610 / entries 2,931,539, load pole 324.6s, peak RSS **72.6GB** (the
same +15GB panel class as L1-L3; the cut path runs unguarded — the
batch shape carries it, per §4). The four-cut table (§10.1) now
reads: L1 563.4s / L2 542.5s / L3 504.2s / L4 554.2s / pilot d=0
485.0s.

**THE MEASURED SERIAL ROUTES WALL — the §10.1 "~30min panel-side,
once" estimate is WRONG by ~15x, and the correction is a lever in its
own right.** The whole-panel derive-routes (in flight, started
15:46Z Oct 17) passed its panel load and is in the per-lane source
hash: `panel_routes::lookup` computes each lane's `sequence_fnv1a64`
by reading the lane's FULL source bytes through the AGC
(core_bp=65536 crops), SERIALLY, one core at 100% — a linear pass
over the whole **1.403 Tbp / 38,790-lane** panel at a steady **~55
Mbp/s** (measured over the per-lane `source-ports/{id}.bin` creation
timestamps: ~500 lanes / 7.3 min, steady from lane 500 through lane
2200; ~6.1% of the panel bp at 16:25Z). Projected completion from
the 15:46Z start: **~23:30Z — a ~7.5h panel-side wall.** The run
STAYS ALIVE (setsid; supervisor-ratified): it is now the free
full-panel SERIAL byte-identity receipt generator — its graph.json
is the cmp target for the parallel artifact below.

**THE PILOT'S COMMITTED-PATH BASELINE: NO BASELINE EXISTS, BY PRIOR
RECEIPT — the §10.1 rerun recipe cannot produce one.** Rerunning
`run-smoke-cli-cut.sh` with the substitute routes
(`writeup-codification/substrate-pilot/routes`, the v2 rebuild built
over the same smoke-cut-pilot census; manifest verified) fails at
the CLI bind: `census occurrence start has no panel step` (exit 1,
wall 7s; receipts `smoke-cli-cut-v2routes.*` in hg002-projection).
This is NOT the reaped-routes issue the §10.1 resume point assumed:
the substitute routes loaded fine — the failure is the census-panel
bind, and it is the ALREADY-RECEIPTED fragmentary-record class
(d1cc84c: "the smoke-cut census is the autopsy's named
fragmentary-record bind rejection"; `scripts/reproduce-
writeup-receipts.sh` routes the pilot chain over the corrected cut3
census instead — and the cut3-based chain-pilot leg itself failed
fail-closed at census record 7479, the named 2-of-4 key ambiguity,
its manifest receipts status "failed"). So the smoke-shard cut
census is NOT a callable committed-path input in EITHER frame's
terms, and no committed-path baseline calls.jsonl exists to cmp
against. THE HONEST PILOT GATE FORM (supervisor-ratified): run the
global-frame pilot call over global-cut-pilot-d0; `cmp` against
`hg002-locus-pilot/cli-product-hap1-alignfold/calls.jsonl` — the
on-disk committed receipt, A DIFFERENT EVIDENCE CHAIN (the
`cosine-multi-census-hap1` derive-cache census, not a cut census;
named as such in the receipt) — with the delta receipted per-field
honestly as the evidence-set margin (the cut-L1.calls-diff.txt
precedent), and the committed-path bind failure carried verbatim.

**THE RATIFIED PARALLEL LANE-HASH (the serial wall's own lever;
supervisor-ratified before implementation, gated before any big
run).** The wall's entire content is `panel_routes::lookup`'s
per-lane source hash, so THAT loop is now parallel across lanes:
rayon `par_iter` over the lane vector with `collect` preserving lane
order; the per-lane crop loop is UNCHANGED; the fetch path is
thread-safe by construction (`Sources::fetch` is `&self`, the AGC
keeps per-rayon-thread decompressor pools, and core_bp crops exceed
the small-crop fetch memo — no shared mutable state anywhere in the
body). NO OUTPUT CHANGES, gated three ways: (a) the new unit test
`lookup_lane_hashes_are_deterministic_and_crop_independent` — the
lane rows are identical across core_bp crop partitionings (37 / 64 /
65536) and across reruns, in lane order (the hash is a pure function
of the lane's own bytes, the thread schedule invisible in the
artifact); (b) THE MINI-PANEL SERIAL-VS-PARALLEL CMP: the pilot panel
(458 lanes, the same inputs as the committed substitute artifact)
derived by a rebuilt HEAD-serial binary and by the parallel binary —
graph.json **BYTE-IDENTICAL, including compiler_identity** (the
parallelization contributes zero bytes); the Oct-14 committed
substitute artifact differs from BOTH only at compiler_identity
(byte 665, line 39) — the named code-drift class since Oct 14 (the
identity chain's own sources, vendor/syng's submodule state; the lane
rows — all 458 sequence_fnv1a64, families, source_paths, panel
identity, total_source_bp — identical); (c) cargo test --release
--lib **559/0** (+1), the sample_projection example proof **2/2**.
The serial whole-panel run STAYS ALIVE: its binary is the same HEAD
identity, so its graph.json (~23:30Z) is the free full-panel identity
receipt — expected BYTE-IDENTICAL to the parallel artifact byte for
byte. The parallel form serves lever 1's whole-panel routes inventory
AND lever 2's shared-index cut lanes (the same per-lane source-hash
shape at the cut's territory index build).

## 10.3 The gates' receipts (the byte-identity cmps, lane scaling,
and the 20-minute budget restatement — landing as they run)

**THE FIRST L1 GATE (the span-scoped LaneCache bug + the ELSEWHERE
FINDING — the headline).** The first complete global-panel call
(the parallel routes inventory, the L1 d=1 cut, run-global-call.sh L1;
wall 1489s solo: inputs 238.2s + binding 1047.7s single-threaded +
loci 162.1s; 1,475,421 AGC fetches / 358.9 MB; peak RSS 39.0GB)
needed one engine fix to complete at all and then RECEIPTED THE DIVERGENCE THE GATES EXIST TO DECIDE:

1. **The span-scoped LaneCache mis-slice (fixed, committed with this
increment):** d843da1's span-scoped fetch sliced the union-grown
span with ABSOLUTE path coordinates instead of span-relative ones —
the first non-zero-start span fetch panicked ('range start index
4764432 out of range for slice of length 426'); the covered-span
branch would have served wrong-offset bytes silently had the fresh
fetch not panicked first (the mini-panel proof never saw it: its
spans start at 0). The fix returns the span's OWN [lo,hi) bounds
with it and slices at start−lo; every served crop is a slice of the
same Sources::fetch answer, exactly as the policy's own contract
states.
2. **THE ELSEWHERE BRANCH IS PATH-LENGTH-CARRIED (the named structural
divergence, fixed by the supervisor ruling):** the committed local
panel IS the territory's member record set — its name map iterates
exactly the ~251kb `path:start-end` slices, so the committed E is
measured at the TERRITORY'S RECORD context (L1: E −14.418156, 466
strains, mean haploid 251,150.0 bp). The global panel's name map
iterates WHOLE CHROMOSOMES (mean haploid 3,011,203,926.4 bp over the
same 466 strain families) — a different model quantity: the first
gate receipted E −23.810253, a uniform −9.392 per-component additive
shift, QUAL-inflating at 12+ loci (mean Δqual +14,222, max
+161,986), and 8 called-fold-material flips. THE FIX (the ruling's
principled form, not a tuning constant): under --global-frame the
elsewhere branch is computed at the SAME RECORD CONTEXT as the
committed path — the territory's own record statements (the frame
table's `records()` enumeration, the frame-invariant extents
`end−start`, grouped by the same haplotype rule; verified pre-impl:
466 distinct records / 466 strains, mean 251,150.049 ≡ the committed
251,150.0; the per-window-crop alternative — the maps' 11,649 member
rows — measures 250,752.58 and is NOT the committed context). THE
RECEIPT: the fixed L1 call prints E = −14.418156, mean haploid
251,150.0 bp — the committed quantity reproduced EXACTLY; the
elsewhere_log_prob field is now IDENTICAL at every locus. THE
HONEST ALTERNATIVE (named for the owner, never switched silently):
the whole-chromosome context is arguably the MORE TRUTHFUL universe
— the 251kb slices were the materialization-era artifact — a
global-frame-NATIVE answer profile with the −9.392 shift named; a
new native profile is the OWNER's product decision, not ours.

**THE L1-FIX GATE (the bar: byte-identity vs cut-L1-cli; verdict:
DELTA, measured per-field, the remainder NAMED).** With the record
context in place: **3/25 loci FULLY IDENTICAL** (byte-equal call
records), 22/25 differ. The remaining divergence is exactly the two
named §6 risk classes, now measured:

- the census read-mass margin (§6 risk #2, the §10.1-known open
  question: the global census keeps 2.8x the reads — reads_kept
  304,793 vs 108,230): record_count/unit_count at 19/25 loci,
  best_log_likelihood at 22, alternative_log_gap at 20, qual at 20
  (marginal — mean Δqual dominated by 6 loci with |Δ|>0.5, max
  +33,630; the other 19 sit within ±0.21);
- the record-edge superset expressibility class (§6 risk #1: +2.2%
  identity anchors, +2.0% territory entries): class_count at 12/25,
  qual_alternative_fold_pairs at 9, called_classes/diplotype at 8,
  and the CALLED FOLD MATERIAL at **5/25 loci (2, 11, 13, 18, 22)**
  — the product-level divergence; 20/25 loci call the same fold
  material (id-identical or superset-renumbered).

Per the ruling: measured, named, and STOPPED — no further tuning;
the §9 precedent (a 13% census margin with byte-identical calls) is
NOT confirmed at L1's 2.8x margin: the census margin moves the
QUAL spectrum and flips 5 loci' fold material. The stale pre-fix
N=2 wave (L2-n2/L3-n2, walls 524s/451s, both DELTA) ran the
pre-fix binary and is superseded; its walls are warm-cache shapes,
not the scaling receipt.

**THE FULL GATE TABLE (all gates run against the parallel routes
inventory `global-routes-par`, the fixed binary, `cmp` against the
committed cut-L*-cli receipts; per-field receipts in
$WORK/gate-*.calls-diff.txt and the analysis tables below):**

| gate | loci | FULLY IDENTICAL | elsewhere E | fold-material flips |
|---|---:|---:|---|---:|
| L1 | 25 | **3** | EXACT (−14.418156) | 5 (2, 11, 13, 18, 22) |
| L2 | 25 | **4** | EXACT | 4 (5, 6, 16, 17) |
| L3 | 24 | **2** | EXACT | 5 (0, 1, 17, 19, 23) |
| L4 | 25 | **0** | −14.413791 vs committed −14.414212 | 2 (6, 23) |
| pilot | 22 | **0** (vs a DIFFERENT evidence chain) | −14.825332 vs −14.827120 | 1 (19) |

The gate verdict is **DELTA at every locality — never byte-identity**,
with the divergence exactly the two named §6 risk classes at every
locality (the census read-mass margin moving record_count/unit_count/
best_log_likelihood/alternative_log_gap/qual at 15-25 loci each, the
QUAL deltas marginal outside a handful of loci; the record-edge
superset expressibility class moving class_count/called_classes and
flipping the called FOLD MATERIAL at 2-5 loci per locality). The
§9 precedent (a 13% census margin with byte-identical calls) is NOT
confirmed at the 2.8x margin: L1-L4 flip fold material at 2-5 loci
each. THE DETERMINISM CROSS-CHECK: the N=2 and N=4 reruns of L2 and
L3 produced BYTE-IDENTICAL calls.jsonl (cmp clean) — the machinery
is thread-schedule-independent and the waves' walls measure real
concurrency, not noise.

**THE PANEL-EXTRA-RECORD CLASS (L4 + the pilot's small E mismatch,
measured and named).** L1-L3's frame tables reproduce their panels'
record sets exactly (E exact); L4's record-context E is −14.413791 vs
the committed −14.414212 (mean haploid 250,056.4 vs 250,161.7 — a
~49kb total difference over 466 strains) and the pilot's is −14.825332
vs −14.827120 (265,636.3 vs 266,111.6). THE PILOT'S CAUSE, measured
exactly (its panel is durable): the committed pilot panel carries
458 records, the called component's locality maps carry 457 — the
panel contains ONE record the maps do not (`NA20809#2#JBIREA010000087.
1:0-210118`, 210,118bp; 117,621,348 − 210,118 = 117,411,230 exactly
the maps' total), so the committed E was computed over a record the
global frame's only record statement (the maps) never states. L4's
committed panel was disk-cleaned (d1cc84c) so its ~49kb difference is
the same class by shape (the committed panels carried records outside
the called component's maps; the L1-L3 panels happened not to). This
is NOT fixable from the maps alone — the global frame's territory
record statement IS the maps — and per the ruling it is receipted,
not patched: restoring it would need a territory-record artifact
beyond the approved interface.

**THE PILOT GATE (the ratified honest form).** The committed-path
baseline does not exist and cannot be produced (§10.2: the
smoke-shard census fails the CLI bind in BOTH frames — the d1cc84c
fragmentary-record class, receipted verbatim; the cut3-based
chain-pilot leg itself failed fail-closed at record 7479). The
global-frame pilot call ran against `global-cut-pilot-d0` (exit 0,
wall 1680s: inputs 252.2s + bind 697.4s + loci 685.1s; 5,828,409
AGC fetches / 696.8 MB) and the cmp target is the on-disk committed
receipt `hg002-locus-pilot/cli-product-hap1-alignfold/calls.jsonl` —
A DIFFERENT EVIDENCE CHAIN (the `cosine-multi-census-hap1`
derive-cache census, not a cut census). THE DELTA, measured per-field
(0/22 loci fully identical): the evidence-chain margin moves
record_count/unit_count/best_log_likelihood/alternative_log_gap/
elsewhere_log_prob/qual at 22/22 loci (the committed run's E was over
its 458-record panel incl. the panel-extra-record class, §10.3);
the CALLED FOLD MATERIAL agrees at 21/22 loci (one flip: locus 19;
called_classes differ at 3) — the cut-L1.calls-diff precedent shape:
the called answers are robust to the evidence-chain margin while
every margin/QUAL field moves. Receipted, never tuned.

**THE LANE-SCALING RECEIPT (N concurrent call lanes against the ONE
shared index, RAYON_NUM_THREADS=1 per lane, the live campaign a
stated co-tenant throughout; per-lane wall + peak RSS/Pss/shared from
the poller logs):**

| lane | wave | wall | peak RSS | peak Pss (private) | shared |
|---|---|---:|---:|---:|---:|
| L1 (solo) | N=1 | 1445s | 39.1GB | 38.2GB | 1.7GB |
| L2 | N=2 | 525s | 35.5GB | 34.2GB | 1.7GB |
| L3 | N=2 | 452s | 38.0GB | 36.8GB | 1.7GB |
| L1 | N=4 | 1454s | 39.5GB | 38.7GB | 1.7GB |
| L2 | N=4 | 542s | 35.5GB | 34.0GB | 1.7GB |
| L3 | N=4 | 464s | 38.0GB | 36.3GB | 1.7GB |
| L4 | N=4 | 522s | 37.3GB | 35.8GB | 1.7GB |
| pilot (5th lane, largest census) | N=4+1 | 1680s | 42.1GB | 41.3GB | 1.9GB |

The scaling shape is NEARLY FLAT: L2/L3 at N=4 vs N=2 cost +3.2%/+2.7%
wall at double concurrency, and L1 at N=4 (1454s) vs its solo run
(1445s) costs **+0.6%** — with the N=1 and N=4 outputs
BYTE-IDENTICAL (cmp clean), the determinism receipt at both ends of
the curve (the bind is single-threaded per lane; the box absorbs 4
lanes trivially). THE MEMORY-MODEL CORRECTION (the
design doc's §4 must be corrected by this receipt): the "~15GB
private per call process" understates the measured per-lane private
footprint ~2.4x — **peak Pss 34-41GB private per lane** (the panel
sidecars' ~15GB PLUS the call's own working set; the committed local
call peaks at ~0.86GB for the same machinery, and the global census
is 2.8x the committed one — the superset class inflates the working
set: the pilot, the largest census of the set, peaks at 41.3GB Pss),
with only **~1.7-1.9GB shared** per lane (the touched pstep/spos
spans + the AGC — the page-cache-shared component of the model
holds, small). The lane ceiling at ~1TB therefore corrects from the
design doc's 40-60 to **~24-29 concurrent global-panel call lanes** —
memory-bounded by the private side, exactly as the doc's own
argument says, at the measured constant.

**THE 20-MINUTE BUDGET RESTATEMENT (from receipts only):**

- **map ~3 min** (the 12-lane projection basis, receipted in the
  closed runtime-verification campaign: 4.74M reads/s aggregate,
  byte-identical outputs) — unchanged.
- **prep ZERO per sample** (the global panel deletes the
  per-territory chain; §1's table). The panel-side ONE-TIME cost,
  now measured honestly: the whole-panel routes inventory is
  **1505s (~25 min) parallel** (the ratified parallel lane-hash, run
  once per panel, amortized over every sample and every territory;
  the serial form's ~7.5h wall is receipted in §10.2) — plus the
  window-bed regeneration (~12s/territory, the three-step partition
  receipt). THIS COST IS PER-PANEL, NOT PER-SAMPLE: at fleet scale
  it is ~0 per sample.
- **cuts**: measured 504-563s per sampled locality at 64 threads on
  the co-tenant box, of which **325s is the ONE global-panel load**
  (the amortizable pole — lever 2's shared-index cut lanes exist to
  share it): the marginal per-locality cut work is ~180-240s. For
  the whole genome (the live campaign's own census: 5,844 accepted
  territories), the per-territory marginal cut at the campaign's
  measured width-4 rate (459-712s per 4-territory batch) is the
  same minutes class; the whole-genome cut wall is then
  fleet-parallelism-bound, not machinery-bound.
  (The lane ceiling at the corrected ~24-29 lanes: a 5,844-territory
  genome at 25 lanes × ~9 min/lane ≈ 35 lane-hours; the 20-min line
  needs the per-lane wall down ~10x.)
- **calls**: the measured per-lane wall at N=4 is **~8.7-9.0 min**
  (522-542s) per sampled locality with the bind run
  single-threaded per lane (L1's bind pole: 1007s of its 1445s solo
  wall over 1.45M AGC fetches). The phase-3 bind seam is
  rayon-parallel in the committed code — the walls above pin it to
  1 thread per lane for lane comparability; releasing it is a
  per-lane thread-count decision, not new machinery. At the measured
  N=4 shape and the corrected ~24-29-lane memory ceiling, a
  5,844-territory genome at 25 lanes × ~9 min/lane is ~35 lane-hours
  — the 20-minute line is NOT reachable by lane parallelism alone at
  this per-lane wall: it needs the per-lane wall down ~10x (the
  parallel bind + the census-superset resolution both named) AND the
  marathon tail bounded (lever 4: the campaign's own deep-lane
  receipts — T03641 11,084s/71.1GB, T03203 2,884s/70.6GB — the
  cross_lcp_merge segdup core, still open).
- **The honest bottom line**: at the measured shape, the global-panel
  mode delivers map ~3min + prep ~0 + panel-side ~25min ONCE, and
  per-locality cuts/calls in the ~8.7-9.5 min/lane band at N=4 —
  the 20-minute whole-genome line remains BEYOND the current
  receipts, conditional on the named levers: lever 2 (shared-index
  cut lanes), the parallel bind seam, the census-superset question
  (the 2.8x margin — an accuracy-relevant divergence, §10.3's gate
  table), and lever 4 (the marathon tail). No number in this
  restatement is extrapolated beyond a measured receipt.

## 10.4 The exact resume point (if the budget ends before the last
receipt lands)

1. The serial whole-panel derive-routes (started 15:46Z Oct 17,
   setsid, ALIVE at 21:24Z at ~22.6k/38,790 lanes, ~55 Mbp/s serial,
   ETA ~02:00Z under co-tenancy) — on completion: `cmp` its
   `graph.json` against `$WORK/global-routes-par/graph.json` (the
   free full-panel byte-identity receipt; the same HEAD compiler
   identity — expected byte-identical; if it differs, name it
   per-field and re-examine the gates' verdicts against the serial
   artifact). THE DETACHED WATCHER (10.4 item 3) performs this
   unattended and writes the receipt.
2. The L1 N=4 lane RESOLVED (wall 1454s, DELTA, its calls BYTE-IDENTICAL
   to the solo run — §10.3's table updated) AND THE PILOT GATE LANE
   RESOLVED (exit 0, wall 1680s, DELTA vs the different-evidence-
   chain receipt, per-field in §10.3's pilot paragraph) — all five
   gates and all scaling waves are COMPLETE.
3. THE SERIAL FULL-PANEL CMP WATCHER (duty 5 of the 10x-path
   increment): a detached nohup loop
   (`routes-cmp-watcher.sh`, setsid, PID receipted in
   `routes-cmp-watcher.pid`) polls every 5 min for the serial run's
   `graph.json`; on completion it cmps the serial and parallel
   artifacts and writes the receipt to
   `serial-vs-par.routes-cmp-receipt.txt` in the work dir. Never
   polled in-session by design.
4. THE NEXT INCREMENT'S EXACT RESUME POINT (SUPERSEDED by §10.6 — kept
   for the record; the 10x-path increment-2 state): (a) the
   shared-index cut+call seam LANDED (§10.6.2, commit 2a88482, the
   combined gate PASS two-sided); (b) the span-cache delta-extend
   LANDED (§10.6.1, commit 44be7b0, the 508x fetch collapse,
   byte-identical); (c) the census-superset ANSWER question — the
   owner's, not ours (§10.5.4's shadow: ~6% wall at saturated
   territories, accuracy-dominated); (d) lever 4's marathon tail —
   OPENED (§10.6.3: the profile receipted, candidate A landed
   byte-identical, candidates B/C named; the marathon substrates are
   reaped, so the marathon gate form at the next campaign window is
   determinism + the wall receipt); (e) the exact NEXT poles are
   §10.6.4's list — the cross-process panel share (the private
   arenas), the call's own wall (the per-fetch AGC seeks, the loci
   core's B/C), the memory ceiling, and the marathon tail; the
   cheapest named step is the concurrent-in-process calls (the
   arithmetic at §10.6.4's close).

## 10.5 The 10x-path receipts (the bind-seam profile +
parallelization, the shared-index cut lanes, the census shadow)

### 10.5.1 THE BIND-SEAM PROFILE (measured before any code change)

**Method.** `perf(1)` is BLOCKED on this box (kernel
`perf_event_paranoid=4` — no unprivileged sampling; `samply` uses the
same `perf_event_open`), so the method is the committed phase
instrumentation — the `Instant` phase timers and atomic counters the
binary already prints to stderr (overhead: sub-1%, they wrap the same
code the committed receipts ride) — plus `/proc/PID/smaps`
per-mapping Pss snapshots every 5 min (no target overhead) and the
rss/pss poller. THE FRESH INSTRUMENTED RUN: the exact L1-fix baseline
recipe re-run serial (same binary, same inputs), exit 0, wall 1450s,
`calls.jsonl` BYTE-IDENTICAL to the committed baseline — the profile
run doubles as a free determinism re-gate.

**The serial-bind fraction of the per-lane wall** (the five lanes'
committed phase summaries + the fresh profile run):

| lane | wall | inputs | bind (serial) | loci | bind share |
|---|---:|---:|---:|---:|---:|
| L1 (baseline) | 1445s | 237.9s | 1006.6s | 159.6s | **70%** |
| L1 (profile re-run) | 1450s | 232.0s | 1013.8s | 163.5s | 70% |
| L2 (N=4) | 542s | 248.5s | 203.6s | 47.8s | 38% |
| L3 (N=4) | 464s | 245.3s | 25.1s | 150.6s | 5% |
| L4 (N=4) | 522s | 242.0s | 121.3s | 119.0s | 23% |
| pilot (N=4+1) | 1680s | 252.2s | 697.4s | 685.1s | 42% |

**What the bind is actually doing.** The bind is NOT the
panel-record pairing (the candidate lookup: `steps` 2.3-13.6s per
lane) and NOT the occurrence assembly compute (`shapes` ~0.0s) — it
is 99.3-99.7% the per-occurrence AGC sequence verification, and the
verification is ~100% **span-load time** (the `LaneCache`'s own nanos
counter ≈ the verify seconds at every lane):

| lane | records | occurrences | bind s | verify s | lane loads | load MB | load-time s |
|---|---:|---:|---:|---:|---:|---:|---:|
| L1 | 3,775 | 1,448,220 | 1013.8 | 1001.2 | 277,796 | 56,671.6 | 1016.3 |
| L2 | 3,542 | 1,387,508 | 203.6 | 201.2 | 47,704 | 10,713.2 | 209.0 |
| L3 | 3,660 | 1,405,637 | 25.1 | 22.0 | 7,211 | 995.3 | 29.0 |
| L4 | 3,698 | 1,415,097 | 121.3 | 118.6 | 30,291 | 6,638.0 | 125.6 |
| pilot | 12,286 | 5,793,360 | 697.4 | 683.4 | 159,540 | 38,468.3 | 685.6 |

THE REDUNDANCY DISEASE the counters name: mean load size 209KB ≈ the
territory extent (L1's 251kb) — every span-cache miss re-decompresses
the WHOLE grown span, and the union-growth policy misses 19.2% of the
1,448,348 range fetches (277,796 loads). L1 fetched **56.7GB for a
~117.6MB unique territory envelope set** (the 466 records' 117,621,348
bp; ~482x re-decompression), at the same ~3.0-4.3ms per load and
~46-58MB/s per-fetch AGC throughput the routes campaign receipted
(~55Mbp/s). THE LOAD VARIANCE IS THE ORDER, NOT THE MASS: L1-L4 all
carry ~1.39-1.45M occurrences, but their loads vary 38x (7,211 to
277,796) — the per-lane union-regrowth count is the occurrence VISIT
ORDER's geometry, not census mass. (The mean-load redundancy and its
fix — delta-extend the span instead of regrowing it whole — is the
named follow-on pole, receipted here, NOT implemented in this
increment: the mandated lever was the seam parallelization, and the
loads parallelize across lanes under it.)

**Where the 34-41GB working set comes from** (smaps per-mapping Pss
at the bind plateau, Pss 34.4GB): the working set is reached AT THE
INPUTS PHASE (rss inputs 36.5GB — the bind adds only ~0.3GB, rss
binding 36.9GB serial): the largest anon arena Pss **12.13GB** plus a
second **4.31GB** arena and the heap (1.26GB) and a tail of smaller
arenas — the panel's in-RAM index sidecars (the loader's own printed
material: the 9.5GB `.1gbwt` as 272,213,874 vertices / 683,539,246
edges; attribution by the loader's counts, not per-structure
measurement) — plus the census (3,775 records / 1.45M occurrences), the
routes inventory (graph.json 34.7MB, 38,790 lanes), the derive cache
(304,793 reads), and the reads. The FILE-BACKED touched sidecars
(pstep/spos class) are only ~0.8GB Pss — the page-cache-shared
component, matching §10.3's ~1.7-1.9GB shared receipt. The 34-41GB
private lane footprint is therefore the GLOBAL PANEL'S OWN IN-RAM
SIDELOADS, not call state — the structural fact the shared-index
lever (below) is built on.

### 10.5.2 THE BIND-SEAM PARALLELIZATION (the committed phase-3 seam,
released under the global frame) + the AGC per-thread decompressor
fix (the enabling change)

The seam machinery itself was already committed (phase 3:
`IMPG_REALIGN_PARALLEL_SEAMS=1 IMPG_REALIGN_SEAM_WIDTH=N` runs the
per-record binding loop and the per-locus loop on dedicated scoped
std::threads); the global-frame gates had simply never set it. THE
PROFILE'S BLOCKING FINDING FIRST: releasing the seam naively would
have been a NO-OP on the measured pole — `AgcIndex::fetch_sequence`
indexed its decompressor pools by `rayon::current_thread_index()` and
routed EVERY non-rayon thread (all seam threads) to the ONE fallback
slot behind a single Mutex (`unwrap_or(len-1)`), serializing every
seam-thread AGC fetch. THE FIX (src/agc_index.rs, committed with this
increment): every rayon thread keeps its own lazily-initialized slot
(exactly as before); every NON-rayon OS thread now claims its OWN pool
lazily, keyed by thread id — decompressors are independent readers of
immutable files, so which pool serves a fetch never changes the
bytes.

**THE GATE (byte-identity vs the serial bind on L1 — the bar: the
L1-fix baseline receipts on disk):**

| run | wall | bind | loci | peak Pss | calls.jsonl | sequences |
|---|---:|---:|---:|---:|---|---|
| serial baseline (L1-fix) | 1445s | 1006.6s | 159.6s | 38.2GB | — (the baseline) | — |
| serial profile re-run | 1450s | 1013.8s | 163.5s | 38.2GB | BYTE-IDENTICAL | — |
| seam W=16, rayon pinned to 1 | 551s | **110.4s (9.2x)** | 1750.5s | 47.7GB | **BYTE-IDENTICAL** | **BYTE-IDENTICAL** |
| seam W=16, rayon pool 16 | 419s | **118.1s (8.6x)** | 107.6s (1.5x) | 47.5GB | **BYTE-IDENTICAL** | **BYTE-IDENTICAL** |

Verdict: the bind seam is GATED at both rungs — `calls.jsonl` and the
sequences sidecar byte-identical to the serial bind; the machinery
is thread-schedule-independent at the seam too (the span cache's
union-regrowth LOAD COUNT varies with the schedule — 125,573 and
100,266 vs the serial 277,796 — but every served crop is a slice of
the same Sources::fetch answer, so the OUTPUTS are invariant; the
load counters are stderr diagnostics, and their sums differ from the
serial run's by construction, receipted here). The two honest
cautions: (1) the per-thread sub-timers now SUM across threads
(steps 363.9s + verify 1402.6s over a 110.4s phase) — timer sums are
no longer walls, only the phase seconds are; (2) THE LOCUS SEAM AT
RAYON-PINNED-1 IS A MEASURED 10.7x REGRESSION (16 seam threads each
enqueueing inner scoring par_iter jobs onto the same ONE-worker
global rayon pool: loci 1750.5s vs 163.5s serial) — the released-pool
rung (W=16, pool 16) is the receipted shape: loci 107.6s (the inner
par_iters get the 16-worker pool), bind 118.1s, wall **419s ≈ 7 min
vs 1450s ≈ 24 min — 3.5x per lane**. The memory cost of the released
seam: peak Pss 47.5GB vs 38.2GB serial (+9.3GB — the in-flight loads
and per-thread decompressors), which corrects the lane ceiling at
this rung to ~21 concurrent lanes at 1TB.

**The per-lane wall after the seam** (L1, the pole lane): inputs
244.7s / bind 118.1s / rows 23.1s / loci 107.6s / quality+elsewhere+
derive ~3s — THE NEW POLE IS THE INPUTS PHASE (the global panel's
own load, 244.7s of 419s = 58%), which is exactly the shared-index
lever's territory.

### 10.5.3 THE SHARED-INDEX CUT LANES (lever 2: the batch-amortized
loads) + the byte-identity gate

`sample_projection` gains `cut-batch-global`: ONE global syng + AGC
load serving a BATCH of global-frame territories (the projection
campaign's region-batching shape, at the global frame; manifest: one
JSONL row per territory {name, axis, bed_directory, component,
out_dir, record_origin}). Per territory the setup reproduces the
single global cut's own Global-frame branch verbatim (the translated
universe, the territory index over the shared AGC lane fetch, the
identity expressibility map); the ONE shard scan routes each read
through an inverted index (global node → territory ids — the local
fleet batch's own routing shape) and runs every routed territory's
committed identity remap; the per-(file, territory) states fold in
read order and merge in file order exactly as the single cut's do,
so each territory's outputs are the single cut's, byte for byte.

**THE GATE (the four-territory L1-L4 batch vs the committed single
cut receipts, all 12 files):**

| territory | cut-census.jsonl | cut-cache.bin | cut-reads.fasta |
|---|---|---|---|
| L1 | BYTE-IDENTICAL | BYTE-IDENTICAL | BYTE-IDENTICAL |
| L2 | BYTE-IDENTICAL | BYTE-IDENTICAL | BYTE-IDENTICAL |
| L3 | BYTE-IDENTICAL | BYTE-IDENTICAL | BYTE-IDENTICAL |
| L4 | BYTE-IDENTICAL | BYTE-IDENTICAL | BYTE-IDENTICAL |

Every semantic cut-stats.json field identical per territory; only the
timing fields differ (the amortized-load receipt below). **THE
MEASUREMENTS:** the 4-territory batch: wall 593s (613s process), of
which the ONE global load 344.7s (58%) — shared across all four —
and the scan 150.4s for ALL FOUR TOGETHER (vs 146.7-197.3s per SINGLE
territory: the shard decode dominates the scan and the per-territory
routing tests are cheap; the marginal scan cost is sub-linear in
batch width at this rung), setup + census emission ~100s. Against
the four single cuts' walls (563.4 / 542.5 / 504.2 / 554.2s =
2,164s total) the batch is **3.7x — the per-territory marginal cut
cost drops from ~541s to ~62s** once the load amortizes. THE MEMORY
MODEL AT BATCH WIDTH: peak RSS 76.3GB at width 4 vs 72.6GB for a
single cut — only +3.7GB for 3 extra territories (the shared ~40GB
panel dominates; per-territory indexes ~1-3GB), so the batch width
is memory-loose: at the same shape, width-16 batches stay in the
~110GB class, far inside 1TB — the cuts' fleet-parallelism ceiling is
process count, not memory.

**THE CUT+CALL SHARING QUESTION (measured facts, machinery NOT built
this increment):** the cut's loaded global SyngIndex and the call's
inputs-phase load are THE SAME PANEL (~36.5GB in RAM, 244.7s to load
at L1) — a combined cut+call process (or an injected-index call API)
would amortize the load across both phases; the memory model ALLOWS
it (a cut process already holds the panel at 76.3GB peak, and the
call's own working set beyond inputs is ~2-9GB), but the call's
`Options` accepts panel PATHS, not a loaded index — the injection
seam is the named next step, not silently built.

### 10.5.4 THE CENSUS-SUPERSET PERFORMANCE SHADOW (duty 4 — the
owner's accuracy decision UNTOUCHED; only its cost measured)

The question is the owner's (§10.3: the global census keeps 2.8x the
committed read mass — reads_kept 304,793 vs 108,230 at L1); this
increment measures only what that mass COSTS. THE INSTRUMENT
(receipt-side, default OFF): `IMPG_CUT_READSET_FASTA` restricts the
cut scan to the reads in a stated FASTA — here the committed
cut-L1/cut-reads.fasta (the committed 1x evidence mass, the natural
anchor). THE SHADOW CUT (exit 0, wall 523s): reads_kept **95,558**
(the committed set ∩ the global routing — 9,393 of the committed
distinct sequences sit in the global frame's inexpressible class and
do not route), routed_patterns 3,586 (vs the 2.8x run's 3,775; the
committed cut's 3,514), census_mult **141,114 vs 384,056 — the shadow
census is 0.37x the mass**. THE SHADOW CALL (the same global-frame
call machinery, serial, over the shadow census/cache; no answer gate —
the accuracy question is out of scope):

| phase | shadow (1x mass) | the 2.8x run (serial profile) |
|---|---:|---:|
| inputs | 247.1s | 232.0s |
| binding | 989.7s | 1013.8s |
| — lane loads | 278,591 / 57.0GB / 990.9s | 277,796 / 56.7GB / 1016.3s |
| — occurrences | 1,415,645 | 1,448,220 |
| loci (scoring) | **85.7s** | **163.5s** |
| wall | 1364s | 1450s |
| peak Pss | 37.3GB | 38.2GB |

**THE SHADOW, honestly: at a SATURATED-pattern territory the 2.8x
read-mass costs the call almost nothing.** The mass arrives as
MULTIPLICITY on the same pattern set (3,586 vs 3,775 records), and
the call's poles are both read-mass-invariant: the inputs (the panel
load — rss inputs identical at 36.5GB either way; the derive cache
95,558 vs 304,793 reads moves nothing measurable) and the bind (the
per-record occurrence structure: 1.42M vs 1.45M occurrences, the same
~278k span loads — the load work is order-geometry, not mass). THE
ONE REAL SHADOW: the loci/scoring phase, which scales with census
mult — **~78s per lane at L1 (85.7s vs 163.5s, 1.9x)** and ~1GB of
peak working set (37.3 vs 38.2GB) — a ~6% wall effect net. THE CUT
SIDE: the shadow cut's wall 523s vs the 2.8x cut's 563-586s (~8% —
the scan decodes the same shards either way). THE OWNER'S CAVEAT,
measured at the other end: where the pattern set is NOT saturated the
mass DOES bite — the pilot, the largest census (12,286 records /
5.79M occurrences), paid bind 697s with the largest working set of
the set (41.3GB Pss). The superset question is therefore ACCURACY-
dominated at saturated territories (§10.3's fold flips) with a small
performance shadow — the owner's decision needs the accuracy case,
not a performance case.

### 10.5.5 The 20-minute budget, restated from the new receipts

- **map ~3 min** — unchanged (the 12-lane projection receipt).
- **prep ZERO per sample**; the panel-side ONE-TIME routes inventory
  1505s parallel, once per panel — unchanged (§10.3).
- **cuts — THE BATCH RECEIPT (§10.5.3):** the 4-territory batch cut
  runs 593s of which the ONE 345s global load is shared by all four —
  the per-territory marginal cut is **~62s vs ~541s single** (8.7x);
  the scan's shard decode is shared per batch (150s for all four
  together, sub-linear in width at the measured rung). Peak RSS
  76.3GB at width 4 (+3.7GB for 3 extra territories — the shared
  ~40GB panel dominates): the batch width is memory-loose, the fleet
  ceiling is process count and scan CPU.
- **calls — THE SEAM RECEIPT (§10.5.2):** the per-lane wall at the
  receipted rung (seam W=16, pool 16) is **419s ≈ 7 min at L1** (vs
  1450s serial = 3.5x), gate-PASS byte-identical; peak Pss 47.5GB
  corrects the concurrent-lane ceiling to **~21 lanes at 1TB**.
- **The census-superset shadow (§10.5.4):** at saturated territories
  ~6% wall / ~1GB / the owner's accuracy decision, not a performance
  lever.
- **THE REMAINING POLES, each named with its receipt:** (1) the
  per-lane inputs phase — **244.7s = 58% of the 419s wall** — the
  global panel's own load, the same panel the cut process already
  holds (the combined cut+call process, or an injected-index call
  API: the memory model ALLOWS it, ~76-85GB at width 4; the seam is
  the named next machinery step); (2) the span cache's union-regrowth
  re-decompression — the ~482x receipt (56.7GB fetched for ~117.6MB
  unique at L1; a delta-extend span policy is the named follow-on
  pole, byte-identity by the same served-crop argument); (3) lever
  4's marathon tail — unchanged (T03641 11,084s/71.1GB,
  T03203 2,884s/70.6GB — the cross_lcp_merge segdup core, still open).
- **The honest bottom line:** the per-locality cost at the new
  receipts is ~62s (batch-marginal cut) + ~419s (call) per territory;
  at the ~21-lane ceiling a 5,844-territory genome is ~30 lane-hours
  of calls alone — the 20-minute line still needs the per-lane wall
  down ~10x, and the receipts now name exactly where it lives: the
  inputs pole (the shared-index cut+call seam), the load redundancy
  (the delta-extend fix), and the marathon tail. No number above is
  extrapolated beyond its receipt; the whole-genome cut/call wall at
  fleet scale is fleet-parallelism-bound (process count × scan CPU ×
  the 1TB memory model), not machinery-bound at the measured margins.

## 10.6 The increment-2 receipts (the two named fixes landed + the
marathon-tail profile opened + the budget restated; every gate
byte-identity, every number a receipt)

### 10.6.1 THE SPAN-CACHE DELTA-EXTEND (commit 44be7b0 — the ~482x
re-decompression, collapsed)

`LaneCache::span` on a miss now fetches ONLY the missing sub-ranges
([new_lo, lo) and/or [hi, new_hi)) and stitches them to the cached
interior — the covered interior is never re-decompressed. Byte-identity
by the served-crop argument restated for the stitch: a crop of a
Sources::fetch answer equals the direct fetch of the crop, so
fetch(new_lo, lo) is fetch(new_lo, new_hi)'s own head crop and
fetch(hi, new_hi) its tail crop, and the cached span is fetch(lo, hi)
by the maintained invariant — the concatenation IS fetch(new_lo,
new_hi) byte for byte, the buffer the union-regrowth form held.
LaneCache moved to module scope and the contract is unit-proved
(lane_cache_span_delta_extend_serves_direct_fetch_bytes: left/right/
both-side extension orders, every served crop byte-equal to the direct
fetch, the counters receipt the DELTA volume — 34,500B vs the regrowth
form's 92,000B on the same request order).

THE GATES (seam W=16, rayon pool 16, the receipted rung; cmp vs the
on-disk L1-fix / L2-fix-n4 receipts):

| lane | wall | bind | lane loads (count / bytes) | vs the serial form | calls + sequences |
|---|---:|---:|---|---:|---|
| L1 | 390s | 55.8s | 126,966 / 111.6 MB | 277,796 / 56.7 GB = **508x collapse** | BYTE-IDENTICAL |
| L2 | 345s | 50.9s | 21,288 / 111.4 MB | 47,704 / 10.7 GB = **96x** | BYTE-IDENTICAL |

The fetched volume is now ~1x the unique territory envelope (~117.6MB
at L1) — the re-decompression disease is gone; the load count is
schedule-visible (receipted as diagnostics only). Peak Pss 47.1 /
45.4GB (the seam rung's class, unchanged). THE NEW PER-LANE POLE the
gate names: the inputs phase (L1: 274.8s of 390s = 70%) — the global
panel's own load, exactly the shared-index seam's territory.

### 10.6.2 THE SHARED-INDEX CUT+CALL SEAM (commit 2a88482 — the
inputs pole deleted)

THE INJECTED-INDEX FORM: `realign_score::run` stays the committed path
verbatim; `run_injected(options, Option<&SyngIndex>)` runs the SAME
machinery over a caller-held panel. THE COMBINED PROCESS
(`sample_projection cut-call-batch-global`): manifest rows gain an
optional `call` object {reads, partition_graphs, loci, out_dir}; ONE
panel+AGC load, the cut phase = cut-batch-global's own core
(cut_batch_global_with, now taking the caller-held panel), then each
call row runs IN-PROCESS over the SAME loaded panel with the committed
CLI's own options construction (the instrument receipts + the product
emission + the sequences sidecar under with_output_model; the per-call
RSS guard defaults 0 = unguarded — the process holds the shared panel
by design).

THE COMBINED-PROCESS CALL BOUNDARY FIX the first gate forced (the
named divergence, receipted in the work dir's first-run-divergence.txt):
READ_LENGTH is a process-global stored from each call's derive cache,
but E is computed BEFORE that store — every committed receipt's E is
measured at the fresh-process 150, and in the combined process L2/L3/
L4's E read L1's stale 148 (a uniform +0.0063 E shift at every locus,
ZERO fold-material flips, verified numerically to 1e-9: E(148, a_L2,
501689.5) = −14.287349, exactly the observed delta, at the identical
mean haploid). `run_injected` now resets the global to the committed
fresh-process value (READ_LENGTH_FRESH = 150, the static's own initial
value) at every call entry — each call sees exactly the fresh-process
trajectory. The honest alternative (E at the cache's declared read
length) would move EVERY receipt's E — an answer-profile change, the
owner's alone, never switched silently.

THE GATE (two-sided, byte-identity; the 4-territory L1-L4 batch, one
process): the CUT side — ALL TWELVE FILES (cut-census/cache/reads at
L1-L4) byte-identical to the committed single-cut receipts; the CALL
side — ALL EIGHT FILES (calls.jsonl + the sequences sidecar at every
territory) byte-identical to the committed-path call receipts
(global-call-L1-fix, the fix-n4 receipts at L2-L4; E −14.418156 /
−14.293680 / −14.290571 / −14.413791, the receipts' own values).

THE MEASUREMENTS (the combined batch, 64 threads, the co-tenant box
stated): wall **1230s for 4 territories** = ONE **241.1s** load
(receipted in every territory's cut-stats.json) + cut phase 451s (the
153.7s scan shared by all four) + calls **147 / 112 / 126 / 125s**.
The per-territory call's inputs residual is **45-47s** (graph.json +
the AGC Sources open + census + maps) — the 274.8s panel-load pole is
GONE. Per-territory marginal cut+call at the combined form: **~240s**
(cut-share ~113s + call ~127s) vs the pre-increment ~481s (62s batch
cut + 419s call) and the campaign's original ~1,991s (541s single cut
+ 1450s serial call). Peak RSS **90.4GB / Pss 89.5GB at width 4**
(inside the 1TB model; the cut phase peaked 76.6GB, the call phase's
working set adds ~14GB at width 4).

### 10.6.3 LEVER 4 OPENED — THE MARATHON-TAIL PROFILE (the receipts
mined; the roster table in the work dir's marathon-tail-profile.txt)

Method: the marathon substrates were REAPED at campaign close (no
rerun is possible), so the profile is the campaign's own retained
per-lane instrumentation (prep/<T>/call.err per-locus [score] lines,
call.rss, call.wall) mined across the roster, plus the code anatomy.
THE ALGORITHMIC SHAPE: the marathon wall is the per-locus SCORING
MATRIX — the dense (folds × units) cross product walked per (fold,
unit) pair, O(folds × units × per-pair), with per-pair work whose
FOLD-INDEPENDENT parts were recomputed at every pair. Folds ≈ the
segdup copy count (the worst loci's member→fold ratios are nearly 1:1:
509→464, 985→884, 764→532) and units ≈ the routed read mass — so at
fixed per-copy depth BOTH factors grow with the copy count: **the
wall is quadratic in the segdup copy count**. Not a merge blowup (the
fold merge is seconds). THE MEMORY TAIL is the dense 13B/pair
ll/pinned/mismatches matrix (~33GB at the worst measured locus) — the
64GiB guard-kill class's direct cause (134 kills; completions 54-71GB).
THE PER-PAIR RATE BAND (the derived constant): pinned pairs 1.3-2.7
μs/pair; no-pin pairs 0.2-0.3 μs/pair — measured at the two ≥98%-no-pin
loci (T03803: 717.4s / 2.36e9 pairs; T03915: 434.8s / 2.07e9), and the
no-pin floor IS the fold-independent recompute class.

THE NAMED CANDIDATES: **(A) THE PER-RECORD HOIST — LANDED (commit
b80a687)**: the canonical decode, the own anchor positions, and the
merged anchor-window hull are pure functions of the RECORD; the
committed pair loop recomputed them per (fold, unit) pair (up to
~2.4e9 recomputes at the worst receipted locus). Now computed once per
record (7,665 decodes at T03641's worst locus, not 1.49e9); the pair
loop reads the table. GATE: L1 byte-identical (calls + sequences, the
committed receipt bar); the L1-scale wall effect is inside the
co-tenant band (the untouched bind phase itself moved 55.8→64.8s
between runs — the noise ruler) and the marathon effect is the
arithmetic, honestly a PROJECTION (no marathon rerun possible; the
marathon gate form at the next campaign window is determinism + the
wall receipt). **(B) THE ALWAYS-ON EXACTNESS ASSERTION** (direct_score
beside placement_score on every pinned pair — ~½ the pinned alignment
work in production; a receipt-shape decision, the owner's). **(C) THE
DENSE MATRIX** (row-compress the constant elsewhere rows + stream per
fold — the guard-kill class's fix; not small, named not attempted).
The irreducible core after A-C: anchor_correspondence per (record,
orientation, fold) + one alignment per pinned pair. THE HONEST
MARATHON BOUND stands: p90 2,912s to deepest 11,084s per lane
single-threaded, 15-25 members per genome — the tail alone can eat
the budget until B/C land.

### 10.6.4 THE 20-MINUTE BUDGET, RESTATED FROM THE NEW RECEIPTS

- **map ~3 min** — unchanged (the 12-lane projection receipt).
- **prep ZERO per sample**; the panel-side ONE-TIME routes inventory
  1505s parallel, once per panel — unchanged.
- **cuts + calls, THE COMBINED FORM (the receipted shape)**: ONE
  241.1s load per width-4 process shared by the batch's cuts AND
  calls; per-territory marginal cut+call **~240s** (cut-share ~113s +
  call 112-147s), all-in ~307s/territory at width 4. The call's own
  wall after both fixes: bind 50-73s (seam W=16), loci 19-137s,
  inputs residual 45-47s, rows 10-18s.
- **THE FLEET MATH AT THE CORRECTED CEILING (the honest arithmetic,
  no term extrapolated)**: a width-4 combined process peaks at 90.4GB
  → ~11 concurrent processes (~44 territories in flight) at the 1TB
  model → a 5,844-territory genome at the receipted shape is
  **~45-48 hours of wall** (1,461 processes × 1230s at 11-concurrent
  waves). The alternative receipted shape — standalone seam+delta-
  extend lanes at 47.5GB Pss, ~21 concurrent lanes, 390s/call — is
  **~30 hours** (the calls pole; the cut batches pipeline behind at
  18.5h). BOTH are far beyond the 20-minute line: 5,844 territories
  in 1,200s needs ~4.9 territories/s sustained; at the 240s
  marginal that is ~1,170 concurrent territory slots — **~26x the
  44-slot ceiling**.
- **THE REMAINING POLES, each named with its receipt**: (1) THE LOAD
  (241s per process — the last per-process pole; amortizes with
  batch width, but the structural fix is the panel shared ACROSS
  processes: the .1gbwt/.1khash are private in-RAM C-struct arenas,
  so cross-process sharing needs a shared-memory sidecar form — the
  named next memory lever); (2) THE CALL'S OWN WALL (112-147s: bind
  ~50-73s — the ~1.45M per-fetch AGC seeks are the floor after the
  span-volume collapse; loci ~19-137s — the marathon core's ordinary
  face, candidates B/C; the inputs residual ~46s — graph.json + AGC
  open + census, per-process cacheable); (3) THE MEMORY CEILING
  (90.4GB/width-4 process vs 47.5GB/single-call lane — the private
  arenas, not the working set); (4) THE MARATHON TAIL (§10.6.3's
  bound, unchanged until B/C land).
- **THE CHEAPEST NAMED NEXT STEP (arithmetic, not built)**: the
  combined process's calls run SEQUENTIALLY — the batch's K calls
  could run CONCURRENTLY inside the same process envelope (K call
  working sets ~10GB each + the shared ~40GB panel ≈ the same 90GB
  at K=4), and the machinery is thread-schedule-invariant
  (gate-proven at every rung so far), so outputs stay byte-identical;
  at K=4 that takes the width-4 process from 1230s to ~840s and the
  genome projection from ~45-48h to ~31h. It is a batch-shape
  decision beyond this increment's mandated seam — named, not built.
- **THE HONEST BOTTOM LINE**: the two named fixes landed and are
  gate-proven byte-identical — the per-territory cut+call marginal is
  down ~8.3x from the campaign's original ~1,991s to ~240s (the
  combined form), and the fetch-volume disease is gone (508x). The
  20-minute line remains BEYOND the current receipts: it needs the
  per-territory wall down ~10x MORE (240s → ~24s) and/or the memory
  ceiling up ~26x (the shared-panel lever), with the marathon tail
  bounded on top (B/C). No number above is extrapolated beyond its
  receipt; the walls carry the box's stated co-tenancy (±20%).

### 10.6.5 THE SERIAL FULL-PANEL CMP (duty 5, checked once at this
increment's close): the detached watcher (PID 2483558,
routes-cmp-watcher.sh) is ALIVE at 4h22m; the serial whole-panel
derive-routes is still IN FLIGHT (~495 min CPU, the ~55Mbp/s serial
lane-hash, ETA ~02:00Z) — the receipt file
(serial-vs-par.routes-cmp-receipt.txt) has NOT landed yet, exactly as
the watcher's design expects (it cmps the moment the serial graph.json
lands, unattended; never polled in-session). The receipt is owed to
the next increment's close.

### 10.6.6 THE EXACT RESUME POINT (if the budget ends here)

All three commits landed and gate-proven (44be7b0 the delta-extend,
2a88482 the combined cut+call seam + the READ_LENGTH call boundary,
b80a687 the per-record hoist); cargo test --release --lib 560/0; the
example proof 2/2; tree clean except vendor/syng (uncommitted by
design); NO PR PUSH. NEXT, in the §10.6.4 pole order:
(1) the serial-vs-par receipt at the watcher's completion (§10.6.5);
(2) the concurrent-in-process calls (the §10.6.4 arithmetic — the
cheapest named step, a batch-shape decision needing the owner's nod);
(3) the cross-process panel share (the private .1gbwt/.1khash arenas
— a shared-memory sidecar form, the memory-ceiling lever);
(4) lever 4's candidates B/C (the exactness assertion's receipt-shape
decision is the owner's; the dense-matrix row-compression/streaming is
the guard-kill class's fix); (5) the marathon re-receipt at the next
campaign window (determinism + wall, no committed receipts exist).

### 10.6.7 THE CONCURRENT-IN-PROCESS CALLS (increment 3, duty 1 — the
§10.6.4 cheapest-named step, LANDED and gate-PASSed)

**THE GLOBALS AUDIT (the boundary lesson applied to every
process-global the call path can reach, BEFORE any concurrency):**
the sweep found exactly ONE mutable process-global — the read
length (the §10.6.2 boundary bug's own static). Everything else the
call touches is per-call by construction (its own LaneCache, its own
Sources handles and per-OS-thread AGC decompressor pools, its own
receipt buffers and out dir) or already seam-thread-shared read-only
(the loaded panel — `unsafe impl Sync`, read-only C queries,
gate-proven schedule-invariant at 16 seam threads since §10.5.2; the
thread-local aligners; the env-derived OnceLocks are read-only after
init). **THE FIX: the read length is now PER-CALL STATE, not a
process-global** — a plain local in `run_injected`, initialized at the
fresh-process 150, stored once at the derive-cache load (the declared
length), shared immutably with the seam threads afterwards. The
committed trajectory is reproduced BY CONSTRUCTION at any call
concurrency: E reads the fresh value before the store (the constant
`read_len()`), every later site (the skipped-positions paths, the
unit likelihood, the scoring receipt) reads the declared value —
no reset discipline to forget, no stale-148 class possible. The free
functions take the length as a parameter (`unit_log_likelihood`,
`Scoring::floor`); the static is deleted.

**THE DRIVER:** the combined process's K calls now run CONCURRENTLY
— one scoped thread per call row (the manifest's own width, never
tuned) — inside the ONE loaded process envelope. Error discipline:
the first error in ROW ORDER fails the process (the serial form's
first-failure semantics preserved).

**THE GATE (the same two-sided bar, the committed receipts):** ALL
TWELVE cut files + ALL EIGHT call files (calls.jsonl + the sequences
sidecar at L1-L4) **BYTE-IDENTICAL** — the strongest determinism
receipt yet: four calls sharing one process, one panel, and one
rayon pool, byte-identical to the committed serial-per-territory
receipts (E exact at every territory, the per-call read-length
boundary included).

**THE MEASUREMENTS (the co-tenant box stated, 256 cores):**

| shape | wall | call phase | peak RSS / Pss |
|---|---:|---:|---|
| serial calls (the §10.6.2 receipt) | 1230s | 147/112/126/125s (510s total) | 90.4GB / 89.5GB |
| **4 concurrent calls (this)** | **923s** | **243/254/269/269s → 269s phase** | 118.0GB / 117.2GB |

The batch wall drops 1230s → 923s (**−25%**); the call PHASE drops
510s → 269s (**1.9x**). THE HONEST CONTENTION RECEIPT: each
individual call got SLOWER (243-269s vs 112-147s serial) — the bind
phase contends ~3x under 4-way concurrency (binding 162-186s per
call vs 50-73s serial; the seam's 4x16=64 threads and the per-call
AGC decompression compete for memory bandwidth on the co-tenant
box), so the concurrency wins the WALL by overlapping the
contended binds, not by making any call faster. THE MEMORY RECEIPT:
the peak grows 90.4 → 118.0GB (the four concurrent calls' working
sets — span caches, decompressors, receipt buffers — ~9GB each
beside the shared panel).

**THE FLEET-LEVEL HONESTY (the width-4 arithmetic, both terms
receipted):** the wall win and the memory cost CANCEL at the fleet
level at this width: serial 1,461 processes × 1230s at an 11-process
ceiling = **~45.4h**; concurrent 1,461 × 923s at an 8-process ceiling
(118GB each) = **~46.8h**. The concurrent calls are therefore NOT a
fleet-throughput lever at width 4 by themselves — they are a wall
lever that the memory ceiling currently absorbs; converting it into
fleet throughput is exactly the cross-process panel share's
territory (§10.6.9). The receipted value stands regardless: the
per-territory all-in cost at the combined form is now
923s/4 ≈ **231s**, and the width-4 process envelope is 118GB.

### 10.6.8 CANDIDATE C — THE DENSE-MATRIX ROW-COMPRESSION/STREAMING
(increment 3, duty 2 — the guard-kill fix, DESIGNED then IMPLEMENTED;
gate PASS at every receipt file; the bound receipted at the profiled
worst-locus shape)

**THE DESIGN (what the matrix's consumers actually read, the
exactness argument, and the memory bound it buys):**

*WHAT THE CONSUMERS READ (the full survey, every site):* the dense
(folds x units) ll/pinned/mismatches triple (13 B/pair, ~33GB at the
receipted worst locus T03411) has exactly seven consumer classes:
(1) THE EXHAUSTIVE CLASS SCAN — for each of the folds^2/2 class
pairs, every unit in order: `count * mix_logsumexp(ll_a[u], ll_b[u])`
— the accumulation ORDER is the byte-identity invariant (a sparse
iteration that skips the constant entries would change the float
summation order and BREAK byte-identity; the compressed scan
therefore walks the same units in the same order reading the same
values); (2) THE MERGE PASS — the ll columns' bit-exact equality
(hash-bucketed, exactly verified) and the representative-row clones
— the predicate is over the LL VALUES ALONE (pinned/mismatches ride
the representative); (3) the exactness sample and the named-fold
strata — the pinned pairs enumerated in (fold, unit) order;
(4) the anatomy — winner/truth-pair pinned probes; (5) the records
sidecar — per-unit fold_ll for the named pairs; (6) THE INGREDIENTS
SIDECAR — the ONLY full-matrix reader: serde_json of the whole
ll_matrix, ~20 B/pair of JSON text buffered IN RAM per locus (a
second, larger matrix beside the first — at T03411's shape ~50GB,
the co-pole of the same kill); (7) the scoring receipt's rss probe.

*THE EXACTNESS ARGUMENT (bit-for-bit, NO approximation):* the pair
loop's two paths PARTITION the value space — a pinned pair's local
branch is finite, so its marginalization logsumexp2(local, E) is
STRICTLY greater than E; an unpinned pair's mismatches stay 0 and
its ll is EXACTLY E (the empty local branch is -inf). So an entry
equals the no-pin constant (E, false, 0) IFF the unit is unpinned
on that fold: compressing the constant entries away is a
VALUE-IDENTICAL encoding, not an approximation. Every consumer reads
the same bits, booleans, u32s, in the same iteration order, so
every downstream float operation is unchanged — byte-identity by
construction, proven again by the gate.

*THE ROW ENCODING (a DERIVED per-row choice, never tuned):* each
fold row takes the SMALLER of the two forms — the sparse form costs
16 B/pinned entry (unit + ll + mismatches, structure-of-arrays),
the dense 13 B/unit — so the pinned-majority marathon loci
(T03641 at 84% pinned, T03499 at 86%) keep the committed size
EXACTLY (never larger), and the no-pin-majority loci (T03803 at 2.1%
pinned, T03915 at 1.3%) collapse. The bound is per row:
min(13 x n_units, 16 x pinned_in_row) bytes — the whole-matrix
bound rides the roster's own pinned shares. THE CLASS SCAN reads
through a monotone cursor (O(1) amortized per unit, sparse or dense
— no binary search on the scan path, and FEWER bytes streamed at the
no-pin-majority loci than the dense form).

*THE STREAMING POLE:* the ingredients sidecar's ll_matrix now
serializes ROW BY ROW from the compressed form — serde_json's own
per-value bytes in the committed separator layout, streamed to a
per-locus spill file and copied into the ordered ingredients file at
the write phase (the seam threads never share the file; the final
bytes are the committed form's exactly — the gate cmps the file).
The dense JSON text never materializes in RAM.

*THE MEMORY BOUND IT BUYS (the roster's own numbers):* the matrix at
T03803's shape (532 x 4,438,575 = 2.36e9 pairs at the roster's own
2.1% pinned share): **750MB live (0.33 B/pair) vs the dense form's
29.3GB (13 B/pair) = a 39x collapse**; at T03641's pinned-majority
shape the bound is the committed 13 B/pair exactly (the derived
choice keeps dense); the ingredients JSON (the ~50GB co-pole at
T03411's shape) is streamed at EVERY locus. The guard-kill class's
arithmetic: the worst roster case's peak becomes the pinned-majority
matrix (~19-24GB) + the streamed receipts + the lane base (~6GB) —
under the 64GiB guard at every roster locus, where the committed
form's matrix + JSON buffer exceeded it.

*WHAT IS HONESTLY NOT CLAIMED:* the marathon WALL effect (the class
scan's fewer bytes, the scoring loop's unchanged work) is a
PROJECTION — the marathon substrates are reaped, no rerun is
possible; the marathon re-receipt at the next campaign window is
determinism + the wall + this bound. The 2.5e9-pinned-pairs class
keeps its ~20-24GB matrix — the exact-value floor of THIS machinery;
removing it would need the class scan's algebraic decomposition
(mix_logsumexp terms precomputed per fold), which changes the float
accumulation and BREAKS byte-identity — named and refused, not
attempted. (SUPERSEDED BY THE OWNER'S RULING 3, §10.6.15: the
sparse-skip form now exists as the NAMED fast mode behind
IMPG_CLASS_SCAN_ORDER=sparse-skip — never the default, its float
deltas measured and labeled in every output; the exact unit-order
accumulation remains the receipts' bit-reproducible form.)

**THE IMPLEMENTATION:** `ScoreMatrix`/`ScoreRow` (module scope, with
the partition and encoding arguments as definition comments), the
pair loop building rows via `ScoreRow::from_pair_loop`, the merge
over the compressed rows (the canonical value hash buckets, the
entrywise bit-exact verify), the class scan via
`ScoreMatrix::class_pair_sum` (the committed accumulation verbatim),
the consumers on probes, the ingredients streaming with spill
cleanup on the error path, and `live_bytes` (the per-locus stderr
receipt: every locus prints its pairs, live MB, and the dense-form
MB). The dense clone-collect (the committed form's SECOND full copy
of the matrix, transiently two) is deleted with the rework.

**THE GATES:**
- *The unit proof* (`score_matrix_row_compression_is_value_identical`):
every probe equals the dense form's own bits, the pinned
enumeration is the dense filter's order, the class pair sum is
bit-for-bit the committed dense loop at every fold pair, the merge
predicates agree with the dense bit comparisons ACROSS ENCODINGS (a
dense row and a sparse row with the same values hash and compare
equal), the streamed JSON is serde_json's own bytes for the dense
form, and the encoding choice matches the derived crossover — for
no-pin-heavy, pinned-heavy, and all-no-pin shapes.
- *The L1 gate (the non-marathon territory, the committed receipts):*
the combined-batch run with candidate C — **every receipt file
byte-identical vs BOTH baselines** (the committed
global-call-L{1..4}-fix receipts AND the duty-1 concurrent run):
calls.jsonl + sequences + **the ingredients file (the streamed
ll_matrix IS the committed buffered bytes)** + records + skeleton +
exactness; the instrument receipt field-identical mod the walls and
rss probes; the cut side 12/12 byte-identical. Wall 891s, call phase
260s, peak RSS 112.6GB (the co-tenant band around duty 1's 923s /
118GB).
- *The marathon-class receipt (no committed marathon receipts exist
— the honest gate is determinism + the bound):*
`marathon_worst_locus_shape_memory_receipt` builds the profiled
worst-locus shape (T03803: 532 x 4,438,575, the roster's own pinned
share) and measures: live 750MB / dense 29.3GB = **39.0x**,
measured RSS delta 754MB (the live bytes' own class, allocator
overhead included); the pinned-majority control asserts the derived
choice keeps the committed dense size exactly. Determinism: the L1
gate's byte-identity at both baselines + the unit proof's
bit-equality.

### 10.6.9 THE CROSS-PROCESS PANEL SHARE (increment 3, duty 3 — the
slot-ceiling lever, OPENED: the design + the measurements; the
implementation is its own increment with full gates, NOT attempted
this window)

**WHAT IS LOAD-TIME-DERIVED vs INTRINSICALLY-PRIVATE (the receipts):**

| sidecar | on-disk | how it loads | sharing today |
|---|---:|---|---|
| `.pstep` | 3.24GB | memmap2 | **already shareable** — page-cache; only touched pages materialize (~0.8GB Pss at the bind plateau) |
| `.spos` | 5.09GB | memmap2 | **already shareable** — same |
| `.names` | 2.3MB | mmap | shared |
| `.1gbwt` | 8.85GB | the C reader (`syngBWTread`) deserializes into a **malloc'd anon arena** (the loader's own print: 272,213,874 vertices / 683,539,246 edges — the 12.13GB anon arena of the §10.5.1 anatomy) | **intrinsically private TODAY** — not because the data is private but because the C reader materializes it in process memory |
| `.1khash` | 5.07GB | `kmerHashReadOneFile` (ONEcode) into the **4.31GB anon arena** + heap | **intrinsically private TODAY** — same |
| the AGC | 3.3GB | page cache via per-call `Sources` | shared |

THE LIVE RECEIPT (the serial derive-routes process, PID 2121938,
holding the whole panel at 43.3GB RSS): **anon private 39.72GB;
file-backed touched 1.60GB** — the file-backed class is the shareable
one and it is TINY; the mass is the C-arena class. The bytes of both
arenas are a PURE FUNCTION of the sidecar files — the privacy is an
artifact of the loader, not of the data.

**THE THREE SHAREABLE FORMS, and the ceiling each buys (at the ~1TB
model, arithmetic from the receipts):**

(a) **FORK-AFTER-LOAD COW (the natural Linux form):** the parent
loads the panel once (the receipted 241s load), children fork per
batch/territory; the arenas' pages stay shared while unwritten
(COW). The shared base is the parent's ~45GB (panel + AGC); each
child's private is its own per-territory work — the receipted
per-territory call working set (~9GB: the width-4 concurrent peaks
112.6-118GB over the 76.6GB cut-phase base, /4), the cut territory
index (1-3GB), and the per-territory inputs (census/cache/reads —
per-territory data, private under every form). **Ceiling: ~
(1,007 − 45) / 13 ≈ 74 territory slots ≈ 2.3-2.5x the current
32-slot ceiling — AND the 241s load pole amortizes to ~0 per
territory** (the parent's lifetime is the fleet run). THE NAMED RISK
the implementation gate must measure FIRST: any post-load WRITE
into the arenas (lazy C-side initialization, in-place counters)
breaks COW page by page — a forked child's Pss growth over a full
call is the gate; the seam-thread era already relies on the same
read-only-post-load claim for thread safety, but fork makes it a
memory claim. The prefork hygiene is standard: rayon pools and
thread-locals re-initialize per child (post-fork single-threaded
start); the C globals are the §10.6.7 audit's own read-only class.

(b) **SHARED-MMAP ARENAS (the sidecar restructure):** convert the
two C arenas to a fixed-layout mmap-able form (a new sidecar class,
or an attach entry point that consumes MAP_SHARED memory), so every
process maps the same pages — the same ceiling as (a) without the
parent lifecycle. The cost: a C-side attach API (`kmerHashReadOneFile`
populates malloc'd memory today) or a Rust-side reader
re-implementing the ONEcode parse into an mmap-backed struct — the
khash is a whole-structure hash (the §4 named floor: no per-lane
span loads without an index restructure). Same ~74-slot ceiling
class.

(c) **A PANEL-SERVER PROCESS — REJECTED by the receipts, honestly:**
the bind's measured pole is the AGC span loads (per-call LOCAL
`Sources` — a server does not move them), and the panel queries are
C calls at the sub-microsecond class; an IPC round trip (1-10us,
even shared-memory) would INVERT the bind's cost. Only viable if the
query surface is coarsened — a different machinery, named as the
fallback if (a) and (b) both fail their gates.

**THE HONEST CEILING STATEMENT:** the panel share is a ~2.3-2.5x
SLOT lever plus the deletion of the per-process load pole — it does
NOT touch the per-territory wall (the bind's AGC-seek floor, the
loci core). The 20-minute line needs BOTH this lever AND the
per-territory wall down (§10.6.10).

### 10.6.10 THE 20-MINUTE BUDGET, RESTATED AT INCREMENT 3 (duty 4 —
every number a receipt; the fleet math at the receipted shape)

- **map ~3 min** — unchanged (the 12-lane projection receipt).
- **prep ZERO per sample**; the panel-side ONE-TIME routes inventory
  1505s parallel (the serial form still IN FLIGHT at this
  increment's close — §10.6.11) — unchanged.
- **cuts + calls, the combined+concurrent+candidate-C form (all
  gate-PASSed this increment):** the width-4 process runs **891s**
  (the candidate-C gate run; the duty-1 concurrent run 923s — the
  co-tenant band) = ONE **~220s load** + cut phase **408-436s** +
  call phase **260-269s for ALL FOUR calls concurrent** (the
  per-call walls 201-269s: the binds contend ~3x under 4-way
  concurrency — §10.6.7's honest contention receipt). Per-territory
  all-in **~223s**; marginal ex-load **~167s**. Peak **112.6-118.0GB
  / Pss 111.4-117.2GB** at width 4.
- **THE FLEET MATH AT THE RECEIPTED CEILING:** 112.6-118GB per
  width-4 process → **8 concurrent processes (32 territory slots)**
  at the 1TB model → a 5,844-territory genome is 1,461 processes x
  891-923s at 8-concurrent waves = **~45-47 hours** — the concurrent
  calls' wall win (1230 -> 891-923s) is exactly absorbed by the
  memory ceiling it costs (90.4 -> 112.6-118GB, 11 -> 8 processes):
  at width 4 the concurrency is a wash at the fleet level (§10.6.7's
  arithmetic, restated here as the budget's own term).
- **WITH THE PANEL-SHARE LEVER (the §10.6.9 design, UNGATED — a
  projection, named as such):** the fork-COW form's ~74 slots at
  ~290s/territory (cut-share 113s + call 130-150s + inputs 46s, the
  receipted walls) is **~6-7 hours** genome wall IF the COW sharing
  holds (the fork-probe gate named in §10.6.9) — the single biggest
  structural lever now open, and the load pole gone.
- **THE REMAINING POLES, each named with its receipt:** (1) THE LOAD
  (220-241s per process — amortizes with width; the structural fix
  is the cross-process share, §10.6.9); (2) THE CALL'S OWN WALL —
  bind 50-73s solo / 137-186s under 4-way contention (the ~1.4-1.5M
  per-fetch AGC seeks are the floor after the §10.6.1 volume
  collapse; the LaneCache's load-time counter receipted 773.7s of
  summed thread-time under concurrency — the contention is
  decompression bandwidth, not volume), loci 16-67s at L1-class
  (the marathon core's ordinary face; candidate B and the
  census-superset are the owner's decisions), inputs residual ~46s;
  (3) THE MEMORY CEILING (112.6-118GB/width-4 — the private C
  arenas + the concurrent working sets; §10.6.9 is the lever);
  (4) THE MARATHON TAIL — the WALL bound unchanged (p90 2,912s to
  11,084s single-threaded, 15-25 members/genome), but the KILL class
  collapses by candidate C's arithmetic at every roster locus (the
  worst roster peak ~24GB matrix + streamed receipts + ~6GB base vs
  the committed 33GB + ~50GB JSON over the 64GiB guard) — the
  marathon re-receipt at the next campaign window is determinism +
  wall + the bound.
- **THE HONEST BOTTOM LINE:** the increment-3 receipts put the
  per-territory all-in cost at **~223s and the fleet at ~45-47h** at
  the receipted 32-slot ceiling. The 20-minute line (5,844
  territories in 1,200s = 4.87 territories/s) needs **~44x the
  current throughput**: the panel share (a ~2.3-2.5x SLOT lever,
  designed ungated) AND the per-territory wall down ~an order
  (the bind's seek floor, the loci core, candidate B — the owner's,
  the census-superset — the owner's) AND the marathon wall bounded —
  or slot counts beyond every named form (the line at 290s/territory
  needs ~1,415 slots). No number above is extrapolated beyond its
  receipt except the two named projections (the fork-COW sharing,
  the marathon wall effect), each labeled and gated.

### 10.6.11 THE SERIAL FULL-PANEL CMP — LANDED AT THIS INCREMENT'S
CLOSE (duty 5, checked once; the watcher's unattended cmp)

The serial whole-panel derive-routes (started 15:46Z Oct 17) COMPLETED
at 01:42Z Oct 18 (exit 0, wall **35,778s ≈ 9.9h** at the ~55 Mbp/s
serial lane-hash under the campaign's co-tenancy — the §10.2-measured
~7.5h class, honestly restated at its actual wall), and the detached
watcher's unattended cmp wrote the receipt the moment it landed:

- **graph.json: BYTE-IDENTICAL** — the free full-panel byte-identity
  receipt **PASSES**: the parallel lane-hash contributes ZERO bytes
  over the whole panel (38,790 lanes / 1.403 Tbp; serial 35,778s vs
  parallel 1,505s — the ~23.8x panel-side ratio, one-time per panel).
- registry.jsonl / ports.bin / ownership-input.json / ownership.jsonl
  identical; source-ports 38,790 = 38,790; native dirs equal.
- manifest.json DIFFERS — the run-metadata class (the manifest's own
  timing/argv fields, the §10.2 compiler-identity precedent's
  neighbor; the artifact's semantic members are the identical ones).

**THE EXACT RESUME POINT (if the budget ends here):** all four duties
landed and gate-proven (3e9f9b8 the concurrent calls + the per-call
read length; c4b7e86 candidate C's row-compression + streaming;
603a226 the panel-share opening + the budget restatement; the serial
cmp receipt above); cargo test --release --lib 561/0 (+2: the
compression value-identity proof, the marathon-shape memory receipt —
the latter #[ignore]d, run explicitly: 39.0x at T03803's shape); the
example proof 2/2; tree clean except vendor/syng (uncommitted by
design); NO PR PUSH. NEXT, in the §10.6.10 pole order:
(1) THE PANEL SHARE'S IMPLEMENTATION INCREMENT (the fork-COW probe
FIRST — a forked child's Pss growth over a full call, the §10.6.9
named risk gate — then the prefork batch driver, full byte-identity
gates at every step); (2) candidate B (the always-on exactness
assertion — the owner's receipt-shape decision) and the
census-superset answer (the owner's); (3) the marathon re-receipt at
the next campaign window (determinism + wall + candidate C's bound);
(4) the bind's per-fetch AGC-seek floor (the span-volume collapse
went 508x; the seek COUNT is the named floor — a bulk-prefetch form
is the candidate, undesigned).

### 10.6.12 THE FORK-AFTER-LOAD COW PANEL SHARE — LANDED (increment
4; commit 3e161eb the driver + the COW audit + the gates
below as run, this commit the receipts and the corrected budget): the §10.6.9 form (a), implemented and gated; the COW
holds at every measured concurrency; the ceiling CORRECTED from the
receipts; the budget restated

**THE COW WRITE AUDIT (duty 1, both halves):**

*The code half — every post-load writer into the panel's arenas,
named:* (1) the GBWT query paths are PURE READERS —
`syngBWTmatchStart/Next`, `syngBWTpathStartOld/pathNext`,
`syngBWTincomingRank`, `syngBWTadvanceRank`, `syngBWTlocFind`, and
the `rsDirRankSyng`/`rsDirSum`/`rsFindSyng` rskip readers never write
`sb->node/status/path/loc` (all mutable query state — the
`SyngBWTpath`, the `SeqhashIterator` — is a fresh per-call
allocation); `startCount(sb,k,false)` is a hashFind, and the lazy
`startHash` create is the add path only. (2) the KmerHash query path
(`kmerHashFindThreadSafe`, every bind-phase syncmer lookup) writes
ONLY the `kh->finds`/`kh->deltas` diagnostic counters — two U64s in
the small-heap KmerHash *header chunk*, a benign cross-thread race
since the seam-thread era: the COW cost is ONE page, and it is not in
either big arena; `kmerHashFind`'s scratch write into `pack[max+1]`
is the index-build path only. (3) the Seqhash is read-only
post-create; the iterators are per-query allocations. (4) the
Rust-side lazy structures do not fire: `locate_checkpoint_index` is
never built (the panel's `.spos` serves locate), and the C
diagnostics (`pathCount`) are set once at load — pre-fork.

*THE PREFORK HAZARD THE FIRST PROBE RECEIPTED (the run hung — CPU
frozen at 05:09:26, all threads S-state, ~20 minutes past the last
locus-open line):* rayon's GLOBAL registry is set ONCE per process
(`THE_REGISTRY_SET` — never replaceable in-process), and a registry
initialized in the PARENT holds worker threads that DO NOT EXIST in
a forked child; the call's seam threads run inner par_iters that
install the GLOBAL registry (the committed inner behavior), enqueue
to dead workers, and hang forever. **THE FIX BY CONSTRUCTION:** the
parent NEVER touches rayon before the fork — `SyngIndex::load` is
verified rayon-free (the C reader's own joined pthreads, mmap'd
sidecars, plain file reads) and the parent builds NO AgcIndex (whose
build sizes its slots with `rayon::current_num_threads()`); every
child builds its OWN AGC index (the struct is small — name maps +
per-thread decompressor slots over the page-cache-shared AGC file;
the per-child build is <1s, receipted), so the child's global
registry initializes at the same sequence point as the committed
fresh process, spawning the child's own 64 threads at the inherited
`RAYON_NUM_THREADS`. Every par_iter in the child — the main thread,
the cut, the call's seam threads — joins the child's own LIVE
registry: the fresh-process trajectory by construction, no
`install()` indirection. This is the read-length-class lesson in its
fork form: the audit names every process-global a child inherits
(driver doc comment, verbatim list), and the one that mattered was
found by the gate, not assumed away.

**THE DRIVER (`sample_projection fork-cut-call-batch-global`):** the
parent loads the panel ONCE (the load owner — 197-221s, the co-tenant
band of the committed 241s panel-load pole) and forks;
`--rows-per-child` chunks the manifest (1 = one child per territory,
the maximal-share form; the manifest's own batch width amortizes the
cut scan within a child — the fleet-shape parameter, never an
algorithm constant). Each child builds its own AGC, runs its chunk's
batch cut (the batch-composition-invariance receipt: single ==
batch member) plus its chunk's calls (the combined process's own
`run_combined_call_rows`, extracted verbatim), then `_exit`s; the
parent waits and fails at the first failing chunk in order. The smaps
receipts print at parent-post-load, child-open, post-cut
`malloc_trim` (before/after private), and child-close with
per-mapping deltas and COW-BREAK marking over every >=1GiB mapping.

**THE GATES (all PASS, the two-sided committed bar):**

| form | wall | byte-identity | COW-hold (per child) |
|---|---:|---|---|
| N=1 probe (L1, rows-per-child 1) | load 213s + child 496s (cut 354s incl. scan 153s; call 129s) | **PASS** vs committed | open private 0MB shared 32.36GB; close: the 12.87GB GBWT arena **private 0.00GB shared 12.06GB**, the 4.31GB khash arena **private 0.00GB shared 4.31GB**; own private +40.24GB |
| N=2 probe (L1 x2, 2 concurrent children) | load 198s + children 541s | **PASS** (both children) | private +36.12GB / +37.13GB (own work; arenas shared) |
| L1-L4 at rows-per-child 1 (4 concurrent children) | load 221s + children 841s (cuts 572-693s, scans contended 378-485s; calls 117-137s) | **PASS** (12 cut + 8 call files + all sidecars vs BOTH baselines; instrument field-identical mod timing) | arenas private 0.00GB shared 12.13/4.31GB (Pss 3.03/1.08GB = the /4 share); per-child private +35.3-38.8GB |
| L1-L4 in ONE width-4 child (rows-per-child 4 — the combined form minus the load pole) | load 197s + child 666s (cut 415s incl. scan 141s; 4 concurrent calls 190-229s) | **PASS** (the same two-sided bar) | arenas private 0.00GB shared (Pss /2 with the parent); child private peak ~72.3GB, +52.81GB at close |

**THE VERDICT ON THE §10.6.9 NAMED RISK:** the COW holds — over a
whole cut+call, at every measured concurrency, the panel's arenas
break ZERO pages (the GBWT and khash arenas read `private 0.00GB,
shared 12.13/4.31GB` at close; the only private growth on an
inherited mapping is the `.pstep` mmap's +0.24GB file-backed
accounting class). The one-page khash-counter break predicted by the
code audit is under the print precision of the arena lines — the
design's viability number is the measured per-child private below,
not COW.

**THE HONEST NEW FINDING (the receipts the §10.6.9 projection
missed):** the child's own private growth is +35-40GB per territory
at rows-per-child 1 — NOT COW, and NOT the projected ~13GB: it is the
CUT SCAN-SIDE working set (the shard decode/remap stream through 64
rayon-thread arenas — 734M reads scanned per child, WIDTH-INSENSITIVE:
the width-4 combined form held the same ~42GB class), retained by the
allocator after the cut phase. The `malloc_trim(0)` hygiene receipt
measures it honestly: 35.2-38.95GB -> 32.16-32.75GB (returns only
3.1-6.2GB — the free TAIL; the ~32GB remainder is fragmentation/live
class, not returnable). The width-4 child amortizes the scan side
over its 4 territories (private peak ~72.3GB, ~18GB/territory).

**THE BIND-CONTENTION RECEIPT (the §10.6.7 question answered):** at
the fork form with per-process decompression pools, the contention
VANISHES — 4 concurrent children (one call each) bind at **46.1-55.6s
each**, the serial-run class (50-73s), versus the combined process's
contended **137.2-180.0s** (~3x): the 773.7s summed-thread-time
contention was the in-process shared-bandwidth class, and per-process
pools delete it at 4-way process concurrency. (Within ONE width-4
child the 4 concurrent calls contend as before: binds 116.7-152.0s —
the in-process class is unchanged, as expected.)

**THE CEILING, CORRECTED FROM THE RECEIPTS (the ~74-slot projection
corrected):** at the 1TB model with the 900GB guard, the arithmetic is
parent ~34GB private + per-child private peak: at rows-per-child 1,
(900-34)/40 -> **~21 territory slots** (WORSE than the combined
form's 32 — the scan side re-binds per child); at rows-per-child 4,
(900-34)/72.3 -> **11-12 children x 4 = 44-48 territory slots** — the
honest receipted ceiling, ~1.4x the combined form's 32, NOT 2.3-2.5x:
the scan-side working set is the term the §10.6.9 model missed. THE
CPU CAVEAT the receipts also force: 4 concurrent children x 64
threads ~ the box's 256 cores — the 4-way scans contended 378-485s vs
141-153s solo (2.5-3.2x) — so the fleet at >~4 children needs child
pool sizing (an operational fleet parameter) or the walls inflate;
the 44-48-slot memory ceiling is not reachable at the committed
64-thread child shape.

**THE GENOME WALL, RESTATED FROM THE RECEIPTS (the §10.6.10 fleet
convention — solo walls, memory-capped slots, no term extrapolated
beyond its receipt):** at rows-per-child 4, 44-48 slots at a 666s
wave (the receipted width-4 child, load excluded — the parent's
197-221s load amortizes to ~0 at fleet scale under the same
convention the 241s pole was charged per process): 5,844 territories /
44-48 slots x 666s = **~22.5-24.6h**, versus the combined form's
~45-47h — **~1.9-2.0x**, the load pole deleted and the slot ceiling
44-48 vs 32. (The parent-level admission control — forking at most
the ceiling's children per wave and looping the manifest — is the
named remaining driver increment; the current form's manifest IS the
wave.) Per-territory marginal at the fork form: 666/4 = **~167s
ex-load**, unchanged from the combined marginal — the panel share is
a SLOT and LOAD lever, not a per-territory wall lever, exactly as
§10.6.9 stated.

**THE MARATHON RE-RECEIPT (duty 6 — the measurable parts receipted,
the unmeasurable named, no marathon wall number fabricated):** the
substrates are reaped, so the re-receipt is the retained-shape
machine: (a) THE MEMORY BOUND re-run, `marathon_worst_locus_shape_
memory_receipt` — at T03803's profiled shape (532 x 4,438,575 at the
roster's 2.1% pinned share): compressed live 750MB (0.33 B/pair) vs
the dense form's 29,275MB (13 B/pair) = **39.0x**, measured RSS delta
1,333MB — candidate C's committed bound reproduces. (b) THE
SCAN-WALL CORE, the new `marathon_worst_locus_shape_scan_wall_
receipt` — the merge's equal-column hash pass (532 rows, 2.36e9
(fold,unit) hash steps) in 25.4s = **10.75 ns/(fold,unit)**
single-threaded; the FULL exhaustive class scan at the all-distinct
worst shape (141,778 fold pairs x 4,438,575 units = 6.29e11
unit-pairs) in 204.7s at 256 threads = **0.33 ns/unit-pair aggregate
(~84ns single-thread equivalent)**. NOT MEASURED, because the reaped
substrates carried them: the pair loop's PINNED alignment class (the
1.3-2.7us/pair band needs real records/reads), the bind phase, the
records/census phases — and therefore any whole-locus marathon wall.
The honest bound the new receipt buys: the class-scan core at the
no-collapse worst shape is 204.7s at full-box width — and the
committed single-threaded marathon receipts (T03803's whole locus
717.4s) imply the real loci's class scans run over far fewer than
532 unmerged classes (the merge collapses first), so the synthetic
all-distinct scan is an UPPER-shape bound, never a wall claim.

**THE 20-MINUTE BUDGET, RESTATED AT INCREMENT 4 (supersedes
§10.6.10's fleet terms):**

- **map ~3 min**; **prep ZERO**; the panel-side routes inventory
  1,505s once — unchanged.
- **cuts + calls at the fork form (rows-per-child 4, all gates
  PASS):** per-territory marginal **~167s ex-load** (the width-4
  child: cut share ~104s + call phase ~62s at 4-concurrent), the
  load pole DEAD (197-221s once per parent), the bind at the serial
  class 46-56s under 4-way process concurrency.
- **THE FLEET MATH:** 44-48 memory slots at 666s waves = **~22.5-24.6h
  genome wall** — with the CPU caveat: >~4 concurrent children at
  the committed 64-thread shape oversubscribe the 256-core box (the
  4-way scan contention receipt: 2.5-3.2x), so the slot ceiling needs
  child pool sizing to be USABLE, an operational fleet decision.
- **THE REMAINING GAP:** 5,844 territories in 1,200s needs 4.87
  territories/s; the fork form receipts 44-48 slots / 666s = 0.066/s —
  **a ~73x gap**, every pole named with its receipt: (1) THE LOAD —
  dead; (2) THE CUT SCAN — the NEW pole (141s scan + ~104s cut share
  per width-4 child, plus the ~32-40GB scan-side working set that
  binds the memory ceiling at ~18GB/territory; the same 65GB shard
  gzip is re-decoded per child — a shared-scan/bulk-prefetch form is
  the undesigned candidate); (3) THE CALL'S OWN WALL — bind 46-56s
  (the ~1.4-1.5M per-fetch AGC seeks, the named floor), the LOCI CORE
  — **candidate B (the always-on exactness assertion) and the
  census-superset are THE OWNER'S decisions**, the receipts say the
  path to 20 runs through them; (4) THE MEMORY CEILING — 44-48 slots
  (corrected from the 74 projection); (5) THE MARATHON TAIL — the
  kill class collapsed by C (re-receipted 39.0x), the WALL not
  re-receiptable (substrates reaped; the measurable scan-core rates
  above); (6) THE CPU SATURATION — the honest new term: the 64-thread
  child shape saturates at ~4 children.

**THE EXACT RESUME POINT:** the fork-after-load panel share is
landed, gated, and receipted (this section); the remaining increment
order: (1) the parent-level wave admission control (fork at most the
ceiling's children per wave, loop the manifest — the fleet driver's
last piece) + child pool sizing as the operational fleet decision;
(2) the cut-scan pole (the shared-scan/bulk-prefetch form,
undesigned); (3) candidate B + the census-superset (the owner's);
(4) the bind's per-fetch AGC-seek floor (a bulk-prefetch form,
undesigned); (5) the marathon re-receipt's full form at the next
campaign window (determinism + wall, if substrates exist).

### 10.6.13 THE UNCLAIMED ENGINEERING LEVERS, MEASURED — THE
WIDTH-AMORTIZATION CURVE, THE SEGMENT-ALIGNED SEEK FIX (GATE
PASS), THE THREAD-RUNG TABLE, AND THE CPU-WORK-PER-TERRITORY
ARITHMETIC THAT CORRECTS THE GENOME WALL (increment 5; commits
f259c47 the instruments + 511834e the fix; every gate PASS at the
committed two-sided bar)

**DUTY 1 — THE SCAN-AMORTIZATION CURVE (the biggest unclaimed CPU
item, measured at widths 8/16/32; gates PASS byte-identity at every
rung, work dirs forkcutcall-width-W{8,16,32}):** only 4 distinct
gate territories exist, so the width-W manifest carries W/4 COPIES
of each (the predecessor's own N=2 COW-probe form, L1x2 —
the predecessor's own gate-PASS precedent) and only the FIRST L1
copy carries a call (the cut side is the measurement; the call stays
the committed 1-call shape). THE HONESTY LABEL the receipts force:
the copies route every touching read to W/4 copies of its territory,
so the scan's REMAP term is measured at its per-copy cost (an upper
bound on the genome-scale marginal for distinct territories — the
committed receipts already showed 1x→4x flat, and 4x→32x adds only
+28-38s); the DECODE share is measured exactly (it is
width-insensitive by construction).

| width (rows/child) | scan | cut total | cut/row | child wall | child/row | pre-trim private | post-trim private |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 (committed) | 135-153s | 354s | 354s | 496s | 496s | 35.2-38.95GB | 32.16-32.75GB |
| 4 (committed) | 141s | 415s | 104s | 666s | 167s | 38.95GB | 32.75GB |
| 8 | 151s | 519s | 65s | 659s | 82s | 44.98GB | 32.72GB |
| 16 | 147s | 662s | 41s | 808s | 50.5s | 57.28GB | 32.69GB |
| 32 | 179s | 968s | 31s | 1116s | 34.9s | 81.75GB | 32.66GB |

THE CURVE: the scan's decode share is WIDTH-INSENSITIVE (141-179s
across 1x→8x copies — the remap term costs ~4-5s per added copy)
and the per-territory cut share falls 354s → 104s → 65s → 41s →
31s; the ~32.7GB post-trim live/fragmentation class is FLAT at
every width (the scan-side retained set does not grow with rows;
the PRE-trim peak grows ~+1.5GB/copy — the per-copy cut-state
arenas — so the child's own memory model at width W is
~32.7GB live + ~1.5GB/row transient + one call's ~6-10GB). The
CEILING arithmetic at width 32: (900-34)/~82GB pre-trim ≈ 10
children x 32 = ~320 territory slots — but the CPU receipt below
says the box holds 4 children at the committed 64-thread shape,
so MEMORY IS NO LONGER THE FIRST-BINDING TERM at width >= 8; the
CPU is.

**THE SCAN-DECODE SERIAL RECEIPT (duty 1's second half, answered by
measurement):** the scan decodes 19 shard gzip streams, ONE flate2
reader thread per file (shard_packet_stream's spawned reader), and
zcat of the biggest shard alone (5.25GB compressed -> 12.9GB) takes
80s single-streamed (~160MB/s) — the scan wall 141-151s is the
slowest files' GUNZIP wall, not parallel compute: the decode uses
at most 19 cores of the 256-core box and IS SERIAL WITHIN EACH
FILE (the same per-stream decompressor disease as the bind's slot
class, in its format form). THE FIX CANDIDATE, NAMED, NOT
IMPLEMENTED: a blocked/frame-indexed shard format (bgzip/zstd
frames at emit time) would parallelize the decode within a file —
an OUTPUT-FORMAT change to the projection's own product machinery,
beyond this increment's seam; the honest receipt is that the scan's
decode pole is gzip-format-bound at ~80-100s per wave no matter the
width or the thread shape.

**DUTY 2 — THE BIND SEEK-FLOOR ANATOMY, MEASURED, AND THE FIX
(GATE PASS, commit 511834e):**

*THE ANATOMY (the f259c47 probe, bind_seek_floor_anatomy_receipt,
explicit-run at the panel AGC):* lane CHM13#0#chr1 = 248,387,328bp
over 4,101 segments (~60.6KB mean, kmer 31), and the per-fetch cost
of fetch_sequence is a FIXED FLOOR, NOT byte-proportional: 200
random 900B crops = 0.4 ms/fetch (2.3 MB/s), 200 random 6KB crops
= the SAME 0.4 ms/fetch (15.2 MB/s) — 6.7x the bytes for the same
wall; 64 sequential 64KB crops = 0.24 ms/fetch (316 MB/s); one
250KB crop = 1.3 ms. THE MECHANISM, named in ragc-core's own code:
get_contig_range pays PER OVERLAPPING SEGMENT, EVERY FETCH — a pack
part re-read + whole-pack zstd re-decompress (get_part_by_id has
no cache; a delta-stream pack holds 50 LZ-encoded segments) + an
LZDiff REFERENCE PREPARE (a reference clone + hash-table build) +
the LZ decode. So the bind's ~6KB-mean lane loads each paid the same
floor a whole 64KB segment would: the seek floor was the FETCH
COUNT times ~0.4ms, never the byte volume. THE PRODUCTION RECEIPT
(the pre-fix IMPG_LANE_TRACE at the T32 rung): p50 per-load = 0.48ms
(the floor, confirmed in production), 7,688-111,212 bind loads per
call over exactly 466 lanes (the schedule-dependent union geometry
the committed receipts already named), sizes mostly <2KB.

*THE FIX — SEGMENT-ALIGNED SPAN GROWTH (a layout change, derived
from the archive's own decompression unit, never a tuned
constant):* LaneCache::span now grows the cached span to the
touched SEGMENTS' own boundaries on every fetch (initial AND both
delta-extend sides) — the boundaries come from ragc's PUBLIC
descriptors (get_contig_segments_desc, the same contribution
formula get_contig_range itself uses) through
AgcSegmentLayout/Sources::aligned_crop (identity for
non-segment-addressed FASTA sources, so the committed local frame
is untouched). BYTE-IDENTITY, the delta-extend argument verbatim:
the covered span is still a contiguous [lo,hi) whose buffer is
exactly fetch(lo,hi) by induction (a crop of a fetch answer IS the
fetch of the crop — the Sources contract), and every served crop is
a slice of it; the alignment changes only WHICH superset range is
fetched, never a served byte.

*THE FIX GATE (the committed fork shape, rpc=1, T=64, 4 children;
work dirs forkcutcall-threads-T64 + forkcutcall-seekfix-T64):* ALL
4 cut + ALL 4 call receipt files BYTE-IDENTICAL vs BOTH baselines —
GATE PASS; the served-crop counters identical (1.41-1.48M fetches,
349.9-358.9MB — the SERVED bytes unchanged, exactly the argument).

*THE FIX RECEIPTS:* the bind-phase lane loads 18,626-107,333 (the
committed class) → **1,468-2,200 per call — the projected 30-50x
collapse — AND THE SCHEDULE VARIANCE ITSELF DELETED** (the count is
now the lane's touched-segment hull, not the thread interleaving:
loads/lane 3-5, vs the pre-fix 33-225); per-load size 63-95KB (the
whole-segment class, p50 1.0-1.7ms — the same fixed floor, now
amortized over the full segment); fetched span volume 111MB →
133-160MB/call (the segment hull, +20-45% — the memory cost,
receipted: child rss-at-binding 75.7-78.1GB vs the committed
77.5-78.3GB class, invisible at scale); binds **36.6-51.2s vs the
committed 46.1-55.6s band** at the same shape and co-tenancy; the
load thread-time fell 100-165s → 44-86s summed. THE RESIDUAL,
NAMED (pre-existing, NOT the fix's — the pre-fix traces carry the
same class): a handful of multi-second fetches per call (5-21s
each, ~3-8s of bind wall) — the fresh decompressor pool's FIRST
fetch pays ragc's whole-archive contig-batch metadata load + cold
page cache, under the per-lane mutex; the POOL-WARM form at the
child's own AGC build (one get_contig_range touch per pool, off
the bind's critical path) is the named next candidate, NOT
implemented.

**DUTY 3 — THE LOC-PARALLEL THREAD-RUNG TABLE (the per-child thread
width audit; gates PASS at every rung, work dirs
forkcutcall-threads-T{32,16,8}; the committed T=64 rung is the
SS10.6.12 receipt):** the fork form, rpc=1 (4 concurrent children,
one territory + one call each), re-run at RAYON_NUM_THREADS 32/16/8
— the audit's question is territory throughput per CORE-SECOND,
not per-process:

| child threads T | box threads (4xT + 4x16 seam) | 4-way scans | cut phases | children | binds | territories/core-second |
|---|---|---|---|---|---|---|
| 64 (committed) | 320 (OVERSUBSCRIBED) | 378-485s | 572-693s | 725-841s | 46.1-55.6s | 4/(841x320) = 1.5e-5 |
| 32 | 192 | 264-286s | 408-454s | 540-604s | 49.1-59.2s | 4/(611x192) = 3.4e-5 |
| 16 | 128 | 147-164s | 279-340s | 411-492s | 47.1-60.8s | 4/(498x128) = 6.3e-5 |
| 8 | 96 | 119-159s | 268-414s | 401-580s | 47.1-61.6s | 4/(585x96) = 7.1e-5 |

THE FINDINGS THE RECEIPTS FORCE: (a) **THE SS10.6.12 SCAN
CONTENTION (2.5-3.2x AT 4-WAY) WAS OVERSUBSCRIPTION, NOT A
MEMORY-BANDWIDTH WALL** — at T=16 the same 4-way scans run at the
SOLO class (147-164s vs 135-153s solo), and even T=8 holds it
(119-159s): 320 threads on 256 cores was the whole disease; (b)
the scan is thread-width-INSENSITIVE below 32 (it is gunzip-IO-
bound, at most 19 streams — the duty-1 decode receipt); (c) the
bind is SEAM-width-bound (16 scoped threads per call), insensitive
to the rayon width (47-62s at every T); (d) the per-core-second
throughput RISES as the width falls (1.5e-5 → 7.1e-5, 4.8x), because
the oversubscribed rung burns cores waiting. THE DERIVED SHAPE at
the box: children x T <= 256 with T in the 8-16 band avoid the
oversubscription class entirely — but the WIDTH receipts above say
the territory count per child should be wide (the scan share), so
the honest fleet shape is WIDE children at NARROW threads; the
W-x-T interaction (e.g. width 16 at T=16) is the named unmeasured
cell of the matrix.

**DUTY 4 — THE BUDGET, RESTATED AT THE CPU-WORK-PER-TERRITORY
ARITHMETIC (the correction this increment's receipts force on the
SS10.6.12 headline):** the box-saturation arithmetic made explicit
— a fleet cannot beat WORK/CORES, and the SS10.6.12 genome-wall
number (44-48 slots x 666s waves = ~22.5-24.6h) charged the SOLO
666s wall to 11-12 CONCURRENT 64-thread children (704-768 threads
on 256 cores, the receipts' own caveat) — it was a memory-ceiling
extrapolation, not a CPU-honest one.

*The CPU-work-per-territory table (wall-charged: phase wall x
threads configured; the cut at 64 rayon, the call at 64 rayon + 16
seam per concurrent call):*

| form | cut core-s/row | call core-s/call | per-territory (real shape: 1 call each) |
|---|---:|---:|---:|
| width 1 (solo) | 354x64 = 22,656 | 129x80 = 10,320 | ~32,976 |
| width 4 (committed) | 415x64/4 = 6,640 | 229x128/4 = 7,328 | ~13,968 |
| width 32 (this increment) | 968x64/32 = 1,936 | ~5,000-7,300 (the bracket: phase-timer CPU ~3,200-4,700 true + inputs/context; wall-charged 7,328 at the 4-concurrent packing) | **~6,900-9,200** |

THE HONEST HEADLINE: **the width amortization collapsed the CUT
(from ~6,640 to ~1,936 core-s/territory — 3.4x) but the CALL never
amortized, and at wide rows the CALL IS NOW ~75-80% OF THE
PER-TERRITORY CPU WORK.** The genome-wall floor at the receipts:
5,844 territories x ~6,900-9,200 core-s / 256 cores = **~44-60h at
the width-32 shape** (vs the width-4 committed shape's own honest
floor, 5,844 x ~13,968 / 256 = ~88.6h — the SS10.6.12 headline of
22.5-24.6h was ~3.7x optimistic, exactly what its own caveat
said); the wall-charged wave arithmetic agrees (width-32 children:
cut 968s + a 16-32-call in-child phase, 4 concurrent children,
~45-60h). NO number in this table is extrapolated beyond its
receipt except the two labeled terms (the 16-32-call in-child
concurrency wall, unmeasured — the committed 4-concurrent receipt
is the rung; and the whole-child 4-way contention at width 32,
the width-4 receipt's 1.26x class).

*THE 20-MINUTE FLOOR:* 5,844 territories x X core-s / 256 cores =
1,200s -> **X = 52.6 core-s/territory**. Today's honest
per-territory work at the best receipted shape is ~6,900-9,200
core-s — **the remaining gap is ~130-175x**, and the earlier
scratch-table figure (~2,237 core-s/territory, ~43x) is CORRECTED
HERE as an artifact: it divided ONE call's cost by 32 rows; every
territory needs its own call, and the call is the dominant term.

*THE POLES, each labeled with its receipt:*

- **THE CALL (~5,000-7,300 core-s/territory, the dominant term) —
  split:** the LOCI CORE (loci walls 30-86s at 64 threads ≈
  2,000-5,500 core-s — **candidate B (the always-on exactness
  assertion) and the census-superset are THE OWNER'S decisions**,
  OWNER-DECISION); the BIND (~880-820 core-s of true summed-thread
  work: steps 246-461s + verify 418-539s at 16 seams — the seek fix
  landed, the pool-warm residual ~3-8s named, ENGINEERING); the
  INPUTS phase (41-50s, mostly single-threaded census/maps/AGC
  loads — an IO class, ENGINEERING); the context/rows assembly
  (11-19s, serial — ENGINEERING).
- **THE CUT at width 32 (~1,936 core-s/territory, wall-charged):**
  the scan decode (~358 core-s — gzip-format-bound, the
  blocked-shard-format candidate, ENGINEERING), the per-territory
  setup/emission/cut-state classes (~1,578 — the remap/records
  classes, ENGINEERING).
- **THE BOX-LIMIT:** 256 cores — the floor arithmetic's own
denominator; the oversubscription receipts say the fleet must run
children x T <= 256 (the T-rung table), and at width >= 8 the
memory ceiling no longer binds first (pre-trim peaks 45-82GB/child:
4-10 children fit the 1TB model).

**THE EXACT RESUME POINT (increment 5 closes complete):** the four
duties landed and gated — (1) the width-amortization curve + the
scan-decode serial receipt (the blocked-shard format the named
unimplemented candidate); (2) the seek-floor anatomy + the
segment-aligned fix (GATE PASS, the load collapse receipted, the
pool-warm residual named); (3) the thread-rung table (the
oversubscription finding; wide-children-at-narrow-threads the
derived shape; the W-x-T cell the named unmeasured term); (4) this
budget restatement. THE REMAINING ORDER: (a) the W-x-T fleet cell
(width 16-32 children at T 8-16) — one measured wave closes the
last labeled extrapolation; (b) the in-child call-phase packing at
wide rows (16-32 concurrent calls, the committed 4-concurrent
receipt's rung); (c) the pool-warm form; (d) the blocked-shard
format (an output-format change, the owner's seam); (e) candidate B
+ the census-superset (the owner's — the receipts now say the path
to 20 minutes runs THROUGH the call, which is 75-80% of the honest
CPU work); (f) the marathon re-receipt's full form at the next
campaign window.

### 10.6.14 THE CALL-CORE AUTOPSY — the per-class, per-locus CPU-work
census, measured (the owner's "why is the call core so slow — is it
just some slow loci?" answered at the arithmetic the receipts force;
the §10.6.13 WALL-CHARGED call term corrected by direct serial
measurement)

**THE ANSWER IN THE FIRST PARAGRAPH.** At TRUE CPU — the serial rung
(RAYON_NUM_THREADS=1, no seam), where every wall second IS a
core-second — an ordinary L1-class territory's whole call costs
**~203 core-seconds ex-load, ~8.1 per locus over its 25 loci**, and
the classes are: **the scoring pair loop 84.3 core-s (42% of the
call; inside it the ALWAYS-ON EXACTNESS ASSERTION is 32.9 core-s =
16% of the whole call and 58% of the pinned alignment work — MORE
than the ~½ receipted in §10.6.3)**, **the exhaustive class scan
68.4 core-s (34%; folds²/2 × units at a constant ~18
ns/unit-pair)**, the fold-build/units/qual classes 12.2 (6%), and
the serial bind+rows+IO classes 37.9 (19%; the bind itself is now
**6.6s** — the seek fix collapsed the serial bind 1013.8s → 6.6s,
153x). AND YES, THE TAIL IS THE STORY INSIDE THE LOCI CORE: the
worst 5% of loci (2 of 25) own **64.9% of L1's loci phase** (65.9%
at the committed serial receipts, 68.5% at L3, **88.1% at L4**),
because the fold×units cross product is 100-360x bigger at the two
or three cross-product-heavy loci than at the median locus. BUT the
130-175x headline itself was MOSTLY THE WALL-CHARGING CONVENTION:
the fleet shape's call phase runs at ~1-4 cores of true work
against 128+ configured threads (the §10.6.13 T-rung table's own
idle finding), so charging phase walls × configured threads
overcounted the call's true work ~25-36x. **The honest call-core
gap to the 20-min line is ~3.9x today, ~2.0x after the owner's
candidate-B + census-restriction decisions**, and the residual is
owned by the class scan (OWNER-DECISION: the byte-identity
accumulation), the pair loop's per-base comparison
(IRREDUCIBLE-FOR-THE-MODEL at a ~0.2 µs/pair floor vs the 0.69 µs
measured placement, plus the re-anchorable anchor class), and the
serial IO/context classes (ENGINEERING).

**THE INSTRUMENT AND ITS GATES (receipt-side, DEFAULT OFF, no
behavior change).** `IMPG_LOCUS_SCORING_CENSUS=1` wraps the per-locus
scoring pair loop's four work classes in per-call Instant timers,
atomically accumulated across the fold par_iter, emitted as a
`scoring_census` member of the per-locus receipt + a per-locus
stderr line: the per-orientation ANCHOR/PIN BUILD
(`anchor_correspondence` + the merged backbone), `placement_score`
(the pinned pair's own alignment), the ALWAYS-ON `direct_score`
exactness assertion, and the marginalization + row build. THE RUN
(work dir `global-call-L1-autopsy`, the exact committed L1-fix
serial recipe at the current binary): exit 0, wall 476s, and
`calls.jsonl` + sequences **BYTE-IDENTICAL to the committed L1-fix
receipt** — the no-behavior-change gate. THE CONTROL (the same
recipe, env unset: `global-call-L1-autopsy-control`): exit 0, wall
439s, **byte-identical again** (a free determinism re-gate), loci
phase 164.8s vs the census run's 182.1s — **the instrument's own
timer overhead is +17.3s (+10.5%), receipted**; every census class
below is scaled by 164.8/182.1 = 0.905 to the no-timer walls (the
per-locus raw census lines stay upper bounds; the receipts are in
the work dir's `autopsy-census.txt` / `autopsy-final.txt` /
`autopsy-mining.txt`).

**(1) THE PER-LOCUS CENSUS AT AN ORDINARY LOCALITY (L1, serial =
true CPU, current binary, timers scaled out).** The call's phase
summary (the control run): inputs 236.0s (≈ a ~222s panel load +
~14s residual IO), bind 6.6s, quality 2.7s, variants 0.2s, rows/
context 14.3s, **loci 164.8s** — the loci phase is now **81% of the
ex-load call CPU**. The loci phase's anatomy, summed over the 25
loci and decomposed to the pair loop's own classes:

| class | core-s | share of loci | the measured rate |
|---|---:|---:|---|
| scoring pair loop | 84.3 | 51% | see the four rows below |
| — anchor/pin build | 16.8 | 10% | 54,236,474 calls × **0.31 µs** |
| — placement_score | 24.0 | 15% | 34,728,402 calls × **0.69 µs** |
| — the direct_score assertion | 32.9 | 20% | 34,697,264 calls × **0.95 µs** |
| — marginalization + row build | 3.0 | 2% | over 48.8M pairs × 0.06 µs |
| — pair-loop residual (per-unit/no-pin pushes) | 7.6 | 5% | |
| the exhaustive class scan | 68.4 | 41% | 3.711e9 unit-pairs × **~18.4 ns** |
| fold build + units classification | 8.3 | 5% | |
| qual | 2.9 | 2% | |
| exactness (the sampled receipt) | 1.0 | 1% | |
| (anatomy: env-off, 0.0s at every locus) | 0.0 | — | |

48,833,152 (fold × unit) pairs at L1, 71% pinned (34.7M
`factorized_checked` — the assertion's own count). THE BIND
DECOMPOSITION at the current binary: steps 3.4s (3,775 record
extractions) + verify 3.0s (4.35M probes) + lane loads 2,318 /
133.5 MB / 2.3s of load time, phase 6.6s — the seek fix's serial
form; THE FLEET-SHAPE BIND WALLS (36.6-60.8s at the T-rungs) are
WAITING, not work: 16 seam threads × their own per-thread AGC
decompressor pools each pay the pool-warm cold-start class (the
§10.6.13-named candidate, 3-8s each, under the per-lane mutex) —
the POOL-WARM form deletes them. The rows/context class (14.3s)
and the inputs residual (~14s warm-cache serial; 41-50s at the
fork-child rung with cold sidecars) are serial IO classes.

**(2) THE "SOME SLOW LOCI" QUESTION, answered quantitatively (the
per-locus WALL = serial CPU distribution, all five lanes' committed
serial receipts + the current binary).** The marathon class is
15-25 LANES per genome at the roster's worst shapes (§10.6.3's
committed bound); the WITHIN-territory tail
is ALSO concentrated — but it is NOT the marathon class: these are
ordinary-territory loci whose fold×units cross product is large:

| lane | n | median locus | p90 | max | worst-5% share of the loci phase |
|---|---:|---:|---:|---:|---:|
| L1 (committed serial receipts) | 25 | 0.60s | 14.54s | 69.83s | **65.9%** (top-1 42.7%) |
| L1 (current binary, autopsy) | 25 | 0.66s | 16.03s | 76.91s | **64.9%** (top-1 42.2%; = 53% of the ex-load CALL CPU) |
| L2 | 25 | 0.59s | 6.78s | 9.00s | 35.1% |
| L3 | 24 | 0.51s | 3.50s | 64.15s | **68.5%** |
| L4 | 25 | 0.45s | 1.41s | 101.22s | **88.1%** (top-1 85.0%) |
| pilot | 22 | 14.93s | 37.36s | 187.44s | 48.6% |

So the direct answer: **yes — at an ordinary territory, 2 of 25
loci own ~2/3 (L1, L3) to ~9/10 (L4) of the loci-phase CPU, and the
loci phase is 81% of the call's ex-load CPU — the worst-5% loci own
~53% of L1's whole ex-load call core.** The tail loci are the
cross-product-heavy ones (L1 locus 18: 111f × 180,179u = 20.0M
pairs; L4 locus 23: 155f × 191,360u = 29.7M; the median locus is
157f × 530u = 83,210 pairs — a 240x pair-count spread), NOT the
marathon class (whose roster shape is 532f × 4.44M u = 2.36e9 pairs,
~100x L1's worst locus).

**(3) THE MULTIPLIER STACK, each with its receipt.**

**(a) THE ALWAYS-ON EXACTNESS ASSERTION (candidate B's target).**
Measured directly: **32.9 core-s of L1's 164.8s loci phase = 20% of
the loci phase, 16% of the ex-load call, and 58% of the pinned
alignment work (placement 24.0 + assertion 32.9)** — MORE than the
§10.6.3 "~½" receipt (per-call: 0.95 µs assertion vs 0.69 µs
placement — `direct_score` re-verifies the backbone base-by-base
with the per-base cover search while `placement_score` trusts node
identity). B (assertion as a receipt-only sample) is an OWNER
receipt-shape decision worth **~33 core-s/territory (−16% of the
call)** at L1-class, and proportionally more at pinned-heavy loci
(locus 18: the assertion is 17.8s of the ~69.6s locus, scaled).
(EXECUTED — §10.6.15: the owner's Ruling 1 removed it from the
production path, byte-identity gated at L1-L4, the returned
per-locus budget measured at 6.66 vs the 6.79 projected here.)

**(b) THE CENSUS-SUPERSET 2.8x MASS, reconciled per-locality.** The
§10.5.4 shadow receipt (1364s vs 1450s, "~6% wall") was measured
against the OLD 1013.8s serial bind — at the CURRENT bind the same
shadow arithmetic says the mass is now **a ~1.6x multiplier on the
call's ex-load CPU at L1** (202.7 → 124.3 core-s, arithmetic below).
WHERE IT BITES, per-locality (the shadow's per-locus receipts, L1
serial): the mass bites exactly at the loci where the added reads
carry NEW routed units — locus 8: units 14,711 vs 1,490 (9.9x) →
scoring 8.4x; locus 18: 180,179 vs 40,092 (4.5x) → scoring 4.0x;
locus 11: 6.0x units → 5.4x scoring; locus 13: 8.4x → 2.7x — and
at the SATURATED loci (units 1.0x: loci 3, 6, 9, 10, 14, 17, 20,
...) the scoring is 1.0x, exactly the §10.5.4 saturation statement
at per-locality resolution. The per-class shadow ratios (both
serial, same binary era): scoring 0.520, class scan 0.496, fold
residual 0.666, qual 0.909, whole loci phase **0.524** (85.7s vs
163.5s); the bind/rows/inputs are mass-invariant at saturation.

**(c) THE FOLD×UNITS CROSS-PRODUCT at ordinary vs segdup loci.**
Ordinary (median) locus: 157f × 530u = 83,210 pairs — the whole
locus is 0.66s. The cross-product tail loci: 111-231f ×
14,711-191,360u = 3.4-29.7M pairs. The measured per-pair rates
(census run, scaled): the pair loop runs **2.06-2.18 µs/pair at the
pinned-heavy loci (73-105% pinned)** and **0.18 µs/pair at the
1%-pinned locus 12** — the no-pin floor class, matching the
marathon band's 0.2-0.3 µs/pair receipt; the class scan is a
CONSTANT **16.4-20.2 ns/unit-pair** at every locus (3.711e9
unit-pairs at L1), quadratic in folds (pairs = folds²/2 × units).
So the cross-product's CPU bites through BOTH terms: the pair loop
scales with PINNED pairs (units × folds × pin rate), the class scan
with ALL pairs.

**(d) THE ALIGNMENT CLASS ITSELF and its floor.** Per ordinary
territory: 54.2M anchor-correspondence calls (1.11 per (unit, fold)
pair — the per-orientation pin build), 34.7M placement calls, 34.7M
assertion calls, at 0.31/0.69/0.95 µs. WHAT THE LIKELIHOOD ACTUALLY
NEEDS: one (matches, mismatches) per (unit, fold, orientation) —
the read's ~150 bases compared against the fold's own sequence at
the anchored coordinates. The IRREDUCIBLE-FOR-THE-MODEL floor is
that per-base comparison: ~150 bases × ~1-2 ns ≈ **0.15-0.3 µs per
pinned pair** (the no-pin floor's 0.18 µs/pair receipts the
per-pair loop overhead). The measured placement (0.69 µs) is 2-4x
over that floor (the `serving_anchor` searches and the per-pair
skipped-position recomputation — ENGINEERING); the ANCHOR class
(0.31 µs × every pair) is per-(unit, fold, orientation)
re-derivation of what the projection's occurrence structure already
holds — the read's node path is known, the fold's node positions are
precomputed (`positions_of`), so the correspondence is a COORDINATE
LOOKUP class, not an alignment (THE RE-ANCHOR PRINCIPLE —
ENGINEERING); and beyond the per-base floor, the only genuinely
per-fold work is the comparison of the bases where the fold's
sequence DIFFERS from its neighbors — the variant-pocket
delta-comparison form (named here as the candidate mechanism; its
per-locus savings are UN-RECEIPTED — a fold-sequence divergence
count is the one instrument this autopsy did not spend).

**(4) THE HONEST FLOOR (the budget table; arithmetic projections
from the receipts above, every row labeled; NO tuning constants).**
Per-locus core-seconds at L1-class (25 loci), the call ex-load
(panel load amortized at the fleet's fork form; the marathon tail
and the cut are separate terms below):

| row | loci phase | call ex-load | per-locus | gap to the 2.10 line | label |
|---|---:|---:|---:|---:|---|
| today (serial true CPU) | 164.8s | 202.7 | 8.11 | **3.9x** | — |
| under B (assertion receipt-only) | 131.9s | 169.8 | 6.79 | 3.2x | **OWNER-DECISION** |
| under census-restriction | 86.4s | 124.3 | 4.97 | 2.4x | **OWNER-DECISION** |
| under both | 69.2s | 107.1 | 4.28 | **2.0x** | (both owner decisions) |
| the 20-min line (52.6/25) | — | 52.6 | 2.10 | — | — |

(census-restriction = the measured 0.524 loci-phase ratio, applied
to the classes the shadow receipts as mass-scaling; B = the measured
32.9 assertion core-s; both are arithmetic projections from THIS
autopsy's receipts, labeled as such. The 202.7 charges the serial
warm-cache inputs residual (~14s); the fork-child rung's receipted
inputs class is 41-50s (cold sidecars), which puts today at ~215-224
core-s / ~8.6-9.0 per locus — the row is the lower receipt, the
spread is receipted, and the class is cacheable IO either way.) **THE RESIDUAL'S OWN
DECOMPOSITION under both** (4.28 core-s/locus, vs the 2.10 line):
the class scan 33.9 core-s = **1.36/locus — OWNER-DECISION** (the
exhaustive folds²/2 × units accumulation IS the committed answer
profile's byte-identity invariant; the algebraic per-fold
decomposition that would remove it changes the float accumulation
and is named and refused, §10.6.8); the pair loop 26.7 core-s =
**1.07/locus — IRREDUCIBLE-FOR-THE-MODEL core ~0.4-0.6** (the
per-base comparison at the 0.15-0.3 µs/pair floor over the pinned
pairs) **+ ENGINEERING** (the re-anchorable anchor class 8.7 core-s
and the 2-4x placement-walk residue over the raw comparison floor);
the serial bind+rows+IO classes 37.9 core-s = **1.52/locus —
ENGINEERING** (rows/context 14.3s serial; the inputs residual
14-50s, cacheable IO; bind 6.6s incl. the seek-floor residue, with
the fleet's pool-warm waiting deleted by the named candidate); the
fold-build/qual/exact classes 8.5 core-s = 0.34/locus, mixed.

**THE WALL-CHARGED RECONCILIATION (the correction this autopsy
forces on §10.6.13's headline).** §10.6.13 charged the call at
5,000-7,300 core-s/territory (phase wall × configured threads) and
bracketed "phase-timer CPU ~3,200-4,700 true"; the direct serial
measurement says the call's true CPU is **~203 core-s/territory
ex-load at L1-class** — the wall-charged convention overcounted the
call ~25-36x, exactly the idle-thread class the §10.6.13 T-rung
table receipted (per-core-second throughput 4.8x better at T8 than
T64: the fleet's call phases run ~1-4 cores of true work against
128+ configured threads — the small loci's T16 walls were ~5.6-6.0s
each vs 0.5-1.6s serial: queue waits, not work). THE HONEST
GENOME TERMS AT TRUE CPU, each with its receipt and its open term:
the CALL at L1-class = 5,844 × 202.7 / 256 ≈ **~1.3h** (after B +
census: 5,844 × 107.1 / 256 ≈ ~0.7h), plus the MARATHON TAIL
(§10.6.3's committed bound: 15-25 members per genome at the p90
2,912s to deepest 11,084s serial lane walls ≈ 44k-277k core-s ≈
~0.05-0.3h of box time at 256 cores — and those walls carry the
OLD bind era, so the marathon's true-CPU mass is the re-receipt's
owed number; the kill class is collapsed by C) — and THE CUT, whose TRUE CPU remains the
honest open term (wall-charged upper 1,936 at width 32; the scan's
gunzip class ~800-1,500 core-s/wave is real and width-amortized,
gzip-format-bound per §10.6.13's decode receipt; the remap/emission
classes' true CPU un-receipted). Whether true CPU converts to wall
is the fleet-discipline question the T-rung table already governs
(children × T ≤ 256). **The 130-175x remaining gap of §10.6.13 is
therefore restated: at TRUE CPU the call-core's own gap is ~3.9x
today and ~2.0x after the owner's two decisions; the bigger honest
unknowns are the cut's true CPU and the marathon wall, not the
call core.**

**THE EXACT RESUME POINT (if the budget ends here):** the census
instrument is committed (env-gated, default OFF, byte-identity
gated twice at the committed receipt bar; cargo test --release
--lib 562/0, 5 ignored); the receipts are in the work dir
(`global-call-L1-autopsy{,-control}`, `autopsy-*.txt`,
`run-locus-census-L1*.sh`). THE OWNER'S DECISIONS THE RECEIPTS NOW
PRICE EXACTLY: (1) candidate B — the always-on assertion, 58% of
the pinned alignment work, 16% of the call, receipted per-locus;
(2) the census-superset — a ~1.6x call-core multiplier at
L1-class, biting per-locality exactly where patterns don't
saturate (the per-locus table above). THE ENGINEERING ORDER THE
AUTOPSY NAMES: the pool-warm form (deletes the fleet bind walls),
the rows/context + inputs serial classes, the re-anchor (coordinate
lookup for the anchor class), the placement-walk residue over the
raw comparison floor, the variant-pocket delta comparison (needs
its own divergence-count receipt). STILL OWED: the cut's true-CPU
receipt (the last un-receipted term of the genome floor), the
marathon wall re-receipt at the next window, the W-x-T fleet cell.
NO PR PUSH (the owner pushes).

### 10.6.15 THE OWNER'S RULINGS, EXECUTED — the exactness assertion
out of the production path (Ruling 1) and the class-scan
accumulation-order toggle with the sanctioned sparse-skip fast mode
(Ruling 3); the census-profile A/B measured as the owner's decision
input (Ruling 2, UNDECIDED — §10.6.16's pointer)

**RULING 1 — THE ALWAYS-ON EXACTNESS ASSERTION IS OUT OF THE
PRODUCTION PATH (receipt/gate/test verification only).** The
always-on `direct_score` re-verification (§10.6.14's candidate B:
0.95 µs/pair, 32.9 core-s at L1 = 16% of the ex-load call, 58% of
the pinned alignment work) is now env-gated:
`IMPG_SCORE_EXACTNESS_ASSERT=1` restores the committed always-on
form VERBATIM (every placement re-verified base-by-base, the
per-locus and whole-run FACTORIZATION gates enforced); the
production default (unset) trusts the factorized (m, c) exactly as
the committed receipts' assertion proved it. THE GATES (the
committed serial recipe, the committed receipts-of-record, all six
answer files — calls.jsonl + sequences + exactness + ingredients +
records + skeleton): **L1-L4 under the production default ALL
BYTE-IDENTICAL** (work dirs `rulings-L{1..4}-prod` vs
`global-call-L1-fix` / `global-call-L{2,3,4}-fix-n4`) — the
assertion emitted no output bit, the gate is the proof of the
removal; the same gate PASSES for the **census form** (after
timers, `rulings-L1-census`), the **gate form** (`rulings-L1-assert`:
IMPG_SCORE_EXACTNESS_ASSERT=1 — the whole-domain proof still
enforces, 34,697,264/34,697,264 placements exact at L1, the exact
§10.6.14 count), and the **compare form** (`rulings-L1-compare`,
Ruling 3's delta receipt — exact downstream, byte-identical). THE
MEASURED RETURN (serial recipes; the eight gate runs co-tenant ×8,
walls 324-481s, the per-phase receipts are the honest numbers):

| lane | loci phase before (committed receipt, assertion ON) | after (production, OFF) | returned |
|---|---:|---:|---:|
| L1 | 162.3s | 128.0s | **−34.3s (−21%)** |
| L2 | 47.8s | 39.2s | −8.6s (−18%) |
| L3 | 150.6s | 115.4s | −35.2s (−23%) |
| L4 | 119.0s | 96.3s | −22.7s (−19%) |

The per-class census (timer-on walls, before = the §10.6.14 autopsy
census run, after = `rulings-L1-census`): anchor 18.6→17.5s,
placement 26.5→22.7s, **assertion 36.4→0.0s (0 calls)**,
marginal 3.3→3.2s — the assertion's whole measured mass is gone and
every other class is stable in the co-tenant band. THE RETURNED
PER-LOCUS BUDGET at L1-class: loci 128.0/25 = **5.12**; the
ex-load call ≈ 166.4 core-s = **6.66/locus** (the §10.6.14
projection "under B" was 6.79 — measured slightly better, inside
the co-tenant band), vs the autopsy's 8.11 today → **the executed
Ruling 1 alone is −18% of the ex-load call** (the projected −16%,
receipted).

**RULING 3 — THE ACCUMULATION-ORDER TOGGLE + THE SANCTIONED
SPARSE-SKIP (the previously REFUSED form, §10.6.8, now the owner's
performance path).** `IMPG_CLASS_SCAN_ORDER` selects the class
scan's accumulation: `exact` (the default, unset — the committed
deterministic unit-order accumulation, bit-reproducible, the
receipts' form) or `sparse-skip` (the fast mode — the §10.6.8
"named and refused" form: the constant entries, the (fold, unit)
pairs where BOTH rows hold the no-pin constant — where
mix_logsumexp(E, E) = E exactly — leave the per-unit loop and
return as ONE `elsewhere × constant_mass` term; the float
accumulation order CHANGES by construction). `IMPG_CLASS_SCAN_COMPARE=1`
(receipt-side, default OFF) sums BOTH ways per class pair and uses
the EXACT value downstream — a delta receipt that changes no output
byte. THE FAST MODE IS LABELED: the model string carries a
`-classscan-sparse` suffix in every emitted record (calls.jsonl,
sequences, the per-locus receipts) and in manifest.json — a
fast-mode output can never be mistaken for an exact-mode one.
THE GATES: (a) **exact mode byte-identity** — the default IS the
committed accumulation; the L1-L4 production gate fleet above
(byte-identical, all six files) and the compare run
(`rulings-L1-compare`, exact downstream, byte-identical) prove the
toggle's default side; (b) **the fast-mode float deltas, measured
and named per locus at L1** (`rulings-L1-compare`'s per-locus
`class_scan` receipt member): **336,879 of 350,365 class-pair sums
differ bitwise; the max ulp delta is 20,858** (the constant-mass
regrouping at the sums' own magnitudes — a ~1e-11 relative class);
**the ANSWER level: ZERO fold-material flips** (the called diplotype
is identical at all 25 loci; the sequences' fold records
identical), `called_classes` identical at every locus, and the
QUAL-visible deltas are **17/25 loci, max |ΔQUAL| = 2.7e-6** —
`best_log_likelihood` differs at 23/25 (the sums themselves);
(c) **the speed receipt at L1** (`rulings-L1-fast` vs
`rulings-L1-prod`, same co-tenant fleet): the class scan **71.6s →
53.0s (−26%)**, the whole loci phase 128.0s → 107.9s (−16%). THE
HONEST CAVEAT: L1 is a pinned-heavy territory (71% of pairs
pinned), so the constant share the skip removes is 1.098e9 of
3.711e9 unit-pairs (**29.6%**) and the union walk pays a merge per
kept entry — the −26% is the real L1-class gain; at the
marathon-class no-pin-heavy loci (1.3-2.1% pinned) the skippable
share is ~96-99%, but those substrates are reaped, so no marathon
speed number is claimed (the receipted rate is the L1-class one,
and the sparse-skip is exactly the form the compressed rows were
built for there).

**THE PER-LOCUS BUDGET, RESTATED UNDER THE TWO EXECUTED RULINGS**
(L1-class, 25 loci, arithmetic from this increment's receipts,
every row labeled):

| row | loci phase | call ex-load | per-locus | gap to the 2.10 line |
|---|---:|---:|---:|---:|
| §10.6.14 today (assertion on, exact scan) | 164.8s | 202.7 | 8.11 | 3.9x |
| **+ Ruling 1 executed (measured)** | 128.0s | 166.4 | 6.66 | 3.2x |
| **+ Ruling 1 + the fast scan at L1 (measured)** | 107.9s | ~146.3 | ~5.85 | ~2.8x |
| + census-restriction (§10.6.14's 0.524 ratio; SUPERSEDED — Ruling 2 executed, MEASURED at 3.96 exact / 3.95 stacked, §10.6.17) | ~56.5s | ~75 | ~3.0 | ~1.4x |

(the fast-scan row charges the L1-class skip share honestly; the
census-restriction row remains the §10.6.14 arithmetic projection,
now on the post-Ruling-1 base, labeled as such — the owner's
Ruling 2 is NOT decided, and §10.6.16's A/B table is the decision
input.) THE RESIDUAL after both executed rulings at L1-class: the
pair loop 42.9 core-s (anchor 17.5 + placement 22.7 + marginal
3.2, timer-on — the re-anchor and placement-walk ENGINEERING classes
of §10.6.14(d) unchanged), the class scan 53.0s under the fast mode
(the exact form's 71.6s remains the default), the serial
bind+rows+IO classes 37.9 (ENGINEERING, unchanged), the
fold-build/qual/exact classes ~12. cargo test --release --lib
**563/0** (5 ignored; the new
`class_scan_sparse_skip_union_walk` proves the union walk's skipped
count against the brute-force complement, the full-union pairs'
bit-identity, and the single-term all-no-pin form). The work-dir
artifacts: `run-rulings-gates.sh` + `rulings-*` (the eight-run
gate fleet), `run-rulings-ab.sh` (Ruling 2's A/B, below).

**RULING 2 — THE CENSUS-PROFILE A/B, MEASURED AS THE OWNER'S
DECISION INPUT (NOT A DECISION).** The two profiles ARE: (A) the
committed-slice restricted read set — the global cut restricted to
the committed 1x kept-reads via IMPG_CUT_READSET_FASTA (the §10.5.4
shadow form; census_mult 141,114 vs 384,056 at L1, and 3.6x/10.6x/
15.4x at L2/L3/L4 — the pilot's two profiles nearly coincide,
receipted in the table); (B) the global-native full mass (the
committed global cut census/cache — the receipts-of-record's own).
Both assessed with the OFFICIAL QV machinery (--truth-qv,
IMPG_TRUTH_QV_ALLWAVE=1, the likegt max-qv convention, the QV-50
cap) over all five localities: **the table is
`docs/syng-gmem-bwt/census-profile-ab.md` (§10.6.16's form)** — per
locality, per profile: rank-1 rate, the QV distribution (median,
≥QV30, at-cap), the call-core cost (the §10.6.14 census-superset
multiplier), and the NAMED loci where the profiles differ (which
fold flips, which profile was right per truth). The owner decides.

### 10.6.16 THE CENSUS-PROFILE A/B RECEIPT (Ruling 2 — the owner's
decision INPUT, measured; the decision itself UNMADE at this
increment — DECIDED AND EXECUTED at §10.6.17)

The five localities (L1-L4 + the pilot) ran under BOTH census
profiles with the official QV machinery (--truth-qv,
IMPG_TRUTH_QV_ALLWAVE=1, the likegt max-qv convention, the QV-50
cap) — profile A the committed-slice restricted read set (the
IMPG_CUT_READSET_FASTA shadow cuts, census_mult 141,114 / 48,641 /
67,359 / 18,121 / 512,929 at L1/L2/L3/L4/pilot), profile B the
global-native full mass (the committed global cuts, 384,056 /
156,976 / 567,243 / 275,405 / 513,836). THE HEADLINE: **a wash on
accuracy, not on cost** — rank-1 rates IDENTICAL at L1/L3/L4/pilot
(16%/8%/12%/18% both profiles) and B buys exactly ONE truth locus
at L2 (24%→28%, locus 17: B calls the truth pair, A does not, QV
50 vs 42); the median official QV within 1.5 everywhere (A higher
at L1 42.6 vs 41.1 and L4 45.0 vs 43.5; equal at L2/L3/pilot); of
the QV-visible fold flips A is closer to truth at 9 loci and B at
3; the rest of the 7-9 fold-material differences per locality are
sequence-equivalent copies (equal QV). THE COST, measured per
locality as the loci-phase ratio B/A: **1.77x (L1), 2.49x (L2),
5.02x (L3), 7.68x (L4), 1.008x (pilot)** — the pilot's two profiles
coincide (0.2% mass apart: identical answers, the A/B's own
control row), and the whole-call ex-load ratios are 1.57x / 1.75x /
3.07x / 3.67x / 1.008x — the §10.6.14 "~1.6x at L1" confirmed and
WORSE away from pattern saturation (the full-mass censuses are
8.4x/15.2x the slice at L3/L4). **THE FULL TABLE, the named
differing loci, and the machinery statement are
`docs/syng-gmem-bwt/census-profile-ab.md` — the owner's decision
input; this increment does not decide.** The runs: `rulings-ab-
{A,B}-<L>` (exit 0 ×10, serial recipes, the executed-Ruling-1
production path, co-tenant ×10), the analysis
`rulings-ab-analysis.txt`, the shadow cuts `global-cut-
{L2,L3,L4,pilot}-shadow1x`.

### 10.6.17 THE OWNER'S RULING 2, EXECUTED — RESTRICT IS THE PRODUCT DEFAULT; the owner's routing model verified against the machinery; the L2-17 exception autopsied (the conserved-elsewhere dictionary class); the gates; the per-locus budget restated under the stacked rulings

**THE RULING + THE DESIGN PRINCIPLE (the owner's own words, now the
stated principle of the census):** RESTRICT is the product default —
the global-frame cut's default census profile is the committed routed
read set (the committed 1x kept-reads, the evidence mass the original
local-frame receipts banked), NOT the global-native full mass; the
extra 2.8x-15.2x mass is edge-grazers + repeat-routed reads —
non-discriminating, it multiplies work without moving the argmax.
*"we just need to just worry about the overlapping ends — the MEM
matches at the ends — and we used to decide where to put each read,
so it's basically being routed by its best MEMs. If they go into
different bins unless they're basically adjacent, we don't want to be
thinking about them."* The one exception the A/B measured (L2 locus
17: the single clean global-native win — B calls the truth pair at
QV 50, the restricted profile calls [0,0] at 42) was autopsied to the
READ before finalizing; the verdict is below, and it is NOT the
owner's boundary case.

**THE ROUTING-MODEL VERIFICATION (the owner's model vs the actual
machinery, code receipts).** The machine IS the model's boundary half,
exactly: a read record routes by its BEST MEM — `route_record`
(examples/sample_projection.rs) seeds placements from the anchor with
the FEWEST territory-index entries (the most discriminative k-mer),
extends each placement through the verified walk
(`occurrences_from_anchor`), and the occurrence's touched bins are
every window whose interval overlaps any anchor's k-mer
(`partitions_at`: `bp + k > start && bp < end`) — a
boundary-straddling read lands in BOTH adjacent windows by pure
interval overlap: no tolerance, no snapping, nothing to tune (the
window grid is contiguous — L2's 25 windows have zero gaps, verified).
The call side keeps the convention (`locality_touch`, realign_score:
an occurrence votes at a locus iff `occ.partitions.contains(&locus)`,
the committed partition mode), and the census's own pair-bins taxonomy
(`same_window`/`adjacent_window`/`same_contig_far`/`other_contig`)
already names the owner's bins. `SortedRanges::insert_with_boundary_snap`
(src/impg.rs, consumed by multi_impg) is NOT read routing — it is the
territory-expansion range union (expansion ranges snap onto adjacent
ranges within a caller-stated min_distance to absorb alignment
boundary jitter; user-requested ranges never snap). THE HONEST
ADDITION: the model's other half — "if they go into different bins
unless they're basically adjacent, we don't want to be thinking about
them" — is NOT implemented as a drop: the machine keeps a
non-adjacent multi-window record at EVERY far bin it touches. The
measured mass share of that class is real at some territories: L1's
restricted census carries **31.1% of its mass in 14 non-adjacent
multi-window records — ONE single-anchor pattern of multiplicity
43,686 routed to windows [12,16]**; L3 21.7% (56 records); L4 0.02%;
L2 0.0% (every multi-window record is adjacent-only); the pilot 2.6%
non-adjacent + 14.0% adjacent. The owner's rule is therefore a REAL,
un-implemented census lever — NAMED here, not built (the ruled-in
restrict profile is the committed routed read set, which contains
these records; an adjacency-filtered profile would be a separate
owner decision). And the mass the ruling actually removes — the B−A
difference — is NOT far-routed reads: at L1, 97% of the
difference mass (274,719 of 283,030) is SINGLE-window single-anchor
repeat mass, dropped by restrict through the local expressibility
basis (the autopsy's class, next).

**THE L2-17 AUTOPSY (the one exception, traced to the read — the
deciding evidence).** The truth variant is one extra (AAAT) unit (a
4bp insertion) in the array at window-local position ~1307 of window
17 ([170035,180062) territory-local) of HG002#1#chr13:28850052-29100914
— MID-window (1.3kb from the left edge, 8.7kb from the right) and
mid-territory: no boundary is involved. In the full-mass census the
decisive record is census record 193 of global-cut-L2-d1: a
SINGLE-anchor pattern (node −6756230, multiplicity 13, 178 panel
occurrences ALL inside window 17 — a single-locus, per-strain
orthologous repeat anchor), whose representative read (cutread 43820)
is an EXACT 148bp match of the truth haplotype HG002#2 across the
junction; its fold6-vs-fold0 log-likelihood gap is **+4.005 — the
largest pro-truth term at the locus** (receipt: the
rulings-ab-B-L2 records receipt; the next three records are
ANTI-fold6 at ≈−3.95 each). The second truth-allele junction read
(43853: junction-exact, one sequencing error elsewhere) shares the
record. THE VERDICT: the lost read is **(c) SOMETHING ELSE — the
conserved-elsewhere dictionary class** (§6 risk #1's own named second
sub-class: "a conserved-elsewhere k-mer that the global dictionary
resolves but the local one lacks"), NEITHER (a) a boundary-straddling
read (its only MEM routes entirely into ONE window — it is already
"kept at its bin") NOR (b) a repeat-routed far read (zero occurrences
outside the locus). The mechanism, receipted end to end: the
committed LOCAL cut dropped the record because its anchor k-mer was
never selected within the locality's own content — the locality
dictionary lacks it while the global dictionary resolves it
(membership earned by selection elsewhere in the panel); the local
node map (12,284 pairs vs the global identity map's 12,548) therefore
lacks the pair, the strict all-or-nothing remap dropped the record,
and the reads whose ONLY record it was died with it — including both
truth-allele junction reads. Reads with other records survived (the
56 array-core reads of the committed set all kept, 55 via their other
records), which is why the restricted census still scores the locus —
but with zero insertion-carrying evidence (its own single-anchor
array record, pattern (192229540,), multiplicity 19, has a
fold-0-favoring representative: gap −3.953) — and calls [0,0]. THE
CODE RECEIPT (examples/census_autopsy_probe.rs, outputs in the work
dir census-autopsy-L2-17/): the canonical anchor walk (the
expressibility map's own rule) over the window-17 slices against BOTH
panels — the GLOBAL walk has the anchors the locality walk lacks, at
the exact positions: 1237 (node −48057385, both haplotype phases) and
1255 (node −6756230, the truth-allele phase); the locality dictionary
does not contain either k-mer (its own array members are the
shifted-phase nodes 11263-11266). **No adjacency rule recovers
L2-17 under restrict — the cost is intrinsic to the restriction
basis, and it is the measured price of the ruling: exactly one truth
locus.** The recovery path, NAMED NOT IMPLEMENTED (a separate gated
increment if the owner ever orders it): a restricted read set whose
expressibility basis is the locality map EXTENDED by the global
dictionary's conserved-elsewhere class — a dictionary-level
widening, decided by its own accuracy/cost receipts, never by a
tuned threshold.

**RESTRICT AS DEFAULT (the implementation).** The global-frame cut
gains the three-step interface's own product flag:
`--census-profile <restricted|global-native>`, DEFAULT RESTRICTED
(the codebase's product doctrine: product behavior is a visible CLI
flag; env vars are test-mode machinery and non-default restores).
The restricted profile keeps only the reads of the territory's
committed kept-read set, stated by the IMPG_CUT_READSET_FASTA
instrument (the restriction basis override, untouched — a stated
basis restricts in either profile, loudly); the restricted profile
WITHOUT a basis fails closed before the panel load, naming both
forms; the global-native profile is the full routed mass,
byte-identical to the committed global-frame cut (the record form,
behind the explicit flag). The LOCAL frame is untouched
byte-for-byte (the flag is refused there; the read-set instrument
still restricts as before). The manifest names the effective profile:
cut-stats.json gains `census_profile` (global frame only — the
local frame keeps its committed stats bytes; the batch forms' stats
name their actual full-mass profile, and carrying the profile into
the batch manifests the way record_origin rides them is the named
follow-up). The resolution is a pure tested function
(`CensusProfile::resolve`): five states receipted — local untouched,
global default fail-closed, basis-required restricted, global-native
record form, basis override.

**THE GATES (all at the byte-identity bar).**
(a) **The default-restrict cuts, L1-L4 + the pilot: ALL PASS** —
census/cache/reads BYTE-IDENTICAL to the A/B's profile-A shadow cuts
(census_mult 141,114 / 48,641 / 67,359 / 18,121 / 512,929 reproduced
exactly; cut-stats.json equal modulo the new census_profile field and
the timing walls; walls 623/606/968/960/556s, co-tenant ×2 waves).
Because the call is a pure function of the census/cache/reads, the
A/B's profile-A calls are the restrict-default ANSWER receipts.
(b) **The end-to-end default form at L2: PASS** — the committed
serial call over the default-restrict L2 cut is byte-identical to
the A/B profile-A receipt (rulings-ab-A-L2; wall 409s), and **locus
17 stays [0,0] — the receipted one-locus cost**, per the autopsy: not
a routing gap, so the fixed form does not and cannot recover it.
(c) **The record form: PASS** — `--census-profile global-native`
reproduces the committed global cut byte-identically (L4 vs
global-cut-L4-d1: census/cache/reads cmp-clean; cut-stats equal
modulo census_profile + timing; census_mult 275,405 reproduced).
(d) **The budget run (restrict + Ruling 1 + the sanctioned
sparse-skip, L1): PASS** — the fast mode over the restricted census:
loci phase 70.7s, ZERO fold-material flips (25/25 identical called
folds), identical official QVs, the -classscan-sparse label present
(wall 415s serial, co-tenant).

**THE PER-LOCUS BUDGET, RESTATED UNDER THE STACKED RULINGS (L1-class,
25 loci, every row a receipt; the ex-load convention of §10.6.14 —
the serial panel load excluded):**

| row | loci phase | ex-load call | per-locus | gap to the 2.10 line |
|---|---:|---:|---:|---:|
| §10.6.14 today (assertion on, exact scan, full mass) | 164.8s | 202.7 | 8.11 | 3.9x |
| + Ruling 1 executed (measured) | 128.0s | 166.4 | 6.66 | 3.2x |
| + Ruling 1 + the fast scan, full mass (measured) | 107.9s | ~146.3 | ~5.85 | ~2.8x |
| **+ Ruling 2: restrict-default, exact scan (MEASURED, this increment)** | **75.0s** | **99.0** | **3.96** | **1.9x** |
| **+ all three rulings stacked (restrict + R1 + fast scan; MEASURED)** | **70.7s** | **98.7** | **3.95** | **1.9x** |

The §10.6.15 projection row ("~3.0 under census-restriction") scaled
the whole ex-load by the 0.524 mass ratio and UNDERSTATED the
mass-invariant classes: the measured 3.96 shows the bind+rows+IO
classes (~1.1/locus) do not scale with census mass — exactly the
§10.5.4 finding, now at the per-locus grain. The sparse-skip's
honest receipt at the RESTRICT form: loci 75.0 → 70.7s (−5.7%),
far under the full-mass −16% — the restricted census is pinned-heavier
(fewer constant-mass pairs to skip); the marathon-class skip share
remains unmeasured (reaped substrates), nothing claimed there.

The per-locality measured shape at the restrict-default form (exact
scan, post-Ruling-1 production path, serial co-tenant runs; walls
include the ~260-290s panel load):

| locality | loci phase | ex-load | per-locus | wall |
|---|---:|---:|---:|---:|
| L1 (25 loci) | 75.0s | 99.0 | 3.96 | 393s |
| L2 (25) | 16.9s | 34.5 | 1.38 | 301s |
| L3 (24) | 24.2s | 47.0 | 1.96 | 327s |
| L4 (25) | 13.6s | 34.4 | 1.38 | 306s |
| pilot (22; the degenerate row — both profiles coincide) | 630.9s | 668.4 | 30.4 | 998s |
| L1 + the fast scan (the budget run) | 70.7s | 98.7 | 3.95 | 415s |

**THE RESIDUAL POLES at the restrict-default form (L1-class, each
labeled):** (a) the loci core 2.83/locus — the re-anchorable anchor
class, the placement walk, and the class scan under the fast mode (the
§10.6.14(d) ENGINEERING classes, unchanged by the census ruling);
(b) the mass-invariant serial classes 1.12/locus — bind 7.9s +
rows/context 16.6s + quality 3.4s (ENGINEERING; the §10.5.4
mass-invariance finding — the census ruling does NOT shrink these,
which is why the measured 3.96 sits above the old ~3.0 projection);
(c) THE FAR-BIN PRUNE LEVER — the owner's own model, UN-IMPLEMENTED:
31.1% of L1's restricted census mass (the one 43,686-mult pattern at
windows [12,16]), 21.7% at L3 — an adjacency-filtered census profile
would remove it and is the owner's call, not built here;
(d) the cut's true-CPU census — still OWED (the last un-receipted term
of the genome floor); (e) the whole-call 20-minute line at the
measured shape: the ex-load gap to the 2.10/locus line is now
**1.9x at L1-class** — the closest the receipts have ever stood to
it — with the remaining distance owned by the loci core's
engineering classes and the far-bin lever, not by the census mass.

The work-dir receipts: `census-autopsy-L2-17/` (the probe outputs,
both panels × both haplotype phases), `census-default-{L1,L2,L3,L4,
pilot}/` + `census-native-L4/` (the gate cuts), `census-default-call-
L2/` (the end-to-end default-form gate), `census-default-budget-L1-
fast/` (the stacked-rulings budget run), `run-census-default-gates.sh`
+ `run-default-call-L2.sh` + `run-budget-L1-fast.sh` (the recipes).
NO PR PUSH (the owner pushes).

### 10.6.18 THE FAR-BIN PRUNE, ENFORCED — the owner's routing principle at the call's census view; the anchor-class hoist; the stacked budget at ~1.2x the 2.10 line; the cut's true-CPU census IN FLIGHT (the resume point)

**THE PRUNE, EXECUTED.** The owner's principle — *"if they go into
different bins unless they're basically adjacent, we don't want to
be thinking about them"* — is now an enforced call-mode flag:
`--far-bin-prune` (requires `--global-frame`; the committed local
frame stays byte-frozen). The rule is pure window geometry, nothing
tuned: a census record's own routing keeps it at the windows its own
walk placements touch — ONE window, or TWO ADJACENT windows when
the anchors straddle a boundary (both by the routing's own interval
overlap; boundary-straddlers keep BOTH bins). A pattern whose
placements scatter into NON-ADJACENT windows — the measured far-bin
class, e.g. L1's one 43,686-multiplicity single-anchor pattern at
windows [12,16] — is the "different bins, not basically adjacent"
class: there is NO choice-free home bin to keep it at, so it is
pruned from the call's census view ENTIRELY (pilot selection, bind,
variants, and every locus's scoring). The test is the touched-window
DIAMETER against the census's own pair-bins adjacency constant 1
(`census_record_far_bin_scattered`; unit-proven). The cut files are
untouched (the census receipt keeps its honest full routing + dense
ids); the run labels itself (`-farbin-prune` on the model string)
and receipts its own mass drop to stderr. Default OFF: the committed
path byte-identical.

**THE GATES (two-sided, all PASS).** (a) BYTE-IDENTITY at every
locus the prune does not touch, across the committed gate-file set
(calls/sequences/exactness/ingredients/records/skeleton + the
instrument receipt modulo its timing fields): **L2 — 0 records
pruned — ALL 25 loci byte-identical**; at the other localities every
locus OUTSIDE the pruned records' own windows is byte-identical (L1
19/25, L3 20/24, L4 22/25, pilot 14/22), and the changed loci are
EXACTLY the pruned records' windows: L1 [2,4,11,12,13,16], L3
[1,13,17,19], L4 [9,17,23], pilot [1,4,7,8,9,10,11,12]. (b) THE
TRUTH ASSESSMENT at the changed loci (official QV machinery, the
census-profile-ab form): **ZERO truth losses** — rank-1 rates
IDENTICAL at L1/L2/L3/L4 (16%/24%/8%/12%), median QVs and >=QV30
counts identical — and at the PILOT the prune WINS FOUR TRUTH LOCI:
rank-1 4/22 (18%) -> **8/22 (36%)**, median official QV 36.3 ->
**45.0**, >=QV30 17 -> 20 (loci 8/9/10/11: the pruned form calls
the truth pair at QV 50 where the restrict-default called
non-truth fold material at QV 0.8-27.3; locus 12 equal-class QV
shift only). The named fold-material changes with NO truth stake
(neither call matches truth, equal or better QV): L1 loci 4, 11
(locus 11's QV RISES 30.4 -> 39.6); L3 loci 1, 13, 19. The owner's
model is not merely cost-free — at the pilot it is accuracy-
POSITIVE: the far-bin chain records ({8,9,10,11,12}-class repeat
patterns, 1,167 records, 10.7% of the pilot's census mass) were
drowning the local evidence at five consecutive windows.

**THE MASS/WALL RECEIPTS.** Census mult pruned per locality (the
run's own receipt line): **L1 14 records / 43,868 of 141,114
(31.1%); L2 0; L3 52 / 29,889 of 67,359 (44.4%); L4 3 / 11 of
18,121 (0.06%); pilot 1,167 / 54,846 of 512,929 (10.7%)**. Loci
phase (serial recipe, RAYON_NUM_THREADS=1, walls with co-tenancy
stated — the box carried a load average of ~180 and another
campaign all evening): L1 75.0s (committed receipt) -> **34.7s**
at the stacked form, measured as a CLEAN SAME-WAVE PAIR (unpruned
new binary 74.8s vs pruned new binary 34.7s, the pair's own wave) =
**−53.6% of the loci core**; the same-wave unpruned references are
byte-identical to the committed receipts (determinism re-gate
PASS at L1+L3); pilot 630.9 -> 288.3s (−54%); L3 24.2 -> 15.2s.

**THE RE-ANCHORABLE ANCHOR CLASS — profiled, verdict, implemented.**
What it computes: per (unit, fold, orientation) the committed form
rebuilt the pin skeleton — `anchor_correspondence`'s monotone DP
over the fold's node positions, the merged pinned-window backbone,
the windows-inside test — which is a PURE function of (record,
fold, orientation): the canonical walk is record-level state, the
node positions fold-level, and NO variant input enters it. The
§10.6.14 receipt: 16.8 core-s at L1 (54.2M calls at 0.31 µs; 26.3M
calls at the restricted census). THE VERDICT: the projection's
occurrence structure ALREADY serves this class — the node->position
lookup IS the k-independent panel coordinate the re-anchor
principle names (conversion-not-alignment); the waste was the
per-VARIANT recomputation of a per-RECORD result. Implemented as
the anchor-class hoist: the pins are built ONCE per (record, fold)
before the unit loop, every per-unit term (the integer
abstained-hull accumulation, the no-pin elsewhere branch, the
per-(unit, pin) placement) kept in the committed unit order with
the committed values — bit-identical by construction. THE GATE:
byte-identity PASS at L1 and L2 (the new binary WITHOUT the prune,
all artifacts, vs the committed receipts). THE RECEIPT (the
scoring-census instrument, serial, restricted census, both forms
raw): anchor class **11.63s / 26,330,418 calls -> 1.08s / 626,701
calls (−90.7%, a 42x call collapse)**; the placement class
unchanged (13.4-13.9s, identical call counts — the hoist touches
only the lookup class). The hoist's whole-loci-phase effect at the
ordinary localities sits inside the evening's contention noise
(74.8-94.5s across waves vs the receipt's 75.0s); the class census
is the solid number.

**THE STACKED BUDGET (L1-class, 25 loci, ex-load = the phase
summary's non-inputs sum, serial recipe, each row's wave stated):**

| row | loci phase | ex-load | per-locus | gap to 2.10 |
|---|---:|---:|---:|---:|
| restrict-default (the §10.6.17 receipt, ×10 co-tenant) | 75.0s | 99.0 | 3.96 | 1.9x |
| + the far-bin prune (pre-hoist binary, exact, ×5 wave) | 39.7s | 66.4 | 2.66 | 1.3x |
| + prune + the anchor hoist (exact; CLEAN PAIR wave) | 34.7s | 61.1 | **2.44** | **1.16x** |
| + prune + hoist + the sparse-skip fast scan | 34.4s | 66.2 | ~2.65* | ~1.3x |

(*the budget run's wave was the evening's most contended; its loci
phase 34.4s is within noise of the exact form's 34.7s CLEAN-PAIR
34.7s — at the pruned form the sparse-skip buys nothing measurable,
honestly receipted: the pruned L1 census is pinned-heavy and the
skip's −5.7% at the unpruned restrict form vanishes under it. The
mass-invariant serial classes — bind 6.5s + rows/context 17.0s +
quality 2.8s ≈ 1.06/locus — are unchanged and now own ~44% of the
budget.) **THE HONEST HEADLINE: the L1-class gap to the 2.10 line
is now ~1.16x** — the 20-minute line's call pole is, within the
contention noise of any single wave, essentially AT the line; the
remaining distance is the mass-invariant serial classes (ENGINEER-
ING) and the placement/class-scan core (the class scan's exact
accumulation is an OWNER byte-identity decision, the placement
class's per-pair comparison is the IRREDUCIBLE floor ~0.4-0.6
core-s/locus of §10.6.14).

**THE CUT'S TRUE-CPU CENSUS — IN FLIGHT (the owed duty-3 receipt).**
The serial-census method (the call autopsy's own: RAYON_NUM_THREADS=1,
walls = true CPU, /usr/bin/time user+sys beside) is running at the
width-amortization shape: ONE serial width-4 batch cut (one copy of
each gate territory) and ONE serial width-32 batch cut (8 copies
each; the genome-scale marginal = the flat decode share + ONE
copy's remap, the §10.6.13 honesty label), both full-mass scans.
AT COMMIT TIME both runs are ~45 minutes into the per-territory
setup phase (the [cut-batch-global] territory/setup lines; the
scan+emission phases follow). **THE RESUME POINT: (1) let
`serial-cut-w4` + `serial-cut-w32` finish (work dir
`run-cut-serial-census.sh`, receipts `serial-cut-{w4,w32}.{log,err,
wall,exit,time}`), decompose decode-vs-per-territory from the
phase lines, and state the cut pole's honest per-territory true CPU
in the 20-min ledger; (2) restate the whole-genome 20-minute
budget table with every pole receipted (the call pole now ~2.44/
locus at L1-class; the cut pole pending (1); the load/scan/ memory
poles unchanged from §10.6.13-14) and each residual labeled
ENGINEERING / OWNER-DECISION / IRREDUCIBLE; (3) if the owner wants
the prune as the product DEFAULT (it is gate-clean: zero truth
loss, four pilot truth gains), flip the flag's default by its own
gated increment — the flag default OFF keeps the committed
receipts byte-frozen until then.**

The work-dir receipts: `farbin-A-{L1,L2,L3,L4,pilot}/` (the pruned
calls), `farbin-ref-A-{L1,L3}/` + `stackref-L1/` (the same-wave
unpruned references, byte-identical to the committed receipts),
`hoistgate-{L1,L2}/` (the anchor-hoist byte-identity gates),
`prehoist-census-L1/` + `hoistgate-census-L1/` (the anchor class's
A/B census receipts, old form vs hoisted), `farbin-hoist-A-L1/` +
`farbin-budget-L1-fast/` (the stacked exact + sparse-skip runs),
`impg-prehoist` (the preserved pre-hoist binary the wave ran),
`run-farbin-gates.sh` + `farbin-compare.py` + `run-cut-serial-census.sh`
(the recipes). NO PR PUSH (the owner pushes).
