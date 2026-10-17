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
   setsid, ALIVE at 19:13Z at ~12.9k/38,790 lanes, ~55 Mbp/s serial,
   ETA past midnight under co-tenancy) — on completion: `cmp` its
   `graph.json` against `$WORK/global-routes-par/graph.json` (the
   free full-panel byte-identity receipt; the same HEAD compiler
   identity — expected byte-identical; if it differs, name it
   per-field and re-examine the gates' verdicts against the serial
   artifact).
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
