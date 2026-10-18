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
attempted.

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

### 10.6.11 THE SERIAL FULL-PANEL CMP + THE EXACT RESUME POINT

The serial whole-panel derive-routes (PID 2121938) is ALIVE at this
increment's close (9h34m CPU, ~574 min; the watcher PID 2483558
performs the cmp the moment its graph.json lands and writes
serial-vs-par.routes-cmp-receipt.txt — never polled in-session). The
receipt is owed to the next increment's close per the watcher's
design.

**THE EXACT RESUME POINT (if the budget ends here):** all three
duties landed and gate-proven (3e9f9b8 the concurrent calls + the
per-call read length; c4b7e86 candidate C's row-compression +
streaming); cargo test --release --lib 561/0 (+2: the compression
value-identity proof, the marathon-shape memory receipt — the
latter #[ignore]d, run explicitly: 39.0x at T03803's shape); the
example proof 2/2; tree clean except vendor/syng (uncommitted by
design); NO PR PUSH. NEXT, in the §10.6.10 pole order: (1) the
serial-vs-par receipt at the watcher's completion; (2) THE PANEL
SHARE'S IMPLEMENTATION INCREMENT (the fork-COW probe FIRST — a
forked child's Pss growth over a full call, the §10.6.9 named risk
gate — then the prefork batch driver, full byte-identity gates at
every step); (3) candidate B (the always-on exactness assertion —
the owner's receipt-shape decision) and the census-superset answer
(the owner's); (4) the marathon re-receipt at the next campaign
window (determinism + wall + candidate C's bound); (5) the bind's
per-fetch AGC-seek floor (the span-volume collapse went 508x; the
seek COUNT is the named floor — a bulk-prefetch form is the
candidate, undesigned).
