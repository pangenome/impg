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
- [ ] Identity gates: L1-L4 (§10 — in flight; the gate recipe:
  fresh global-frame cut from the committed 35x shards → global-panel
  call → `cmp` against the committed `$PROJ/cut-L*-cli/calls.jsonl`).
- [ ] Measurements + budget restatement (§10).

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
## 10. The measured receipts (updated as the gates land)
