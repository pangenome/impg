# THE EVIDENCE STORE, THE UNION SCAN, AND THE BATCHED CUT — receipts of the performance night

> **WHAT THIS IS.** Three changes, each gated on the pipeline's own
> artifacts, plus the measurements that say what the per-territory cost is
> made of. Everything below is measured; the projections are labelled as
> projections. The target the work is measured against is **likegt-class
> cost** — a caller that genotypes by a single dot product per candidate —
> not a wish for "minutes".

## 1. THE DIAGNOSIS THE NIGHT STARTED FROM

Per territory, measured, before any of this:

| stage | wall | note |
|---|---:|---|
| map leg (`infer map`, one 250 kb territory) | **695 s** | 96% of it ONE scan |
| call-region's cut | 317 s | its own scan of the same evidence |
| the model (24 loci) | 103 s | ~10% of the total |

and the map leg's scan, from its own receipt:

```
[cut] scan done: 734,443,930 reads scanned, 276,042 kept, 5,988 touching patterns
```

**The entire sample's reads decompressed per territory to keep 0.04%.** The
same shape on the call side. So 73% of the per-territory cost was two
full-sample scans, and the model — the part that looks expensive — was 10%.

## 2. THE STORE (what was missing)

The projection emits reads in ARRIVAL order, so nothing can be seeked. The
store makes **position order equal file order** — BAM's own trick — and a
sparse manifest is then all the index needed.

Two artifacts, two lifetimes:

| artifact | lifetime | size |
|---|---|---|
| `node -> (path, offset)` (`examples/panel_positions.rs`) | **once per panel, forever** | 4.1 GB, 271,000,143 nodes |
| the evidence store (`IMPG_EVIDENCE_STORE`) | **once per sample** | **22.4 GB** (168 GB raw, zstd per bucket) |

Two facts that decided the design, both measured:

1. **The syng node id space is NOT positionally ordered.** A median 3-node
   walk spans 4,697 ten-thousand-node buckets, and 63% of reads span >100.
   Keying by node id is keying by a hash; the key must be `(path, offset)`.
2. **A vertex's position is orientation-independent.** Keying on the SIGNED
   node id silently dropped every reverse-strand read: `+118354468 FOUND,
   -118354468 ABSENT` — 45% placement before the fix, 83% after, and then
   **0 walk-bearing reads unplaced** (the completeness property; the 17%
   that do not place are reads with no MEM records at all).

Exactness is structural: a read is stored in the bucket of EVERY walk
node's position, so "does a walk touch node N" is answered by reading N's
bucket — no window geometry, no rewind. (The bed-window query is the wrong
predicate and fails as a strict SUBSET: 488 of 5,678 records, zero extra.
That failure shape is what made the diagnosis immediate.)

## 3. THE GATES

**Store as the map leg's packet source** (`IMPG_EVIDENCE_STORE` +
`IMPG_PANEL_POSITIONS`), against the full-shard run:

```
cut-census.jsonl   BYTE-IDENTICAL (148,292,612 B)
cut-cache.bin      BYTE-IDENTICAL ( 54,141,088 B)
cut-reads.fasta    BYTE-IDENTICAL ( 45,435,820 B)
scan counters      identical: 276,042 kept, 5,988 touching patterns
fetch              145,441 nodes -> 8,543 buckets (30 unresolved),
                   14,949,866 packets, 668 MB read, 23 s
wall               371.30 s  vs  637-670 s   (1.75x)
```

**Batched cut** (`cut-batch-seq`: ONE panel load, the committed single cut
per territory, byte-identity by construction):

```
panel loaded once in 184.3 s; 3 territories
  T00048 155.0 s   T00049 67.6 s   T00055 44.3 s      WALL 463.29 s = 154 s/territory
GATE: T00048 and T00049 both BYTE-IDENTICAL on cut-census, cut-cache and
      cut-reads against their own standalone map legs
standalone comparison: T00049 264.51 s alone -> 67.6 s in the batch (3.9x)
```

The panel load was **51% of the 371 s** (190 s, paid per territory for the
same panel). `172 + N x ~120` instead of `N x 371`.

**The map leg: 637-670 s -> 154 s/territory, byte-exact.** Per sample at
16 lanes: **66 h -> ~21 h**.

## 4. THE UNION SCAN (honest version)

The ll matrix is 96.9% the elsewhere constant at the tail locus, so in
principle a class-pair sum needs ~2,700 corrections, not 45,423 — 11.1x
territory-wide IN UNIT-PAIR COUNT. `IMPG_CLASS_SCAN_UNION` wires that into
the bounded scan's real work sites (stage one) with certification: the
winner's float window (measured: ONE pair), the top 32, and the truth pair
are recomputed in committed order.

**Realized speedup: 2.1x on the class phase, not 11x** — because the
excess-mass bound keeps exactly the pairs with HIGH excess mass (the
pinned-heavy folds), so the two optimizations work against each other. And
the class phase is ~5 s of the call's ~103 s: the poles there are the
spectrum (16 s) and the certificate (14.8 s).

Accuracy: **decision-exact at 24/24 loci** — identical `diplotype`,
`called_classes`, `qual_cluster_k`, `qual_cluster_size`, `qual_knee_distance`,
`qual_spectrum_shape`, `qual_unbounded`. Only float digits move, bounded at
9e-9 relative (the gap fields amplify an ulp by cancellation). NOT
byte-identical, and never could be: the union walk reassociates the sum.

## 5. WHAT THE SEARCH QUESTION COST TO ANSWER (three negatives)

The tail locus owns **73% of a territory's compute** (locus 8: F=429,
U=45,423, 4.19e9 unit-pairs; the bottom half of the loci own 1%). So every
screen was tested against it:

| approach | recall | pairs saved |
|---|---:|---:|
| hierarchical centroid drill-down (k=4, top-1 descent) | **67%** | 115x |
| per-fold haploid screen | 48% @ m=16, 83% @ m=32 | 10-25x |
| per-fold likegt-style cosine screen | 13-17% | 10-25x |
| cheap per-pair max-dominance proxy | 22% top-1, 65% top-8 | none |

**The diplotype objective is interaction-dominated**: two individually
mediocre haplotypes make the best diplotype, so nothing that ranks folds
solo, clusters them, or approximates the per-pair sum recovers the winner.
That is why the model enumerates pairs, and it is why these three ideas
are closed with numbers rather than with opinions.

## 6. WHAT IS LEFT, MEASURED

| term | per territory | note |
|---|---:|---|
| map leg (now) | 154 s | 155/68/44 s marginal for three adjacent territories — the tail shape again |
| call-region cut | ~231 s per CHUNK | already amortized over 32 territories (~7 s/territory) — do NOT wire the store here |
| the model's poles | ~31 s | spectrum 16 s + certificate 14.8 s at the degenerate locus |
| the certificate's O(f^4) | 17.5 s at f=124 | 46.5M entries; 10^10 at f=450 (~40 min) — a plateau, the exact bound prunes only 36% there |
| times the panel is loaded | 2 per territory | the map leg now amortizes; the call leg already COW-amortizes per chunk |

## 7. TRAPS FOUND (all the author's, all caught by gating)

1. `payload_base = 40` copied from the single-file prototype — the
   whole-store payload is its own file at offset 0. The gate surfaced it in
   seconds as a zstd `Unknown frame descriptor`.
2. The positions table loaded into a Vec and RE-SORTED per cut child —
   4.1 GB and a 271M-element sort on every territory, which would have
   eaten the entire win. Now mmapped.
3. **`impg infer map` execs a sibling binary** (`examples/sample_projection
   cut`): a library change compiles, runs in the parent, and then measures
   a STALE child doing the full scan for 11 minutes. "It built" does not
   mean "it ran". The CLI seam removes this class of trap.
4. The union-walk patch landed on stage two when the work is in stage one —
   the census line said `0 of 0 kept-pair sums`, which is what caught it.
5. An inverted compare branch (compare mode returned the FAST value) —
   caught because `calls.jsonl` diffed.
6. A CLI rename (`infer` -> `genotype-infer`) broke 16 integration
   invocations while `cargo test --lib` stayed green. The gate is
   `cargo test --release` with no filter.

**AND THE METHOD ONE:** the exactness question that cost four 190 s
panel-load runs was answerable OFFLINE in one minute (the store's fetched
read set vs the scan's own `cut-reads.fasta`: 239,124 kept reads, 0
missing). Artifact-vs-artifact for design; full runs only for the closing
gate.

NO PR PUSH (the owner pushes).

## 8. THE OBJECTIVE: THE LIKELIHOOD STAYS (the owner's ruling)

Asked whether the pipeline should move to likegt's objective (cosine over
node-coverage vectors — one dot product per candidate, cheap by
construction) the owner's answer was that **the likelihood is the stronger
object, and node-coverage cosine is a soft feature by comparison**
("the cosine distance of the node coverage ... is kind of a soft object").

So the expensive objective stays, and the cost work is to make IT cheaper —
which is what sections 2-4 above did (and what the interaction-dominated
evidence in section 5 explains: the pair structure is the signal, which is
also why no cheap screen recovered the winner).

Consequence for the cost target: likegt-class wall is NOT reachable while
keeping this objective. The reachable envelope is what the measurements say
it is — the map leg at 154 s/territory, the model's poles (spectrum 16 s +
certificate) to be attacked on their own terms, and the certificate's
enumeration removed (the O(f^4)) without touching the answer.

## 9. THE O(f^4) ATTEMPT: RIGHT ANSWER, WRONG COST (a recorded negative)

The certificate's enumeration (`|called| x |spectrum|`, 46.5M entries,
14.8 s at the degenerate locus) exists to feed two questions: `k` (the
cluster count) and the `excluded` union. `k` is exactly the number of
connected components of the graph on {winner} + tied classes with an edge
where `d(L,R) <= cut` (the lists encode that relation, and
`differing_observed_mass` is symmetric), and `excluded` is the union of the
bands. Both are answerable from the same bounds without ever building the
lists — no threshold, no approximation.

IMPLEMENTED (`tie_cluster_fast`, `IMPG_TIED_FOLD_FAST`, default OFF) and
GATED:

  **calls.jsonl and calls.jsonl.sequences.jsonl BYTE-IDENTICAL** to the
  committed path at T00048 (24 loci, the degenerate locus included), and
  both band-oracle tests now compare the fast path against the naive band
  form at both tie shapes.

**AND IT IS 20x SLOWER — the attempt is a recorded negative, not a fix:**

| phase | exact | fast |
|---|---:|---:|
| qual_spectrum (contains the enumeration) | 16.00 s | 1.25 s |
| qual_certificate | 14.79 s | **0.00 s** |
| qual_cluster | 0.95 s | **299.32 s** |

The enumeration already shares its per-class work: one multiset merge per
tied class, reused across all spectrum positions, with the bounds deciding
almost every position. A direct pairwise component sweep gives that sharing
up — it pays a per-pair distance decision instead of a per-(class,
position) one, and the surviving undecided-band walks dominate. So the
enumeration is the cheaper way to answer these questions **at this shape**,
and the O(f^4) stays open.

The next attempt should profile the cluster step first (the phase census
already brackets it) and keep the per-class sharing: a per-class *sorted
band scan* that stops at the cut rather than materializing the set, or a
bound-refined partition of the spectrum that lets one scan serve many
classes.

## 10. THE PER-LANE MEMORY: IT IS THE PANEL, NOT THE WORKING SET

The lane ceiling (which sets the whole-genome wall) is per-lane RSS. Measured
on one map leg, same territory, same box:

| configuration | peak RSS |
|---|---:|
| the fleet parent holding the panel alone (the COW receipt) | **33.88 GB** |
| map leg, shard scan (pre-store) | 36.40 GB |
| map leg, evidence store, `IMPG_AGC_DECOMPRESSOR_SLOTS=1` | 39.70 GB |
| map leg, evidence store, + `drop(agc)` at the frame setup's end + slot affinity | 39.30 GB |

**The global panel is 86% of it.** Everything the cut and the census hold -
the evidence packets, the packet arena, the territory index, `FileCut`, the
census JSON, the local syng, `locus.fa` - is ~2-6 GB.

A static audit had argued the dominant term was eight AGC decompressor
copies (~30 GB at ~3.7 GB each, the pool's round-robin opening one per
fetch). **The measurement refuses that**: forcing one slot changed the peak
by nothing (39.70 GB), and slot affinity changed it by 0.4 GB (39.30 GB).
The reason is visible in the numbers: the panel is already 33.88 GB, so a
few decompressors cannot be 30 GB on top of it. The two fixes were kept -
they are correct hygiene (one archive copy instead of up to eight; the AGC
released after its last reader, which is what the compiler confirms) - but
they are NOT the lever, and the dashboard above is why.

**THE LEVER IS SHARING THE PANEL, NOT TRIMMING THE WORKING SET:**

| lanes | panel paid privately (33.88 GB each) | panel shared (fork-COW) |
|---:|---:|---:|
| 16 | 542 GB - fits | 34 + 16x5.4 = 120 GB |
| 28 | 949 GB - barely | 185 GB |
| 64 | **2,168 GB - impossible** | **~380 GB - fits** |

So "64 threads / ~5 hours" is reachable, but only as a **fork-COW map fleet**
- one panel load, one child per territory, children inheriting the panel the
way the call fleet already does (`fork-cut-call-batch-global`). Thread count
per lane buys nothing: a map leg already saturates its cores (153 of 192
when alone). The present `cut-batch-seq` amortizes the *load* by running
territories serially inside one process; it does not add cross-territory
parallelism, which is what a COW child-per-territory fleet would.

The floor under all of this is the panel's own 32.36 GB of shared arenas
(GBWT 12.87 + khash 4.31 + the mmap'd sidecars) plus a ~185 s load. Only an
index-format change moves that.

## 11. TWO MEASUREMENT CORRECTIONS THAT CHANGE THE ARITHMETIC

Both from one map leg, same territory, same store, same output.

### 11.1 THE DEFAULT RAYON POOL WASTES 98% OF ITS CPU

| configuration | core-s | wall | cores avg |
|---|---:|---:|---:|
| default rayon (all 192 threads) | **17,154** | 395 s | 43.5 |
| `RAYON_NUM_THREADS=12` | **381** | 289 s | 1.3 |

The same work, the same output, **45x the CPU and 1.37x the wall** at the
default pool. The wall is worse *and* the CPU is absurd: the extra threads
are contending, not working. So **core-seconds measured at the default pool
are not a measure of work** - and any earlier "core-s" figure in this
document (or in the campaign's ledgers) taken at the default pool inherits
that inflation. The routes derivation's own 163,159-vs-62,169 core-s swing
between thread counts is the same phenomenon.

The consequences:
* **The campaign's measured throughput was accidentally right for the wrong
  reason.** 88 territories/hour at 16 lanes is consistent with ~7,800
  core-s/leg, not with the 106,496 that a leg reports when it runs alone
  over the whole pool.
* **Wall at a fixed concurrency is the only trustworthy cost metric here.**
  Quote (lanes x rayon-per-lane, wall/territory, throughput/hour); do not
  quote core-s unless the thread count is stated and sane.

### 11.2 THE STORE-FED SCAN IS SERIAL - and that is now the opportunity

With the store the cut reads ONE packet stream (the territory's fetched
evidence) instead of nineteen sharded ones, so the scan no longer
parallelizes: the leg averages 1.3 cores. That is *why* it is cheap in core-s
and why the wall did not fall further.

**So the honest current arithmetic, at a sane thread count:**

| | value |
|---|---:|
| work per territory (map leg) | **381 core-s** |
| whole sample, map only | 2.2M core-s |
| floor at 192 cores | **3.2 h** |
| lanes that fit at 41 GB/lane | ~24 -> ~20 h |
| lanes that fit if the panel is fork-COW shared (~7 GB private/lane) | **~140 -> ~3-4 h** |

**That is the 5-hour target, and it is reachable** - not by more threads and
not by more lanes on a private panel, but by (a) sharing the panel across
lanes (COW, the call fleet's proven form) and (b) giving the store-fed scan
its parallelism back, since today it leaves ~190 of the 192 cores idle.

## 12. NODE->READ POSTINGS: THE MEASURED GATE (and what it is not)

`.n2r` is a per-sample sibling artifact: node -> the read ordinals whose walk
touches it. Built once in **409 s** from the store (22.4 GB read, 168.4 GB
decompressed, 349 MB peak RSS via a hash-partitioned spill), **3.41 GB**, and
it leaves `sample.evd*` byte-untouched. Ordinals are the store's own domain,
which IS globally unique in the whole-sample store (`evidence_index.rs:874-881`
increments one counter across all shards) - unlike the per-shard prototypes,
whose ordinals are per shard.

The cut then refuses non-candidate records BEFORE decode/derive. The per-read
predicate is untouched, so this is exact by construction, and measured so on
T00048 (same thread count, RAYON=12):

| | without postings | with postings |
|---|---:|---:|
| packets examined | 14,949,866 | 14,949,866 |
| candidate reads | 14,949,866 | **7,157,270** (47.9%) |
| kept | 276,042 | **276,042** |
| cut phase | 68.2 s | **49.9 s** (-27%) |
| map leg wall | 289 s | **228 s** (-21%) |
| map leg core-s | 381 | **294** (-23%) |
| cut-census / cut-cache / cut-reads | | **BYTE-IDENTICAL** |

### 12.1 WHAT IT IS NOT: the bucket sweep was never the problem

The premise that motivated this work - "the fetch is 54:1, so an exact
node->read index makes it 1:1" - is WRONG, and the measurement says why: the
territory's node set is **145,441 nodes**, and fully half the reads sitting in
its locus-granular buckets genuinely touch one of them. The postings can only
discharge the other half; the remaining 3.9% keep-rate is the per-read
predicate's own work, which no index above it can do.

So: **-27% on the cut, not 47x.** Bounded by construction, because the bytes of
a candidate read only exist inside its own bucket and the bucket is a ~10 kb
locus, not a node.

### 12.2 AND NOW THE LOAD DOMINATES

With the cut at 50 s, the map leg is **183 s of panel load + 50 s of cut**: the
load is 80% of the leg. That is why the batch behaviour (section 9) matters more
than the cut does - and why any further cut-side work should be justified
against 50 s, not 191 s.

### 12.3 THE BUILD-CHAIN TRAP, AGAIN (it cost a full run)

`impg infer map` does not run the cut itself: it execs
`target-frontier/release/examples/sample_projection`. Rebuilding `--bin impg`
therefore measures a STALE child, and `strings .../impg | grep node-postings`
returns 0 even when the code is current, because `run_cut` is only reachable
from the example binary and the linker strips it from `impg`. The only valid
staleness check on this path is:

    strings target-frontier/release/examples/sample_projection | grep node-postings

Rebuild `--bin impg` AND `--example sample_projection` before believing any
cut-side measurement.

## 13. THE DEGENERATE CERTIFICATE, THREE MEASURED CUTS

`qual_certificate` at T00048 locus 2 (f=124 folds, 7,749 tied classes) was
19.68 s. Its structure, measured rather than assumed:

    60,052,001 pair tests
    46,456,131  TRUE predicates, admitted by the O(1) upper bound  <- this IS the answer
    13,592,xxx  needed the exact distance walk (1.45 us each)      <- 16.3 s of the 16.7
    ...and the emitted predicate lists collapse to 9 distinct member sets
       (6,786 + 582 + 232 + 116 + 14 + 10 + 5 + 3 + 1 = 7,749: a partition into cliques)

| cut | mechanism | effect | identity |
|---|---|---|---|
| spectrum shape cache | build each entry's merged multisets once, index them | 19.80 -> 16.84 s (-15.0%) | byte-identical |
| observed-mass arrays | cache `observed(key)` per entry; the walk multiplies | 16.67 -> 13.30 s (-20.2%) | byte-identical |
| unit-column fold (flag) | drop constant ll columns, merge bit-identical ones | see section 4 | not bit-identical (1,362 ulp) |

Cumulative: **19.68 -> 13.30 s (-32%)**, with `calls.jsonl`,
`calls.jsonl.sequences.jsonl`, `cut-census.jsonl` and all four instrument
receipt parts byte-identical at every step.

### 13.1 THREE HYPOTHESES OF MINE DIED HERE - record them

1. "423 ns/pair means the loop is 10-100x off; it must be allocation churn."
   The cache removed only 0.22 us of the 1.45 us per test. The other ~1.2 us
   was the walk itself (a HashMap probe per merged key), i.e. the WORK.
2. "Caching the probes should be 2-4x." It was 20%.
3. "The component-only path (`tie_cluster_fast`) is 20x slower because of a
   different reason." It carried the SAME per-test construction defect as the
   default path. The original 299 s negative was measured against the buggy
   code and must not be treated as a property of the algorithm until the
   cached variant is timed.

### 13.2 A FLAKY TEST IN THE SUITE IS A GATE RISK

`resolution::tests::abpoa_ordered_sequence_helper_preserves_exact_paths_in_header_order`
is order-dependent: the same tree gives

    run 1: 602 passed; 1 failed
    run 2: 603 passed; 0 failed

It passes in isolation. Every byte-identity claim in this document rests on
`cargo test --release --lib`, so a flaky member means a real regression can
hide in a green-looking run. Treat a single green run as insufficient.
