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
