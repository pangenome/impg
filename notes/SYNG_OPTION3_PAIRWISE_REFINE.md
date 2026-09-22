# Option 3: syng-for-discovery + per-homolog pairwise refine

**Status:** Design study.
**Context:** `notes/SYNG_VS_PAF_VALIDATION.md` documented a ~1500bp
systematic undershoot on `-` strand RC homologs — syng's anchor grid
is too sparse under RC to support accurate linear projection at the
boundary, and my edge-realignment step is capped too tightly to
recover the real extent. This note studies the alternative: use syng
*only* to identify which target paths are homologous, then do a fresh
pairwise alignment per homolog to recover PAF-quality boundaries.

---

## 1. The fundamental limit syng alone cannot cross

Syng's anchor graph encodes homology at syncmer resolution. Shared
syncmers between query and target paths are:
- Dense on collinear forward-strand homology in closely-related
  genomes (~35/kb).
- Sparse (~3–5/kb) on RC homology, because closed syncmers aren't
  RC-symmetric — the minimizer k-mer of a window depends on
  orientation.
- Sparse in divergent regions (each SNP that lands in a syncmer
  disrupts it).
- Absent in transposon-rich or subtelomeric regions with few
  informative k-mers.

**Crucial property:** indels DISRUPT syncmers (they change k-mer
content) rather than sit between shared ones. So the shared-anchor
grid has zero indel information *within* it — between any two
shared anchors on the same path, syng cannot tell you whether there's
a 10bp insertion or a perfect match.

Linear projection from anchors assumes zero indels between anchor
and edge. That's fine when the anchor sits <30bp from the edge
(indel burden negligible). It fails when the anchor is 500bp+ from
the edge (typical on RC paths), because any indel in that gap shifts
the real edge somewhere the anchor grid can't show.

Edge realignment (current impl) is a bounded BiWFA around the
outermost anchor. It extends linear projection by the size of the
fetch window; beyond that, it degenerates. We observe ~1500bp
undershoots on RC paths, meaning the real edge is 1500bp beyond where
both linear projection AND our bounded edge realignment reach.

**PAF-based queries don't have this problem**: FastGA's stored CIGAR
describes every base of the alignment between query and target. The
indel content is in the CIGAR directly; impg's query projects through
it trivially.

To get PAF-quality boundaries from a syng-seeded query, we have to
produce a PAF-quality alignment at query time. That's what this study
is about.

---

## 2. Proposal

**Syng is a discovery oracle only.** It tells us: for the query
region, these N paths have at least one shared syncmer, and the
approximate syncmer-resolution extent on each path is `[homolog.start,
homolog.end]` with strand `s`.

**Per-homolog pairwise alignment** then takes that seed and produces a
full CIGAR between `query[qs, qe]` and `target[homolog.start - pad,
homolog.end + pad]` (RC'd if `s == '-'`). The CIGAR tells us exactly
where qs and qe land on target. No linear projection. No anchor-
dependent approximation.

### Per-homolog pipeline

```
For each merged homolog (produced by query_region_with_anchors +
                        distance_merge + strand_dedupe):

    Fetch query bytes:   q = query[qs, qe]                     (width W_q)
    Fetch target bytes:  t = target[homolog.start - pad,       (forward strand
                                    homolog.end   + pad]        coords, both strands)

    If homolog.strand == '-':
        t = reverse_complement(t)

    Align q against t:
        semi-global, one of:
          (a) EndsFree on both target ends, pattern end2end
          (b) End2End if syng bounds are tight
        Pick the one that matches the biology best — syng bounds are
        approximate (padded syncmer resolution), so (a) is safer.

    Read the CIGAR:
      - Number of leading 'I' ops (target consumed before query starts)
        → offset of qs on target forward strand.
      - Number of trailing 'I' ops (target consumed after query ends)
        → complement for qe.
      - For '-' strand, the CIGAR is against RC'd target; mirror the
        offset back to forward-strand coords.

    Output refined (target_start, target_end, strand) with bp-precise
    boundaries.
```

### Key algorithmic choices

1. **Target fetch size.** Syng's `homolog.start/end` bounds the
   region but can be off by hundreds of bp (especially under RC).
   Pad the target fetch by, say, `max(1000bp, W_q)` on each side
   beyond syng's bounds, so the real alignment endpoint is within
   the fetched window. Bounded cost because we clamp to sequence
   length.

2. **Alignment form.** Two-sided EndsFree: both left and right ends
   of target are free (target_begin_free = target_end_free = pad).
   Query is End2End (pattern_begin_free = pattern_end_free = 0).
   This lets BiWFA find the best position for the query inside the
   fetched target window without penalty for the flanking context.

3. **Aligner.** lib_wfa2 edit-distance, MemoryMode::High (small
   windows, not a concern). Could move to GapAffine for more
   biologically meaningful alignments but edit-distance matches
   PAF's typical CIGAR semantics closely enough for boundary
   projection.

4. **Score thresholding.** If the alignment score is too low
   (e.g., edit distance > 10% of query length), the "homolog" is
   probably syncmer noise rather than real homology. Fall back to
   syng's padded bounds, or drop the homolog entirely with a flag.

5. **CIGAR → target bounds projection.**
   - Scan from left: count consumed target bases until query starts
     (leading 'I'). That's `target_offset_for_qs`.
   - Scan from right: count trailing 'I's. Complement gives
     `target_offset_for_qe_exclusive`.
   - Target forward-strand bounds:
     - `+` strand: `[t_fetch_start + offset_for_qs, t_fetch_end -
       trailing_I_count]`.
     - `-` strand: target was RC'd for alignment. The fetched
       forward-coord window is `[t_fetch_start, t_fetch_end)`. In
       oriented coords, position 0 corresponds to `t_fetch_end - 1`.
       So `qs_offset_in_oriented → t_fetch_end - qs_offset` (forward-
       strand coord of qs). Mirror for qe.

---

## 3. Corner cases

### 3a. Partial homology
Query `[qs, qe)` may extend past the homology on one or both sides.
The alignment will show trailing/leading gaps (high edit distance at
those ends). Two handling options:
- **Truncate**: the refined interval covers only the aligned portion,
  not the full query. Set refined_start to the first matched base's
  target position, refined_end to the last matched base's.
- **Extend to the padded bound**: treat the query's outside-homology
  portion as "no projection data" and fall back to syng's padded
  edge. Easier but less accurate.

Prefer *truncate*. Output is cleaner for downstream GFA construction
(no ragged edges).

### 3b. Divergent / low-identity homology
If the aligned edit distance is > some threshold (proposed 20% of
query length), the "homology" is probably syncmer noise or a very
distant paralog. Options:
- Drop the homolog.
- Keep it but flag it.
- Log it and let the user decide.

Prefer *drop* for graph-quality outputs; *keep-and-flag* for bed
output where the user is exploring.

### 3c. Repeats / multimapping
The target region may contain several copies of a query-similar
sequence. BiWFA will pick one alignment and stick with it. If the
true homolog is elsewhere, this could mispoint. Mitigation:
- Syng's merge already splits paralogs into separate homologs by
  the co-linearity filter. Each one gets its own alignment, so this
  case is typically handled at the homolog level.
- For repeats WITHIN a single homolog (e.g., a tandem duplication
  in one strain), alignment picks the best single mapping. That's
  correct for boundary projection, just not ideal for graph
  representation. Acceptable for this phase.

### 3d. Multihop
At hop k+1, syng seeds come from hop k's refined (tight) intervals.
Because hop k's boundaries are now PAF-quality, the syng query at
hop k+1 operates on tight ranges, and slop doesn't compound across
hops. This actually gets simpler, not more complex, with pairwise
refinement.

---

## 4. Cost analysis

Per query region:
- N = number of merged homologs after syng + strand_dedupe.
  On yeast235 2kb query: N ≈ 270.
- Per homolog: 2 AGC fetches (query stays constant across homologs
  so can be cached → only 1 target fetch per homolog) + 1 BiWFA
  alignment on ~(W_q + 2×pad) × W_q bases.
  For W_q = 2000, pad = 1000: 4000 × 2000 BiWFA. Edit-distance
  BiWFA runs in O(ns) time where s is the edit distance. For ~5%
  divergence on 2kb, s ≈ 100, n = 4000 → ~400k wavefront ops ≈
  a few ms.
- Total: 270 × (fetch + ~5ms BiWFA) ≈ 10–30 seconds.

AGC fetch is likely the dominant cost (similar to current edge
realignment on target fetches). BiWFA itself is cheap.

vs. current edge realignment on the same query: 5 seconds with 540
small-window fetches. Option 3 does 270 larger fetches. Probably
trades to ~10-15s wall total. Acceptable.

---

## 5. Accuracy expectations

### Forward homology, low divergence
Current syng+linear+edge-realign gets |start Δ| ≤ 5bp for ~60% of
these paths. Option 3 should get it to ≥95%, matching PAF precisely
for the collinear cases.

### RC homology
Current: ~1500bp undershoot on outer edge for 10+ paths. Option 3
should eliminate this entirely — the alignment covers the full
homology regardless of anchor density.

### Divergent regions
Current: edge realignment runs off its bounded window, returns
padded bound.
Option 3: alignment either finds the real edge (if divergence is
handleable by BiWFA) or reports low score → drop/flag (if not).

### Expected PAF-vs-syng diff after option 3
- ≥ 95% of common paths within 5bp on both edges
- ≥ 99% within 20bp
- Remaining disagreement should be in:
  - Truly divergent regions (alignment quality-limited)
  - Transposon-rich regions (repeats cause localized wrong mappings)
  - Sample-coverage gaps (orthogonal — different source of
    disagreement)

---

## 6. Implementation plan

Not coding yet. Just the sequence of steps when we're ready:

1. **Extract the per-homolog alignment helper into its own function.**
   `fn refine_homolog_by_alignment(hit, query_bytes, sequence_index,
   query_edges) -> Option<(u64, u64)>`. Returns refined forward-
   strand (start, end) or None if alignment fails / score too low.

2. **Replace the body of `refine_boundaries`.** Keep the signature
   so call sites in `one_hop` don't need to change. Internally, if
   `sequence_index` is provided, call the new helper. If it fails
   or returns low score, fall back to linear projection (current
   behavior) with a debug log.

3. **Cache the query bytes per-hop.** Fetch `query[qs, qe)` once at
   the top of `one_hop`, pass into the refiner. Avoids 270 redundant
   AGC fetches of the same query region.

4. **Decide the score threshold.** Start with "reject if edit
   distance > 30% of query length" (very permissive; only filters
   obvious junk). Tune on the yeast235 validation data.

5. **Re-run validation.** Expect:
   - 95%+ paths within 5bp on both edges.
   - Max |Δ| drops from 1583bp to single digits.
   - Some syng-only paths may drop (scored below threshold).

6. **Integration test.** `test_syng_query_reconstructs_homology_with_diffs`
   already exists for the simple case. Add a specifically-RC fixture
   similar to `test_syng_rc_homolog_end_to_end` but verify exact
   coordinate agreement now (not just base-content overlap).

7. **Remove or demote current edge realignment.** Once option 3 is
   the refinement path, `align_edge_to_anchor` and the
   ANCHOR_FLANK/EDGE_ALIGN_CAP constants are dead code for the
   default path. Either delete them or keep them behind a debug flag
   for comparison.

---

## 7. Open questions before coding

1. **Scoring / alignment form.** EndsFree with both target sides
   free is the safest bet. But might cost us precision if the
   padding is very asymmetric. Empirical — try on yeast data and see.

2. **Query cache granularity.** Should we cache query bytes at the
   `one_hop` level, or at `query_transitive` level (across hops)?
   Given hop count is usually 1, one_hop is fine.

3. **Gap-affine vs edit.** Edit distance is simpler but punishes
   indels linearly. Gap-affine is more biological but needs
   scoring parameters. For boundary projection (not biological
   interpretation), edit distance should be adequate.

4. **Target fetch padding default.** `max(1000bp, W_q)` is generous.
   Tune based on the yeast RC undershoot: we saw ~1500bp, so pad
   should comfortably exceed that. Proposal: `max(2000bp, W_q × 0.5)`.
   Revisit after first run.

5. **Multithreading.** 270 homologs is embarrassingly parallel.
   Rayon-ize the loop? BiWFA aligner is thread-local already. Would
   probably cut wall time to ~3-5s instead of 10-15s.

6. **Semantic of "refined homolog bounds"**. Should refined bounds
   represent the full length of aligned query on target, or the
   target portion aligning to the user's `[qs, qe)` specifically?
   They're almost the same when query fully lies inside the homology;
   they differ when query extends past. Decide up-front for
   consistency.

---

## 8. Relationship to the existing roadmap

`notes/SYNG_NEXT_STEPS.md` lists three items:
- #1 PAF validation — done, produced the data that motivates this study.
- #2 smoothxg inflation — unaffected by this work.
- #3 partitioned-POA engine — unaffected; this is about the syng
  *query* side, not the graph-engine side.

Option 3 is effectively a refinement of #1's output. Once it lands,
re-run validation and it should tighten from ~60% within 5bp to ~95%.
