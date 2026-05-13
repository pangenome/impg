# Syng vs PAF validation — run #1

**Date:** 2026-04-16
**Query:** `S288C#0#chrIV:408000-410000 -d 10000`
**Data:** yeast235 (235 haplotypes)
**Binaries:** commit `d48afc8`

## Setup

```bash
# Build sparse all-vs-S288C PAF via sweepga --pairs
agc listset yeast235.agc | awk '$1 != "S288C#0" {print "S288C#0\t"$1}' \
    > /tmp/pairs_S288C.txt
sweepga yeast235.agc --pairs /tmp/pairs_S288C.txt -t 8 \
    --output-file yeast235_vs_S288C.paf
# 234 pairs, ~880k unfiltered alignments → 88,762 after filtering, 245MB PAF

impg index -a yeast235_vs_S288C.paf

# PAF-based query
impg query -a yeast235_vs_S288C.paf --sequence-files yeast235.agc \
    -r 'S288C#0#chrIV:408000-410000' --transitive -d 10000 -o bed \
    > paf_query.bed

# Syng-based query
impg query -a yeast235 --sequence-files yeast235.agc \
    -r 'S288C#0#chrIV:408000-410000' -d 10000 -o bed \
    > syng_query.bed

# Diff
tests/validation/compare_syng_vs_paf.sh syng_query.bed paf_query.bed
```

## Headline results

| Metric | Value |
| ------ | ----- |
| Syng rows | 269 |
| PAF rows | 228 |
| Common paths | 209 |
| Syng-only paths | 58 |
| PAF-only paths | 19 |
| Strand agreement | **209/209 (100%)** |
| \|start Δ\| ≤ 5bp | 124/209 (59.3%) |
| \|start Δ\| ≤ 20bp | 128/209 (61.2%) |
| \|end Δ\| ≤ 5bp | 155/209 (74.2%) |
| \|end Δ\| ≤ 20bp | 156/209 (74.6%) |
| max \|start Δ\| | 1583bp |
| max \|end Δ\| | 254bp |

## Interpretation

**What agrees:** strand assignment is 100% consistent between the two
methods on the 209 shared paths. This validates the RC detection
(`query_region_with_anchors` both-orientation GBWT query) and the
majority-strand dedupe. No phantom inversions from either side.

**What mostly agrees:** ~60% of common paths have `|start delta| ≤
5bp` and ~74% have `|end delta| ≤ 5bp`. The tight cases are the
well-behaved collinear homologies where both pipelines converge on
the same coordinates.

## Residual disagreement — four observed patterns

The non-zero deltas cluster in specific patterns, suggesting
systematic (not random) differences:

### 1. `+` strand paths, small shifts (≤ 300bp)

Common patterns:
- `start=+97 end=0`    ~20+ paths (CBK, CBM, SK1, W303, Y55, …)
- `start=0   end=-78`   ~10+ paths (ABA, ABH, AEH, AKR, BHH, CDG_1a, …)
- `start=0   end=-192`  ~15+ paths (BAI_1a, BAK_1a, BAL_1a, CCN, HN1, …)

These are small systematic shifts across groups of closely-related
strains that likely share identical indel structures in the probed
region. Plausible reading: syng and PAF each make slightly different
choices about where homology ends at the boundary of a
strain-specific indel block. Not a bug in either, but resolvable with
a follow-up if we want sub-100bp agreement.

### 2. `+` strand paths, moderate shifts (~115, ~245)

- `start=+115 end=-254`  8 paths (ALI, BCE_3a, BCN, CPI_1c, CQS_1a)
- `start=+245 end=-192`  4 paths (AMH_1a, BBM_1a, CEI_1a, EM14S013B,
  XXYS14)

Similar story — specific indel blocks handled slightly differently.

### 3. `-` strand paths, **LARGE** start shifts (~1450–1600bp)

Here the disagreement is large:

```
BAD#4#block15_contig1  syng=12473..12890(-)   paf=10890..12890(-)    ΔS=+1583
BAD#3#block8_contig1   syng=93907..94324(-)   paf=92324..94324(-)    ΔS=+1583
BTE#4#block86_contig2  syng=57343..57894(-)   paf=55891..57894(-)    ΔS=+1452
BBT#4#block55_contig1  syng=1064653..1065087(-) paf=1063201..1065201(-) ΔS=+1452
CGH_3#4#block18_contig1 syng=575619..575922(-) paf=574036..576036(-) ΔS=+1583
...
```

Consistent pattern on `-` strand: syng's target intervals are
~1500bp **NARROWER** than PAF's on the left (= the projection of
`qe=410000` under RC). End coordinate (= projection of `qs=408000`)
agrees. So syng correctly anchors the right edge but trims ~1500bp
off the left.

PAF intervals are consistently ~2000bp (matching the query width).
Syng intervals are ~400–550bp.

**Hypothesis:** under `-` strand, the innermost anchor on query axis
(`rightmost_anchor` for qe projection) is probably sitting in the
middle of the homology, not near its left edge, because syng's
shared-syncmer coverage is sparser on the RC side (known: closed
syncmers aren't RC-symmetric). Linear projection from that anchor
under RC lands in the middle of target; edge realignment doesn't
extend outward far enough (capped at 4kb but with a default buffer
sized for forward homology). PAF's stored CIGAR captures the full
extent of the RC alignment directly.

This is a real issue worth fixing, and it's a specific fix: for
`-` strand edges, extend the realignment window OR relax the
co-linearity tolerance so far-edge anchors aren't filtered as
paralogs when they're actually valid terminal anchors of the same RC
homolog.

### 4. Path-set disagreement (not boundary)

- **Syng-only (58 paths)**: dominated by oddly-named contigs like
  `ANL#3#block107_contig3`, `AIF#1#chrIV`, etc. These are
  non-PanSN-named contigs that syng finds via syncmer sharing but
  that didn't get aligned by FastGA in the sweepga run — possibly
  because the contig is too short, too divergent, or in a repeat
  region that FastGA filters.
- **PAF-only (19 paths)**: CBS432, CBS7001, CL216, CL450, DG1768,
  FM1318, HNe8, and similar well-named chrIV contigs that sweepga
  aligned but syng's query didn't discover at this query region.
  Likely explanation: these are strains where the chrIV homolog of
  `S288C:408000-410000` has enough sequence divergence that shared
  syncmers are too sparse for syng to report above the padding
  threshold.

Neither is necessarily a bug — they're each picking up homologies the
other missed, for different reasons.

## Actionable findings

1. **[High] Fix `-` strand edge realignment.** The ~1500bp systematic
   left-edge shortfall on 10+ RC paths is the largest disagreement
   and the one most likely to affect graph quality. Plan:
   - Probe one of the affected paths (BAD#3#block8_contig1) with
     `syng_anchor_probe` to see where the anchors actually sit.
   - If the innermost anchor on the RC side is far from the query
     boundary, extend `EDGE_ALIGN_CAP_BP` (currently 4096) for `-`
     strand OR enlarge the target buffer for RC.
   - Re-run validation and check that the max |start Δ| on `-`
     strand drops.

2. **[Medium] Understand sample coverage gaps.** Write small scripts
   to spot-check one syng-only and one PAF-only case — confirm the
   explanations above. Not likely to change the algorithm but worth
   understanding.

3. **[Low] Resolve sub-100bp `+` strand shifts.** The `97`, `78`,
   `115`, `192`, `254` clusters probably reveal whether edge
   realignment is properly extending into small indels. Compare on
   a known-indel fixture.

## Conclusion

Strand assignment: correct.
Forward-strand boundaries: good (~60-75% within 5bp).
RC-strand boundaries: systematically undershoot by ~1500bp on the
outer edge for a specific class of paths. This is a real bug worth
fixing before declaring the syng query pipeline validated.
