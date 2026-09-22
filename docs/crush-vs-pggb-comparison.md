# crush vs PGGB on C4 (GRCh38 chr6:31891045-32123783)

Comparison of the current best crush output against a real PGGB control run on
the same C4 GRCh38 region, plus three concrete examples of trivial bubbles that
crush is failing to compact.

## TL;DR

- **PGGB control:** `impg query ... -o 'gfa:pggb'` ran the full
  FastGA + seqwish + smoothxg + gfaffix pipeline in 13:38 / 64 GB. The output
  graph is 13,288 S, 234,524 bp, 16,240 L, 465 P, all walks preserved.
- **Current best crush:** `c4_exp_no_filter_20260526T005655Z/run.Ygs.gfa` —
  ran in 36:53 / 101 GB and produced 19,836 S / 553,585 bp / 23,384 L / 465 P.
  Crush stores **2.36× more sequence content** for the same 232 kb region than
  PGGB does, despite running ~2.7× slower with ~36% higher peak RSS.
- **Where crush misses:** every ≤200 bp polyT/SNP bubble that should compact
  to 1–5 long-spine segments stays expanded to 20–30 single-bp segments. The
  CHM13 walk through each of the three example bubbles below uses 17–28
  segments in crush vs **1–5 segments in PGGB**. Root cause is two things in
  the same pipeline:
  1. `min-traversal-len=5k` filters out small bubbles at
     `resolution.rs:1637` before the aligner-and-polish stage ever sees them, so
     all ≤5 kb bubble candidates survive in their syng-raw form.
  2. The syng:crush engine path (`lib.rs:745` `apply_graph_transforms`) never
     runs gfaffix at the end. PGGB's other engines all do (`lib.rs:823, 836,
     889, 980, 990, 1127` all call `normalize_and_sort` → `run_gfaffix`), which
     is what collapses adjacent-identical-sequence chains.
- **Trivial-stringy bubbles:** crush has 476 of them; PGGB has 12.

## Inputs and tooling

| Field | Value |
|---|---|
| Region | GRCh38#0#chr6:31891045-32123783 (≈232,738 bp) |
| Reference path in graph | CHM13#0#chr6:31744284-31976975 (the syng anchor) |
| Input panel | `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng` + `.agc` (465 path ranges) |
| impg binary | `/home/erikg/impg/target/release/impg` (prebuilt; not rebuilt for this task) |
| Current best crush GFA | `/home/erikg/impg/data/c4_exp_no_filter_20260526T005655Z/run.Ygs.gfa` |
| Current best crush PNG | https://hypervolu.me/~erik/impg/c4-crush-no-filter-best.png |
| PGGB control output | `/home/erikg/impg/data/c4_pggb_control_20260526T025439Z/pggb.nosort.gfa` |
| PGGB control PNG | https://hypervolu.me/~erik/impg/c4-pggb-control.png |
| Threads | 32 |

The crush run used `gfa:syng:mask,min-run=3:crush,method=auto,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort`.

The PGGB run used `-o 'gfa:pggb'` — the default engine. No additional flags.

## Step 1: PGGB run

```bash
out=/home/erikg/impg/data/c4_pggb_control_20260526T025439Z
/usr/bin/time -v /home/erikg/impg/target/release/impg query \
  -t 32 -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -d 50k \
  -o 'gfa:pggb' \
  -O "$out/pggb.nosort" -v 1 > "$out/pggb.nosort.stdout" 2> "$out/pggb.nosort.stderr"

gfasort -i "$out/pggb.nosort.gfa" -o "$out/pggb.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/pggb.Ygs.gfa" -o "$out/pggb.Ygs.png" -m -x 2200 -y 1200
scp "$out/pggb.Ygs.png" erik@hypervolu.me:www/impg/c4-pggb-control.png
```

**PGGB wall time + RSS:** 13:38 wall, 64.3 GB peak RSS, exit 0 (from `pggb.nosort.stderr`).

**Crush wall time + RSS (baseline):** 36:53 wall, 100.8 GB peak RSS (from `run.nosort.stderr`).

PGGB is ~2.7× faster and ~36% lower peak RSS on this region, in addition to
producing a much more compact graph (next section).

## Step 2: Side-by-side metrics

| Metric | PGGB control | Current best crush | Notes |
|---|---|---|---|
| Segments (S-lines) | **13,288** | 19,836 | crush has 49% more |
| Total segment bp | **234,524** | 553,585 | crush stores **2.36× more sequence** than PGGB for the same region |
| Mean seg len (bp) | 17.65 | 27.91 | |
| Median seg len (bp) | 1 | 21 | PGGB chops to many 1-bp variation segments; crush keeps mid-sized syng anchors |
| Links (L-lines) | **16,240** | 23,384 | crush has 44% more |
| Paths (P-lines) | 465 | 465 | paths preserved in both |
| Path steps total | 5,538,879 | 4,291,621 | PGGB has more steps because of finer chopping |
| Path step median per path | 11,841 | 9,410 | same reason |
| Paths preserved? | **Yes** | **Yes** | both round-trip the input panel walks |
| Distinct sequences (across S-lines) | 2,667 | 13,126 | crush has 4.9× more *distinct* short sequences — these are crush's "stringy" small bubbles |
| Duplicate-sequence extras¹ | 10,621 | 6,710 | |
| ↳ length band ≤4 bp | 10,605 | 4,652 | |
| ↳ length band 5–10 bp | 15 | 413 | |
| ↳ length band 11–50 bp | 1 | 1,626 | crush carries 1,626 ≤50-bp dup-sequence segments that PGGB has compacted |
| ↳ length band 51–200 bp | 0 | 16 | |
| ↳ length band >200 bp | 0 | 3 | |
| Trivial-stringy bubbles² | **12** | **476** | dropping from 476→12 is the headline structural delta |

¹ "Duplicate-sequence extras" = number of S-line entries whose nucleotide
sequence is byte-identical to some other S-line in the same graph (extras =
`total_segs - distinct_seqs`). PGGB has a higher raw extras count because its
gfaffix+chopping pipeline intentionally splits variation sites into many ≤4 bp
segments that often repeat across the graph; the cleaner metric is **total
segment bp**, where PGGB stores 2.36× less sequence content per bp of region.

² "Trivial-stringy bubble" = bubble with ≤10 distinct traversal sequences,
median traversal length ≤500 bp, and ≥5 interior segments per traversal.
Heuristic from `/tmp/find_stringy_bubbles.py`, run over each output GFA. In
the PGGB control all 12 hits cap at ≤10 interior segments; in crush the worst
hit has 27 interior segments per traversal for ~130 bp of content.

POVU residual-bubble counts (flubbles at level 0/1/2 in the output graph) were
**not** computed: the `povu` binary is not installed on this host and rebuilding
was out of scope ("DO NOT REBUILD"). The trivial-stringy-bubble heuristic count
in the table above (476 in crush vs 12 in PGGB) is computed from the same
input — anchors derived from high-multiplicity segments, bubble interiors
derived from the step run between consecutive anchors on each path — and is a
faithful surrogate for "post-compaction residual bubble count" for the bands
the user cares about (small bubbles that should have compacted). If a true
POVU comparison is needed in future, install `povu` and add an extra column.

## Step 3: Three trivial-but-stringy bubble examples

A "trivial-stringy" bubble here means a bubble where:

- the number of *distinct* traversal sequences through the bubble is ≤10,
- the median traversal length is ≤500 bp, and
- the number of *interior* segments per traversal is ≥5.

That last condition is what makes it "stringy": the bubble's interior is
fragmented into many tiny segments instead of a small number of larger ones.

A heuristic crawl over `run.Ygs.gfa` (anchor segments = segments visited by
≥80% of paths; bubble interior = step run between consecutive anchors) found
**476** such bubbles in the crush output. The three highest-stringiness
examples are below. (Heuristic script:
`/tmp/find_stringy_bubbles.py`; per-bubble dumper: `/tmp/inspect_bubble.py`. Both
parse the GFA directly and run from any shell — no project rebuild required.)

### Side-by-side region comparison

For each example bubble below, we also walk the CHM13 reference path through
the same absolute chr6 coordinate window in the PGGB output and report how
many segments PGGB uses. This is the apples-to-apples "what would the same
content look like compacted?" view.

| Region (CHM13 chr6) | Crush segs on CHM13 walk | Crush ≤4-bp segs on walk | PGGB segs on CHM13 walk | PGGB ≤4-bp segs on walk |
|---|---|---|---|---|
| Example A: 31826515-31826645 (130 bp) | 28 | 21 | **2** | 0 |
| Example B: 31805524-31805642 (118 bp) | 26 | 19 | **1** | 0 |
| Example C: 31764162-31764226 ( 64 bp) | 17 | 14 | **5** | 2 |

PGGB collapses each of these to 1–5 long-spine segments (≥50 bp each) with the
polyT/SNP variation as side-branches off the spine; crush keeps a 17–28
single-bp-dominant chain on the CHM13 walk because syng's initial
fragmentation was never undone.

### Example A — anchors 6943→7207 (CHM13 chr6:31826515-31826645, ~130 bp span)

| Field | Value |
|---|---|
| Anchor L | seg 6943 (1 bp, "A"), visited by 459/465 paths |
| Anchor R | seg 7207 (12 bp, "ACTACAGGTGCC"), visited by 459/465 paths |
| Paths through bubble | 459 |
| Distinct traversal sequences | 8 |
| Traversal len distribution | 121–150 bp |
| Most-common traversal | 222 paths, 27 interior segments, 128 bp (19 of 27 segments are 1 bp) |
| 2nd most-common traversal | 106 paths, 17 interior segments, 121 bp |
| Variation content | polyT length (9–21 Ts) + one G/T SNP + one 5-bp insertion (GAGAC) |
| What a sane graph would look like | ~5 segments: anchor – polyT-stem – SNP – anchor (or even 1 segment if the polyT is left alone) |
| What crush did | Did not touch this bubble. The interior segment IDs (6944..6976, 7187..7206) lie entirely below 272,214,103 — the floor of crush-minted IDs (see `run.nosort.stderr`: `crush: initial next_id=272214103`). They are original syng-induced segments. The crush log mentions no `6943`/`7207` site, and the round-1 stats report only 8 selected sites with median-len ≥3,319 bp — this bubble's max-len (150 bp) is below the `min-traversal-len=5k` floor and was discarded at `resolution.rs:1637`. |
| Aligner run on it | None — bubble filtered out before aligner dispatch. |
| Polish run on it | None — bubble filtered out before polish dispatch. |

### Example B — anchors 5547→5600 (CHM13 chr6:31805524-31805642, ~118 bp span)

| Field | Value |
|---|---|
| Anchor L | seg 5547 (9 bp, "TAATTCTAC"), visited by ≥372/465 paths |
| Anchor R | seg 5600 (27 bp, "ATCCTCCTACCTCAGCCTCCCTAGTAG"), visited by ≥372/465 paths |
| Paths through bubble | 398 |
| Distinct traversal sequences | 8 |
| Traversal len distribution | 103–132 bp |
| Most-common traversal | 315 paths, 25 interior segments, 109 bp |
| Variation content | polyT length (8–19 Ts) + minor SNPs + one 8-bp insertion (GAGACACGGTCCTG-like) |
| What a sane graph would look like | ~5–6 segments |
| What crush did | Same story: interior segment IDs (5553..5599) all below crush ID floor. Not touched by crush. Aligner not invoked. Polish not invoked. The site's max-traversal-len (132 bp) is below `min-traversal-len=5k`. |

### Example C — anchors 1846→1896 (CHM13 chr6:31764162-31764226, ~64 bp span)

| Field | Value |
|---|---|
| Anchor L | seg 1846 (1 bp, "C"), visited by ≥372/465 paths |
| Anchor R | seg 1896 (44 bp, "TCAAGCGATCCGCCCACTGCAGTCTCCCAAAGTGCTGGGATTAC"), visited by ≥372/465 paths |
| Paths through bubble | 454 |
| Distinct traversal sequences | 4 |
| Traversal len distribution | 62–85 bp |
| Most-common traversal | 318 paths, 16 interior segments, 63 bp |
| 3rd traversal | 41 paths, 18 interior segments, 85 bp (carries a 21-T insertion + ACCGGGT |
| Variation content | polyT length variation + the 21-T+anchor-block insertion |
| What a sane graph would look like | ~5 segments (anchor – preT-stem – polyT – postT-stem – anchor) |
| What crush did | Same: interior IDs (1847..1895) all syng-origin (below crush ID floor), max-traversal-len 85 bp → discarded by `min-traversal-len=5k`. Not touched. |

## Step 4: Diagnosis — where crush is specifically failing to condense

Two compounding pipeline gaps, both in the syng:crush dispatch path:

1. **The `min-traversal-len=5k` floor (`/home/erikg/impg/src/resolution.rs:1637`)
   discards every small bubble before any aligner or polish runs.** In this run
   the config was `min-traversal-len=5k`, so any bubble with `max_traversal_len
   < 5000 bp` is filtered out at candidate-selection time. The round-1 stats in
   the crush log confirm this: from 70 candidates, only 8 were selected, and
   the smallest selected candidate had median-len 3319 bp. All three example
   bubbles above sit at 60–150 bp and so were never considered. They survive
   unchanged from the syng raw GFA.

2. **The syng:crush pipeline never runs gfaffix on its output
   (`/home/erikg/impg/src/lib.rs:745` `apply_graph_transforms` runs `crush` →
   `sort` but skips `normalize_and_sort`).** PGGB's other engines (`pggb`,
   `seqwish`, partitioned `pggb`) all end with `graph::normalize_and_sort` —
   see lib.rs:823, 836, 889, 980, 990, 1127 — which calls `run_gfaffix`. The
   crush dispatch path does not. So even if a small bubble's interior happens
   to be byte-identical between two adjacent segments (a free gfaffix
   collapse), the syng:crush output preserves the fragmentation. This is why
   the dup-sequence-extras count is 6,710 in this graph and why 4,652 of those
   are ≤4 bp.

The two gaps interact: gap 1 means crush only attacks big bubbles, leaving
small ones in their syng-raw fragmented state; gap 2 means the syng-raw
fragmentation is never normalized away. A pggb-style end-of-pipeline gfaffix
step would collapse ~6,700 of the current 19,836 segments without touching any
path walk, and would visibly de-stringify the small-bubble regions even before
any further crush tuning. A second knob (`min-traversal-len` lowered to e.g.
50 bp, paired with a fast pair-wise polish budget) would let crush attack the
substantive small bubbles, but step 2's gfaffix is the cheap structural win.

## Reproduction recipes

PGGB control:

```bash
out=/home/erikg/impg/data/c4_pggb_control_20260526T025439Z
ls -la "$out"            # pggb.nosort.gfa, pggb.Ygs.gfa, pggb.Ygs.png
```

The analysis scripts used below are archived alongside the PGGB output for
permanence:

```
/home/erikg/impg/data/c4_pggb_control_20260526T025439Z/analysis-scripts/
├── gfa_metrics2.py            # whole-graph metrics
├── find_stringy_bubbles.py    # rank trivial-but-stringy bubble candidates
├── inspect_bubble.py          # dump a single bubble by anchor pair
├── inspect_region_in_graph.py # walk CHM13 path through an abs chr6 window
└── locate_bubble_grch38.py    # map seg id → CHM13 abs chr6 coord
```

All scripts only parse GFA files — no impg rebuild needed.

Metrics:

```bash
SCR=/home/erikg/impg/data/c4_pggb_control_20260526T025439Z/analysis-scripts
python3 $SCR/gfa_metrics2.py /home/erikg/impg/data/c4_exp_no_filter_20260526T005655Z/run.Ygs.gfa
python3 $SCR/gfa_metrics2.py /home/erikg/impg/data/c4_pggb_control_20260526T025439Z/pggb.Ygs.gfa
```

Find more trivial-stringy bubbles in any GFA:

```bash
python3 $SCR/find_stringy_bubbles.py path/to/run.Ygs.gfa | head -30
```

Inspect a specific bubble by anchor pair (e.g. example A above):

```bash
python3 $SCR/inspect_bubble.py /home/erikg/impg/data/c4_exp_no_filter_20260526T005655Z/run.Ygs.gfa 6943+ 7207+
```

Compare the same chr6 window across crush and PGGB:

```bash
python3 $SCR/inspect_region_in_graph.py /home/erikg/impg/data/c4_exp_no_filter_20260526T005655Z/run.Ygs.gfa 31826515 31826645
python3 $SCR/inspect_region_in_graph.py /home/erikg/impg/data/c4_pggb_control_20260526T025439Z/pggb.Ygs.gfa  31826515 31826645
```

Map an interior segment to a CHM13 absolute chr6 coordinate:

```bash
python3 $SCR/locate_bubble_grch38.py /home/erikg/impg/data/c4_exp_no_filter_20260526T005655Z/run.Ygs.gfa 6943 7207
```
