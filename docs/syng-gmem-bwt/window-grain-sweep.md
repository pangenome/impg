# THE WINDOW-GRAIN SWEEP (2026-10-09) — the owner's ruling tested: 10kb / 20kb / 50kb windows x deep / 35x, five HG002 localities

The owner's ruling: the committed 10kb windows are too small; test 20kb and 50kb
(rationale: a 50kb window at 35x carries ~10kb-at-175x discriminating mass; coarser
windows also reduce the DP spine's chain boundaries). The sweep re-ran the committed
per-locality truth chain VERBATIM in stage shape and options at each new partition
grain (`impg partition -w <W> -d <W/10>`, the committed recipe's own 10:1 w:d ratio —
the grains are the owner's specified points, not tuning choices), over the SAME five
localities (the chr1 pilot + the four ordinary-genome localities L1-L4), the SAME
reads (the committed deep fetches AND the 35x corrected cuts,
`IMPG_CUT_LOCAL_DERIVE=1`), and the SAME machinery (fold rules, territory images at
the new window sizes, the CLI with `--truth-qv`).

## The substrate answer (measured at L1-w20k first, held everywhere)

The grain is GRAIN-LOCKED into the partition substrate, not a calling-layer
parameter: everything downstream of `impg partition -w` is per-grain (beds, catalog,
export graphs, pggb graphs, projection, locality windows, census, derive cache);
there is no knob in `call-diplotypes` that coarsens windows over an existing
substrate. The panel syng, the sample membwt and the truth-routes lanes are
grain-independent and REUSED; the routes are PARTIALLY grain-dependent (lanes, ports,
native/ and source-ports/ byte-identical across grains; ownership + registry track
the catalog's occurrence tiling — receipt `grain-sweep/routes-identity.txt`) and are
rebuilt per grain from the new-grain catalog. The one-time re-graining cost per
locality-grain, measured: ~9-11 min substrate chain + census 194-414s + derive-cache
~11s + CLI 34-125s + cut 408-484s (the full per-run walls are in the tables below).

The sweep also surfaced a new input class — degenerate non-axis SLIVER partitions
(all rows a few bp) that the coarse tilings emit — and one of them DETERMINISTICALLY
CRASHES the FastGA index builder (the engine-bug receipt at the end of this report;
the sweep's receipt-side containment builds exactly the provably-consumed axis set).

## THE COMPLETE GRAIN x DEPTH TABLES

(The generated artifact of record: `/home/erikg/hprcv2/grain-sweep/grain-tables.txt`;
`grain-tables.py` re-derives every cell from the CLI truth-QV receipts. Assessment-side
truth from `calls.jsonl.truth-qv.jsonl` ONLY — the CLI test mode; the 34-43 QV band is
the receipts' own named near-twin class, 34/60 the datasheet's named QUAL bands; no
thresholds, no tuning constants.)

### DEPTH: DEEP (committed per-locality fetches, 274-320x)

| locality | grain | loci | expr | rank1 | perfect | rank2 | neartwin-r2 (34-43) | QV med/p10/p90 | rank1-QV med | QUAL med(p25/p75;>=60/34-60/<34) | classes sum/max | cli_wall_s | census_wall_s | census_RSS_GiB | cut_wall_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pilot | 10kb | 22 | 22 | 4 | 4 | 3 | 3 | 25.77/2.45/60.0 | 60.0 | 156.53559774527022(53.56850658137531/1170.2344291317388;15/3/4) | 753649/205120 | - | - | - | - |
| pilot | 20kb | 12 | 11 | 2 | 2 | 1 | 1 | 26.34/5.44/60.0 | 60.0 | 703.4840662543083(268.7196581134482/9995.415752865536;10/1/1) | 530235/119805 | 69 | 936 | 40.1 | 442 |
| pilot | 50kb | 4 | 3 | 1 | 1 | 0 | 0 | 25.50/8.19/60.0 | 60.0 | 965.8172635565556(481.5632394038305/965.8172635565556;4/0/0) | 358663/96580 | 67 | 640 | 19.1 | 433 |
| L1 | 10kb | 25 | 25 | 20 | 20 | 2 | 2 | 60.0/33.98/60.0 | 60.0 | 180.50828630293756(67.80911046598506/796.0668948025335;20/2/3) | 415147/60726 | - | - | - | - |
| L1 | 20kb | 12 | 12 | 9 | 9 | 2 | 0 | 60.0/18.87/60.0 | 60.0 | 1358.968430384015(334.67556263325804/1992.783990735497;10/0/2) | 517154/75855 | 40 | 229 | 18.5 | 453 |
| L1 | 50kb | 5 | 5 | 4 | 4 | 0 | 0 | 60.0/6.16/60.0 | 60.0 | 3949.309349862773(3301.325801476777/7193.769829340598;5/0/0) | 433831/96580 | 65 | 276 | 16.2 | 460 |
| L2 | 10kb | 25 | 25 | 22 | 22 | 1 | 1 | 60.0/43.35/60.0 | 60.0 | 317.80713703280605(121.22719995910502/652.5769900207001;21/0/4) | 277007/28920 | - | - | - | - |
| L2 | 20kb | 12 | 11 | 10 | 10 | 0 | 0 | 60.0/60.0/60.0 | 60.0 | 183.66317426943203(29.358999797516407/740.3599784585583;7/1/4) | 389178/61075 | 39 | 259 | 23.2 | 484 |
| L2 | 50kb | 5 | 5 | 3 | 3 | 0 | 0 | 60.0/6.81/60.0 | 60.0 | 785.2041073450605(382.6292371236128/1100.142612531535;4/1/0) | 349338/92665 | 58 | 344 | 21.2 | 474 |
| L3 | 10kb | 24 | 24 | 18 | 18 | 4 | 2 | 60.0/33.99/60.0 | 60.0 | 249.19292598836032(90.0658458703721/452.8967043726997;22/2/0) | 181049/29890 | - | - | - | - |
| L3 | 20kb | 12 | 12 | 9 | 9 | 3 | 0 | 60.0/28.69/60.0 | 60.0 | 283.30388789651516(164.3721080938579/518.0491327954085;10/2/0) | 253307/45753 | 34 | 194 | 12.5 | 453 |
| L3 | 50kb | 4 | 4 | 3 | 3 | 1 | 0 | 60.0/50.00/60.0 | 60.0 | 50912.81792621389(834.20063776453/50912.81792621389;4/0/0) | 320486/89253 | 52 | 219 | 7.3 | 468 |
| L4 | 10kb | 25 | 25 | 17 | 17 | 3 | 3 | 60.0/40.00/60.0 | 60.0 | 176.14758758305183(8.183451172447494/260.1072443492265;15/0/10) | 279693/62835 | - | - | - | - |
| L4 | 20kb | 12 | 12 | 7 | 7 | 4 | 1 | 60.0/43.02/60.0 | 60.0 | 170.91017989396863(9.277811527918097/218.84989129836032;7/0/5) | 354008/70125 | 42 | 288 | 29.3 | 408 |
| L4 | 50kb | 5 | 5 | 3 | 3 | 2 | 0 | 60.0/46.99/60.0 | 60.0 | 1734.489699316841(100.50422983173442/58248.76985825295;4/0/1) | 387851/100576 | 54 | 414 | 23.7 | 427 |

Every rank-1 call at every grain and depth is PERFECT (zero-edit, QV 60) — the
sweep's rank-1s are exactly the perfect calls (the committed instrument's class
identity, the red-team-audited `--truth-qv` construction).

### DEPTH: 35X (corrected cuts, IMPG_CUT_LOCAL_DERIVE=1)

| locality | grain | loci | expr | rank1 | perfect | rank2 | neartwin-r2 (34-43) | QV med/p10/p90 | rank1-QV med | QUAL med(p25/p75;>=60/34-60/<34) | classes sum/max | cli_wall_s | census_wall_s | census_RSS_GiB | cut_wall_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pilot | 10kb | 22 | 22 | 1 | 1 | 3 | 3 | 28.09/2.44/40.02 | 60.0 | 55.97173842258253(4.4509147141124386/700.4409155814743;11/1/10) | 768252/204480 | - | - | - | - |
| pilot | 20kb | 12 | 11 | 1 | 1 | 0 | 0 | 25.35/5.44/43.02 | 60.0 | 386.81881492093214(34.088273239765066/681.7932114007174;8/1/3) | 542444/119316 | 99 | 936 | 40.1 | 442 |
| pilot | 50kb | 4 | 3 | 1 | 1 | 0 | 0 | 13.74/9.97/60.0 | 60.0 | 402.8256688457699(267.30604273862434/402.8256688457699;4/0/0) | 358235/96580 | 125 | 640 | 19.1 | 433 |
| L1 | 10kb | 25 | 25 | 12 | 12 | 6 | 6 | 43.03/17.24/60.0 | 60.0 | 14.782666905822328(3.681431238773913/82.65354400976987;7/3/15) | 409934/60726 | - | - | - | - |
| L1 | 20kb | 12 | 12 | 8 | 8 | 0 | 0 | 60.0/16.77/60.0 | 60.0 | 72.94640902331435(29.38449727476882/220.9543832322099;6/1/5) | 515088/75855 | 41 | 229 | 18.5 | 453 |
| L1 | 50kb | 5 | 5 | 1 | 1 | 0 | 0 | 22.29/6.16/60.0 | 60.0 | 359.8916438968936(137.72746182246232/542.8090030280175;4/0/1) | 433046/96580 | 75 | 276 | 16.2 | 460 |
| L2 | 10kb | 25 | 25 | 15 | 15 | 3 | 3 | 60.0/15.37/60.0 | 60.0 | 24.26596486912991(6.730036947353809/112.32136841400968;7/3/15) | 274358/28920 | - | - | - | - |
| L2 | 20kb | 12 | 11 | 4 | 4 | 2 | 1 | 46.02/18.26/60.0 | 60.0 | 78.2065418090888(38.757484524405406/374.6487034099961;7/2/3) | 385621/60726 | 26 | 259 | 23.2 | 484 |
| L2 | 50kb | 5 | 5 | 2 | 2 | 0 | 0 | 21.72/6.81/60.0 | 60.0 | 66.06566703174876(64.36337076465409/341.72947692736943;4/0/1) | 347681/92665 | 47 | 344 | 21.2 | 474 |
| L3 | 10kb | 24 | 24 | 10 | 10 | 7 | 3 | 43.02/24.54/60.0 | 60.0 | 16.816920897613027(6.499184825754122/39.00742741784913;5/2/17) | 176224/29890 | - | - | - | - |
| L3 | 20kb | 12 | 12 | 5 | 5 | 2 | 0 | 47.77/21.24/60.0 | 60.0 | 35.30766429802632(24.780315144786925/46.81723593617627;3/3/6) | 247786/45753 | 51 | 194 | 12.5 | 453 |
| L3 | 50kb | 4 | 4 | 2 | 2 | 0 | 0 | 60.0/35.68/60.0 | 60.0 | 520.7603903304403(337.75862448998726/520.7603903304403;3/0/1) | 315288/87990 | 71 | 219 | 7.3 | 468 |
| L4 | 10kb | 25 | 25 | 11 | 11 | 6 | 5 | 43.02/36.02/60.0 | 60.0 | 8.479584812745413(4.78051317136137/60.54874248661784;7/0/18) | 278046/62481 | - | - | - | - |
| L4 | 20kb | 12 | 12 | 4 | 4 | 2 | 0 | 46.02/40.00/60.0 | 60.0 | 6.272638014838547(3.1736916132862736/13.571770750054046;2/0/10) | 352713/69751 | 35 | 288 | 29.3 | 408 |
| L4 | 50kb | 5 | 5 | 3 | 3 | 1 | 0 | 60.0/29.28/60.0 | 60.0 | 103.51159572088534(13.674026910694872/3692.121094666208;3/0/2) | 387424/100576 | 38 | 414 | 23.7 | 427 |

## THE ACCURACY-VS-GRAIN CURVE (five localities pooled; window count and bp-weighted)

| grain | DEEP rank1/loci | DEEP bp-weighted | 35X rank1/loci | 35X bp-weighted | near-twin r2 deep / 35x |
|---|---|---|---|---|---|
| 10kb | 81/121 = 66.9% | 821,134/1,252,016 = 65.6% | 49/121 = 40.5% | 490,848/1,252,016 = 39.2% | 11 / 20 |
| 20kb | 37/60 = 61.7% | 760,126/1,252,016 = 60.7% | 22/60 = 36.7% | 440,103/1,252,016 = 35.2% | 2 / 1 |
| 50kb | 14/23 = 60.9% | 750,101/1,252,016 = 59.9% | 9/23 = 39.1% | 450,115/1,252,016 = 36.0% | 0 / 0 |

(The bp-weighted view counts each window's axis-row bp; the five territories total
1,252,016 bp of axis sequence at every grain — the denominators are identical, the
window mix is all that changes.)

## The near-twin rank-2 class's fate (the 10kb QV34-43 windows traced into their enclosing coarser windows)

The class NEVER entrenches as a decisive rank-2 anywhere. At 50kb it leaves the
34-43 band entirely (0 windows at both depths), by two mechanisms:

- **ABSORBED into rank-1** (the coarser window's clean flanking material drowns the
  twin): DEEP — L1's two near-twins (rank2/QV43 -> rank1/QV60), L3's two (one of
  which had worsened to QV27.6 at 20kb), pilot's [30010,40047); 35X — L1's
  [50058,60059) and L2's two. Absorption is dilution-by-clean-flanks, not twin
  separation: it works exactly where the surrounding windows are already rank-1.
- **DEGRADED to rival-dominant** (the coarser window swallows the repeat/segdup
  domain and a rival wins outright): DEEP — L2's [240047,250862) rank2/QV43.4 ->
  rank332 at 20kb -> rank112 at 50kb (monotone worsening with grain); pilot's
  [50013,80029) cluster -> rank 8667 at 50kb, and the pilot segdup-block window
  -> rank 21705 at 50kb. 35X — L1's four near-twins -> rank 600/769/870/16372 at
  50kb. Degradation is the segdup-swallowing cost the sweep's design named as a
  watch item: it MATERIALIZED, at every coarser grain, at the repeat-adjacent
  localities.
- A residue at L4-50kb sits at rank-2 with sequence QV 50/47 — outside the 34-43
  band only because the called near-twin pair shares MORE bp with truth in the
  bigger window (the QV of a rank-2 call rises mechanically with window size);
  the rank and the LL gap are the near-twin class, unconverted.

## The per-bp invariance test (is accuracy grain-invariant with only the window mix changing?)

NOT as a law — as a locality property. Where the underlying window verdicts are
stable, the bp-weighted accuracy is grain-invariant to within a window's worth of bp
(L1-DEEP 80.1% -> 71.7% -> 79.7% — the named "L1-w50k 80% = 10kb's 80%" signal is
real; L3-DEEP 76.0% -> 76.0% -> 80.0%; L3-35X 40.0% at ALL THREE grains; pilot-DEEP
16.0% -> 16.0% -> 20.0%). Where the coarser window's domain swallows repeat material,
the bp view moves as much as the window view (L2-DEEP 87.7% -> 79.7% -> 59.8%;
L1-35X 47.8% -> 63.7% -> 19.9%; L4-35X 44.0% -> 32.0% -> 60.0%). POOLED, the bp view
declines with grain at both depths (65.6% -> 60.7% -> 59.9% deep; 39.2% -> 35.2% ->
36.0% at 35x): the window mix is not neutral — coarse windows trade absorbed
near-twins for swallowed segdup domains, and the trade is net-negative.

## The honest costs, measured

- **Segdup swallowing at coarse grains** (the design's named watch item, materialized):
  five L4 10kb rank-1 PERFECT windows fall to rank-2 at 20kb; L1's last two 10kb
  windows (rank-1 at 10kb) fall to rank 5,907 at 20kb; L2's near-twin window degrades
  monotonically (rank 2 -> 332 -> 112); pilot's segdup block goes deep rival-dominant
  (rank 8667 / 21705 at 50kb).
- **Expressibility losses** (the coarser fold universes drop truth pairs): the pilot's
  segdup-block window is INEXPRESSIBLE at 20kb AND 50kb at BOTH depths (expressible
  but rival-dominant at 10kb); L2 loses one window at 20kb. Every other window keeps
  the truth pair in the universe.
- **The census/monster growth: REFUTED by measurement** — the census peak SHRINKS
  with window count (pilot 40.1 GiB at 20kb -> 19.1 GiB at 50kb; L3 12.5 -> 7.3;
  L1 18.5 -> 16.2). The exhaustive-class domain grows per window (max classes
  L1 60,726 -> 75,855 -> 96,580; L2 28,920 -> 61,075 -> 92,665; sub-linear in window
  size) but the SUM stays flat or shrinks (fewer windows): no universe blow-up.
- **The chain-boundary reduction: REAL but accuracy-negative** — the DP spine's
  windows drop 22-25 -> 12 -> 4-5 per locality (boundaries 24 -> 11 -> 3-4); the
  census wall is flat-to-slightly-up (229-414s), and the accuracy pays for it at
  every coarser point on both pooled curves.
- **The in-window switch check**: HG002's truth haplotypes stay contiguous HPRC
  paths at every grain (one fragmented-truth window at pilot-20kb, zero everywhere
  else) — the cross-window frame class does not arise.

## THE ENGINE-BUG RECEIPT: the FastGA sliver-partition heap corruption (new, real, deterministic; NO fix implemented without the owner's gate)

The coarse tilings emit degenerate non-axis SLIVER partitions (all member rows a few
bp: L3-w20k partition17 = 48 rows of 8bp; L4-w20k partition12 = 14 rows of 11bp).
Such a partition DETERMINISTICALLY CRASHES the pggb engine's FastGA index builder —
the 10kb tilings never emitted the class at these localities. Reproduced from the
preserved L4-w20k partition12 members (grain-sweep/crash-probe2/): the committed
release `GIXmake` SIGSEGVs (or, per heap layout, glibc aborts "corrupted
double-linked list (not small)" — the original receipts; the predecessor's -T1 probe
failed as "Rewind of file .post.2661645.0.idx failed"). ROOT CAUSE, line-precise
(pinned engine rev pangenome/fastga-rs 5216a157f761545c9c3160ea1501ef1beab9c663,
`deps/fastga/GIXmake.c`, the copy impg's build.rs pulls): with every member sequence
shorter than the index k-mer (KMER=40), the syncmer scan samples ZERO posts; then
`distribute()`'s part-split training loop (GIXmake.c:675-691) runs on an all-zero
bucket distribution — `t = buck[NUM_BUCK-1]/NPARTS = 0`, every one of the 1024
buckets satisfies `buck[i] >= t`, `n` increments per bucket, and `Ksplit[n]` is
written for n = 1..1024 into the (NPARTS+1)-int allocation made at GIXmake.c:1918
(9 ints at -T8): a ~4KB out-of-bounds write train that clobbers the adjacent heap
(measured post-mortem: `Buckets[0] = 0x1d0000001c`, `Select[i] >= 15` at NPARTS=8).
ASan names the first OOB access at GIXmake.c:677 (the i=0 iteration reads
`t - buck[i-1]` = `buck[-1]`, 8 bytes before the Buckets array). THE FIX, stated for
the owner's gate: in `distribute()`, guard the zero-total case before the training
loop — if `buck[NUM_BUCK-1] == 0`, set the uniform split directly (`Ksplit[p] =
p*NUM_BUCK/NPARTS` for p <= NPARTS, `Select[i] = min(i*NPARTS/NUM_BUCK, NPARTS-1)`)
and skip the loop (or bail with a named empty-index error); the loop's invariant
n <= NPARTS is guaranteed only by a non-zero total. Optional impg-side hardening:
the `graph` command could pre-reject members sets with no sequence >= KMER with a
named receipt instead of an abort. The sweep's receipt-side containment (the
axis-only, provably-consumed build set) already avoids the class; `short_GDB_fix`
(GIXmake.c:1605-1620) does not cover it (it pads too-FEW contigs, not too-SHORT
ones, and its fake contigs are scan-skipped anyway).

## THE FINAL VERDICT on the owner's grain hypothesis

**NOT SUPPORTED at either depth; the recommended grain derived from the curve is
the committed 10kb.** The pooled accuracy-vs-grain curves are monotonically
non-improving: DEEP 66.9% -> 61.7% -> 60.9% (window) and 65.6% -> 60.7% -> 59.9%
(bp); 35X 40.5% -> 36.7% -> 39.1% (window) and 39.2% -> 35.2% -> 36.0% (bp) — the
50kb@35x "more discriminating mass" prediction does not rescue the 35x operating
point (its partial recovery over 20kb still lands below 10kb). The near-twin
QV34-43 rank-2 class never becomes decisive: it is absorbed by clean flanks or
degraded by segdup swallowing, never separated. Every rank-1 call at every grain is
a PERFECT zero-edit call, so the accuracy differences are pure window-mix effects
— and the measured window mix says the coarser trades lose more (segdup domains
swallowed, expressibility lost at the pilot segdup block) than they win (near-twins
absorbed where the flanks were already clean). The coarser grain's real measured
benefits — chain boundaries 24 -> 11 -> 3-4, census peaks halved (the feared
monster growth refuted), per-window domains sub-linear — do not offset the
accuracy paid. The 10kb grain stands as the accuracy-optimal point of the sweep;
the grain axis is measured to exhaustion at these five localities.

## Standing rules (held throughout)

Assessment-side truth only via the test mode; no thresholds; no tuning constants
(the grains are the owner's specified points, the w:d ratio the committed recipe's
own); scoreboard machinery unmodified; zero src/ changes (the sweep is receipt-side:
the only repo-adjacent artifacts are the out-of-repo grain-sweep scripts and this
docs section); the 64GiB guard clean everywhere (census peaks 7.3-40.1 GiB); disk
floor discipline (the 18G v2 floor; the disk-exhaustion casualties of Oct-9 all
recovered by idempotent relaunch); NO PR PUSH (the owner pushes).
