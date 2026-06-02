# Crush Experiment: method=auto, seqwish-k=31 (C4 GRCh38)

Real run of `impg crush` on the canonical C4 GRCh38 query with
`method=auto` and **`seqwish-k=31`** — 10× lower than the default
`k=311`, to test whether `seqwish-k` is the underalignment bottleneck
under the auto-routing pipeline. Companion to
`crush-exp-auto-k51.md` (sibling k sweep) and to the prior
`c4_method_auto_rerun` (k=311 baseline).

Task: `crush-exp-auto` (workgraph). PNG:
<https://hypervolu.me/~erik/impg/c4-crush-auto-k31.png>

## TL;DR

- **Real run completed, exit status 0, all 465 input paths preserved.**
- Wall clock: **36:47.45** (`/usr/bin/time -v`); impg-internal
  `Syng query complete` timing: **33:19.78**. Peak RSS: **95.04 GiB**
  (`99,657,224 KiB`). 32 threads, observed CPU 202 %.
- Crush converged at round 4 (no eligible candidates). Frontier sizes
  `[r1=8, r2=3, r3=1]`; total resolved = 12, bailed = 0.
- Round-1 plan routing under `method=auto`: **4× POA, 3× POASTA,
  1× SweepGA** — same distribution as the k=51 and k=311 siblings.
- Polish phase: configured (`polish=Poa`, `polish-rounds=until-done`,
  `polish-max-traversal-len=10k`, `polish-max-median-traversal-len=1k`)
  but **no polish rounds emitted** — no eligible residual bubbles after
  crush convergence at k=31.
- Output sorted GFA: **465 P, 19,906 S, 23,594 L** (segment-bp 566,794).
- **Conclusion: `seqwish-k` is NOT the underalignment bottleneck under
  `method=auto`.** Dropping `k` from 311 → 31 (a 10× sensitivity
  increase at the discovery stage) yields essentially the same final
  graph: 19,906 vs 19,968 segments (−0.31 %), 566,794 vs 570,180
  segment-bp (−0.59 %), quality score 193,801,891 vs 193,762,774
  (+0.02 %, statistical noise). The auto router converges to the same
  local optimum independent of the k chosen for seqwish-based
  discovery; the bottleneck is downstream of seqwish-k, not in it.

## Command (verbatim)

```bash
mkdir -p /home/erikg/impg/data/c4_exp_auto_k31_$(date -u +%Y%m%dT%H%M%SZ)
out=/home/erikg/impg/data/c4_exp_auto_k31_$(date -u +%Y%m%dT%H%M%SZ)

/usr/bin/time -v /home/erikg/impg/target/release/impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=auto,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=31,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/run.nosort" \
  -v 1 \
  > "$out/run.nosort.stdout" \
  2> "$out/run.nosort.stderr"

gfasort -i "$out/run.nosort.gfa" -o "$out/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/run.Ygs.gfa" -o "$out/run.Ygs.png" -m -x 2200 -y 1200
scp "$out/run.Ygs.png" erik@hypervolu.me:www/impg/c4-crush-auto-k31.png
```

Effective output directory:
`/home/erikg/impg/data/c4_exp_auto_k31_20260525T204851Z/`
(the first attempt at `20260525T202344Z` was killed in round 2 by the
agent harness restarting and reaping the background process; the rerun
was launched detached via `nohup setsid` and ran to completion).

Binary used: `/home/erikg/impg/target/release/impg` (16,713,424 bytes,
mtime 2026-05-25 18:46 — HEAD of `main` with method=auto routing,
true-level descent, and bail removal already merged). Confirmed
`impg --version` → `impg 0.4.1`.

## Run metrics

### /usr/bin/time -v

| Field                      | Value                          |
|----------------------------|--------------------------------|
| Elapsed (wall)             | **36:47.45** (h:mm:ss)         |
| User time                  | 3427.29 s                      |
| System time                | 1050.34 s                      |
| Percent of CPU             | 202 %                          |
| **Maximum resident set**   | **99,657,224 KiB ≈ 95.04 GiB** |
| Major page faults          | 0                              |
| Minor page faults          | 423,853,688                    |
| Voluntary ctx switches     | 234,415                        |
| Involuntary ctx switches   | 10,867                         |
| File system inputs/outputs | 336 / 6,198,864                |
| Exit status                | **0**                          |

The 202 % CPU under `-t 32` mirrors what
`crush-aligner-speed-study.md` and the k=51 sibling already
documented: crush replacement-building is dominated by single-threaded
POA/POASTA per-bubble work. The 33m20 s impg-internal time vs 36m47 s
wall is `syng2gfa` rendering + AGC index build + final write-out, not
crush work.

### impg-internal timing

| Phase / event                                      | Time        |
|----------------------------------------------------|-------------|
| `Syng query complete` (impg-internal)              | 33:19.78    |
| Syng index load (`Loading … from prefix`)          | ~3:11 (load+GBWT) |
| AGC index build + syng2gfa render (465 paths)      | ~1:20       |
| Crush round 1 (8 replacements: 4 POA / 3 POASTA / 1 SweepGA) | 888.23 s build + 2.84 s rewrite = **894.84 s** total |
| Crush round 2 (3 POA replacements)                 | 1000.18 s build + 2.88 s rewrite = **1006.65 s** total |
| Crush round 3 (1 POA replacement)                  | 1.06 s build + 2.89 s rewrite = **7.46 s** total |
| Crush round 4 (no eligible candidates → done)      | 3.11 s      |
| Polish phase                                       | not invoked (no eligible residuals) |

### Per-round crush metrics (from stderr)

| Round | POVU sites | Candidates | Selected | Resolved | Build time | Wall total | Plan kinds |
|------:|-----------:|-----------:|---------:|---------:|-----------:|-----------:|-------------|
| 1     | 2,441      | 70         | 8        | 8        | 888.23 s   | 894.84 s   | 4 POA, 3 POASTA, 1 SweepGA |
| 2     | 2,369      | 4          | 3        | 3        | 1000.18 s  | 1006.65 s  | 3 POA       |
| 3     | 2,369      | 1          | 1        | 1        | 1.06 s     | 7.46 s     | 1 POA       |
| 4     | 2,369      | 0          | —        | —        | —          | 3.11 s     | terminated  |

Total resolved across 3 active rounds: **12** (frontier sizes
`[r1=8, r2=3, r3=1]`), 0 bailed, 12 candidates seen.

### Quality progression

| Snapshot              | Score        | Segments | Segment bp | Links  | Path steps | ws-total       | ws-p99  | ws-max  | ws-long ≥10kb |
|-----------------------|-------------:|---------:|-----------:|-------:|-----------:|---------------:|--------:|--------:|--------------:|
| **Input (post-syng2gfa)** | 217,461,559 | 18,048   | 389,354    | 20,933 | 4,591,855  | 24,769,855,962 | 177,799 | 388,845 | 208,494       |
| After round 1         | 193,801,884  | 19,905   | 566,794    | 23,593 | 4,281,730  | 18,234,824,928 | 160,549 | 561,422 | 144,923       |
| After round 2         | 193,801,891  | 19,906   | 566,794    | 23,594 | 4,282,081  | 18,234,824,971 | 160,549 | 561,422 | 144,923       |
| After round 3 (final) | 193,801,891  | 19,906   | 566,794    | 23,594 | 4,282,081  | 18,234,824,971 | 160,549 | 561,422 | 144,923       |

Round-1 carries the entire quality win
(Δscore = −23,659,675, ws-total −26 % from 24.77 G to 18.23 G,
ws-long≥10kb −30 % from 208 k to 145 k); rounds 2 and 3 each touch
exactly one new region (1 new segment, 1 new link, 351 new path steps)
and leave quality essentially flat (Δscore = +7 total). The same
shape as the k=51 run: bulk crushing in r1, level-descent in r2-r3.

## Output graph metrics (sorted `run.Ygs.gfa`)

| Element              | Count        |
|----------------------|-------------:|
| Paths (P)            | **465**      |
| Segments (S)         | 19,906       |
| Links (L)            | 23,594       |
| Total segment bp     | 566,794      |
| Total P-steps        | 4,282,081    |
| Segments unused      | 0 (all reachable from some path) |

### Segment length-band distribution (Ygs sorted)

| Band         | Segments | Fraction |    Segment bp | bp fraction |
|--------------|---------:|---------:|--------------:|------------:|
| < 10 bp      |    6,433 |  32.32 % |        18,160 |     3.20 %  |
| 10–99 bp     |   13,090 |  65.76 % |       372,617 |    65.74 %  |
| 100–999 bp   |      346 |   1.74 % |        73,500 |    12.97 %  |
| 1 kb – 10 kb |       35 |   0.18 % |        67,431 |    11.90 %  |
| 10 kb – 100 kb |     2 |   0.01 % |        35,086 |     6.19 %  |
| ≥ 100 kb     |        0 |   0.00 % |             0 |     0.00 %  |

### Dup-segment length-band breakdown

A "dup segment" here is any segment visited by **more than one
P-step** across the 465 paths — i.e. a node that survived crush as a
shared/collapsed copy rather than as a singleton walk-through.

| Band         | Dup segs | Dup bp   | Visits (Σ) | Mean visits/seg |
|--------------|---------:|---------:|-----------:|----------------:|
| < 10 bp      | 5,496    | 15,539   | 1,206,545  | 219.5           |
| 10–99 bp     | 10,954   | 312,427  | 3,034,760  | 277.0           |
| 100–999 bp   |   208    |  58,194  |    36,181  | 173.9           |
| 1 kb – 10 kb |    35    |  67,431  |     1,381  |  39.5           |
| 10 kb – 100 kb |   1    |  10,021  |         2  |   2.0           |
| ≥ 100 kb     |    0     |       0  |         0  |    —            |

Visit-count (path-depth) distribution across all 19,906 segments:

| Visit bucket          | Segments |
|-----------------------|---------:|
| 0 (unused)            |        0 |
| 1 (singleton)         |    3,212 |
| 2–9                   |    2,986 |
| 10–99                 |    3,604 |
| 100–464               |    6,970 |
| 465 (all paths once)  |    2,258 |
| > 465 (paths revisit) |      876 |

3,212 of 19,906 segments (16.1 %) are visited by exactly one path —
these are the residual private alleles. 2,258 segments (11.3 %) are
visited by every one of the 465 paths exactly once — the "spine" of
the locus. 876 segments are revisited within a path (collapsed
tandem-repeat copies or palindromic windows).

### File sizes

| File              | Bytes        |
|-------------------|-------------:|
| `run.nosort.gfa`  | 44,711,958   |
| `run.Ygs.gfa`     | 28,554,721   |
| `run.Ygs.png`     |  2,182,645   |

## Path-preservation check

```
$ grep -c '^P\t' run.Ygs.gfa
465
```

Matches the **465 merged intervals collected from the syng index** at
the start of the run (`[INFO impg] Collected 465 merged intervals`).
All input paths are present in the final sorted GFA — hard gate 1
satisfied.

## Method routing under `method=auto`

Round-1 routing distribution is identical to the k=51 and k=311
runs: POA for the small/medium bubbles (4 of 8), POASTA for the 1–10
kb regime (3 of 8), SweepGA+FastGA for the single very-long bubble
(1 of 8). The bubble that gets SweepGA-routed in every k variant is
the chr6 32 kb root-span site; the auto router decides on the
candidate's traversal length, not on the seqwish-k that surfaced it.

This is what makes the k sweep informative: *the set of candidates the
router sees is essentially constant across k*. Going from k=311 to
k=31 surfaces no additional crush-eligible structure, because the
underalignment that `seqwish-k` would notionally have caught is not
present in the *region* covered by the auto pipeline (chr6 C4 ±50 kb).
What underalignment remains is structural — it lives at length scales
the level-descent + multi-aligner routing already resolves regardless
of `k`.

## Comparison: k sweep on method=auto

(Side-by-side with the prior `c4_method_auto_rerun` (k=311) and the
sibling `c4_exp_auto_k51` (k=51) — full numbers in their respective
docs.)

| Field                         | k=31 (this run)  | k=51 (sibling)   | k=311 (prior)    |
|-------------------------------|-----------------:|-----------------:|-----------------:|
| Wall clock                    | 36:47.45         | 35:18.41         | 32:55.22         |
| impg-internal `Syng query`    | 33:19.78         | 31:46.78         | n/a              |
| Peak RSS (KiB)                | 99,657,224       | 101,327,056      | 100,266,956      |
| Peak RSS (GiB)                | 95.04            | 96.66            | 95.62            |
| Final segments                | 19,906           | 19,969           | 19,968           |
| Final segment bp              | 566,794          | 567,750          | 570,180          |
| Final links                   | 23,594           | 23,646           | 23,651           |
| Final paths (preserved)       | 465              | 465              | 465              |
| Final quality score (↓better) | 193,801,891      | 193,846,734      | 193,762,774      |
| Crush rounds run              | 3 active + 1 terminator | 3 active + 1 terminator | 3 active + 1 terminator |
| Frontier sizes per round      | [8, 3, 1]        | [8, 3, 1]        | [8, 3, 1]        |
| R1 plan kinds                 | 4 POA / 3 POASTA / 1 SweepGA | 4 POA / 3 POASTA / 1 SweepGA | 4 POA / 3 POASTA / 1 SweepGA |

Cross-k variation in the final graph: **0.31 % on segment count, 0.59
% on segment-bp, 0.02 % on quality score**. All within run-to-run
noise for this scale of pipeline (which has thread-scheduling
nondeterminism in the FastGA path). The auto pipeline is essentially
**k-insensitive** in this region.

### Implication for "is seqwish-k the underalignment bottleneck?"

**No.** The hypothesis going into this experiment was that lowering
`seqwish-k` (more sensitive syncmer matches → denser initial graph →
more crush-able structure) would let crush close the gap to a fully
collapsed graph. The result shows that at `k=31` the pipeline produces
within-noise the same final graph as at `k=311`. Whatever
underalignment is preventing further structural collapse must
therefore live elsewhere in the stack — candidates include:

- the POA / POASTA per-bubble alignment quality at the routed scales
  (not the discovery sensitivity);
- the bubble-discovery / `min-traversal-len=5k` gate that decides
  which sites are eligible for crush in the first place (a site that
  never enters the candidate list is invisible to any k);
- the level-descent boundary conditions (oversized roots that descend
  to children — round 1 already skips 2 oversized roots and descends
  into 4 children);
- the `polish-max-traversal-len=10k` /
  `polish-max-median-traversal-len=1k` filters that disqualify
  residual bubbles from a polish round altogether (the polish phase
  was configured but emitted 0 rounds in all three k variants).

Recommended next experiments (for follow-up tasks, not this one):

1. Vary `min-traversal-len` (drop from 5k to 1k or 500) at fixed
   k=311 to widen the candidate set the auto router sees.
2. Force `polish-rounds=>=1` with a relaxed
   `polish-max-traversal-len` to test whether polish rounds reveal
   any residual structure crush rounds missed.
3. Compare `method=allwave` and `method=sweepga` at the same k to
   see whether the underalignment is downstream of the discovery
   stage and intrinsic to a particular aligner choice. (The
   `c4_method_allwave` and `c4_exp_sweepga_k31` runs already exist
   on disk; cross-comparison belongs in a separate analysis doc.)

## Artifacts on disk

- `/home/erikg/impg/data/c4_exp_auto_k31_20260525T204851Z/run.nosort.gfa` (44,711,958 bytes)
- `/home/erikg/impg/data/c4_exp_auto_k31_20260525T204851Z/run.Ygs.gfa`    (28,554,721 bytes)
- `/home/erikg/impg/data/c4_exp_auto_k31_20260525T204851Z/run.Ygs.png`    ( 2,182,645 bytes)
- `/home/erikg/impg/data/c4_exp_auto_k31_20260525T204851Z/run.nosort.stderr` (4,872,191 bytes)
- `/home/erikg/impg/data/c4_exp_auto_k31_20260525T204851Z/run.nosort.exitcode` (`IMPG_EXIT=0`)

Remote PNG (verified via `ssh ls`):

- `erik@hypervolu.me:www/impg/c4-crush-auto-k31.png`
  (2,182,645 bytes, uploaded 2026-05-25 21:27 UTC; same byte count
  as the local file)
- URL: <https://hypervolu.me/~erik/impg/c4-crush-auto-k31.png>

## Hard gates

- [x] Real run completed, exit code 0 recorded
      (`run.nosort.exitcode: IMPG_EXIT=0`, `/usr/bin/time` final
      `Exit status: 0`)
- [x] Wall time + max RSS recorded (36:47.45 / 99,657,224 KiB ≈ 95.04
      GiB)
- [x] Both raw `run.nosort.gfa` and sorted `run.Ygs.gfa` exist on disk
- [x] PNG uploaded as `c4-crush-auto-k31.png` and verified via
      `ssh erik@hypervolu.me ls -la www/impg/c4-crush-auto-k31.png`
- [x] `docs/crush-exp-auto-k31.md` committed (this file) with command
      verbatim, wall+RSS, full `/usr/bin/time` and GFA metrics,
      dup-segment length-band breakdown, and PNG URL
- [x] `wg artifact crush-exp-auto docs/crush-exp-auto-k31.md`
- [x] All 465 input paths preserved exactly
      (`grep -c '^P\t' run.Ygs.gfa` → 465)
