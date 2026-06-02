# crush + gfaffix workflow — slice run + recommended pipeline

**Task:** `crush-add-gfaffix-2` (opus retry). The prior sonnet attempt produced
zero code or run evidence; this run records actual numbers and PNG URLs.

**Date:** 2026-05-24
**Worktree:** `wg/agent-79/crush-add-gfaffix-2`
**Author:** agent-79 (Programmer role, opus tier)

## Option chosen

**Option A — external post-process workflow.**

Rationale:

- `gfaffix` is already vendored in this repo (`vendor/gfaffix/`) and is built
  into a top-level binary via the second `[[bin]]` entry in `Cargo.toml`. The
  binary lives at `target/release/gfaffix` after `cargo build --release`.
- `docs/crush-spec-audit.md` §Phase-8 is explicit that the post-emit gfaffix
  pass is a **band-aid** that "should not substitute for the Phase 3/4 fix".
  Adding an in-pipeline `:gfaffix` operator at this point hides the upstream
  duplicate-sequence emission bug rather than fixing it.
- Option A is reversible (users opt in by piping through the binary); Option B
  changes the default behaviour of every `gfa:syng:...:nosort` consumer.

The recommended workflow is documented in [§Recommended workflow](#recommended-workflow)
below. A follow-up task can promote it to a `:gfaffix` operator after the
Phase 3/4 fixes land.

## Slice run — verbatim commands and results

Slice fixture: `tests/test_data/crush/c4_slice_1500_3000.gfa` (7,169,407 bytes,
2,942 segments, 465 paths, 64,348 bp). Output dir: `/tmp/crush-gfaffix-run/`.

### Step 1 — impg crush on the slice

```bash
mkdir -p /tmp/crush-gfaffix-run
cd /tmp/crush-gfaffix-run

/usr/bin/time -v /home/erikg/impg/target/release/impg crush \
  -g /home/erikg/impg/tests/test_data/crush/c4_slice_1500_3000.gfa \
  --method auto --max-iterations 3 \
  > slice.raw.gfa 2> slice.stderr
```

Result:

- exit 0, wall **20.83 s**, max RSS 732 MB
- `slice.raw.gfa` = 7,169,407 bytes, 6727 lines, **byte-identical to input**
- Reason from stderr: `crush round 1: rejecting 33/33 replacement(s); score grew 15.3% (threshold 10.0%)`
  — every candidate replacement was rejected by the regression-control gate
  in `apply_replacement_frontier`, so the working graph passes through
  unchanged. The "raw" file we feed to gfaffix is therefore the post-syng
  blunt graph in its as-extracted form.

### Step 2 — gfaffix on the crush output

```bash
cd /tmp/crush-gfaffix-run

/usr/bin/time -v /home/erikg/impg/target/release/gfaffix \
  slice.raw.gfa -o slice.gfaffix.gfa -c \
  > gfaffix.log 2> gfaffix.stderr
```

The `-c` flag enables gfaffix's built-in correctness check
("verifies that the transformed parts of the graphs spell out the identical
sequence as in the original graph") and the run logs `all correct!`.

Result:

- exit 0, wall **0.12 s**, max RSS 18 MB
- gfaffix reported: `found and collapsed 775 shared prefixes, 458 of which are overlapping, and 182 of which are bubbles`
- `slice.gfaffix.gfa` = 8,979,672 bytes, 6981 lines

### Step 3 — metrics

Computed by `/tmp/crush-gfaffix-run/metrics.py` (a small Python pass over both
GFAs that counts unique canonical-sequence segments, splits dup-segments by
size bucket, and walks every `P`-line through the segment table to compare
path-spelled sequences byte-for-byte).

| metric                                       | slice.raw.gfa | slice.gfaffix.gfa | delta            |
|---|---:|---:|---:|
| total segments                               | 2,942         | 3,070             | +128 (+4.4%)     |
| total segment bp                             | 64,348        | 53,653            | −10,695 (−16.6%) |
| duplicate-canonical-sequence segments        | 974           | 932               | −42              |
| duplicate-segment fraction                   | 33.1%         | 30.4%             | −2.7 pp          |
| dup-seg size buckets — ≤8 bp / 9-32 / 33-128 | 558 / 318 / 98 | 811 / 95 / 26    | small bucket grows, larger buckets collapse 70-73% |
| paths preserved                              | 465           | 465               | 0                |
| sum of path-spelled bp                       | 15,221,205    | 15,221,205        | 0                |
| path-sequence mismatches                     | 0             | 0                 | byte-identical   |

### Interpretation

- **gfaffix reduced segment-bp by 16.6%** (64,348 → 53,653) — this is the
  payload that downstream genotyping cares about, because each duplicate byte
  is a position where read mapping has to pick between equivalent landing
  sites.
- **Path-sequence preservation is exact** for all 465 paths (15.2 Mbp
  walked, 0 mismatches). gfaffix's `-c` correctness check passed; our
  independent byte-equality pass agrees.
- **Segment count rose +4.4%** because gfaffix splits segments to factor out
  shared affixes — the new short segments are visible in the dup-size buckets
  (≤8 bp: 558 → 811). This is expected gfaffix behaviour and is what makes
  the bp shrink possible.
- **Duplicate-segment count is not zero (932 ≠ 0).** The task description
  asks to either reach 0 or document the residual. The residual is dominated
  by **≤8 bp segments** (811 of 932 = 87% of residual dups). With only
  `4^8 = 65,536` possible canonical 8-mers, short segments collide on
  canonical sequence by chance even when their graph contexts are completely
  different — these are *not* gfaffix-collapsible affixes (their walks are
  not shared) and removing them would mis-merge unrelated paths. The
  9-32 bp bucket shrank by 70% (318 → 95) and the 33-128 bp bucket by 73%
  (98 → 26); the >128 bp bucket was already empty in the input. That is
  the expected gfaffix outcome on a well-formed graph: collapse what shares
  walks, leave random k-mer collisions alone.

### Visualizations

Both PNGs were rendered with `gfalook` (`-m` = colour by mean coverage,
2200x1200):

```bash
gfalook -i slice.raw.gfa -o slice.raw.png -m -x 2200 -y 1200
gfalook -i slice.gfaffix.gfa -o slice.gfaffix.png -m -x 2200 -y 1200
```

Uploaded to:

- https://hypervolu.me/~erik/impg/c4_slice_1500_3000.crush.raw.png (1,234,162 bytes)
- https://hypervolu.me/~erik/impg/c4_slice_1500_3000.crush.gfaffix.png (3,523,996 bytes)

The gfaffix render is materially larger (more segment edges to draw) and
shows a tighter, less redundant graph — consistent with the −16.6% bp number.

## Full C4 GRCh38 — attempted with 15-min wall timeout

Canonical command from `docs/c4-crush-handoff.md`:

```bash
timeout 900 /usr/bin/time -v /home/erikg/impg/target/release/impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=sweepga,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O /tmp/crush-gfaffix-run/c4full/C4A.parent5k.crush_full.nosort \
  -v 1
```

Result (15-minute wall timeout enforced via `timeout 900`): **completed
successfully within the timeout, exit 0**.

- wall **12:00.99**, max RSS **65.3 GB**, user time 9,321 s on 32 threads
- `C4A.parent5k.crush_full.nosort.gfa` = 47,182,519 bytes
- last crush stderr lines:
  ```
  [impg::resolution] crush round 1: replacement build progress 1/1 (accepted)
  [impg::resolution] crush round 1: rejecting 1/1 replacement(s); score grew 84.2% (threshold 10.0%); ws-p99 177799 -> 189428, ws-max 388845 -> 464223; total 440.69s
  [impg] crush: 0 resolved, 1 bailed, 1 candidates seen across 1 rounds
  [impg] Syng query complete: GRCh38#0#chr6:31891045-32123783 in 8m44.438s
  ```
- The known-good baseline in `docs/c4-crush-handoff.md` (line 64) is wall
  **4:39.95** at max RSS **55.1 GB**; this run is **2.6× slower** and
  **1.18× heavier** at peak. The crush bailed at the regression-control
  gate (84.2% score growth on a single candidate), which is the same
  symptom — different surface — as the >30 min hang the audit describes
  in "Failure 1 + cascade". It completes here because the gate rejects the
  candidate rather than spinning further, but it makes zero crush progress.
- Built artifact = the post-syng blunt graph in raw form (same situation as
  the slice case, see Step 1 above). gfaffix on it does real work because
  the upstream graph has 41.2% duplicate-sequence segments out of the box.

#### gfaffix on the full C4 output

```bash
cd /tmp/crush-gfaffix-run/c4full

/usr/bin/time -v /home/erikg/impg/target/release/gfaffix \
  C4A.parent5k.crush_full.nosort.gfa \
  -o C4A.parent5k.crush_full.gfaffix.gfa -c \
  > c4full.gfaffix.log 2> c4full.gfaffix.stderr
```

Result:

- exit 0, wall **0.95 s**, max RSS 98 MB
- gfaffix reported: `found and collapsed 5540 shared prefixes, 3307 of which are overlapping, and 1526 of which are bubbles`
- `C4A.parent5k.crush_full.gfaffix.gfa` = 60,957,026 bytes (file grew because
  paths now reference more, shorter segments — segment-bp shrank, the
  serialised P-lines grew)

| metric                                       | crush_full.nosort | crush_full.gfaffix | delta            |
|---|---:|---:|---:|
| total segments                               | 18,048            | 18,652             | +604 (+3.3%)     |
| total segment bp                             | 389,354           | 307,751            | −81,603 (−21.0%) |
| duplicate-canonical-sequence segments        | 7,430             | 6,972              | −458             |
| duplicate-segment fraction                   | 41.2%             | 37.4%              | −3.8 pp          |
| dup-seg size buckets — ≤8 bp / 9-32 / 33-128 | 3,910 / 2,690 / 830 | 5,989 / 777 / 206 | 9-32 bp shrank 71%, 33-128 bp shrank 75% |
| paths preserved                              | 465               | 465                | 0                |
| sum of path-spelled bp                       | 101,633,470       | 101,633,470        | 0                |
| path-sequence mismatches                     | 0                 | 0                  | byte-identical   |

Same pattern as the slice: the long-segment duplicates collapse, the residual
is dominated by short k-mer collisions (≤8 bp = 86% of the gfaffix residual),
all 465 paths preserved byte-for-byte.

PNGs uploaded:

- https://hypervolu.me/~erik/impg/C4A.parent5k.crush_full.raw.png (1,229,084 bytes)
- https://hypervolu.me/~erik/impg/C4A.parent5k.crush_full.gfaffix.png (2,599,369 bytes)

## Recommended workflow

Until Phases 3/4 of `docs/crush-spec-audit.md` are fixed and the
post-emit dedup is provably a no-op, downstream genotyping pipelines should
run gfaffix as an explicit second step:

```bash
# 1. Extract + crush the locus (any :nosort variant).
impg query ... \
  -o 'gfa:syng:mask,...:crush,...:nosort' \
  -O work.nosort

# 2. Collapse walk-preserving affixes.
gfaffix work.nosort.gfa -o work.nosort.gfaffix.gfa -c

# 3. Downstream sort / index / genotyping consumes work.nosort.gfaffix.gfa.
```

The `gfaffix` binary is the same one shipped in this repo (built from
`vendor/gfaffix/` via the `[[bin]] name = "gfaffix"` entry in
`Cargo.toml:9-11`), so `cargo build --release` produces both binaries
side-by-side under `target/release/`.

## Hard validation gate status

| gate | status |
|---|---|
| slice.raw.gfa exists at /tmp/crush-gfaffix-run/ | yes — 7,169,407 bytes |
| slice.gfaffix.gfa exists at /tmp/crush-gfaffix-run/ | yes — 8,979,672 bytes |
| wall time recorded | crush 20.83 s, gfaffix 0.12 s |
| duplicate-sequence segment count = 0 in gfaffix output | **932**, documented as ≤8 bp k-mer collisions (87% of residual); 9-32 bp bucket shrank 70%, 33-128 bp bucket shrank 73% |
| path sequences in raw == gfaffix (byte equality) | 465/465 paths match, 15,221,205 bp on each side |
| both PNGs uploaded to hypervolu.me with confirmed URLs | yes, see Visualizations section |
| docs/crush-gfaffix-run.md committed | this file |
| Option A or B documented | Option A; rationale at top |
| full C4 GRCh38 attempted with 15-min wall timeout | **succeeded** in 12:00.99 wall / 65.3 GB RSS; crush bailed (84.2% score growth), gfaffix collapsed 5,540 prefixes in 0.95 s; segment-bp −21.0%, dup-seg fraction 41.2% → 37.4%, 465/465 paths byte-identical; PNGs uploaded |
