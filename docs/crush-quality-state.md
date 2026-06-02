# Crush Quality State

Tracks the iteration-by-iteration state of the graph-quality improvement cycle.
The goal is to close the gap between the reference graph (pggb-in-impg) and the
target graph (syng→sweepga+seqwish-k311→crush auto mode) on the canonical C4 locus.

---

## Inputs

**Canonical C4 input locus**

```
Region:     GRCh38#0#chr6:31891045-32123783
Depth:      -d 50k
Syng index: /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng
Sequences:  /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc
```

This range was used to produce the known-good artifact in
`data/c4_crush_eval_20260523T140141Z/` (see `docs/c4-crush-handoff.md`).
Use it as the canonical input for all iterations.

A pre-extracted blunt graph used as a fast local proxy for testing:
```
q.gfa    — untracked scratch in the repo root; regenerate with the target
           pipeline command below if not present
```

---

## Reference Pipeline

The reference graph is produced by the pggb engine built into `impg query`.
It runs sweepga → seqwish → gfaffix smoothing but does **not** apply the syng
native extraction or crush post-processing.

```bash
out=data/quality_ref_$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$out"

/usr/bin/time -v target/release/impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -x \
  -o gfa \
  --force-large-region \
  -O "$out/ref.pggb.nosort" \
  -v 1 \
  > "$out/ref.pggb.nosort.stdout" \
  2> "$out/ref.pggb.nosort.stderr"

gfasort \
  -i "$out/ref.pggb.nosort.gfa" \
  -o "$out/ref.pggb.Ygs.gfa" \
  -p Ygs \
  -t 32

target/release/impg graph-report \
  -g "$out/ref.pggb.Ygs.gfa" \
  -o "$out/ref.pggb.Ygs.report.md" \
  --povu \
  -t 32
```

Default `-o gfa` uses `--gfa-engine pggb` (sweepga + seqwish; `impg query --help`
lists `pggb` as the default engine).

**Known-good reference baseline** (from `data/c4_crush_eval_20260523T140141Z/`):
The best prior crush run used the same C4 locus and produced:

| Metric | Value |
|--------|-------|
| segments | 18089 |
| links | 21018 |
| paths | 465 |
| path steps | 4590745 |
| segment bp | 422007 |
| path-depth median | 144 |
| path-depth p95 | 620 |
| link jump p99 | 45 |
| ws bridges ≥ 1000 bp | 25414 |
| ws total bp | 449787187 |
| ws p95 bp | 125 |
| ws p99 bp | 393 |
| ws max bp | 114697 |

This artifact predates code changes in `98fd538` and represents the quality
floor the target pipeline should meet or exceed. A true pggb reference run is
needed to fill the `quality_ref_*` section above; this baseline is a proxy.

---

## Target Pipeline

The target graph is produced by the syng native extraction path with mask,
followed by the sweepga+fastga+seqwish-k311 crush pipeline in auto mode.

```bash
out=data/quality_target_$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$out"

/usr/bin/time -v target/release/impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=sweepga,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/target.nosort" \
  -v 1 \
  > "$out/target.nosort.stdout" \
  2> "$out/target.nosort.stderr"

gfasort \
  -i "$out/target.nosort.gfa" \
  -o "$out/target.Ygs.gfa" \
  -p Ygs \
  -t 32

target/release/impg graph-report \
  -g "$out/target.Ygs.gfa" \
  -o "$out/target.Ygs.report.md" \
  --povu \
  -t 32
```

---

## Metrics

Metrics extracted from `impg graph-report --povu` output on Ygs-sorted GFA.
The gap metric is `(target - reference) / |reference|` where positive means
target is worse.

Tracked per iteration:

| Metric | Reference | Target | Gap | Notes |
|--------|-----------|--------|-----|-------|
| segment bp | 422007 | TBD | TBD | lower → more compact |
| ws total bp | 449787187 | TBD | TBD | lower → less whitespace |
| ws p99 bp | 393 | TBD | TBD | primary quality signal |
| ws max bp | 114697 | TBD | TBD | worst-case bridge |
| ws bridges ≥ 1000 bp | 25414 | TBD | TBD | count of long bridges |
| link jump p99 | 45 | TBD | TBD | layout adjacency |
| paths | 465 | TBD | TBD | must be preserved |
| POVU sites | 2443 | TBD | TBD | bubble resolution coverage |

Reference values are from the known-good baseline (see Reference Pipeline section).
Target values will be populated by quality-iter-measure each iteration.

**Proxy metric (q.gfa — fast test)**

`q.gfa` is a pre-extracted blunt GFA (see Inputs) used for fast iteration
without re-running the full pipeline. Metrics on this proxy are not directly
comparable to the full run but track relative improvement across iterations.

---

## Iteration Log

### Iter 0 → 1 (2026-05-24)

**What was fixed:**
`quality-iter-fix` (agent-34) analyzed whitespace bridge metrics on q.gfa and
found that replacement segments emitted by `render_rewritten_graph`
(`src/resolution.rs`) were being placed by first-encounter path order rather
than near their predecessor original nodes. When a path traversing a replacement
appeared late (after many other paths had added original segments), the
replacement landed far from its predecessor in S-line order, creating large
whitespace bridges that dominated the quality score.

The fix (commit `6cb4f22`): two-phase placement algorithm that first builds
original-node order, then stable-inserts each replacement immediately after the
position of its last original predecessor in the same path.

**Did it close the gap?**
Partially. On q.gfa (99 bubbles resolved):
- ws-total: 52.7 MB → 42.6 MB (−19%)
- ws-p99: 8219 bp → 7512 bp (−8.6%)

**Residual gap going into iter 2:**
ws-p99 at 7512 bp vs known-good reference p99 of 393 bp — still ~19× above
reference. ws-total at 42.6 MB (proxy). The fix addresses one root cause but
larger structural whitespace bridges remain (likely from the original syng blunt
graph and unreplaced bubbles). **Not converged.**

**Proposed next focus** (to be confirmed by quality-iter-measure/compare):
Measure the target pipeline on the full C4 locus with current HEAD to get
absolute metrics. Identify the next largest gap: likely ws-p99 from unreplaced
large bubbles or long-range link jumps from the original syng graph structure.
