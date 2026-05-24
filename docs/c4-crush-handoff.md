# C4 Crush Handoff

This note records the current C4 graph-crushing state so another pass can pick
up from the useful artifacts instead of repeating the failed experiments.

## Current State

- Branch: `eg/c4-crush-resolution-controls`
- PR: https://github.com/pangenome/impg/pull/222
- Latest pushed commit: `98fd538` (`Fix SweepGA crush local seed frequency`)
- Local tests after that commit: `cargo fmt --check && cargo test` passed.
- Do not merge this PR solely on that commit. The commit fixes one observed
  local FastGA seed-frequency starvation case, but it was validated on a
  reduced CHM13/frozen-blunt experiment, not on the known-good GRCh38 C4
  command below.

The worktree also has unrelated local state that should not be committed as
part of the C4 fix:

- `vendor/gfaffix` is dirty as a submodule/worktree.
- `q.gfa`, `q.s.gfa`, `q.s.gfa.png`, and `renderings/` are untracked scratch.

## Known-Good C4 Artifact

The best current C4 result is:

```text
data/c4_crush_eval_20260523T140141Z/C4A.parent5k.sweepga_allvsall_fastga.k311.done.Ygs.mean-depth.paths.png
```

Generated graph/report:

```text
data/c4_crush_eval_20260523T140141Z/C4A.parent5k.sweepga_allvsall_fastga.k311.done.Ygs.gfa
data/c4_crush_eval_20260523T140141Z/C4A.parent5k.sweepga_allvsall_fastga.k311.done.Ygs.report.md
data/c4_crush_eval_20260523T140141Z/C4A.parent5k.sweepga_allvsall_fastga.k311.done.nosort.gfa
data/c4_crush_eval_20260523T140141Z/C4A.parent5k.sweepga_allvsall_fastga.k311.done.nosort.stderr
```

The command recorded by `/usr/bin/time` was:

```bash
target/release/impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=sweepga,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O data/c4_crush_eval_20260523T140141Z/C4A.parent5k.sweepga_allvsall_fastga.k311.done.nosort \
  -v 1
```

Observed resource use:

```text
wall: 4:39.95
max RSS: 55,106,276 KB
```

The reported graph is not perfect, but it is the right baseline:

```text
18089 segments
21018 links
465 paths
4590745 path steps
422007 segment bp
path-depth median 144, p95 620, max 1345
link jumps p95 11, p99 45, max 3270
path white-space bridges >=1000 bp: 25414
```

Interpretation from the report: focused connected locus, common start/end
nodes, still underaligned in places, still has repeated local contexts and a few
long link/path-jump artifacts.

## Important Correction

Do not compare the known-good image above to these reduced artifacts as if they
were the same experiment:

```text
data/c4_crush_rnd_20260523T235304Z/C4A.syng_blunt.mask_min3.gfa
data/c4_crush_rnd_20260523T235304Z/default_after_fix/*
data/c4_crush_rnd_20260523T235304Z/default_after_fix_filter/*
```

Those were based on a frozen CHM13-range blunt graph:

```bash
-r CHM13#0#chr6:31744284-31976975
```

The known-good artifact was generated from the GRCh38 range:

```bash
-r GRCh38#0#chr6:31891045-32123783
```

They are both C4-related, but they are not identical test inputs. I used the
CHM13 frozen graph for fast debugging and then compared it too loosely against
the GRCh38 known-good render. That was a bad comparison.

## What Was Actually Diagnosed

In the frozen CHM13 standalone crush experiment, the bad first pass inflated the
graph from roughly `394 kb` of segment sequence to roughly `15 Mb`. Debug dumps
under:

```text
data/c4_crush_rnd_20260523T235304Z/debug_fastga/
```

showed:

- `seqwish.gfa` and `unchopped.gfa` were byte-identical for replacements.
- The explosion was already present in seqwish output.
- Raw PAF covered only a small fraction of traversal names in several large
  candidates, so seqwish was being asked to place many paths without alignment
  support.

Raising the local FastGA seed frequency made that reduced experiment no longer
explode. That motivated commit `98fd538`, which makes `--method sweepga`
replacement frequency auto-scale instead of defaulting to the historical
whole-genome `10` seed cap.

This is useful, but it is not a complete explanation of the known-good GRCh38
C4 behavior because the known-good run used the full query path and already
looked nearly correct.

## Baseline To Reproduce First

The next person should start by reproducing the known-good command exactly on
current `HEAD`, then diff graph metrics against the existing artifact:

```bash
out=data/c4_handoff_repro_$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$out"

/usr/bin/time -v target/release/impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=sweepga,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/C4A.parent5k.sweepga_allvsall_fastga.k311.done.nosort" \
  -v 1 \
  > "$out/C4A.parent5k.sweepga_allvsall_fastga.k311.done.nosort.stdout" \
  2> "$out/C4A.parent5k.sweepga_allvsall_fastga.k311.done.nosort.stderr"

gfasort \
  -i "$out/C4A.parent5k.sweepga_allvsall_fastga.k311.done.nosort.gfa" \
  -o "$out/C4A.parent5k.sweepga_allvsall_fastga.k311.done.Ygs.gfa" \
  -p Ygs \
  -t 32

gfalook \
  -i "$out/C4A.parent5k.sweepga_allvsall_fastga.k311.done.Ygs.gfa" \
  -o "$out/C4A.parent5k.sweepga_allvsall_fastga.k311.done.Ygs.mean-depth.paths.png" \
  -m \
  -x 2200 \
  -y 1200
```

After that, run the graph report on the new sorted GFA and compare against:

```text
data/c4_crush_eval_20260523T140141Z/C4A.parent5k.sweepga_allvsall_fastga.k311.done.Ygs.report.md
```

## Questions To Resolve

1. Does current `HEAD` reproduce the known-good GRCh38 image and metrics?
2. If not, which commit changed behavior between the known-good artifact and
   current `HEAD`?
3. Is `98fd538` beneficial, neutral, or harmful on the exact GRCh38 command?
4. Are the remaining long links in the known-good graph from:
   - original syng/blunt graph glue,
   - SweepGA/seqwish replacement alignment,
   - polish/POA cleanup,
   - or sorting/layout only?
5. Does the second round need to be rejected by the visual-tail score, or is the
   score rejecting a biologically/visually better cleanup?

The report says the final known-good run attempted a second round, accepted
three replacements locally, then rejected the round because visual-tail score
grew by `0.0618` over the allowed `0.0200`. That round-level rejection may be
right, but it should be checked against the image rather than assumed.

## Recommended Next Step

Freeze the known-good GRCh38 command as an integration/regression test fixture
at the graph-stat level before changing the algorithm again. The first pass
should assert broad metrics, not byte-for-byte graph identity:

- paths remain `465`
- segment bp stays in the same order as `422007`, not multi-megabase explosion
- link/path jump p99 remains bounded near the existing report
- no large increase in sparse coverage / white-space bridges
- path sequences validate exactly

Only after that should the algorithm be changed to reduce the remaining
white-space and long-link artifacts.
