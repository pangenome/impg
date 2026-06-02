# Graph Quality Validation

This note captures the current validation pattern for local GFA graphs emitted by
`impg query -o gfa:syng:crush`.

## What "Good" Means

A rendered local graph is considered healthy when these checks hold together:

- Query succeeds within a bounded time on boring loci and known SV loci.
- The extracted paths are preserved exactly by `crush`; path validation is part
  of the resolver.
- POVU decomposition finds many eligible local sites in variable regions, not
  just a handful of candidates.
- The selected bubbles resolve without sending very dense high-copy tangles into
  allwave/seqwish.
- Sorted-order jumps are mostly local. Large jumps are acceptable in hard loci,
  but they should be visible and explainable in renders.
- 1D and 2D gfalook renders are non-empty. Stress should improve in 2D layouts;
  renderer NaN cases are tracked separately from GFA construction failures.

## Default Pipeline

The default syng GFA shorthand:

```bash
impg query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -b regions.bed \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 100k \
  -o gfa:syng:crush \
  -O output_graphs \
  -t 32 -v 1
```

`gfa:syng:crush` currently means syng native extraction, blunt GFA conversion,
one conservative POVU-guided crush round, and `gfasort` with pipeline `Ygs`.
During syng native extraction, the top 0.05% most frequent local syncmer nodes
are cloned by default, and rare repeated-copy local syncmer contexts are split.
This keeps repetitive syncmers as sequence evidence without letting them act as
global repeat-glue nodes.
`method=auto` routes very small budget-safe bubbles to global SPOA and sends
larger, wider, or high-copy bounded bubbles through pairwise graph induction.
The replacement seqwish step uses `seqwish-k=311` by default in crush so human
repeats require long exact matches before graph induction; residual bounded
bubbles are then handled by global SPOA polish. Direct SPOA/POASTA replacements
are accepted after exact path-sequence validation, while pairwise-induced
replacements remain guarded by local and round graph-quality checks.

## Current HPRCv2 Validation Panel

The current validation run lives in:

```text
data/graph_quality_sweep_20260521T162851Z/
```

Small control loci:

```text
small_cmrg.bed
small_graphs/
small_metrics.tsv
small_log_metrics.tsv
small_renders_1d/
small_renders_2d/
```

SV loci:

```text
sv_panel.bed
sv_graphs_full/
sv_metrics.tsv
sv_log_metrics.tsv
sv_renders_1d/
sv_renders_2d/
```

The SV panel contains AMY1A, KIR2DL1, CYP2D7, HLA-A, and C4A. The five-locus
batch completed in 8m43s wall time with a 45.9 GB peak RSS. The syng index load
is a large fixed cost; per-locus graph construction after load is much faster.

## Metrics

The metric TSVs record:

- segment/link/path counts
- total segment bases
- total path steps
- link and path jump percentiles after sorting
- duplicate segment sequence counts
- query intervals, raw syncmer nodes/edges, crush candidates/resolutions, sort
  time, and per-locus query elapsed time

These metrics are not a formal correctness proof. They are smoke tests that
catch obvious failures: missing paths, runaway graph size, no POVU candidates,
quadratic resolver behavior, bad sort ordering, blank renders, and unstable 2D
layouts.

## Clean Graph QC

Use `scripts/graph-clean-qc.py` for a deterministic "does this look like a
focused local graph?" check:

```bash
scripts/graph-clean-qc.py \
  data/graph_quality_sweep_20260521T162851Z/frequency_mask_validation/C4A.masked.gfa \
  --tsv data/graph_quality_sweep_20260521T162851Z/frequency_mask_validation/clean-qc.tsv
```

The checker reports and threshold-tests:

- connected components and largest-component fraction
- degree-0/1 tips, separating path-endpoint tips from unsupported internal tips
- most common path start and end nodes, with path support fractions
- sorted-order link jumps and consecutive path-step jumps
- duplicate segment sequence counts

The default profile is intentionally tuned for a focused SVR bubble with common
head/tail structure. Hard loci can fail this profile in informative ways. For
example, the current masked C4A graph passes, while the current masked AMY1A
graph fails on split start/end support plus one very long sorted-order link.

## Current Findings

- Small controls render cleanly. Their 2D stress improves strongly; MC1R is the
  noisiest small control and remains worth visual inspection.
- AMY1A, KIR2DL1, CYP2D7, HLA-A, and C4A all build successfully from the HPRCv2
  syng index.
- AMY1A exposed an allwave/seqwish stall when high-copy bubbles with hundreds
  of traversals were routed to allwave. The default auto router now sends
  bubbles above `auto-allwave-max-traversals=128` or
  `auto-allwave-max-total-seq=200k` to direct SPOA.
- AMY1A 2D gfalook stress currently becomes `NaN` and the 2D PNG is effectively
  blank even with bounded learning rate and Hilbert initialization. The 1D
  render and GFA parse are OK; this is tracked as a gfalook layout robustness
  issue for huge repeated path walks.
