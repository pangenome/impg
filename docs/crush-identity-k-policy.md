# Identity-Driven Crush `seqwish-k`

Crush replacement graph induction currently has an explicit `seqwish-k`
minimum exact-match length. A high value such as `311` is useful in human
repetitive sequence because it is just above Alu length and prevents local graph
induction from gluing through short off-diagonal repeat matches. It should not
remain a magic constant forever.

The automatic policy should be driven by a global identity estimate for the
sequence set being crushed.

## Model

For two homologous sequences with per-base identity `p` and effective aligned
length `L`, a rough independent-site model gives the expected number of exact
runs of length at least `k` as:

```text
E[runs >= k] ~= L * (1 - p) * p^k
```

Solving for `k` at a target number of exact anchors `A` gives:

```text
k_supported ~= ln(A / (L * (1 - p))) / ln(p)
```

Because `ln(p)` is negative, lower identity or shorter sequence produces a
smaller supported exact-match length. This is the key reason a fixed `311` is
not universal: it is reasonable for very high identity human haplotypes, but it
can underalign lower-identity bubbles or short local regions.

## Policy

The automatic replacement `seqwish-k` should use:

```text
k_auto = min(repeat_guard_k, k_supported)
```

then clamp to practical bounds, for example:

```text
min_k <= k_auto <= max_k
```

For human high-identity whole-genome panels, `repeat_guard_k=311` is a sensible
default because it is just above Alu length. For lower identity estimates,
`k_supported` should pull the value down so the replacement graph still has
enough homologous anchors. In that regime, repeat control must come mostly from
the pairwise alignment filter: SweepGA scaffold/plane-sweep filtering or
AllWave-selected pairs, not from seqwish transitive closure alone.

## Implementation Target

Add an automatic mode alongside the existing explicit integer:

```text
seqwish-k=auto
identity=0.99
target-anchors=3
repeat-guard-k=311
min-seqwish-k=31
max-seqwish-k=511
min-map-length=0
min-identity=0.0
```

The explicit integer remains exact and reproducible:

```text
seqwish-k=311
```

When `seqwish-k=auto`, crush should compute the value per candidate from the
candidate traversal statistics and a global identity estimate. The first
implementation can take `identity` explicitly. Later, impg can estimate it from
a sample of accepted pairwise alignments and log the estimate once per run.

The pairwise mapping filter should track the same scale. In current CLI terms,
`min-map-length=0` means "use the resolved `seqwish-k`", while an explicit
positive value overrides it. `min-identity` is independent: it filters
inappropriate chains before seqwish sees them, which is important when large
parent bubbles contain repetitive local homology that is chainable but not
useful for local graph induction.

## Logging

Every replacement run should log the resolved policy:

```text
crush seqwish-k auto: identity=0.9900, L=10000, target-anchors=3,
repeat-guard-k=311 -> k=311
```

This makes the behavior reproducible and prevents hidden changes in graph
condensation.
