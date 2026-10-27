# THE AGC DECOMPRESSOR POOL — the K-sweep of record

> **WHAT THIS IS.** The measurement that prices the memory fix. `AgcIndex`
> used to open one AGC `Decompressor` per ambient thread slot; a
> decompressor holds the archive's own decompression state in RAM, and
> that state measures about the archive's own size. On a 192-thread box
> that was 192 private copies of a 3.3 GB archive — the routes
> derivation's measured **711 GB** peak. The pool is now bounded by
> `IMPG_AGC_DECOMPRESSOR_SLOTS`. This note is the five-point sweep that
> says what bounding it costs in wall, and where the knee is.

## The sweep (one run, one method, one binary)

Fix `RAYON_NUM_THREADS=192` and vary **K alone** — the pool size. Every
point must emit a byte-identical `graph.json`, because sharing a
decompressor may not change a byte.

| K | peak RSS | wall | core-s | eff. concurrency | `graph.json` |
|---|---:|---:|---:|---:|---|
| 8 | 88.58 GB | 4729 s (78.8 min) | 54,222 | 11.5 of 192 | `ef7f39ea…` |
| 16 | 134.96 GB | 2710 s (45.2 min) | 53,909 | 19.9 | **identical** |
| 32 | 219.64 GB | 1536 s (25.6 min) | 51,391 | 33.5 | **identical** |
| **64** | **367.09 GB** | **1129 s (18.8 min)** | **50,802** | 45.0 | **identical** |
| 192 | 861.71 GB | 1318 s (22.0 min) | 62,169 | 47.2 | **identical** |

- **The identity gate holds at every point**: `graph.json` is
  `ef7f39ea d80fe8bd af6c77e9 91b3fd31 c88dd0ba 327d835c bd8c9c05
  1d9022ff`, 34,682,464 bytes, at K = 8, 16, 32, 64 and 192. The bounded
  pool is answer-preserving across a 24× range of slot counts.
- **Effective concurrency** = core-s ÷ wall: how many cores were actually
  busy. It tracks K to 32 (11.5 → 19.9 → 33.5) and then **saturates at
  ~45 of 192**. The loop is not thread-bound past K≈32; something else
  (per-lane decompression serialization inside AGC) is the limit.

## The memory model

A five-point power fit over the same run:

> **peak ≈ 18.8 GB × K^0.720**

| K | predicted | measured |
|---:|---:|---:|
| 32 | 228 GB | 219.64 GB |
| 64 | 376 GB | 367.09 GB |
| 192 | 829 GB | 861.71 GB |

Within 4% at every point. The marginal cost per slot *falls* with K
(5.80 → 5.29 → 4.61 → 3.86 GB/slot), i.e. the relation is genuinely
sublinear — an earlier two-point linear fit predicted 1,090 GB at K=192
and was wrong.

## The verdict: the committed default is dominated

**K=192 — one decompressor per thread, the shipped behaviour — is worse
than K=64 in every dimension simultaneously:** +135% memory, +22% CPU,
**+17% wall**. It buys nothing and costs everything. The knee is K=32–64.

- **K=64** is the wall optimum (18.8 min) and the CPU optimum (50,802).
- **K=32** is the citizenship optimum: 25.6 min for **40% less memory**
  (220 GB vs 367 GB) on a step that runs **once per panel, forever**.
- **K=8** (the shipped default for the *call* paths) is 78.8 min at
  88.6 GB.

## The measurement-era caveat (stated, not hidden)

The campaign's earlier receipt for this step recorded **163,159 core-s**
at 192 threads. This sweep measures **62,169 core-s** for the same
nominal configuration. The two are not reconcilable, and that era has a
**known core-s capture bug** (the `TIMEFORMAT`-with-redirect fix is in
this repo's own notes). The consequences:

- The *robust* claims here are the same-run ones: flat CPU to K=64, the
  ~45-core concurrency ceiling, and K=192 being dominated by K=64.
- The earlier inference that "~109,000 core-s were burned on contention
  at 192 threads" is **RETRACTED** — it rested on that single
  unreproducible number.

## The operating point, and which default belongs where

The same slot count serves two very different jobs, and they want
opposite things:

- **the whole-panel routes derivation** — a once-per-panel, whole-box job
  over 38,790 lanes: wants many slots (32–64).
- **a per-territory call** — 16 threads, and per-lane memory is *the*
  whole-genome residency constraint: wants few slots (8).

So the shipped default stays at **8** (it bounds the call path's
per-lane appetite), and the derivation opts in. Today that opt-in is the
env knob:

```bash
IMPG_AGC_DECOMPRESSOR_SLOTS=32 impg infer derive-routes …
```

The named follow-up is to thread the hint through `Sources::open` so the
derivation sets its own default (~min(threads, 32)) and the knob is not
required for a good number — the plumbing is the reason it is not in this
change.

## Reproduction

```bash
KS="8 16 32 64 192" bash scripts/complexity-ledger/run-agc-pool-ksweep.sh
```

Marker-gated per K, self-draining (each point's artifact is removed after
its receipt is mined), gated on the sweep's own disk headroom, and gated
on the `graph.json` identity above. Receipt:
`/mnt/nvme2n1/erikg/agc-pool-ksweep/ksweep-receipt.txt`. Total run 3.2 h.

No PR PUSH (the owner pushes).
