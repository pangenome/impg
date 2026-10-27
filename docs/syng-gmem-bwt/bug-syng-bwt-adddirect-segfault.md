# BUG: `impg syng` segfaults on a single sequence (syngBWTadd → addDirect)

**STATUS: OPEN, reproduced under the pinned syng (054acfc), deterministic.**

Found while re-running the campaign's `exit 139` terminal territory (T00024,
`chr1:9,750,000-10,000,000`) against the *corrected* pin — i.e. after the
worktree's stale submodule (`68ac197`, which predates `b092a4c "Fix
uninitialised rs.dynamic->start in buildDynamic when nRun == 0"`, the fix
that landed with #247) was replaced by the pinned `054acfc`. **The crash is
therefore not that class** — it survives the fix.

## THE REPRO (321 KB, one sequence)

```bash
# save one FASTA record:
#   >HG03927#1#CM086638.1:130277163-130594662   (321,510 bp)

impg syng -f one-sequence.fa \
  --syncmer-length 63 --smer-length 8 --syncmer-seed 7 \
  --position-sample-rate 64 \
  -o out.syng
# -> Segmentation fault (core dumped), exit 139
```

## CHARACTERISATION (measured)

| probe | result |
|---|---|
| the same sequence, run twice | **139, 139** (deterministic) |
| `--position-sample-rate` 64 / absent / 8 | **139 in all three** — the flag is not implicated |
| the current binary (`target-frontier/release/impg`) | 139 |
| the committed substrate binary (`target/release/impg`) | 139 |
| the full 478 MB territory `locus.fa` (56 sequences) | 139 on the **5th** sequence |
| 5 sequences from that file (1.6 MB) | 139 |
| sequences 4+5 together (708 KB) | **hangs** (>200 s), does not fault |

So the trigger is content-shaped, not size-, flag- or binary-dependent, and
there are TWO symptoms: a segfault for one input and an unbounded hang for a
nearby one. Neither is a routing or threshold question; both are in the BWT
insertion.

## THE STACK (gdb, current binary)

```
Thread 1 "impg" received signal SIGSEGV, Segmentation fault.
0x000055555662fb43 in addDirect ()
0x0000555556629baf in syngBWTadd ()
0x000055555662a50d in syngBWTpathAdd ()
0x0000555556519033 in <impg::syng::SyngIndex>::add_sequence_internal ()
0x00005555565151d5 in <impg::syng::SyngIndex>::add_sequence ()
0x0000555555af16b8 in impg::run ()
```

`addDirect` is the direct-insert path of the dynamic BWT
(`syngbwt3.c`, called through `syngBWTadd`). The Rust side is only the
caller: the sequence is added one record at a time and the fault is inside
the C insertion.

## WHY IT MATTERS BEYOND THE CAMPAIGN

`impg syng` is the panel-side index builder for every route: a substrate
whose `syng` step faults cannot be mapped or called at all. In the campaign
this class appears as the fail-closed **SUBSTRATE-TERMINAL POISON** entries
(`T00024`, `T00197`, both `exit 139`) — and, per the hang probe above, the
`low-complexity hang` class shares the entry point.

EIGHT further territories were poisoned on the *extract* assertion
(`len N != bed N+1` at AGC contig boundaries), which is a different class
and a substrate-ledger finding.

## THE HONEST CAVEAT FOR THE CAMPAIGN RECEIPTS

The campaign's earlier receipts were produced by a worktree whose syng
predated `b092a4c`. Main's pin has that fix; the crash above happens *with*
it. So:

* the terminals are **not** artefacts of the stale submodule — they
  reproduce on the pin, and
* any campaign run made before the submodule was corrected describes a syng
  that main does not have. Both statements belong in the datasheet.

## NEXT

The repro is small enough to hand upstream as-is. A fix needs
`addDirect`'s insert path understood (the campaign's `IMPG_PATCHES.md`
records the impg-side divergences in this tree, which is the place to
check first for a local patch).
