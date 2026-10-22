# THE STREAM DISCIPLINE — the owner's standing ruling (2026-10-21/22; law for every future worker, including the merge fix when it resumes)

**THE RULING (verbatim intent, from the tranche-2 mine-and-reap
order):** extract the answer verdicts from the completed legs, then the
1.1T goes; "there shouldn't be any writing." A workspace that only
exists to produce ANSWERS must never accumulate raw intermediates —
the 318G `callA` dir of tranche 2 held ~300M PER TERRORY of
instrument-receipt streams beside a 64K calls.jsonl, and the verdicts
lived in a 68K truth-qv sidecar. That shape is banned going forward.

## THE RULES (standing law)

1. **RECEIPTS COMPUTED AT WRITE TIME.** Every instrumented path
   computes its receipts (counts, distributions, verdict tables,
   denominators) as it writes, in the same pass — never "dump now,
   analyze later." If a number is worth keeping, it is kept as a
   number with its n and its status label; the raw stream it was
   computed from is not kept.
2. **GATE COMPARISONS ON DIGESTS, IN-PROCESS.** Determinism and
   cross-form gates compare digests/hashes in-process at the moment
   both sides exist — never by materializing both sides to disk and
   diffing later. The tranche-2 determinism gate (57 IDENTICAL /
   9 GATE-DIFFERS at rvcut) is the shape: a one-line verdict per
   member, no banked copies.
3. **INSTRUMENT PATHS EMIT COUNTERS + BOUNDED WITNESS SAMPLES, NEVER
   RAW DUMPS.** An instrument may emit aggregate counters and a
   bounded, named witness sample (the first k + the worst-j items,
   k small and stated). Unbounded per-record/per-ingredient streams
   are the banned class — they were the bulk of the 1.1T.
4. **TRANSIENT FASTAS TO BOUNDED SCRATCH, REAPED AT EXIT.** Extraction
   FASTAs, window FASTAs, and other regenerable transients go to a
   named scratch with a hard size budget and are deleted by the
   runner's own exit path (success OR failure OR kill — trap-based),
   never left for a later session to reap.
5. **HARD WORKDIR BUDGET ENFORCED INSIDE EVERY RUNNER.** Every runner
   (stage script, fleet driver, call lane) enforces its own workdir
   budget — check before write, kill or fail-closed at the bound —
   the same discipline as the RSS guards. A stage that would grow the
   workdir past its stated budget stops and says so; it never silently
   fills the disk.
6. **MATERIALIZING RAW INTERMEDIATES REQUIRES AN EXPLICIT EXCEPTION.**
   If a genuine consumer needs a raw intermediate on disk (a
   restriction basis a later stage reads, a cache a call re-loads),
   the exception is written down AT CREATION naming (a) its consumer,
   (b) its reap time or reap condition, and (c) its budget. No
   exception, no materialization. Exceptions do not survive their reap
   condition.
7. **MINE BEFORE REAP, NEVER THE REVERSE.** When a workspace is
   closed, the answer receipts are computed and COMMITTED first; mass
   deletion happens only after. Partial legs are receipted as
   "partial: n loci complete of m" and never contribute answer loci.

## WHY (the receipts, not opinion)

- Tranche 2's answers fit in **376KB** (per-locus verdict tables,
  walls, counters, the 69 named prune truth-stakes, the marathon
  memory-tail resolution) — against **1.13T** of workspace mass, of
   which the product output (calls.jsonl + truth-qv) was ~10MB total.
- The panel asset question is the doctrine's own lesson: the
  once-forever routes artifact survived only because it lived outside
  the workspace; everything regenerable in the workspace was
  regenerable, and everything definitional (11MB of ledgers, axis
  rows, refined beds) was kept and is enough to re-run any leg.

## SCOPE

This is the standing rule for every future worker in this session —
the merge fix when it resumes, any tranche 3, any re-run of the
lost-to-stop legs (tranche 2's TB form-B legs; the 25/20 guard-killed
call legs, re-runnable under a quieter box). New runners copy the
tranche-2 mining shape (`scripts/genome-tranche/tranche2-mine.py`):
compute at write time, label every number with its n and status, reap
what you no longer need at your own exit.
