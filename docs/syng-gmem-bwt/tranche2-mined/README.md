# THE TRANCHE-2 MINED RECEIPTS (the stream doctrine's own shape: computed at mining time, before the reap)

**THE OWNER'S ORDER, executed:** mine the answer verdicts out of the
completed legs into small receipts, keep the definitional smalls, reap
the workspace mass. Nothing raw is copied — these files ARE the
tranche-2 answers. Workspace:
`/mnt/nvme2n1/erikg/hprcv2/global-panel-call-mode/genome-tranche2`
(reaped after these receipts were committed; see the reap receipt at
the bottom).

## THE INVENTORY VERDICT (the surprise: only ONE leg was frozen mid-write)

109 members drawn (S_EXP 22, S_BIND 71, S_MAR 16); 108 runnable
(T00762 excluded: extract-length-mismatch, `tranche-failed.tsv`).

| leg | terminal (answers) | terminal-fail (named class) | partial |
|---|---:|---:|---:|
| form A (restricted) | **82** | 25 (21 rss-guard, 2 verify-key, 1 qv-biwfa-helper, 1 overlapping-records) | **1** |
| form C (restricted+prune) | **88** | 20 (15 rss-guard, 2 verify-key, 2 qv-biwfa-helper, 1 overlapping-records) | 0 |
| form B (global-native) | 0 | 0 | 0 — **NEVER STARTED** (the master was stopped mid-TA; TB never opened) |

The ONE partial leg: `callA/T04063` — **partial: 0 loci complete of 14**
(killed by the stop's external SIGTERM mid-run at 89.4GB RSS, 13.5h
wall). Every other non-answer leg has a terminal exit receipt with a
named failure class (see `metrics.txt` and `inventory.tsv`); the
rss-guard kills are honest guard-kills (the memory tail), re-runnable
later under a quieter box; T04066 is the deterministic
overlapping-records fail-closed class (same family as tranche-1's
T06027/T05603); T00717/T01909 are the deterministic verify-key
fail-closed class; T00157/T02597 the qv-biwfa helper fail-closed class.

The failed legs' `calls.jsonl` are all EMPTY (0 bytes) — zero loci
landed in any of them; no partial loci are counted anywhere in these
receipts. Verdict loci come ONLY from `call{A,C}.done` legs.

## THE MINED DENOMINATORS (every number labeled with its n and status; full tables in metrics.txt)

- **The wash (rank-1 of truth-expressible, TERMINAL legs):**
  S_EXP A 89/391 (23%, 20 terr) vs C 87/391 (22%, 20 terr);
  S_BIND A 275/1251 (22%, 61 terr) vs C 267/1270 (21%, 62 terr);
  S_MAR A 5/20 (1 terr) vs C 36/129 (28%, 6 terr).
  QV median 43.5 everywhere; >=QV30 92-93%.
- **The prune's own verdict (1,643 paired C-vs-A loci, both forms
  TERMINAL): 28 WINS, 41 LOSSES — the zero-loss gate bar does NOT hold
  at tranche-2 scale**; both sides concentrated at the heavy-share
  members (all 69 named in `metrics.txt` and rowed in
  `prune-pairs-C-vs-A.tsv`). QV-visible changes at 210 paired loci
  (76 higher under C, 134 higher under A).
- **The prune's wall value: median C/A 0.74x over 81 paired legs
  (0.63x at far-bin share > 1%, n=66; 1.03x at share ~0, n=7)** — the
  heavy tail pairs drop 3-9x (e.g. T00824 84.0% share: 231.9s to 26.5s).
- **The marathon roster (16, outside every stratified estimate): the
  prune RESOLVES most of the memory tail** — 6 of 16 terminal under C
  where only 1 landed under A, with the guard-killed A peaks at
  67-133GB vs the surviving C peaks at 43-63GB (e.g. T00129: A peak
  129GB guard-kill vs C 63GB terminal). Full per-member receipt in
  `metrics.txt`.
- **The census counters (all 108 runnable members, both cuts
  TERMINAL):** census-mult B/A median 1.11-1.12x (S_EXP/S_BIND);
  far-bin share median S_EXP 30.1% (max 85.4%) / S_BIND 13.7% (max
  88.9%) / S_MAR 23.4%; **tranche-2 aggregate 36.17%** of the
  restricted census mass far-bin scattered; form-C machinery prune
  receipts cross-checked: 107 agree, 0 disagree.
- **The universe-wall fragment tranche 2 alone can state** (S_EXP +
  S_BIND only; tranche 1 owns S_ORD/S_SEG/S_REP and 8 of the 79
  S_BIND): form C 3.20e6 core-s (+/- 3.04e6, 95% CI) = 3.5h +/- 3.3h at
  256 cores — a FRAGMENT, never a full-universe total.
- **Form B: ZERO denominators** — TB never started; the
  conserved-elsewhere rate is NOT computable from this workspace
  (named lost-to-stop; tranche 1's 1-in-981 remains the standing
  receipt).

## WHAT WAS KEPT (and why) — moved to /home/erikg/hprcv2/panel/tranche2/ (11MB)

- `tranche.tsv` (the draw ledger), `truth2.tsv` (the truth2 bridge
  rows), `tranche-failed.tsv` — the tranche's own definitional
  ledgers, not regenerable without re-running the bridge.
- `bridge/target.specs`, `bridge/window.specs`,
  `bridge/truth2-bridge.paf` — the truth2 bridge's own definitional
  receipts (what was extracted and aligned; the FASTAs themselves are
  re-derivable from the AGC and were reaped).
- `members/<T>/{axis-path.txt, locus.refined.bed, reference-axis.json,
  extract.specs.txt, locality-domain-structure.txt}` for all 109
  members — the per-territory definitional beds/ledgers (the axis row,
  the refined-bed rows the truth pair is drawn from, the locus
  census). The heavy per-territory mass (locus.fa, syng indexes,
  export graphs/walks, partition beds, locality graphs) is REGENERABLE
  via the committed prep chain and was reaped.
- All stage scripts and the bridge/extract scripts were already
  committed in the repo (`scripts/genome-tranche/`,
  `scripts/runtime-verification/`); the workspace copies were verbatim
  (verified: `run-truthcalls.sh` differs only in its WT default).
  The mining script itself: `scripts/genome-tranche/tranche2-mine.py`.

## THE PANEL ASSET — NOT HERE (receipted as ordered)

The three-step doctrine's once-forever panel asset is NOT in this
workspace: the syng panel itself lives at
`/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng` (outside the
workspace, untouched), and the routes artifact (`global-routes-par`)
lived in the ROOT-DISK copy (`/home/erikg/hprcv2/global-panel-call-mode/`),
which is already reaped — that path no longer exists. It is
re-derivable (routes ~25min once/panel, per the owner). Nothing in
this workspace was the panel asset; nothing panel-definitional was
lost with the reap.

## FILES

- `metrics.txt` — the computed decision receipts (inventory verdict,
  wash, prune verdict, walls, marathon receipt, census counters, rvcut
  gate, fleet/prod receipts)
- `inventory.tsv` — per member: every leg's status/exit/wall/RSS
  peak/class
- `verdicts-A.tsv`, `verdicts-C.tsv` — the per-territory per-locus
  verdict tables (TERMINAL legs only)
- `prune-pairs-C-vs-A.tsv` — the prune's paired verdicts (69 truth
  stakes + the differing no-stake loci)
- `territory-walls.tsv` — the serial phase walls per territory per form
- `census-counters.tsv` — the cut-side counters per territory
- `receipts/` — verbatim small terminal receipts (rvcut determinism
  gate + campaign, stage.log, co-tenant.log, tranche-failed.tsv,
  bridge specs + paf)

## THE REAP RECEIPT

(to be appended after the deletion — receipts were committed first,
deletion second, never the reverse order.)
