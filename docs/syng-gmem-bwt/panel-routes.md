# Automatic all-source panel routing (experimental v1)

This route advances beyond hand-authored finite alternatives. It automatically
retains every panel source lane, compiles native source-start profiles, declares
all panel-native identity topology families, and searches genuine successive donor
switches. Neither successful generation nor exact evaluation proves genome coverage
repair. No uncertified sequence is emitted on this route.

**Fixed depth/background, composite likelihood.** Native exposure and complete
assignment evaluation are exact under the error-free read/MEM observation operator.
The Poisson objective is a working **composite** likelihood: overlapping MEM
features are correlated, so this is not an exact probabilistic genome likelihood,
posterior or calibrated confidence. Depth and background remain fixed positive
inputs. No annealing, learned background, family-frequency/switch priors, frequency
masking, reference-member restriction or threshold tuning is introduced.

## Callable routes

```
impg genome-infer build-panel-routes --panel PREFIX --catalog catalog.json \
  --sources panel.agc --read-lengths 150 --out-dir panel-routes

# Automatic complete native assignment for an exact identity; no hand layout.
impg genome-infer evaluate-panel-routes --panel PREFIX --routes panel-routes \
  --sample sample.membwt --haploid-depth 10 --background 0.1 \
  --native-family 'A#0' --out-dir native-A

impg genome-infer search-panel-routes --panel PREFIX --routes panel-routes \
  --sample sample.membwt --haploid-depth 10 --background 0.1 \
  --max-work 100000 --max-evaluations 10000 --max-frontier 1024 \
  --max-optima 1000 --tie-epsilon 1e-9 --out-dir searched

impg genome-infer evaluate-panel-routes --panel PREFIX --routes panel-routes \
  --sample sample.membwt --haploid-depth 10 --background 0.1 \
  --assignment searched/incumbent-assignment.json --out-dir evaluated
```

All output directories must be new. Running/failed manifests and per-lane
checkpoints survive failure. Only a succeeded, checksummed graph can be loaded.
A succeeded search command may report budget exhaustion: inspect its result,
not just its process status or manifest. Legacy/v2/f7 calling, threading and
reconstruction formats are separate and reject these dedicated route artifacts.

## Generation rules and physical capacity

`panel-route-relative-poisson-v1` freezes the rule string, actual source-bound
operator identity, panel dictionary/sidecars, source file content fingerprints,
and actual source-access sidecars (`.fai`, plus `.gzi` for BGZF FASTA).
The builder binds these after build-time index creation; evaluators verify them
before opening any native reader. Missing/modified indexes are rejected, not
regenerated and adopted. Auxiliary bindings are separate from sequence input paths;
AGC inputs have no external access sidecars. Pre-fix graphs lack these bindings and
must be rebuilt, not silently upgraded. Read lengths and all topology families are
fixed. No existing sample or count-policy change;
the f7 explicit engine and its source identity remain unchanged.

Every full panel source ID/name/length and native continuity path survives,
including zero-anchor paths, partition-label mismatch, unknown ownership orientation
and internal inversions. PanSN sample/haplotype identity supplies a topology family;
its **entire native path inventory** supplies one start/end-paired linear molecule
per path. Family sizes are not forced to 17 or any biological chromosome count.
There is no preference for a particular family. Donors inside routes may be any
lane, regardless of family or ownership label.

The declared domain is **all panel-native identity path inventories and their
native start/end pairings**. These are assembly-path hypotheses, not assertions
of biological chromosome ends. Novel/unprovided endpoint or karyotype families
are excluded. Unknown families, endpoint changes and unknown JSON fields fail;
none are coerced into native pairings. A later versioned topology provider may
supply additional inventories/pairings; that extension seam is documented, not a
placeholder input pretending to support undeclared biology in v1.

Ownership is exact-once source annotation, not an allelic route chain. Duplicate
same-owner descriptions remain in `ownership-input.json`; normalized disjoint
coverage is streamed to `ownership.jsonl`. Conflicting owners/uncovered bases fail.
Large legacy feature-location and occurrence-multiplicity values are ignored while
ALL feature definitions and original IDs are streamed to `registry.jsonl`. The
legacy envelope header is explicitly recorded as unverified; a separate complete
raw-file fingerprint binds the import. No 41GB location matrix is loaded.

The builder extracts both independent raw native anchor views in bounded source
cores with word halos. It deduplicates physical `(source, anchor-position)` views,
validates the complete unambiguous DNA word, and emits two oriented ports. For
syncmer word width k and fixed cut d=floor(k/2):

- forward word: source cut `anchor+d`;
- reverse-complement word: source cut `anchor+k-d`.

Equal **oriented full DNA words** define switch hubs. They are declared sequence-cut
hypotheses, not inferred homology. Repeated words keep every physical occurrence.
Opposite orientations are views of the same capacity, never extra copies. Native
continuation does not depend on ownership labels or anchor polarity voting.

`ports.bin` is externally sorted fixed-width hub membership (`k DNA bytes`, little-
endian u64 source, u64 anchor, u8 reverse). Binary-search word bounds are its implicit
bucket index. Per-source files are sorted by `(anchor, reverse)` for traversal and
cut lookup. There is **no all-pairs edge/orientation matrix**. Sort buffers cap both
record count and word bytes; merge fan-in is 32.

Every traversed segment has positive source length, including between successive
switches. Thus zero-length switch cycles are impossible. Capacity is canonical
source-forward span overlap across the **whole selected genome**, shared between
orientations and every storage/ownership split. Distinct source paths remain
physical copies even when DNA-identical. Disjoint spans of one source may be used
in different molecules; overlapping spans may not. Harmless same-source contiguous
cuts normalize without adding rates, copies or a switch cost.

## Native baseline and exact lazy corrections

For every source/read length, compilation writes:

- `native/S-L.events.jsonl`: exact native start runs with sparse per-start counts;
- `native/S-L.index`: fixed-width `(start, byte-offset)` records;
- `native/S-L.totals.jsonl`: integrated exact totals for the full native domain;
- a lane checkpoint and content seals for all three files.

Cores own at most `--core-bp` starts (default 65,536, max 1,048,576), with complete
read-only halos. Anchor enter/leave events partition all admitted starts; the
unchanged native MEM collector preserves overlapping/count-two events. Counts are
RC invariant because reversing swaps the two independent raw collection passes;
this is orientation averaging, not adding two physical copies. Read lengths are
1..1,048,576 bp, at most 32 per graph. L > source has zero native starts.

For a realized route, the evaluator:

1. Sums native profiles only for starts owned by selected oriented source spans.
   Full domains use persisted totals; restricted domains seek through the start
   index and integrate intersecting event runs exactly.
2. Unions all seam-crossing windows `[b-L+1,b)` and removed molecule-terminal
   windows. Nearby switches are corrected **once**, not pairwise.
3. Subtracts the native contribution on that union, with checked integers.
4. Spells/adds actual admitted route contexts across arbitrarily many pieces.

For reverse traversal of `[a,b)`, local start t corresponds to native forward read
start `b-L-t`, not `b-1-t`. Signed/clipped transforms handle native source extension,
newly admitted starts, removed native starts at a declared molecule end and reads
longer than a piece/source/molecule. Unknown context is never substituted with an
old whole incidence total. Fresh singleton/public-BWT and frozen-f7 replay are
independent correctness oracles.

A bounded cache retains computed per-route integer contributions. It charges both
feature cells and segment bindings, including zero-feature routes. Cache hits save
computation; they do not deduplicate selected physical multiplicity. Each selected
route is added to the genome rate. Shards are content-verified on first use, before
being marked verified; failures cannot be bypassed by retrying the same evaluator.

## Fixed implicit universe and reported objective

The universe is independent of visited routes: complete sample-positive support,
all imported registry definitions, and all two-node features possible under the
frozen route/read-length rules. Instantiating all possible short route contexts is
neither required nor attempted. Unmaterialized zero-count/zero-selected-signal
features have symbolic relative contribution zero. **Every realized zero-count
feature with positive signal pays the full signal penalty.** Every sample-positive
factor is materialized exactly once, including unsupported positive residuals.

For read histogram H and integer orientation-averaged totals q:

```
e_f = sum_L H(L) * q_f(L) / sum_L L*H(L)
s_f = fixed_depth * sum_selected_physical_exposures
R   = sum_f [s_f - observed_f * ln1p(s_f / fixed_background)]
```

This is **background-relative NLL only**. No absolute NLL is reported from a
currently materialized subset. New zero features cannot change an earlier route's
baseline. Against f7, compare the complete explicit-union NLL minus **that same
union's background constant**. The removed constant depends on background: this
relative objective **cannot be reused unchanged to optimize beta**. Depth/beta
learning or a prior would require a separate declared model decision.

`result.json` binds graph/sample checksums, compiler/count policy, fixed parameters
and the objective kind. Evaluations include the complete realized/sample-positive
factor ledger, per-length integer counts, residual flags, exact admitted starts,
correction/cache diagnostics and explicit context/topology limits. Registry-only
symbolic zeros retain definitions/IDs in the sealed registry rather than millions
of redundant zero rows per evaluation.

## Actual search and honest guarantees

Search scores every native family first (unless its explicit budget is too small,
which is reported). It then performs deterministic, indexed successive-donor-switch
DFS: choose positive traversal to a DNA port, enter any matching donor port, and
continue until the original paired endpoint is reached. Complete routes are combined
across all slots with global source-span feasibility. This explores multi-donor
chains directly; a preexisting single-donor return bridge is not required. It does
not enumerate chromosome alternatives eagerly or stop at whole-template ranking.

Work counts native initializations and individual source-port/hub/DFS steps,
including skipped/infeasible repeat proposals. `--max-evaluations` bounds scored
complete assignments, `--max-frontier` bounds active traversal frames, and
`--max-optima` bounds correlated epsilon-optimum storage. None changes generation
rules or masks repeats. Budget exhaustion saves exact traversal cursors plus a
not-yet-pushed child in `frontier.json`; resume is not implemented. A streamed
`evaluations.jsonl` preserves every completed assignment score, including repeated
native initialization/DFS visits. Correlated alternatives are whole family/routes
assignments, not independently combinable marginals.

Only complete exhaustion of every family's globally feasible traversal/product
search certifies the finite rule-domain optimum. Alternative-storage exhaustion
withholds complete support, conservatively even if a later improvement could make
discarded ties irrelevant. A donor neighborhood exhaustion or good incumbent is
not a global certificate. Separate fields report native-initialization completeness,
exact native evaluation, frozen generation-rule completeness, search exhaustion,
objective/support certification and false biological topology completeness.

The always-valid nonnegative-signal lower bound uses an intentionally conservative
numerical relaxation. For c>beta, `ln(c/beta) <= c/beta-1` gives the lower bound
`-(c-beta)^2/beta` beneath the exact independent-factor minimum. Integer conversion,
products/division and sums are outward rounded; null means negative-infinity
arithmetic overflow, not an absent unknown claim. Unseen zeros cannot lower it.
It is weaker than evaluating the logarithmic relaxation but avoids claiming a
rigorous bound from unspecified libm logarithm rounding. Search does not certify
from approximate bound equality. Certification is exhaustive under its stated
f64 objective/tie tolerance, not probabilistic confidence.

## Synthetic validation gates

Default native/CLI tests require neither a private f7 executable nor `sha256sum`,
Python, or the Linux parent driver. They compare automatic profiles/search against
the independently implemented finite native engine in the current build, plus
public singleton-BWT unit checks. **This is not a frozen-f7 proof.** Separately
configured frozen-f7 gates retain independent bounded enumeration of 513 complete
routes through 52 oriented ports, matching every integer profile and relative
objective against the frozen binary's complete explicit feature union. A second
coupled gate uses 64 alternatives per molecule and atomic canonical source spans:
f7's 4,096-assignment Cartesian search and route DFS agree on exactly 64 correlated
feasible optima, not the invalid Cartesian product of their marginals. It exercises
simultaneous physical copy exchanges across two molecules. The parent driver also
passes a synthetic end-to-end execution. Independent review identified four
provenance/validation issues; all were reproduced, corrected and re-reviewed with
no remaining findings. The parent release gate passed **615 ordinary Rust tests**
across 13 suites, then explicitly ran and passed the three normally ignored
frozen-oracle/driver tests, plus **six Python tests**. Source hashes matched before
and after the gate. Real-panel performance and coverage remain separate acceptance
gates. Sidecar mutation/missing-index, same-source repeat versus actual
source/identity mixing, and malformed native stream regressions cover the review
corrections. Harmless contiguous segmentation is normalized before switch metrics.

Portable gates (use the normal offline native build environment):

```sh
cargo test --offline --locked --release -j4 --lib genome_inference::panel_routes -- --test-threads=1
cargo test --offline --locked --release -j4 --test test_panel_routes --test test_panel_route_coupling -- --test-threads=1
python3 tests/test_panel_route_driver.py
```

The private gates are explicitly ignored by default. On the configured Linux
parent host, with CPUs 252–255 and nice +10, opt in explicitly:

```sh
export IMPG_TEST_F7_FINITE=/path/to/sha256-pinned/frozen-impg
# Requires sha256sum; a missing/wrong executable is a failure, never a fallback.
taskset -c 252-255 nice -n 10 cargo test --offline --locked --release -j4 --test test_panel_routes --test test_panel_route_coupling frozen_ -- --ignored --test-threads=1
IMPG_TEST_PANEL_ROUTE_DRIVER=1 taskset -c 252-255 nice -n 10 cargo test --offline --locked --release -j4 --test test_panel_routes linux_parent_panel_route_driver -- --ignored --exact --test-threads=1
```

Set `IMPG_TEST_PANEL_ROUTE_OUTPUT` to a fresh absolute directory to retain the
synthetic driver fixture, all successful outputs, and the expected missing-FAI-SHA
rejection logs/status. Existing output directories are never reused.

The driver gate additionally requires Python 3, `/usr/bin/time`, `taskset`, `nice`,
and the explicitly configured parent native runtime environment. No private path
is substituted when configuration is absent.

## Native diagnostic startup and first full-panel attempt

The first real all-panel build was stopped after 578 of 9,901 source-path
checkpoints: unintended native per-match tracing had written 256,361,382,017 bytes
of diagnostics. Native `pathCount` and `PATH_DEBUG` both default to zero; raw replay
could reach their equality guard before an existing path-walking API suppressed
it. This was a logging/resource failure, not an inference or coverage result.
Search had not started. Checkpoints and bounded diagnostic excerpts were preserved;
the full runaway log was removed with user authorization.

`SyngIndex` now invokes the existing native debug-suppression helper once at
construction/loading, before publishing an index to callers. No vendor/gitlink,
counting, routing or likelihood changes are involved. Fresh-process route CLI
regressions reject native match/path dumps and enforce a 1 MiB diagnostic budget
for each small fixture command (not a cap on biological output files). Updating
`syng.rs` changes source-bound compiler identities: use fresh artifacts, not a
spoofed identity or an implied resume of the interrupted graph. Preserved frozen
executables remain independent count oracles.

## Parent-run driver and residual scalability limits

`python3 scripts/panel-route-gate.py` is a parent-only fresh output driver. `--help` describes
inputs; without `--execute` it prints a plan and opens no biological inputs. Actual
execution requires an exact candidate binary SHA256, the frozen f7 binary SHA256,
a parent inventory and a JSON mapping of all input/panel/frozen-artifact absolute
paths to SHA256, including each actual source-access `.fai` and applicable `.gzi`.
Prepare and freeze these indexes before `--execute`; the driver must not create
and adopt undeclared source-access inputs. It checks those inputs before/after, serializes pinned CPU/nice
commands, records `/usr/bin/time -v` logs, validates all native inventories/end
pairings, exact-once ownership and physical port/reverse membership. Independently
of Rust replay, it streams every actual native index/event pair: exact record-boundary
offsets, starts, contiguous nonoverlapping `[0, max(0,n-L+1))` coverage, positive run
lengths, run counts, byte sizes, EOF/truncation/trailing-data and zero-start checks.
It also exactly reconciles integrated sparse totals, using at most one profile's
`--max-profile-terms` token map (cap exhaustion fails rather than truncates the gate).
`generation-gate.json` records verified profile/run/start counts and totals verification.
It then runs all-native plus mixed search and re-evaluates the incumbent in a fresh
process. No truth option or sequence-emission command exists. Operational success
requires complete native initialization and at least one complete route crossing
source path IDs, not merely a same-source repeat deletion/inversion or two different
chromosomes in an assignment. `mixed_assignments_evaluated` counts those cross-source
assignments; `non_native_assignments_evaluated` separately counts rearranged
assignments. The retained `donor_transitions_examined` name is a hub-successor work
counter (including same-source jumps), not evidence of mixing.
`mixed_identity_assignments_evaluated` requires an actual
`sample#haplotype` change **within one route**. Different source chromosomes of one
identity are not mixed donor identities; zero identity mixing is reported honestly,
not treated as donor-identity success. Complete correlated assignments remain the
scoring/support unit. Operational success still does not imply complete search or
repaired coverage.

The driver is tested on synthetic data only. The real parent inventory's 9,901
paths / 3,336,976,856 bp / 235 families is an input validation target, **not a worker
run result**. Real full-panel generation/search and biological controls remain
parent-owned after independent review.

Remaining costs are explicit: native event/index/totals storage can be large; shard
verification performs full sequential reads on first use; compact ownership metadata,
complete sample support, realized sparse counts and the bounded route cache occupy
RAM. `--max-profile-terms`, `--max-feature-terms` and `--cache-terms` control resources,
never feature scope. Failures retain partial checkpoints and do not claim complete
generation/evaluation. Search is exponential, deterministic-order/budget dependent,
rechecks span feasibility and accumulates full-genome factors for each complete
assignment; it is not yet an incremental optimizer or a scalable proof method.
Sealed source/shard inputs must remain immutable within an evaluator session;
the parent driver independently fingerprints all declared frozen inputs before
and after its run. FNV seals are accidental-corruption checks, not adversarial
authentication. More efficient cumulative indexes, genome-factor deltas, frontier scheduling/resume
and tighter conservative bounds are future engineering, not silently implemented
capabilities. All source candidates and physical copies remain in the declared
domain despite those exploration limits. No resolved sequence is authorized even
by this version's exhaustive certificate; emission/probabilistic uncertainty needs
its own reviewed contract.
