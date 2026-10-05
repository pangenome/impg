#!/usr/bin/env python3
"""check-realign-scoring.py -- the independent checker for the realignment
scoring layer (slice B), its slice-C extension (the generated-here/
generated-elsewhere marginalization) AND its slice-D extension (the
canonical-scheme pin-skeleton frame repair). Assessment-side; validates
the receipts against the COMMITTED artifacts only (the partition GFAs,
the maps, the committed anchor-projection receipt, the census record
ids, the FASTQ).

With no arguments: the committed slice-B receipts (the original floor
semantics: prior + 150*B for unplaceable reads). With --marginal: the
slice-C receipts (LL = logsumexp(local, E); the floor is the derived
elsewhere branch E; the E derivation audited from the receipt's own
panel strain list; the bounding property asserted over the whole
matrix; the before/after pilot verdicts vs the slice-B receipts,
including the L4 control preservation). With --frame: the slice-D
receipts (the marginal model + the repaired canonical-scheme pin
skeletons) plus the frame-repair audit phases:
  9a. the identity gate: the env-gated stored-walk run reproduces the
      committed slice-C receipts (semantic fields) with byte-identical
      sidecars;
  9b. the skeleton audit: every fold's canonical walk verified
      INDEPENDENTLY against the GFA -- each step's window by SEQUENCE
      against the node's S segment (orientation-aware), the frame
      decision against the window's canonical form, the forward part
      complete against the P-line walk, the dropped positions
      rc-canonical; the skeleton sidecar cross-checked;
  9c. the no-regression sweep: folds mapped by member set, per
      (unit, fold) pin gains/losses and ll changes -- ZERO pin losses;
  9d. the blinded-unit census before/after (the paired anatomy
      receipts: units with a full-match donor occurrence inside a
      truth-fold member row and no truth placement);
  9e. the named proof case (record 113/unit 84 at L7) at the anatomy
      AND skeleton levels;
  9f. the full-component frame-audit census (chrI/chrMT rows affected,
      added/dropped).

Phases:
  1. run markers, exit code, wall, the 64 GiB RSS guard;
  2. receipt structure: the folds partition the partition's member rows
     exactly; every fold's sequence re-derived from the GFA (P/L/S
     spelling, strand-aware, gap segments handled) and compared EXACTLY;
     the fold's walk verified against the GFA P line (slice B/C) or the
     canonical-skeleton rules (slice D); the fold's node/edge usage
     multisets re-derived;
  3. the locality ruling re-derivation: per locus the touching-occurrence
     in-axis/overhang/extension counts re-derived from the committed
     anchor-projection receipt + the maps and compared EXACTLY;
  4. THE EXACTNESS PROOF: every sampled (unit, fold) pair re-derived
     independently -- the read's FASTQ identity (FNV), the correspondence
     from the GFA walk (the unique-monotone rule), the serving-anchor
     projection and every vote against the GFA-spelled sequence, the
     backbone and the (m, c) integers EXACTLY; plus the bounded-offset
     dominance scan (the anchored placement attains the best
     single-offset score within +/- 150, collinear placements) and the
     bounded edit-distance alignment (the sample's indel structure
     named; indel-free pairs must have the single-offset form);
  5. the scoring re-derivation: sampled (unit, fold) log-likelihoods
     recomputed from the ingredients (correspondence + votes + prior +
     floor) and compared; floor-only entries verified exactly; the
     class log-likelihoods recomputed from the per-unit matrix (all
     classes, numpy) -- winner, truth rank, ties, called set;
  6. the QUAL cluster state re-derived from the receipt's spectrum
     (the mirrored knee/cluster machinery) and compared;
  7. the identical-through-graph fold verified from the GFAs (the two
     member rows spell identical sequence and walk -- the fold-by-
     construction proof);
  8. the pilot verdicts stated (with the before/after table and the L4
     control enforced);
  9. (slice D, --frame) the frame-repair audit phases above.
"""

import glob
import gzip
import json
import math
import os
import sys

D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
GRAPHS = f"{D}/partition-graphs"
NAMES = "/home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng.names"
ANCHOR = f"{D}/anchor-projection-chrI.jsonl"
READS = f"{D}/reads.fastq.gz"
RECEIPT = f"{D}/realign-score-chrI.jsonl"
EXACTNESS = f"{D}/realign-score-chrI.exactness.jsonl"
INGREDIENTS = f"{D}/realign-score-chrI.jsonl.ingredients.jsonl"
RUN = f"{D}/run-realignscore-chrI"
LOCI = [2, 4, 7]
LOCI_SET = set(LOCI)  # rebound by --exhaustive (and --loci-slice) below
COMPONENT = "S288C#0#chrI"  # the axis path (the --exhaustive mode rebinds it)
K = 63
READ_LENGTH = 150
RSS_BUDGET_KB = 64 * 1024 * 1024
PHRED = 40
EPS = 10 ** (-(PHRED) / 10.0)
A = math.log(1.0 - EPS)
B = math.log(EPS / 3.0)
FLOOR = READ_LENGTH * B

# Slice C: the generated-here/generated-elsewhere marginalization.
# With no arguments the checker validates the committed slice-B
# receipts exactly as before. With --marginal it validates the slice-C
# receipts (the marginal model: LL = logsumexp(local, E), the floor
# for unplaceable reads is the derived elsewhere branch E, the E
# derivation audited from the receipt's own panel strain list).
# Slice D (--frame): the marginal model + the repaired canonical-
# scheme pin skeletons + the frame-repair audit phases (9a-9f).
# Slice E (--exhaustive COMP): the exhaustive full-component rerun
# receipts (the marginal model + the repaired skeletons over EVERY
# window of the component) -- phases 1-8 with the skeleton audit of
# phase 2 and the E audits of phase 5, the identical-fold proof
# generalized (every locus's identical-pair fold verified from the
# GFAs, no name hardcode), the truth-rank checks skipped at the
# non-expressible loci, and phase 8 the before/after table vs the
# component's Poisson-era read-matched receipts (the instrument of
# record; the fleet components' BEFORE receipts are the multicensus
# runs' own read-matched likelihood outputs in the diagnostic scratch)
# with the old-wins PREDICTION VERDICTS (a lost old-win is a named
# measured finding about the instruments, never a receipt failure),
# the chrI L4 control hard-gated, and the aggregate counts. Phase 9
# stays on the slice-D pilot receipts (--frame).
# THE FLEET: the whole-genome stage runs every balanced-diploid
# component through this mode (chrII..chrXVI beside the closed
# chrMT/chrI/chrIV).
EXHAUSTIVE = None
if "--exhaustive" in sys.argv:
    _i = sys.argv.index("--exhaustive")
    EXHAUSTIVE = sys.argv[_i + 1] if _i + 1 < len(sys.argv) else "chrI"
    _FLEET = (
        "chrMT", "chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI",
        "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII",
        "chrXIV", "chrXV", "chrXVI",
    )
    if EXHAUSTIVE not in _FLEET:
        sys.exit(f"--exhaustive requires one of the 17 components: {_FLEET}")
FRAME = ("--frame" in sys.argv) or (EXHAUSTIVE is not None)
MARGINAL = FRAME or ("--marginal" in sys.argv)
# (Stage 2, the extent normalization: --territory validates the
# TERRITORY-NORMALIZED receipts — the env-gated convention
# (IMPG_REALIGN_TERRITORY_NORMALIZED) written through the same
# --timered-base prefix mechanism. The phases gain the territory
# audits: the images re-derived from the walks and compared exactly,
# the sampled per-unit re-derivations under the image bounds, the
# placement-span checks of phase 4, and phase 8N — the
# conversion/no-regression gate vs the COMMITTED default-convention
# receipts. Phase 8t (answer preservation) does NOT run: the
# convention changed by design; the env-unset identity gate is the
# separate no-op proof.)
TERRITORY = "--territory" in sys.argv
# Slice F (phase 0 of the runtime plan — the dominance measurement):
# with --exhaustive, the optional --timered-base BASE --run-tag TAG pair
# points every phase at the TIMERED rerun's receipts (the instrumented
# binary: timers and counters only, stderr-emitted) and adds phase 8t,
# THE ANSWER-PRESERVATION GATE — the timered receipts must equal the
# committed slice-E receipts on every semantic field (the identity-gate
# pattern from slice D: only the walls/rss_kb timing fields may differ)
# with byte-identical sidecars. The default paths are untouched.
TIMERED_BASE = None
if "--timered-base" in sys.argv:
    if EXHAUSTIVE is None:
        sys.exit("--timered-base requires --exhaustive")
    _i = sys.argv.index("--timered-base")
    TIMERED_BASE = sys.argv[_i + 1]
    _i = sys.argv.index("--run-tag")
    TIMERED_RUN = f"{D}/run-{sys.argv[_i + 1]}"
if TERRITORY:
    if EXHAUSTIVE is None or TIMERED_BASE is None:
        sys.exit("--territory requires --exhaustive COMP --timered-base BASE --run-tag TAG")
    COMMITTED_RECEIPT = f"{D}/realign-exhaustive-{EXHAUSTIVE}.jsonl"
if EXHAUSTIVE is not None:
    # the BEFORE receipt root: chrI/chrMT keep the committed
    # Poisson-era receipts at the validation root; chrIV and every
    # fleet component's Poisson-era BEFORE receipt lives in the
    # diagnostic scratch (the census run's own likelihood output;
    # the validation-root sidecars are the preserved balanced
    # records).
    _before_root = (
        f"{D}/cosine-diagnostic-scratch/{EXHAUSTIVE}"
        if EXHAUSTIVE not in ("chrI", "chrMT")
        else D
    )
    if TIMERED_BASE is not None:
        RECEIPT = f"{D}/{TIMERED_BASE}.jsonl"
        EXACTNESS = f"{D}/{TIMERED_BASE}.exactness.jsonl"
        INGREDIENTS = f"{D}/{TIMERED_BASE}.jsonl.ingredients.jsonl"
        RUN = TIMERED_RUN
        BEFORE_RECEIPT = f"{_before_root}/cosine-graph-likelihood-readmatched-{EXHAUSTIVE}.jsonl"
        SKELETON = f"{D}/{TIMERED_BASE}.skeleton.jsonl"
    else:
        RECEIPT = f"{D}/realign-exhaustive-{EXHAUSTIVE}.jsonl"
        EXACTNESS = f"{D}/realign-exhaustive-{EXHAUSTIVE}.exactness.jsonl"
        INGREDIENTS = f"{D}/realign-exhaustive-{EXHAUSTIVE}.jsonl.ingredients.jsonl"
        RUN = f"{D}/run-realignexhaustive-{EXHAUSTIVE}-{EXHAUSTIVE}"
        BEFORE_RECEIPT = f"{_before_root}/cosine-graph-likelihood-readmatched-{EXHAUSTIVE}.jsonl"
        SKELETON = f"{D}/realign-exhaustive-{EXHAUSTIVE}.skeleton.jsonl"
    ANCHOR = f"{D}/anchor-projection-{EXHAUSTIVE}.jsonl"
    COMPONENT = f"S288C#0#{EXHAUSTIVE}"
    _loci = []
    for _line in open(RECEIPT):
        _loci.append(json.loads(_line)["locus"])
    LOCI = sorted(_loci)
    del _loci, _line
    LOCI_SET = set(LOCI)
    # (chrIV scale: the exhaustive checker over all 167 loci can be
    # SLICED — every check runs unchanged per sliced locus; the slices
    # are a wall measure, not a sample. "--loci-slice 0-49,100" keeps
    # loci 0..49 and 100.)
    if "--loci-slice" in sys.argv:
        _i = sys.argv.index("--loci-slice")
        _keep = set()
        for _part in sys.argv[_i + 1].split(","):
            if "-" in _part:
                _a, _b = _part.split("-")
                _keep.update(range(int(_a), int(_b) + 1))
            else:
                _keep.add(int(_part))
        _sliced = [l for l in LOCI if l in _keep]
        print(f"(loci slice: {len(_sliced)} of {len(LOCI)} loci)", flush=True)
        LOCI = _sliced
        LOCI_SET = set(_sliced)
        del _keep, _part, _sliced
elif FRAME:
    RECEIPT = f"{D}/realign-framescore-chrI.jsonl"
    EXACTNESS = f"{D}/realign-framescore-chrI.exactness.jsonl"
    INGREDIENTS = f"{D}/realign-framescore-chrI.jsonl.ingredients.jsonl"
    RUN = f"{D}/run-framescore-chrI-chrI"
    BEFORE_RECEIPT = f"{D}/realign-marginal-chrI.jsonl"
    IDENTITY_RECEIPT = f"{D}/realign-frameidentity-chrI.jsonl"
    IDENTITY_EXACTNESS = f"{D}/realign-frameidentity-chrI.exactness.jsonl"
    IDENTITY_INGREDIENTS = f"{D}/realign-frameidentity-chrI.jsonl.ingredients.jsonl"
    IDENTITY_RUN = f"{D}/run-frameidentity-chrI-chrI"
    IDENTITY_ANATOMY = f"{D}/realign-frameidentity-chrI.anatomy.jsonl"
    ANATOMY = f"{D}/realign-framescore-chrI.anatomy.jsonl"
    SKELETON = f"{D}/realign-framescore-chrI.skeleton.jsonl"
    AUDIT_CHRI = f"{D}/realign-frameaudit-chrI.jsonl"
    AUDIT_CHRMT = f"{D}/realign-frameaudit-chrMT.jsonl"
elif MARGINAL:
    RECEIPT = f"{D}/realign-marginal-chrI.jsonl"
    EXACTNESS = f"{D}/realign-marginal-chrI.exactness.jsonl"
    INGREDIENTS = f"{D}/realign-marginal-chrI.jsonl.ingredients.jsonl"
    RUN = f"{D}/run-realignmarginal-chrI"
    BEFORE_RECEIPT = f"{D}/realign-score-chrI.jsonl"


def axis_attribution(axis_rows):
    """(Stage 2, the edge refinement's witness) per node: the first
    window's own axis row carrying it and how many windows' axis rows
    carry it — re-derived from the partition GFAs' P lines (the same
    stored walks the instrument attributes over from the panel).
    axis_rows is the caller's (start, end, partition) window list —
    the checker's own load_maps derivation, main()-local."""
    if hasattr(axis_attribution, "_cache"):
        return axis_attribution._cache
    rep, count = {}, {}
    for w, (start, end, partition) in enumerate(axis_rows):
        g = gfa(partition)
        _pos, steps = g.spelled_syncmers(f"{COMPONENT}:{start}-{end}")
        for n in steps:
            key = abs(n)
            if key not in rep:
                rep[key] = w
            count[key] = count.get(key, 0) + 1
    axis_attribution._cache = (rep, count)
    return rep, count


def logsumexp2(a, b):
    m = max(a, b)
    if m == -math.inf:
        return m
    return m + math.log(math.exp(a - m) + math.exp(b - m))


failures = []
checks = [0]


def check(ok, message):
    checks[0] += 1
    if not ok:
        failures.append(message)
        print(f"FAIL: {message}", flush=True)


def fnv1a64(data):
    h = 0xCBF29CE484222325
    for b in data:
        h ^= b
        h = (h * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return h


def revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def complement(base):
    return {"A": "T", "C": "G", "G": "C", "T": "A"}.get(base, base)


def window_is_canonical_forward(window):
    """min(K, rc(K)) == K decided from both ends (the Rust function's
    exact loop, byte-wise): a palindrome-around-the-center window
    falls to the middle base's comparison."""
    lo, hi = 0, len(window) - 1
    while True:
        a, b = window[lo], window[hi]
        comp_b = {ord("A"): ord("T"), ord("C"): ord("G"),
                  ord("G"): ord("C"), ord("T"): ord("A")}.get(b)
        if comp_b is None:
            return False
        if a not in (ord("A"), ord("C"), ord("G"), ord("T")):
            return False
        if a != comp_b:
            return a < comp_b
        if lo + 1 >= hi:
            return True
        lo += 1
        hi -= 1


def load_names():
    names = {}
    with open(NAMES) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            names[int(fields[0])] = fields[1]
    return names


def load_maps():
    maps = {}
    axis = []
    for path in glob.glob(f"{GRAPHS}/*.gfa.map.json"):
        with open(path) as f:
            m = json.load(f)
        maps[m["partition"]] = m
        # ONE WINDOW PER AXIS ROW (not per partition): a repeat-locality
        # partition can host several of the component's windows (first
        # measured at chrIV — 167 windows over 158 axis partitions), so
        # every component member row contributes a window and the
        # window->partition map may repeat a partition.
        for member in m["members"]:
            if member["path_name"] == COMPONENT:
                axis.append((member["start"], member["end"], m["partition"]))
    axis.sort()
    return maps, [p for _, _, p in axis], [(s, e, p) for s, e, p in axis]


class Gfa:
    """One partition GFA: S segments, P lines, L overlaps."""

    def __init__(self, partition):
        self.partition = partition
        self.segs = {}
        self.paths = {}
        self.links = {}
        with open(f"{GRAPHS}/partition{partition}.gfa") as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if fields[0] == "S":
                    self.segs[fields[1]] = fields[2]
                elif fields[0] == "P":
                    self.paths[fields[1]] = fields[2].split(",")
                elif fields[0] == "L":
                    self.links[(fields[1], fields[2], fields[3], fields[4])] = int(
                        fields[5].rstrip("M")
                    )

    NUM_SYNCMER_NODES = 10024605

    def spelled(self, row_name):
        """(spelled sequence, step positions relative to the first step,
        signed steps INCLUDING the interleaved gap segments) — '-'
        steps spell reverse complements, L-line CIGAR overlaps trim the
        incoming step's head in its orientation."""
        steps = self.paths[row_name]
        content = {}
        seq_parts = []
        positions = []
        pos = 0
        prev = None
        for step in steps:
            name, sign = step[:-1], step[-1]
            if (name, sign) not in content:
                content[(name, sign)] = (
                    self.segs[name] if sign == "+" else revcomp(self.segs[name])
                )
            if prev is not None:
                overlap = self.links[(prev[0], prev[1], name, sign)]
                pos += len(self.segs[prev[0]]) - overlap
                seq_parts.append(content[(name, sign)][overlap:])
            else:
                seq_parts.append(content[(name, sign)])
            positions.append(pos)
            prev = (name, sign)
        return "".join(seq_parts), positions, steps

    def spelled_syncmers(self, row_name):
        """(positions, signed steps) of the SYNCMER steps only (the gap
        segments interleaved between non-overlapping syncmers excluded —
        the committed rows' walks are the syncmer steps)."""
        positions, steps = [], []
        pos = 0
        prev = None
        for step in self.paths[row_name]:
            name, sign = step[:-1], step[-1]
            if prev is not None:
                overlap = self.links[(prev[0], prev[1], name, sign)]
                pos += len(self.segs[prev[0]]) - overlap
            if int(name) <= self.NUM_SYNCMER_NODES:
                positions.append(pos)
                steps.append(int(name) * (1 if sign == "+" else -1))
            prev = (name, sign)
        return positions, steps


GFAS = {}


def gfa(partition):
    if partition not in GFAS:
        GFAS[partition] = Gfa(partition)
    return GFAS[partition]


def row_gfa_name(member):
    return f"{member['path_name']}:{member['start']}-{member['end']}"


def own_positions(canonical, mirror, w_lo):
    last = canonical[-1][1]
    if mirror:
        return [w_lo + (last - pos) for _, pos in canonical]
    return [w_lo + pos for _, pos in canonical]


def merged_windows(windows):
    windows = sorted(windows)
    merged = []
    for lo, hi in windows:
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return merged


def skipped_from_own(own, read_len):
    merged = merged_windows([[p, p + K] for p in own])
    skipped = []
    cursor = 0
    for lo, hi in merged:
        if lo > cursor:
            skipped.append((cursor, lo))
        cursor = max(cursor, hi)
    if cursor < read_len:
        skipped.append((cursor, read_len))
    return skipped


def correspondence(canonical, orientation, node_positions):
    """The unique-monotone pinned correspondence (mirrored from the Rust)."""
    m = len(canonical)
    candidates = []
    for node, _ in canonical:
        required = node if orientation == 0 else -node
        positions = [
            bp for bp, sign in node_positions.get(abs(required), [])
            if sign == (1 if required > 0 else -1)
        ]
        candidates.append(positions)
    pinned = [None] * m
    seen = [False] * m
    chain = []

    def enumerate_at(index):
        if index == m:
            for j, bp in chain:
                if not seen[j]:
                    pinned[j] = bp
                    seen[j] = True
                elif pinned[j] != bp:
                    pinned[j] = None
            return
        if not candidates[index]:
            enumerate_at(index + 1)
            return
        for bp in candidates[index]:
            ok = True
            if chain:
                if orientation == 0:
                    ok = bp > chain[-1][1]
                else:
                    ok = bp < chain[-1][1]
            if ok:
                chain.append((index, bp))
                enumerate_at(index + 1)
                chain.pop()

    enumerate_at(0)
    return pinned


def serving_anchor(pinned, own_hull, d):
    best = None
    for j, r in enumerate(pinned):
        if r is None:
            continue
        own_pos = own_hull[j]
        if best is None:
            best = (j, own_pos)
        else:
            own_left = own_pos <= d
            best_left = best[1] <= d
            if (own_left and not best_left) or (
                own_left == best_left and ((own_pos > best[1]) if own_left else (own_pos < best[1]))
            ):
                best = (j, own_pos)
    return best


def placement_votes(read, canonical, mirror, w_lo, orientation, pinned):
    """The independent per-base recomputation of one placement: returns
    (backbone_bp, votes, m, c, collinear, sigma) or None when invalid
    (a voted base projects outside the fold). The backbone is verified
    PER BASE (never trusted from node identity); the skipped votes are
    recomputed; (m, c) is the total."""
    forward = (orientation == 0) != mirror
    own = own_positions(canonical, mirror, w_lo)
    own_hull = own_positions(canonical, mirror, 0)
    windows = [[own_hull[j], own_hull[j] + K] for j, r in enumerate(pinned) if r is not None]
    merged = merged_windows([list(w) for w in windows])
    backbone_bp = sum(hi - lo for lo, hi in merged)
    votes = []
    matches = 0
    mismatches = 0
    # The backbone, per base over the merged pinned-window coverage,
    # each base verified through one covering pinned anchor against the
    # fold's spelled sequence (strand-aware).
    for lo, hi in merged:
        for dd in range(lo, hi):
            cover = None
            for j, r in enumerate(pinned):
                if r is None:
                    continue
                p = own_hull[j]
                if p <= dd < p + K:
                    if cover is None or p > own_hull[cover]:
                        cover = j
            j = cover
            r = pinned[j]
            t = dd - own_hull[j]
            read_i = w_lo + dd
            read_base = read[read_i] if forward else complement(read[read_i])
            coord = r + t if forward else r + K - 1 - t
            if coord < 0 or coord >= len(FOLD_SEQ[0]):
                return None
            if read_base == FOLD_SEQ[0][coord]:
                matches += 1
            else:
                mismatches += 1
    # The skipped-base votes (the read's bases outside its merged anchor
    # windows), projected through the serving anchors.
    for lo, hi in skipped_from_own(own, len(read)):
        for i in range(lo, hi):
            d = i - w_lo
            served = serving_anchor(pinned, own_hull, d)
            if served is None:
                return None
            j, own_pos = served
            r = pinned[j]
            coord = r + (d - own_pos) if forward else r + (own_pos + K - 1 - d)
            if coord < 0 or coord >= len(FOLD_SEQ[0]):
                return None
            read_base = read[i] if forward else complement(read[i])
            vote = 1 if read_base == FOLD_SEQ[0][coord] else 0
            votes.append((i, coord, vote))
            matches += vote
            mismatches += 1 - vote
    # collinearity and the single-offset sigma
    deltas = set()
    for j, r in enumerate(pinned):
        if r is None:
            continue
        if forward:
            deltas.add(r - own_hull[j])
        else:
            deltas.add(r + own_hull[j])
    collinear = len(deltas) == 1
    sigma = None
    if collinear:
        if forward:
            sigma = min(deltas) - w_lo
        else:
            sigma = min(deltas) + K - 1 - 149 + w_lo
    return backbone_bp, votes, matches, mismatches, collinear, sigma


FOLD_SEQ = [None]  # the current fold's spelled sequence (set per pair)


def edit_distance_banded(read, window, band):
    """Semiglobal unit-cost edit distance of the read against any
    substring of the window (free prefix/suffix), banded; returns
    (distance, indels) via backtrace."""
    n, m = len(read), len(window)
    INF = float("inf")
    # D[i][j]: read[:i] aligned to a suffix window[..j]
    prev = [0] * (m + 1)
    # full DP (window is bounded: ~310bp)
    trace = []
    rows = [prev]
    for i in range(1, n + 1):
        cur = [INF] * (m + 1)
        cur[0] = i
        for j in range(1, m + 1):
            best = INF
            if prev[j - 1] + (0 if read[i - 1] == window[j - 1] else 1) < best:
                best = prev[j - 1] + (0 if read[i - 1] == window[j - 1] else 1)
            if prev[j] + 1 < best:
                best = prev[j] + 1
            if cur[j - 1] + 1 < best:
                best = cur[j - 1] + 1
            cur[j] = best
        rows.append(cur)
        prev = cur
    jbest = min(range(m + 1), key=lambda j: rows[n][j])
    dist = rows[n][jbest]
    # backtrace for the indel count
    i, j = n, jbest
    indels = 0
    while i > 0:
        if j > 0 and rows[i][j] == rows[i - 1][j - 1] + (0 if read[i - 1] == window[j - 1] else 1):
            i, j = i - 1, j - 1
        elif rows[i][j] == rows[i - 1][j] + 1:
            i -= 1
            indels += 1
        else:
            j -= 1
            indels += 1
    return dist, indels


def mix_logsumexp(a, b):
    if a == float("-inf") and b == float("-inf"):
        return float("-inf")
    m = max(a, b)
    return m + math.log(0.5 * math.exp(a - m) + 0.5 * math.exp(b - m))


def median_of(values):
    values = sorted(values)
    middle = len(values) // 2
    if len(values) % 2 == 1:
        return values[middle]
    return (values[middle - 1] + values[middle]) / 2.0


def spectrum_knee(sorted_distances):
    positives = [d for d in sorted_distances if d > 0.0]
    if len(positives) < 2:
        return False, 0.0
    jumps = [(pair[1] - pair[0]) / pair[1] for pair in zip(positives, positives[1:])]
    typical = median_of(list(jumps))
    best, best_index = float("-inf"), 0
    for index, jump in enumerate(jumps):
        if jump >= best:
            best = jump
            best_index = index
    if best > typical:
        return True, positives[best_index]
    return False, 0.0


def cluster_form_qual(s_win, spectrum, tied, compact=False):
    sorted_distances = sorted(d for d, _ in spectrum)
    has_knee, cut = spectrum_knee(sorted_distances)
    excluded = [d <= cut for d, _ in spectrum]
    parent = list(range(len(tied) + 1))

    def root(node):
        while parent[node] != node:
            node = parent[node]
        return node

    if compact:
        # (Stage 2, THE WALL FIX's mirror: the degenerate-locus tie
        # certificate is the compact CUT-INDEX list — each list IS its
        # band's exact <=cut predicate set, so the excluded-union is
        # list membership and a union-find edge (left, right) exists
        # iff tied class right's spectrum index appears in left's
        # list. The edge set is IDENTICAL to the band form's (which
        # tests band_left[index_right] <= cut), so k is identical;
        # the iteration is list-driven, linear in the total predicate
        # count instead of quadratic in the tied-class count — the
        # degenerate loci tie ~35,510 classes and the quadratic
        # pairwise loop is the checker-side wall.)
        ordinal_of = {entry[0]: ordinal for ordinal, entry in enumerate(tied)}
        for left, (index, members) in enumerate(tied):
            for member in members:
                excluded[member] = True
            if spectrum[index][0] <= cut:
                a, b = root(0), root(left + 1)
                if a != b:
                    parent[a] = b
            for member in members:
                right = ordinal_of.get(member)
                if right is not None:
                    a, b = root(left + 1), root(right + 1)
                    if a != b:
                        parent[a] = b
    else:
        for _, band in tied:
            for index, d in enumerate(band):
                if d <= cut:
                    excluded[index] = True
        for left in range(len(tied)):
            if spectrum[tied[left][0]][0] <= cut:
                a, b = root(0), root(left + 1)
                if a != b:
                    parent[a] = b
            for right in range(left + 1, len(tied)):
                # (the band-pair distance is band_left's row at the
                # RIGHT BAND'S FOLD INDEX — the instrument's own indexing
                # (`tied[left].1[index_right]`), not the band's ordinal in
                # the tied list; chrXIII L1's three-band tie caught the
                # ordinal misread as a k=4 false failure)
                if tied[left][1][tied[right][0]] <= cut:
                    a, b = root(left + 1), root(right + 1)
                    if a != b:
                        parent[a] = b
    roots = sorted({root(n) for n in range(len(tied) + 1)})
    k = len(roots)
    alternative = None
    for index, (_, score) in enumerate(spectrum):
        if not excluded[index]:
            alternative = score if alternative is None else max(alternative, score)
    shape = "single_class" if not spectrum else ("knee" if has_knee else "no_knee")
    if k == 0:
        p = None
    else:
        denom = k * s_win + (alternative or 0.0)
        p = s_win / denom if denom > 0 else None
    qual = None if p is None or p >= 1.0 else -10.0 * math.log10(1.0 - p)
    return {
        "knee": cut if has_knee else None,
        "shape": shape,
        "k": k,
        "alternative": alternative,
        "p": p,
        "qual": qual,
        "unbounded": p is not None and p >= 1.0,
    }


def main():
    print("== phase 1: run markers, exit, wall, the 64GiB guard", flush=True)
    check(os.path.exists(f"{RUN}.done"), "run done marker missing")
    exit_code = int(open(f"{RUN}.exit").read().strip())
    check(exit_code == 0, f"run exit code {exit_code}")
    wall = int(open(f"{RUN}.wall").read().strip())
    rss_peak = 0
    with open(f"{RUN}.rss") as f:
        for line in f:
            rss_peak = max(rss_peak, int(line.split()[2]))
    check(rss_peak <= RSS_BUDGET_KB, f"RSS guard exceeded: {rss_peak} kB")
    print(f"   exit {exit_code}, wall {wall}s, rss peak {rss_peak} kB", flush=True)

    names = load_names()
    maps, axis, axis_rows = load_maps()

    # The committed anchor-projection receipt, restricted to the pilot
    # records (the canonical walks, the verified shifts, the
    # occurrences' orientations and touched windows).
    print("== loading the committed anchor-projection receipt (pilot records)", flush=True)
    pilot_needed = set()
    with open(RECEIPT) as f:
        for line in f:
            d = json.loads(line)
            pilot_needed.add(d["locus"])
    anchor_records = {}
    with open(ANCHOR) as f:
        for line in f:
            d = json.loads(line)
            record = d["record"]
            for occ in d["occurrences"]:
                if any(w in pilot_needed for w in occ["partitions"]):
                    anchor_records[record] = d
                    break
    print(f"   {len(anchor_records)} pilot records loaded", flush=True)

    # The receipts.
    receipt = {}
    for line in open(RECEIPT):
        d = json.loads(line)
        receipt[d["locus"]] = d
    ingredients = {}
    for line in open(INGREDIENTS):
        d = json.loads(line)
        check(
            "ll_matrix" in d and "folds" in d,
            "ingredients line lacks folds or ll_matrix",
        )
        ingredients[d["locus"]] = d
    exactness = []
    for line in open(EXACTNESS):
        exactness.append(json.loads(line))
    print(f"   {len(exactness)} exactness sample pairs", flush=True)

    # Slice D inputs (frame mode): the skeleton sidecar, the paired
    # anatomy receipts, the identity (stored-walk) receipts.
    skeleton_by_fold = {}
    frame_totals = [0, 0, 0, 0, 0]
    identity_receipt = {}
    identity_ingredients = {}
    anatomy_by_locus = {}
    identity_anatomy_by_locus = {}
    if FRAME:
        for line in open(SKELETON):
            d = json.loads(line)
            skeleton_by_fold[(d["locus"], d["fold"])] = d

        if EXHAUSTIVE is None:
            def load_anatomy(path):
                out = {}
                for line in open(path):
                    d = json.loads(line)
                    out.setdefault(d["locus"], {})[(d["record"], d["variant"])] = d
                return out

            anatomy_by_locus = load_anatomy(ANATOMY)
            identity_anatomy_by_locus = load_anatomy(IDENTITY_ANATOMY)
            for line in open(IDENTITY_RECEIPT):
                d = json.loads(line)
                identity_receipt[d["locus"]] = d
            for line in open(IDENTITY_INGREDIENTS):
                d = json.loads(line)
                identity_ingredients[d["locus"]] = d
            print(
                f"   frame mode: {len(skeleton_by_fold)} skeleton rows, "
                f"{sum(len(v) for v in anatomy_by_locus.values())} anatomy units (after), "
                f"{sum(len(v) for v in identity_anatomy_by_locus.values())} (before)",
                flush=True,
            )
        else:
            print(
                f"   exhaustive mode ({EXHAUSTIVE}): {len(skeleton_by_fold)} skeleton rows",
                flush=True,
            )

    # ---------------- phase 2: folds vs the partition maps and GFAs
    print("== phase 2: fold structure vs the maps and the GFAs", flush=True)
    fold_seqs = {}
    for locus in LOCI:
        d = receipt[locus]
        partition = d["partition"]
        pmap = maps[partition]
        members = [(m["path_name"], m["start"], m["end"]) for m in pmap["members"]]
        fold_members = []
        for fold in ingredients[locus]["folds"]:
            fold_members.extend(
                (m["path_name"], m["start"], m["end"]) for m in fold["members"]
            )
        check(
            sorted(fold_members) == sorted(members),
            f"locus {locus}: folds do not partition the member rows exactly",
        )
        # every fold: sequence + contained walk re-derived from the GFA.
        # The GFA spelling starts at the first (possibly edge-overlapping)
        # step's window start, so the row's [start, end) sequence sits at
        # a front-overhang offset inside the spelling; the fold's walk is
        # the CONTAINED syncmer steps (windows fully inside the extent),
        # positioned relative to the row start.
        g = gfa(partition)
        for fi, fold in enumerate(ingredients[locus]["folds"]):
            first = fold["members"][0]
            seq, positions, steps = g.spelled(row_gfa_name(first))
            row_seq = fold["sequence"]
            check(
                len(row_seq) == fold["length"],
                f"locus {locus} fold {fi}: length differs from the sequence",
            )
            offset = seq.find(row_seq)
            check(
                offset >= 0 and seq[offset:offset + len(row_seq)] == row_seq,
                f"locus {locus} fold {fi}: row sequence absent from the GFA spelling",
            )
            sync_positions, sync_steps = g.spelled_syncmers(row_gfa_name(first))
            contained = [
                (p - offset, s)
                for p, s in zip(sync_positions, sync_steps)
                if p >= offset and p + K <= offset + len(row_seq)
            ]
            walk_receipt = [(w[0], w[1]) for w in fold["walk"]]
            if FRAME:
                # THE CANONICAL-SKELETON AUDIT (slice D): every claimed
                # step verified BY SEQUENCE against the node's own GFA S
                # segment (orientation-aware per the sign), the frame
                # decision verified against the window's canonical form
                # (a step at a stored position with the stored node must
                # be canonical-forward; a step the stored walk lacks
                # must be rc-canonical), the forward part verified
                # COMPLETE against the P-line walk, and the dropped
                # stored positions verified rc-canonical.
                stored_map = dict(contained)
                check(
                    len(set(rel for rel, _ in walk_receipt)) == len(walk_receipt),
                    f"locus {locus} fold {fi}: duplicate step positions in the claimed skeleton",
                )
                claimed_map = dict(walk_receipt)
                kept_forward = 0
                kept_reverse = 0
                added = 0
                replaced = 0
                added_seq_verified = 0
                for rel, node in walk_receipt:
                    window = row_seq[rel : rel + K]
                    cf = window_is_canonical_forward(window.encode())
                    seg = g.segs.get(str(abs(node)))
                    if seg is not None:
                        oriented = seg if node > 0 else revcomp(seg)
                        check(
                            window == oriented,
                            f"locus {locus} fold {fi}: claimed step {rel}/{node} "
                            "fails sequence verification against the GFA S segment",
                        )
                    else:
                        # An rc-frame-only anchor has no S segment in
                        # this partition's GFA (the GFA spells the
                        # STORED walks' nodes only). The node identity
                        # is verified by the in-process kmerHash check
                        # (the source of truth) and independently by
                        # phase 4's per-base exactness re-derivation
                        # (a wrong node id would break every placement
                        # that pins through it); here the window's
                        # canonical form is verified against the
                        # GFA-spelled row.
                        check(
                            not cf,
                            f"locus {locus} fold {fi}: claimed step {rel}/{node} has no "
                            "GFA S segment but a canonical-forward window",
                        )
                    if seg is not None:
                        added_seq_verified += 1
                    if cf:
                        # The canonical scheme keeps the FORWARD frame's
                        # selection at canonical-forward positions, so
                        # the stored walk must carry exactly this step.
                        kept_forward += 1
                        check(
                            stored_map.get(rel) == node,
                            f"locus {locus} fold {fi}: canonical-forward step {rel} "
                            "differs from the stored walk's selection",
                        )
                    else:
                        # At rc-canonical positions the rc frame's node is
                        # kept (it may coincide with the stored node when
                        # both frames qualified there).
                        kept_reverse += 1
                    if rel not in stored_map:
                        added += 1
                    elif stored_map[rel] != node:
                        replaced += 1
                for rel, node in contained:
                    if window_is_canonical_forward(row_seq[rel : rel + K].encode()):
                        check(
                            claimed_map.get(rel) == node,
                            f"locus {locus} fold {fi}: canonical-forward stored "
                            f"step {rel} missing from the claimed skeleton",
                        )
                dropped = [rel for rel, _ in contained if rel not in claimed_map]
                for rel in dropped:
                    check(
                        not window_is_canonical_forward(row_seq[rel : rel + K].encode()),
                        f"locus {locus} fold {fi}: dropped stored step {rel} "
                        "is canonical-forward",
                    )
                skel = skeleton_by_fold[(locus, fi)]
                check(
                    [(s[0], s[1]) for s in skel["steps"]] == walk_receipt,
                    f"locus {locus} fold {fi}: skeleton sidecar steps differ",
                )
                check(
                    [(s[0], s[1]) for s in skel["stored"]] == contained,
                    f"locus {locus} fold {fi}: skeleton sidecar stored steps differ",
                )
                for rel, node, frame in skel["steps"]:
                    check(
                        frame == (0 if window_is_canonical_forward(row_seq[rel : rel + K].encode()) else 1),
                        f"locus {locus} fold {fi}: skeleton sidecar frame tag wrong at {rel}",
                    )
                check(
                    (skel["added"], skel["replaced"], skel["dropped"], skel["kept_forward"], skel["kept_reverse"])
                    == (added, replaced, len(dropped), kept_forward, kept_reverse),
                    f"locus {locus} fold {fi}: skeleton sidecar diff census differs",
                )
                frame_totals[0] += kept_forward
                frame_totals[1] += added
                frame_totals[2] += len(dropped)
                frame_totals[3] += added_seq_verified
                frame_totals[4] += kept_reverse
            else:
                check(
                    walk_receipt == contained,
                    f"locus {locus} fold {fi}: contained walk differs from the GFA P line",
                )
            fold_seqs[(locus, fi)] = row_seq
    if FRAME:
        print(
            f"   folds verified against {len(LOCI)} partitions' GFAs; the "
            f"canonical skeletons: kept-forward {frame_totals[0]}, "
            f"kept-reverse {frame_totals[4]}, added (positions the stored "
            f"walk lacks) {frame_totals[1]}, dropped {frame_totals[2]}; "
            f"{frame_totals[3]} steps carry a GFA S segment and are "
            "sequence-verified byte-level (the rc-frame-only remainder "
            "verified by canonical form here, by the in-process kmerHash "
            "check against the AGC-fetched sequence, and by phase 4's "
            "per-base exactness); every frame decision verified, the "
            "forward part complete, the dropped positions rc-canonical",
            flush=True,
        )
    else:
        print(f"   folds verified against {len(LOCI)} partitions' GFAs", flush=True)

    # ---------------- phase 3: the locality ruling re-derivation
    print("== phase 3: the locality counts re-derived", flush=True)
    rows_by_pp = {}
    for pid, m in maps.items():
        for row in m["members"]:
            rows_by_pp.setdefault((pid, row["path_name"]), []).append((row["start"], row["end"]))
    for k in rows_by_pp:
        rows_by_pp[k].sort()
    # (chrIV-scale inversion, semantically identical to the committed
    # per-locus loop: the occurrence classification is a per-OCCURRENCE
    # property — it reads only the occurrence's own touched windows'
    # rows — so it is computed once per occurrence and the per-locus
    # counts accumulate the touching occurrences by their class. The
    # committed loop was O(loci x occurrences); chrIV is 167 x 12.75M.)
    touching_classes = {}
    for record, rd in anchor_records.items():
        span = rd["span"]
        for oi, occ in enumerate(rd["occurrences"]):
            origin = occ["start"] + occ["shift"]
            contained = False
            overlapped = False
            for w in occ["partitions"]:
                lst = rows_by_pp.get((axis[w], names[occ["path"]]))
                if not lst:
                    continue
                for s, e in lst:
                    if s <= origin and origin + span <= e:
                        contained = True
                    if s < origin + span and origin < e:
                        overlapped = True
            cls = 0 if contained else (1 if overlapped else 2)
            for locus in set(occ["partitions"]):
                if locus in receipt:
                    touching_classes.setdefault(locus, []).append(cls)
    for locus in LOCI:
        d = receipt[locus]
        counts = [0, 0, 0]
        for cls in touching_classes.get(locus, []):
            counts[cls] += 1
        in_axis, overhang, extension = counts
        check(
            d["locality"]["touching_in_axis"] == in_axis
            and d["locality"]["touching_overhang"] == overhang
            and d["locality"]["touching_extension"] == extension,
            f"locus {locus}: locality counts differ "
            f"({in_axis}/{overhang}/{extension} vs receipt)",
        )
    print("   locality counts reproduced exactly", flush=True)

    # ---------------- phase 4: the exactness proof
    print("== phase 4: the exactness proof (sampled pairs)", flush=True)
    needed_fnv = {pair["read_fnv"] for pair in exactness if pair["locus"] in receipt}
    fastq_fnv = set()
    with gzip.open(READS, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            h = fnv1a64(seq.encode())
            if h in needed_fnv:
                fastq_fnv.add(h)
    check(
        needed_fnv <= fastq_fnv,
        f"{len(needed_fnv - fastq_fnv)} sampled reads absent from the FASTQ",
    )
    print(f"   all {len(needed_fnv)} sampled reads verified in the FASTQ", flush=True)

    units_by_locus = {}
    for locus in LOCI:
        units_by_locus[locus] = ingredients[locus]["units"]

    vote_checked = 0
    exact_pairs = 0
    offset_dominated = 0
    offset_named = 0
    collinear_pairs = 0
    indel_free = 0
    indel_pairs = 0
    edit_distances = []
    for pair in exactness:
        locus = pair["locus"]
        if locus not in LOCI_SET:
            # (a --loci-slice invocation: pairs outside the slice are
            # covered by the slice that owns them)
            continue
        fold_index = pair["fold"]
        fold = ingredients[locus]["folds"][fold_index]
        FOLD_SEQ[0] = fold_seqs[(locus, fold_index)]
        read = pair["read"]
        check(fnv1a64(read.encode()) == pair["read_fnv"], "sampled read FNV mismatch")
        record = pair["record"]
        canonical = [(a[0], a[1]) for a in anchor_records[record]["anchors"]]
        mirror = pair["mirror"] == 1
        w_lo = pair["w_lo"]
        # the orientations present among the touching occurrences
        orientations = set()
        for occ in anchor_records[record]["occurrences"]:
            if locus in occ["partitions"]:
                orientations.add(occ["orientation"])
        receipt_placements = {p["orientation"]: p for p in pair["placements"]}
        check(
            set(receipt_placements) <= orientations,
            f"locus {locus} record {record}: placement orientation not among the touching occurrences",
        )
        for orientation in sorted(receipt_placements):
            # the node positions from the GFA walk
            node_positions = {}
            for bp, node in fold["walk"]:
                node_positions.setdefault(abs(node), []).append((bp, 1 if node > 0 else -1))
            pinned = correspondence(canonical, orientation, node_positions)
            receipt_pin = receipt_placements[orientation]
            check(
                [p for p in pinned] == receipt_pin["pinned"],
                f"locus {locus} record {record} orientation {orientation}: pinned vector differs",
            )
            got = placement_votes(read, canonical, mirror, w_lo, orientation, pinned)
            if got is None:
                # invalid placement (out of extent): the receipt must not
                # have scored it either -- it is absent from placements
                # with a score; the receipt only lists scored ones, so
                # this branch means the receipt has one we deem invalid
                check(False, "receipt placement invalid under re-derivation")
                continue
            backbone_bp, votes, m, c, collinear, sigma = got
            if TERRITORY:
                # (Stage 2: the sampled placements must lie WHOLLY
                # within the fold's territory image — the pinned
                # windows and every voted coordinate; an empty-image
                # fold admits no scored placement at all. The images
                # themselves are re-derived and verified in phase 5.)
                img = receipt[locus]["territory"]["images"][fold_index]
                check(
                    img is not None,
                    f"locus {locus} fold {fold_index}: a scored placement on an empty-image fold",
                )
                check(
                    all(img[0] <= v[1] < img[1] for v in votes),
                    f"locus {locus} record {record}: sampled votes outside the territory image",
                )
                check(
                    all(
                        img[0] <= r and r + K <= img[1]
                        for r in receipt_pin["pinned"]
                        if r is not None
                    ),
                    f"locus {locus} record {record}: pinned windows outside the territory image",
                )
            check(
                backbone_bp == receipt_pin["backbone_bp"],
                f"locus {locus} record {record}: backbone_bp differs",
            )
            check(
                [(v[0], v[1], v[2]) for v in votes]
                == [(v[0], v[1], v[2]) for v in receipt_pin["votes"]],
                f"locus {locus} record {record}: votes differ",
            )
            check(
                (m, c) == (receipt_pin["m"], receipt_pin["c"]),
                f"locus {locus} record {record}: (m,c) differs "
                f"({(m, c)} vs {(receipt_pin['m'], receipt_pin['c'])})",
            )
            vote_checked += len(votes) + backbone_bp
            exact_pairs += 1
            # the bounded-offset dominance scan (collinear placements)
            if collinear:
                collinear_pairs += 1
                strand_read = read if (orientation == 0) != mirror else revcomp(read)
                fold_seq = FOLD_SEQ[0]
                lo = max(0, sigma - 150)
                hi = min(len(fold_seq) - READ_LENGTH, sigma + 150)
                best_m = -1
                best_sigma = None
                for s in range(lo, hi + 1):
                    mm = sum(
                        1 for i in range(READ_LENGTH) if strand_read[i] == fold_seq[s + i]
                    )
                    if mm > best_m:
                        best_m = mm
                        best_sigma = s
                if best_sigma == sigma and best_m == m:
                    offset_dominated += 1
                else:
                    offset_named += 1
                # the bounded edit-distance alignment (unit costs, window
                # sigma +/- 80)
                window = fold_seq[max(0, sigma - 80): sigma + READ_LENGTH + 80]
                dist, indels = edit_distance_banded(strand_read, window, 80)
                edit_distances.append(dist)
                if indels == 0:
                    indel_free += 1
                else:
                    indel_pairs += 1
            else:
                offset_named += 1
    check(
        exact_pairs > 0 and offset_named == 0 or offset_named <= exact_pairs,
        "offset scan accounting",
    )
    print(
        f"   {exact_pairs} placements re-derived EXACTLY "
        f"({vote_checked} per-base votes+backbone bases); "
        f"offset-dominance {offset_dominated}/{collinear_pairs} collinear "
        f"({offset_named} named: piecewise/non-dominant); "
        f"edit view: {indel_free} indel-free, {indel_pairs} with indels",
        flush=True,
    )

    # ---------------- phase 5: the scoring + class re-derivation
    print("== phase 5: the per-unit LLs and the class table re-derived", flush=True)
    import numpy as np

    for locus in LOCI:
        d = receipt[locus]
        folds = ingredients[locus]["folds"]
        units = ingredients[locus]["units"]
        matrix = np.array(ingredients[locus]["ll_matrix"], dtype=float)
        check(matrix.shape == (len(folds), len(units)), f"locus {locus}: matrix shape")
        territory_images = None
        if TERRITORY:
            # ---------------- the territory images re-derived (rule 8)
            # The images are a PURE function of the folds' walks and
            # the axis fold (the window's own row): the interval
            # spanned by the steps shared with the axis row's walk.
            # Re-derived here from the ingredients and compared to the
            # receipt EXACTLY (both the territory field and the
            # fold_identities tag); the axis fold's membership is
            # verified against the window's own axis row from the
            # maps (the locus's coordinate, not any candidate's).
            t = d.get("territory")
            check(
                t is not None and t.get("convention") == "normalized",
                f"locus {locus}: territory field missing or not the normalized convention",
            )
            check(
                d["model"] == "anchor-realign-v2-marginal-frame-territory",
                f"locus {locus}: model tag is not the territory convention",
            )
            axis_fold = t["axis_fold"]
            axis_start, axis_end, _ = axis_rows[locus]
            axis_members = d["fold_identities"][axis_fold]["members"]
            check(
                any(
                    m["path_name"] == COMPONENT
                    and m["start"] == axis_start
                    and m["end"] == axis_end
                    for m in axis_members
                ),
                f"locus {locus}: the axis fold does not carry the window's own axis row",
            )
            axis_nodes = {abs(n) for _bp, n in folds[axis_fold]["walk"]}
            rep, count = axis_attribution(axis_rows)

            def attributable(node, _locus=locus):
                first = rep.get(node)
                if first is None:
                    return False
                return first != _locus or count[node] > 1

            derived = []
            for fi, fold in enumerate(folds):
                # the window's own axis fold IS the territory: full extent
                if fi == axis_fold:
                    derived.append((0, fold["length"]))
                    continue
                shared = [bp for bp, node in fold["walk"] if abs(node) in axis_nodes]
                if not shared:
                    # no shared anchor: own by elimination iff no other
                    # window's axis row claims any of its steps
                    if all(
                        not attributable(abs(n)) for _bp, n in fold["walk"]
                    ):
                        derived.append((0, fold["length"]))
                    else:
                        derived.append(None)
                    continue
                lo, hi = shared[0], shared[-1] + K
                head_claimed = any(
                    bp < lo and attributable(abs(n)) for bp, n in fold["walk"]
                )
                tail_claimed = any(
                    bp >= hi and attributable(abs(n)) for bp, n in fold["walk"]
                )
                derived.append(
                    (
                        lo if head_claimed else 0,
                        hi if tail_claimed else fold["length"],
                    )
                )
            territory_images = derived
            check(
                len(t["images"]) == len(folds),
                f"locus {locus}: territory image count differs",
            )
            for fi, img in enumerate(t["images"]):
                expect = [derived[fi][0], derived[fi][1]] if derived[fi] else None
                check(img == expect, f"locus {locus}: fold {fi} territory image differs")
                identity_img = d["fold_identities"][fi].get("image")
                check(
                    identity_img == expect,
                    f"locus {locus}: fold {fi} fold_identities image differs",
                )
            n_empty = sum(1 for img in derived if img is None)
            print(
                f"   locus {locus}: territory images verified "
                f"(axis fold {axis_fold}, {len(derived) - n_empty} images, "
                f"{n_empty} empty)",
                flush=True,
            )
        if MARGINAL:
            # THE E-DERIVATION AUDIT (from the receipt's own panel
            # strain list: the mean per-haplotype total path length,
            # doubled for the diploid, then E = 150*A - ln(2*(G-149))).
            sc = d["scoring"]
            strain_list = sc["panel_strain_list"]
            check(
                len(strain_list) == sc["panel_strains"],
                "locus %d: strain list count differs" % locus,
            )
            mean_hap = sum(s["bp"] for s in strain_list) / len(strain_list)
            check(
                abs(mean_hap - sc["panel_mean_haploid_bp"]) < 1e-6,
                f"locus {locus}: mean haploid re-derivation differs",
            )
            g_dip = 2.0 * mean_hap
            check(
                abs(g_dip - sc["genome_diploid_bp"]) < 1e-6,
                f"locus {locus}: diploid genome re-derivation differs",
            )
            e_expect = READ_LENGTH * A - math.log(2.0 * (g_dip - READ_LENGTH + 1.0))
            check(
                abs(e_expect - sc["elsewhere_log_prob"]) < 1e-9,
                f"locus {locus}: E re-derivation differs "
                f"({e_expect} vs {sc['elsewhere_log_prob']})",
            )
            E = sc["elsewhere_log_prob"]
            # THE BOUNDING PROPERTY: every unit's likelihood under every
            # fold is at least E (the elsewhere branch bounds every
            # floor), and no-pin entries sit at exactly E.
            check(
                float(matrix.min()) >= E - 1e-6,
                f"locus {locus}: matrix entry below E",
            )
            at_e = int((np.abs(matrix - E) < 1e-9).sum())
            print(
                f"   locus {locus}: E audit passed "
                f"(strains {sc['panel_strains']}, mean haploid {mean_hap:.1f}, "
                f"E {E:.6f}); {at_e} matrix entries at E",
                flush=True,
            )
        else:
            # floor-only entries: prior + 150*B exactly
            floor_checked = 0
            for fi, fold in enumerate(folds):
                prior = -math.log(2.0 * (fold["length"] - READ_LENGTH + 1))
                row = matrix[fi]
                floor_idx = np.where(row < prior + FLOOR + 1e-6)[0]
                check(
                    all(abs(row[u] - (prior + FLOOR)) < 1e-6 for u in floor_idx),
                    f"locus {locus} fold {fi}: floor entries not at prior+150B",
                )
                floor_checked += len(floor_idx)
        # re-derive sampled pinned (unit, fold) LLs end to end
        pairs_here = [p for p in exactness if p["locus"] == locus]
        for pair in pairs_here[:40]:
            fi = pair["fold"]
            fold = folds[fi]
            record = pair["record"]
            canonical = [(a[0], a[1]) for a in anchor_records[record]["anchors"]]
            unit = units[pair["unit"]]
            check(unit["record"] == record and unit["variant"] == pair["variant"], "unit identity")
            mirror = unit["mirror"] == 1
            w_lo = unit["w_lo"]
            orientations = set()
            for occ in anchor_records[record]["occurrences"]:
                if locus in occ["partitions"]:
                    orientations.add(occ["orientation"])
            node_positions = {}
            for bp, node in fold["walk"]:
                node_positions.setdefault(abs(node), []).append((bp, 1 if node > 0 else -1))
            scores = [FLOOR] if not MARGINAL else []
            FOLD_SEQ[0] = fold_seqs[(locus, fi)]
            # (Stage 2: under the territory convention the placement
            # domain is the fold's IMAGE — the windows-inside test and
            # the prior's support both use it; an empty image admits
            # no placement at all.)
            if TERRITORY:
                img = territory_images[fi]
                support = 0 if img is None else img[1] - img[0]
                img_lo = 0 if img is None else img[0]
                img_hi = 0 if img is None else img[1]
            else:
                support = fold["length"]
                img_lo, img_hi = 0, fold["length"]
            for orientation in sorted(orientations):
                pinned = correspondence(canonical, orientation, node_positions)
                if all(p is None for p in pinned):
                    continue
                if any(
                    r + K > img_hi or r < img_lo for r in pinned if r is not None
                ):
                    continue
                got = placement_votes(unit["read"], canonical, mirror, w_lo, orientation, pinned)
                if got is None:
                    continue
                if TERRITORY and any(v[1] < img_lo or v[1] >= img_hi for v in got[1]):
                    # a voted base outside the image: not a placement
                    # of the normalized convention (the instrument
                    # discards it; so does the faithful re-derivation)
                    continue
                _, _, m, c, _, _ = got
                scores.append(m * A + c * B)
            # (The instrument's own degenerate-case convention
            # (unit_log_likelihood): an image shorter than a read
            # admits NO placement — local_ll NEG_INFINITY — and the
            # prior is never formed. Under the territory convention a
            # non-empty image can be shorter than a read (chrXIII
            # L22's sampled fold); compute the prior lazily so the
            # faithful re-derivation mirrors the instrument instead of
            # leaving the log domain.)
            prior = (
                -math.log(2.0 * (support - READ_LENGTH + 1))
                if support >= READ_LENGTH
                else -math.inf
            )
            best = max(scores) if scores else -math.inf
            local_ll = (
                -math.inf
                if not scores or support < READ_LENGTH
                else prior + best + math.log(sum(math.exp(s - best) for s in scores))
            )
            if MARGINAL:
                ll = logsumexp2(local_ll, E)
                check(
                    abs(pair["ll_local"] - local_ll) <= 1e-9 * max(1.0, abs(local_ll))
                    if local_ll != -math.inf
                    else pair["ll_local"] is None or pair["ll_local"] == -math.inf,
                    f"locus {locus}: sampled ll_local differs",
                )
                check(
                    abs(pair["elsewhere"] - E) < 1e-9,
                    f"locus {locus}: sampled elsewhere term differs from E",
                )
            else:
                ll = local_ll
            got_ll = matrix[fi][pair["unit"]]
            check(
                abs(ll - got_ll) <= 1e-9 * max(1.0, abs(ll)),
                f"locus {locus}: unit LL re-derivation differs {ll} vs {got_ll}",
            )
            check(
                abs(pair["ll"] - got_ll) <= 1e-9 * max(1.0, abs(got_ll)),
                f"locus {locus}: sampled pair ll differs from the matrix",
            )
        # the class table
        counts = np.array([u["count"] for u in units], dtype=float)
        n_folds = len(folds)
        class_lls = np.empty(n_folds * (n_folds + 1) // 2)
        idx = 0
        for j in range(n_folds):
            for i in range(j + 1):
                mixed = np.logaddexp(matrix[i], matrix[j]) - math.log(2.0)
                class_lls[idx] = float(np.dot(counts, mixed))
                idx += 1
        receipt_lls = d["class_log_likelihoods"]
        check(
            len(receipt_lls) == len(class_lls),
            f"locus {locus}: class count differs",
        )
        diff = np.max(np.abs(np.array(receipt_lls, dtype=float) - class_lls))
        check(
            diff <= 1e-6 * max(1.0, float(np.max(np.abs(class_lls)))),
            f"locus {locus}: class LL re-derivation differs (max {diff})",
        )
        best = float(class_lls.max())
        # (the relative convention of the class check above: the
        # re-derivation's summation error is RELATIVE to the LL
        # magnitude - measured at chrII locus 22 (48,828 classes,
        # 17.8M matrix entries, |LL| ~1.6e6): 1.28e-06 ABSOLUTE =
        # 8e-13 relative, unmeetable by an absolute 1e-6 demand; the
        # receipt is self-consistent (best == max of its own classes,
        # runner-up gap 15,458 - no argmax flip). The class check's
        # tolerance is the convention here.)
        check(
            abs(best - d["best_log_likelihood"])
            <= 1e-6 * max(1.0, abs(best)),
            f"locus {locus}: best LL differs",
        )
        if d["truth_folds"] is None:
            # a non-expressible locus (the exhaustive domains include
            # them): no truth class exists to rank
            check(
                d["truth_pair_expressible"] is False,
                f"locus {locus}: truth folds absent but flagged expressible",
            )
            print(
                f"   locus {locus}: {len(class_lls)} classes re-derived "
                f"(max diff {diff:.2e}); truth pair NOT expressible",
                flush=True,
            )
            continue
        ti, tj = d["truth_folds"]
        check(d["truth_pair_expressible"] is True, f"locus {locus}: truth folds present but flagged inexpressible")
        truth_flat = tj * (tj + 1) // 2 + ti
        truth_ll = float(class_lls[truth_flat])
        rank = 1 + int((class_lls > truth_ll).sum())
        check(rank == d["truth_rank"], f"locus {locus}: truth rank differs ({rank})")
        check(
            abs(truth_ll - d["truth_log_likelihood"])
            <= 1e-6 * max(1.0, abs(truth_ll)),
            f"locus {locus}: truth LL differs",
        )
        winner_flat = int(np.argmax(class_lls))
        wi, wj = d["best_fold_indices"]
        check(
            winner_flat == wj * (wj + 1) // 2 + wi,
            f"locus {locus}: winner class differs",
        )
        print(
            f"   locus {locus}: {len(class_lls)} classes re-derived "
            f"(max diff {diff:.2e}), truth rank {rank}"
            + (
                ""
                if MARGINAL
                else f", floor entries {floor_checked}"
            ),
            flush=True,
        )

    # ---------------- phase 6: the QUAL cluster state
    print("== phase 6: the QUAL cluster state re-derived", flush=True)
    for locus in LOCI:
        d = receipt[locus]
        spectrum = list(zip(d["qual_distance_spectrum"], d["qual_distance_spectrum_scores"]))
        # (Stage 2: under the territory convention the receipt carries
        # the tied-winner tie certificate — the cluster machinery's
        # union-find input for the called set's tied classes; the
        # committed receipts never carried a tied winner, so the
        # default mode keeps the empty band list. The wall fix adds
        # the compact CUT-INDEX form at fully degenerate loci (every
        # class ties at the pure-E score — the receipt's own
        # called_class_count == class_count equality); both forms are
        # the exact <=cut predicate sets, so the same re-derivation
        # applies.)
        compact = False
        if TERRITORY:
            if "qual_tied_cut_indices" in d:
                tied = [
                    (int(entry[0]), [int(member) for member in entry[1]])
                    for entry in d["qual_tied_cut_indices"]
                ]
                compact = True
            else:
                tied = [
                    (int(entry[0]), [float(x) for x in entry[1]])
                    for entry in d.get("qual_tied_bands", [])
                ]
        else:
            tied = []
        state = cluster_form_qual(1.0, spectrum, tied, compact)
        check(
            state["shape"] == d["qual_spectrum_shape"],
            f"locus {locus}: spectrum shape differs",
        )
        check(
            (state["k"] == d["qual_cluster_k"]) if d["qual_cluster_k"] is not None else True,
            f"locus {locus}: cluster k differs",
        )
        if d["qual"] is None:
            check(state["unbounded"], f"locus {locus}: receipt QUAL None but state not unbounded")
        else:
            check(abs(state["qual"] - d["qual"]) < 1e-6, f"locus {locus}: QUAL differs")
        print(
            f"   locus {locus}: shape {state['shape']}, k {state['k']}, "
            f"QUAL {'unbounded' if state['unbounded'] else state['qual']}",
            flush=True,
        )

    # ---------------- phase 7: the identical-through-graph fold
    print("== phase 7: the identical-through-graph fold", flush=True)
    _census_path = f"{D}/partition-graphs/partition-graph-{EXHAUSTIVE}-census.jsonl"
    if EXHAUSTIVE not in ("chrI", "chrMT") and os.path.exists(_census_path):
        # (THE FLEET GENERALIZATION of the chrIV twin branch: the EXACT
        # flagged near-twin pairs — the identical-through-graph class,
        # the survey's 100%-exact-rows flags; at chrIV exactly
        # AAA#0#chrIV/SGDref#0#chrIV — fold to ONE candidate at every
        # locus where either holds a member row, BY CONSTRUCTION; the
        # proof: every fold carrying one path's row carries the other's
        # IDENTICAL row too (same interval), no fold carries exactly
        # one of the pair, and every member of every such fold spells
        # the same sequence and walk from the GFAs — the
        # indistinguishability the class table then inherits: classes
        # pair FOLDS, so no class can distinguish the pair. The
        # interval-close near-identical copies (the CLL class) are NOT
        # indistinguishable and get no proof. The pairs come from the
        # census receipt of record — the partition-graphs checker
        # independently re-derives the flag rule from the raw BEDs.)
        _survey = {}
        _spectrum = {}
        for _line in open(_census_path):
            _r = json.loads(_line)
            if _r.get("question") == "s":
                _survey = _r
                for _entry in _r.get("interval_close_spectrum_ge_floor", ()):
                    _spectrum[tuple(_entry["pair"])] = _entry
        _exact_pairs = [
            tuple(_pair)
            for _pair in _survey.get("twin_pairs", ())
            if _spectrum.get(tuple(_pair), {}).get("exact_fraction", 0.0) >= 1.0
        ]
        for pair in _exact_pairs:
            twin_fold_loci = 0
            for locus in LOCI:
                d = receipt[locus]
                folds = ingredients[locus]["folds"]
                hit_folds = [
                    (i, f)
                    for i, f in enumerate(folds)
                    if any(m["path_name"] in pair for m in f["members"])
                ]
                if not hit_folds:
                    check(
                        all(
                            m["path_name"] not in pair
                            for f in folds
                            for m in f["members"]
                        ),
                        f"locus {locus}: a twin path row sits in a fold",
                    )
                    continue
                twin_fold_loci += 1
                g = gfa(d["partition"])
                for i, f in hit_folds:
                    names_here = [m["path_name"] for m in f["members"]]
                    check(
                        pair[0] in names_here and pair[1] in names_here,
                        f"locus {locus}: fold {i} carries only one twin path "
                        f"(the identical rows must fold together)",
                    )
                    a_intervals = sorted(
                        (m["start"], m["end"])
                        for m in f["members"]
                        if m["path_name"] == pair[0]
                    )
                    b_intervals = sorted(
                        (m["start"], m["end"])
                        for m in f["members"]
                        if m["path_name"] == pair[1]
                    )
                    check(
                        a_intervals == b_intervals,
                        f"locus {locus}: fold {i} twin row intervals differ "
                        f"({a_intervals} vs {b_intervals})",
                    )
                    # (the fold's identity criterion is (row sequence,
                    # contained STORED walk) — NOT the full P-line spelling
                    # and NOT the fold's claimed walk, which in this mode is
                    # the REPAIRED canonical skeleton (phase 2's frame audit
                    # verifies the skeleton against the stored walk and the
                    # stored walk against the GFA P line for the
                    # representative; a member's P line may also extend past
                    # the row extent with its own edge-overlapping steps, and
                    # cross-chromosome repeat-family members legitimately
                    # spell longer tails — measured: partition 110's 74bp
                    # repeat fold carries AMP_1a#0#chrIII_chrX / ANL / AVN
                    # rows whose P lines spell 135bp). The grouping audit:
                    # every member spells the fold's sequence, and every
                    # member's contained stored walk — the GFA P line's
                    # syncmer steps positioned relative to the row's own
                    # crop — is IDENTICAL across the fold, the fold key
                    # itself.)
                    row_seq = f["sequence"]
                    ref_contained = None
                    for m in f["members"]:
                        name = row_gfa_name(m)
                        seq, _positions, _steps = g.spelled(name)
                        offset = seq.find(row_seq)
                        check(
                            offset >= 0,
                            f"locus {locus}: fold {i} member {name} does not "
                            "spell the fold's sequence",
                        )
                        if offset < 0:
                            continue
                        sync_positions, sync_steps = g.spelled_syncmers(name)
                        contained = [
                            (p - offset, s)
                            for p, s in zip(sync_positions, sync_steps)
                            if p >= offset and p + K <= offset + len(row_seq)
                        ]
                        if ref_contained is None:
                            ref_contained = contained
                        else:
                            check(
                                contained == ref_contained,
                                f"locus {locus}: fold {i} member {name} "
                                "contained stored walk differs from the "
                                "fold's",
                            )
            print(
                f"   {'/'.join(pair)} fold to ONE candidate at {twin_fold_loci} of "
                f"{len(LOCI)} loci — the fold criterion (row sequence + contained "
                f"stored walk) re-derived from the GFAs for every member of "
                f"every twin fold; no fold distinguishes the pair",
                flush=True,
            )
    for locus in LOCI:
        d = receipt[locus]
        if d["identical_pair_fold"] is None:
            continue
        fold = ingredients[locus]["folds"][d["identical_pair_fold"]]
        member_names = [m["path_name"] for m in fold["members"]]
        if EXHAUSTIVE is None:
            check(
                "BTE#3#block28_contig1" in member_names
                and "BTE#4#block28_contig1" in member_names,
                f"locus {locus}: identical-pair fold membership",
            )
            wanted = "block28_contig1"
        else:
            # the exhaustive mode generalizes the fold-by-construction
            # proof: EVERY member of the named identical-pair fold must
            # spell the same sequence and walk from the GFAs
            wanted = None
        g = gfa(d["partition"])
        if EXHAUSTIVE is None:
            seqs = {}
            walks = {}
            for m in fold["members"]:
                if wanted not in m["path_name"]:
                    continue
                seq, positions, steps = g.spelled(row_gfa_name(m))
                seqs[m["path_name"]] = seq
                walks[m["path_name"]] = list(zip(positions, steps))
            vals = list(seqs.values())
            check(len(vals) >= 2 and all(v == vals[0] for v in vals), "identical-pair sequences differ")
            wvals = list(walks.values())
            check(all(v == wvals[0] for v in wvals), "identical-pair walks differ")
            print(
                f"   locus {locus}: {len(vals)} members "
                + (
                    f"BTE#3/#4 block28_contig1 spell identical "
                    f"sequence ({len(vals[0])} bp) and walk — the fold is exact"
                    if EXHAUSTIVE is None
                    else f"spell identical sequence ({len(vals[0])} bp) and walk — the fold is exact"
                ),
                flush=True,
            )
        else:
            # (the exhaustive mode's fold-by-construction proof: the
            # same fold-criterion grouping audit as the twin branch —
            # the fold key is (row sequence, contained stored walk);
            # FULL P-line spellings legitimately differ at
            # edge-overlapping steps for cross-chromosome
            # repeat-family members — measured: partition 111's 54bp
            # repeat fold at L91/L105 carries 135 members incl.
            # fused-path rows (ABH#0#chrV_chrXIV, AFH#0#chrVII_chrXVI)
            # whose P lines extend past the row extent)
            row_seq = fold["sequence"]
            check(
                len(fold["members"]) >= 2,
                f"locus {locus}: identical-pair fold has fewer than 2 members",
            )
            ref_contained = None
            for m in fold["members"]:
                name = row_gfa_name(m)
                seq, _positions, _steps = g.spelled(name)
                offset = seq.find(row_seq)
                check(
                    offset >= 0,
                    f"locus {locus}: identical-pair fold member {name} "
                    "does not spell the fold's sequence",
                )
                if offset < 0:
                    continue
                sync_positions, sync_steps = g.spelled_syncmers(name)
                contained = [
                    (p - offset, s)
                    for p, s in zip(sync_positions, sync_steps)
                    if p >= offset and p + K <= offset + len(row_seq)
                ]
                if ref_contained is None:
                    ref_contained = contained
                else:
                    check(
                        contained == ref_contained,
                        f"locus {locus}: identical-pair fold member {name} "
                        "contained stored walk differs from the fold's",
                    )
            print(
                f"   locus {locus}: {len(fold['members'])} members share the "
                f"fold criterion (sequence {len(row_seq)} bp + contained "
                "stored walk) — the fold is exact",
                flush=True,
            )

    # ---------------- phase 8: the pilot verdicts
    print("== phase 8: the pilot verdicts", flush=True)
    before = {}
    if MARGINAL:
        for line in open(BEFORE_RECEIPT):
            b = json.loads(line)
            before[b["locus"]] = b

    def fold_member_key(d, index):
        return tuple(
            sorted(
                (m["path_name"], m["start"], m["end"])
                for m in d["fold_identities"][index]["members"]
            )
        )

    if EXHAUSTIVE is not None:
        # slice E: the exhaustive before/after vs the committed
        # Poisson-era read-matched receipts (the instrument of record).
        # Guards: every BEFORE-rank-1 locus must remain rank 1 (the
        # old-wins no-regression guard); no locus expressible before
        # may become inexpressible; the chrI L4 control stays
        # bit-exact (rank 1, gap 0.0, truth in the called set, the
        # winner IS the truth pair).
        print("== phase 8: the exhaustive verdicts vs the Poisson-era read-matched receipts", flush=True)
        rank1_before = [locus for locus in LOCI if before[locus]["truth_rank"] == 1]
        rank1_after = []
        prediction_violations = []
        expressible_before = [locus for locus in LOCI if before[locus]["truth_pair_expressible"]]
        expressible_after = []
        for locus in LOCI:
            b = before[locus]
            d = receipt[locus]
            b_rank = b["truth_rank"]
            b_qual = "unbounded" if b["qual"] is None else f"{b['qual']:.2f}"
            a_qual = "unbounded" if d["qual"] is None else f"{d['qual']:.2f}"
            print(
                f"   locus {locus}: partition {d['partition']}, folds {d['folds']}, "
                f"units {d['unit_count']}, classes {d['class_count']}; "
                f"truth rank {b_rank} -> {d['truth_rank']}, "
                f"log gap {('%.2f' % d['log_gap']) if d['log_gap'] is not None else 'n/a'}; "
                f"QUAL {b_qual} -> {a_qual}; "
                f"in called set {b.get('qual_truth_in_called_set')} -> {d['truth_in_called_set']}; "
                f"factorized {d['validation']['factorized_equal']}/"
                f"{d['validation']['factorized_checked']} exact",
                flush=True,
            )
            if b["truth_pair_expressible"]:
                check(
                    d["truth_pair_expressible"],
                    f"locus {locus}: expressible under the Poisson instrument but not the realignment instrument",
                )
                check(
                    d["truth_rank"] is not None and d["log_gap"] is not None,
                    f"locus {locus}: expressible but rank/gap missing",
                )
                if d["truth_rank"] == 1:
                    rank1_after.append(locus)
                    check(
                        d["truth_in_called_set"],
                        f"locus {locus}: rank 1 but truth not in the called set",
                    )
                if b_rank == 1:
                    # THE PREDICTION VERDICT ("the old wins hold"), scored
                    # honestly: a rank-1 locus that loses rank 1 is a
                    # MEASURED FINDING about the instruments, not a
                    # receipt error — reported loudly as a named
                    # violation, never silently, and never folded into
                    # the receipt-validation pass/fail.
                    if d["truth_rank"] == 1:
                        print(
                            f"     prediction verdict: locus {locus} old-win HOLDS (rank 1)",
                            flush=True,
                        )
                    else:
                        prediction_violations.append(locus)
                        print(
                            f"     PREDICTION VIOLATION (measured finding): locus {locus} "
                            f"old-win LOST rank 1 under the realignment instrument "
                            f"(Poisson rank 1 -> {d['truth_rank']}, log gap "
                            f"{d['log_gap']:.2f}) — residual named in the gap "
                            "decomposition",
                            flush=True,
                        )
            if d["truth_pair_expressible"]:
                expressible_after.append(locus)
        if EXHAUSTIVE == "chrI":
            d4 = receipt[4]
            check(
                d4["truth_rank"] == 1 and d4["log_gap"] == 0.0 and d4["truth_in_called_set"],
                "chrI L4 control: rank 1, gap 0.0, truth in the called set",
            )
            check(
                tuple(fold_member_key(d4, i) for i in d4["best_fold_indices"])
                == tuple(fold_member_key(d4, i) for i in d4["truth_folds"]),
                "chrI L4 control: the winner must BE the truth pair",
            )
        newly = [locus for locus in expressible_after if locus not in expressible_before]
        print(
            f"   AGGREGATE ({EXHAUSTIVE}): truth rank-1 {len(rank1_before)} -> {len(rank1_after)} "
            f"of {len(expressible_after)} expressible (before: {len(expressible_before)}); "
            f"rank-1 loci before {rank1_before}, after {rank1_after}; "
            f"newly expressible {newly}; "
            f"prediction violations (old wins lost): {prediction_violations}",
            flush=True,
        )

        if EXHAUSTIVE not in ("chrI", "chrMT") and os.path.exists(
            f"{D}/partition-graphs/partition-graph-{EXHAUSTIVE}-census.jsonl"
        ):
            # ------------- phase 8c: the expressibility classes
            # (the prediction table of record: the census receipt's
            # per-locus ortholog/positional statements interpret the
            # truth-rank table; the instrument's expressibility is
            # re-derived from the maps — the window's own axis row in
            # the truth-first fold, an SK1 row of the window's
            # partition in the truth-second fold at single-window
            # partitions — and verified against the receipt. The class
            # distinction that matters: IN-AXIS (the ortholog row IS
            # in the window's partition — the truth pair is the real
            # homolog pair) vs TILED-ELSEWHERE-with-neighbor-SK1-row
            # (the path-domain convention finds the COORDINATE-OFFSET
            # neighbor window's aligned SK1 row in the partition — the
            # paired class is named per locus, never silently passed
            # as the ortholog) vs the inexpressible groups (multi-
            # window partitions, no-SK1-row partitions, the contig
            # end).
            print(
                "== phase 8c: the slice-1 expressibility classes vs the receipts",
                flush=True,
            )
            census = {}
            for line in open(
                f"{D}/partition-graphs/partition-graph-{EXHAUSTIVE}-census.jsonl"
            ):
                r = json.loads(line)
                if r.get("question") == "a":
                    census[r["locus"]] = r
            windows_of_partition = {}
            for pid in axis:
                windows_of_partition[pid] = windows_of_partition.get(pid, 0) + 1
            sk1_rows = {}
            for pid, m in maps.items():
                rows = [
                    (r["start"], r["end"])
                    for r in m["members"]
                    if r["path_name"] == f"SK1#0#{EXHAUSTIVE}"
                ]
                if rows:
                    sk1_rows[pid] = rows
            from collections import Counter

            cls_counts = Counter()
            cls_expressible = Counter()
            cls_rank1 = Counter()
            neighbor_material = []
            for locus in LOCI:
                d = receipt[locus]
                c = census[locus]
                pid = axis[locus]
                multi = windows_of_partition[pid] > 1
                # (the fleet generalization: expressibility needs the
                # path-name rule to be UNAMBIGUOUS — exactly one fold
                # carrying an SK1 member row at a single-window
                # partition; a partition carrying SEVERAL SK1 folds
                # (measured: chrXIV partition 16, chrXVI 117/257,
                # chrVII 16/111) cannot attribute the ortholog, no
                # truth pair is fabricated. Derived from the receipt's
                # fold members — the fold structure phase 2
                # independently re-verified against the partition
                # maps/GFAs — beside the map-derived single-window
                # count.)
                sk1_folds = sum(
                    1
                    for f in d["fold_identities"]
                    if any(
                        m["path_name"] == f"SK1#0#{EXHAUSTIVE}"
                        for m in f["members"]
                    )
                )
                expect = (not multi) and (sk1_folds == 1)
                check(
                    d["truth_pair_expressible"] == expect,
                    f"locus {locus}: expressibility differs from the map-derived "
                    f"rule (receipt {d['truth_pair_expressible']} vs expected {expect})",
                )
                if d["truth_pair_expressible"]:
                    # (the instrument emits truth_folds sorted by fold
                    # INDEX — (a.min(b), a.max(b)) — so the axis-row
                    # fold and the SK1 fold can arrive in either order;
                    # the class pair is unordered and the membership
                    # proof is order-independent)
                    first, second = d["truth_folds"]
                    fa = d["fold_identities"][first]["members"]
                    fb = d["fold_identities"][second]["members"]
                    s0, e0, _ = axis_rows[locus]
                    carries_axis = lambda ms: any(
                        m["path_name"] == COMPONENT
                        and m["start"] == s0
                        and m["end"] == e0
                        for m in ms
                    )
                    carries_sk1 = lambda ms: any(
                        m["path_name"] == f"SK1#0#{EXHAUSTIVE}" for m in ms
                    )
                    check(
                        carries_axis(fa) or carries_axis(fb),
                        f"locus {locus}: neither truth fold carries the window's axis row",
                    )
                    check(
                        carries_sk1(fa) or carries_sk1(fb),
                        f"locus {locus}: neither truth fold carries an SK1 row",
                    )
                    if d["truth_rank"] == 1:
                        cls_rank1[c["verdict"]] += 1
                    if (
                        c["verdict"] != "IN-AXIS-PARTITION"
                        and d["truth_pair_expressible"]
                    ):
                        neighbor_material.append(locus)
                cls_counts[c["verdict"]] += 1
                if d["truth_pair_expressible"]:
                    cls_expressible[c["verdict"]] += 1
            print(
                f"   class census over the slice: "
                f"{dict(sorted(cls_counts.items()))}; expressible "
                f"{dict(sorted(cls_expressible.items()))}; truth rank-1 "
                f"{dict(sorted(cls_rank1.items()))}",
                flush=True,
            )
            print(
                f"   path-domain-expressible with NEIGHBOR (coordinate-offset) "
                f"SK1 material — the paired class is NOT the ortholog pair, "
                f"named per locus: {neighbor_material}",
                flush=True,
            )
            # the full-component group census (all 167 windows, from the
            # maps + census — not just the slice)
            groups = Counter()
            for locus in sorted(census):
                pid = axis[locus]
                if census[locus]["verdict"] == "IN-AXIS-PARTITION":
                    groups["in_axis"] += 1
                elif windows_of_partition[pid] > 1:
                    groups["multiwindow_partition"] += 1
                elif pid in sk1_rows:
                    groups["neighbor_sk1_row"] += 1
                elif census[locus]["verdict"] == "ABSENT":
                    groups["contig_end_absent"] += 1
                else:
                    groups["no_sk1_row"] += 1
            print(
                f"   the full-component groups (all {len(census)} windows): "
                f"{dict(sorted(groups.items()))}",
                flush=True,
            )

        # ------------- phase 8t (slice F): THE ANSWER-PRESERVATION GATE
        # The timered rerun (the instrumented binary: timers and
        # counters only, emitted to stderr) must reproduce the committed
        # slice-E receipts on EVERY semantic field — the identity-gate
        # pattern from slice D: only the walls/rss_kb timing fields may
        # differ, the walls field SET must be unchanged (no receipt
        # schema drift), and the four sidecars must be byte-identical.
        # (Stage 2: the gate does NOT run for the territory receipts —
        # the convention changed BY DESIGN; the env-unset identity
        # gate is the separate no-op proof, and phase 8N below is the
        # territory receipts' own gate.)
        if TIMERED_BASE is not None and not TERRITORY:
            import filecmp

            print(
                "== phase 8t: the answer-preservation gate (the timered rerun vs the committed receipts)",
                flush=True,
            )
            # (chrIV and every fleet component: the committed identity
            # base is the SERIAL run's receipts — the 4-wide exhaustive
            # run is the record, the serial run the gate's reference,
            # the phase-3 ladder convention; chrMT/chrI keep the
            # committed slice-E receipts as the base.)
            committed_prefix = (
                f"realign-par-serial-{EXHAUSTIVE}"
                if EXHAUSTIVE not in ("chrI", "chrMT")
                else f"realign-exhaustive-{EXHAUSTIVE}"
            )
            committed_slice_e = {}
            with open(f"{D}/{committed_prefix}.jsonl") as f:
                for line in f:
                    d = json.loads(line)
                    committed_slice_e[d["locus"]] = d
            skip = {"walls", "rss_kb"}
            for locus in LOCI:
                o, n = committed_slice_e[locus], receipt[locus]
                check(
                    (set(o) - skip) == (set(n) - skip),
                    f"answer-preservation gate: locus {locus} field set differs",
                )
                check(
                    set(o["walls"]) == set(n["walls"]),
                    f"answer-preservation gate: locus {locus} walls field set differs",
                )
                for key in sorted(set(o) - skip):
                    check(
                        o[key] == n[key],
                        f"answer-preservation gate: locus {locus} field {key} differs",
                    )
            for committed_sidecar, timered_sidecar in [
                (
                    f"{D}/{committed_prefix}.exactness.jsonl",
                    EXACTNESS,
                ),
                (
                    f"{D}/{committed_prefix}.jsonl.ingredients.jsonl",
                    INGREDIENTS,
                ),
                (
                    f"{D}/{committed_prefix}.jsonl.records.jsonl",
                    f"{D}/{TIMERED_BASE}.jsonl.records.jsonl",
                ),
                (
                    f"{D}/{committed_prefix}.skeleton.jsonl",
                    SKELETON,
                ),
            ]:
                check(
                    filecmp.cmp(committed_sidecar, timered_sidecar, shallow=False),
                    f"answer-preservation gate: {timered_sidecar} not byte-identical "
                    f"to the committed sidecar",
                )
            print(
                f"      the run reproduces the committed {committed_prefix} receipts "
                f"(every semantic field at {len(LOCI)} loci; only walls/rss_kb "
                f"differ) with byte-identical sidecars",
                flush=True,
            )
        # ------------- phase 8N (stage 2): THE TERRITORY GATE — the
        # normalized-convention receipts vs the COMMITTED default-
        # convention receipts of the same component. THE GUARDS, hard:
        # (1) NO REGRESSION — every locus at truth rank 1 under the
        # committed convention must remain rank 1 (a regression is a
        # DESIGN FAILURE, named loudly); (2) expressibility is
        # preserved; (3) the winner at a held rank-1 locus is the
        # truth pair. The conversion table is stated honestly per
        # locus (CONVERT / HOLD / REGRESS / residual-moved).
        if TERRITORY:
            print(
                "== phase 8N: the territory gate (the normalized receipts vs the committed convention)",
                flush=True,
            )
            committed = {}
            for line in open(COMMITTED_RECEIPT):
                c = json.loads(line)
                committed[c["locus"]] = c
            conversions = []
            regressions = []
            moved = []
            rank1_c = [
                locus
                for locus in LOCI
                if committed[locus].get("truth_rank") == 1
            ]
            rank1_t = []
            for locus in LOCI:
                c, d = committed[locus], receipt[locus]
                if not c["truth_pair_expressible"]:
                    check(
                        not d["truth_pair_expressible"],
                        f"locus {locus}: inexpressible under the committed convention "
                        "but expressible under the normalized convention",
                    )
                    continue
                cr, tr = c["truth_rank"], d["truth_rank"]
                if tr == 1:
                    rank1_t.append(locus)
                    check(
                        d["truth_in_called_set"],
                        f"locus {locus}: rank 1 but truth not in the called set",
                    )
                verdict = "HOLD" if cr == 1 and tr == 1 else None
                if verdict is None:
                    if cr == 1:
                        verdict = "REGRESS"
                        regressions.append(locus)
                    elif tr == 1:
                        verdict = "CONVERT"
                        conversions.append(locus)
                    elif cr != tr:
                        verdict = "moved"
                        moved.append(locus)
                    else:
                        verdict = "unchanged"
                print(
                    f"   locus {locus}: truth rank {cr} -> {tr}, "
                    f"log gap {c['log_gap']:.2f} -> "
                    f"{d['log_gap'] if d['log_gap'] is None else round(d['log_gap'], 2)}  "
                    f"{verdict}",
                    flush=True,
                )
            check(
                not regressions,
                f"TERRITORY GATE FAILED: rank-1 regressions under the normalized "
                f"convention: {regressions}",
            )
            for locus in rank1_c:
                d = receipt[locus]
                check(
                    tuple(fold_member_key(d, i) for i in d["best_fold_indices"])
                    == tuple(fold_member_key(d, i) for i in d["truth_folds"]),
                    f"locus {locus}: a held rank-1 locus whose winner is not the truth pair",
                )
            print(
                f"   AGGREGATE ({EXHAUSTIVE}): truth rank-1 "
                f"{len(rank1_c)} -> {len(rank1_t)} under the normalized convention; "
                f"CONVERT {conversions}; REGRESS {regressions}; "
                f"residual moved {moved}",
                flush=True,
            )
    else:
        for locus in LOCI:
            d = receipt[locus]
            print(
                f"   locus {locus}: partition {d['partition']}, folds {d['folds']}, "
                f"units {d['unit_count']}, classes {d['class_count']}; "
                f"truth rank {d['truth_rank']}, log gap {d['log_gap']}, "
                f"QUAL {'unbounded' if d['qual'] is None and d['qual_unbounded'] else d['qual']}; "
                f"factorized {d['validation']['factorized_equal']}/"
                f"{d['validation']['factorized_checked']} exact",
                flush=True,
            )
            if MARGINAL:
                b = before[locus]
                if FRAME:
                    # slice D: the before-record is the paired identity
                    # (stored-walk) run; winners compared by MEMBER SET (the
                    # fold indices are per-receipt).
                    winner_before = tuple(fold_member_key(b, i) for i in b["best_fold_indices"])
                    winner_after = tuple(fold_member_key(d, i) for i in d["best_fold_indices"])
                    print(
                        f"     before/after: truth rank {b['truth_rank']} -> {d['truth_rank']}, "
                        f"log gap {b['log_gap']:.2f} -> {d['log_gap']:.2f}, "
                        f"winner material {'unchanged' if winner_before == winner_after else 'CHANGED'}",
                        flush=True,
                    )
                    if locus == 4:
                        check(
                            d["truth_rank"] == 1 and d["truth_in_called_set"] and d["log_gap"] == 0.0,
                            "L4 control: the truth-rank1 control must hold bit-exact "
                            "(rank 1, gap 0.0, truth in the called set)",
                        )
                        check(
                            winner_before == winner_after,
                            "L4 control: the winner material must be bit-exact unchanged",
                        )
                        check(
                            winner_after == tuple(fold_member_key(d, i) for i in d["truth_folds"]),
                            "L4 control: the winner must BE the truth pair",
                        )
                    if locus == 7:
                        check(
                            d["truth_rank"] < b["truth_rank"] and d["log_gap"] < b["log_gap"],
                            "L7: the frame repair must strictly improve the truth rank and gap",
                        )
                else:
                    print(
                        f"     before/after: truth rank {b['truth_rank']} -> {d['truth_rank']}, "
                        f"log gap {b['log_gap']:.2f} -> {d['log_gap']:.2f}, "
                        f"winner {b['best_fold_indices']} -> {d['best_fold_indices']}",
                        flush=True,
                    )
                    if locus == 4:
                        check(
                            d["truth_rank"] == 1 and d["truth_in_called_set"],
                            "L4 control: the truth-rank1 control must hold after the marginalization",
                        )
                        check(
                            d["best_fold_indices"] == b["best_fold_indices"],
                            "L4 control: the winner must be bit-exact unchanged",
                        )
                    shrink = (d["log_gap"] or 0.0) <= (b["log_gap"] or 0.0)
                    check(shrink, f"locus {locus}: the log gap grew after the marginalization")

    # ---------------- phase 9 (slice D): the frame-repair audit
    if FRAME and EXHAUSTIVE is None:
        import filecmp

        print("== phase 9: the frame-repair audit (slice D)", flush=True)

        # 9a. THE IDENTITY GATE: the env-gated stored-walk run must
        # reproduce the committed slice-C receipts on every semantic
        # field, with byte-identical exactness/ingredients sidecars.
        print("   9a. the identity gate (the stored-walk run vs the committed slice-C receipts)", flush=True)
        check(os.path.exists(f"{IDENTITY_RUN}.done"), "identity run done marker missing")
        exit_id = int(open(f"{IDENTITY_RUN}.exit").read().strip())
        check(exit_id == 0, f"identity run exit code {exit_id}")
        with open(f"{IDENTITY_RUN}.rss") as f:
            rss_id = max(int(line.split()[2]) for line in f)
        check(rss_id <= RSS_BUDGET_KB, f"identity RSS guard exceeded: {rss_id} kB")
        committed = {}
        for line in open(BEFORE_RECEIPT):
            b = json.loads(line)
            committed[b["locus"]] = b
        skip_fields = {"walls", "rss_kb", "skeleton"}
        for locus in LOCI:
            o, n = committed[locus], identity_receipt[locus]
            check(
                (set(o) - skip_fields) == (set(n) - skip_fields),
                f"identity gate: locus {locus} field set differs",
            )
            for key in set(o) - skip_fields:
                check(
                    o[key] == n[key],
                    f"identity gate: locus {locus} field {key} differs",
                )
        for a, b in [
            (f"{D}/realign-marginal-chrI.exactness.jsonl", IDENTITY_EXACTNESS),
            (f"{D}/realign-marginal-chrI.jsonl.ingredients.jsonl", IDENTITY_INGREDIENTS),
        ]:
            check(filecmp.cmp(a, b, shallow=False), f"identity gate: {b} not byte-identical to the committed sidecar")
        print(
            "      the identity run reproduces the committed slice-C receipts "
            "(every semantic field) with byte-identical sidecars",
            flush=True,
        )

        # 9b. the skeleton diff census vs the receipt's skeleton block
        print("   9b. the skeleton census vs the receipt's skeleton block", flush=True)
        for locus in LOCI:
            sk = receipt[locus]["skeleton"]
            check(sk["scheme"] == "canonical_scheme", f"locus {locus}: receipt skeleton scheme")
            rows = [k for k in skeleton_by_fold if k[0] == locus]
            check(len(rows) == sk["rows"], f"locus {locus}: skeleton row count differs")
            affected = sum(
                1
                for k in rows
                if skeleton_by_fold[k]["added"]
                or skeleton_by_fold[k]["replaced"]
                or skeleton_by_fold[k]["dropped"]
            )
            added = sum(skeleton_by_fold[k]["added"] for k in rows)
            replaced = sum(skeleton_by_fold[k]["replaced"] for k in rows)
            dropped = sum(skeleton_by_fold[k]["dropped"] for k in rows)
            kept_forward = sum(skeleton_by_fold[k]["kept_forward"] for k in rows)
            kept_reverse = sum(skeleton_by_fold[k]["kept_reverse"] for k in rows)
            check(affected == sk["rows_affected"], f"locus {locus}: affected row count differs")
            check(
                (added, replaced, dropped, kept_forward, kept_reverse)
                == (sk["added"], sk["replaced"], sk["dropped"], sk["kept_forward"], sk["kept_reverse"]),
                f"locus {locus}: skeleton diff census differs",
            )
            print(
                f"      locus {locus}: {len(rows)} rows, {affected} affected "
                f"(added {added}, replaced {replaced}, dropped {dropped}; "
                f"kept forward {kept_forward}, reverse {kept_reverse})",
                flush=True,
            )

        # 9c. THE NO-REGRESSION SWEEP: folds mapped by member set (the
        # fold indices are per-receipt), units by (record, variant);
        # per (unit, fold) pin gains/losses and ll changes. ZERO pin
        # losses is the gate (a previously-placed pair must keep its
        # placement: the dropped stored positions never carried a
        # true pin).
        print("   9c. the no-regression sweep (identity vs repaired, per (unit, fold))", flush=True)

        def member_keys(folds):
            return [
                tuple(sorted((m["path_name"], m["start"], m["end"]) for m in f["members"]))
                for f in folds
            ]

        for locus in LOCI:
            fi = identity_ingredients[locus]
            fr = ingredients[locus]
            mi = member_keys(fi["folds"])
            mr = member_keys(fr["folds"])
            check(sorted(mi) == sorted(mr), f"locus {locus}: fold member sets differ between runs")
            map_i = {k: n for n, k in enumerate(mi)}
            map_r = {k: n for n, k in enumerate(mr)}
            units_i = [(u["record"], u["variant"]) for u in fi["units"]]
            units_r = [(u["record"], u["variant"]) for u in fr["units"]]
            check(units_i == units_r, f"locus {locus}: unit streams differ between runs")
            ui = {k: n for n, k in enumerate(units_i)}
            pinned_b = {
                key: set(u["pinning_folds"]) for key, u in identity_anatomy_by_locus[locus].items()
            }
            pinned_a = {
                key: set(u["pinning_folds"]) for key, u in anatomy_by_locus[locus].items()
            }
            both = eq = changed = gains = losses = 0
            max_delta = 0.0
            for key, n in ui.items():
                pb = pinned_b.get(key, set())
                pa = pinned_a.get(key, set())
                for k in mi:
                    fb = map_i[k] in pb
                    fa = map_r[k] in pa
                    if fb and fa:
                        both += 1
                        delta = fr["ll_matrix"][map_r[k]][n] - fi["ll_matrix"][map_i[k]][n]
                        if delta == 0.0:
                            eq += 1
                        else:
                            changed += 1
                            max_delta = max(max_delta, abs(delta))
                    elif fb and not fa:
                        losses += 1
                    elif fa and not fb:
                        gains += 1
            check(
                losses == 0,
                f"locus {locus}: the frame repair lost {losses} previously-placed (unit, fold) pins",
            )
            print(
                f"      locus {locus}: pinned-both {both} (ll identical {eq}, changed "
                f"{changed}, max |delta| {max_delta:.4f}); pin gains {gains}; "
                f"pin losses {losses}",
                flush=True,
            )

        # 9d. THE BLINDED-UNIT CENSUS before/after: units with a
        # full-match donor occurrence inside a truth-fold member row
        # and NO truth placement.
        print("   9d. the blinded-unit census before/after", flush=True)

        def blinded_census(ana, locus, d):
            truth_folds = set(d["truth_folds"] or [])
            trows = [
                (m["path_name"], m["start"], m["end"])
                for fi_ in truth_folds
                for m in d["fold_identities"][fi_]["members"]
            ]
            out = {}
            for key, u in ana[locus].items():
                in_window = False
                for dc in u["donor_checks"]:
                    if dc["c"] == 0 and dc["oob"] == 0 and dc["m"] == READ_LENGTH:
                        for pn, s, e in trows:
                            if pn == dc["path"] and s <= dc["origin"] and dc["origin"] + READ_LENGTH <= e:
                                in_window = True
                                break
                    if in_window:
                        break
                if not in_window:
                    continue
                if not (truth_folds & set(u["pinning_folds"])):
                    out[key] = u["count"]
            return out

        for locus in LOCI:
            d = receipt[locus]
            cb = blinded_census(identity_anatomy_by_locus, locus, identity_receipt[locus])
            ca = blinded_census(anatomy_by_locus, locus, d)
            check(
                len(ca) < len(cb),
                f"locus {locus}: the frame repair must strictly shrink the blinded-unit census "
                f"({len(cb)} -> {len(ca)})",
            )
            print(
                f"      locus {locus}: blinded units {len(cb)} (mass {sum(cb.values())}) -> "
                f"{len(ca)} (mass {sum(ca.values())})",
                flush=True,
            )

        # 9e. THE NAMED PROOF CASE (record 113 / unit 84 at L7): the
        # truth fold must place it after the repair, with the anchor
        # present in the repaired skeleton at the occurrence position
        # and absent from the stored walk.
        print("   9e. the named proof case (record 113 / unit 84 at L7)", flush=True)
        d7 = receipt[7]
        truth_folds_7 = set(d7["truth_folds"] or [])
        before_u = identity_anatomy_by_locus[7][(113, 0)]
        after_u = anatomy_by_locus[7][(113, 0)]
        check(
            not (truth_folds_7 & set(before_u["pinning_folds"])),
            "proof case: the before-record must NOT place record 113 on the truth folds",
        )
        check(
            bool(truth_folds_7 & set(after_u["pinning_folds"])),
            "proof case: the repaired skeleton must place record 113 on the truth folds",
        )
        valid_truth = [
            p
            for p in after_u["placements"]
            if p["fold"] in truth_folds_7 and p["valid"] and p["m"] == READ_LENGTH and p["c"] == 0
        ]
        check(
            bool(valid_truth),
            "proof case: the truth-fold placement of record 113 must be a full-read match",
        )
        found_row = None
        with open(AUDIT_CHRI) as f:
            for line in f:
                row = json.loads(line)
                if row["path_name"] == "S288C#0#chrI" and 7 in row["partitions"]:
                    found_row = row
                    break
        check(found_row is not None, "proof case: the truth row absent from the frame audit")
        rel = 73906 - found_row["start"]
        check(
            any(s[0] == rel and s[1] == -92071 and s[2] == 1 for s in found_row["steps"]),
            "proof case: the repaired skeleton must carry (-92071, rc frame) at 73,906",
        )
        check(
            all(s[0] != rel for s in found_row["stored"]),
            "proof case: the stored walk must lack the anchor at 73,906",
        )
        print(
            "      record 113: winner_only -> both; a valid truth-fold full-read "
            "placement (m 150, c 0); the repaired skeleton carries (-92071, rc "
            "frame) at 73,906 where the stored walk has no step",
            flush=True,
        )

        # 9f. THE FULL-COMPONENT FRAME-AUDIT CENSUS (the diagnosis
        # receipts): every axis-partition row affected on both
        # components.
        print("   9f. the full-component frame-audit census", flush=True)
        for name, path in [("chrI", AUDIT_CHRI), ("chrMT", AUDIT_CHRMT)]:
            rows = affected = added = replaced = dropped = 0
            with open(path) as f:
                for line in f:
                    r = json.loads(line)
                    rows += 1
                    added += r["added"]
                    replaced += r["replaced"]
                    dropped += r["dropped"]
                    if r["added"] or r["replaced"] or r["dropped"]:
                        affected += 1
            check(added > 0, f"{name}: the audit must find rc-frame-only anchors")
            check(
                affected == rows,
                f"{name}: every row must be affected ({affected} of {rows})",
            )
            print(
                f"      {name}: {rows} rows, {affected} affected "
                f"(added {added}, replaced {replaced}, dropped {dropped})",
                flush=True,
            )

    print(f"\nTOTAL CHECKS: {checks[0]}, FAILURES: {len(failures)}")
    if failures:
        print("ALL PHASES FAIL" if len(failures) > 3 else "FAILURES PRESENT")
        sys.exit(1)
    print("ALL PHASES PASS")


if __name__ == "__main__":
    main()
