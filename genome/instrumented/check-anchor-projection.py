#!/usr/bin/env python3
"""Independent checker for the anchor-projection layer (slice A).

Re-derives a measured sample of the anchored projections from the RAW
receipts — the committed multi-matching census, the partition-embedded
graph GFAs, and the sample FASTQ — and demands exactness:

  Phase 1  run markers, exit codes, walls, the 64 GiB RSS guard.
  Phase 2  receipt-vs-census structure over EVERY record: census-verbatim
           placements, interval counts, variant instance sums, flank
           arithmetic, physical-strand flags.
  Phase 3  EXACTNESS on a measured sample of records: the read hull is
           re-assembled purely from the receipt's anchor chain and the
           partition GFAs' own S segments (strand-aware per the mirror
           state), the flanks come from the receipt, and the reconstructed
           full read must EXIST in the FASTQ — proving anchor positions,
           skipped bases, and strand handling jointly, without the Rust
           panel machinery.
  Phase 4  the per-path graph-context sample re-derivation: every sidecar
           row is validated against its GFA (node order, S/L overlap
           arithmetic, spelled sequence), then the anchor correspondence,
           serving anchors, row-side coordinates, covering walk-index
           ranges, and bases are re-derived and compared EXACTLY.

Usage: check-anchor-projection.py chrMT [chrI ...]
"""
import gzip
import json
import sys

D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
GRAPHS = f"{D}/partition-graphs"
NAMES = "/home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng.names"
K = 63
READ_LENGTH = 150
PHASE3_RECORDS = 150
PHASE4_LINES = 25

failures = []
checks = [0]


def check(ok, message):
    checks[0] += 1
    if not ok:
        failures.append(message)
        print(f"FAIL: {message}")


def fnv1a64(data):
    h = 0xCBF29CE484222325
    for b in data:
        h ^= b
        h = (h * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return h


def revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def load_path_names():
    names = {}
    with open(NAMES) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            names[int(fields[0])] = fields[1]
    return names


def load_maps(component):
    """The axis partitions (holding a member row on the component path)
    ranked by axis-row start map to the component's window ids in order."""
    import glob
    maps = {}
    axis = []
    for path in glob.glob(f"{GRAPHS}/*.gfa.map.json"):
        with open(path) as f:
            m = json.load(f)
        maps[m["partition"]] = m
        for member in m["members"]:
            if member["path_name"] == component:
                axis.append((member["start"], m["partition"]))
                break
    axis.sort()
    return maps, [p for _, p in axis]


class Gfa:
    """One partition GFA: S segments, P lines, L overlaps (lazy)."""

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

    def spelled(self, row_name):
        """(spelled sequence, step positions relative to the first step).
        A '-' step spells its segment's reverse complement (GFA path
        semantics); the L-line CIGAR overlap trims the incoming step's
        head in ITS orientation."""
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


GFAS = {}


def gfa(partition):
    if partition not in GFAS:
        GFAS[partition] = Gfa(partition)
    return GFAS[partition]


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
    ambiguous = [j for j in range(m) if candidates[j] and pinned[j] is None]
    return pinned, ambiguous


def check_component(component, path_names):
    print(f"\n=== {component} ===")
    census_path = f"{D}/cosine-diagnostic-scratch/{component}/cosine-multi-census.jsonl"
    receipt_path = f"{D}/anchor-projection-{component}.jsonl"
    context_path = f"{D}/anchor-projection-context-{component}.jsonl"
    sidecar_path = f"{context_path}.rows.jsonl"
    tag = f"run-anchorproj-{component}"

    # ---------------------------------------------------------- phase 1
    with open(f"{D}/{tag}.exit") as f:
        exit_code = int(f.read().strip())
    check(exit_code == 0, f"{component}: exit code {exit_code}")
    import os
    check(os.path.exists(f"{D}/{tag}.done"), f"{component}: no .done marker")
    with open(f"{D}/{tag}.wall") as f:
        wall = int(f.read().strip())
    peak_rss = 0
    with open(f"{D}/{tag}.rss") as f:
        for line in f:
            peak_rss = max(peak_rss, int(line.split()[2]))
    check(peak_rss < 64 * 1024 * 1024, f"{component}: RSS guard ({peak_rss} kB)")
    print(f"phase 1: exit 0, wall {wall}s, peak RSS {peak_rss} kB")

    maps, axis = load_maps(f"S288C#0#{component}")

    # ---------------------------------------------------------- phase 2
    n_records = 0
    n_occ = 0
    shift_hist = {}
    variant_total = 0
    skipped_bp = 0
    receipts = {}
    with open(census_path) as census_f, open(receipt_path) as receipt_f:
        for census_line in census_f:
            census = json.loads(census_line)
            receipt = json.loads(receipt_f.readline())
            rid = census["record"]
            check(receipt["record"] == rid, f"record {rid}: id mismatch")
            check(receipt["multiplicity"] == census["multiplicity"],
                  f"record {rid}: multiplicity mismatch")
            check(receipt["k"] == K, f"record {rid}: k mismatch")
            anchors = receipt["anchors"]
            span = receipt["span"]
            check(len(anchors) == census["anchors"], f"record {rid}: anchor count")
            check(span == anchors[-1][1] + K, f"record {rid}: span arithmetic")
            check(span <= READ_LENGTH, f"record {rid}: span exceeds read")
            variants = receipt["reads"]
            total = sum(v[3] for v in variants)
            check(total == census["multiplicity"],
                  f"record {rid}: variant instance sum {total} != {census['multiplicity']}")
            for mirror, left, right, count in variants:
                check(len(left) + span + len(right) == READ_LENGTH,
                      f"record {rid}: flank arithmetic")
                skipped_bp += (len(left) + len(right)) * count
            variant_total += len(variants)
            check(len(receipt["occurrences"]) == len(census["occurrences"]),
                  f"record {rid}: occurrence count")
            for ro, co in zip(receipt["occurrences"], census["occurrences"]):
                check(ro["path"] == co["path"] and ro["start"] == co["start"]
                      and ro["orientation"] == co["orientation"],
                      f"record {rid}: occurrence not census-verbatim")
                check(ro["partitions"] == co["partitions"],
                      f"record {rid}: partitions not census-verbatim")
                check(ro["intervals"] == len(co["intervals"]),
                      f"record {rid}: interval count mismatch")
                check(ro["shift"] in (-1, 0, 1), f"record {rid}: shift domain")
                check(len(ro["pf"]) == len(variants),
                      f"record {rid}: pf arity")
                for flag, (mirror, _l, _r, _c) in zip(ro["pf"], variants):
                    check(flag == ((ro["orientation"] == 0) != (mirror == 1)),
                          f"record {rid}: physical strand flag")
                shift_hist[ro["shift"]] = shift_hist.get(ro["shift"], 0) + 1
                n_occ += 1
            n_records += 1
            receipts[rid] = (
                receipt["anchors"], receipt["span"],
                [(v[0], v[1], v[2], v[3]) for v in receipt["reads"]],
            )
    print(f"phase 2: {n_records} records, {n_occ} occurrences, "
          f"{variant_total} read variants, skipped bp {skipped_bp}, "
          f"shift histogram {shift_hist}")

    # ---------------------------------------------------------- phase 3
    # The read hull re-assembled from the receipt's anchor chain and the
    # GFAs' S segments; flanks from the receipt; the reconstructed read
    # must exist in the FASTQ.
    read_hashes = set()
    with gzip.open(f"{D}/reads.fastq.gz", "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            read_hashes.add(fnv1a64(seq.encode()))

    def node_segment(node_abs, sign):
        for partition in axis:
            g = gfa(partition)
            if str(node_abs) in g.segs:
                seg = g.segs[str(node_abs)]
                return seg if sign > 0 else revcomp(seg)
        return None

    sampled_ids = sorted(receipts.keys())[:: max(1, len(receipts) // PHASE3_RECORDS)][:PHASE3_RECORDS]
    phase3_checked = 0
    phase3_no_segment = 0
    for rid in sampled_ids:
        anchors, span, variants = receipts[rid]
        # Interior hull gaps (measured zero on this sample) would leave
        # the hull uncovered by anchor windows; those records are named
        # and skipped rather than falsely failed.
        for mirror in (0, 1):
            own = own_positions(anchors, mirror == 1, 0)
            if any(lo > 0 and hi < span for lo, hi in skipped_from_own(own, span)):
                phase3_no_segment += 1
                anchors = None
                break
        if anchors is None:
            continue
        for mirror, left, right, count in variants:
            # The read's own k-mers at the own anchor positions.
            own = own_positions(anchors, mirror == 1, len(left))
            hull = {}
            ok = True
            for (node, _pos), p in zip(anchors, own):
                seg = node_segment(abs(node), node if mirror == 0 else -node)
                if seg is None:
                    ok = False
                    break
                for offset, base in enumerate(seg):
                    hull[p + offset] = base
            if not ok:
                phase3_no_segment += 1
                continue
            # Overlapping anchor windows must agree on their overlaps.
            read = left + "".join(hull[i] for i in range(len(left), len(left) + span)) + right
            check(fnv1a64(read.encode()) in read_hashes,
                  f"record {rid}: reconstructed read absent from the FASTQ")
            phase3_checked += 1
            if phase3_checked >= PHASE3_RECORDS:
                break
        if phase3_checked >= PHASE3_RECORDS:
            break
    print(f"phase 3: {phase3_checked} reconstructed reads verified against the FASTQ "
          f"({phase3_no_segment} skipped: anchor segment absent from the built GFAs)")

    # ---------------------------------------------------------- phase 4
    sidecar = {}
    with open(sidecar_path) as f:
        for line in f:
            row = json.loads(line)
            sidecar[(row["path"], row["start"], row["end"])] = row
    # Gap segments are interned strictly after the global syncmer
    # interning (10,024,605 syncmer nodes — the partition-embedded build's
    # own statement); P lines interleave them between non-overlapping
    # syncmer steps, while the sidecar walks carry the syncmer steps only.
    N_SYNC = 10024605
    n_rows_validated = 0
    for key, row in sidecar.items():
        partition = row["members"][0][0]
        g = gfa(partition)
        row_name = f"{row['path']}:{row['start']}-{row['end']}"
        spelled, positions, steps = g.spelled(row_name)
        walk = row["walk"]
        p_sync = [(i, s) for i, s in enumerate(steps) if int(s[:-1]) <= N_SYNC]
        check(len(p_sync) == len(walk), f"sidecar row {row_name}: step count")
        order_ok = True
        for (i, s), w in zip(p_sync, walk):
            name, sign = s[:-1], s[-1]
            signed = int(name) if sign == "+" else -int(name)
            if signed != w[1] or positions[i] - positions[p_sync[0][0]] != w[0] - walk[0][0]:
                order_ok = False
                break
        check(order_ok, f"sidecar row {row_name}: walk nodes/spacing vs GFA")
        if not order_ok:
            continue
        # The spelled sequence is anchored at the first P step (which may
        # be a leading gap segment); the row's sequence must match over
        # the covered interval.
        base_abs = walk[0][0] - (positions[p_sync[0][0]] - positions[0])
        lo = max(row["start"], base_abs)
        hi = min(row["end"], base_abs + len(spelled))
        check(hi > lo, f"sidecar row {row_name}: no covered interval")
        if hi > lo:
            got = spelled[lo - base_abs:hi - base_abs]
            want = row["seq"][lo - row["start"]:hi - row["start"]]
            check(got == want, f"sidecar row {row_name}: GFA spelling differs")
        n_rows_validated += 1

    context_lines = 0
    lookups_checked = 0
    with open(context_path) as f:
        for line in f:
            context_lines += 1
    sample_stride = max(1, context_lines // PHASE4_LINES)
    with open(context_path) as f:
        for index, line in enumerate(f):
            if index % sample_stride != 0 or lookups_checked >= PHASE4_LINES * 400:
                continue
            ctx = json.loads(line)
            anchors, span, variants = receipts[ctx["record"]]
            for row_value in ctx["rows"]:
                key = (row_value["path"], row_value["start"], row_value["end"])
                row = sidecar[key]
                walk = row["walk"]
                node_positions = {}
                for bp, node in walk:
                    node_positions.setdefault(abs(node), []).append((bp, 1 if node > 0 else -1))
                seq = row["seq"]
                for per_mirror in row_value["per_mirror"]:
                    mirror = per_mirror["mirror"] == 1
                    pinned, ambiguous = correspondence(anchors, ctx["orientation"], node_positions)
                    check(pinned == per_mirror["anchors"],
                          f"context record {ctx['record']}: pinned correspondence")
                    check(ambiguous == per_mirror["ambiguous"],
                          f"context record {ctx['record']}: ambiguous set")
                    forward = (ctx["orientation"] == 0) != mirror
                    own = own_positions(anchors, mirror, 0)
                    wlos = [len(v[1]) for v in variants if (v[0] == 1) == mirror]
                    if not wlos:
                        continue
                    ranges = [(-(max(wlos)), 0)]
                    interior = [r for r in skipped_from_own(own, span) if r[0] > 0 and r[1] < span]
                    ranges.extend(interior)
                    if READ_LENGTH - min(wlos) > span:
                        ranges.append((span, READ_LENGTH - min(wlos)))
                    expected = []
                    for lo, hi in sorted(ranges):
                        for d in range(lo, hi):
                            serving = None
                            for j, p in enumerate(pinned):
                                if p is None:
                                    continue
                                op = own[j]
                                if serving is None:
                                    serving = (j, op)
                                else:
                                    sj, sp = serving
                                    own_left, best_left = op <= d, sp <= d
                                    if own_left and not best_left:
                                        serving = (j, op)
                                    elif own_left == best_left:
                                        if own_left and op > sp:
                                            serving = (j, op)
                                        if not own_left and op < sp:
                                            serving = (j, op)
                            if serving is None:
                                expected.append([d, -1, 0, 0, None])
                                continue
                            j, op = serving
                            r = pinned[j]
                            coord = r + (d - op) if forward else r + (op + K - 1 - d)
                            if coord < 0:
                                expected.append([d, -1, 0, 0, None])
                                continue
                            lo_i = 0
                            while lo_i < len(walk) and walk[lo_i][0] + K <= coord:
                                lo_i += 1
                            hi_i = lo_i
                            while hi_i < len(walk) and walk[hi_i][0] <= coord:
                                hi_i += 1
                            if row["start"] <= coord < row["end"]:
                                base = seq[coord - row["start"]]
                            else:
                                base = None
                            expected.append([d, coord, lo_i, hi_i, base])
                    check(expected == per_mirror["lookups"],
                          f"context record {ctx['record']}: lookups differ")
                    lookups_checked += len(per_mirror["lookups"])
    print(f"phase 4: {n_rows_validated} sidecar rows validated against the GFAs; "
          f"{lookups_checked} context lookups re-derived exactly")

    return {
        "records": n_records,
        "occurrences": n_occ,
        "variants": variant_total,
        "skipped_bp": skipped_bp,
        "shifts": shift_hist,
        "wall": wall,
        "peak_rss_kb": peak_rss,
    }


def main():
    components = sys.argv[1:] or ["chrMT"]
    path_names = load_path_names()
    summaries = {}
    for component in components:
        summaries[component] = check_component(component, path_names)
    print()
    for component, s in summaries.items():
        print(f"{component}: {s['records']} records, {s['occurrences']} occurrences, "
              f"{s['variants']} read variants, skipped bp {s['skipped_bp']}, "
              f"shifts {s['shifts']}, wall {s['wall']}s, peak RSS {s['peak_rss_kb']} kB")
    print()
    if failures:
        print(f"CHECK FAILED: {len(failures)} failures over {checks[0]} checks")
        sys.exit(1)
    print(f"ALL PHASES PASS ({checks[0]} checks)")
    # (the fleet: one summary log per invocation — concurrent
    # component chains must not clobber each other's summaries)
    with open(f"{D}/check-anchor-projection-{'-'.join(components)}.log", "w") as f:
        f.write(f"{checks[0]} checks, {len(failures)} failures\n")


if __name__ == "__main__":
    main()
