#!/usr/bin/env python3
"""THE IDENTITY GATE (the campaign's absolute gate).

Compares shard packet streams as RECORD/MULTISET identities:

  mode committed:  slice-run shard  vs  the committed projection's shard
                   PREFIX (first N packets, N = the slice's per-stream
                   read count). The committed whole-sample shard's first
                   N packets are exactly the same stream's first N reads
                   (the slice extraction replays the committed samtools
                   stream command verbatim), so the multisets must agree
                   record-for-record.

  mode runs:       two runs' shard files (baseline vs optimized), packet
                   multiset AND (if present) the per-shard .occ.bin
                   pattern-multiset artifacts.

Byte-identity of the packet ORDER is NOT the standard: the committed
send-in-par scheduling sends packets to the shard writer from inside the
rayon parallel body, so within-batch order is the scheduler's own
nondeterminism (receipted at commit 51a730a). The records themselves
(seq + canonical tokens + (node,pos) walk) are the census source and must
be identical.
"""
import gzip
import hashlib
import struct
import sys
from collections import Counter


def packet_stream(path, limit=None):
    """Yield raw packet byte strings from a shard .bin.gz, up to limit."""
    n = 0
    with gzip.open(path, "rb") as f:
        while limit is None or n < limit:
            head = f.read(4)
            if len(head) < 4:
                return
            (plen,) = struct.unpack("<I", head)
            payload = f.read(plen)
            if len(payload) < plen:
                raise SystemExit(f"truncated packet in {path}")
            yield head + payload
            n += 1


def packet_counts_window(path, offset, count):
    """Counts over the packet window [offset, offset+count)."""
    counts = Counter()
    reads = 0
    bases = 0
    records = 0
    for i, pkt in enumerate(packet_stream(path, limit=offset + count)):
        if i < offset:
            continue
        h = hashlib.blake2b(pkt, digest_size=16).digest()
        counts[h] += 1
        (seq_len,) = struct.unpack_from("<I", pkt, 4)
        (nrec,) = struct.unpack_from("<I", pkt, 8 + seq_len)
        reads += 1
        bases += seq_len
        records += nrec
    return counts, reads, bases, records


def packet_window_sets(path, offsets, count):
    """One streaming pass to max(offsets)+count: the packet multiset +
    counters of each disjoint window [off, off+count), for full-vs-full
    shard comparisons (the optimized shard's window vs the committed
    shard's SAME window - not the whole file vs a window). Windows past
    the end of the file come back empty, identically for equal files."""
    res = {o: [Counter(), 0, 0, 0] for o in offsets}
    hi = max(offsets) + count
    ordered = sorted(offsets)
    cur = 0  # index into ordered of the window containing packet i
    for i, pkt in enumerate(packet_stream(path, limit=hi)):
        while cur < len(ordered) and i >= ordered[cur] + count:
            cur += 1
        if cur >= len(ordered):
            break
        o = ordered[cur]
        if i < o:
            continue
        r = res[o]
        h = hashlib.blake2b(pkt, digest_size=16).digest()
        r[0][h] += 1
        (seq_len,) = struct.unpack_from("<I", pkt, 4)
        (nrec,) = struct.unpack_from("<I", pkt, 8 + seq_len)
        r[1] += 1
        r[2] += seq_len
        r[3] += nrec
    return res


def packet_window_multiset(path, offset, count):
    """Multiset of packet hashes + counters over [offset, offset+count)."""
    return packet_counts_window(path, offset, count)


def packet_counts(path, limit=None):
    counts = Counter()
    reads = 0
    bases = 0
    records = 0
    for pkt in packet_stream(path, limit):
        h = hashlib.blake2b(pkt, digest_size=16).digest()
        counts[h] += 1
        (seq_len,) = struct.unpack_from("<I", pkt, 4)
        (nrec,) = struct.unpack_from("<I", pkt, 8 + seq_len)
        reads += 1
        bases += seq_len
        records += nrec
    return counts, reads, bases, records


def occ_counts(path):
    counts = Counter()
    with open(path, "rb") as f:
        while True:
            head = f.read(4)
            if len(head) < 4:
                break
            (n,) = struct.unpack("<I", head)
            toks = f.read(8 * n)
            if len(toks) < 8 * n:
                raise SystemExit(f"truncated occ entry in {path}")
            (mult_b,) = struct.unpack("<Q", f.read(8))
            counts[toks] += mult_b
    return counts


def diff_label(a, b, name):
    only_a = a - b
    only_b = b - a
    if not only_a and not only_b:
        print(f"  {name}: IDENTICAL ({sum(a.values())} items)")
        return True
    print(f"  {name}: MISMATCH  only-left={sum(only_a.values())} "
          f"distinct {len(only_a)}; only-right={sum(only_b.values())} "
          f"distinct {len(only_b)}")
    for h, c in list(only_a.items())[:3]:
        print(f"    e.g. left-only {h.hex()} x{c}")
    for h, c in list(only_b.items())[:3]:
        print(f"    e.g. right-only {h.hex()} x{c}")
    return False


def guarded_gate(slice_shard, committed_shard, off, n, J, guard_shards):
    """The slice committed-identity gate with the edge-jitter guard.

    The committed send-in-par scheduling scrambles packet order within
    each 16,384-read arrival batch, so a committed shard's packet-window
    edge set is an arbitrary subset of the boundary batches' stream
    reads. Legitimacy of a symmetric difference requires BOTH:
      (a) every window-extra packet (in the committed window, not among
          the slice's core reads) is a packet of a GUARD read - a stream
          read from outside the window, within J of an edge;
      (b) every core read whose packet is missing from the window has its
          packet in the committed shard's own guard bands [off-J, off) and
          [off+n, off+n+J).
    (b) alone proves the difference is the jitter class from the
    committed side; (a) closes it from the read side.
    """
    a, ra, ba, reca = packet_counts(slice_shard)
    w, rb, bb, recb = packet_counts_window(committed_shard, off, n)
    print(f"  slice: reads {ra} bases {ba} records {reca}")
    print(f"  committed window [{off}, {off + n}): reads {rb} bases {bb} records {recb}")
    if rb != n:
        raise SystemExit(f"committed shard window is not {n} packets (got {rb})")
    only_a = a - w
    only_b = w - a
    if not only_a and not only_b:
        print("  packet multiset: IDENTICAL (no edge jitter)")
        return True
    gl = packet_counts_window(committed_shard, max(0, off - J), min(off, J))[0]
    gr = packet_counts_window(committed_shard, off + n, J)[0]
    b_ok = sum((only_a - (gl + gr)).values()) == 0
    g = Counter()
    for path in guard_shards:
        g += packet_counts(path)[0]
    a_ok = sum((only_b - g).values()) == 0
    print(f"  window extras {sum(only_b.values())} distinct {len(only_b)}; "
          f"guard-run confined: {a_ok}")
    print(f"  core reads missing from window {sum(only_a.values())} distinct {len(only_a)}; "
          f"committed-band confined: {b_ok}")
    return a_ok and b_ok


def main():
    mode = sys.argv[1]
    if mode == "committed":
        # args: <slice-shard.bin.gz> <committed-shard.bin.gz> <offset> <N>
        # (N = 0 -> the whole window to end of file)
        #
        # THE EDGE-JITTER GUARD: the committed send-in-par scheduling
        # sends packets to the shard writer from inside the rayon body,
        # so a committed shard's packet order is scrambled within each
        # 16,384-read arrival batch; a packet window's edge set is an
        # arbitrary subset of the boundary batch's stream reads. One
        # arrival batch holds at most 16,384 reads of one stream, so any
        # legitimate symmetric-difference packet must lie within
        # +/-16,384 packets of the window edges (the guard bands).
        sl, com, off, n = sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
        a, ra, ba, reca = packet_counts(sl)
        whole = n == 0
        if whole:
            n = 1 << 62
        b, rb, bb, recb = packet_counts_window(com, off, n)
        print(f"  slice: reads {ra} bases {ba} records {reca}")
        print(f"  committed window [{off}, {off + n}): reads {rb} bases {bb} records {recb}")
        if rb < n and n != 1 << 62:
            raise SystemExit(f"committed shard has fewer than {off + n} packets")
        if whole:
            ok = diff_label(a, b, "packet multiset")
            sys.exit(0 if ok else 1)
        only_a = a - b
        only_b = b - a
        if not only_a and not only_b:
            print("  packet multiset: IDENTICAL (no edge jitter)")
            sys.exit(0)
        J = 16384
        lg, *_ = packet_counts_window(com, max(0, off - J), min(off, J))
        rg, *_ = packet_counts_window(com, off + n, J)
        guard = lg + rg
        confined = sum((only_a - guard).values()) == 0 and sum((only_b - guard).values()) == 0
        print(f"  symmetric difference: {sum(only_a.values())} + {sum(only_b.values())} "
              f"packets; guard bands [+/-{J}]: {confined and 'CONFINED' or 'NOT CONFINED'}")
        print("  (record/multiset identity modulo the committed run's own "
              "within-batch send jitter at the window edges)" if confined else "")
        sys.exit(0 if confined else 1)
    elif mode == "guarded":
        # args: <slice-shard> <committed-shard> <off> <n> <J> <guard-shard>...
        sl, com, off, n, J = sys.argv[2:7]
        guards = sys.argv[7:]
        ok = guarded_gate(sl, com, int(off), int(n), int(J), guards)
        sys.exit(0 if ok else 1)
    elif mode == "windows":
        # args: <shardA.bin.gz> <shardB.bin.gz> <count> <offset> ...
        # FULL-VS-FULL window comparison: the two whole-shard files must
        # carry equal packet counts (the file-level identity is the
        # commutative hash standard); this mode samples equal windows of
        # BOTH files and requires multiset + counter identity per window.
        a_path, b_path, count = sys.argv[2], sys.argv[3], int(sys.argv[4])
        offsets = [int(x) for x in sys.argv[5:]]
        ok = True
        a = packet_window_sets(a_path, offsets, count)
        b = packet_window_sets(b_path, offsets, count)
        for o in offsets:
            ca, ra, ba, reca = a[o]
            cb, rb, bb, recb = b[o]
            if ca == cb and (ra, ba, reca) == (rb, bb, recb):
                print(f"  window [{o}, {o + count}): IDENTICAL "
                      f"({ra} reads {reca} records)")
            else:
                ok = False
                print(f"  window [{o}, {o + count}): MISMATCH "
                      f"reads {ra}/{rb} bases {ba}/{bb} records {reca}/{recb}")
        sys.exit(0 if ok else 1)
    elif mode == "runs":
        # args: <shardA.bin.gz> <shardB.bin.gz> [<occA.bin> <occB.bin>] ...
        ok = True
        it = iter(sys.argv[2:])
        for a_path, b_path in zip(it, it):
            if a_path.endswith(".occ.bin") and b_path.endswith(".occ.bin"):
                continue  # raw occ artifacts: gated in the second loop below
            a, ra, ba, reca = packet_counts(a_path)
            b, rb, bb, recb = packet_counts(b_path)
            print(f"  {a_path} vs {b_path}: reads {ra}/{rb} records {reca}/{recb}")
            if ra != rb or ba != bb or reca != recb:
                ok = False
                print("    READ/RECORD COUNT MISMATCH")
                continue
            ok = diff_label(a, b, "packet multiset") and ok
        it = iter(sys.argv[2:])
        for a_path, b_path in zip(it, it):
            if a_path.endswith(".occ.bin") and b_path.endswith(".occ.bin"):
                ok = diff_label(occ_counts(a_path), occ_counts(b_path),
                                f"occ multiset {a_path.rsplit('/', 1)[-1]}") and ok
        sys.exit(0 if ok else 1)
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
