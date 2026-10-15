#!/usr/bin/env python3
"""Full-shard packet multiset hash: the commutative sum (mod 2^128) of
blake2b-128 digests over every packet of a shard file. Order-independent,
so the send-in-par scheduling's within-batch nondeterminism cancels;
equal sums over equal packet counts = multiset identity (with negligible
collision probability, stated). One process per shard file; run under
GNU parallel/xargs for the fleet.

usage: shard-hash.py <shard.bin.gz>...
prints: <path> <packets> <records> <sum-hi-hex> <sum-lo-hex>
"""
import gzip
import hashlib
import struct
import sys


def shard_hash(path):
    lo = 0
    hi = 0
    packets = 0
    records = 0
    with gzip.open(path, "rb") as f:
        while True:
            head = f.read(4)
            if len(head) < 4:
                break
            (plen,) = struct.unpack("<I", head)
            pkt = f.read(plen)
            if len(pkt) < plen:
                raise SystemExit(f"truncated packet in {path}")
            (seq_len,) = struct.unpack_from("<I", pkt, 0)
            (nrec,) = struct.unpack_from("<I", pkt, 4 + seq_len)
            d = hashlib.blake2b(head + pkt, digest_size=16).digest()
            a = int.from_bytes(d[:8], "little")
            b = int.from_bytes(d[8:], "little")
            lo = (lo + a) & 0xFFFFFFFFFFFFFFFF
            hi = (hi + b) & 0xFFFFFFFFFFFFFFFF
            packets += 1
            records += nrec
    return packets, records, hi, lo


for path in sys.argv[1:]:
    packets, records, hi, lo = shard_hash(path)
    print(f"{path} packets {packets} records {records} "
          f"sum {hi:016x}{lo:016x}")
