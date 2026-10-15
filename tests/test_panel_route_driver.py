"""Portable stdlib regressions for the parent's independent stream validator."""
import importlib.util
import json
from pathlib import Path
import struct
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("panel_route_gate", Path(__file__).resolve().parents[1] / "scripts/panel-route-gate.py")
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


class NativeProfileStreamTest(unittest.TestCase):
    def fixture(self, root, zero=False):
        rows = [] if zero else [{"start": 0, "end": 2, "counts": [[[6, 21, 10], 2]]},
                               {"start": 2, "end": 5, "counts": [[[6, 21, 10], 2]]}]
        events = b""
        index = b""
        for row in rows:
            index += struct.pack("<QQ", row["start"], len(events))
            events += json.dumps(row).encode() + b"\n"
        totals = b"" if zero else b"[[6, 21, 10], 10]\n"
        profile = {"length": 6, "starts": 0 if zero else 5, "runs": len(rows)}
        for key, data in [("events", events), ("index", index), ("totals", totals)]:
            (root / key).write_bytes(data)
            profile[key] = {"path": key, "bytes": len(data), "hash": "unused-by-independent-check"}
        return profile

    def test_valid_complete_and_zero_start_streams(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            gate.validate_native_profile(root, 10, 6, self.fixture(root))
            gate.validate_native_profile(root, 5, 6, self.fixture(root, zero=True))

    def test_malformed_streams_are_not_certified_by_consistent_metadata(self):
        mutations = ["gap", "overlap", "bad-index-start", "interior-offset", "wrong-offset",
                     "truncated-index", "truncated-event", "trailing-index", "trailing-event",
                     "trailing-blank", "wrong-run-count", "missing-final-newline", "bad-total",
                     "missing-total", "duplicate-total", "empty-nonzero", "short-domain",
                     "beyond-domain", "zero-length-run"]
        for mutation in mutations:
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as d:
                root = Path(d)
                profile = self.fixture(root)
                rows = (root / "events").read_bytes().splitlines(keepends=True)
                if mutation in ("gap", "overlap"):
                    row = json.loads(rows[1]); row["start"] = 3 if mutation == "gap" else 1
                    rows[1] = json.dumps(row).encode() + b"\n"
                    (root / "events").write_bytes(b"".join(rows))
                    (root / "index").write_bytes(struct.pack("<QQQQ", 0, 0, row["start"], len(rows[0])))
                elif mutation in ("short-domain", "beyond-domain"):
                    row = json.loads(rows[1]); row["end"] = 4 if mutation == "short-domain" else 6
                    rows[1] = json.dumps(row).encode() + b"\n"
                    (root / "events").write_bytes(b"".join(rows))
                    (root / "totals").write_bytes(json.dumps([[6, 21, 10], row["end"] * 2]).encode() + b"\n")
                elif mutation == "zero-length-run":
                    zero = b'{"start":0,"end":0,"counts":[]}\n'
                    (root / "events").write_bytes(zero + b"".join(rows))
                    (root / "index").write_bytes(struct.pack("<QQQQQQ", 0, 0, 0, len(zero), 2, len(zero) + len(rows[0])))
                    profile["runs"] = 3; profile["index"]["bytes"] = 48
                elif mutation == "bad-index-start":
                    (root / "index").write_bytes(struct.pack("<QQQQ", 0, 0, 3, len(rows[0])))
                elif mutation in ("interior-offset", "wrong-offset"):
                    (root / "index").write_bytes(struct.pack("<QQQQ", 0, 0, 2, 1 if mutation == "interior-offset" else 0))
                elif mutation == "truncated-index": (root / "index").write_bytes((root / "index").read_bytes()[:-1])
                elif mutation == "truncated-event": (root / "events").write_bytes(b"".join(rows)[:-8])
                elif mutation == "trailing-index": (root / "index").write_bytes((root / "index").read_bytes() + b"x")
                elif mutation == "trailing-event": (root / "events").write_bytes(b"".join(rows) + rows[-1])
                elif mutation == "trailing-blank": (root / "events").write_bytes(b"".join(rows) + b"\n")
                elif mutation == "wrong-run-count": profile["runs"] = 1; profile["index"]["bytes"] = 16
                elif mutation == "missing-final-newline": (root / "events").write_bytes(b"".join(rows)[:-1])
                elif mutation == "bad-total": (root / "totals").write_bytes(b"[[6, 21, 10], 11]\n")
                elif mutation == "missing-total": (root / "totals").write_bytes(b"")
                elif mutation == "duplicate-total": (root / "totals").write_bytes((root / "totals").read_bytes() * 2)
                else: (root / "events").write_bytes(b"")
                # Keep byte-size metadata consistent for JSON mutations, so the
                # regression exercises contents, not just file-length sanity.
                for key in ("events", "totals"):
                    profile[key]["bytes"] = (root / key).stat().st_size
                with self.assertRaises((RuntimeError, ValueError, EOFError)):
                    gate.validate_native_profile(root, 10, 6, profile)

    def test_zero_start_domain_rejects_any_records(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            profile = self.fixture(root, zero=True)
            (root / "events").write_bytes(b'{"start":0,"end":1,"counts":[]}\n')
            profile["events"]["bytes"] = (root / "events").stat().st_size
            with self.assertRaises(RuntimeError):
                gate.validate_native_profile(root, 5, 6, profile)

    def test_exact_sparse_totals_verification_fails_instead_of_truncating_at_cap(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            profile = self.fixture(root)
            events, index = b"", b""
            for line in (root / "events").read_bytes().splitlines():
                row = json.loads(line)
                row["counts"].append([[6, 23, 10], 1])
                index += struct.pack("<QQ", row["start"], len(events))
                events += json.dumps(row).encode() + b"\n"
            for key, data in [("events", events), ("index", index),
                              ("totals", b"[[6,21,10],10]\n[[6,23,10],5]\n")]:
                (root / key).write_bytes(data)
                profile[key]["bytes"] = len(data)
            with self.assertRaisesRegex(RuntimeError, "resource cap"):
                gate.validate_native_profile(root, 10, 6, profile, max_terms=1)
            gate.validate_native_profile(root, 10, 6, profile, max_terms=2)

    def test_gate_distinguishes_same_source_cross_source_and_cross_identity(self):
        result = {"native_initializations_complete": True, "sequence_emission_authorized": False,
                  "non_native_assignments_evaluated": 5, "mixed_assignments_evaluated": 0,
                  "mixed_identity_assignments_evaluated": 0}
        with self.assertRaises(RuntimeError):
            gate.validate_search_gate(result)
        result["mixed_assignments_evaluated"] = 1
        gate.validate_search_gate(result)  # source paths may share an identity
        self.assertEqual(result["mixed_identity_assignments_evaluated"], 0)
        result["mixed_identity_assignments_evaluated"] = 1
        gate.validate_search_gate(result)

    def test_actual_source_access_dependencies_are_separate_inputs(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            fasta = root / "plain.fa"
            bgzf = root / "bgzf.fa.gz"
            fasta.write_bytes(b">s\nACGT\n")
            bgzf.write_bytes(b"\x1f\x8b" + b"fixture")
            self.assertEqual(gate.source_access_paths([fasta, bgzf, root / "archive.agc"]),
                             [str(fasta) + ".fai", str(bgzf) + ".fai", str(bgzf) + ".gzi"])
            self.assertFalse(Path(str(fasta) + ".fai").exists())


if __name__ == "__main__":
    unittest.main()
