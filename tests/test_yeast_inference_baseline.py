"""Tiny-fixture tests; never open an AGC or invoke a real impg binary."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "yeast_inference_baseline.py"
spec = importlib.util.spec_from_file_location("yeast_inference_baseline", SCRIPT)
baseline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(baseline)


class BaselineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.names = self.write("panel.syng.names", "0\tS#1#contig\t25\t1\t0\t2\t3\n1\tplain path with description\t40\n")
        self.records = baseline.load_names(self.names)
        self.seeds = baseline.make_seeds(self.records, "S#1#contig", 0, 10, 3)

    def write(self, name, text):
        path = self.root / name
        path.write_text(text, encoding="utf-8")
        return path

    def cli(self, argv):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return baseline.main(argv)

    def test_names_formats_exact_names_and_pansn(self):
        names = self.write("names", "4\tS#1#ctg\t9\n7\tplain name\t12\t0\t0\t0\n"
                           "8\tS#1#ctg description\t3\t0\t0\t0\t0\n"
                           "9\tS##ctg\t4\n10\tS#1#ctg#extra\t8\n11\tS#2\t6\n")
        records = baseline.load_names(names)
        self.assertEqual(records[0]["source_id"], 4)
        self.assertEqual(records[0]["pansn"], {"sample": "S", "haplotype": "1", "contig": "ctg"})
        self.assertTrue(all(r["pansn"] is None for r in records[1:]))
        self.assertEqual(records[2]["name"], "S#1#ctg description")
        inventory = baseline.inventory_data(records)
        self.assertEqual(inventory["source_bases"], 42)
        self.assertEqual(inventory["parsed_sample_count"], 1)
        self.assertEqual(inventory["legacy_anchor_offset_names"], ["plain name"])
        self.assertFalse(inventory["source_ids_match_row_order"])
        self.assertEqual(inventory["validation_checks"]["ploidy"], "pending")

    def test_names_reject_duplicates_and_invalid_fields(self):
        for text in ("0\ta\t1\n0\tb\t2\n", "0\ta\t1\n1\ta\t2\n", "00\ta\t1\n0\tb\t2\n",
                     "0\ta\t-1\n", "0\ta\t0\n", "0\ta\t1.5\n", "0\ta\tNaN\n",
                     "0\ta\t18446744073709551616\n", "4294967296\ta\t2\n", "-1\ta\t2\n",
                     "0\t\t3\n", "0\ta\n", "", "0\ta\t3\textra\n"):
            with self.subTest(text=text), self.assertRaises(baseline.BaselineError):
                baseline.load_names(self.write("bad.names", text))

    def test_tail_windows_and_explicit_bound(self):
        self.assertEqual([(s["start"], s["end"]) for s in self.seeds], [(0, 10), (10, 20), (20, 25)])
        self.assertEqual(len(baseline.make_seeds(self.records, "S#1#contig", 0, 10, 2)), 2)
        self.assertEqual([(s["start"], s["end"]) for s in baseline.make_seeds(self.records, "S#1#contig", 24, 10, 100)], [(24, 25)])
        for ref, start, size, count in (("S", 0, 10, 1), ("S#1#contig", -1, 10, 1),
                                        ("S#1#contig", 25, 10, 1), ("S#1#contig", 0, 0, 1),
                                        ("S#1#contig", 0, 10, 0)):
            with self.subTest(start=start, size=size, count=count), self.assertRaises(baseline.BaselineError):
                baseline.make_seeds(self.records, ref, start, size, count)

    def test_repeated_occurrences_labels_strands_and_coverage_retained(self):
        bed = self.write("hits.bed", "plain path with description\t3\t13\tseed_000000\t.\t-\n" * 2
                         + "plain path with description\t20\t30\tseed_000000\t0\t+\n"
                         + "plain path with description\t3\t13\tseed_000001\t.\t-\n")
        summary, rows = baseline.summarize_data(bed, self.records, self.seeds, 10)
        self.assertEqual(len(rows), 4)
        self.assertEqual(len({r["occurrence_id"] for r in rows}), 4)
        self.assertEqual([r["start"] for r in rows], [3, 3, 20, 3])
        first = summary["per_seed"][0]
        self.assertEqual(first["candidate_count"], 3)
        self.assertEqual(first["unique_path_count"], 1)
        self.assertEqual(first["strand_counts"], {"+": 1, "-": 2})
        self.assertEqual(first["length_sum"], 30)
        self.assertEqual(first["length_min"], 10)
        self.assertEqual(summary["source_coverage"][1]["candidate_union_bases"], 20)
        self.assertEqual(summary["source_coverage"][1]["candidate_uncovered_bases"], 20)
        self.assertEqual(summary["source_coverage"][0]["seed_union_bases"], 25)
        self.assertEqual(summary["per_seed"][2]["failures"], ["no_candidates"])

    def test_refined_bed_without_seed_is_valid_and_does_not_invent_seed_hit(self):
        bed = self.write("refined.bed", "plain path with description\t1\t11\tseed_000000\t.\t+\n")
        summary, rows = baseline.summarize_data(bed, self.records, self.seeds[:1], 10)
        seed = summary["per_seed"][0]
        self.assertEqual(seed["failures"], [])
        self.assertEqual(seed["exact_seed_interval_count"], 0)
        self.assertEqual((seed["name"], seed["start"], seed["end"]), ("S#1#contig", 0, 10))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["name"], "plain path with description")
        self.assertIn("not required", summary["seed_membership_policy"])
        self.assertEqual(summary["validation_checks"]["spanning_support"], "pending")

    def test_malformed_bed_unknown_labels_names_and_bounds(self):
        bad_rows = ["diagnostic on stdout", "S#1#contig\t0\t1\tseed_000000",
                    "S#1#contig\t0\t1\tunknown\t.\t+", "unknown\t0\t1\tseed_000000\t.\t+",
                    "S#1#contig\t0\t26\tseed_000000\t.\t+", "S#1#contig\t-1\t1\tseed_000000\t.\t+",
                    "S#1#contig\t1\t1\tseed_000000\t.\t+", "S#1#contig\t2\t1\tseed_000000\t.\t+",
                    "S#1#contig\tx\t1\tseed_000000\t.\t+", "S#1#contig\t0\t1\tseed_000000\t.\t?",
                    "S#1#contig\t0\t1\tseed_000000\t1001\t+", "S#1#contig\t0\t1\tseed_000000\t.\t+\textra"]
        for row in bad_rows:
            with self.subTest(row=row), self.assertRaises(baseline.BaselineError):
                baseline.summarize_data(self.write("bad.bed", row + "\n"), self.records, self.seeds, 10)

    def test_hash_prefixed_exact_source_name_is_not_a_comment(self):
        names = self.write("hash.names", "0\t#literal\t10\n")
        records = baseline.load_names(names)
        seeds = baseline.make_seeds(records, "#literal", 0, 5, 1)
        bed = self.write("hash.bed", "# ordinary comment\n#literal\t0\t5\tseed_000000\t.\t+\n")
        summary, rows = baseline.summarize_data(bed, records, seeds, 10)
        self.assertEqual(summary["candidate_count"], 1)
        self.assertEqual(rows[0]["name"], "#literal")

    def test_candidate_row_limit(self):
        bed = self.write("hits.bed", "S#1#contig\t0\t1\tseed_000000\t.\t+\n" * 2)
        with self.assertRaisesRegex(baseline.BaselineError, "limit exceeded"):
            baseline.summarize_data(bed, self.records, self.seeds, 1)

    def test_no_hits_and_seed_only_diagnostics(self):
        summary, rows = baseline.summarize_data(self.write("empty.bed", "# comment\n\n"), self.records, self.seeds, 10)
        self.assertEqual(rows, [])
        self.assertTrue(all(s["failures"] == ["no_candidates"] for s in summary["per_seed"]))
        self.assertIsNone(summary["per_seed"][0]["length_min"])
        summary, _ = baseline.summarize_data(self.write("self.bed", "S#1#contig\t0\t10\tseed_000000\t.\t+\n"), self.records, self.seeds, 10)
        self.assertEqual(summary["per_seed"][0]["failures"], ["only_exact_seed_intervals"])

    def test_seed_bed_duplicate_labels_rejected_occurrences_allowed(self):
        seed_bed = self.write("seeds.bed", "S#1#contig\t0\t10\tone\nS#1#contig\t0\t10\ttwo\n")
        self.assertEqual(len(baseline.load_seeds(seed_bed, self.records)), 2)
        seed_bed.write_text("S#1#contig\t0\t10\tone\nS#1#contig\t10\t20\tone\n")
        with self.assertRaises(baseline.BaselineError):
            baseline.load_seeds(seed_bed, self.records)

    def fake_pilot(self, mode="ok", timeout="2"):
        for suffix in (".1gbwt", ".1khash", ".spos", ".pstep", ".meta"):
            self.write("panel.syng" + suffix, "fixture, not a real index\n")
        agc = self.write("panel.agc", "fixture, never parsed\n")
        executable = self.write("fake impg", f'''#!{sys.executable}
import json, pathlib, sys, time
root = pathlib.Path(__file__).parent
with (root / "calls.jsonl").open("a") as f:
    f.write(json.dumps(sys.argv[1:]) + "\\n")
mode = {mode!r}
if sys.argv[1:] == ["--version"]:
    print("impg fake-test-version")
    sys.exit(3 if mode == "version_failure" else 0)
assert sys.argv[1] == "query"
print("fixture query log", file=sys.stderr, flush=True)
if mode == "timeout":
    time.sleep(10)
if mode == "failure":
    print("partial invalid output", flush=True)
    sys.exit(7)
if mode == "malformed":
    print("not BED")
elif mode != "empty":
    seeds = pathlib.Path(sys.argv[sys.argv.index("-b") + 1]).read_text().splitlines()
    for line in seeds:
        print(line + "\\t.\\t+")
    print("plain path with description\\t1\\t11\\tseed_000000\\t.\\t-")
    print("plain path with description\\t1\\t11\\tseed_000000\\t.\\t-")
''')
        executable.chmod(0o755)
        return ["pilot", "--index-prefix", str(self.root / "panel.syng"), "--agc", str(agc),
                "--impg", str(executable), "--reference", "S#1#contig", "--start", "0",
                "--window-size", "10", "--max-windows", "3", "--threads", "2",
                "--merge-distance", "0", "--timeout", timeout, "--out-dir", str(self.root / "out")]

    def test_pilot_one_batched_query_and_saved_reproducibility(self):
        self.assertEqual(self.cli(self.fake_pilot()), 0)
        out = self.root / "out"
        calls = [json.loads(line) for line in (self.root / "calls.jsonl").read_text().splitlines()]
        self.assertEqual([c[0] for c in calls], ["--version", "query"])
        self.assertEqual(calls[1][calls[1].index("-o") + 1], "bed")
        self.assertEqual(calls[1][calls[1].index("-t") + 1], "2")
        self.assertIn("--consider-strandness", calls[1])
        manifest = json.loads((out / "manifest.json").read_text())
        self.assertEqual(manifest["argv"][1:], calls[1])
        self.assertEqual(manifest["status"], "succeeded")
        self.assertEqual(manifest["version"], "impg fake-test-version")
        self.assertEqual(manifest["validation_checks"]["spanning_support"], "pending")
        self.assertEqual(len((out / "seeds.bed").read_text().splitlines()), 3)
        self.assertIn("S#1#contig\t20\t25\tseed_000002", (out / "seeds.bed").read_text())
        self.assertEqual(len((out / "candidates.tsv").read_text().splitlines()), 6)
        for file in ("argv.json", "command.sh", "version.txt", "version.stderr.log", "query.stderr.log", "inventory.json"):
            self.assertTrue((out / file).is_file(), file)

    def test_failed_query_rejects_partial_output(self):
        self.assertEqual(self.cli(self.fake_pilot("failure")), 2)
        out = self.root / "out"
        manifest = json.loads((out / "manifest.json").read_text())
        self.assertEqual(manifest["query_returncode"], 7)
        self.assertEqual(manifest["status"], "failed")
        self.assertFalse((out / "summary.json").exists())
        failure = json.loads((out / "failure.json").read_text())
        self.assertEqual(len(failure["per_seed"]), 3)
        self.assertIsNone(failure["per_seed"][0]["candidate_count"])
        self.assertIn("partial", (out / "candidates.bed").read_text())

    def test_timed_out_query_is_failure(self):
        self.assertEqual(self.cli(self.fake_pilot("timeout", "0.3")), 2)
        failure = json.loads((self.root / "out/failure.json").read_text())
        self.assertIn("timed out", failure["error"])
        self.assertFalse((self.root / "out/summary.json").exists())

    def test_version_failure_prevents_query(self):
        self.assertEqual(self.cli(self.fake_pilot("version_failure")), 2)
        self.assertEqual(len((self.root / "calls.jsonl").read_text().splitlines()), 1)

    def test_successful_command_with_malformed_bed_is_failure(self):
        self.assertEqual(self.cli(self.fake_pilot("malformed")), 2)
        manifest = json.loads((self.root / "out/manifest.json").read_text())
        self.assertEqual(manifest["status"], "failed")
        self.assertEqual(manifest["stage"], "output_validation")
        self.assertFalse((self.root / "out/candidates.tsv").exists())

    def test_empty_successful_query_reports_no_hits_not_process_failure(self):
        self.assertEqual(self.cli(self.fake_pilot("empty")), 0)
        summary = json.loads((self.root / "out/summary.json").read_text())
        self.assertEqual(summary["candidate_count"], 0)
        self.assertEqual(summary["per_seed"][0]["failures"], ["no_candidates"])

    def test_nonempty_directory_refused_before_executable(self):
        args = self.fake_pilot()
        (self.root / "out").mkdir()
        sentinel = self.write("out/keep.txt", "unchanged")
        self.assertEqual(self.cli(args), 2)
        self.assertEqual(sentinel.read_text(), "unchanged")
        self.assertFalse((self.root / "calls.jsonl").exists())

    def test_missing_executable_and_invalid_numeric_options(self):
        args = self.fake_pilot()
        args[args.index("--impg") + 1] = str(self.root / "missing")
        self.assertEqual(self.cli(args), 2)
        for option, value in (("--threads", "0"), ("--merge-distance", "-1"),
                              ("--timeout", "nan"), ("--timeout", "0")):
            args = self.fake_pilot()
            args[args.index(option) + 1] = value
            with self.subTest(option=option, value=value):
                self.assertEqual(self.cli(args), 2)
        self.assertFalse((self.root / "out").exists())

    def test_inventory_and_external_summary_cli_and_overwrite(self):
        out = self.root / "inventory"
        args = ["inventory", "--names", str(self.names), "--out-dir", str(out)]
        self.assertEqual(self.cli(args), 0)
        self.assertTrue((out / "sequences.tsv").is_file())
        self.assertEqual(self.cli(args), 2)
        seeds = self.write("seeds.bed", "S#1#contig\t0\t10\tone\n")
        bed = self.write("hits.bed", "plain path with description\t2\t7\tone\t.\t-\n")
        args = ["summarize", "--names", str(self.names), "--seeds", str(seeds), "--bed", str(bed),
                "--out-dir", str(self.root / "summary")]
        self.assertEqual(self.cli(args), 0)
        summary = json.loads((self.root / "summary/summary.json").read_text())
        self.assertEqual(summary["command_status"], "not_verified_external_bed")
        self.assertEqual(summary["per_seed"][0]["length_mean"], 5)
        self.assertEqual(self.cli(args), 2)

    def test_external_summary_validation_failure_records_all_seed_failures(self):
        seeds = self.write("seeds.bed", "S#1#contig\t0\t10\tone\n")
        bed = self.write("bad.bed", "S#1#contig\t0\t26\tone\t.\t+\n")
        args = ["summarize", "--names", str(self.names), "--seeds", str(seeds), "--bed", str(bed),
                "--out-dir", str(self.root / "summary")]
        self.assertEqual(self.cli(args), 2)
        failure = json.loads((self.root / "summary/failure.json").read_text())
        self.assertEqual(failure["per_seed"][0]["failures"], ["invalid_output"])
        self.assertFalse((self.root / "summary/summary.json").exists())

    def test_pilot_rejects_ids_that_impg_would_remap(self):
        args = self.fake_pilot()
        self.names.write_text("4\tS#1#contig\t25\n")
        self.assertEqual(self.cli(args), 2)
        self.assertFalse((self.root / "calls.jsonl").exists())

    def test_empty_existing_output_directory_is_allowed(self):
        out = self.root / "empty-out"
        out.mkdir()
        self.assertEqual(self.cli(["inventory", "--names", str(self.names), "--out-dir", str(out)]), 0)
        self.assertTrue((out / "inventory.json").exists())

    def test_interleaved_directory_acquisition_has_only_one_owner(self):
        out = self.root / "race"
        mkdir = Path.mkdir
        acquired = []

        def interleave(path, *args, **kwargs):
            # Both callers observed the directory absent before this mkdir.
            with mock.patch.object(Path, "mkdir", mkdir):
                acquired.append(baseline.fresh_directory(out))
            return mkdir(path, *args, **kwargs)

        with mock.patch.object(Path, "mkdir", interleave):
            with self.assertRaisesRegex(baseline.BaselineError, "reserved"):
                baseline.fresh_directory(out)
        self.assertEqual(acquired, [out])
        self.assertTrue((out / ".baseline-reserved").is_file())

    def test_partial_candidate_write_never_publishes_summary(self):
        args = self.fake_pilot()
        original = baseline.write_tsv

        def fail_candidates(path, columns, rows):
            if path.name == "candidates.tsv.tmp":
                path.write_text("incomplete candidate data")
                raise OSError("injected TSV write failure")
            return original(path, columns, rows)

        with mock.patch.object(baseline, "write_tsv", fail_candidates):
            self.assertEqual(self.cli(args), 2)
        out = self.root / "out"
        self.assertEqual(json.loads((out / "manifest.json").read_text())["status"], "failed")
        self.assertTrue((out / "failure.json").exists())
        for name in ("summary.json", "candidates.tsv", "candidates.tsv.tmp"):
            self.assertFalse((out / name).exists(), name)

    def test_external_summary_candidate_write_failure_is_not_published(self):
        seeds = self.write("seeds.bed", "S#1#contig\t0\t10\tone\n")
        bed = self.write("hits.bed", "S#1#contig\t1\t9\tone\t.\t+\n")
        out = self.root / "summary"
        args = ["summarize", "--names", str(self.names), "--seeds", str(seeds),
                "--bed", str(bed), "--out-dir", str(out)]
        with mock.patch.object(baseline, "write_tsv", side_effect=OSError("injected TSV failure")):
            self.assertEqual(self.cli(args), 2)
        self.assertTrue((out / "failure.json").exists())
        self.assertFalse((out / "summary.json").exists())
        self.assertFalse((out / "candidates.tsv").exists())

    def test_inventory_write_failure_records_preparation_stage(self):
        args = self.fake_pilot()
        with mock.patch.object(baseline, "write_inventory", side_effect=OSError("injected inventory failure")):
            self.assertEqual(self.cli(args), 2)
        failure = json.loads((self.root / "out/failure.json").read_text())
        self.assertEqual(failure["stage"], "preparation")
        self.assertIn("inventory failure", failure["error"])

    def test_final_manifest_failure_discards_accepted_artifacts(self):
        args = self.fake_pilot()
        original = baseline.write_json

        def fail_final_manifest(path, data):
            if path.name == "manifest.json" and data["status"] == "succeeded":
                raise OSError("injected final manifest failure")
            return original(path, data)

        with mock.patch.object(baseline, "write_json", fail_final_manifest):
            self.assertEqual(self.cli(args), 2)
        self.assertFalse((self.root / "out/summary.json").exists())
        self.assertFalse((self.root / "out/candidates.tsv").exists())

    def test_json_write_failure_preserves_previous_complete_version(self):
        output = self.write("record.json", '{"old": true}\n')
        with mock.patch.object(Path, "replace", side_effect=OSError("injected rename failure")):
            with self.assertRaises(OSError):
                baseline.write_json(output, {"new": True})
        self.assertEqual(json.loads(output.read_text()), {"old": True})
        self.assertFalse((self.root / "record.json.tmp").exists())

    def test_symlink_output_refused(self):
        target = self.root / "target"
        target.mkdir()
        link = self.root / "link"
        link.symlink_to(target, target_is_directory=True)
        with self.assertRaises(baseline.BaselineError):
            baseline.fresh_directory(link)


if __name__ == "__main__":
    unittest.main()
