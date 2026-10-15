#!/usr/bin/env python3
"""Parent-owned, fresh-only all-panel route gate. No truth or sequence emission.

Without --execute this prints a plan and opens no biological inputs. Execution
requires frozen SHA256 declarations and a new output directory. Python stdlib only.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import time

F7_SHA = "a79200bd1ea1c4b3eacbe71effbe10adb6d61ec9904b381250c0d4052f8d8be4"
MODEL = "panel-route-relative-poisson-v1"


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    temporary = path.with_suffix(path.suffix + ".incomplete")
    with open(temporary, "x") as f:
        json.dump(value, f, indent=2)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    temporary.replace(path)


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def source_access_paths(sources):
    """Mirror actual HTSlib lookup dependencies without creating/adopting indexes."""
    result = []
    for source in sources:
        path = Path(source).resolve()
        if path.suffix == ".agc":
            continue
        require(str(path).endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")), "unsupported source access format")
        result.append(str(path) + ".fai")
        with open(path, "rb") as f:
            if f.read(2) == b"\x1f\x8b":
                result.append(str(path) + ".gzi")
    return result


def validate_search_gate(result):
    require(result["native_initializations_complete"], "native family evaluation gate incomplete; increase only declared work resources")
    require(result["mixed_assignments_evaluated"] > 0, "search gate incomplete: no complete route crosses source path IDs")
    require(not result["sequence_emission_authorized"], "unexpected sequence authorization")


def validate_native_profile(root, source_length, length, profile, max_terms=2000000):
    """Independently read every event/index pair and exactly reconcile sparse totals.

    Memory is bounded to one profile's max_terms sparse tokens, not all sources or
    read locations. Limits fail the gate rather than truncate verification.
    """
    u64_max = (1 << 64) - 1
    integer = lambda x: type(x) is int and 0 <= x <= u64_max
    require(integer(length) and length > 0 and integer(source_length) and max_terms > 0, "invalid profile validation parameters")
    starts = max(0, source_length - length + 1)
    require(integer(profile["length"]) and profile["length"] == length
            and integer(profile["starts"]) and profile["starts"] == starts
            and integer(profile["runs"]) and profile["index"]["bytes"] == profile["runs"] * 16,
            "incomplete native start-domain index")
    files = {}
    for key in ("events", "index", "totals"):
        path = Path(profile[key]["path"])
        require(not path.is_absolute() and path.parts and all(p not in (".", "..") for p in path.parts), "invalid native profile path")
        files[key] = Path(root) / path
        require(files[key].stat().st_size == profile[key]["bytes"], "native stream byte size mismatch")

    def unique_object(pairs):
        require(len({key for key, _ in pairs}) == len(pairs), "duplicate native JSON field")
        return dict(pairs)

    def read_row(f):
        line = f.readline(67108865)
        require(line and len(line) <= 67108864 and line.endswith(b"\n"), "truncated/oversized native JSON record")
        return json.loads(line, object_pairs_hook=unique_object)

    def token_count(row):
        require(isinstance(row, list) and len(row) == 2, "invalid native count row")
        token, count = row
        require(isinstance(token, list) and len(token) == 3 and all(integer(t) for t in token)
                and 4 <= token[0] <= 8589934590 and 4 <= token[2] <= 8589934590
                and token[0] % 2 == token[2] % 2 == 0
                and 3 <= token[1] <= 8589934591 and token[1] % 2 == 1 and integer(count) and count > 0,
                "invalid native token/count")
        def flip(t):
            z = (t - 2) // 2
            return 2 + 2 * (z + 1 if z % 2 else z - 1)
        require(token <= [flip(token[2]), token[1], flip(token[0])], "noncanonical native token")
        return tuple(token), count

    totals = {}
    cursor = 0
    with open(files["index"], "rb") as index, open(files["events"], "rb") as events:
        for _ in range(profile["runs"]):
            record = index.read(16)
            require(len(record) == 16, "truncated native index record")
            start, offset = struct.unpack("<QQ", record)
            require(start == cursor and offset == events.tell(), "native index start/record boundary offset mismatch")
            event = read_row(events)
            require(isinstance(event, dict) and set(event) == {"start", "end", "counts"}
                    and integer(event["start"]) and integer(event["end"])
                    and event["start"] == cursor and cursor < event["end"] <= starts
                    and isinstance(event["counts"], list), "native event gap/overlap/domain mismatch")
            previous = None
            for row in event["counts"]:
                token, count = token_count(row)
                require(previous is None or previous < token, "duplicate/unsorted event tokens")
                previous = token
                if token not in totals:
                    require(len(totals) < max_terms, "independent profile verification resource cap exhausted")
                total = totals.get(token, 0) + count * (event["end"] - cursor)
                require(total <= u64_max, "native total overflow")
                totals[token] = total
            cursor = event["end"]
        require(cursor == starts and index.read(1) == b"" and events.read(1) == b"", "incomplete/trailing native event or index data")
    with open(files["totals"], "rb") as f:
        previous = None
        while f.peek(1):
            token, count = token_count(read_row(f))
            require(previous is None or previous < token, "duplicate/unsorted native totals")
            require(totals.pop(token, None) == count, "native integrated totals mismatch")
            previous = token
    require(not totals, "missing native integrated totals")
    return {"runs_verified": profile["runs"], "starts_verified": starts, "totals_verified": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--impg", required=True)
    parser.add_argument("--expected-binary-sha256", required=True)
    parser.add_argument("--f7-finite", required=True)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--sources", nargs="+", required=True)
    parser.add_argument("--inventory", required=True)
    parser.add_argument("--frozen-inputs", required=True,
                        help="JSON object mapping absolute input/frozen artifact paths to SHA256")
    parser.add_argument("--sample", nargs=3, action="append", required=True,
                        metavar=("LABEL", "MEMBWT", "FIXED_DEPTH"))
    parser.add_argument("--read-lengths", nargs="+", type=int, required=True)
    parser.add_argument("--background", type=float, default=0.1)
    parser.add_argument("--max-work", type=int, default=2000000)
    parser.add_argument("--max-evaluations", type=int, default=10000)
    parser.add_argument("--max-frontier", type=int, default=1024)
    parser.add_argument("--max-optima", type=int, default=1000)
    parser.add_argument("--max-profile-terms", type=int, default=2000000)
    parser.add_argument("--max-feature-terms", type=int, default=50000000)
    parser.add_argument("--cache-terms", type=int, default=1000000)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    prefix = ["taskset", "-c", "252-255", "nice", "-n", "10", "/usr/bin/time", "-v"]
    plan = {"scope": "parent-only all-panel generation and fixed-parameter evaluation/search",
            "binary": args.impg, "panel": args.panel, "out_dir": args.out_dir,
            "read_lengths": args.read_lengths, "samples": args.sample,
            "prefix": prefix, "truth": "not accepted", "sequence_emission": "never",
            "relative_objective": "fixed beta/depth composite Poisson, not an exact probabilistic likelihood"}
    if not args.execute:
        print(json.dumps(plan, indent=2))
        return
    root = Path(args.out_dir).absolute()
    root.mkdir(exist_ok=False)
    state = {"status": "running", "plan": plan, "stages": []}
    write(root / "driver-status.json", state)
    try:
        require(sha(args.impg) == args.expected_binary_sha256, "worker binary SHA mismatch")
        require(sha(args.f7_finite) == F7_SHA, "frozen independent f7 binary mismatch")
        with open(args.frozen_inputs) as f:
            frozen = json.load(f)
        require(isinstance(frozen, dict) and frozen, "empty frozen input declarations")
        needed = [args.catalog, args.inventory, args.impg, args.f7_finite, *args.sources]
        needed += [row[1] for row in args.sample]
        # Auxiliary bindings must already be frozen. Never create an index here
        # and silently adopt it as a new input to cached native computations.
        access_paths = source_access_paths(args.sources)
        needed += access_paths
        # Resolve exactly the primary/legacy panel sidecar convention.
        for suffix in ("1gbwt", "1khash", "names", "meta", "spos", "pstep"):
            primary = f"{args.panel}.{suffix}" if suffix.startswith("1") or args.panel.endswith(".syng") else f"{args.panel}.syng.{suffix}"
            needed.append(primary if Path(primary).exists() else f"{args.panel}.syng.{suffix}")
        require(all(str(Path(p).resolve()) in frozen for p in needed), "missing frozen input SHA declaration")
        for path, digest in frozen.items():
            require(Path(path).is_absolute() and sha(path) == digest, f"frozen input mismatch: {path}")
        write(root / "frozen-inputs.json", frozen)
        with open(args.inventory) as f:
            inventory = json.load(f)
        expected = {entry["source"]: (family["identity"], entry["path"], entry["length"])
                    for family in inventory["families"] for entry in family["paths"]}
        require(len(expected) == inventory["path_count"] and not inventory["zero_length_paths"], "invalid parent topology inventory")
        env = os.environ.copy()
        env.update(RAYON_NUM_THREADS="4", CARGO_BUILD_JOBS="4",
                   LD_LIBRARY_PATH="/home/erikg/.cache/impg/native-build-deps/root/usr/lib/x86_64-linux-gnu",
                   IMPG_TEST_WFMASH="/home/erikg/.cargo/bin/wfmash")

        def run(name, command):
            begin = time.time()
            state["stages"].append({"name": name, "command": prefix + command, "status": "running"})
            write(root / "driver-status.json", state)
            with open(root / f"{name}.log", "xb") as log:
                result = subprocess.run(prefix + command, stdout=log, stderr=subprocess.STDOUT, env=env)
            state["stages"][-1].update(status="succeeded" if result.returncode == 0 else "failed",
                                      seconds=time.time() - begin, returncode=result.returncode)
            write(root / "driver-status.json", state)
            require(result.returncode == 0, f"stage failed; preserved {name}.log")

        graph_dir = root / "routes"
        run("build", [args.impg, "genome-infer", "build-panel-routes", "--panel", args.panel,
                      "--catalog", args.catalog, "--sources", *args.sources,
                      "--read-lengths", ",".join(map(str, args.read_lengths)),
                      "--max-profile-terms", str(args.max_profile_terms), "--out-dir", str(graph_dir)])
        with open(graph_dir / "manifest.json") as f:
            manifest = json.load(f)
        require(manifest["status"] == "succeeded" and manifest["model"] == MODEL, "incomplete route graph")
        with open(graph_dir / "graph.json") as f:
            graph = json.load(f)
        require([b["path"] for b in graph["source_access"]] == access_paths, "source-access bindings differ from frozen inputs")
        for path in access_paths:
            require(sha(path) == frozen[str(Path(path).resolve())], f"source-access sidecar changed during build: {path}")
        actual = {lane["id"]: (graph["families"][lane["family"]]["identity"], lane["name"], lane["length"])
                  for lane in graph["lanes"]}
        require(actual == expected, "native paths/identities/lengths missing or relabeled")
        require(graph["total_source_bp"] == inventory["source_bp"] and len(graph["families"]) == inventory["family_count"], "family/bp inventory mismatch")
        expected_families = {f["identity"]: [p["source"] for p in f["paths"]] for f in inventory["families"]}
        require({f["identity"]: f["paths"] for f in graph["families"]} == expected_families, "endpoint path inventory/pairing mismatch")
        require(graph["generation_complete"] and all(len(l["sequence_fnv1a64"]) == 16 for l in graph["lanes"]), "incomplete native source spelling scan")
        # Stream exact-once normalized ownership; never load original location maps.
        ends = [0] * len(actual)
        with open(graph_dir / graph["ownership"]["path"]) as f:
            for line in f:
                row = json.loads(line)
                require(row["start"] == ends[row["source"]] and row["end"] > row["start"], "ownership gap/overlap")
                ends[row["source"]] = row["end"]
        require(ends == [actual[i][2] for i in range(len(actual))], "uncovered native source bases")
        k = graph["k"]
        width = k + 17
        def read_port(f):
            data = f.read(width)
            if not data:
                return None
            require(len(data) == width, "partial port record")
            return (data[:k], *struct.unpack("<QQB", data[k:]))
        # Independent streaming membership identity: sorted global hub stream is
        # ordered/unique; per-lane pairs retain every opposite traversal description.
        total, per_source, previous = 0, [0] * len(actual), None
        with open(graph_dir / graph["ports"]["path"], "rb") as f:
            while (port := read_port(f)) is not None:
                word, source, anchor, reverse = port
                require(source in actual and anchor + k <= actual[source][2] and reverse in (0, 1), "invalid physical port")
                require(previous is None or previous < port, "duplicate/unsorted hub membership")
                previous = port
                per_source[source] += 1
                total += 1
        complement = bytes.maketrans(b"ACGT", b"TGCA")
        profiles_verified, runs_verified, starts_verified = 0, 0, 0
        for lane in graph["lanes"]:
            count, previous_anchor = 0, -1
            with open(graph_dir / lane["ports"]["path"], "rb") as f:
                while (forward := read_port(f)) is not None:
                    reverse = read_port(f)
                    require(reverse is not None and forward[1:] == (lane["id"], forward[2], 0)
                            and reverse[1:] == (lane["id"], forward[2], 1)
                            and forward[2] > previous_anchor and reverse[0] == forward[0].translate(complement)[::-1], "physical-view dedup/reverse port mismatch")
                    previous_anchor = forward[2]
                    count += 2
            require(count == lane["port_count"] == per_source[lane["id"]], "incomplete lane hub membership")
            require(len(lane["profiles"]) == len(graph["read_lengths"]), "incomplete native length inventory")
            for length, profile in zip(graph["read_lengths"], lane["profiles"]):
                verified = validate_native_profile(graph_dir, lane["length"], length, profile, args.max_profile_terms)
                profiles_verified += 1
                runs_verified += verified["runs_verified"]
                starts_verified += verified["starts_verified"]
        require(total == graph["port_count"], "global port enumeration mismatch")
        write(root / "generation-gate.json", {"native_paths": len(actual), "source_bp": sum(ends),
              "families": len(graph["families"]), "oriented_ports": total,
              "zero_anchor_paths": sum(l["port_count"] == 0 for l in graph["lanes"]),
              "native_profiles_verified": profiles_verified, "native_event_runs_verified": runs_verified,
              "native_starts_verified": starts_verified, "native_totals_independently_verified": True,
              "generation_rule_complete": True, "biological_topology_complete": False})
        labels = set()
        for label, sample_path, depth in args.sample:
            require(label and label not in labels and all(c.isalnum() or c in "-_" for c in label), "invalid/duplicate sample label")
            labels.add(label)
            output = root / f"search-{label}"
            run(f"search-{label}", [args.impg, "genome-infer", "search-panel-routes", "--panel", args.panel,
                "--routes", str(graph_dir), "--sample", sample_path, "--haploid-depth", depth,
                "--background", str(args.background), "--max-work", str(args.max_work),
                "--max-evaluations", str(args.max_evaluations), "--max-frontier", str(args.max_frontier),
                "--max-optima", str(args.max_optima), "--max-feature-terms", str(args.max_feature_terms),
                "--cache-terms", str(args.cache_terms), "--out-dir", str(output)])
            with open(output / "result.json") as f:
                result = json.load(f)["result"]
            validate_search_gate(result)
            run(f"evaluate-{label}", [args.impg, "genome-infer", "evaluate-panel-routes", "--panel", args.panel,
                "--routes", str(graph_dir), "--sample", sample_path, "--haploid-depth", depth,
                "--background", str(args.background), "--assignment", str(output / "incumbent-assignment.json"),
                "--max-feature-terms", str(args.max_feature_terms), "--out-dir", str(root / f"evaluate-{label}")])
            with open(root / f"evaluate-{label}" / "result.json") as f:
                evaluation = json.load(f)["result"]
            require(abs(evaluation["relative_objective"] - result["incumbent"]["relative_objective"]) <= 1e-8 * max(1, abs(evaluation["relative_objective"])), "independent process route evaluation mismatch")
            # Bound support remains honest even when this operational gate succeeds.
            state.setdefault("search_guarantees", {})[label] = {k: result[k] for k in (
                "status", "search_exhausted", "global_optimum_certified", "correlated_optima_complete",
                "native_initializations_complete", "mixed_assignments_evaluated",
                "mixed_identity_assignments_evaluated", "non_native_assignments_evaluated", "lower_bound")}
        for path, digest in frozen.items():
            require(sha(path) == digest, f"frozen input changed during gate: {path}")
        state.update(status="succeeded-operational-gates-only", retained_bytes=sum(p.stat().st_size for p in root.rglob("*") if p.is_file()),
                     genome_coverage_repair_proven=False, sequence_emission_authorized=False)
    except Exception as exc:
        state.update(status="failed-or-incomplete", error=str(exc))
        write(root / "driver-status.json", state)
        raise
    write(root / "driver-status.json", state)


if __name__ == "__main__":
    main()
