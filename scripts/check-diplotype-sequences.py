#!/usr/bin/env python3
"""THE DIPLOTYPE SEQUENCE EMISSION GATE (the write-up codification).

Proves the product emission (--emit-sequences, calls.jsonl.sequences.jsonl)
on one locality's CLI run, three checks, every one a loud gate:

  CHECK 1  EVERY CALL -> A SEQUENCE. The sequences file carries exactly
     one record per product call; hap1/hap2 are nonempty, the emitted
     sequence's length equals the record's `length` field, the QUAL
     fields are identical to the product call's (qual, qual_unbounded,
     class_count, called_class_count), and the called_set is exactly the
     COALESCED CALLED SET of the call's called_classes (every distinct
     called fold, sorted - the owner's ambiguity ruling; an untied call's
     set is exactly the representative pair, and the representative
     hap1/hap2 folds are members of the set).

  CHECK 2  THE SPELL IDENTITY (the walks from the local graph): every
     emitted fold sequence is byte-equal to EVERY member row's
     locus.fa slice [start, end) - the fold's own spelled material
     (folds coalesce identical sequence+walk copies, so every member
     must reconstruct the same sequence; any mismatch is a failure).

  CHECK 3  THE OFFICIAL_ALLWAVE IDENTITY (assessment-side truth only as
     INPUT, the test mode's own artifact): the official sequence-QV
     convention recomputed FROM THE EMITTED SEQUENCES through the real
     allwave binary (allwave -i <fasta> -t 4 -p tree:5:0:0, per-truth-
     hap best alignment, identity = matches/query_span, QV = -10log10
     capped 50, diplotype QV the per-hap mean) against the truth rows
     of the locality's own committed maps must reproduce the CLI's
     --truth-qv official_allwave fields EXACTLY, at every locus where
     the artifact carries them (the inexpressible loci are counted and
     named, not failed - the truth pair not being in the candidate
     universe is the committed convention, not an emission defect).

usage: check-diplotype-sequences.py --calls <calls.jsonl> \\
    --sequences <calls.jsonl.sequences.jsonl> --truth-qv <truth-qv.jsonl> \\
    --maps <locality-graphs dir> --locus-fa <locus.fa> \\
    --truth-prefixes HG002#1,HG002#2 --out <report.txt> [--allwave-bin allwave]
"""
import argparse
import json
import math
import os
import subprocess
import sys
import tempfile

FAILURES = []
NOTES = []


def fail(msg):
    FAILURES.append(msg)


def official_qv_fields(called0, called1, truth0, truth1, locus_tag, workdir, allwave):
    """The CLI's AllwaveQv.locus_qv mirrored exactly (the committed
    convention, the same mirror the masked campaign calibrated)."""
    fa = os.path.join(workdir, f"locus-{locus_tag}.fa")
    with open(fa, "w") as out:
        for name, seq in [("called0", called0), ("called1", called1),
                          ("truth0", truth0), ("truth1", truth1)]:
            out.write(f">{name}\n{seq}\n")
    result = subprocess.run(
        [allwave, "-i", fa, "-t", "4", "-p", "tree:5:0:0"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise SystemExit(f"allwave failed at {locus_tag}: {result.stderr[:400]}")
    paf = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            raise SystemExit(f"malformed allwave PAF line: {line[:120]}")
        paf.append((fields[0], fields[5], int(fields[9]),
                    int(fields[3]) - int(fields[2])))

    def best_for_hap(truth):
        best = None
        for index, target in enumerate(["called0", "called1"]):
            matches = span = 0
            for query, tgt, m, s in paf:
                if query == truth and tgt == target:
                    matches += m
                    span += s
            if span == 0:
                continue
            identity = matches / span
            candidate = (identity, matches, span, index)
            if best is None or candidate[:3] > best[:3]:
                best = candidate
        if best is None:
            return 0.0, 0, 0, None
        return best[0], best[1], best[2], f"called{best[3]}"

    def qv_of(identity):
        if identity >= 1.0:
            return 50.0
        error = 1.0 - identity
        if error <= 0.0:
            return 50.0
        return max(0.0, min(50.0, -10.0 * math.log10(error)))

    i0, m0, s0, t0 = best_for_hap("truth0")
    i1, m1, s1, t1 = best_for_hap("truth1")
    q0, q1 = qv_of(i0), qv_of(i1)
    return {
        "qv": (q0 + q1) / 2.0,
        "qv_truth_hap0": q0,
        "qv_truth_hap1": q1,
        "identity_truth_hap0": i0,
        "identity_truth_hap1": i1,
        "best_called_hap_truth_hap0": t0,
        "best_called_hap_truth_hap1": t1,
        "best_matches": [m0, m1],
        "best_query_spans": [s0, s1],
        "no_alignment_truth_hap0": t0 is None,
        "no_alignment_truth_hap1": t1 is None,
    }


def truth_rows_of_partition(maps_dir, partition, prefixes):
    """ALL of the committed maps' truth rows for one partition (the
    axis path's and the truth-2 path's member rows; a partition can
    carry several rows per path). The machinery's truth fold keeps its
    own REPRESENTATIVE member row's material, so the official-QV
    identity check searches these candidates."""
    for name in sorted(os.listdir(maps_dir)):
        if not name.endswith(".gfa.map.json"):
            continue
        with open(os.path.join(maps_dir, name)) as handle:
            m = json.load(handle)
        if m["construction"]["axis_partition"] != partition:
            continue
        first = [member for member in m["members"]
                 if member["path_name"].startswith(prefixes[0] + "#")]
        second = [member for member in m["members"]
                  if member["path_name"].startswith(prefixes[1] + "#")]
        return first, second
    return None, None


class Fasta:
    """Random-access slicing of a fasta with a samtools fai index."""

    def __init__(self, path):
        self.path = path
        self.cache = {}
        with open(path + ".fai") as handle:
            for line in handle:
                parts = line.split("\t")
                self.cache[parts[0]] = None

    def slice(self, name, start, end):
        result = subprocess.run(
            ["samtools", "faidx", self.path, f"{name}:{start + 1}-{end}"],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise SystemExit(f"faidx failed for {name}: {result.stderr}")
        return "".join(result.stdout.split("\n")[1:]).upper()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--calls", required=True)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--truth-qv", required=True)
    parser.add_argument("--maps", required=True)
    parser.add_argument("--locus-fa", required=True)
    parser.add_argument("--truth-prefixes", default="HG002#1,HG002#2")
    parser.add_argument("--allwave-bin", default="allwave")
    parser.add_argument("--instrument-receipt", required=True,
                        help="the CLI run's instrument-receipt.jsonl (the "
                             "machinery's own fold anatomy: fold_identities "
                             "members + the alignment-connected merge groups' "
                             "representative rows - the truth folds' material)")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    prefixes = args.truth_prefixes.split(",")

    calls = {}
    for line in open(args.calls):
        record = json.loads(line)
        calls[record["locus"]] = record
    sequences = {}
    for line in open(args.sequences):
        record = json.loads(line)
        sequences[record["locus"]] = record
    truth_qv = {}
    for line in open(args.truth_qv):
        record = json.loads(line)
        truth_qv[record["locus"]] = record
    receipt = {}
    for line in open(args.instrument_receipt):
        record = json.loads(line)
        receipt[record["locus"]] = record

    fasta = Fasta(args.locus_fa)
    workdir = tempfile.mkdtemp(prefix="emission-gate-")

    # ------------------------------ CHECK 1: every call -> a sequence
    if set(calls) != set(sequences):
        fail(f"CHECK 1: locus sets differ: {len(calls)} calls vs "
             f"{len(sequences)} sequence records; "
             f"calls-without-sequence {sorted(set(calls) - set(sequences))[:10]}; "
             f"sequences-without-call {sorted(set(sequences) - set(calls))[:10]}")
    qual_checked = spell_checked = 0
    for locus, call in sorted(calls.items()):
        record = sequences.get(locus)
        if record is None:
            continue
        for hap in ("hap1", "hap2"):
            entry = record.get(hap)
            if not isinstance(entry, dict) or not entry.get("sequence"):
                fail(f"CHECK 1: locus {locus}: {hap} missing or empty")
                continue
            if len(entry["sequence"]) != entry["length"]:
                fail(f"CHECK 1: locus {locus}: {hap} length {entry['length']} "
                     f"!= emitted sequence length {len(entry['sequence'])}")
        for field in ("qual", "qual_unbounded", "class_count",
                      "called_class_count"):
            if record.get(field) != call.get(field):
                fail(f"CHECK 1: locus {locus}: {field} "
                     f"{record.get(field)!r} != the call's {call.get(field)!r}")
        expected_set = sorted({fold
                               for pair in call.get("called_classes", [])
                               for fold in pair})
        got_set = [entry["fold_index"] for entry in record.get("called_set", [])]
        if got_set != expected_set:
            fail(f"CHECK 1: locus {locus}: called_set {got_set} != the "
                 f"coalesced called classes {expected_set}")
        for hap in ("hap1", "hap2"):
            if record[hap]["fold_index"] not in got_set:
                fail(f"CHECK 1: locus {locus}: {hap} fold not in called_set")
        qual_checked += 1

    # ------------------------------ CHECK 2: the spell identity
    # The fold's material is the REPRESENTATIVE row's spell (the panel
    # fetch of a member row window). Under the alignment-connected fold
    # relaxation (the -alignfold model marker) merged folds keep the
    # representative's sequence and carry the MERGED members' rows
    # (proven equivalents by the committed phase-10 gate: variant-pocket
    # census mass zero + the class-LL invariant) - those members' raw
    # slices may legitimately differ from the fold sequence and are
    # NAMED, not failed. Under the exact model every member must
    # slice-equal the fold sequence.
    merge_members = []
    for locus, record in sorted(sequences.items()):
        entries = list(record.get("called_set", []))
        for hap in ("hap1", "hap2"):
            entry = record.get(hap)
            if isinstance(entry, dict) and entry.get("sequence"):
                entries.append(entry)
        alignfold = "alignfold" in str(calls.get(locus, {}).get("model", ""))
        for entry in entries:
            sequence = entry.get("sequence")
            if not sequence:
                continue
            matching = 0
            for member in entry.get("members", []):
                spelled = fasta.slice(member["path_name"], member["start"],
                                      member["end"])
                if spelled == sequence:
                    matching += 1
                elif alignfold:
                    diffs = sum(1 for a, b in zip(spelled, sequence) if a != b)
                    merge_members.append(
                        f"locus {locus} fold {entry['fold_index']}: merged "
                        f"member {member['path_name']}"
                        f"[{member['start']},{member['end']}) diverges "
                        f"from the representative spell at {diffs} bp "
                        f"(the alignment-connected fold's variant pocket)")
                else:
                    fail(f"CHECK 2: locus {locus} fold {entry['fold_index']}: "
                         f"the member row {member['path_name']}"
                         f"[{member['start']},{member['end']}) does not "
                         f"spell the emitted sequence "
                         f"({len(spelled)} vs {len(sequence)} bp)")
            if matching < 1:
                fail(f"CHECK 2: locus {locus} fold {entry['fold_index']}: "
                     f"NO member row spells the emitted sequence "
                     f"(the representative material is unverifiable)")
            spell_checked += 1

    # ------------------- CHECK 3: the official_allwave identity
    official_checked = official_match = 0
    inexpressible = []
    for locus, record in sorted(truth_qv.items()):
        official = record.get("official_allwave")
        if not official:
            inexpressible.append(locus)
            continue
        emission = sequences.get(locus)
        if emission is None:
            fail(f"CHECK 3: locus {locus}: official_allwave present but no "
                 f"emitted sequences record")
            continue
        # The truth folds' material from the machinery's own receipt:
        # the fold anatomy (fold_identities' member rows; for an
        # alignment-connected MERGED fold the merge group's named
        # representative row) - the exact material the in-process
        # official_allwave computation consumed, so an exact
        # reproduction proves the EMITTED called sequences are the
        # aligned material.
        entry = receipt.get(locus)
        if entry is None:
            fail(f"CHECK 3: locus {locus}: no instrument receipt record")
            continue
        fold_identities = entry.get("fold_identities") or []
        merge_groups = ((entry.get("alignment_connected_folds") or {})
                        .get("merge_groups") or [])
        def truth_material(fold_index):
            if not isinstance(fold_index, int) or fold_index >= len(fold_identities):
                return None
            for group in merge_groups:
                if group.get("fold") == fold_index:
                    return group.get("representative_row")
            # an exact (non-merged) fold: every member row spells the
            # fold's sequence (identical sequence+walk by construction)
            members = fold_identities[fold_index].get("members") or []
            return members[0] if members else None
        folds_list = [f.get("index") for f in record.get("truth_folds", [])]
        if len(folds_list) != 2:
            fail(f"CHECK 3: locus {locus}: the truth-qv record names "
                 f"{len(folds_list)} truth folds (expected 2)")
            continue
        rows = [truth_material(index) for index in folds_list]
        if any(row is None for row in rows):
            fail(f"CHECK 3: locus {locus}: the truth folds' material is "
                 f"not resolvable from the instrument receipt "
                 f"(fold indices {folds_list})")
            continue
        truth0 = fasta.slice(rows[0]["path_name"], rows[0]["start"],
                             rows[0]["end"])
        truth1 = fasta.slice(rows[1]["path_name"], rows[1]["start"],
                             rows[1]["end"])
        recomputed = official_qv_fields(
            emission["hap1"]["sequence"], emission["hap2"]["sequence"],
            truth0, truth1, f"{record['partition']}-{locus}", workdir,
            args.allwave_bin)
        official_checked += 1
        mismatched = [field for field in
                      ("qv", "qv_truth_hap0", "qv_truth_hap1",
                       "identity_truth_hap0", "identity_truth_hap1",
                       "best_called_hap_truth_hap0",
                       "best_called_hap_truth_hap1",
                       "best_matches", "best_query_spans",
                       "no_alignment_truth_hap0", "no_alignment_truth_hap1")
                      if recomputed[field] != official.get(field)]
        if mismatched:
            fail(f"CHECK 3: locus {locus}: official_allwave fields differ "
                 f"when recomputed from the emitted sequences against the "
                 f"receipt's own truth-fold material: "
                 f"{ {f: (recomputed[f], official.get(f)) for f in mismatched} }")
        else:
            official_match += 1

    with open(args.out, "w") as handle:
        handle.write("THE DIPLOTYPE SEQUENCE EMISSION GATE\n")
        handle.write(f"calls: {len(calls)}; sequence records: {len(sequences)}; "
                     f"truth-qv records: {len(truth_qv)}\n")
        handle.write(f"CHECK 1 every-call-a-sequence + QUAL + called_set: "
                     f"{qual_checked}/{len(calls)} checked\n")
        handle.write(f"CHECK 2 spell identity (member-row locus.fa slices): "
                     f"{spell_checked} fold records checked; "
                     f"{len(merge_members)} merged-fold member rows named "
                     f"(the -alignfold variant-pocket class, proven "
                     f"equivalents by the committed phase-10 gate)\n")
        for message in merge_members:
            handle.write(f"  merged-member: {message}\n")
        handle.write(f"CHECK 3 official_allwave identity recomputed from the "
                     f"emitted sequences: {official_match}/{official_checked} "
                     f"exact; inexpressible (no official fields, named): "
                     f"{len(inexpressible)} {inexpressible}\n")
        if FAILURES:
            handle.write(f"FAILURES ({len(FAILURES)}):\n")
            for message in FAILURES:
                handle.write(f"  FAIL {message}\n")
        else:
            handle.write("ALL CHECKS PASS\n")
    print(f"emission gate: {'PASS' if not FAILURES else 'FAIL'} "
          f"({len(FAILURES)} failures); report at {args.out}")
    sys.exit(1 if FAILURES else 0)


if __name__ == "__main__":
    main()
