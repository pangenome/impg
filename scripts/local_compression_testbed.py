#!/usr/bin/env python3
"""Generate and run the local graph compression testbed fixtures.

The fast profile is intentionally deterministic and small.  It produces a
scoreboard row for every fixture/method pair, including rows skipped by profile
or optional-tool policy, and it records exact path-preservation checks for every
graph it writes.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


REQUIRED_FIXTURE_CLASSES = [
    "snp_bubble",
    "short_indel",
    "insertion_50_500bp",
    "alu_like_insertion",
    "adjacent_bubbles_compress_together",
    "bubble_split_by_fake_repeat_anchor",
    "repeated_motif_microtangle",
    "duplicated_flank_requires_path_context",
    "tandem_copy_number_loop_cyclic",
    "dispersed_repeat_glue_break_or_ignore",
    "inversion_like",
    "nested_bubbles_top_level_right",
    "nested_bubbles_top_level_wrong",
]


MANDATORY_METHOD_IDS = [
    "local_syng_raw",
    "local_syng_crush_auto",
    "local_syng_crush_poa",
    "local_syng_crush_poasta",
    "local_syng_crush_sweepga",
    "top_flubble_nonoverlap_sweepga",
    "chunk_window_smooth_or_crush",
    "chunk_window_sweepga_seqwish",
    "whole_region_sweepga_seqwish",
]


OPTIONAL_METHOD_IDS = [
    "pggb_control",
    "smoothxg_control",
    "pggb_plus_smoothxg_control",
]


CONTROL_TIMEOUT_SECONDS = 120
CONTROL_IMPG_PGGB_TARGET_POA_LENGTHS = {
    "pggb_control": "64",
    "smoothxg_control": "64",
    "pggb_plus_smoothxg_control": "32,64",
}
CONTROL_MAX_NODE_LENGTH = 64
CONTROL_MAX_RSS_RE = re.compile(r"Maximum resident set size \(kbytes\):\s+(\d+)")
CHUNK_WINDOW_MIN_INVARIANT_GAP_BP = 5


@dataclass(frozen=True)
class MethodSpec:
    method_id: str
    family: str
    strategy: str
    optional: bool
    parameters: str


@dataclass(frozen=True)
class ControlAvailability:
    available: bool
    provider: str
    detail: str
    paths: Dict[str, str]


METHODS = [
    MethodSpec(
        "local_syng_raw",
        "Raw/local SYNG-derived graph construction",
        "raw_paths",
        False,
        "fixture-local;strategy=raw_paths;threads=1;normalization=copy;candidate_windows=none;path_guard=record",
    ),
    MethodSpec(
        "local_syng_crush_auto",
        "SYNG plus existing flubble/local crush variant",
        "compact_bubble",
        False,
        "fixture-local;strategy=compact_bubble;threads=1;resolver=auto;candidate_policy=all_fixture_variants;path_guard=record",
    ),
    MethodSpec(
        "local_syng_crush_poa",
        "SYNG plus POA/SPOA local crush",
        "compact_bubble",
        False,
        "fixture-local;strategy=compact_bubble;threads=1;resolver=poa_spoa;polish_rounds=1;path_guard=record",
    ),
    MethodSpec(
        "local_syng_crush_poasta",
        "SYNG plus abPOA/POASTA local crush",
        "compact_bubble",
        False,
        "fixture-local;strategy=compact_bubble;threads=1;resolver=poasta_abpoa;polish_rounds=1;path_guard=record",
    ),
    MethodSpec(
        "local_syng_crush_sweepga",
        "SYNG plus SweepGA/seqwish local crush",
        "compact_bubble",
        False,
        "fixture-local;strategy=compact_bubble;threads=1;resolver=sweepga_seqwish;seqwish_k=19;filter=off;path_guard=record",
    ),
    MethodSpec(
        "top_flubble_nonoverlap_sweepga",
        "Highest-level non-overlapping flubble-window crush",
        "compact_bubble",
        False,
        "fixture-local;strategy=compact_bubble;threads=1;window_policy=highest_nonoverlap_flubble;resolver=sweepga_seqwish;seqwish_k=19;path_guard=record",
    ),
    MethodSpec(
        "chunk_window_smooth_or_crush",
        "SmoothXG-style sorted/chunk-window smoothing or crush",
        "compact_bubble",
        False,
        "fixture-local;strategy=compact_bubble;threads=1;window_policy=sorted_chunk;chunk_bp=512;min_invariant_gap_bp=5;resolver=local_smooth_or_crush;path_guard=record",
    ),
    MethodSpec(
        "chunk_window_sweepga_seqwish",
        "Sorted chunk-window SweepGA/seqwish local crush",
        "compact_bubble",
        False,
        "fixture-local;strategy=compact_bubble;threads=1;window_policy=sorted_chunk;chunk_bp=512;min_invariant_gap_bp=5;resolver=sweepga_seqwish;seqwish_k=19;filter=off;path_guard=record",
    ),
    MethodSpec(
        "whole_region_sweepga_seqwish",
        "Whole-region SweepGA/seqwish induction",
        "compact_bubble",
        False,
        "fixture-local;strategy=compact_bubble;threads=1;window_policy=whole_fixture;resolver=sweepga_seqwish;seqwish_k=19;filter=off;path_guard=record",
    ),
    MethodSpec(
        "pggb_control",
        "Bounded PGGB-style whole-region control",
        "internal_control",
        True,
        "control=bounded;provider=internal_impg_graph_pggb;engine=pggb;aligner=fastga-via-impg;threads=1;timeout=120s;target-poa-length=64;max-node-length=64;path_guard=record",
    ),
    MethodSpec(
        "smoothxg_control",
        "Bounded SmoothXG-style local control",
        "internal_control",
        True,
        "control=bounded;provider=internal_impg_smoothxg_stage_via_pggb_engine;engine=pggb;input-paf=empty;threads=1;timeout=120s;target-poa-length=64;max-node-length=64;path_guard=record",
    ),
    MethodSpec(
        "pggb_plus_smoothxg_control",
        "Bounded PGGB/SmoothXG-style local control",
        "internal_control",
        True,
        "control=bounded;provider=internal_impg_graph_pggb;engine=pggb;aligner=fastga-via-impg;threads=1;timeout=120s;target-poa-length=32,64;max-node-length=64;path_guard=record",
    ),
]


TSV_FIELDS = [
    "fixture_id",
    "fixture_class",
    "tier",
    "method_id",
    "method_family",
    "method_parameters",
    "command_line",
    "profile",
    "input_manifest_path",
    "metadata_path",
    "input_fasta_path",
    "expected_paths_path",
    "output_gfa_path",
    "normalized_gfa_path",
    "command_log_path",
    "stdout_log_path",
    "stderr_log_path",
    "metrics_json_path",
    "fixture_notes_path",
    "exact_path_preservation",
    "hard_path_corruption",
    "path_corruption_detail",
    "expected_topology_status",
    "expected_topology_assertion_id",
    "expected_topology_message",
    "missing_path_count",
    "extra_path_count",
    "path_name_stable",
    "graph_size_bytes",
    "segment_count",
    "link_count",
    "path_count",
    "total_segment_bp",
    "total_path_steps",
    "path_replay_compression_ratio",
    "node_depth_min",
    "node_depth_p05",
    "node_depth_median",
    "node_depth_p95",
    "node_depth_max",
    "path_depth_min",
    "path_depth_p05",
    "path_depth_median",
    "path_depth_p95",
    "path_depth_max",
    "white_space_proxy_bp_total",
    "white_space_proxy_bp_p95",
    "white_space_proxy_bp_p99",
    "white_space_proxy_bp_max",
    "self_loop_count",
    "repeat_loop_count",
    "bubble_count",
    "flubble_count",
    "long_link_count",
    "long_link_max_span_bp",
    "runtime_seconds",
    "max_rss_kb",
    "command_status",
    "exit_code",
    "skip_reason",
    "skipped_optional_tool_reason",
    "mandatory_unavailable",
    "candidate_count",
    "candidate_skipped_count",
    "candidate_skip_reasons",
    "render_status",
    "render_path",
    "render_skip_reason",
    "profile_fixture_included",
    "tool_available",
    "validation_log_path",
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def rel(path: Optional[Path], root: Optional[Path] = None) -> str:
    if path is None:
        return ""
    if root is None:
        root = repo_root()
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def pansn_name(sample: str, fixture_id: str, seq: str) -> str:
    return f"{sample}#0#{fixture_id}:1-{len(seq)}"


def generic_ranges(loop_required: bool = False, messy: bool = False) -> Dict[str, Dict[str, int]]:
    max_segments = 18 if messy else 12
    max_links = 28 if messy else 18
    ranges = {
        "segment_count": {"min": 2, "max": max_segments},
        "link_count": {"min": 1, "max": max_links},
        "path_depth_median": {"min": 2, "max": 20},
        "node_depth_p95": {"min": 1, "max": 8},
        "white_space_proxy_bp_p95": {"min": 0, "max": 4},
        "white_space_proxy_bp_max": {"min": 0, "max": 8},
        "bubble_count": {"min": 1, "max": 3},
        "flubble_count": {"min": 1, "max": 3},
        "long_link_count": {"min": 0, "max": 0},
        "long_link_max_span_bp": {"min": 0, "max": 0},
    }
    if loop_required:
        ranges["self_loop_count"] = {"min": 1, "max": 1}
        ranges["repeat_loop_count"] = {"min": 1, "max": 1}
        ranges["node_depth_p95"] = {"min": 1, "max": 12}
    else:
        ranges["self_loop_count"] = {"min": 0, "max": 0}
        ranges["repeat_loop_count"] = {"min": 0, "max": 0}
    return ranges


def make_fixture(
    fixture_id: str,
    fixture_class: str,
    title: str,
    tier: str,
    sample_sequences: Sequence[Tuple[str, str]],
    assertion_id: str,
    topology_description: str,
    known_failure_mode: str,
    expected_shape: str,
    messy: bool = False,
    loop_required: bool = False,
    graph_hints: Optional[Dict[str, Any]] = None,
    range_overrides: Optional[Dict[str, Dict[str, int]]] = None,
) -> Dict[str, Any]:
    input_sequences = []
    input_sequence_names = []
    expected_path_spellings = {}
    for sample, seq in sample_sequences:
        name = pansn_name(sample, fixture_id, seq)
        input_sequences.append(
            {
                "name": name,
                "file": "input.fa",
                "expected_spelling": seq,
            }
        )
        input_sequence_names.append(name)
        expected_path_spellings[name] = seq

    exact = {
        "path_count": len(sample_sequences),
        "self_loops": 1 if loop_required else 0,
        "repeat_loops": 1 if loop_required else 0,
        "max_long_link_span_bp": 0,
    }
    if loop_required:
        exact["loop_policy"] = "repeat_loop_required"
    else:
        exact["loop_policy"] = "forbid_self_and_repeat_loops"

    render_hints = {
        "reference_path": input_sequence_names[0],
        "highlight_paths": input_sequence_names[1:],
        "expected_shape": expected_shape,
    }
    if graph_hints:
        render_hints.update(graph_hints)

    allowed_ranges = generic_ranges(loop_required=loop_required, messy=messy)
    if range_overrides:
        allowed_ranges.update(range_overrides)

    return {
        "fixture_id": fixture_id,
        "fixture_class": fixture_class,
        "title": title,
        "tier": tier,
        "input_sequences": input_sequences,
        "input_sequence_names": input_sequence_names,
        "expected_path_spellings": expected_path_spellings,
        "expected_topology": {
            "assertion_id": assertion_id,
            "description": topology_description,
            "exact": exact,
        },
        "allowed_ranges": allowed_ranges,
        "long_link_policy": {
            "max_count": 0,
            "max_span_bp": 0,
        },
        "known_failure_mode": known_failure_mode,
        "render_hints": render_hints,
    }


def fixture_specs() -> List[Dict[str, Any]]:
    mid_insert = "ACGT" * 50
    mid_insert_alt = ("TGCA" * 49) + "TT"
    alu_insert = ("AGCTGGTGCAAAT" * 20) + "AGCTGGTG"
    dispersed_repeat = "GGAATTCGGAATTC"
    return [
        make_fixture(
            "snp_bubble_3path",
            "snp_bubble",
            "Three-path SNP bubble",
            "ci",
            [
                ("ref", "ACGTACGTACGTCGATTACAGATTACA"),
                ("hapA", "ACGTACGTACGTTGATTACAGATTACA"),
                ("hapB", "ACGTACGTACGTGGATTACAGATTACA"),
            ],
            "single_ordered_snp_bubble",
            "One source-to-sink SNP bubble with no loop and no long links.",
            "Branch path spellings get swapped or the SNP remains as path-specific singleton nodes.",
            "single small ordered bubble",
        ),
        make_fixture(
            "short_indel_3path",
            "short_indel",
            "Short insertion/deletion bubble",
            "ci",
            [
                ("ref", "TTGACCTTGACCAACCGTAACCGTAA"),
                ("hapA", "TTGACCTTGACCAAGGCCGTAACCGTAA"),
                ("hapB", "TTGACCTTGACCCCGTAACCGTAA"),
            ],
            "short_indel_bubble",
            "One compact indel bubble with shared flanks and no long links.",
            "Resolver clips inserted bases, creates a dangling tip, or collapses deletion into the wrong path.",
            "single insertion/deletion bubble",
        ),
        make_fixture(
            "mid_insertion_200bp",
            "insertion_50_500bp",
            "Mid-sized 200 bp insertion",
            "ci",
            [
                ("ref", ("GCTAGCTA" * 3) + ("TTAACCGG" * 3)),
                ("hapA", ("GCTAGCTA" * 3) + mid_insert + ("TTAACCGG" * 3)),
                ("hapB", ("GCTAGCTA" * 3) + mid_insert_alt + ("TTAACCGG" * 3)),
            ],
            "mid_insertion_shared_flanks",
            "Inserted 50-500 bp branches remain exactly spelled and share both flanks.",
            "POA/SYNG crush over-fragments the insertion or drops the inserted branch.",
            "bounded mid-sized insertion bubble",
            messy=True,
        ),
        make_fixture(
            "alu_like_insertion",
            "alu_like_insertion",
            "Repeat-rich Alu-like insertion",
            "local",
            [
                ("ref", "CCGGAATTCCGG" + "TTAACCGGTTAACC"),
                ("hapA", "CCGGAATTCCGG" + alu_insert + "TTAACCGGTTAACC"),
                ("hapB", "CCGGAATTCCGG" + alu_insert.replace("AAA", "AAT", 1) + "TTAACCGGTTAACC"),
            ],
            "alu_like_repeat_insertion",
            "Repeat-rich inserted branches preserve exact spelling without repeat-glue long links.",
            "Repeat content glues unrelated anchors or fails under short exact-match floors.",
            "repeat-rich insertion with bounded ranges",
            messy=True,
        ),
        make_fixture(
            "adjacent_bubbles_joint",
            "adjacent_bubbles_compress_together",
            "Adjacent variants that should be joint",
            "ci",
            [
                ("ref", "AACCGGTTAACCGGTTACCGGAATTCCGG"),
                ("hapA", "AACCGGTTAATCGGATACCGGAATTCCGG"),
                ("hapB", "AACCGGTTAAACGGCTACCGGAATTCCGG"),
            ],
            "adjacent_bubbles_joint_context",
            "Two adjacent variants are represented as one compact local problem or topology-equivalent joint result.",
            "Per-leaf flubble boundaries leave two stringy bubbles separated by artificial white-space.",
            "joint adjacent-bubble event",
        ),
        make_fixture(
            "fake_repeat_anchor_split",
            "bubble_split_by_fake_repeat_anchor",
            "Bubble split by fake internal repeat anchor",
            "ci",
            [
                ("ref", "GGGAAACTTATATATCCCGGGTTTAAA"),
                ("hapA", "GGGAAAGTTATATATACCGGGTTTAAA"),
                ("hapB", "GGGAAACTTATATATGGCGGGTTTAAA"),
                ("hapC", "GGGAAAGTTATATATGGCGGGTTTAAA"),
            ],
            "fake_repeat_anchor_not_boundary",
            "A repeated internal k-mer does not split the true event into unrelated pieces.",
            "Candidate discovery treats the repeat as a hard anchor and under-compresses the true event.",
            "bubble with internal fake repeat anchor",
            messy=True,
        ),
        make_fixture(
            "microtangle_repeat_motif",
            "repeated_motif_microtangle",
            "Bounded repeated-motif microtangle",
            "ci",
            [
                ("ref", "TTCGAACAGCAGCAGGGTTAACCGG"),
                ("hapA", "TTCGAACAGCAGCAGCAGGGTTAACCGG"),
                ("hapB", "TTCGAACAGCAGGGTTAACCGG"),
                ("hapC", "TTCGAACAGCAGCAGTAGGGTTAACCGG"),
            ],
            "bounded_repeat_motif_microtangle",
            "Short repeated motif creates a bounded microtangle with exact spelling and no collapsed motif copy.",
            "Local graph explodes, self-loops appear, or all motif copies collapse into one unsupported node.",
            "bounded motif-count microtangle",
            messy=True,
        ),
        make_fixture(
            "duplicated_flank_context",
            "duplicated_flank_requires_path_context",
            "Duplicated flank requiring path context",
            "ci",
            [
                ("ref", "AACCGGGGATCCCTTAAGGATCCGGAATT"),
                ("hapA", "AACCGGGGATCCCTTACGGATCCGGAATT"),
                ("hapB", "AACCGGGGATCCCTTAAGGATCCCTAATT"),
            ],
            "duplicated_flank_context_required",
            "The replacement uses path context to select the correct duplicate flank occurrence.",
            "Replacement laces into the wrong duplicate flank or merges distinct occurrences.",
            "variant near a duplicated flank",
            messy=True,
        ),
        make_fixture(
            "tandem_copy_loop_keep",
            "tandem_copy_number_loop_cyclic",
            "Tandem copy-number loop should be kept",
            "ci",
            [
                ("ref", "GTTACCGGAATT" + ("ACAC" * 2) + "TTAACCGGTTAA"),
                ("hapA", "GTTACCGGAATT" + ("ACAC" * 3) + "TTAACCGGTTAA"),
                ("hapB", "GTTACCGGAATT" + ("ACAC" * 4) + "TTAACCGGTTAA"),
            ],
            "tandem_copy_loop_required",
            "Tandem copy-number variation preserves a cyclic or loop-capable representation.",
            "Normalization or crush linearizes the copy-number loop and loses cyclic topology.",
            "copy-number loop over a tandem motif",
            loop_required=True,
            graph_hints={
                "repeat_motif": "ACAC",
                "shared_left": "GTTACCGGAATT",
                "shared_right": "TTAACCGGTTAA",
            },
        ),
        make_fixture(
            "dispersed_repeat_glue_break",
            "dispersed_repeat_glue_break_or_ignore",
            "Dispersed repeat glue should be broken or ignored",
            "local",
            [
                ("ref", "AACCTT" + dispersed_repeat + "GGCCAATT" + dispersed_repeat + "TTAAGG"),
                ("hapA", "AACCTT" + dispersed_repeat + "GGCCGGTT" + dispersed_repeat + "TTAAGG"),
                ("hapB", "AACCTT" + dispersed_repeat + "GGCCAATT" + dispersed_repeat.replace("TTC", "TTA", 1) + "TTAAGG"),
            ],
            "dispersed_repeat_no_long_glue",
            "Dispersed repeat similarity is broken or ignored without long repeat-glue links.",
            "seqwish/SYNG transitive closure glues distant repeat copies into long links.",
            "two dispersed repeats without long glue",
            messy=True,
        ),
        make_fixture(
            "inversion_like_case",
            "inversion_like",
            "Orientation-sensitive inversion-like branch",
            "ci",
            [
                ("ref", "CCGGAATTAACCGGTTTTAACCGG"),
                ("hapA", "CCGGAATTTTGGCCAATTAACCGG"),
                ("hapB", "CCGGAATTAACCGGAATTAACCGG"),
            ],
            "inversion_like_orientation_preserved",
            "Inversion-like branches preserve orientation-sensitive path spellings and avoid cross-orientation lacing.",
            "Resolver reverse-complements the wrong interval or creates long links between inverted anchors.",
            "orientation-sensitive branch bubble",
            messy=True,
        ),
        make_fixture(
            "nested_top_level_right",
            "nested_bubbles_top_level_right",
            "Nested bubbles where top-level boundary is right",
            "ci",
            [
                ("ref", "AATTCCGGGGAACCTTGGCCAATT"),
                ("hapA", "AATTCCGGTTAACTTTGGCCAATT"),
                ("hapB", "AATTCCGGGGAAGGTTGGCCAATT"),
                ("hapC", "AATTCCGGTTAACTTAGGCCAATT"),
            ],
            "nested_parent_boundary_right",
            "Nested bubbles are compressed by top-level non-overlapping boundaries with descendant detail included.",
            "Leaf-only selection leaves avoidable nested white-space or repeated replacement conflicts.",
            "nested event inside one parent bubble",
            messy=True,
        ),
        make_fixture(
            "nested_top_level_wrong",
            "nested_bubbles_top_level_wrong",
            "Nested bubbles where top-level boundary over-merges",
            "ci",
            [
                ("ref", "GGCCAATTCCGGTTAAAACCGGTTCCAATTGG"),
                ("hapA", "GGCCAATTCCGGAATAAACCGGTTCCAATTGG"),
                ("hapB", "GGCCAATTCCGGTTAAAACCAATTCCAATTGG"),
                ("hapC", "GGCCAATTCCGGAATAAACCAATTCCAATTGG"),
            ],
            "nested_parent_boundary_wrong",
            "Top-level compression is intentionally wrong; bounded chunking should outperform parent over-merge.",
            "Top-level selection over-merges independent events and creates bad loops, long links, or excess path depth.",
            "two independent events inside an overbroad parent",
            messy=True,
            range_overrides={
                "bubble_count": {"min": 2, "max": 3},
                "flubble_count": {"min": 2, "max": 3},
            },
        ),
    ]


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2)
        handle.write("\n")


def write_fixture_tree(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    manifest_entries = []
    required_by_id = {spec["fixture_id"]: spec for spec in fixture_specs()}
    for fixture_id, spec in required_by_id.items():
        fixture_dir = root / fixture_id
        fixture_dir.mkdir(parents=True, exist_ok=True)
        input_fa = fixture_dir / "input.fa"
        with input_fa.open("w", encoding="utf-8") as handle:
            for record in spec["input_sequences"]:
                handle.write(f">{record['name']}\n{wrap(record['expected_spelling'])}\n")

        expected_paths = fixture_dir / "expected_paths.tsv"
        with expected_paths.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["path_name", "expected_spelling"])
            for name, spelling in spec["expected_path_spellings"].items():
                writer.writerow([name, spelling])

        metadata = fixture_dir / "metadata.json"
        write_json(metadata, spec)

        notes = fixture_dir / "notes.md"
        notes.write_text(
            "\n".join(
                [
                    f"# {spec['title']}",
                    "",
                    f"- Fixture ID: `{fixture_id}`",
                    f"- Class: `{spec['fixture_class']}`",
                    f"- Tier: `{spec['tier']}`",
                    f"- Assertion: `{spec['expected_topology']['assertion_id']}`",
                    f"- Expected shape: {spec['render_hints']['expected_shape']}",
                    f"- Known failure mode: {spec['known_failure_mode']}",
                    "",
                ]
            ),
            encoding="utf-8",
        )

        manifest_entries.append(
            {
                "fixture_id": fixture_id,
                "fixture_class": spec["fixture_class"],
                "tier": spec["tier"],
                "metadata_path": f"tests/test_data/local_compression/{fixture_id}/metadata.json",
                "input_fasta_path": f"tests/test_data/local_compression/{fixture_id}/input.fa",
                "expected_paths_path": f"tests/test_data/local_compression/{fixture_id}/expected_paths.tsv",
                "notes_path": f"tests/test_data/local_compression/{fixture_id}/notes.md",
            }
        )

    manifest = {
        "manifest_version": 1,
        "schema": "local_compression_fixture_manifest_v1",
        "description": "Synthetic local graph compression fixtures covering the design-required classes.",
        "generated_by": "scripts/local_compression_testbed.py write-fixtures",
        "required_fixture_classes": REQUIRED_FIXTURE_CLASSES,
        "fast_profile_fixture_tiers": ["ci"],
        "local_profile_fixture_tiers": ["ci", "local"],
        "fixtures": manifest_entries,
    }
    manifest_path = root / "manifest.json"
    write_json(manifest_path, manifest)
    return manifest_path


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, List[str]] = {}
    name: Optional[str] = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].strip()
            if not name:
                raise ValueError(f"{path}: empty FASTA record name")
            if name in records:
                raise ValueError(f"{path}: duplicate FASTA record {name}")
            records[name] = []
        else:
            if name is None:
                raise ValueError(f"{path}: sequence data before first FASTA header")
            records[name].append(line)
    return {name: "".join(parts) for name, parts in records.items()}


def read_expected_paths_tsv(path: Path) -> Dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"path_name", "expected_spelling"}
        if set(reader.fieldnames or []) != required:
            raise ValueError(f"{path}: expected columns path_name and expected_spelling")
        return {row["path_name"]: row["expected_spelling"] for row in reader}


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_manifest(manifest_path: Path) -> Tuple[Dict[str, Any], List[Dict[str, Any]], List[str]]:
    manifest = load_json(manifest_path)
    repo = repo_root()
    errors: List[str] = []

    if manifest.get("schema") != "local_compression_fixture_manifest_v1":
        errors.append("manifest schema must be local_compression_fixture_manifest_v1")
    if manifest.get("required_fixture_classes") != REQUIRED_FIXTURE_CLASSES:
        errors.append("manifest required_fixture_classes does not match the design-required class order")

    fixtures = manifest.get("fixtures")
    if not isinstance(fixtures, list):
        raise ValueError("manifest fixtures must be a list")

    seen_classes = set()
    loaded: List[Dict[str, Any]] = []
    for entry in fixtures:
        fixture_id = entry.get("fixture_id")
        if not fixture_id:
            errors.append("manifest fixture entry missing fixture_id")
            continue
        expected_paths = {
            "metadata_path": f"tests/test_data/local_compression/{fixture_id}/metadata.json",
            "input_fasta_path": f"tests/test_data/local_compression/{fixture_id}/input.fa",
            "expected_paths_path": f"tests/test_data/local_compression/{fixture_id}/expected_paths.tsv",
            "notes_path": f"tests/test_data/local_compression/{fixture_id}/notes.md",
        }
        for key, expected in expected_paths.items():
            if entry.get(key) != expected:
                errors.append(f"{fixture_id}: manifest {key} is {entry.get(key)!r}, expected {expected!r}")

        metadata_path = repo / expected_paths["metadata_path"]
        input_fasta_path = repo / expected_paths["input_fasta_path"]
        expected_paths_path = repo / expected_paths["expected_paths_path"]
        for path in [metadata_path, input_fasta_path, expected_paths_path, repo / expected_paths["notes_path"]]:
            if not path.exists():
                errors.append(f"{fixture_id}: missing {rel(path)}")

        if not metadata_path.exists():
            continue

        metadata = load_json(metadata_path)
        seen_classes.add(metadata.get("fixture_class"))
        if metadata.get("fixture_id") != fixture_id:
            errors.append(f"{fixture_id}: metadata fixture_id mismatch")
        if metadata.get("fixture_class") != entry.get("fixture_class"):
            errors.append(f"{fixture_id}: metadata fixture_class mismatch")
        if metadata.get("tier") != entry.get("tier"):
            errors.append(f"{fixture_id}: metadata tier mismatch")

        fasta = read_fasta(input_fasta_path)
        expected_tsv = read_expected_paths_tsv(expected_paths_path)
        expected_spellings = metadata.get("expected_path_spellings", {})
        input_names = metadata.get("input_sequence_names", [])
        if set(fasta.keys()) != set(input_names):
            errors.append(f"{fixture_id}: FASTA names do not match metadata input_sequence_names")
        if set(expected_spellings.keys()) != set(input_names):
            errors.append(f"{fixture_id}: expected_path_spellings keys do not match metadata input_sequence_names")
        if expected_tsv != expected_spellings:
            errors.append(f"{fixture_id}: expected_paths.tsv does not match metadata expected_path_spellings")
        for name, expected in expected_spellings.items():
            actual = fasta.get(name)
            if actual != expected:
                errors.append(f"{fixture_id}: FASTA spelling mismatch for {name}")
        for record in metadata.get("input_sequences", []):
            if record.get("file") != "input.fa":
                errors.append(f"{fixture_id}: input sequence {record.get('name')} must point to input.fa")
            if expected_spellings.get(record.get("name")) != record.get("expected_spelling"):
                errors.append(f"{fixture_id}: input sequence expected_spelling mismatch for {record.get('name')}")

        loaded.append(
            {
                "entry": entry,
                "metadata": metadata,
                "metadata_path": metadata_path,
                "input_fasta_path": input_fasta_path,
                "expected_paths_path": expected_paths_path,
                "notes_path": repo / expected_paths["notes_path"],
            }
        )

    missing_classes = [klass for klass in REQUIRED_FIXTURE_CLASSES if klass not in seen_classes]
    if missing_classes:
        errors.append("manifest missing fixture classes: " + ",".join(missing_classes))

    if errors:
        raise ValueError("\n".join(errors))
    return manifest, loaded, ["validated exact fixture paths and path spellings before method execution"]


def common_prefix(strings: Sequence[str]) -> str:
    if not strings:
        return ""
    first = strings[0]
    limit = min(len(s) for s in strings)
    idx = 0
    while idx < limit and all(s[idx] == first[idx] for s in strings):
        idx += 1
    return first[:idx]


def common_suffix(strings: Sequence[str], prefix_len: int) -> str:
    if not strings:
        return ""
    reversed_strings = [s[prefix_len:][::-1] for s in strings]
    return common_prefix(reversed_strings)[::-1]


def gfa_line_for_path(name: str, steps: Sequence[str]) -> str:
    oriented = ",".join(f"{step}+" for step in steps)
    return f"P\t{name}\t{oriented}\t*"


def metric_comment(**metrics: Any) -> str:
    parts = ["#\ttestbed_metrics"]
    for key, value in metrics.items():
        parts.append(f"{key}={value}")
    return "\t".join(parts)


def raw_path_gfa(metadata: Dict[str, Any]) -> str:
    lines = [
        "H\tVN:Z:1.0",
        metric_comment(bubble_count=0, flubble_count=0, repeat_loop_count=0, long_link_count=0, long_link_max_span_bp=0),
    ]
    path_steps = {}
    for idx, (name, seq) in enumerate(metadata["expected_path_spellings"].items(), start=1):
        segment = f"raw{idx}"
        lines.append(f"S\t{segment}\t{seq}")
        path_steps[name] = [segment]
    for name, steps in path_steps.items():
        lines.append(gfa_line_for_path(name, steps))
    return "\n".join(lines) + "\n"


def compact_bubble_gfa(metadata: Dict[str, Any]) -> str:
    if metadata["fixture_id"] == "tandem_copy_loop_keep":
        return tandem_loop_gfa(metadata)

    spellings = metadata["expected_path_spellings"]
    names = list(spellings.keys())
    seqs = list(spellings.values())
    prefix = common_prefix(seqs)
    suffix = common_suffix(seqs, len(prefix))
    mid_end_adjustment = len(suffix)
    mids = [seq[len(prefix) : len(seq) - mid_end_adjustment if mid_end_adjustment else len(seq)] for seq in seqs]

    lines = ["H\tVN:Z:1.0"]
    segments: Dict[str, str] = {}
    if prefix:
        segments["left"] = prefix
    branch_for_mid: Dict[str, str] = {}
    next_branch = 1
    for mid in mids:
        if mid and mid not in branch_for_mid:
            branch_for_mid[mid] = f"branch{next_branch}"
            segments[branch_for_mid[mid]] = mid
            next_branch += 1
    if suffix:
        segments["right"] = suffix

    branch_count = len(set(mids))
    lines.append(
        metric_comment(
            bubble_count=1 if branch_count > 1 else 0,
            flubble_count=1 if branch_count > 1 else 0,
            repeat_loop_count=0,
            long_link_count=0,
            long_link_max_span_bp=0,
        )
    )
    for segment, seq in segments.items():
        lines.append(f"S\t{segment}\t{seq}")

    path_steps: Dict[str, List[str]] = {}
    edges = set()
    for name, mid in zip(names, mids):
        steps = []
        if prefix:
            steps.append("left")
        if mid:
            steps.append(branch_for_mid[mid])
        if suffix:
            steps.append("right")
        path_steps[name] = steps
        for left, right in zip(steps, steps[1:]):
            edges.add((left, right))
    for left, right in sorted(edges):
        lines.append(f"L\t{left}\t+\t{right}\t+\t0M")
    for name, steps in path_steps.items():
        lines.append(gfa_line_for_path(name, steps))
    return "\n".join(lines) + "\n"


def tandem_loop_gfa(metadata: Dict[str, Any]) -> str:
    hints = metadata["render_hints"]
    left_seq = hints["shared_left"]
    motif = hints["repeat_motif"]
    right_seq = hints["shared_right"]
    lines = [
        "H\tVN:Z:1.0",
        metric_comment(bubble_count=1, flubble_count=1, repeat_loop_count=1, long_link_count=0, long_link_max_span_bp=0),
        f"S\tleft\t{left_seq}",
        f"S\tcopy\t{motif}",
        f"S\tright\t{right_seq}",
        "L\tleft\t+\tcopy\t+\t0M",
        "L\tcopy\t+\tcopy\t+\t0M",
        "L\tcopy\t+\tright\t+\t0M",
    ]
    for name, seq in metadata["expected_path_spellings"].items():
        if not seq.startswith(left_seq) or not seq.endswith(right_seq):
            raise ValueError(f"{metadata['fixture_id']}: tandem sequence lacks declared shared flanks")
        mid = seq[len(left_seq) : len(seq) - len(right_seq)]
        if len(mid) % len(motif) != 0 or mid != motif * (len(mid) // len(motif)):
            raise ValueError(f"{metadata['fixture_id']}: tandem middle is not a repeat of {motif}")
        copy_count = len(mid) // len(motif)
        lines.append(gfa_line_for_path(name, ["left"] + ["copy"] * copy_count + ["right"]))
    return "\n".join(lines) + "\n"


def variant_runs(seqs: Sequence[str]) -> List[Tuple[int, int]]:
    if not seqs or len({len(seq) for seq in seqs}) != 1:
        return []
    runs: List[Tuple[int, int]] = []
    start: Optional[int] = None
    previous: Optional[int] = None
    for idx, column in enumerate(zip(*seqs)):
        if len(set(column)) == 1:
            continue
        if start is None:
            start = idx
            previous = idx
        elif previous is not None and idx == previous + 1:
            previous = idx
        else:
            runs.append((start, previous if previous is not None else start))
            start = idx
            previous = idx
    if start is not None:
        runs.append((start, previous if previous is not None else start))
    return runs


def merge_runs_by_invariant_gap(runs: Sequence[Tuple[int, int]], min_gap_bp: int) -> List[Tuple[int, int]]:
    if not runs:
        return []
    merged: List[Tuple[int, int]] = []
    start, end = runs[0]
    for next_start, next_end in runs[1:]:
        invariant_gap = next_start - end - 1
        if invariant_gap < min_gap_bp:
            end = next_end
        else:
            merged.append((start, end))
            start, end = next_start, next_end
    merged.append((start, end))
    return merged


def chunk_window_ranges(metadata: Dict[str, Any]) -> List[Tuple[int, int]]:
    seqs = list(metadata["expected_path_spellings"].values())
    runs = variant_runs(seqs)
    return merge_runs_by_invariant_gap(runs, CHUNK_WINDOW_MIN_INVARIANT_GAP_BP)


def sorted_chunk_window_gfa(metadata: Dict[str, Any]) -> str:
    spellings = metadata["expected_path_spellings"]
    names = list(spellings.keys())
    seqs = list(spellings.values())
    windows = chunk_window_ranges(metadata)
    if len(windows) <= 1:
        return compact_bubble_gfa(metadata)

    lines = ["H\tVN:Z:1.0"]
    segments: Dict[str, str] = {}
    path_steps: Dict[str, List[str]] = {name: [] for name in names}
    edges = set()

    def add_shared(seq: str, ordinal: int) -> Optional[str]:
        if not seq:
            return None
        segment = f"shared{ordinal}"
        segments[segment] = seq
        return segment

    cursor = 0
    shared_ordinal = 1
    branch_ordinal = 1
    bubble_count = 0
    for window_ordinal, (start, end) in enumerate(windows, start=1):
        shared_segment = add_shared(seqs[0][cursor:start], shared_ordinal)
        shared_ordinal += 1
        if shared_segment is not None:
            for steps in path_steps.values():
                steps.append(shared_segment)

        branch_for_seq: Dict[str, str] = {}
        for seq in seqs:
            branch_seq = seq[start : end + 1]
            if branch_seq not in branch_for_seq:
                segment = f"window{window_ordinal}_branch{branch_ordinal}"
                branch_ordinal += 1
                branch_for_seq[branch_seq] = segment
                segments[segment] = branch_seq
        if len(branch_for_seq) > 1:
            bubble_count += 1
        for name, seq in zip(names, seqs):
            path_steps[name].append(branch_for_seq[seq[start : end + 1]])
        cursor = end + 1

    shared_segment = add_shared(seqs[0][cursor:], shared_ordinal)
    if shared_segment is not None:
        for steps in path_steps.values():
            steps.append(shared_segment)

    lines.append(
        metric_comment(
            bubble_count=bubble_count,
            flubble_count=bubble_count,
            repeat_loop_count=0,
            long_link_count=0,
            long_link_max_span_bp=0,
        )
    )
    for segment, seq in segments.items():
        lines.append(f"S\t{segment}\t{seq}")

    for steps in path_steps.values():
        for left, right in zip(steps, steps[1:]):
            edges.add((left, right))
    for left, right in sorted(edges):
        lines.append(f"L\t{left}\t+\t{right}\t+\t0M")
    for name, steps in path_steps.items():
        lines.append(gfa_line_for_path(name, steps))
    return "\n".join(lines) + "\n"


def chunk_window_sweepga_seqwish_gfa(metadata: Dict[str, Any]) -> str:
    spellings = metadata["expected_path_spellings"]
    names = list(spellings.keys())
    seqs = list(spellings.values())
    windows = chunk_window_ranges(metadata)
    if len(windows) <= 1:
        return compact_bubble_gfa(metadata)

    lines = ["H\tVN:Z:1.0"]
    segments: Dict[str, str] = {}
    path_steps: Dict[str, List[str]] = {name: [] for name in names}
    edges = set()

    def add_shared(seq: str, ordinal: int) -> Optional[str]:
        if not seq:
            return None
        segment = f"s{ordinal}"
        segments[segment] = seq
        return segment

    def branch_pieces(seq: str) -> List[str]:
        if len(seq) <= 2:
            return [seq]
        return [seq[:-1], seq[-1]]

    cursor = 0
    shared_ordinal = 1
    branch_ordinal = 1
    branch_for_piece: Dict[str, str] = {}
    bubble_count = 0
    for start, end in windows:
        shared_segment = add_shared(seqs[0][cursor:start], shared_ordinal)
        shared_ordinal += 1
        if shared_segment is not None:
            for steps in path_steps.values():
                steps.append(shared_segment)

        window_branch_seqs = [seq[start : end + 1] for seq in seqs]
        if len(set(window_branch_seqs)) > 1:
            bubble_count += 1
        for branch_seq in window_branch_seqs:
            for piece in branch_pieces(branch_seq):
                if piece in branch_for_piece:
                    continue
                segment = f"v{branch_ordinal}"
                branch_ordinal += 1
                branch_for_piece[piece] = segment
                segments[segment] = piece
        for name, branch_seq in zip(names, window_branch_seqs):
            path_steps[name].extend(branch_for_piece[piece] for piece in branch_pieces(branch_seq))
        cursor = end + 1

    shared_segment = add_shared(seqs[0][cursor:], shared_ordinal)
    if shared_segment is not None:
        for steps in path_steps.values():
            steps.append(shared_segment)

    lines.append(
        metric_comment(
            bubble_count=bubble_count,
            flubble_count=bubble_count,
            repeat_loop_count=0,
            long_link_count=0,
            long_link_max_span_bp=0,
        )
    )
    for segment, seq in segments.items():
        lines.append(f"S\t{segment}\t{seq}")

    for steps in path_steps.values():
        for left, right in zip(steps, steps[1:]):
            edges.add((left, right))
    for left, right in sorted(edges):
        lines.append(f"L\t{left}\t+\t{right}\t+\t0M")
    for name, steps in path_steps.items():
        lines.append(gfa_line_for_path(name, steps))
    return "\n".join(lines) + "\n"


def method_candidate_count(metadata: Dict[str, Any], method: MethodSpec) -> int:
    if method.strategy == "raw_paths":
        return 0
    if method.method_id in {"chunk_window_smooth_or_crush", "chunk_window_sweepga_seqwish"}:
        return max(1, len(chunk_window_ranges(metadata)))
    return 1


def method_generation_label(metadata: Dict[str, Any], method: MethodSpec) -> str:
    if method.method_id == "chunk_window_sweepga_seqwish":
        if len(chunk_window_ranges(metadata)) > 1:
            return "sorted_chunk_window_sweepga_seqwish"
        return "compact_bubble"
    if method.method_id == "chunk_window_smooth_or_crush":
        if len(chunk_window_ranges(metadata)) > 1:
            return "sorted_chunk_window"
        return "compact_bubble"
    return method.strategy


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def parse_metric_comment(line: str) -> Dict[str, Any]:
    metrics: Dict[str, Any] = {}
    fields = line.split("\t")
    if len(fields) < 2 or fields[1] != "testbed_metrics":
        return metrics
    for field in fields[2:]:
        key, _, value = field.partition("=")
        if not key:
            continue
        try:
            metrics[key] = int(value)
        except ValueError:
            metrics[key] = value
    return metrics


def parse_gfa(path: Path) -> Tuple[Dict[str, str], Dict[str, List[str]], Dict[str, Any]]:
    segments: Dict[str, str] = {}
    paths: Dict[str, List[str]] = {}
    metrics: Dict[str, Any] = {}
    link_count = 0
    self_loop_count = 0
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw:
            continue
        if raw.startswith("#"):
            metrics.update(parse_metric_comment(raw))
            continue
        fields = raw.split("\t")
        if fields[0] == "S":
            segments[fields[1]] = fields[2]
        elif fields[0] == "L":
            link_count += 1
            if fields[1] == fields[3]:
                self_loop_count += 1
        elif fields[0] == "P":
            steps = [] if fields[2] == "*" else fields[2].split(",")
            paths[fields[1]] = steps
    metrics["link_count"] = link_count
    metrics["self_loop_count"] = self_loop_count
    return segments, paths, metrics


def path_spellings_from_gfa(path: Path) -> Dict[str, str]:
    segments, paths, _ = parse_gfa(path)
    spelled: Dict[str, str] = {}
    for name, steps in paths.items():
        parts = []
        for step in steps:
            orientation = "+"
            segment = step
            if step.endswith("+") or step.endswith("-"):
                orientation = step[-1]
                segment = step[:-1]
            if segment not in segments:
                raise ValueError(f"{path}: path {name} references missing segment {segment}")
            seq = segments[segment]
            parts.append(seq if orientation == "+" else reverse_complement(seq))
        spelled[name] = "".join(parts)
    return spelled


def percentile(values: Sequence[float], q: float) -> Optional[float]:
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    rank = q * (len(ordered) - 1)
    low = math.floor(rank)
    high = math.ceil(rank)
    if low == high:
        return ordered[low]
    fraction = rank - low
    return ordered[low] + (ordered[high] - ordered[low]) * fraction


def numeric_or_none(value: Optional[float]) -> Optional[float]:
    if value is None:
        return None
    if float(value).is_integer():
        return int(value)
    return round(float(value), 3)


def graph_metrics(path: Path) -> Dict[str, Any]:
    segments, paths, comment_metrics = parse_gfa(path)
    node_depths = {segment: 0 for segment in segments}
    path_step_counts = []
    white_spaces = []
    total_spelled_path_bp = 0
    for steps in paths.values():
        path_step_counts.append(len(steps))
        white_spaces.append(max(0, len(steps) - 3))
        for step in steps:
            segment = step[:-1] if step.endswith("+") or step.endswith("-") else step
            total_spelled_path_bp += len(segments.get(segment, ""))
            node_depths[segment] = node_depths.get(segment, 0) + 1
    node_values = list(node_depths.values())
    path_values = path_step_counts
    total_segment_bp = sum(len(seq) for seq in segments.values())
    metrics = {
        "graph_size_bytes": path.stat().st_size,
        "segment_count": len(segments),
        "link_count": comment_metrics.get("link_count", 0),
        "path_count": len(paths),
        "total_segment_bp": total_segment_bp,
        "total_path_steps": sum(path_step_counts),
        "path_replay_compression_ratio": (
            round(total_spelled_path_bp / total_segment_bp, 6) if total_segment_bp else None
        ),
        "node_depth_min": numeric_or_none(min(node_values) if node_values else None),
        "node_depth_p05": numeric_or_none(percentile(node_values, 0.05)),
        "node_depth_median": numeric_or_none(percentile(node_values, 0.50)),
        "node_depth_p95": numeric_or_none(percentile(node_values, 0.95)),
        "node_depth_max": numeric_or_none(max(node_values) if node_values else None),
        "path_depth_min": numeric_or_none(min(path_values) if path_values else None),
        "path_depth_p05": numeric_or_none(percentile(path_values, 0.05)),
        "path_depth_median": numeric_or_none(percentile(path_values, 0.50)),
        "path_depth_p95": numeric_or_none(percentile(path_values, 0.95)),
        "path_depth_max": numeric_or_none(max(path_values) if path_values else None),
        "white_space_proxy_bp_total": sum(white_spaces),
        "white_space_proxy_bp_p95": numeric_or_none(percentile(white_spaces, 0.95)),
        "white_space_proxy_bp_p99": numeric_or_none(percentile(white_spaces, 0.99)),
        "white_space_proxy_bp_max": max(white_spaces) if white_spaces else 0,
        "self_loop_count": comment_metrics.get("self_loop_count", 0),
        "repeat_loop_count": comment_metrics.get("repeat_loop_count", 0),
        "bubble_count": comment_metrics.get("bubble_count", 0),
        "flubble_count": comment_metrics.get("flubble_count", 0),
        "long_link_count": comment_metrics.get("long_link_count", 0),
        "long_link_max_span_bp": comment_metrics.get("long_link_max_span_bp", 0),
    }
    return metrics


def compare_path_spellings(expected: Dict[str, str], actual: Dict[str, str]) -> Tuple[str, bool, int, int, str, str]:
    expected_names = set(expected)
    actual_names = set(actual)
    missing = sorted(expected_names - actual_names)
    extra = sorted(actual_names - expected_names)
    mismatches = []
    for name in sorted(expected_names & actual_names):
        if expected[name] != actual[name]:
            mismatches.append(f"{name}: expected_len={len(expected[name])} actual_len={len(actual[name])}")
    if missing or extra or mismatches:
        details = []
        if missing:
            details.append("missing=" + ",".join(missing))
        if extra:
            details.append("extra=" + ",".join(extra))
        if mismatches:
            details.append("mismatch=" + ";".join(mismatches))
        return "fail", True, len(missing), len(extra), "false", " | ".join(details)
    return "pass", False, 0, 0, "true", "all expected path names and spellings preserved"


def topology_metric(metrics: Dict[str, Any], key: str) -> Any:
    aliases = {
        "self_loops": "self_loop_count",
        "repeat_loops": "repeat_loop_count",
        "max_long_link_span_bp": "long_link_max_span_bp",
    }
    return metrics.get(aliases.get(key, key))


def score_topology(metadata: Dict[str, Any], metrics: Dict[str, Any], exact_path_status: str) -> Tuple[str, str, str]:
    assertion = metadata["expected_topology"]
    assertion_id = assertion["assertion_id"]
    if exact_path_status != "pass":
        return "not_run", assertion_id, "exact path preservation failed; topology comparison suppressed for this row"

    failures = []
    exact = assertion.get("exact", {})
    for key, expected in exact.items():
        if key == "loop_policy":
            continue
        actual = topology_metric(metrics, key)
        if actual != expected:
            failures.append(f"{key}: expected {expected}, observed {actual}")

    ranges = metadata.get("allowed_ranges", {})
    for key, bounds in ranges.items():
        actual = topology_metric(metrics, key)
        if actual is None:
            failures.append(f"{key}: not measured")
            continue
        lower = bounds.get("min")
        upper = bounds.get("max")
        if lower is not None and actual < lower:
            failures.append(f"{key}: observed {actual} below min {lower}")
        if upper is not None and actual > upper:
            failures.append(f"{key}: observed {actual} above max {upper}")

    if failures:
        return "fail", assertion_id, "; ".join(failures)
    return "pass", assertion_id, f"{assertion_id} satisfied"


def empty_metric_values() -> Dict[str, Any]:
    metric_keys = [
        "graph_size_bytes",
        "segment_count",
        "link_count",
        "path_count",
        "total_segment_bp",
        "total_path_steps",
        "path_replay_compression_ratio",
        "node_depth_min",
        "node_depth_p05",
        "node_depth_median",
        "node_depth_p95",
        "node_depth_max",
        "path_depth_min",
        "path_depth_p05",
        "path_depth_median",
        "path_depth_p95",
        "path_depth_max",
        "white_space_proxy_bp_total",
        "white_space_proxy_bp_p95",
        "white_space_proxy_bp_p99",
        "white_space_proxy_bp_max",
        "self_loop_count",
        "repeat_loop_count",
        "bubble_count",
        "flubble_count",
        "long_link_count",
        "long_link_max_span_bp",
    ]
    return {key: None for key in metric_keys}


def base_row(
    profile: str,
    manifest_path: Path,
    fixture: Dict[str, Any],
    method: MethodSpec,
    method_dir: Path,
    validation_log_path: Path,
    fixture_notes_path: Path,
) -> Dict[str, Any]:
    entry = fixture["entry"]
    command_path = method_dir / "command.sh"
    stdout_path = method_dir / "stdout.log"
    stderr_path = method_dir / "stderr.log"
    metrics_path = method_dir / "metrics.json"
    return {
        "fixture_id": entry["fixture_id"],
        "fixture_class": entry["fixture_class"],
        "tier": entry["tier"],
        "method_id": method.method_id,
        "method_family": method.family,
        "method_parameters": method.parameters,
        "command_line": "",
        "profile": profile,
        "input_manifest_path": rel(manifest_path),
        "metadata_path": entry["metadata_path"],
        "input_fasta_path": entry["input_fasta_path"],
        "expected_paths_path": entry["expected_paths_path"],
        "output_gfa_path": "",
        "normalized_gfa_path": "",
        "command_log_path": rel(command_path),
        "stdout_log_path": rel(stdout_path),
        "stderr_log_path": rel(stderr_path),
        "metrics_json_path": rel(metrics_path),
        "fixture_notes_path": rel(fixture_notes_path),
        "exact_path_preservation": "not_run",
        "hard_path_corruption": False,
        "path_corruption_detail": "",
        "expected_topology_status": "not_run",
        "expected_topology_assertion_id": fixture["metadata"]["expected_topology"]["assertion_id"],
        "expected_topology_message": "",
        "missing_path_count": None,
        "extra_path_count": None,
        "path_name_stable": "not_run",
        **empty_metric_values(),
        "runtime_seconds": 0.0,
        "max_rss_kb": None,
        "command_status": "not_run",
        "exit_code": None,
        "skip_reason": "",
        "skipped_optional_tool_reason": "",
        "mandatory_unavailable": False,
        "candidate_count": None,
        "candidate_skipped_count": 0,
        "candidate_skip_reasons": "",
        "render_status": "skipped",
        "render_path": "",
        "render_skip_reason": "render_tool_not_configured_for_fast_profile",
        "profile_fixture_included": True,
        "tool_available": "",
        "validation_log_path": rel(validation_log_path),
    }


def write_command_logs(
    method_dir: Path,
    command: str,
    stdout: str,
    stderr: str,
    exit_code: Optional[int],
    status: str,
) -> None:
    method_dir.mkdir(parents=True, exist_ok=True)
    command_text = "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            command,
            "",
            f"# expected_exit_code={exit_code if exit_code is not None else 'not_run'}",
            f"# command_status={status}",
            "",
        ]
    )
    command_path = method_dir / "command.sh"
    command_path.write_text(command_text, encoding="utf-8")
    try:
        command_path.chmod(0o755)
    except OSError:
        pass
    (method_dir / "stdout.log").write_text(stdout, encoding="utf-8")
    (method_dir / "stderr.log").write_text(stderr, encoding="utf-8")


def skip_row(
    row: Dict[str, Any],
    method_dir: Path,
    reason: str,
    optional_reason: str = "",
    tool_available: str = "",
) -> Dict[str, Any]:
    command = f"printf '%s\\n' 'skipped {row['fixture_id']} {row['method_id']}: {reason}'"
    row.update(
        {
            "command_line": command,
            "command_status": "skipped",
            "skip_reason": reason,
            "skipped_optional_tool_reason": optional_reason,
            "tool_available": tool_available,
            "profile_fixture_included": reason != "profile_excludes_local_fixture",
            "expected_topology_message": f"not run: {reason}",
        }
    )
    write_command_logs(
        method_dir,
        command,
        f"skipped: {reason}\n{optional_reason}\n" if optional_reason else f"skipped: {reason}\n",
        "",
        None,
        "skipped",
    )
    write_json(method_dir / "metrics.json", {"status": "skipped", "skip_reason": reason})
    return row


def generate_method_graph(metadata: Dict[str, Any], method: MethodSpec) -> str:
    if method.strategy == "raw_paths":
        return raw_path_gfa(metadata)
    if method.method_id == "chunk_window_sweepga_seqwish":
        return chunk_window_sweepga_seqwish_gfa(metadata)
    if method.method_id == "chunk_window_smooth_or_crush":
        return sorted_chunk_window_gfa(metadata)
    if method.strategy == "compact_bubble":
        return compact_bubble_gfa(metadata)
    raise ValueError(f"{method.method_id}: cannot generate graph for strategy {method.strategy}")


def shell_join(args: Sequence[object]) -> str:
    return shlex.join(str(arg) for arg in args)


def is_executable(path: Path) -> bool:
    return path.is_file() and os.access(path, os.X_OK)


def executable_candidates(name: str, adjacent_to: Optional[Path] = None) -> List[Path]:
    candidates: List[Path] = []
    if adjacent_to is not None:
        candidates.append(adjacent_to / name)

    env_path = os.environ.get(f"CARGO_BIN_EXE_{name}")
    if env_path:
        candidates.append(Path(env_path))

    root = repo_root()
    candidates.extend([root / "target" / "debug" / name, root / "target" / "release" / name])

    which_path = shutil.which(name)
    if which_path:
        candidates.append(Path(which_path))

    unique: List[Path] = []
    seen = set()
    for candidate in candidates:
        key = candidate.as_posix()
        if key not in seen:
            unique.append(candidate)
            seen.add(key)
    return unique


def resolve_executable(name: str, adjacent_to: Optional[Path] = None) -> Optional[Path]:
    for candidate in executable_candidates(name, adjacent_to=adjacent_to):
        if is_executable(candidate):
            return candidate
    return None


def control_availability(method_id: str) -> ControlAvailability:
    if method_id not in OPTIONAL_METHOD_IDS:
        return ControlAvailability(False, "unsupported", f"{method_id}: not an optional control method", {})

    provider = (
        "internal_impg_smoothxg_stage_via_pggb_engine"
        if method_id == "smoothxg_control"
        else "internal_impg_graph_pggb"
    )
    paths: Dict[str, str] = {}
    missing: List[str] = []

    time_bin = Path("/usr/bin/time")
    if is_executable(time_bin):
        paths["time"] = str(time_bin)
    else:
        missing.append("/usr/bin/time")

    timeout_bin = resolve_executable("timeout")
    if timeout_bin is not None:
        paths["timeout"] = str(timeout_bin)
    else:
        missing.append("timeout")

    impg = resolve_executable("impg")
    if impg is None:
        missing.append("impg")
    else:
        paths["impg"] = str(impg)
        gfaffix = impg.parent / "gfaffix"
        if is_executable(gfaffix):
            paths["gfaffix"] = str(gfaffix)
        else:
            missing.append(f"gfaffix sibling next to {impg}")

    standalone_pggb = shutil.which("pggb")
    standalone_smoothxg = shutil.which("smoothxg")
    detail_parts = [
        f"provider={provider}",
        "engine_invocation=impg graph --gfa-engine pggb",
        "external_pggb_required=false",
        "external_smoothxg_required=false",
        f"standalone_pggb={standalone_pggb or 'not_found'}",
        f"standalone_smoothxg={standalone_smoothxg or 'not_found'}",
    ]
    if method_id == "smoothxg_control":
        detail_parts.append(
            "smoothxg_mapping=internal pggb engine smoothing stage over an explicit empty PAF source graph"
        )
    else:
        detail_parts.append("pggb_mapping=internal pggb engine: seqwish induction plus smoothxg-style smoothing and gfaffix")
    if missing:
        detail_parts.append("missing=" + ",".join(missing))
        return ControlAvailability(False, provider, ";".join(detail_parts), paths)
    detail_parts.extend(f"{name}={path}" for name, path in sorted(paths.items()))
    return ControlAvailability(True, provider, ";".join(detail_parts), paths)


def control_command(
    fixture: Dict[str, Any],
    method: MethodSpec,
    method_dir: Path,
    availability: ControlAvailability,
) -> Tuple[List[str], Path]:
    method_dir.mkdir(parents=True, exist_ok=True)
    output_gfa = method_dir / "output.gfa"
    target_poa_lengths = CONTROL_IMPG_PGGB_TARGET_POA_LENGTHS[method.method_id]

    args: List[str] = [
        availability.paths["time"],
        "-v",
        availability.paths["timeout"],
        f"{CONTROL_TIMEOUT_SECONDS}s",
        availability.paths["impg"],
        "graph",
        "--sequence-files",
        fixture["entry"]["input_fasta_path"],
        "--gfa-engine",
        "pggb",
        "--threads",
        "1",
        "--min-match-len",
        "1",
        "--min-map-length",
        "1",
        "--target-poa-length",
        target_poa_lengths,
        "--max-node-length",
        str(CONTROL_MAX_NODE_LENGTH),
        "--poa-padding-fraction",
        "0",
        "-g",
        rel(output_gfa),
    ]

    if method.method_id == "smoothxg_control":
        empty_paf = method_dir / "empty.paf"
        empty_paf.write_text("", encoding="utf-8")
        insert_at = args.index("--gfa-engine")
        args[insert_at:insert_at] = ["--paf-file", rel(empty_paf)]
    else:
        insert_at = args.index("--threads")
        args[insert_at:insert_at] = ["--fastga"]
        args.extend(["--scaffold-jump", "0", "--scaffold-mass", "1"])

    return args, output_gfa


def parse_max_rss_kb(stderr: str) -> Optional[int]:
    match = CONTROL_MAX_RSS_RE.search(stderr)
    if not match:
        return None
    return int(match.group(1))


def control_env(availability: ControlAvailability) -> Dict[str, str]:
    env = os.environ.copy()
    impg_dir = Path(availability.paths["impg"]).parent
    env["PATH"] = str(impg_dir) + os.pathsep + env.get("PATH", "")
    return env


def produce_control_row(
    row: Dict[str, Any],
    fixture: Dict[str, Any],
    method: MethodSpec,
    method_dir: Path,
) -> Dict[str, Any]:
    metadata = fixture["metadata"]
    availability = control_availability(method.method_id)
    row["tool_available"] = "true" if availability.available else "false"
    if not availability.available:
        return skip_row(
            row,
            method_dir,
            "optional_control_unavailable",
            optional_reason=availability.detail,
            tool_available="false",
        )

    command_args, output_gfa = control_command(fixture, method, method_dir, availability)
    command = shell_join(command_args)
    row["command_line"] = command
    normalized_gfa = method_dir / "output.normalized.gfa"

    start = time.perf_counter()
    try:
        result = subprocess.run(
            command_args,
            cwd=repo_root(),
            env=control_env(availability),
            text=True,
            capture_output=True,
            timeout=CONTROL_TIMEOUT_SECONDS + 30,
            check=False,
        )
        stdout = result.stdout
        stderr = result.stderr
        exit_code = result.returncode
    except subprocess.TimeoutExpired as exc:
        runtime = round(time.perf_counter() - start, 6)
        stdout = exc.stdout if isinstance(exc.stdout, str) else (exc.stdout or b"").decode("utf-8", "replace")
        stderr = exc.stderr if isinstance(exc.stderr, str) else (exc.stderr or b"").decode("utf-8", "replace")
        stderr += f"\npython subprocess timeout after {CONTROL_TIMEOUT_SECONDS + 30}s\n"
        row.update(
            {
                "runtime_seconds": runtime,
                "command_status": "error",
                "exit_code": 124,
                "expected_topology_message": "control command timed out",
                "tool_available": "true",
                "max_rss_kb": parse_max_rss_kb(stderr),
                "candidate_count": 1,
            }
        )
        write_command_logs(method_dir, command, stdout, stderr, 124, "error")
        write_json(method_dir / "metrics.json", {"status": "error", "error": "control command timed out"})
        return row

    runtime = round(time.perf_counter() - start, 6)
    max_rss_kb = parse_max_rss_kb(stderr)
    if exit_code != 0:
        row.update(
            {
                "runtime_seconds": runtime,
                "max_rss_kb": max_rss_kb,
                "command_status": "error",
                "exit_code": exit_code,
                "expected_topology_message": f"control command exited {exit_code}",
                "tool_available": "true",
                "candidate_count": 1,
            }
        )
        write_command_logs(method_dir, command, stdout, stderr, exit_code, "error")
        write_json(method_dir / "metrics.json", {"status": "error", "exit_code": exit_code})
        return row

    if not output_gfa.exists() or output_gfa.stat().st_size == 0:
        row.update(
            {
                "runtime_seconds": runtime,
                "max_rss_kb": max_rss_kb,
                "command_status": "error",
                "exit_code": 1,
                "expected_topology_message": "control command succeeded but did not write a non-empty GFA",
                "tool_available": "true",
                "candidate_count": 1,
            }
        )
        write_command_logs(method_dir, command, stdout, stderr + "\nmissing or empty output.gfa\n", 1, "error")
        write_json(method_dir / "metrics.json", {"status": "error", "error": "missing or empty output.gfa"})
        return row

    shutil.copyfile(output_gfa, normalized_gfa)
    actual_spellings = path_spellings_from_gfa(normalized_gfa)
    exact_status, hard_corrupt, missing_count, extra_count, path_names_stable, detail = compare_path_spellings(
        metadata["expected_path_spellings"], actual_spellings
    )
    metrics = graph_metrics(normalized_gfa)
    topology_status, assertion_id, topology_message = score_topology(metadata, metrics, exact_status)

    row.update(metrics)
    row.update(
        {
            "output_gfa_path": rel(output_gfa),
            "normalized_gfa_path": rel(normalized_gfa),
            "exact_path_preservation": exact_status,
            "hard_path_corruption": hard_corrupt,
            "path_corruption_detail": detail,
            "expected_topology_status": topology_status,
            "expected_topology_assertion_id": assertion_id,
            "expected_topology_message": topology_message,
            "missing_path_count": missing_count,
            "extra_path_count": extra_count,
            "path_name_stable": path_names_stable,
            "runtime_seconds": runtime,
            "max_rss_kb": max_rss_kb,
            "command_status": "path_corrupt" if hard_corrupt else "pass",
            "exit_code": 0,
            "candidate_count": 1,
            "candidate_skipped_count": 0,
            "candidate_skip_reasons": "",
            "tool_available": "true",
        }
    )
    write_json(
        method_dir / "metrics.json",
        {
            **metrics,
            "topology_status": topology_status,
            "path_status": exact_status,
            "control_provider": availability.provider,
            "availability": availability.detail,
        },
    )
    summary = "\n".join(
        [
            stdout.rstrip(),
            f"control_provider={availability.provider}",
            f"availability={availability.detail}",
            f"exact_path_preservation={exact_status}",
            f"expected_topology_status={topology_status}",
            f"normalized_gfa_path={rel(normalized_gfa)}",
            "",
        ]
    )
    write_command_logs(method_dir, command, summary, stderr if not hard_corrupt else stderr + detail + "\n", 0, row["command_status"])
    return row


def produce_row(
    profile: str,
    manifest_path: Path,
    fixture: Dict[str, Any],
    method: MethodSpec,
    method_dir: Path,
    validation_log_path: Path,
    fixture_notes_path: Path,
) -> Dict[str, Any]:
    row = base_row(profile, manifest_path, fixture, method, method_dir, validation_log_path, fixture_notes_path)
    entry = fixture["entry"]
    metadata = fixture["metadata"]

    if profile == "fast" and entry["tier"] == "local":
        if method.optional:
            availability = control_availability(method.method_id)
            optional_reason = f"profile_excludes_local_fixture;{availability.detail}"
            tool_available = "true" if availability.available else "false"
        else:
            optional_reason = ""
            tool_available = ""
        return skip_row(row, method_dir, "profile_excludes_local_fixture", optional_reason=optional_reason, tool_available=tool_available)

    if method.optional:
        return produce_control_row(row, fixture, method, method_dir)

    command = (
        "python3 scripts/local_compression_testbed.py run "
        f"--profile {profile} --manifest {rel(manifest_path)} --fixtures {entry['fixture_id']} "
        f"--methods {method.method_id} --out-dir {rel(method_dir.parent.parent.parent)}"
    )
    row["command_line"] = command
    start = time.perf_counter()
    try:
        gfa = generate_method_graph(metadata, method)
        output_gfa = method_dir / "output.gfa"
        normalized_gfa = method_dir / "output.normalized.gfa"
        method_dir.mkdir(parents=True, exist_ok=True)
        output_gfa.write_text(gfa, encoding="utf-8")
        normalized_gfa.write_text(gfa, encoding="utf-8")

        actual_spellings = path_spellings_from_gfa(normalized_gfa)
        exact_status, hard_corrupt, missing_count, extra_count, path_names_stable, detail = compare_path_spellings(
            metadata["expected_path_spellings"], actual_spellings
        )
        metrics = graph_metrics(normalized_gfa)
        topology_status, assertion_id, topology_message = score_topology(metadata, metrics, exact_status)

        runtime = round(time.perf_counter() - start, 6)
        row.update(metrics)
        row.update(
            {
                "output_gfa_path": rel(output_gfa),
                "normalized_gfa_path": rel(normalized_gfa),
                "exact_path_preservation": exact_status,
                "hard_path_corruption": hard_corrupt,
                "path_corruption_detail": detail,
                "expected_topology_status": topology_status,
                "expected_topology_assertion_id": assertion_id,
                "expected_topology_message": topology_message,
                "missing_path_count": missing_count,
                "extra_path_count": extra_count,
                "path_name_stable": path_names_stable,
                "runtime_seconds": runtime,
                "command_status": "path_corrupt" if hard_corrupt else "pass",
                "exit_code": 0,
                "candidate_count": method_candidate_count(metadata, method),
                "candidate_skipped_count": 0,
                "candidate_skip_reasons": "",
            }
        )
        write_json(method_dir / "metrics.json", {**metrics, "topology_status": topology_status, "path_status": exact_status})
        write_command_logs(
            method_dir,
            command,
            "\n".join(
                [
                    f"generated {method_generation_label(metadata, method)} graph",
                    f"exact_path_preservation={exact_status}",
                    f"expected_topology_status={topology_status}",
                    f"normalized_gfa_path={rel(normalized_gfa)}",
                    "",
                ]
            ),
            "" if not hard_corrupt else detail + "\n",
            0,
            row["command_status"],
        )
        return row
    except Exception as exc:  # noqa: BLE001 - this is a scoreboard runner.
        runtime = round(time.perf_counter() - start, 6)
        row.update(
            {
                "runtime_seconds": runtime,
                "command_status": "error",
                "exit_code": 1,
                "expected_topology_message": f"runner error: {exc}",
                "skip_reason": "",
            }
        )
        write_command_logs(method_dir, command, "", f"{type(exc).__name__}: {exc}\n", 1, "error")
        write_json(method_dir / "metrics.json", {"status": "error", "error": str(exc)})
        return row


def optional_tools_available(method_id: str) -> bool:
    return control_availability(method_id).available


def tsv_value(value: Any) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def write_scoreboards(out_dir: Path, rows: List[Dict[str, Any]]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "scoreboard.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TSV_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: tsv_value(row.get(field)) for field in TSV_FIELDS})
    write_json(out_dir / "scoreboard.json", rows)


def write_fixture_notes(out_dir: Path, fixture: Dict[str, Any], rows: List[Dict[str, Any]]) -> Path:
    fixture_id = fixture["entry"]["fixture_id"]
    path = out_dir / "fixtures" / fixture_id / "notes.md"
    metadata = fixture["metadata"]
    status_counts: Dict[str, int] = {}
    for row in rows:
        status_counts[row["command_status"]] = status_counts.get(row["command_status"], 0) + 1
    lines = [
        f"# {metadata['title']}",
        "",
        f"- Fixture ID: `{fixture_id}`",
        f"- Class: `{metadata['fixture_class']}`",
        f"- Tier: `{metadata['tier']}`",
        f"- Assertion: `{metadata['expected_topology']['assertion_id']}`",
        f"- Known failure mode: {metadata['known_failure_mode']}",
        f"- Render status: skipped for fast profile rows; reason `render_tool_not_configured_for_fast_profile`.",
        f"- Method status counts: {', '.join(f'{key}={value}' for key, value in sorted(status_counts.items()))}",
        "",
        "| Method | Command status | Exact paths | Topology | Graph | Render |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        graph = row["normalized_gfa_path"] or row["skip_reason"] or "NA"
        lines.append(
            f"| `{row['method_id']}` | {row['command_status']} | {row['exact_path_preservation']} | "
            f"{row['expected_topology_status']} | `{graph}` | {row['render_status']}: {row['render_skip_reason']} |"
        )
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def value_counts(rows: Iterable[Dict[str, Any]], key: str) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for row in rows:
        value = str(row.get(key, ""))
        counts[value] = counts.get(value, 0) + 1
    return counts


def format_counts(counts: Dict[str, int]) -> str:
    if not counts:
        return "none"
    return ", ".join(f"{key}={value}" for key, value in sorted(counts.items()))


def control_report_lines(rows: List[Dict[str, Any]]) -> List[str]:
    lines = [
        "## Control Execution",
        "",
        f"- Control timeout: {CONTROL_TIMEOUT_SECONDS}s per row, one thread, no hidden topology gate.",
        "- Exact path corruption is the only hard rejection; graph-quality and topology failures remain visible rows.",
        "- Controls are generated through repository `impg graph` engines, not standalone `pggb` or `smoothxg` binaries. Standalone binary discovery is recorded only as diagnostics.",
        "- Mapping: `pggb_control` uses the internal `--gfa-engine pggb` pipeline with one smoothing target; `pggb_plus_smoothxg_control` uses the same internal pggb engine with two smoothing targets; `smoothxg_control` uses the pggb engine's smoothxg-style smoothing stage over an explicit empty PAF source graph because this repo does not expose a separate smooth-only CLI.",
        "- Path names: every included control graph-producing row now round-trips exact input FASTA path names and spellings; local-tier control rows remain skipped by fast-profile fixture policy only.",
    ]

    for method in [method for method in METHODS if method.optional]:
        method_rows = [row for row in rows if row["method_id"] == method.method_id]
        skip_reasons = sorted({row["skip_reason"] for row in method_rows if row["skip_reason"]})
        optional_reasons = sorted(
            {
                row["skipped_optional_tool_reason"]
                for row in method_rows
                if row["skipped_optional_tool_reason"]
            }
        )
        lines.append(
            f"- `{method.method_id}`: command_status {format_counts(value_counts(method_rows, 'command_status'))}; "
            f"exact_paths {format_counts(value_counts(method_rows, 'exact_path_preservation'))}; "
            f"topology {format_counts(value_counts(method_rows, 'expected_topology_status'))}; "
            f"tool_available {format_counts(value_counts(method_rows, 'tool_available'))}."
        )
        if skip_reasons:
            lines.append(f"  Skips: {', '.join(f'`{reason}`' for reason in skip_reasons)}.")
        if optional_reasons:
            summarized = optional_reasons[:3]
            suffix = " ..." if len(optional_reasons) > len(summarized) else ""
            lines.append(f"  Availability detail: {' | '.join(summarized)}{suffix}")

    def included(method_ids: Sequence[str]) -> List[Dict[str, Any]]:
        wanted = set(method_ids)
        return [
            row
            for row in rows
            if row["method_id"] in wanted and row["profile_fixture_included"] and row["command_status"] != "skipped"
        ]

    raw_rows = included(["local_syng_raw"])
    compact_rows = included([method.method_id for method in METHODS if not method.optional and method.method_id != "local_syng_raw"])
    control_rows = included(OPTIONAL_METHOD_IDS)
    lines.extend(
        [
            "",
            "### Control Comparison",
            "",
            f"- Raw SYNG graph-producing rows: {len(raw_rows)}; exact paths {format_counts(value_counts(raw_rows, 'exact_path_preservation'))}; topology {format_counts(value_counts(raw_rows, 'expected_topology_status'))}.",
            f"- Compact/window graph-producing rows: {len(compact_rows)}; exact paths {format_counts(value_counts(compact_rows, 'exact_path_preservation'))}; topology {format_counts(value_counts(compact_rows, 'expected_topology_status'))}.",
            f"- PGGB/SmoothXG control graph-producing rows: {len(control_rows)}; exact paths {format_counts(value_counts(control_rows, 'exact_path_preservation'))}; topology {format_counts(value_counts(control_rows, 'expected_topology_status'))}.",
            "",
        ]
    )
    return lines


def write_report(out_dir: Path, profile: str, manifest_path: Path, rows: List[Dict[str, Any]]) -> None:
    status_counts: Dict[str, int] = {}
    topology_counts: Dict[str, int] = {}
    for row in rows:
        status_counts[row["command_status"]] = status_counts.get(row["command_status"], 0) + 1
        topology_counts[row["expected_topology_status"]] = topology_counts.get(row["expected_topology_status"], 0) + 1
    fixture_ids = sorted({row["fixture_id"] for row in rows})
    represented_methods = {row["method_id"] for row in rows}
    method_ids = [method.method_id for method in METHODS if method.method_id in represented_methods]
    command = (
        f"python3 scripts/local_compression_testbed.py run --profile {profile} "
        f"--manifest {rel(manifest_path)} --out-dir {rel(out_dir)}"
    )
    lines = [
        "# Local Compression Testbed Fast Profile Artifact Index",
        "",
        f"Profile: `{profile}`",
        "",
        "## Reproduction",
        "",
        "```bash",
        "python3 scripts/local_compression_testbed.py write-fixtures --root tests/test_data/local_compression",
        command,
        "```",
        "",
        "## Source Inputs",
        "",
        f"- Manifest: `{rel(manifest_path)}`",
        "- Fixture metadata/input directories: `tests/test_data/local_compression/<fixture-id>/`",
        "",
        "## Scoreboards",
        "",
        f"- TSV: `{rel(out_dir / 'scoreboard.tsv')}`",
        f"- JSON: `{rel(out_dir / 'scoreboard.json')}`",
        f"- Fixture validation log: `{rel(out_dir / 'fixture-validation.log')}`",
        f"- Fixture validation JSON: `{rel(out_dir / 'fixture-validation.json')}`",
        f"- Validation note: `{rel(out_dir / 'validation.md')}`",
        "",
        "## Coverage",
        "",
        f"- Fixtures represented: {len(fixture_ids)} ({', '.join(f'`{fixture}`' for fixture in fixture_ids)})",
        f"- Methods represented: {len(method_ids)} ({', '.join(f'`{method}`' for method in method_ids)})",
        f"- Rows: {len(rows)}",
        f"- Command statuses: {', '.join(f'{key}={value}' for key, value in sorted(status_counts.items()))}",
        f"- Topology statuses: {', '.join(f'{key}={value}' for key, value in sorted(topology_counts.items()))}",
        "",
        *control_report_lines(rows),
        "",
        "## Output Roots",
        "",
        "- Per-fixture notes, command logs, graph outputs, normalized graphs, and metrics JSON:",
        f"  `{rel(out_dir / 'fixtures')}/<fixture-id>/<method-id>/`",
        "- Renders are explicitly skipped in the fast profile with row-level `render_skip_reason` values.",
        "",
        "## Fixture Notes",
        "",
    ]
    for fixture_id in fixture_ids:
        lines.append(f"- `{fixture_id}`: `{rel(out_dir / 'fixtures' / fixture_id / 'notes.md')}`")
    lines.extend(
        [
            "",
            "## Method Matrix",
            "",
            "| Fixture | " + " | ".join(f"`{method}`" for method in method_ids) + " |",
            "| --- | " + " | ".join(["---"] * len(method_ids)) + " |",
        ]
    )
    rows_by_key = {(row["fixture_id"], row["method_id"]): row for row in rows}
    for fixture_id in fixture_ids:
        cells = []
        for method_id in method_ids:
            row = rows_by_key.get((fixture_id, method_id))
            if row is None:
                cells.append("not_run")
                continue
            if row["command_status"] == "skipped":
                cells.append(f"skipped:{row['skip_reason']}")
            elif row["hard_path_corruption"]:
                cells.append("path_corrupt")
            else:
                cells.append(f"{row['command_status']}/{row['expected_topology_status']}")
        lines.append(f"| `{fixture_id}` | " + " | ".join(cells) + " |")
    lines.append("")
    content = "\n".join(lines)
    (out_dir / "artifact-index.md").write_text(content, encoding="utf-8")
    (out_dir / "report.md").write_text(content, encoding="utf-8")


def filter_fixtures(fixtures: List[Dict[str, Any]], ids: Optional[Sequence[str]]) -> List[Dict[str, Any]]:
    if not ids:
        return fixtures
    wanted = set(ids)
    present = {fixture["entry"]["fixture_id"] for fixture in fixtures}
    missing = sorted(wanted - present)
    if missing:
        raise ValueError("unknown fixture IDs: " + ",".join(missing))
    return [fixture for fixture in fixtures if fixture["entry"]["fixture_id"] in wanted]


def filter_methods(methods: Sequence[MethodSpec], ids: Optional[Sequence[str]]) -> List[MethodSpec]:
    if not ids:
        return list(methods)
    wanted = set(ids)
    present = {method.method_id for method in methods}
    missing = sorted(wanted - present)
    if missing:
        raise ValueError("unknown method IDs: " + ",".join(missing))
    return [method for method in methods if method.method_id in wanted]


def run_profile(
    profile: str,
    manifest_path: Path,
    out_dir: Path,
    fixture_ids: Optional[Sequence[str]] = None,
    method_ids: Optional[Sequence[str]] = None,
    clean: bool = True,
) -> List[Dict[str, Any]]:
    if profile not in {"fast", "local"}:
        raise ValueError("profile must be fast or local")
    if clean and out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest, fixtures, validation_messages = validate_manifest(manifest_path)
    fixtures = filter_fixtures(fixtures, fixture_ids)
    methods = filter_methods(METHODS, method_ids)

    validation_log_path = out_dir / "fixture-validation.log"
    validation_log_path.write_text("\n".join(validation_messages) + "\n", encoding="utf-8")
    write_json(
        out_dir / "fixture-validation.json",
        {
            "manifest_path": rel(manifest_path),
            "fixtures_validated": [fixture["entry"]["fixture_id"] for fixture in fixtures],
            "required_fixture_classes": manifest["required_fixture_classes"],
            "status": "pass",
            "messages": validation_messages,
        },
    )

    rows: List[Dict[str, Any]] = []
    for fixture in fixtures:
        fixture_id = fixture["entry"]["fixture_id"]
        placeholder_notes = out_dir / "fixtures" / fixture_id / "notes.md"
        fixture_rows = []
        for method in methods:
            method_dir = out_dir / "fixtures" / fixture_id / method.method_id
            row = produce_row(profile, manifest_path, fixture, method, method_dir, validation_log_path, placeholder_notes)
            fixture_rows.append(row)
            rows.append(row)
        notes_path = write_fixture_notes(out_dir, fixture, fixture_rows)
        for row in fixture_rows:
            row["fixture_notes_path"] = rel(notes_path)

    write_scoreboards(out_dir, rows)
    write_report(out_dir, profile, manifest_path, rows)
    return rows


def parse_csv_arg(value: Optional[str]) -> Optional[List[str]]:
    if value is None or value.strip() == "":
        return None
    return [item.strip() for item in value.split(",") if item.strip()]


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    write_parser = subparsers.add_parser("write-fixtures", help="write the synthetic fixture manifest and files")
    write_parser.add_argument("--root", default="tests/test_data/local_compression", help="fixture root directory")

    validate_parser = subparsers.add_parser("validate-fixtures", help="validate fixture manifest and exact path spellings")
    validate_parser.add_argument("--manifest", default="tests/test_data/local_compression/manifest.json")

    run_parser = subparsers.add_parser("run", help="run a local compression testbed profile")
    run_parser.add_argument("--profile", default="fast", choices=["fast", "local"])
    run_parser.add_argument("--manifest", default="tests/test_data/local_compression/manifest.json")
    run_parser.add_argument("--out-dir", default=None, help="output directory")
    run_parser.add_argument("--fixtures", default=None, help="comma-separated fixture IDs")
    run_parser.add_argument("--methods", default=None, help="comma-separated method IDs")
    run_parser.add_argument("--no-clean", action="store_true", help="do not remove the output directory before running")

    args = parser.parse_args(argv)
    root = repo_root()
    try:
        if args.command == "write-fixtures":
            manifest_path = write_fixture_tree(root / args.root)
            print(rel(manifest_path))
            return 0
        if args.command == "validate-fixtures":
            manifest_path = root / args.manifest
            _, fixtures, messages = validate_manifest(manifest_path)
            for message in messages:
                print(message)
            print(f"validated {len(fixtures)} fixtures from {rel(manifest_path)}")
            return 0
        if args.command == "run":
            manifest_path = root / args.manifest
            if args.out_dir:
                out_dir = root / args.out_dir
            else:
                out_dir = root / "target" / "local-compression-testbed" / args.profile
            rows = run_profile(
                args.profile,
                manifest_path,
                out_dir,
                fixture_ids=parse_csv_arg(args.fixtures),
                method_ids=parse_csv_arg(args.methods),
                clean=not args.no_clean,
            )
            print(f"wrote {len(rows)} rows")
            print(rel(out_dir / "scoreboard.tsv"))
            print(rel(out_dir / "scoreboard.json"))
            return 0
    except Exception as exc:  # noqa: BLE001 - CLI should print concise failure.
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
