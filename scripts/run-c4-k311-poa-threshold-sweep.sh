#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-32}"
WIDTH="${WIDTH:-3200}"
HEIGHT="${HEIGHT:-1800}"
MAX_CRUSH_PASSES="${MAX_CRUSH_PASSES:-3}"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
OUT="${1:-/home/erikg/impg/data/c4_k311_poa_threshold_${RUN_ID}}"

IMPG="${IMPG_BIN:-/home/erikg/.cargo/bin/impg}"
GFASORT="${GFASORT_BIN:-/home/erikg/.cargo/bin/gfasort}"
GFALOOK="${GFALOOK_BIN:-/home/erikg/.cargo/bin/gfalook}"
COMPARE="${COMPARE_GFA_PATHS:-/home/erikg/impg/target/release/examples/compare_gfa_paths}"

FASTA="${FASTA:-/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/c4_whole_region.fa}"
THRESHOLDS="${THRESHOLDS:-500 1000 2000 5000 10000}"
REUSE_SEED_GFA="${SEED_GFA:-}"
REUSE_SEED_RUNTIME_TSV="${SEED_RUNTIME_TSV:-}"
REUSE_SEED_RUNTIME_LABEL="${SEED_RUNTIME_LABEL:-}"

GRAPHS="$OUT/graphs"
REPORTS="$OUT/reports"
SORTED="$OUT/sorted"
RENDERS="$OUT/renders"
LOGS="$OUT/logs"
VALIDATION="$OUT/validation"
DEBUG="$OUT/debug"
TMP_ROOT="${TMP_ROOT:-/tmp/c4p_${RUN_ID##*T}}"

if [[ -e "$OUT" ]] && [[ -n "$(find "$OUT" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null || true)" ]]; then
  echo "Output directory already exists and is not empty: $OUT" >&2
  exit 2
fi

mkdir -p "$GRAPHS" "$REPORTS" "$SORTED" "$RENDERS" "$LOGS" "$VALIDATION" "$DEBUG" "$TMP_ROOT"

COMMANDS="$OUT/commands.md"
RUNTIME="$OUT/runtime_summary.tsv"
ROWS="$OUT/report_rows.tsv"
SWEEP="$OUT/poa_sweep_status.tsv"
UPLOADS="$OUT/uploaded_urls.tsv"
SCOREBOARD="$OUT/scoreboard.tsv"
SUMMARY="$OUT/summary.md"

printf '# C4 k311 POA threshold sweep commands\n\n' > "$COMMANDS"
printf 'tool\tpath\tversion_or_note\n' > "$OUT/tool_versions.tsv"
printf 'label\twall_time\tmax_rss_kb\texit_status\tstdout\tstderr\n' > "$RUNTIME"
printf 'row_label\tX\tstage\tgfa\treport_tsv\treport_md\tseed_time_label\tcrush_time_labels\tsort_time_label\trender_time_label\tupload_url\tcompare_stdout\tcompare_stderr\n' > "$ROWS"
printf 'X\tstatus\tresolved\tbailed\tcandidates_seen\trounds\tpasses\tconverged\tnote\tcrush_time_labels\tfinal_gfa\n' > "$SWEEP"
printf 'label\tlocal_png\turl\tremote_path\n' > "$UPLOADS"

require_file() {
  local path="$1"
  if [[ ! -s "$path" ]]; then
    echo "Required file is missing or empty: $path" >&2
    exit 2
  fi
}

require_exec() {
  local path="$1"
  if [[ ! -x "$path" ]]; then
    echo "Required executable is missing or not executable: $path" >&2
    exit 2
  fi
}

append_command() {
  local label="$1"
  shift
  {
    printf '## %s\n\n```bash\n/usr/bin/time -v' "$label"
    local arg
    for arg in "$@"; do
      printf ' %q' "$arg"
    done
    printf '\n```\n\n'
  } >> "$COMMANDS"
}

run_timed() {
  local label="$1"
  shift
  local stdout="$LOGS/${label}.stdout.log"
  local stderr="$LOGS/${label}.stderr.log"

  append_command "$label" "$@"

  set +e
  /usr/bin/time -v "$@" >"$stdout" 2>"$stderr"
  local status=$?
  set -e

  local wall rss exit_status
  wall="$(awk -F': ' '/Elapsed \(wall clock\) time/{v=$2} END{print v}' "$stderr")"
  rss="$(awk -F': ' '/Maximum resident set size/{v=$2} END{print v}' "$stderr")"
  exit_status="$(awk -F': ' '/Exit status/{v=$2} END{print v}' "$stderr")"
  if [[ -z "$exit_status" ]]; then
    exit_status="$status"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$label" "$wall" "$rss" "$exit_status" "$stdout" "$stderr" >> "$RUNTIME"

  if [[ "$status" -ne 0 ]]; then
    echo "Command failed ($label), see $stderr" >&2
    exit "$status"
  fi
}

record_tool_versions() {
  local impg_version
  impg_version="$("$IMPG" --version 2>&1 || true)"
  printf 'impg\t%s\t%s\n' "$IMPG" "$impg_version" >> "$OUT/tool_versions.tsv"
  printf 'gfasort\t%s\t%s\n' "$GFASORT" "no --version support" >> "$OUT/tool_versions.tsv"
  printf 'gfalook\t%s\t%s\n' "$GFALOOK" "no --version support" >> "$OUT/tool_versions.tsv"
  printf 'compare_gfa_paths\t%s\t%s\n' "$COMPARE" "release example binary" >> "$OUT/tool_versions.tsv"
}

report_graph() {
  local row="$1"
  local gfa="$2"
  local tsv="$REPORTS/${row}.graph-report.tsv"
  local md="$REPORTS/${row}.graph-report.md"
  run_timed "${row}.graph-report.tsv" "$IMPG" graph-report -g "$gfa" -o "$tsv" --format tsv --povu --top 20 -t "$THREADS" -v 1
  run_timed "${row}.graph-report.md" "$IMPG" graph-report -g "$gfa" -o "$md" --format markdown --povu --top 20 -t "$THREADS" -v 1
}

LAST_UPLOAD_URL=""
sort_render_upload() {
  local row="$1"
  local gfa="$2"
  local sorted="$SORTED/${row}.Ygs.gfa"
  local png="$RENDERS/${row}.Ygs.mean-depth.png"
  local remote="c4-k311-poa-threshold-${row}.Ygs.mean-depth.png"
  local remote_path="erik@hypervolu.me:www/impg/${remote}"

  run_timed "${row}.gfasort" "$GFASORT" -i "$gfa" -o "$sorted" -p Ygs -t "$THREADS" -v 1
  run_timed "${row}.gfalook" "$GFALOOK" -i "$sorted" -o "$png" -m -x "$WIDTH" -y "$HEIGHT" -a 3 -t "$THREADS" -v 1
  run_timed "${row}.upload" scp -o BatchMode=yes -o ConnectTimeout=20 "$png" "$remote_path"

  LAST_UPLOAD_URL="http://hypervolu.me/~erik/impg/${remote}"
  printf '%s\t%s\t%s\t%s\n' "$row" "$png" "$LAST_UPLOAD_URL" "$remote_path" >> "$UPLOADS"
}

compare_paths() {
  local row="$1"
  local seed="$2"
  local out_gfa="$3"
  local stdout="$VALIDATION/${row}.compare_gfa_paths.stdout.log"
  local stderr="$VALIDATION/${row}.compare_gfa_paths.stderr.log"

  append_command "${row}.compare_gfa_paths" "$COMPARE" "$seed" "$out_gfa"
  set +e
  /usr/bin/time -v "$COMPARE" "$seed" "$out_gfa" >"$stdout" 2>"$stderr"
  local status=$?
  set -e

  local wall rss exit_status
  wall="$(awk -F': ' '/Elapsed \(wall clock\) time/{v=$2} END{print v}' "$stderr")"
  rss="$(awk -F': ' '/Maximum resident set size/{v=$2} END{print v}' "$stderr")"
  exit_status="$(awk -F': ' '/Exit status/{v=$2} END{print v}' "$stderr")"
  if [[ -z "$exit_status" ]]; then
    exit_status="$status"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${row}.compare_gfa_paths" "$wall" "$rss" "$exit_status" "$stdout" "$stderr" >> "$RUNTIME"

  if [[ "$status" -ne 0 ]]; then
    echo "Path comparison failed ($row), see $stdout and $stderr" >&2
    exit "$status"
  fi
}

register_row() {
  local row="$1"
  local x="$2"
  local stage="$3"
  local gfa="$4"
  local seed_label="$5"
  local crush_labels="$6"
  local sort_label="$7"
  local render_label="$8"
  local url="$9"
  local compare_stdout="${10}"
  local compare_stderr="${11}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$row" "$x" "$stage" "$gfa" "$REPORTS/${row}.graph-report.tsv" "$REPORTS/${row}.graph-report.md" \
    "$seed_label" "$crush_labels" "$sort_label" "$render_label" "$url" "$compare_stdout" "$compare_stderr" >> "$ROWS"
}

build_seed() {
  local seed="$GRAPHS/c4.k311.seed.gfa"
  local debug_dir="$DEBUG/c4.k311.seed"
  local tmp_dir="$TMP_ROOT/k311"
  mkdir -p "$debug_dir" "$tmp_dir"

  if [[ -n "$REUSE_SEED_GFA" ]]; then
    require_file "$REUSE_SEED_GFA"
    cp -f "$REUSE_SEED_GFA" "$seed"
    {
      printf '## c4.k311.seed.graph\n\n'
      printf 'Reused completed k311 seed from `%s`.\n\n' "$REUSE_SEED_GFA"
      if [[ -n "$REUSE_SEED_RUNTIME_TSV" ]]; then
        printf 'Runtime provenance copied from `%s` label `%s`.\n\n' "$REUSE_SEED_RUNTIME_TSV" "$REUSE_SEED_RUNTIME_LABEL"
      fi
    } >> "$COMMANDS"
    if [[ -n "$REUSE_SEED_RUNTIME_TSV" && -n "$REUSE_SEED_RUNTIME_LABEL" && -s "$REUSE_SEED_RUNTIME_TSV" ]]; then
      python3 - "$REUSE_SEED_RUNTIME_TSV" "$REUSE_SEED_RUNTIME_LABEL" "$RUNTIME" <<'PY'
import csv
import sys
from pathlib import Path

source = Path(sys.argv[1])
source_label = sys.argv[2]
dest = Path(sys.argv[3])
with source.open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row["label"] == source_label:
            with dest.open("a", newline="") as out:
                writer = csv.writer(out, delimiter="\t", lineterminator="\n")
                writer.writerow([
                    "c4.k311.seed.graph",
                    row.get("wall_time", ""),
                    row.get("max_rss_kb", ""),
                    row.get("exit_status", ""),
                    row.get("stdout", ""),
                    row.get("stderr", ""),
                ])
            break
    else:
        raise SystemExit(f"missing runtime label {source_label!r} in {source}")
PY
    else
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "c4.k311.seed.graph" "" "" "0" "" "" >> "$RUNTIME"
    fi
    printf '%s\n' "$seed"
    return
  fi

  run_timed "c4.k311.seed.graph" "$IMPG" graph \
    --sequence-files "$FASTA" \
    --gfa-engine seqwish \
    --fastga \
    --num-mappings 1:many \
    --scaffold-filter 1:many \
    --scaffold-jump 0 \
    --min-match-len 311 \
    --temp-dir "$tmp_dir" \
    --debug-dir "$debug_dir" \
    -g "$seed" \
    -t "$THREADS" \
    -v 1

  printf '%s\n' "$seed"
}

parse_crush_summary() {
  local stderr="$1"
  python3 - "$stderr" <<'PY'
import re
import sys
from pathlib import Path

text = Path(sys.argv[1]).read_text(errors="replace")
matches = re.findall(r"crush: (\d+) resolved, (\d+) bailed, (\d+) candidates seen across (\d+) rounds", text)
if not matches:
    print("0\t0\t0\t0")
else:
    print("\t".join(matches[-1]))
PY
}

run_crush_threshold() {
  local x="$1"
  local seed="$2"
  local input="$seed"
  local final=""
  local labels=()
  local total_resolved=0
  local total_bailed=0
  local total_candidates=0
  local total_rounds=0
  local converged="no"
  local note="reached MAX_CRUSH_PASSES=${MAX_CRUSH_PASSES}; final pass did not report no eligible candidates"

  for ((pass = 1; pass <= MAX_CRUSH_PASSES; pass++)); do
    local pass_gfa="$GRAPHS/c4.k311.poa${x}.pass${pass}.gfa"
    local label="c4.k311.poa${x}.pass${pass}.crush"
    labels+=("$label")
    run_timed "$label" "$IMPG" crush \
      -g "$input" \
      -o "$pass_gfa" \
      --method poa \
      --max-iterations 5 \
      --max-traversal-len "$x" \
      --max-median-traversal-len "$x" \
      --max-total-sequence 1m \
      --max-traversals 10k \
      -t "$THREADS" \
      -v 1

    final="$pass_gfa"
    IFS=$'\t' read -r resolved bailed candidates rounds < <(parse_crush_summary "$LOGS/${label}.stderr.log")
    total_resolved=$((total_resolved + resolved))
    total_bailed=$((total_bailed + bailed))
    total_candidates=$((total_candidates + candidates))
    total_rounds=$((total_rounds + rounds))

    if grep -q 'no eligible candidates' "$LOGS/${label}.stderr.log"; then
      converged="yes"
      note="crush log reported no eligible candidates"
      break
    fi
    input="$pass_gfa"
  done

  local polished="$GRAPHS/c4.k311.poa${x}.gfa"
  cp -f "$final" "$polished"
  local joined
  joined="$(IFS=';'; printf '%s' "${labels[*]}")"
  local status="resolved"
  if [[ "$converged" != "yes" ]]; then
    status="round_cap"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$x" "$status" "$total_resolved" "$total_bailed" "$total_candidates" "$total_rounds" \
    "${#labels[@]}" "$converged" "$note" "$joined" "$polished" >> "$SWEEP"

  CRUSH_TIME_LABELS="$joined"
  CRUSH_POLISHED="$polished"
}

write_scoreboard() {
  python3 - "$OUT" <<'PY'
import csv
import sys
from pathlib import Path

out = Path(sys.argv[1])

def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

runtime = {row["label"]: row for row in read_tsv(out / "runtime_summary.tsv")}
rows = read_tsv(out / "report_rows.tsv")
sweep = {row["X"]: row for row in read_tsv(out / "poa_sweep_status.tsv")}

def wall_seconds(value):
    if not value:
        return 0.0
    parts = value.split(":")
    try:
        if len(parts) == 3:
            return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
        if len(parts) == 2:
            return int(parts[0]) * 60 + float(parts[1])
        return float(value)
    except ValueError:
        return 0.0

def fmt_seconds(seconds):
    return f"{seconds:.2f}" if seconds else ""

def time_group(labels):
    labels = [label for label in labels.split(";") if label]
    if not labels:
        return "", ""
    seconds = sum(wall_seconds(runtime.get(label, {}).get("wall_time", "")) for label in labels)
    rss_values = []
    for label in labels:
        rss = runtime.get(label, {}).get("max_rss_kb", "")
        if rss:
            try:
                rss_values.append(int(rss))
            except ValueError:
                pass
    return fmt_seconds(seconds), str(max(rss_values)) if rss_values else ""

metric_cols = [
    "segments",
    "links",
    "paths",
    "path_steps",
    "total_segment_bp",
    "direct_self_loop_edges",
    "adjacent_same_node_path_steps",
    "singleton_bp",
    "node_coverage_bp_weighted_mean",
    "path_white_space_bp_p99",
    "path_white_space_bp_max",
    "duplicate_sequence_frac",
    "local_repeat_context_nodes",
    "local_repeat_context_occurrences",
]

fieldnames = [
    "X",
    "status",
    "resolved",
    "bailed",
    "rounds",
    "passes",
    "converged",
    "segments",
    "links",
    "path_steps",
    "segment_bp",
    "direct_self_loop_edges",
    "adjacent_same_node_path_steps",
    "singleton_bp",
    "bp_weighted_node_coverage",
    "path_white_space_p99",
    "path_white_space_max",
    "duplicate_sequence_fraction",
    "local_repeat_contexts",
    "seed_wall_s",
    "seed_max_rss_kb",
    "crush_wall_s",
    "crush_max_rss_kb",
    "sort_wall_s",
    "sort_max_rss_kb",
    "render_wall_s",
    "render_max_rss_kb",
    "paths_preserved",
    "compare_expected_paths",
    "compare_observed_paths",
    "compare_spelling_mismatches",
    "png_url",
    "gfa",
    "note",
]

scoreboard = out / "scoreboard.tsv"
with scoreboard.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        if row["stage"] != "poa":
            continue
        report_path = Path(row["report_tsv"])
        with report_path.open(newline="") as report_handle:
            report_rows = list(csv.DictReader(report_handle, delimiter="\t"))
        if len(report_rows) != 1:
            raise SystemExit(f"expected one data row in {report_path}, found {len(report_rows)}")
        report = report_rows[0]
        compare = {}
        if row["compare_stdout"]:
            with Path(row["compare_stdout"]).open() as compare_handle:
                for line in compare_handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) == 2:
                        compare[parts[0]] = parts[1]
        x = row["X"]
        sweep_row = sweep.get(x, {})
        seed_wall, seed_rss = time_group(row["seed_time_label"])
        crush_wall, crush_rss = time_group(row["crush_time_labels"])
        sort_wall, sort_rss = time_group(row["sort_time_label"])
        render_wall, render_rss = time_group(row["render_time_label"])
        spelling = compare.get("spelling_mismatches", "")
        expected = compare.get("expected_paths", "")
        observed = compare.get("observed_paths", "")
        missing = compare.get("missing_paths", "")
        extra = compare.get("extra_paths", "")
        paths_preserved = (
            expected == "465"
            and observed == "465"
            and missing == "0"
            and extra == "0"
            and spelling == "0"
        )
        writer.writerow({
            "X": x,
            "status": sweep_row.get("status", ""),
            "resolved": sweep_row.get("resolved", ""),
            "bailed": sweep_row.get("bailed", ""),
            "rounds": sweep_row.get("rounds", ""),
            "passes": sweep_row.get("passes", ""),
            "converged": sweep_row.get("converged", ""),
            "segments": report.get("segments", ""),
            "links": report.get("links", ""),
            "path_steps": report.get("path_steps", ""),
            "segment_bp": report.get("total_segment_bp", ""),
            "direct_self_loop_edges": report.get("direct_self_loop_edges", ""),
            "adjacent_same_node_path_steps": report.get("adjacent_same_node_path_steps", ""),
            "singleton_bp": report.get("singleton_bp", ""),
            "bp_weighted_node_coverage": report.get("node_coverage_bp_weighted_mean", ""),
            "path_white_space_p99": report.get("path_white_space_bp_p99", ""),
            "path_white_space_max": report.get("path_white_space_bp_max", ""),
            "duplicate_sequence_fraction": report.get("duplicate_sequence_frac", ""),
            "local_repeat_contexts": f"{report.get('local_repeat_context_nodes', '')}/{report.get('local_repeat_context_occurrences', '')}",
            "seed_wall_s": seed_wall,
            "seed_max_rss_kb": seed_rss,
            "crush_wall_s": crush_wall,
            "crush_max_rss_kb": crush_rss,
            "sort_wall_s": sort_wall,
            "sort_max_rss_kb": sort_rss,
            "render_wall_s": render_wall,
            "render_max_rss_kb": render_rss,
            "paths_preserved": "yes" if paths_preserved else "no",
            "compare_expected_paths": expected,
            "compare_observed_paths": observed,
            "compare_spelling_mismatches": spelling,
            "png_url": row["upload_url"],
            "gfa": row["gfa"],
            "note": sweep_row.get("note", ""),
        })

seed_row = next((row for row in rows if row["stage"] == "seed"), None)
seed_metrics = {}
if seed_row:
    with Path(seed_row["report_tsv"]).open(newline="") as handle:
        seed_metrics = list(csv.DictReader(handle, delimiter="\t"))[0]

score_rows = read_tsv(scoreboard)
summary = out / "summary.md"
with summary.open("w") as handle:
    handle.write("# C4 k311 POA Threshold Sweep\n\n")
    handle.write(f"Output directory: `{out}`\n\n")
    if seed_metrics:
        handle.write("## Seed\n\n")
        handle.write(
            f"k311 seed paths={seed_metrics.get('paths')}, segments={seed_metrics.get('segments')}, "
            f"links={seed_metrics.get('links')}, path_steps={seed_metrics.get('path_steps')}, "
            f"segment_bp={seed_metrics.get('total_segment_bp')}, direct_self_loop_edges={seed_metrics.get('direct_self_loop_edges')}, "
            f"singleton_bp={seed_metrics.get('singleton_bp')}, white_space_p99={seed_metrics.get('path_white_space_bp_p99')}.\n\n"
        )
    handle.write("## Scoreboard\n\n")
    handle.write("| X | status | resolved | bailed | rounds | segs | links | path_steps | segment_bp | self_loops | adjacent_same_node | singleton_bp | bp_cov | ws_p99/max | dup_frac | repeat_ctx | crush_s/RSS_kb | paths | PNG |\n")
    handle.write("|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---|---|---|\n")
    for r in score_rows:
        handle.write(
            f"| {r['X']} | {r['status']} | {r['resolved']} | {r['bailed']} | {r['rounds']} | "
            f"{r['segments']} | {r['links']} | {r['path_steps']} | {r['segment_bp']} | "
            f"{r['direct_self_loop_edges']} | {r['adjacent_same_node_path_steps']} | {r['singleton_bp']} | "
            f"{r['bp_weighted_node_coverage']} | {r['path_white_space_p99']}/{r['path_white_space_max']} | "
            f"{r['duplicate_sequence_fraction']} | {r['local_repeat_contexts']} | "
            f"{r['crush_wall_s']}/{r['crush_max_rss_kb']} | {r['paths_preserved']} | {r['png_url']} |\n"
        )
    handle.write("\n## Answer\n\n")
    if score_rows:
        best = min(
            score_rows,
            key=lambda r: (
                int(r["path_white_space_p99"] or "999999999"),
                int(r["direct_self_loop_edges"] or "999999999"),
                float(r["crush_wall_s"] or "999999999"),
            ),
        )
        handle.write(
            f"The best observed ceiling by white-space p99, with path preservation as the safety gate, was X={best['X']} bp. "
            "Use the full scoreboard above for the speed and fragmentation trade-off.\n\n"
        )
    handle.write("All POA outputs were compared against the k311 seed with `compare_gfa_paths`; `paths` must be `yes` in the scoreboard.\n")
PY
}

require_file "$FASTA"
require_exec "$IMPG"
require_exec "$GFASORT"
require_exec "$GFALOOK"
require_exec "$COMPARE"
record_tool_versions

{
  printf 'Input FASTA: `%s`\n\n' "$FASTA"
  printf 'Seed command follows the requested full-C4 FastGA/seqwish recipe with `--min-match-len 311`, `--num-mappings 1:many`, `--scaffold-filter 1:many`, and `--scaffold-jump 0`.\n\n'
  printf 'POA thresholds: `%s` bp. Each threshold runs `impg crush --method poa --max-iterations 5` with both traversal ceilings set to X; the runner repeats up to `%s` 5-round passes only if the log does not report exhaustion.\n\n' "$THRESHOLDS" "$MAX_CRUSH_PASSES"
} >> "$COMMANDS"

SEED="$(build_seed)"
report_graph "c4.k311.seed" "$SEED"
sort_render_upload "c4.k311.seed" "$SEED"
register_row "c4.k311.seed" "seed" "seed" "$SEED" "c4.k311.seed.graph" "" "c4.k311.seed.gfasort" "c4.k311.seed.gfalook" "$LAST_UPLOAD_URL" "" ""

for x in $THRESHOLDS; do
  run_crush_threshold "$x" "$SEED"
  row="c4.k311.poa${x}"
  report_graph "$row" "$CRUSH_POLISHED"
  compare_paths "$row" "$SEED" "$CRUSH_POLISHED"
  sort_render_upload "$row" "$CRUSH_POLISHED"
  register_row "$row" "$x" "poa" "$CRUSH_POLISHED" "c4.k311.seed.graph" "$CRUSH_TIME_LABELS" "${row}.gfasort" "${row}.gfalook" "$LAST_UPLOAD_URL" "$VALIDATION/${row}.compare_gfa_paths.stdout.log" "$VALIDATION/${row}.compare_gfa_paths.stderr.log"
done

write_scoreboard
printf 'Wrote %s\n' "$OUT"
