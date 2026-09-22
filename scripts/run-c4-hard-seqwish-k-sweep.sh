#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-32}"
WIDTH="${WIDTH:-3200}"
HEIGHT="${HEIGHT:-1800}"
MAX_POA_PASSES="${MAX_POA_PASSES:-3}"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
OUT="${1:-/home/erikg/impg/data/c4_hard_seqwish_k_${RUN_ID}}"

IMPG="${IMPG_BIN:-/home/erikg/.cargo/bin/impg}"
GFASORT="${GFASORT_BIN:-/home/erikg/.cargo/bin/gfasort}"
GFALOOK="${GFALOOK_BIN:-/home/erikg/.cargo/bin/gfalook}"
COMPARE="${COMPARE_GFA_PATHS:-/home/erikg/impg/target/release/examples/compare_gfa_paths}"

FASTA="${FASTA:-/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/c4_whole_region.fa}"
BASELINE_SEED="${BASELINE_SEED:-/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_minmatch1_scaffold0.initial.gfa}"
BASELINE_POA="${BASELINE_POA:-/home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.gfa}"
BASELINE_NORM="${BASELINE_NORM:-/home/erikg/impg/data/c4_self_loop_normalization_20260605T130000Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.selfloop-normalized.gfa}"

GRAPHS="$OUT/graphs"
REPORTS="$OUT/reports"
SORTED="$OUT/sorted"
RENDERS="$OUT/renders"
LOGS="$OUT/logs"
VALIDATION="$OUT/validation"
DEBUG="$OUT/debug"
# FAtoGDB has a small path buffer in this environment. Keep alignment temp
# paths deliberately short; the long-lived debug/output paths stay under OUT.
TMP_ROOT="${TMP_ROOT:-/tmp/c4k_${RUN_ID##*T}}"

if [[ -e "$OUT" ]] && [[ -n "$(find "$OUT" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null || true)" ]]; then
  echo "Output directory already exists and is not empty: $OUT" >&2
  exit 2
fi

mkdir -p "$GRAPHS" "$REPORTS" "$SORTED" "$RENDERS" "$LOGS" "$VALIDATION" "$DEBUG" "$TMP_ROOT"

COMMANDS="$OUT/commands.md"
RUNTIME="$OUT/runtime_summary.tsv"
ROWS="$OUT/report_rows.tsv"
UPLOADS="$OUT/uploaded_urls.tsv"
CONTINUATION="$OUT/continuation_notes.tsv"

printf '# C4 hard seqwish min-match sweep commands\n\n' > "$COMMANDS"
printf 'label\twall_time\tmax_rss_kb\texit_status\tstdout\tstderr\n' > "$RUNTIME"
printf 'row_label\tvariant_k\tstage\tgfa\treport_tsv\treport_md\tgraph_time_labels\tpoa_time_labels\tnormalize_time_labels\tsort_time_labels\trender_time_labels\tupload_url\tcompare_stdout\tcompare_stderr\n' > "$ROWS"
printf 'label\tlocal_png\turl\tremote_path\n' > "$UPLOADS"
printf 'label\tpasses\tconverged\tnote\n' > "$CONTINUATION"

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

require_file "$FASTA"
require_file "$BASELINE_SEED"
require_file "$BASELINE_POA"
require_file "$BASELINE_NORM"
require_exec "$IMPG"
require_exec "$GFASORT"
require_exec "$GFALOOK"
require_exec "$COMPARE"

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
  local remote="c4-hard-seqwish-${row}.Ygs.mean-depth.png"
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
  local k="$2"
  local stage="$3"
  local gfa="$4"
  local graph_labels="$5"
  local poa_labels="$6"
  local norm_labels="$7"
  local sort_labels="$8"
  local render_labels="$9"
  local url="${10}"
  local compare_stdout="${11}"
  local compare_stderr="${12}"
  local report_tsv="$REPORTS/${row}.graph-report.tsv"
  local report_md="$REPORTS/${row}.graph-report.md"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$row" "$k" "$stage" "$gfa" "$report_tsv" "$report_md" "$graph_labels" "$poa_labels" "$norm_labels" "$sort_labels" "$render_labels" "$url" "$compare_stdout" "$compare_stderr" >> "$ROWS"
}

build_seed() {
  local k="$1"
  local raw_paf="$2"
  local label="one_many_minmatch${k}_scaffold0"
  local seed="$GRAPHS/${label}.initial.gfa"
  local debug_dir="$DEBUG/${label}"
  local tmp_dir="$TMP_ROOT/k${k}"
  mkdir -p "$debug_dir" "$tmp_dir"

  local cmd=(
    "$IMPG" graph
    --sequence-files "$FASTA"
    --gfa-engine seqwish
    --fastga
    --num-mappings 1:many
    --scaffold-filter 1:many
    --scaffold-jump 0
    --min-match-len "$k"
    --temp-dir "$tmp_dir"
    --debug-dir "$debug_dir"
    -g "$seed"
    -t "$THREADS"
    -v 1
  )
  if [[ -n "$raw_paf" ]]; then
    cmd+=(--paf-file "$raw_paf")
  fi

  run_timed "${label}.graph" "${cmd[@]}"
  printf '%s\n' "$seed"
}

polish_poa1kb() {
  local label="$1"
  local seed="$2"
  local input="$seed"
  local final=""
  local labels=()
  local converged="no"
  local note="continued until pass limit because no convergence marker was found"

  for ((pass = 1; pass <= MAX_POA_PASSES; pass++)); do
    local outpass="$GRAPHS/${label}.poa1kb.pass${pass}.gfa"
    local time_label="${label}.poa1kb.pass${pass}"
    labels+=("$time_label")
    run_timed "$time_label" "$IMPG" crush \
      -g "$input" \
      -o "$outpass" \
      --method poa \
      --max-iterations 5 \
      --max-traversal-len 1k \
      --max-median-traversal-len 1k \
      --max-total-sequence 1m \
      --max-traversals 10k \
      -t "$THREADS" \
      -v 1
    final="$outpass"
    if grep -q 'no eligible candidates' "$LOGS/${time_label}.stderr.log"; then
      converged="yes"
      note="crush log reported no eligible candidates"
      break
    fi
    input="$outpass"
  done

  local polished="$GRAPHS/${label}.poa1kb.gfa"
  cp -f "$final" "$polished"
  local joined
  joined="$(IFS=';'; printf '%s' "${labels[*]}")"
  printf '%s\t%s\t%s\t%s\n' "$label" "${#labels[@]}" "$converged" "$note" >> "$CONTINUATION"
  POA_TIME_LABELS="$joined"
  POA_POLISHED="$polished"
}

normalize_variant() {
  local row="$1"
  local gfa="$2"
  local norm="$GRAPHS/${row}.selfloop-normalized.gfa"
  local stats="$REPORTS/${row}.selfloop-normalized.stats.json"
  run_timed "${row}.normalize-self-loops" "$IMPG" normalize-self-loops -g "$gfa" -o "$norm" --report "$stats"
  printf '%s\n' "$norm"
}

process_variant() {
  local k="$1"
  local raw_paf="$2"
  local label="one_many_minmatch${k}_scaffold0"
  local seed
  seed="$(build_seed "$k" "$raw_paf")"

  local row_initial="${label}.initial"
  report_graph "$row_initial" "$seed"
  sort_render_upload "$row_initial" "$seed"
  register_row "$row_initial" "$k" "seed" "$seed" "${label}.graph" "" "" "${row_initial}.gfasort" "${row_initial}.gfalook" "$LAST_UPLOAD_URL" "" ""

  polish_poa1kb "$label" "$seed"
  local polished="$POA_POLISHED"
  local poa_labels="$POA_TIME_LABELS"
  local row_poa="${label}.poa1kb"
  report_graph "$row_poa" "$polished"
  compare_paths "$row_poa" "$seed" "$polished"
  sort_render_upload "$row_poa" "$polished"
  register_row "$row_poa" "$k" "poa1kb" "$polished" "${label}.graph" "$poa_labels" "" "${row_poa}.gfasort" "${row_poa}.gfalook" "$LAST_UPLOAD_URL" "$VALIDATION/${row_poa}.compare_gfa_paths.stdout.log" "$VALIDATION/${row_poa}.compare_gfa_paths.stderr.log"

  local norm
  norm="$(normalize_variant "$row_poa" "$polished")"
  local row_norm="${label}.poa1kb.selfloop-normalized"
  report_graph "$row_norm" "$norm"
  compare_paths "$row_norm" "$seed" "$norm"
  sort_render_upload "$row_norm" "$norm"
  register_row "$row_norm" "$k" "poa1kb_selfloop_normalized" "$norm" "${label}.graph" "$poa_labels" "${row_poa}.normalize-self-loops" "${row_norm}.gfasort" "${row_norm}.gfalook" "$LAST_UPLOAD_URL" "$VALIDATION/${row_norm}.compare_gfa_paths.stdout.log" "$VALIDATION/${row_norm}.compare_gfa_paths.stderr.log"
}

register_existing_report() {
  local row="$1"
  local k="$2"
  local stage="$3"
  local gfa="$4"
  local seed="$5"
  local known_url="$6"

  report_graph "$row" "$gfa"
  local compare_stdout=""
  local compare_stderr=""
  if [[ -n "$seed" ]]; then
    compare_paths "$row" "$seed" "$gfa"
    compare_stdout="$VALIDATION/${row}.compare_gfa_paths.stdout.log"
    compare_stderr="$VALIDATION/${row}.compare_gfa_paths.stderr.log"
  fi
  register_row "$row" "$k" "$stage" "$gfa" "" "" "" "" "" "$known_url" "$compare_stdout" "$compare_stderr"
}

printf 'Input FASTA: `%s`\n\n' "$FASTA" >> "$COMMANDS"
printf 'Baseline seed: `%s`\n\n' "$BASELINE_SEED" >> "$COMMANDS"
printf 'PAF reuse policy: build k=311 with FastGA once and save debug/raw.paf; build k=511 and k=1001 from that raw PAF so the only induction change is seqwish min-match-len.\n\n' >> "$COMMANDS"

register_existing_report \
  "baseline_minmatch1.initial" \
  "1" \
  "baseline_seed" \
  "$BASELINE_SEED" \
  "" \
  "https://hypervolu.me/~erik/impg/c4-whole-region-sweepga-seed-one-many-minmatch1-scaffold0-initial.png"

register_existing_report \
  "baseline_minmatch1.poa1kb" \
  "1" \
  "baseline_poa1kb" \
  "$BASELINE_POA" \
  "$BASELINE_SEED" \
  "https://hypervolu.me/~erik/impg/full_c4.one_many_minmatch1_scaffold0.poa1kb.Ygs.mean-depth.png"

register_existing_report \
  "baseline_minmatch1.poa1kb.selfloop-normalized" \
  "1" \
  "baseline_poa1kb_selfloop_normalized" \
  "$BASELINE_NORM" \
  "$BASELINE_SEED" \
  "https://hypervolu.me/~erik/impg/c4-self-loop-normalized-poa1kb.Ygs.mean-depth.png"

process_variant 311 ""
RAW_PAF="$DEBUG/one_many_minmatch311_scaffold0/raw.paf"
require_file "$RAW_PAF"
process_variant 511 "$RAW_PAF"
process_variant 1001 "$RAW_PAF"

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
    if not seconds:
        return ""
    return f"{seconds:.2f}"

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
    "direct_self_loop_nodes",
    "adjacent_same_node_path_steps",
    "adjacent_same_step_path_steps",
    "self_loop_repeat_runs",
    "self_loop_max_repeat_run_len",
    "self_loop_length_buckets",
    "singleton_bp",
    "node_coverage_bp_weighted_mean",
    "path_white_space_bp_p99",
    "path_white_space_bp_max",
    "duplicate_sequence_frac",
    "local_repeat_context_nodes",
    "local_repeat_context_occurrences",
    "povu_sites",
    "povu_leaf_sites",
]

time_cols = [
    ("graph_wall_s", "graph_max_rss_kb", "graph_time_labels"),
    ("poa_wall_s", "poa_max_rss_kb", "poa_time_labels"),
    ("normalize_wall_s", "normalize_max_rss_kb", "normalize_time_labels"),
    ("sort_wall_s", "sort_max_rss_kb", "sort_time_labels"),
    ("render_wall_s", "render_max_rss_kb", "render_time_labels"),
]

scoreboard = out / "metrics.scoreboard.tsv"
with scoreboard.open("w", newline="") as handle:
    fieldnames = [
        "row_label",
        "variant_k",
        "stage",
        "gfa",
        "status",
        "warnings",
        "failures",
        *metric_cols,
        *(name for pair in time_cols for name in pair[:2]),
        "upload_url",
        "compare_expected_paths",
        "compare_observed_paths",
        "compare_missing_paths",
        "compare_extra_paths",
        "compare_spelling_mismatches",
    ]
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        report_path = Path(row["report_tsv"])
        with report_path.open(newline="") as report_handle:
            report_rows = list(csv.DictReader(report_handle, delimiter="\t"))
        if len(report_rows) != 1:
            raise SystemExit(f"expected one data row in {report_path}, found {len(report_rows)}")
        report = report_rows[0]
        out_row = {
            "row_label": row["row_label"],
            "variant_k": row["variant_k"],
            "stage": row["stage"],
            "gfa": row["gfa"],
            "status": report.get("status", ""),
            "warnings": report.get("warnings", ""),
            "failures": report.get("failures", ""),
            "upload_url": row["upload_url"],
        }
        for col in metric_cols:
            out_row[col] = report.get(col, "")
        for wall_col, rss_col, label_col in time_cols:
            wall, rss = time_group(row[label_col])
            out_row[wall_col] = wall
            out_row[rss_col] = rss
        compare_stdout = row.get("compare_stdout", "")
        compare = {}
        if compare_stdout:
            with Path(compare_stdout).open() as compare_handle:
                for line in compare_handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) == 2:
                        compare[parts[0]] = parts[1]
        out_row["compare_expected_paths"] = compare.get("expected_paths", "")
        out_row["compare_observed_paths"] = compare.get("observed_paths", "")
        out_row["compare_missing_paths"] = compare.get("missing_paths", "")
        out_row["compare_extra_paths"] = compare.get("extra_paths", "")
        out_row["compare_spelling_mismatches"] = compare.get("spelling_mismatches", "")
        writer.writerow(out_row)
PY

printf 'Wrote %s\n' "$OUT"
