#!/usr/bin/env bash
# THE WHOLE-GENOME RUNTIME VERIFICATION - STAGE 4: THE STREAMING CALL WAVE.
# Per territory (marker-gated, resumable, killable): the 35x region reads
# fetch (the committed fetch rule: primaries in window +/-50kb, flat fastq)
# -> the three-step PARTITION step (genome-infer partition: the local syng +
# window beds + routes by lookup + the registry; flock - the BED4 append)
# -> call-region (--emit-sequences ALWAYS ON, the export maps) with the RSS
# poller -> THE JIT REAPER (the receipts moved out; the census, the cache,
# the reads fasta/fq, the substrate and the partition precomputes DELETED -
# the disk stays at the in-flight window).
# usage: run-wave.sh   (env: RV, LANES=96, LANE_RAYON=1, NICE=10)
set -uo pipefail
RV=${RV:-/home/erikg/hprcv2/wgs35x-fleet/runtime-verification}
IMPG=${IMPG:-/home/erikg/impg-genome-inference/target/release/impg}
BAM=/home/erikg/hprcv2/HG002.GRCh38.35x.bam
LANES=${LANES:-96}
LANE_RAYON=${LANE_RAYON:-1}
LOCK=$RV/partitions.lock
mkdir -p "$RV/reads" "$RV/receipts" "$RV/wave-logs"

rss_peak() { # rss_peak <rssfile> -> echoes peak kB
  awk '{if ($2+0 > m) m = $2+0} END {print m+0}' "$1" 2>/dev/null
}

one_territory() {
  local T=$1
  local P=$RV/prep/$T
  [ -f "$P/call.done" ] && return 0
  [ -f "$P/cut.done" ] || return 0
  [ -f "$P/call.exit" ] && return 0
  local CHROM S E AXIS
  CHROM=$(awk -F'\t' -v t="$T" '$4==t {print $1}' "$RV/bed4-genome.tsv")
  S=$(awk -F'\t' -v t="$T" '$4==t {print $2}' "$RV/bed4-genome.tsv")
  E=$(awk -F'\t' -v t="$T" '$4==t {print $3}' "$RV/bed4-genome.tsv")
  AXIS=$(cat "$P/axis-path.txt")
  local CUT=$RV/map/$T/$T

  # the no-evidence class: an empty cut is honestly skipped, named
  local KEEP
  KEEP=$(python3 -c "import json;print(json.load(open('$CUT/cut-stats.json')).get('reads_kept',0))" 2>/dev/null || echo -1)
  if [ "$KEEP" = "0" ]; then
    echo "$(date -u +%FT%TZ) $T no-evidence (empty cut); call skipped" >> "$RV/wave-logs/no-evidence.txt"
    touch "$P/call.skipped-no-evidence"
    # reap the empty cut's outputs
    rm -rf "$RV/map/$T"
    return 0
  fi

  # 1. the reads fetch (the committed fetch rule)
  if [ ! -f "$P/reads.done" ]; then
    local start
    start=$(date +%s)
    samtools view -u -F 0x904 "$BAM" "$CHROM:$((S-50000))-$((E+50000))" \
      | samtools fastq - 2> "$P/reads-fq.err" \
      | gzip -1 > "$RV/reads/$T.fq.gz"
    local st=${PIPESTATUS[0]}
    echo "$st" > "$P/reads.exit"
    echo "$(( $(date +%s) - start ))" > "$P/reads.wall"
    [ "$st" -eq 0 ] && touch "$P/reads.done" || { touch "$P/call.exit"; echo "$(date -u +%FT%TZ) $T reads fetch failed ($st)" >> "$RV/wave-logs/failures.txt"; return 1; }
  fi

  # 2. the per-territory map manifest (call-region reads the reads path here)
  mkdir -p "$RV/map/$T"
  printf '{"version": 1, "reads": "%s/reads/%s.fq.gz", "partitions": ["%s"]}\n' \
    "$RV" "$T" "$T" > "$RV/map/$T/manifest.json"

  # 3. the three-step PARTITION step (flock: partitions.bed append)
  if [ ! -f "$P/gpartition.done" ]; then
    local start
    start=$(date +%s)
    flock "$LOCK" env PATH="$PATH" "$IMPG" genome-infer partition \
      --name "$T" --chrom "$CHROM" --start "$S" --end "$E" \
      --sources "$P/locus.fa" \
      --component "$AXIS" \
      --axis "$P/reference-axis.json" \
      --partition-graphs "$P/locality-graphs" \
      --census "$CUT/cut-census.jsonl" \
      --read-lengths 148 \
      --out-dir "$RV/partitions" \
      > "$P/gpartition.log" 2> "$P/gpartition.err" \
      && { touch "$P/gpartition.done"; echo "$(( $(date +%s) - start ))" > "$P/gpartition.wall"; } \
      || { echo 1 > "$P/gpartition.exit"; touch "$P/call.exit"; echo "$(( $(date +%s) - start ))" > "$P/gpartition.wall"; \
           echo "$(date -u +%FT%TZ) $T partition step failed: $(tail -1 "$P/gpartition.err" | head -c 200)" >> "$RV/wave-logs/failures.txt"; return 1; }
  fi

  # 4. the CALL (the three-step call-region; the emission always on)
  if [ ! -f "$P/call.done" ]; then
    rm -rf "$RV/calls/$T"
    local start
    start=$(date +%s)
    RAYON_NUM_THREADS=$LANE_RAYON IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS=1 \
      nice -n 10 stdbuf -oL -eL "$IMPG" genome-infer call-region \
      --map "$RV/map/$T" --partitions "$RV/partitions" --region "$T" \
      --out-dir "$RV/calls/$T" \
      > "$P/call.log" 2> "$P/call.err" &
    local pid=$!
    ( while kill -0 "$pid" 2>/dev/null; do
        rss_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
        echo "$(date -u +%FT%TZ) ${rss_kb:-0}" >> "$P/call.rss"
        if [ "${rss_kb:-0}" -gt $((64*1024*1024)); then
          echo "$(date -u +%FT%TZ) GUARD kill ${rss_kb}kB" >> "$P/call.rss"
          kill -9 "$pid" 2>/dev/null; break
        fi
        sleep 5
      done ) &
    local poll=$!
    local st=0; wait "$pid" || st=$?
    kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
    echo "$st" > "$P/call.exit.code"
    echo "$(( $(date +%s) - start ))" > "$P/call.wall"
    if [ "$st" -ne 0 ]; then
      rm -f "$P/call.done"; touch "$P/call.exit"
      echo "$(date -u +%FT%TZ) $T CALL FAILED ($st): $(tail -1 "$P/call.err" | head -c 200)" >> "$RV/wave-logs/failures.txt"
      return 1
    fi
    touch "$P/call.done"
  fi

  # 5. THE JIT REAPER: receipts out, intermediates deleted
  local R=$RV/receipts/$T
  mkdir -p "$R"
  cp "$RV/calls/$T/$T/calls.jsonl" "$R/calls.jsonl" 2>/dev/null
  cp "$RV/calls/$T/$T/calls.jsonl.sequences.jsonl" "$R/calls.jsonl.sequences.jsonl" 2>/dev/null
  cp "$CUT/cut-stats.json" "$R/cut-stats.json"
  for m in call.wall gpartition.wall reads.wall extract.wall syng.wall beds.wall export.wall walks.wall maps.wall; do
    [ -f "$P/$m" ] && cp "$P/$m" "$R/$m"
  done
  [ -f "$P/call.rss" ] && rss_peak "$P/call.rss" > "$R/call.peak_rss_kb"
  [ -f "$P/call.exit.code" ] && cp "$P/call.exit.code" "$R/call.exit.code"
  # the receipts are real files, not empty
  [ -s "$R/calls.jsonl" ] || { touch "$P/call.exit"; echo "$(date -u +%FT%TZ) $T REAPER REFUSED (empty calls.jsonl)" >> "$RV/wave-logs/failures.txt"; return 1; }
  rm -rf "$RV/map/$T" "$RV/calls/$T" "$RV/partitions/$T" "$RV/reads/$T.fq.gz"
  rm -rf "$P/locus.fa" "$P/locus.fa.fai" "$P"/locus-pos64.syng.* \
         "$P/partition-w10k-d1k" "$P/export-graphs" "$P/export-walks.jsonl" \
         "$P/walk-requests.tsv" "$P/locality-graphs" "$P/reference-axis.json" \
         "$P/locality-domain-structure.txt"
  echo "$(date -u +%FT%TZ) $T DONE (call $(cat "$P/call.wall")s, rss $(cat "$R/call.peak_rss_kb" 2>/dev/null)kB)" >> "$RV/wave-logs/wave-progress.log"
  return 0
}
export -f one_territory
export RV IMPG BAM LANE_RAYON LOCK

# the work loop: every cut-done, not-yet-called territory, at LANES wide
while true; do
  mapfile -t TODO < <(while read -r T; do
      [ -f "$RV/prep/$T/cut.done" ] || continue
      [ -f "$RV/prep/$T/call.done" ] && continue
      [ -f "$RV/prep/$T/call.exit" ] && continue
      [ -f "$RV/prep/$T/call.skipped-no-evidence" ] && continue
      echo "$T"
    done < "$RV/prep/accepted.txt")
  [ "${#TODO[@]}" -eq 0 ] && break
  printf '%s\n' "${TODO[@]}" | xargs -P "$LANES" -I{} bash -c 'one_territory {}'
  sleep 30
done
NCALL=$(while read -r T; do [ -f "$RV/prep/$T/call.done" ] && echo x; done < "$RV/prep/accepted.txt" | wc -l)
echo "=== wave complete: $NCALL territories called $(date -u +%FT%TZ)" >> "$RV/wave-logs/wave-progress.log"
touch "$RV/wave-logs/wave.done"
