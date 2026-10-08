#!/usr/bin/env bash
# Diagnostic-only serial extraction from the saved balanced reads/index/panel.
# The existing component runner supplies its 64GiB guard and external RSS poller.
set -uo pipefail
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
ROOT=/home/erikg/impg-genome-inference
for C in chrI chrVI chrMT chrIII chrIX chrVIII chrV chrXI chrX chrII chrXIV chrXIII chrXVI chrXII chrXV chrVII chrIV; do
  P="$D/run-cosine-diag2-$C"
  A="$D/cosine-probe2-$C.jsonl"
  if [ -e "$P.done" ] && [ "$(cat "$P.exit" 2>/dev/null)" = 0 ] && [ -s "$A" ]; then
    echo "skip $C (clean)" >> "$D/cosine-probe-fleet.log"
    continue
  fi
  echo "start $C $(date -u +%FT%TZ)" >> "$D/cosine-probe-fleet.log"
  if [ "$C" = chrIV ]; then
    IMPG_CORES=0-63 IMPG_COSINE_DIAG_OUTPUT="$A" \
      "$ROOT/genome/instrumented/run-balanced-diploid-component.sh" "$C" cosine-diag2
  else
    IMPG_CORES=0-255 IMPG_COSINE_DIAG_OUTPUT="$A" \
      "$ROOT/genome/instrumented/run-balanced-diploid-component.sh" "$C" cosine-diag2
  fi
  status=$?
  echo "end $C exit=$status wall=$(cat "$P.wall" 2>/dev/null) $(date -u +%FT%TZ)" >> "$D/cosine-probe-fleet.log"
  if [ "$status" != 0 ] || [ ! -s "$A" ]; then
    echo "$C" > "$D/cosine-probe-fleet.failed"
    if [ "$status" = 0 ]; then status=1; fi
    exit "$status"
  fi
done
python3 -B "$ROOT/genome/instrumented/score-cosine-probe.py" \
  "$D"/cosine-probe2-chr*.jsonl > "$D/cosine-score-genome.json.tmp" || exit $?
mv "$D/cosine-score-genome.json.tmp" "$D/cosine-score-genome.json"
touch "$D/cosine-probe-fleet.done"
