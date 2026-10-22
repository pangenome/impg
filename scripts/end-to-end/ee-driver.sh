#!/usr/bin/env bash
# THE END-TO-END RECEIPT - THE PIPELINE DRIVER (the multi-window
# campaign's own master; marker-gated, resumable, killable): four
# stage pools over the full accepted universe (5,844 territories),
# each bounded per iteration and re-checked every wave (the valve
# discipline), under the workdir budget + the worker guards + the
# heartbeat (df, memavail, the terminal counters). The stages:
#   P  the substrate chain (run-ee-prep.sh)        PREP_W wide
#   PP panel-prep (run-ee-panelprep.sh)            PP_W   wide
#   M  the map legs (run-ee-mapcut.sh)             MAP_W  wide
#   C  the call chunks (run-ee-call-chunk.sh)      CALL_W wide (32/chunk)
# The JIT reaper runs inside each call chunk (transients reaped at
# the leg's completion; answers + receipts kept).
# usage: nohup bash ee-driver.sh &   (env: EE, PREP_W, PP_W, MAP_W, CALL_W, VALVE)
source /home/erikg/impg-genome-inference/scripts/end-to-end/ee-env.sh
PREP_W=${PREP_W:-12}
PP_W=${PP_W:-4}
MAP_W=${MAP_W:-5}
CALL_W=${CALL_W:-2}
VALVE=${VALVE:-160}
SCRIPTS=/home/erikg/impg-genome-inference/scripts/end-to-end
CHUNKSZ=${CHUNKSZ:-32}

# ---- the chunk lists (deterministic: bed4 order, adjacent) ----------
mkdir -p "$EE/chunks"
if [ ! -f "$EE/chunks/chunks.built" ]; then
  awk 'NR==FNR{a[$1]=1; next} a[$4]{print $4}' "$RV/prep/accepted.txt" \
    "$RV/bed4-genome.tsv" > "$EE/chunks/universe.order"
  rm -f "$EE/chunks"/chunk-*.list "$EE/chunks"/x??
  split -l "$CHUNKSZ" "$EE/chunks/universe.order" "$EE/chunks/x"
  i=0
  for f in "$EE/chunks"/x??; do
    [ -e "$f" ] || continue
    mv "$f" "$EE/chunks/chunk-$i.list"; i=$((i+1))
  done
  echo "$i" > "$EE/chunks/chunks.built"
  ee_log "chunk lists built: $i chunks of $CHUNKSZ over $(wc -l < "$EE/chunks/universe.order") territories"
fi
NCHUNKS=$(cat "$EE/chunks/chunks.built")
NTERR=$(wc -l < "$EE/chunks/universe.order")

# pending sets by predicate: <outfile> <max> <predicate-expr in $t>
# (also echoes to stdout, so it can be counted)
pending_set() {
  local out=$1 max=$2 pred=$3 t c=0
  : > "$out"
  while read -r t; do
    [ -n "$t" ] || continue
    if eval "$pred"; then
      echo "$t" >> "$out"; echo "$t"
      c=$((c+1))
      [ "$c" -ge "$max" ] && return 0
    fi
  done < "$EE/chunks/universe.order"
  return 0
}

term_census() {
  local sub=0 pp=0 mp=0 cl=0 fail=0 t
  while read -r t; do
    [ -n "$t" ] || continue
    [ -f "$EE/prep/$t/substrate.done" ] && sub=$((sub+1))
    [ -f "$EE/prep/$t/substrate.exit" ] && [ ! -f "$EE/prep/$t/substrate.done" ] && fail=$((fail+1))
    [ -f "$EE/ledger/$t/panelprep.poison" ] && [ ! -f "$EE/ledger/$t/panelprep.done" ] && fail=$((fail+1))
    [ -f "$EE/ledger/$t/panelprep.done" ] && pp=$((pp+1))
    [ -f "$EE/ledger/$t/map.done" ] && mp=$((mp+1))
    if [ -f "$EE/ledger/$t/call.done" ]; then cl=$((cl+1))
    elif [ -f "$EE/ledger/$t/call.exit" ] || [ -f "$EE/ledger/$t/map.exit" ]; then fail=$((fail+1)); fi
  done < "$EE/chunks/universe.order"
  echo "sub=$sub panelprep=$pp map=$mp call=$cl failclosed=$fail of $NTERR"
}

ee_log "driver start: $NTERR territories, $NCHUNKS chunks; widths P=$PREP_W PP=$PP_W M=$MAP_W C=$CALL_W valve=$VALVE budget=${EE_BUDGET_GB}G guard=${GUARD_GIB}GiB"
iter=0
while true; do
  iter=$((iter+1))
  ee_budget_check || { ee_log "BUDGET STOP at iter $iter"; break; }
  census=$(term_census)
  ee_log "iter $iter: $census df=$(df -h /home/erikg | tail -1 | awk '{print $4}') memavail=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)G"
  set -- $census
  c=${4#call=}; f=${5#failclosed=}
  if [ $((c+f)) -ge "$NTERR" ]; then
    ee_log "ALL TERMINAL: $census"
    touch "$EE/ledger/campaign.done"
    break
  fi

  # the valve: prep only while the substrate-done-but-uncut window is bounded
  nsub=$(pending_set /dev/null 999999 '[ -f "$EE/prep/$t/substrate.done" ]' | wc -l)
  nmapdone=$(pending_set /dev/null 999999 '[ -f "$EE/ledger/$t/map.done" ]' | wc -l)
  pend_cut=$(( nsub - nmapdone ))
  if [ "$pend_cut" -lt "$VALVE" ]; then
    pending_set "$EE/scratch/P.set" $((PREP_W*2)) \
      '[ ! -f "$EE/prep/$t/substrate.done" ] && [ ! -f "$EE/prep/$t/substrate.exit" ] && [ ! -f "$EE/prep/$t/substrate.skipped" ]'
  else
    : > "$EE/scratch/P.set"
  fi
  pending_set "$EE/scratch/PP.set" $((PP_W*2)) \
    '[ -f "$EE/prep/$t/substrate.done" ] && [ ! -f "$EE/ledger/$t/panelprep.done" ] && [ ! -f "$EE/ledger/$t/panelprep.poison" ]'
  pending_set "$EE/scratch/M.set" $((MAP_W*2)) \
    '[ -f "$EE/ledger/$t/panelprep.done" ] && [ ! -f "$EE/ledger/$t/map.done" ] && [ ! -f "$EE/ledger/$t/map.exit" ]'
  : > "$EE/scratch/C.set"
  for i in $(seq 0 $((NCHUNKS-1))); do
    [ -f "$EE/ledger/chunk-$i/call.done" ] && continue
    alldone=1
    while read -r t; do
      [ -n "$t" ] || continue
      # the chunk-readiness gate (the fail-closed doctrine at chunk
      # shape): a member is TERMINAL when its map is done OR it is a
      # named poison (map.exit / substrate-terminal) - the chunk
      # runner calls the callable subset and names the poison
      if [ -f "$EE/ledger/$t/map.done" ] || [ -f "$EE/ledger/$t/map.exit" ] \
         || [ -f "$EE/ledger/$t/panelprep.poison" ] \
         || [ -f "$EE/prep/$t/substrate.skipped" ] \
         || { [ -f "$EE/prep/$t/substrate.exit" ] && [ ! -f "$EE/prep/$t/substrate.done" ]; }; then
        continue
      fi
      alldone=0; break
    done < "$EE/chunks/chunk-$i.list"
    [ "$alldone" -eq 1 ] && echo "$i" >> "$EE/scratch/C.set"
    [ "$(wc -l < "$EE/scratch/C.set")" -ge $((CALL_W*2)) ] && break
  done

  ee_log "iter $iter dispatch: P=$(wc -l < "$EE/scratch/P.set") PP=$(wc -l < "$EE/scratch/PP.set") M=$(wc -l < "$EE/scratch/M.set") C=$(wc -l < "$EE/scratch/C.set") (pend_cut=$pend_cut)"

  pids=""
  if [ -s "$EE/scratch/P.set" ]; then
    xargs -P "$PREP_W" -n 1 bash "$SCRIPTS/run-ee-prep.sh" < "$EE/scratch/P.set" &
    pids="$pids $!"
  fi
  if [ -s "$EE/scratch/PP.set" ]; then
    xargs -P "$PP_W" -n 1 bash "$SCRIPTS/run-ee-panelprep.sh" < "$EE/scratch/PP.set" &
    pids="$pids $!"
  fi
  if [ -s "$EE/scratch/M.set" ]; then
    xargs -P "$MAP_W" -n 1 bash "$SCRIPTS/run-ee-mapcut.sh" < "$EE/scratch/M.set" &
    pids="$pids $!"
  fi
  if [ -s "$EE/scratch/C.set" ]; then
    xargs -P "$CALL_W" -n 1 bash "$SCRIPTS/run-ee-call-chunk.sh" < "$EE/scratch/C.set" &
    pids="$pids $!"
  fi
  [ -z "$pids" ] && { ee_log "nothing dispatchable (prep/valve blocked or upstream missing) - sleeping 300s"; sleep 300; }
  for p in $pids; do wait "$p" || true; done
done
ee_log "driver exit after $iter iters: $(term_census)"
bash "$SCRIPTS/ee-ledger.sh" > "$EE/ledger/ledger-final.tsv" 2>&1 || true
