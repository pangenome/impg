#!/usr/bin/env bash
# SHORT POLL (the ride's one-liner): the masters' markers + the widen watcher,
# compact, one line per poll. usage: bash poll-once.sh
RV=${RV:-/home/erikg/hprcv2/wgs35x-fleet/runtime-verification}
cd "$RV" || exit 1
echo "$(date -u +%FT%TZ) split=$(wc -l < prep/split-done.txt) acc=$(wc -l < prep/accepted.txt) sub=$(ls prep/*/substrate.done 2>/dev/null | wc -l) cut=$(ls map/ 2>/dev/null | wc -l) call=$(ls receipts/ 2>/dev/null | wc -l) failcall=$(ls prep/*/call.exit 2>/dev/null | wc -l) df=$(df -BG / | awk 'NR==2{gsub("G","",$4); print $4}')G mem=$(free -g | awk 'NR==2{print $7}')Gib rvsize=$(du -sBG --apparent-size . 2>/dev/null | awk '{gsub("G","",$1); print $1}' 2>/dev/null)G qdone=$([ -f query/queries.done ] && echo Y || echo n) fsplit=$([ -f prep/final-split.done ] && echo Y || echo n) widen=$([ -f widen-done ] && echo Y || echo n) poison=$(cat cuts/poison.tsv 2>/dev/null | grep -c . ) lastkill=$(tail -1 cuts/batch-failures.log 2>/dev/null | cut -c1-90)"
