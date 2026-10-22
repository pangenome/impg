#!/usr/bin/env bash
# THE GENOME TRANCHE - stage 2: the reads fetch (the committed fetch
# rule, run-wave.sh VERBATIM: primaries in window +/-50kb, flat fastq)
# for every tranche territory, marker-gated.
# usage: run-reads.sh   (env: WT, POOL_W=8)
set -uo pipefail
WT=${WT:-/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1}
POOL_W=${POOL_W:-8}
BAM=/home/erikg/hprcv2/HG002.GRCh38.35x.bam
mkdir -p "$WT/reads" "$WT/logs"
export WT BAM
grep -v '^#' "$WT/tranche.tsv" | cut -f1 | \
  xargs -P "$POOL_W" -I{} bash -c '
    T={}
    [ -f "$WT/prep/$T/reads.done" ] && exit 0
    CHROM=$(awk -F"\t" -v t="$T" "\$1==t {print \$3}" "$WT/tranche.tsv")
    S=$(awk -F"\t" -v t="$T" "\$1==t {print \$4}" "$WT/tranche.tsv")
    E=$(awk -F"\t" -v t="$T" "\$1==t {print \$5}" "$WT/tranche.tsv")
    start=$(date +%s)
    samtools view -u -F 0x904 "$BAM" "$CHROM:$((S-50000))-$((E+50000))" \
      | samtools fastq - 2> "$WT/logs/reads-$T-fq.err" \
      | gzip -1 > "$WT/reads/$T.fq.gz"
    st=${PIPESTATUS[0]}
    echo "$st" > "$WT/prep/$T/reads.exit"
    echo "$(( $(date +%s) - start ))" > "$WT/prep/$T/reads.wall"
    [ "$st" -eq 0 ] && touch "$WT/prep/$T/reads.done"
  '
echo "$(date -u +%FT%TZ) reads fetch complete: $(ls $WT/reads/*.fq.gz 2>/dev/null | wc -l) of $(grep -v '^#' $WT/tranche.tsv | wc -l)"
