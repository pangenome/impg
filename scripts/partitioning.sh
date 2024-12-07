#!/bin/bash

# Exit on error, undefined variables, and pipe failures
set -euo pipefail
#set -x # for debugging

# Check if all required arguments are provided
if [ $# -ne 4 ]; then
    echo "Usage: $0 <paf_file> <fasta_index> <window_size> <sample_name>"
    exit 1
fi

# Input parameters
PAF=$1
FASTA_FAI=$2
AVG_WINDOW_SIZE=$3
SAMPLE=$4

# Prepare initial windows
WINDOWS_BED=windows.bed
grep $SAMPLE "$FASTA_FAI" -w | awk -v OFS='\t' '{print($1,"0",$2)}' > "$SAMPLE.bed"
bedtools makewindows -b "$SAMPLE.bed" -w "$AVG_WINDOW_SIZE" > "$WINDOWS_BED"

# Prepare mask file
MASK_BED=mask.bed
cat /dev/null > "$MASK_BED"

# Create initial missing regions file
MISSING_BED=missing.bed
cut -f 1,2 "$PAF" | sort | uniq | awk -v OFS='\t' '{print($1,"0",$2)}' > "$MISSING_BED"

# Initialize partition counter
num=0

# Process windows until no missing regions remain
while [ -s "$WINDOWS_BED" ]; do
    echo "Processing new window set"

    while IFS=$'\t' read -r chrom start end; do
        # Get partition, sort, and merge to avoid duplicates
        REGION_FORMATTED="${chrom}:${start}-${end}"

        echo "-- Querying region $REGION_FORMATTED"
        impg -p "$PAF" -r "$REGION_FORMATTED" -x | bedtools sort | bedtools merge -d 10000 > "partition$num.tmp.bed"

        # Apply mask
        bedtools subtract -a "partition$num.tmp.bed" -b "$MASK_BED" > "partition$num.bed"

        # Check if the partition is not empty
        if [ -s "partition$num.bed" ]; then
            echo "-- Processing partition $num"

            # Update masked regions
            cat "partition$num.bed" "$MASK_BED" | bedtools sort > "$num.mask.bed"
            bedtools merge -i "$num.mask.bed" > "$MASK_BED"

            # Update missing regions
            bedtools subtract -a "$MISSING_BED" -b "partition$num.bed" > "$num.missing.bed"
            cp "$num.missing.bed" "$MISSING_BED"

            num=$((num + 1))
        fi
    done < $WINDOWS_BED

    # Check if there are any missing regions not covered by mask
    bedtools subtract -a "$MISSING_BED" -b "$MASK_BED" -s > "$num.remaining.bed"
    if [ ! -s "$num.remaining.bed" ]; then
        # If no remaining regions, we're done
        break
    fi

    # Create new windows from longest remaining region
    awk -v OFS='\t' 'BEGIN{max_len=0} {len=$3-$2; if(len>max_len){max_len=len; line=$0}} END{print line}' "$num.remaining.bed" > "$num.new.bed"

    bedtools makewindows -b "$num.new.bed" -w "$AVG_WINDOW_SIZE" > "$WINDOWS_BED"
done

# Cleanup
rm -f "$SAMPLE.bed" "$WINDOWS_BED" "$MASK_BED" "$MISSING_BED" partition[0-9]*.tmp.bed [0-9]*.{mask,missing,remaining,new}.bed
