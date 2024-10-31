#!/bin/bash

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

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
grep $SAMPLE $FASTA_FAI -w -m 1 | awk -v OFS='\t' '{print($1,"0",$2)}' > $SAMPLE.bed
bedtools makewindows -b $SAMPLE.bed -w $AVG_WINDOW_SIZE > $WINDOWS_BED

# Prepare mask file
MASK_BED=mask.bed
cat /dev/null > $MASK_BED

# Create initial missing regions file
MISSING_BED=missing.bed
cut -f 1,2 $PAF | sort | uniq | awk -v OFS='\t' '{print($1,"0",$2,".",".","+")}' > $MISSING_BED

# Initialize partition counter
num=0

# Process windows until no missing regions remain
while [ -s $WINDOWS_BED ]; do
    echo "Processing new window set..."

    while IFS=$'\t' read -r chrom start end; do
        # Get partition, sort, and merge to avoid duplicates
        REGION_FORMATTED="${chrom}:${start}-${end}"

        echo "Processing partition $num, region $REGION_FORMATTED..."
        impg -p $PAF -r $REGION_FORMATTED -x | bedtools sort | bedtools merge -s -c 4,5,6 -o distinct -d 5000 > partition$num.tmp.bed

        # Apply mask
        bedtools subtract -a partition$num.tmp.bed -b $MASK_BED -s > partition$num.bed

        # Check if the partition is not empty
        if [ -s partition$num.bed ]; then

            # Update masked regions
            cat partition$num.bed $MASK_BED | bedtools sort > $num.bed
            bedtools merge -i $num.bed -s -c 4,5,6 -o distinct > $num.mask.bed
            cp $num.mask.bed $MASK_BED

            # Update missing regions
            bedtools subtract -a $MISSING_BED -b partition$num.bed -s > $num.missing.bed
            cp $num.missing.bed $MISSING_BED

            # Cleanup
            rm partition$num.tmp.bed $num.bed $num.mask.bed $num.missing.bed

            num=$((num + 1))
        else
            rm partition$num.tmp.bed partition$num.bed
        fi
    done < $WINDOWS_BED

    # Check if there are any missing regions not covered by mask
    bedtools subtract -a $MISSING_BED -b $MASK_BED -s > $num.remaining.bed
    if [ ! -s $num.remaining.bed ]; then
        # If no remaining regions, we're done
        rm $num.remaining.bed
        break
    fi

    # Create new windows from longest remaining region
    awk -v OFS='\t' '{print($0,$3-$2)}' $num.remaining.bed | sort -k 7,7nr | cut -f 1-6 | head -n 1 > $num.new.bed
    bedtools makewindows -b $num.new.bed -w $AVG_WINDOW_SIZE > $WINDOWS_BED
    
    # Cleanup
    rm $num.remaining.bed $num.new.bed
done

rm $SAMPLE.bed $WINDOWS_BED $MASK_BED $MISSING_BED
