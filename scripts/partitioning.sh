#!/bin/bash

PAF=$1
FASTA_FAI=$2
AVG_WINDOW_SIZE=$3
SAMPLE=$4

grep $SAMPLE $FASTA_FAI -w -m 1 | awk -v OFS='\t' '{print($1,"0",$2)}' > $SAMPLE.bed
bedtools makewindows -b $SAMPLE.bed -w $AVG_WINDOW_SIZE > windows.bed

MASK_BED=mask.bed
touch $MASK_BED

MISSING_BED=missing.bed
cut -f 1,2 $PAF | sort | uniq | awk -v OFS='\t' '{print($1,"0",$2,".",".","+")}' > $MISSING_BED

num=0
while [ -s windows.bed ]; do
    echo "Processing new window set..."

    cat windows.bed | while read REGION; do
        REGION2=$(echo $REGION | sed 's/\t/:/' | sed 's/\t/-/')
        echo Partition $num - $REGION

        # Get partition
        impg -p $PAF -r $REGION2 -x | bedtools sort | bedtools merge -s -c 4,5,6 -o distinct > partition$num.tmp.bed

        # Apply mask
        bedtools subtract -a partition$num.tmp.bed -b $MASK_BED -s > partition$num.bed

        # Update mask
        cat partition$num.bed $MASK_BED | bedtools sort > $num.bed
        bedtools merge -i $num.bed -s -c 4,5,6 -o distinct > $num.mask.bed
        cp $num.mask.bed $MASK_BED

        # Update missing
        bedtools subtract -a $MISSING_BED -b partition$num.bed -s > $num.missing.bed
        cp $num.missing.bed $MISSING_BED

        # Cleanup
        rm partition$num.tmp.bed $num.bed $num.merge.bed $num.missing.bed

        num=$((num + 1))
    done

    # Check if there are any missing regions not covered by mask
    bedtools subtract -a $MISSING_BED -b $MASK_BED -s > $num.remaining.bed
    
    if [ ! -s $num.remaining.bed ]; then
        # If no remaining regions, we're done
        ###rm $num.remaining.bed
        break
    fi

    # Take longest missing range and create new windows
    awk -v OFS='\t' '{print($0,$3-$2)}' $num.remaining.bed | sort -k 7,7nr | cut -f 1-6 | head -n 1 > $num.new.bed
    bedtools makewindows -b $num.new.bed -w $AVG_WINDOW_SIZE > windows.bed
    
    # Cleanup
    rm $num.remaining.bed $num.new.bed
done

rm windows.bed
