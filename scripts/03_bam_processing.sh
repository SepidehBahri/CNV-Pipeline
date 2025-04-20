#!/bin/bash
# Step 3: Convert SAM to BAM, sort, and index using SAMtools

INPUT_DIR="../data/processed"
OUTPUT_DIR="../data/processed"
THREADS=4

for samfile in $INPUT_DIR/*.sam
do
    base=$(basename "$samfile" .sam)

    # Convert SAM to BAM
    samtools view -@ $THREADS -bS $samfile > $OUTPUT_DIR/${base}.bam

    # Sort BAM
    samtools sort -@ $THREADS -o $OUTPUT_DIR/${base}_sorted.bam $OUTPUT_DIR/${base}.bam

    # Index BAM
    samtools index $OUTPUT_DIR/${base}_sorted.bam

    # Optionally remove intermediate files
    rm $OUTPUT_DIR/${base}.sam $OUTPUT_DIR/${base}.bam

    echo "Processed: $base"
done

echo "All BAM files sorted and indexed."
