#!/bin/bash
# Step 2: Align reads to the reference genome using BWA

# Directories
READ_DIR="../data/raw"
OUT_DIR="../data/processed"
REF_GENOME="../reference/genome.fa"  # Update with your actual path
THREADS=4

mkdir -p $OUT_DIR

# Index the reference genome if not already done
if [ ! -e "${REF_GENOME}.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index $REF_GENOME
fi

# Align all paired-end FASTQ files
for file in $READ_DIR/*_R1.fastq.gz
do
    base=$(basename "$file" _R1.fastq.gz)

    bwa mem -t $THREADS $REF_GENOME \
        $READ_DIR/${base}_R1.fastq.gz \
        $READ_DIR/${base}_R2.fastq.gz > $OUT_DIR/${base}.sam

    echo "Aligned: $base"
done

echo " All samples aligned. SAM files saved to $OUT_DIR"
