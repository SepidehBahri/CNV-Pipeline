#!/bin/bash
# Step 4: CNV calling using CNVkit

# Directories
INPUT_DIR="../data/processed"
OUTPUT_DIR="../results/cnv_reports"
REF_DIR="../reference"
THREADS=4

mkdir -p $OUTPUT_DIR

# Path to target and antitarget BED files (designed for WES)
TARGETS=$REF_DIR/targets.bed
ANTITARGETS=$REF_DIR/antitargets.bed
REF_GENOME=$REF_DIR/genome.fa

# Optional pooled normal BAMs to build reference (skip if you already have .cnn)
NORMAL_BAMS=$(ls $INPUT_DIR/*_normal*_sorted.bam 2>/dev/null)

# Step 1: Build reference if not available
if [ ! -f "$REF_DIR/cnv_reference.cnn" ] && [ ! -z "$NORMAL_BAMS" ]; then
    cnvkit.py batch $NORMAL_BAMS \
        --normal \
        --targets $TARGETS \
        --antitargets $ANTITARGETS \
        --fasta $REF_GENOME \
        --output-reference $REF_DIR/cnv_reference.cnn \
        --output-dir $OUTPUT_DIR \
        --processes $THREADS
fi

# Step 2: Run CNVkit on each sample
for bamfile in $INPUT_DIR/*_sorted.bam
do
    sample=$(basename "$bamfile" _sorted.bam)

    cnvkit.py batch $bamfile \
        --targets $TARGETS \
        --antitargets $ANTITARGETS \
        --fasta $REF_GENOME \
        --reference $REF_DIR/cnv_reference.cnn \
        --output-dir $OUTPUT_DIR \
        --output-reference $REF_DIR/cnv_reference.cnn \
        --processes $THREADS

    echo "CNV analysis complete for: $sample"
done
