#!/bin/bash
# Step 1: Quality control using FastQC and MultiQC

RAW_DIR="../data/raw"
OUT_DIR="../results/QC_reports"

mkdir -p $OUT_DIR

# Run FastQC
fastqc $RAW_DIR/*.fastq.gz -o $OUT_DIR

# Run MultiQC
multiqc $OUT_DIR -o $OUT_DIR

echo " Quality control completed. Reports saved in $OUT_DIR"
