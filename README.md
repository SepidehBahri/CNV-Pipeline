# CNV-Pipeline

A pipeline for detecting Copy Number Variations (CNVs) from NGS data using BWA, SAMtools, and CNVkit.

## Pipeline Overview

This pipeline processes whole exome or targeted sequencing data to detect CNVs. It includes quality control, alignment, BAM processing, CNV calling, and visualization steps.

### Steps

| Step | Task                                  | Script                      |
|------|----------------------------------------|-----------------------------|
| 1    | Quality control (FastQC, MultiQC)      | scripts/01_quality_control.sh |
| 2    | Alignment to reference genome (BWA)    | scripts/02_alignment.sh     |
| 3    | BAM sorting and indexing (SAMtools)    | scripts/03_bam_processing.sh |
| 4    | CNV calling (CNVkit)                   | scripts/04_cnv_calling.sh   |
| 5    | CNV visualization and reporting (R)    | scripts/05_cnv_plotting.R   |

## Folder Structure

```
CNV-Pipeline/
├── data/
│   ├── raw/                # Raw FASTQ files
│   └── processed/          # Aligned, sorted, indexed BAM files
├── reference/              # Genome FASTA, BED targets, etc.
├── scripts/                # All processing scripts
├── config/                 # Sample metadata
├── results/
│   ├── cnv_reports/        # CNVkit outputs (.cns, .cnr, .cnn)
│   └── cnv_plots/          # Plots generated in R
├── env/                    # Conda environment file
├── logs/                   # Log files for each stage
├── .gitignore
└── README.md
```

## Requirements

Create the conda environment:

```bash
conda create -n cnv_env -c bioconda -c conda-forge \
  fastqc multiqc bwa samtools cnvkit r-base
```

Then activate it:

```bash
conda activate cnv_env
```

## Usage

Run scripts step by step:

```bash
bash scripts/01_quality_control.sh
bash scripts/02_alignment.sh
bash scripts/03_bam_processing.sh
bash scripts/04_cnv_calling.sh
Rscript scripts/05_cnv_plotting.R
```

## Notes

- Reference files (`genome.fa`, `targets.bed`, `antitargets.bed`) must be placed in the `reference/` folder.
- The `sample_info.csv` file in `config/` should list sample names and conditions.
- CNVkit reference `.cnn` will be created from normals if available, or must be provided manually.

## Output

- CNVkit result files: `.cns`, `.cnr`, `.cnn`
- CNV plots: per-sample PDF visualizations of log2 CNV segments
