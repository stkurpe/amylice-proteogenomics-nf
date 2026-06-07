# Reproducibility Guide

This document describes how to reproduce a single-sample run from public SRA accession to mutant proteome and amyloidogenicity predictions.

## 1. Compute Environment

Recommended server:

- Ubuntu 24.04 LTS or similar Linux distribution.
- At least 60 GB disk for one sample at a time.
- 8 CPU threads.
- 32 GB RAM recommended for comfortable STAR/GATK/AmyloGram runs.
- Docker available to the workflow user.
- AWS CLI configured with read access to reference data and write access to the output bucket.

Required command-line tools:

```text
aws
docker
fastqc
fastp
kallisto
STAR
samtools
bcftools
gffread
gatk
python3
Rscript, via Docker image for AmyloGram
```

The workflow is designed for one sample per run. This keeps disk pressure predictable and makes cleanup after successful S3 upload safe.

## 2. Reference Data

The pipeline expects reference files in `REF_DIR`, defaulting to:

```bash
/home/codex/references
```

Default reference filenames are configured in `pipeline_scripts/config.sh`:

```bash
GRCh38.primary_assembly.genome.fa
gencode.v46.primary_assembly.annotation.gtf
common_all_20180418.vcf
gencode.v46.pc_transcripts.fa
```

If references are stored in S3, configure `pipeline_scripts/config_s3.sh` or set environment variables before launch.

## 3. AWS Configuration

Configure an AWS profile outside the repository:

```bash
aws configure --profile codex-sandbox
```

Never commit `.aws`, access key CSV files, tokens, or credentials.

For a custom output bucket:

```bash
export AWS_PROFILE=codex-sandbox
export AWS_REGION=us-east-1
export S3_BUCKET=s3://your-bucket/public-sra-test
```

The final output path for sample `SRR32060234` will be:

```text
s3://your-bucket/public-sra-test/SRR32060234/
```

## 4. Build AmyloGram Docker Image

Build once on the compute server:

```bash
docker build \
  -t codex/amylogram:1.1-r4.3.3 \
  docker/amylogram
```

The optimized AmyloGram runner uses this image and does not install R packages on every run.

## 5. Prepare Samples

Create `samples.txt`:

```text
SRR32060234
```

Only one sample is recommended per run. Multiple samples are supported as separate sequential jobs in the sample file.

## 6. Run End-to-End Workflow

```bash
export AWS_PROFILE=codex-sandbox
export AWS_REGION=us-east-1
export THREADS=8
export CLEANUP_AFTER_UPLOAD=true
export PREPARE_REFS=true
export AMYLOGRAM_CORES=8
export AMYLOGRAM_CHUNK_SIZE=50

./run_full_pipeline.sh \
  --samples samples.txt \
  --s3-bucket "$S3_BUCKET"
```

This runs:

```text
SRA download
  -> FastQC
  -> Kallisto
  -> STAR/GATK
  -> mutant proteome generation
  -> AMYPred-FRL + AmyloGram
  -> S3 upload
```

## 7. Run Proteome Only

```bash
./run_full_pipeline.sh \
  --samples samples.txt \
  --s3-bucket "$S3_BUCKET" \
  --skip-amyloid
```

## 8. Run Amyloid Predictors Only

Use this after `combine_proteome.fasta` exists locally under:

```text
/home/codex/amyloid_runs/<SAMPLE_ID>/input/combine_proteome.fasta
```

or after it is available in:

```text
s3://bucket/prefix/<SAMPLE_ID>/results_proteins/combine_proteome.fasta
```

Command:

```bash
./run_amyloid_predictors.sh \
  --sample SRR32060234 \
  --input-s3 "$S3_BUCKET/SRR32060234/results_proteins/combine_proteome.fasta" \
  --output-s3 "$S3_BUCKET/SRR32060234/results_amyloid"
```

## 9. Integrity Checks

Proteome generation writes:

```text
verification_report.tsv
manifest.txt
local_proteome.events.tsv
local_proteome.log
```

Amyloid prediction writes:

```text
amyloid_predictors_status.tsv
amyloid_predictors_summary.tsv
amyloid_combined_predictions.csv
```

Inspect status files after every run:

```bash
cat /home/codex/amyloid_runs/SRR32060234/results/amyloid_predictors_status.tsv
cat /home/codex/amyloid_runs/SRR32060234/results/amyloid_predictors_summary.tsv
```

## 10. Expected SRR32060234 Amyloid Summary

From the validated run:

```text
input_records:   28327
AMYPred rows:    27978
AmyloGram rows:  28312
combined rows:   28312
```

AmyloGram optimized run:

```text
valid records:      28312
unique sequences:   14199
duplicates saved:   14113
chunks:             284
cores:              8
runtime:            about 52 minutes
errors:             0
```

## 11. Cleanup Policy

Set:

```bash
export CLEANUP_AFTER_UPLOAD=true
```

The pipeline is intended to remove local sample work only after successful S3 upload. Do not manually delete intermediate files unless the uploaded artifacts have been verified.

## 12. Re-running

For amyloid predictors:

```bash
./run_amyloid_predictors.sh --sample SRR32060234
```

If outputs already exist, heavy steps are skipped. To recompute:

```bash
./run_amyloid_predictors.sh --sample SRR32060234 --force
```

## 13. Run Proteome Generation From Prepared S3 Inputs

If expression and variant results already exist in S3, the proteome step can be run without redownloading FASTQ or rerunning STAR/GATK.

Expected input structure:

```text
s3://source-bucket/<SAMPLE_ID>/results_expression/abundance.tsv
s3://source-bucket/<SAMPLE_ID>/results_gatk/variants_filtered.vcf
```

or:

```text
s3://source-bucket/<SAMPLE_ID>/results_gatk/variants_filtered.vcf.gz
```

The sample list is read from:

```text
s3://source-bucket/samples.txt
```

Run:

```bash
export SOURCE_S3=s3://bioinfo-data-amylice-2026
export DEST_S3=s3://codex-test-ngsdata-calculations/prepared-bioinfo-proteome
export AWS_PROFILE=codex-sandbox
export AWS_REGION=us-east-1

./run_prepared_proteomes_from_bioinfo.sh
```

This script processes one sample at a time, uploads proteome outputs to `DEST_S3`, and cleans the local sample directory only after the upload success marker exists.
