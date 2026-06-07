# Amyloid Mutant Proteome Pipeline

Reproducible single-sample bioinformatics workflow for generating a cleaned mutant proteome from public SRA RNA-seq data, scoring amyloidogenicity, and calculating protein physicochemical properties.

- **AMYPred-FRL** for protein-level amyloid probability.
- **AmyloGram** in an optimized Docker runner with sequence deduplication, chunk checkpoints, and S3 upload.
- **ProteinFeatures** R step for protein length, Kyte-Doolittle hydrophobicity, charge, aromaticity, alpha/beta propensity, chameleon score, and pI.

The pipeline was developed for translational bioinformatics and immuno-oncology use cases where each sample is processed independently on a server.

## What It Does

For each SRA accession, the workflow performs:

1. FASTQ download from a public repository/SRA.
2. FastQC read quality control.
3. Kallisto transcript quantification.
4. STAR alignment and GATK variant calling.
5. Mutant proteome generation:
   - expression and mitochondrial filtering;
   - protection against ultra-short CDS records that can crash `gffread`;
   - nonsense mutation truncation at the first stop;
   - frameshift protein reconstruction by surgical indel insertion/deletion;
   - per-gene sequence deduplication;
   - final `combine_proteome.fasta`.
6. Amyloidogenicity prediction:
   - AMYPred-FRL;
   - AmyloGram;
   - merged consensus table.
7. Protein physicochemical feature calculation from `combine_proteome.fasta`.
8. Upload of results, logs, manifests, QC, expression, BAM and variant artifacts to S3.

## Repository Layout

```text
.
├── pipeline_scripts/          # Main bash pipeline steps
├── proteome_pipeline/         # Python library for FASTA/GTF/VCF/proteome logic
├── docker/amylogram/          # Reproducible AmyloGram Docker image
├── tests/                     # Proteome/static tests
├── tests_amyloid/             # Amyloid predictor contract tests
├── run_full_pipeline.sh       # One-command end-to-end runner
├── run_complete_processing_pipeline.sh
│                              # Explicit complete workflow runner including protein features
├── run_prepared_proteomes_from_bioinfo.sh
│                              # Proteome-only batch from prepared S3 inputs
├── run_amyloid_predictors.sh  # AMYPred-FRL + AmyloGram + ProteinFeatures runner
├── predict_amylogram_fast.R   # Optimized AmyloGram batch runner
└── REPRODUCIBILITY.md         # Full reproducibility guide
```

## Quick Start

Prepare a sample list:

```bash
cp samples.txt.example samples.txt
```

Edit `samples.txt`:

```text
SRR32060234
```

Run the complete workflow:

```bash
export AWS_PROFILE=codex-sandbox
export AWS_REGION=us-east-1
export S3_BUCKET=s3://your-bucket/your-prefix
export CLEANUP_AFTER_UPLOAD=true

./run_full_pipeline.sh \
  --samples samples.txt \
  --s3-bucket "$S3_BUCKET"
```

Equivalent explicit complete-processing entrypoint:

```bash
./run_complete_processing_pipeline.sh \
  --samples samples.txt \
  --s3-bucket "$S3_BUCKET"
```

Run amyloid predictors only, after `combine_proteome.fasta` already exists:

```bash
./run_amyloid_predictors.sh \
  --sample SRR32060234 \
  --output-s3 s3://your-bucket/your-prefix/SRR32060234/results_amyloid
```

Force amyloid recomputation:

```bash
./run_amyloid_predictors.sh --sample SRR32060234 --force
```

Run proteome generation from already prepared S3 inputs:

```bash
SOURCE_S3=s3://bioinfo-data-amylice-2026 \
DEST_S3=s3://your-bucket/prepared-bioinfo-proteome \
./run_prepared_proteomes_from_bioinfo.sh
```

## Main Outputs

For each sample:

```text
results_proteins/
  combine_proteome.fasta
  protein.fasta
  frameshift_unique.fasta
  clean_ids.txt
  nonsense_candidates.txt
  verification_report.tsv
  manifest.txt

results_amyloid/
  amypred_frl_prediction.csv
  amylogram_prediction_fast.csv
  protein_features.csv
  amyloid_combined_predictions.csv
  amyloid_predictors_status.tsv
  amyloid_predictors_summary.tsv

results_protein_features/
  protein_features.csv
  protein_features_status.tsv
  protein_features_summary.tsv
  logs/

results_qc/
results_expression/
results_star/
results_gatk/
logs/
```

## Tests

Run static and unit tests:

```bash
python3 -m pytest tests tests_amyloid
bash tests/run_static_tests.sh
bash tests_amyloid/run_amyloid_contract_tests.sh
bash tests_amyloid/test_run_amyloid_predictors_static.sh ./run_amyloid_predictors.sh
```

## Safety

This repository intentionally does **not** include:

- AWS credentials;
- SSH keys;
- `.aws/`;
- FASTQ/BAM/VCF/reference datasets;
- generated results.

Use S3 buckets and AWS profiles through environment variables.
