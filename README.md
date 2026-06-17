# Amyloid Mutant Proteome Nextflow Pipeline

Reproducible Nextflow DSL2 workflow for SRA-to-mutant-proteome processing, amyloidogenicity scoring, and protein sequence feature calculation.

The repository is Nextflow-first. The workflow now includes a native Nextflow full mode for SRA download, QC, expression quantification, alignment, variant calling, proteome generation, AmyloGram-Py scoring, optional AMYPred-FRL scoring, and lightweight feature calculation. Production scientific components such as R `AmyloGram` and full `ProteinFeatures` are retained for Nextflow integration.

## What It Does

In full mode, each sample starts from an SRA accession:

```text
SRA accession
```

In prepared mode, each sample starts from prepared upstream files:

```text
abundance.tsv
variants_filtered.vcf
reference genome FASTA
reference GTF
AmyloGram-Py six-mer table
```

The Nextflow workflow performs:

1. FASTQ download with SRA Toolkit.
2. FastQC read quality reports.
3. Kallisto transcript quantification.
4. STAR alignment.
5. Picard duplicate marking.
6. GATK RNA-seq variant calling and filtration.
7. Expression and CDS filtering.
8. SNP/nonsense proteome reconstruction.
9. Frameshift protein reconstruction.
10. Combined mutant proteome generation and verification.
11. AmyloGram-Py prediction.
12. Consensus-compatible amyloid table normalization.
13. Lightweight protein feature calculation.

`AMYPred-FRL` and full `ProteinFeatures` are kept under `amyloid_predictors/` and are not considered legacy. AMYPred-FRL is available in Nextflow prepared/full branches when its `data/` and `model/` runtime assets are present. Full `ProteinFeatures` is retained for production integration; the current default branch emits light feature outputs.

## Repository Layout

```text
.
├── nextflow/                  # DSL2 workflow and modules
├── nextflow.config            # Profiles and default parameters
├── proteome_pipeline/         # Python library used by Nextflow proteome modules
├── amyloid_predictors/        # AMYPred-FRL and ProteinFeatures production components
├── amylogram-py/              # AmyloGram-Py package and fixtures
├── docker/                    # Nextflow Docker images
├── tests_nextflow/            # Nextflow contract, smoke, and fixture tests
├── docs/                      # Pipeline documentation and figures
├── article_reproduction/      # Article reproduction notebooks and analysis outputs
└── REPRODUCIBILITY.md         # Nextflow reproducibility guide
```

## Quick Start

Run the minimal prepared-input fixture:

```bash
bash tests_nextflow/run_prepared_smoke.sh
```

Run the amyloid-only fixture:

```bash
bash tests_nextflow/run_amyloid_smoke.sh
```

Run Nextflow directly:

```bash
nextflow run nextflow/main.nf \
  -profile test,docker \
  --mode prepared \
  --samples tests_nextflow/fixtures/prepared_minimal/samples.prepared.csv \
  --ref_genome tests_nextflow/fixtures/prepared_minimal/genome.fa \
  --ref_gtf tests_nextflow/fixtures/prepared_minimal/annotation.gtf \
  --sixmer_table amylogram-py/tests/fixtures/amylogram_sixmer_probabilities.bin \
  --outdir results_nextflow
```

For AWS/reference-compatible proteome reproduction, add the `aws_reference`
profile. This sets `--min_cds_bp 1` and preserves ultra-short CDS-derived
peptides that are present in the reference proteome:

```bash
nextflow run nextflow/main.nf \
  -profile docker,aws_reference \
  --mode prepared \
  --samples samples.prepared.csv \
  --ref_genome references/GRCh38.primary_assembly.genome.fa \
  --ref_gtf references/gencode.v46.primary_assembly.annotation.gtf \
  --sixmer_table amylogram-py/tests/fixtures/amylogram_sixmer_probabilities.bin \
  --outdir results_nextflow_reference
```

The default `--min_cds_bp 10` intentionally filters ultra-short CDS records.
Use `aws_reference` when exact comparison against the AWS/reference
`combine_proteome.fasta` is required.

Run prepared mode with AMYPred-FRL enabled:

```bash
nextflow run nextflow/main.nf \
  -profile docker,aws_reference \
  --mode prepared \
  --samples samples.prepared.csv \
  --ref_genome references/GRCh38.primary_assembly.genome.fa \
  --ref_gtf references/gencode.v46.primary_assembly.annotation.gtf \
  --sixmer_table amylogram-py/tests/fixtures/amylogram_sixmer_probabilities.bin \
  --run_amypred true \
  --outdir results_nextflow_amypred
```

This requires a complete `amyloid_predictors/AMYPred-FRL` directory containing
`predict.py`, `data/PAAC.txt`, `data/TR_P_132.fasta`, `data/TR_N_305.fasta`,
and `model/pima.pickle_model_svm_PF.dat`.

Add old R AmyloGram and full R ProteinFeatures when production-style outputs
are needed:

```bash
nextflow run nextflow/main.nf \
  -profile docker,aws_reference \
  --mode prepared \
  --samples samples.prepared.csv \
  --ref_genome references/GRCh38.primary_assembly.genome.fa \
  --ref_gtf references/gencode.v46.primary_assembly.annotation.gtf \
  --sixmer_table amylogram-py/tests/fixtures/amylogram_sixmer_probabilities.bin \
  --run_amypred true \
  --run_amylogram_r true \
  --run_protein_features true \
  --outdir results_nextflow_production
```

Run full SRA-to-results mode:

```bash
nextflow run nextflow/main.nf \
  -profile docker \
  --mode full \
  --samples samples.txt \
  --ref_dir references \
  --ref_genome_name GRCh38.primary_assembly.genome.fa \
  --ref_gtf references/gencode.v46.primary_assembly.annotation.gtf \
  --star_index references/star_index_full \
  --kallisto_index references/transcriptome.idx \
  --sixmer_table amylogram-py/tests/fixtures/amylogram_sixmer_probabilities.bin \
  --outdir results_nextflow
```

## Bash/S3 Compatibility Entrypoints

The repository also retains the original server-oriented bash/S3 entrypoints
for compatibility with existing production runs:

```bash
export AWS_PROFILE=user-sandbox
export AWS_REGION=us-east-1
export S3_BUCKET=s3://your-bucket/your-prefix
export CLEANUP_AFTER_UPLOAD=true
```

Run the complete bash workflow:

```bash
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
SOURCE_S3=s3://prepared-ngsdata \
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
  amypred_frl_prediction.csv          # when --run_amypred true
  amylogram_py_prediction.csv
  amylogram_py_report.json
  amylogram_py_report.md
  amylogram_py_skipped.tsv
  amylogram_py_top_hits.tsv
  amyloid_combined_predictions.csv
  amyloid_predictors_status.tsv
  amyloid_predictors_summary.tsv

results_protein_features/
  protein_features_light.csv
  protein_features_light_status.tsv
  protein_features_light_summary.tsv
  protein_features.csv                # when full ProteinFeatures integration is enabled
```

## Tests

Run the lightweight contract tests:

```bash
python3 -m unittest tests_nextflow.test_contract
python3 -m py_compile tests_nextflow/validate_amyloid_outputs.py tests_nextflow/test_amyloid_outputs.py
bash -n tests_nextflow/run_prepared_smoke.sh
bash -n tests_nextflow/run_amyloid_smoke.sh
bash -n tests_nextflow/test_resume.sh
```

Run executable smoke tests when Docker and Nextflow are available:

```bash
bash tests_nextflow/run_prepared_smoke.sh
bash tests_nextflow/run_amyloid_smoke.sh
```

## Safety

This repository intentionally does **not** include:

- AWS credentials;
- SSH keys;
- `.aws/`;
- FASTQ/BAM/reference production datasets;
- generated results.
