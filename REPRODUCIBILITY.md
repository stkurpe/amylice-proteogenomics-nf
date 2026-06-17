# Reproducibility

This repository is a Nextflow workflow with full SRA-to-results mode and prepared-input mode.

## Requirements

- Nextflow
- Docker
- Python 3
- Java runtime compatible with Nextflow

The smoke tests build local Docker images for:

```text
amyloid-proteome-nextflow:local
amylogram-py-nextflow:local
```

## Prepared-Input Mode

The workflow starts from files that are already available locally:

```text
abundance.tsv
variants_filtered.vcf
reference genome FASTA
reference GTF
AmyloGram-Py six-mer table
```

Minimal fixture run:

```bash
bash tests_nextflow/run_prepared_smoke.sh
```

Direct Nextflow invocation:

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

Reference-compatible prepared invocation:

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

The `aws_reference` profile sets `--min_cds_bp 1`. This retains ultra-short
CDS-derived peptides and is intended for exact comparison against the AWS
reference `combine_proteome.fasta`. The default `--min_cds_bp 10` intentionally
filters those records.

AWS reference regression:

```bash
AWS_PREPARED_ABUNDANCE=/path/to/abundance.tsv \
AWS_PREPARED_VCF=/path/to/variants_filtered.vcf \
AWS_REFERENCE_FASTA=/path/to/reference/combine_proteome.fasta \
REF_GENOME=/path/to/GRCh38.primary_assembly.genome.fa \
REF_GTF=/path/to/gencode.v46.primary_assembly.annotation.gtf \
bash tests_nextflow/run_aws_prepared_reference.sh
```

For the validated `SRR32060234` reference run, the expected
`combine_proteome.fasta` contains 28,343 records and has sorted ID/sequence MD5
`de4d381238adb9c1e72e714e1e72abb1`.

Set `RUN_AMYPRED=1` on the same script to exercise the optional AMYPred-FRL
prepared branch:

```bash
RUN_AMYPRED=1 \
AWS_PREPARED_ABUNDANCE=/path/to/abundance.tsv \
AWS_PREPARED_VCF=/path/to/variants_filtered.vcf \
AWS_REFERENCE_FASTA=/path/to/reference/combine_proteome.fasta \
REF_GENOME=/path/to/GRCh38.primary_assembly.genome.fa \
REF_GTF=/path/to/gencode.v46.primary_assembly.annotation.gtf \
bash tests_nextflow/run_aws_prepared_reference.sh
```

This requires the complete AMYPred-FRL runtime directory under
`amyloid_predictors/AMYPred-FRL`, including `data/` training FASTA/PAAC files
and `model/pima.pickle_model_svm_PF.dat`.

## Amyloid-Only Mode

```bash
bash tests_nextflow/run_amyloid_smoke.sh
```

## Full SRA-To-Results Mode

Full mode requires SRA access plus local reference assets:

```text
references/GRCh38.primary_assembly.genome.fa
references/GRCh38.primary_assembly.genome.fa.fai
references/GRCh38.primary_assembly.genome.dict
references/gencode.v46.primary_assembly.annotation.gtf
references/star_index_full/
references/transcriptome.idx
```

Run:

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

## Resume Check

```bash
bash tests_nextflow/test_resume.sh
```

## Expected Outputs

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
  amylogram_py_prediction.csv
  amyloid_combined_predictions.csv
  amyloid_predictors_status.tsv
  amyloid_predictors_summary.tsv

results_protein_features/
  protein_features_light.csv
  protein_features_light_status.tsv
  protein_features_light_summary.tsv
```
