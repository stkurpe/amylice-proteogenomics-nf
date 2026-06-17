# Pipeline Overview

This document describes the current Nextflow full and prepared-input workflows and the retained production scientific components that are being integrated into them.

Visual biological block scheme: [PIPELINE_BLOCK_SCHEME.md](PIPELINE_BLOCK_SCHEME.md).

## High-Level Flow

```text
SRA accession
  -> FASTQ download
  -> FastQC
  -> Kallisto expression
  -> STAR alignment
  -> GATK variant calling
  -> expression-filtered mutant proteome
  -> proteome verification
  -> AmyloGram-Py prediction
  -> consensus-compatible amyloid table
  -> lightweight protein features
```

## Inputs

Full mode uses a text file with one SRA accession per line and requires prepared reference assets:

```text
reference directory with genome FASTA sidecars
reference GTF annotation
STAR genome index
Kallisto transcriptome index
AmyloGram-Py six-mer table
```

Prepared mode uses local paths from the sample manifest:

```text
sample_id
abundance.tsv
variants_filtered.vcf
```

It also requires:

```text
reference genome FASTA
reference GTF annotation
AmyloGram-Py six-mer table
```

Prepared mode intentionally starts from prepared expression and variant files.

## Upstream Processing

Full mode adds native Nextflow processes for:

```text
DOWNLOAD_FASTQ
FASTQC_READS
KALLISTO_QUANT
STAR_ALIGN
MARK_DUPLICATES
INDEX_DEDUP_BAM
GATK_SPLIT_N_CIGAR
INDEX_SPLIT_BAM
GATK_HAPLOTYPE_CALLER
GATK_VARIANT_FILTRATION
```

The final upstream outputs passed into proteome generation are:

```text
results_expression/abundance.tsv
variants_filtered.vcf.gz
```

## Proteome Generation

The proteome branch keeps expressed, non-mitochondrial transcripts with valid CDS records, then reconstructs sample-specific protein sequences from SNP/nonsense and frameshift variants.

By default, Nextflow uses `--min_cds_bp 10` to remove ultra-short CDS records
that produce 1-3 amino-acid peptides. For exact reproduction of the AWS
reference proteome, use `-profile docker,aws_reference`, which sets
`--min_cds_bp 1` and retains those ultra-short CDS-derived proteins.

Main outputs:

```text
results_proteins/
  clean_ids.txt
  protein.fasta
  frameshift_unique.fasta
  combine_proteome.fasta
  nonsense_candidates.txt
  verification_report.tsv
  manifest.txt
```

`VERIFY_PROTEOME` checks that required outputs exist and that the combined proteome is suitable for downstream prediction.

## AmyloGram-Py

Input:

```text
results_proteins/combine_proteome.fasta
```

Outputs:

```text
results_amyloid/
  amylogram_py_prediction.csv
  amylogram_py_report.json
  amylogram_py_report.md
  amylogram_py_skipped.tsv
  amylogram_py_top_hits.tsv
```

The predictor output is normalized into:

```text
results_amyloid/amyloid_combined_predictions.csv
```

Current columns:

```text
Sequence_ID
AMYPred_Prob
AMYPred_Pred
AmyloGramPy_Prob
AmyloGramPy_Pred
Consensus
```

`AMYPred-FRL` is retained under `amyloid_predictors/AMYPred-FRL` and is wired into the prepared Nextflow branch behind `--run_amypred true`. The process requires the complete AMYPred runtime assets:

```text
amyloid_predictors/AMYPred-FRL/data/PAAC.txt
amyloid_predictors/AMYPred-FRL/data/TR_P_132.fasta
amyloid_predictors/AMYPred-FRL/data/TR_N_305.fasta
amyloid_predictors/AMYPred-FRL/model/pima.pickle_model_svm_PF.dat
```

When `--run_amypred false`, AMYPred fields are emitted as placeholders and `Consensus` is `Partial` when AmyloGram-Py produced a usable class.

Current `Consensus` values:

```text
Partial
NA
```

Summary and status files:

```text
results_amyloid/amyloid_predictors_summary.tsv
results_amyloid/amyloid_predictors_status.tsv
```

## Lightweight Protein Features

Input:

```text
results_proteins/combine_proteome.fasta
```

Outputs:

```text
results_protein_features/
  protein_features_light.csv
  protein_features_light_summary.tsv
  protein_features_light_status.tsv
```

The light feature module is pure Python and computes simple sequence-level descriptors:

```text
protein_length
KD_mean
KD_max
KD_sd
charge
charge_density
aromaticity
alpha_propensity
beta_propensity
chameleon_score
pI
```

`pI` is currently present as a schema field and left blank.

The full `ProteinFeatures` implementation is retained under `amyloid_predictors/ProteinFeatures` and is pending a dedicated Nextflow process.

## Workflow Modes

Prepared mode:

```bash
nextflow run nextflow/main.nf -profile test,docker --mode prepared
```

AWS/reference-compatible prepared mode:

```bash
nextflow run nextflow/main.nf -profile docker,aws_reference --mode prepared
```

AWS/reference-compatible prepared mode with AMYPred-FRL:

```bash
nextflow run nextflow/main.nf -profile docker,aws_reference --mode prepared --run_amypred true
```

Amyloid-only mode:

```bash
nextflow run nextflow/main.nf -profile test,docker --mode amyloid --input_fasta PATH
```
