# Pipeline Overview

This document describes the end-to-end logic of the amyloid mutant proteome pipeline: what each stage does, why it exists, what files it produces, and how the outputs are combined.

## High-Level Flow

```text
FASTQ / SRA
  -> read QC and preprocessing
  -> expression quantification
  -> alignment and variant calling
  -> expression-filtered mutant proteome
  -> amyloidogenicity prediction
  -> combined reports and S3 archive
```

The current system supports two execution modes.

## 1. Input Data

### Full Mode

The full workflow starts from a public sequencing accession:

```text
SRA accession
  -> FASTQ
```

This mode runs download, QC, expression quantification, alignment, variant calling, proteome generation, and amyloidogenicity prediction.

### Prepared-Input Mode

The currently running batch uses prepared upstream outputs:

```text
results_expression/abundance.tsv
results_gatk/variants_filtered.vcf or variants_filtered.vcf.gz
reference genome FASTA
GTF annotation
  -> mutant proteome
```

In this mode, the workflow does not rerun FASTQ download, FastQC, STAR, Kallisto, or GATK. It starts from existing expression and variant files.

## 2. Read QC

Purpose: evaluate sequencing read quality before downstream analysis.

Tools:

```text
FastQC
fastp
```

Typical outputs:

```text
results_qc/
  FastQC HTML reports
  FastQC ZIP archives
  fastp reports
```

Why this matters:

- Detects low-quality reads.
- Detects adapter contamination.
- Helps identify GC/content anomalies.
- Provides evidence that downstream expression and variant results are trustworthy.

## 3. Expression Quantification

Expression is a central part of this pipeline because the final proteome is sample-specific.

Tool:

```text
Kallisto
```

Inputs:

```text
FASTQ
reference transcriptome
```

Outputs:

```text
results_expression/abundance.tsv
results_expression/abundance.h5
```

The key file is:

```text
abundance.tsv
```

It contains fields such as:

```text
target_id
length
eff_length
est_counts
tpm
```

The most important value for proteome generation is:

```text
TPM
```

TPM estimates how strongly a transcript is expressed in the sample.

### Why Expression Filtering Is Used

The pipeline does not blindly translate every transcript in the annotation. Instead, it keeps only transcripts that are likely to be active in the sample.

The default filter is:

```text
TPM > 1
```

This produces:

```text
results_proteins/clean_ids.txt
```

`clean_ids.txt` is the whitelist of transcripts that:

- are expressed;
- are not mitochondrial;
- have a coding sequence long enough to translate safely;
- can be used for mutant protein generation.

In short:

```text
Expression tells the pipeline which transcripts are biologically relevant for this sample.
```

## 4. Alignment and Variant Calling

Tools:

```text
STAR
GATK
samtools
bcftools
```

Inputs:

```text
FASTQ
reference genome
```

Outputs:

```text
results_star/
  aligned BAM files

results_gatk/
  variants_filtered.vcf
  variants_filtered.vcf.gz
```

Why this matters:

- STAR aligns RNA-seq reads to the genome.
- GATK identifies variants supported by the RNA-seq data.
- The filtered VCF describes SNPs and indels that may alter protein sequences.

The key file for proteome generation is:

```text
results_gatk/variants_filtered.vcf
```

or:

```text
results_gatk/variants_filtered.vcf.gz
```

## 5. Transcript Cleaning

Inputs:

```text
GTF annotation
Kallisto abundance.tsv
```

Filtering logic:

```text
keep transcripts with TPM > 1
remove chrM transcripts
remove transcripts with CDS < 10 bp
```

Why:

- Mitochondrial genes use a different genetic code, so standard translation can create noisy protein sequences.
- Ultra-short CDS records can trigger `gffread` failures.
- Non-expressed transcripts should not contribute to a sample-specific proteome.

Output:

```text
results_proteins/clean_ids.txt
```

This file controls which transcripts are allowed into both SNP/nonsense and frameshift proteome generation.

## 6. SNP and Nonsense Proteome

Inputs:

```text
reference genome
variants_filtered.vcf
clean_ids.txt
GTF annotation
```

Processing:

```text
VCF
  -> SNP-only VCF
  -> haplotype 1 consensus genome
  -> haplotype 2 consensus genome
  -> gffread translation using GTF and clean_ids.txt
  -> raw protein FASTA
  -> cleaned protein FASTA
```

Tools:

```text
bcftools
samtools faidx
gffread
proteome_pipeline
```

Intermediate outputs:

```text
prot_h1.fa
prot_h2.fa
protein_raw.fasta
```

Final cleaned output:

```text
results_proteins/protein.fasta
```

Cleaning behavior:

- Removes problematic protein sequences with excessive internal stop codons.
- Truncates nonsense-mutant proteins at the first stop codon.
- Keeps biologically plausible shortened proteins.

This file is the main SNP/nonsense-derived protein FASTA for the sample.

## 7. Frameshift Proteome

Frameshifts are handled separately because indels cannot be reliably translated by simply applying `bcftools consensus` to the whole genome and then running `gffread`.

Inputs:

```text
variants_filtered.vcf
reference genome
GTF annotation
clean_ids.txt
```

Logic:

```text
find coding indels
map each indel to expressed transcripts
extract the relevant coding DNA
apply REF/ALT surgically
translate the altered coding sequence
stop translation at the first stop codon
deduplicate identical proteins per gene
```

Output:

```text
results_proteins/frameshift_unique.fasta
```

Why this matters:

Frameshift mutations can create entirely new C-terminal protein tails. These novel tails may be highly relevant for amyloidogenicity, aggregation, immune recognition, or other downstream analyses.

The pipeline keeps valid frameshift proteins even if the mutant protein becomes longer than the reference protein.

Default minimum length:

```text
min_frameshift_aa = 6
```

## 8. Mutant Proteome Assembly

Inputs:

```text
protein.fasta
frameshift_unique.fasta
```

Combination:

```text
protein.fasta
  + frameshift_unique.fasta
  -> combine_proteome.fasta
```

Final output:

```text
results_proteins/combine_proteome.fasta
```

This is the main sample-specific mutant proteome.

It contains:

- SNP-derived protein sequences;
- nonsense-truncated proteins;
- frameshift-derived proteins;
- deduplicated valid protein records.

This file is the direct input for amyloidogenicity prediction.

## 9. Proteome Verification

After proteome assembly, the workflow validates the result.

Checks include:

```text
clean_ids.txt exists and is non-empty
protein.fasta exists and is non-empty
frameshift_unique.fasta exists and is non-empty
combine_proteome.fasta exists and is non-empty
combine_proteome.fasta contains protein.fasta records
combine_proteome.fasta contains frameshift_unique.fasta records
```

Output:

```text
results_proteins/verification_report.tsv
```

If everything is valid, rows are marked:

```text
OK
```

If something requires attention, rows are marked:

```text
ATTENTION
```

Only verified outputs should be treated as ready for downstream amyloidogenicity prediction.

## 10. Amyloidogenicity Prediction

Input:

```text
results_proteins/combine_proteome.fasta
```

Two predictors are used.

### AMYPred-FRL

Input:

```text
combine_proteome.fasta
```

Output:

```text
results_amyloid/amypred_frl_prediction.csv
```

AMYPred-FRL uses feature representation learning with multiple base classifiers and a final model to estimate amyloid probability.

### AmyloGram

Input:

```text
combine_proteome.fasta
```

Output:

```text
results_amyloid/amylogram_prediction_fast.csv
```

The optimized runner:

- cleans protein sequences;
- removes very short invalid inputs;
- deduplicates identical sequences;
- runs AmyloGram on unique proteins only;
- processes chunks in parallel;
- expands predictions back to all protein IDs.

## 11. Protein Physicochemical Features

Input:

```text
results_proteins/combine_proteome.fasta
```

Output:

```text
results_protein_features/protein_features.csv
```

The `ProteinFeatures` R step uses `Peptides`, `protr`, `Biostrings`, and `zoo`
to calculate:

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

Alpha and beta propensity come from `protr::AAindex` records:

```text
CHOP780201 — normalized alpha-helix frequency
CHOP780202 — normalized beta-sheet frequency
```

The chameleon score is derived from centered sliding-window alpha/beta
propensity similarity. Sequences shorter than the configured window can have
blank sliding-window aggregate fields.

## 12. Combined Amyloid Table

Inputs:

```text
amypred_frl_prediction.csv
amylogram_prediction_fast.csv
```

Output:

```text
results_amyloid/amyloid_combined_predictions.csv
```

Columns:

```text
Sequence_ID
AMYPred_Prob
AMYPred_Pred
AmyloGram_Prob
AmyloGram_Pred
Consensus
```

Possible `Consensus` values:

```text
Amyloid
Non-Amyloid
Discordant
Partial
```

Meaning:

- `Amyloid`: both predictors agree on amyloidogenicity.
- `Non-Amyloid`: both predictors agree on non-amyloidogenicity.
- `Discordant`: predictors disagree.
- `Partial`: only one predictor produced a result for that sequence.

## 13. S3 Archive

All final files are uploaded to S3.

Proteome outputs:

```text
s3://bucket/sample/results_proteins/
  protein.fasta
  frameshift_unique.fasta
  combine_proteome.fasta
  clean_ids.txt
  nonsense_candidates.txt
  verification_report.tsv
  manifest.txt
  local_proteome.log
  local_proteome.events.tsv
  UPLOAD_SUCCESS
```

Amyloid outputs:

```text
s3://bucket/sample/results_amyloid/
  amypred_frl_prediction.csv
  amylogram_prediction_fast.csv
  amyloid_combined_predictions.csv
  amyloid_predictors_status.tsv
  amyloid_predictors_summary.tsv
```

Protein feature outputs:

```text
s3://bucket/sample/results_protein_features/
  protein_features.csv
  protein_features_status.tsv
  protein_features_summary.tsv
  logs/
```

`UPLOAD_SUCCESS` means the proteome stage uploaded its required outputs.

## 14. How the Pieces Fit Together

```text
Expression tells us what is active.
Variants tell us what is changed.
Proteome generation turns active genetic changes into protein sequences.
Amyloid predictors score those protein sequences for aggregation risk.
Protein feature calculation adds interpretable physicochemical descriptors.
```

In practical terms:

```text
abundance.tsv
  -> clean_ids.txt

variants_filtered.vcf
  -> protein.fasta
  -> frameshift_unique.fasta

protein.fasta + frameshift_unique.fasta
  -> combine_proteome.fasta

combine_proteome.fasta
  -> AMYPred-FRL
  -> AmyloGram
  -> amyloid_combined_predictions.csv
```

## 14. Current Prepared-Input Batch

For the prepared-input batch, the workflow starts from:

```text
s3://bioinfo-data-amylice-2026/<SAMPLE_ID>/results_expression/abundance.tsv
s3://bioinfo-data-amylice-2026/<SAMPLE_ID>/results_gatk/variants_filtered.vcf
```

and writes new proteome outputs to:

```text
s3://codex-test-ngsdata-calculations/prepared-bioinfo-proteome/<SAMPLE_ID>/results_proteins/
```

This mode intentionally avoids rerunning upstream FASTQ, QC, STAR, Kallisto, and GATK stages.
