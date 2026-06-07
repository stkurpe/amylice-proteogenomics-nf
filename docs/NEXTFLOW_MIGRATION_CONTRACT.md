# Nextflow Migration Contract

Phase 0 fixes the initial migration contract before replacing bash
orchestration.

## Selected Samples

- `fixture_minimal`: small local prepared-input fixture for fast contract tests.
- `SRR32060234`: validated real sample selected as the legacy bash baseline.

## Prepared Proteome Inputs

Each prepared sample requires:

- `sample_id`
- `results_expression/abundance.tsv`
- `results_gatk/variants_filtered.vcf` or `variants_filtered.vcf.gz`
- reference genome FASTA
- reference GTF annotation

`abundance.tsv` must include:

- `target_id`
- `length`
- `eff_length`
- `est_counts`
- `tpm`

## Prepared Proteome Outputs

The prepared proteome workflow must publish:

- `results_proteins/clean_ids.txt`
- `results_proteins/protein.fasta`
- `results_proteins/frameshift_unique.fasta`
- `results_proteins/combine_proteome.fasta`
- `results_proteins/nonsense_candidates.txt`
- `results_proteins/verification_report.tsv`
- `results_proteins/manifest.txt`

## Phase 0 Checks

- sample manifest is valid CSV and includes both selected samples;
- local fixture input files exist;
- local fixture `abundance.tsv` has required columns;
- local fixture VCF has a valid VCF header;
- local fixture GTF has transcript and CDS rows;
- `proteome_pipeline.cli` imports successfully.
