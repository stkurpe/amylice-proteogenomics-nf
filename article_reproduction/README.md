# Amylice Article Reproduction Bundle

This directory contains the downstream materials used to reproduce manuscript result figures for:

**Amylice: a reproducible RNA-seq proteogenomic workflow for estimating patient-specific amyloidogenic proteome burden in cancer immunotherapy toxicity**

## Contents

- `Amylice_reproduce_article_results.ipynb` - notebook that displays Figures 3-8 with figure text and runnable code cells.
- `technical_validation_figure3_source.tsv` - processed source table for Figure 3 technical-validation metrics.
- `generate_technical_validation_figure3.py` - script that regenerates Figure 3 from `technical_validation_figure3_source.tsv`.
- `analysis_outputs/` - processed downstream tables, reports, and PNG/PDF figure outputs.
- `project.zarr/` - processed discovery-cohort zarr store used by the downstream figure/statistics scripts.
- `validation_project.zarr/` - processed validation-cohort zarr manifest/store.
- `GSE287540_SraRunTable.csv` - sample metadata used for cohort/group selection.
- `pdf_extracted_text.txt` - extracted manuscript text used to verify figure headings.
- `*.py` - downstream analysis and figure-generation scripts.

## Reproducing The Notebook Figures

Open `Amylice_reproduce_article_results.ipynb` from this directory.

By default, the notebook displays already generated local PNG files. To rerun figure-generation scripts before display, set:

```python
RUN_RECOMPUTE = True
```

Figure 3 is regenerated from `technical_validation_figure3_source.tsv`. Figures 4-8 point to downstream generated result panels in `analysis_outputs/`.

## Scope

This bundle reproduces the manuscript's processed analytical layer. It does not rerun the full FASTQ-to-proteome workflow from raw SRA reads. Full end-to-end processing requires the external tools described in the main repository documentation, including SRA Toolkit, FastQC, Kallisto, STAR, GATK, bcftools, gffread, Docker, and/or Nextflow.
