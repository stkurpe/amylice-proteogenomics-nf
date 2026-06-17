# Amyloidogenicity at Baseline BC: control/AC vs chili/TOL

## Project Context

Project hypothesis: amyloidogenic potential of the expressed/personalized proteome may act as a biomarker or mechanistic contributor to immune-related adverse events during immune checkpoint inhibitor treatment. In GSE287540, `chili` represents checkpoint inhibitor-induced liver injury (ChILI) and `control` represents cancer patients without irAE. The most relevant predictive question is whether baseline (`time == BC`) blood already contains an amyloidogenic expression/proteome signature associated with later ChILI/TOL status.

Working labels used in this analysis:

- `control -> AC`
- `chili -> TOL`
- timepoint: `BC`

## Data Analyzed

Local data sources:

- `project.zarr`
- `GSE287540_SraRunTable.csv`
- notebooks in `Article/`
- project plan and progress report

BC samples included:

- AC/control: 6 samples
- TOL/chili: 5 samples

One sample had very low amyloid/protein annotation coverage:

- `SRR32060276`, TOL/chili, protein/proteome annotation coverage = 0.087

Primary statistics therefore use a coverage QC subset excluding samples with annotation fraction `< 0.5`.

## Main Sample-Level Result

The strict consensus amyloid expression burden is directionally higher in baseline ChILI/TOL samples than in controls.

Coverage-QC subset:

| Metric | AC/control mean | TOL/chili mean | Direction | p-value |
|---|---:|---:|---|---:|
| strict_amyloid_TPM | 74,691 | 114,685 | higher in TOL | Welch p = 0.220 |
| strict_amyloid_TPM_fraction | 0.0786 | 0.1325 | higher in TOL | Welch p = 0.160 |
| n strict amyloid expressed transcripts | 1,215 | 2,028 | higher in TOL | Welch p = 0.099 |

Interpretation: no metric reaches statistical significance after multiple testing, but the effect direction is consistent with the project hypothesis: ChILI/TOL baseline blood shows a higher fraction and count of expressed strict-consensus amyloidogenic transcripts.

## Differential Expression

Counts-based exploratory DE used library-size normalized logCPM + Welch tests. Inputs for external DESeq2 were also saved.

Results:

- transcripts tested: 21,045
- FDR < 0.05: 0
- strict amyloid transcripts with raw p < 0.05: 93
- strict amyloid transcripts with FDR < 0.05: 0

Interpretation: at this pilot sample size, transcript-level DE is underpowered. Raw-p candidates can guide validation but should not be claimed as significant biomarkers.

Top raw-p strict amyloid candidates higher in TOL include:

- `CS`
- `CRACR2A`
- `RFC1`
- `ZNF790`
- `CFLAR`
- `GINS1`

Top raw-p strict amyloid candidates higher in AC include:

- `NABP1`
- `CD164`
- `SRRM2`
- `HSD17B7`

## Top Contributors To Amyloid Burden

These are not necessarily DE-significant; they are transcripts that most drive group-level amyloid burden differences because they combine expression and amyloid score.

Higher strict amyloid contribution in AC/control:

- `PTMA`
- `S100A8`
- `RPL9`
- `AKAP13`
- `MTDH`
- `CD52`
- `ATP5MG`
- `RBM39`
- `ANP32B`
- `PRR13`

Higher strict amyloid contribution in TOL/chili:

- `HBA1`
- `OAZ1`
- `S100A9`
- `B2M`
- `RPL21`
- `LYZ`
- `MYL6`
- `RPL37`
- `SLIRP`
- `AP2B1`
- `CLU`

Biological interpretation:

- TOL/chili contributors include inflammation/innate immune and blood-cell-associated genes such as `S100A9`, `B2M`, `LYZ`, `CLU`, and hemoglobin-related signal (`HBA1`).
- AC/control contributors include `S100A8`, `PTMA`, ribosomal/translation-related transcripts, and nuclear/cytosolic proteins.
- The `S100A8/S100A9` split is especially interesting: both are inflammatory alarmins and neutrophil/myeloid markers; their differential contribution may reflect pre-irAE innate immune state rather than amyloidogenicity alone.

## Subcellular Context

Among top 100 absolute strict-contribution transcripts with HPA annotation:

- cytosol: 35
- nucleus: 46
- mitochondria: 9
- membrane: 8
- extracellular: 8

Top TOL/chili contributors include extracellular or secreted/immune-relevant proteins (`S100A9`, `B2M`, `LYZ`, `CLU`), which are biologically plausible as inflammation/proteostasis mediators.

## Nonsense/Truncation Candidates

Nonsense candidate counts were directionally higher in TOL/chili:

- AC/control mean: 36
- TOL/chili mean: 57
- Mann-Whitney p = 0.143
- Cliff's delta AC vs TOL = -0.567

Interpretation: not significant, but directionally consistent with a broader altered personalized proteome in ChILI/TOL baseline samples.

## HLA

HLA calls were summarized but sample size is too small for allele-level statistical inference. HLA should be treated as mechanistic context and covariate/stratification feature in a larger model rather than a standalone result in this pilot analysis.

## Key Limitations

1. GENCODE could not be downloaded in the local run, so gene symbols were taken from the expression layer and strict protein-coding filtering could not be independently revalidated locally. The zarr protein-feature layer still anchors the analysis to translated/predicted protein sequences.
2. Sample size is small: 6 control and 5 ChILI baseline samples; after coverage QC, 6 vs 4.
3. DE is underpowered; no transcript passes FDR < 0.05.
4. Transcript expression does not distinguish H1/H2 personalized protein variants. The analysis collapses variants to transcript level to avoid double-counting TPM.
5. Amyloid predictors are screening tools. Top candidates need validation with slower/more specific tools such as WALTZ, TANGO, AGGRESCAN, PASTA/FoldAmyloid.

## Recommended Main Strategy

The strongest strategy is not to claim individual DE biomarkers yet. Instead, frame the paper around a reproducible proteogenomic workflow and a pilot signal:

> Baseline ChILI/TOL blood shows a directionally higher expression burden of strict-consensus amyloidogenic transcripts, especially among immune/inflammatory and proteostasis-associated genes, supporting amyloidogenic burden as a candidate pre-irAE risk feature.

Primary analysis hierarchy:

1. **Sample-level burden as the main endpoint**
   - `strict_amyloid_TPM_fraction`
   - `strict_amyloid_TPM`
   - `n_strict_amyloid_expressed_TPM_ge_1`

2. **Top burden contributors as biological interpretation**
   - emphasize `S100A9`, `B2M`, `LYZ`, `CLU`, `HBA1` for TOL/chili
   - compare with AC/control contributors such as `PTMA`, `S100A8`, `RPL9`

3. **Counts-based DE as supportive/exploratory**
   - run external DESeq2 using saved `counts_for_external_DESeq2.tsv` and `metadata_for_external_DESeq2.tsv`
   - use DE to prioritize candidates, not as the central result

4. **Personalized proteome/nonsense/HLA as mechanistic layers**
   - nonsense burden: supportive directional evidence
   - HLA: context for immune presentation, not powered for association

5. **Validation plan**
   - validate top strict amyloid contributors with WALTZ/TANGO/AGGRESCAN
   - test whether burden predicts future ChILI in all available baseline samples
   - extend to post-treatment longitudinal changes: BC -> TOL/AC within subject

## Files Generated

- `sample_level_amyloid_burden.tsv`
- `sample_level_burden_statistics.tsv`
- `counts_logCPM_DE_with_amyloid_annotation.tsv`
- `counts_for_external_DESeq2.tsv`
- `metadata_for_external_DESeq2.tsv`
- `strict_amyloid_top_contributors.tsv`
- `strict_amyloid_top_contributors_with_hpa.tsv`
- `nonsense_count_per_sample.tsv`
- `nonsense_count_statistics.tsv`
- `hla_calls_BC.tsv`
- `hla_allele_group_summary.tsv`
- figures in `figures/`
