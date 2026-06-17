# Directional Replication Score: GENCODE Protein-Coding Fixed Panel

## What Was Done

Method: `directional_replication_score_gencode_protein_coding_panel`.
I re-read `project.zarr`, restricted the comparison to baseline `time == "BC"`, mapped `control -> AC` and `chili -> TOL`, and used only GENCODE `biotype == "protein_coding"` transcripts.
Protein variants were collapsed to `sample x transcript` before joining expression, so H1/H2 protein variants cannot double-count transcript-level expression.
Strict amyloid labels were `Consensus_collapsed == "Amyloid"`; continuous burden used `TPM x Amyloid_Index_max`.

Fixed validation panel: ZNF580, STMP1, FXYD5, SEPTIN7, AP2B1, YLPM1, SNX10, TACC1, HBP1, PSMD8.
GENCODE source: `/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/rank_based_amyloid_burden_signature/gencode_v48_transcript_annotation.tsv`.

## Samples

All BC samples: 11 total; AC/control=6, TOL/chili=5.
Coverage-QC subset: 10 total; AC/control=6, TOL/chili=4.

QC-excluded samples:
| sample | Group | raw_Group | protein_coding_annotation_fraction |
| --- | --- | --- | --- |
| SRR32060276 | TOL | chili | 0.0976868 |

## Main Result

All BC samples: 4/10 fixed-panel genes had higher median strict amyloid TPM in TOL than AC; fixed-panel directional replication score=0.400, exact binomial p with ties treated as not replicated=0.8281. Among non-tie genes only: 4/4, sign-test p=0.0625.
Coverage-QC subset: 6/10 fixed-panel genes had higher median strict amyloid TPM in TOL than AC; fixed-panel directional replication score=0.600, exact binomial p with ties treated as not replicated=0.377. Among non-tie genes only: 6/6, sign-test p=0.01562.

Primary panel burden endpoint (`strict_amyloid_TPM_fraction`):
All BC: mean TOL−AC=-0.0276058, median TOL−AC=0.0762542, Cliff's delta=0.067, Mann-Whitney greater p=0.4633.
Coverage-QC: mean TOL−AC=0.0337605, median TOL−AC=0.0897319, Cliff's delta=0.250, Mann-Whitney greater p=0.2965.

## Gene-Level Direction

| analysis_set | gene_name | median_strict_amyloid_TPM_AC | median_strict_amyloid_TPM_TOL | delta_median_TOL_minus_AC | replicated_expected_direction |
| --- | --- | --- | --- | --- | --- |
| all_BC_samples | ZNF580 | 0 | 0 | 0 | False |
| all_BC_samples | STMP1 | 0 | 0 | 0 | False |
| all_BC_samples | FXYD5 | 25.3097 | 37.9834 | 12.6737 | True |
| all_BC_samples | SEPTIN7 | 0 | 0 | 0 | False |
| all_BC_samples | AP2B1 | 0 | 0 | 0 | False |
| all_BC_samples | YLPM1 | 0 | 0 | 0 | False |
| all_BC_samples | SNX10 | 73.9572 | 138.827 | 64.8698 | True |
| all_BC_samples | TACC1 | 0 | 0 | 0 | False |
| all_BC_samples | HBP1 | 0 | 19.9331 | 19.9331 | True |
| all_BC_samples | PSMD8 | 0 | 1.54958 | 1.54958 | True |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | ZNF580 | 0 | 0 | 0 | False |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | STMP1 | 0 | 0 | 0 | False |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | FXYD5 | 25.3097 | 47.4746 | 22.1649 | True |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | SEPTIN7 | 0 | 4.00423 | 4.00423 | True |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | AP2B1 | 0 | 32.8905 | 32.8905 | True |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | YLPM1 | 0 | 0 | 0 | False |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | SNX10 | 73.9572 | 157.049 | 83.0923 | True |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | TACC1 | 0 | 0 | 0 | False |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | HBP1 | 0 | 20.2758 | 20.2758 | True |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | PSMD8 | 0 | 2.91076 | 2.91076 | True |

## Interpretation

This is a small fixed-panel validation-style analysis, not biomarker proof. The relevant signal is direction/stability across predefined contributors, with p-values treated as descriptive because the sample size is small.

## Limitations

- Panel was fixed from prior/discovery contributors, so this should be read as directional validation/exploration rather than independent biomarker discovery.
- Amyloid calls and `Amyloid_Index_max` are computational predictions.
- The cohort is small and includes treatment/batch heterogeneity.
- One TOL sample, `SRR32060276`, has very low protein/amyloid annotation coverage and is excluded in the QC subset.

## Next Best Step

Run the same fixed panel on an independent or held-out baseline cohort with the panel frozen, then report only direction score and effect sizes as the primary validation readout.
