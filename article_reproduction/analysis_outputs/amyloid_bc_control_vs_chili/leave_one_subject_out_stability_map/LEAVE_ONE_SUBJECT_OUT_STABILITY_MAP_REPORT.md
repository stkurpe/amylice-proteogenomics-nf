# Leave-One-Subject-Out Stability Map

## What Was Done

Recomputed baseline (`time == "BC"`) expression-weighted amyloidogenic burden directly from `project.zarr` and `GSE287540_SraRunTable.csv`.
Groups were mapped as `control -> AC` and `chili -> TOL`. Protein H1/H2 variants were collapsed to `sample x transcript` before joining transcript-level TPM, preventing expression double-counting. Strict amyloid label was `Consensus_collapsed == "Amyloid"`; continuous score was `Amyloid_Index_max`.

Primary endpoint: `strict_amyloid_TPM_fraction`.

## Samples

- All BC samples: 6 AC/control and 5 TOL/chili.
- Coverage-QC subset (`protein_coding_annotation_fraction >= 0.5`): 6 AC/control and 4 TOL/chili.
- QC-excluded samples: SRR32060276.

All BC AC/control samples: SRR32060215, SRR32060234, SRR32060240, SRR32060309, SRR32060312, SRR32060326.

All BC TOL/chili samples: SRR32060216, SRR32060276, SRR32060306, SRR32060336, SRR32060341.

Coverage-QC AC/control samples: SRR32060215, SRR32060234, SRR32060240, SRR32060309, SRR32060312, SRR32060326.

Coverage-QC TOL/chili samples: SRR32060216, SRR32060306, SRR32060336, SRR32060341.

## Main Result

Primary endpoint, all BC samples:

- AC/control mean: 0.0765011
- TOL/chili mean: 0.139244
- Delta mean TOL - AC: 0.062743
- Hedges g: 0.972
- Exact empirical one-sided p (`TOL > AC`): 0.05616
- Bootstrap 95% CI for delta: [0.00111376, 0.124392]
- LOSO direction stability: 11/11 (100.0%)

Primary endpoint, coverage-QC subset:

- AC/control mean: 0.0765011
- TOL/chili mean: 0.123834
- Delta mean TOL - AC: 0.0473326
- Hedges g: 0.741
- Exact empirical one-sided p (`TOL > AC`): 0.109
- Bootstrap 95% CI for delta: [-0.0112636, 0.10655]
- LOSO direction stability: 10/10 (100.0%)

## Secondary Endpoints, Coverage-QC Subset

| Metric | Delta mean TOL - AC | Hedges g | Empirical p, one-sided TOL > AC | LOSO direction |
|---|---:|---:|---:|---:|
| `strict_amyloid_TPM_fraction` | 0.0473326 | 0.741 | 0.109 | 10/10 (100.0%) |
| `strict_amyloid_TPM` | 32345.6 | 0.677 | 0.1469 | 10/10 (100.0%) |
| `n_strict_amyloid_expressed_TPM_ge_1` | 493.083 | 0.931 | 0.07109 | 10/10 (100.0%) |
| `strict_weighted_index_sum` | -8209.69 | -0.739 | 0.8863 | 0/10 (0.0%) |
| `continuous_expression_weighted_amyloid_index` | 0.0856186 | 0.653 | 0.1611 | 10/10 (100.0%) |

The top strict amyloid contributor panel was regenerated as `top_strict_amyloid_contributor_panel.tsv`.

## Interpretation

For the primary endpoint, the direction is higher in ChILI/TOL than control/AC in both all-sample and coverage-QC analyses. However, exact permutation p-values and bootstrap confidence intervals remain weak/wide, consistent with a pilot discovery signal rather than a validated biomarker. LOSO stability shows whether the direction survives removal of each subject and is the preferred robustness readout for this small cohort.

## Limitations

- Small sample size: 6 vs 5 all BC; 6 vs 4 after coverage QC.
- `SRR32060276` has very low protein/amyloid annotation coverage and inflates uncertainty in all-sample analyses.
- Amyloid labels are computational predictions and need orthogonal validation.
- This is discovery/exploratory stability analysis; no independent validation cohort is included here.

## Next Best Step

Validate the primary endpoint and the top strict amyloid contributor panel in an independent baseline cohort or by a prespecified cross-validation/permutation framework, then re-score top contributors with slower specialized aggregation predictors.
