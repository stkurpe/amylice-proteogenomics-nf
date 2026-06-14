# Amyloid Burden Quantile Shift: BC control/AC vs chili/TOL

## What Was Done

This analysis tests whether the upper tail of expression-weighted amyloidogenic burden is shifted upward in baseline ChILI/TOL blood RNA-seq samples compared with control/AC samples.

Only `time == "BC"` samples were used. Groups were mapped from `GSE287540_SraRunTable.csv` as `control -> AC` and `chili -> TOL`.

The transcript universe was restricted to GENCODE protein-coding transcripts (`biotype == "protein_coding"`) present in the expression layer. Protein variants were collapsed to `sample x transcript` before joining to transcript-level TPM, avoiding H1/H2 double-counting. Strict amyloid label was `Consensus_collapsed == "Amyloid"`. Continuous burden was `TPM x Amyloid_Index_max`.

## Samples

- All BC samples: AC/control n = 6, TOL/chili n = 5
- Coverage-QC subset (`protein_coding_annotation_fraction >= 0.5`): AC/control n = 6, TOL/chili n = 4
- QC-excluded samples: SRR32060276

## Main Result

Primary summary metric: `strict_burden_q95`, the 95th percentile of `TPM x Amyloid_Index_max` among strict amyloid protein-coding transcripts in each sample.

All BC samples:

- AC/control mean = 13.6532
- TOL/chili mean = 16.9077
- Delta mean TOL - AC = 3.2544
- Exact empirical one-sided p (TOL > AC) = 0.2900
- Exact empirical two-sided p = 0.5887
- Bootstrap 95% CI for delta = [-6.2853, 12.4069]
- Cliff's delta TOL vs AC = 0.200

Coverage-QC subset:

- AC/control mean = 13.6532
- TOL/chili mean = 17.2044
- Delta mean TOL - AC = 3.5512
- Exact empirical one-sided p (TOL > AC) = 0.3048
- Exact empirical two-sided p = 0.6190
- Bootstrap 95% CI for delta = [-6.1924, 12.6973]
- Cliff's delta TOL vs AC = 0.250

## Tail Metrics Direction

Delta mean TOL - AC:

| metric | all_BC_samples | coverage_qc_protein_annotation_fraction_ge_0.5 |
| --- | --- | --- |
| continuous_burden_q90 | -0.9902 | -1.2837 |
| continuous_burden_q95 | -2.4963 | -3.2662 |
| continuous_top_1pct_burden_mass | -25856.1896 | -161.8542 |
| continuous_top_1pct_positive_burden_mass | -25856.1896 | -161.8542 |
| strict_burden_q90 | 0.9790 | 1.2543 |
| strict_burden_q95 | 3.2544 | 3.5512 |
| strict_top_1pct_burden_mass | 3257.5818 | 5863.8904 |
| strict_top_1pct_positive_burden_mass | 3257.5818 | 5863.8904 |

## Interpretation

The primary strict 95th-percentile burden metric is higher in ChILI/TOL than in control/AC in both the all-sample and coverage-QC analyses. The same direction is seen for the strict 90th percentile and strict top-1% burden mass, including leave-one-out stability for `strict_burden_q95` and `strict_top_1pct_burden_mass`.

This supports the hypothesis at the level of discovery-direction/stability for strict amyloid-labelled protein-coding transcripts. However, empirical p-values are not small and bootstrap CIs cross zero, so this should not be described as a validated biomarker.

Continuous all-annotated burden metrics do not show the same direction and should be treated as exploratory. The most defensible interpretation is therefore narrow: a weak but directionally stable strict-label upper-tail signal, not broad evidence that all amyloid-index-weighted annotated transcripts are shifted upward.

## Limitations

- Small sample size: all BC n = 11, QC n = 10.
- `SRR32060276` has low protein/amyloid annotation coverage and is excluded from the QC subset.
- Amyloid labels and `Amyloid_Index_max` are computational predictions.
- The burden score can be negative because the integrative amyloid index is z-score based.
- This is a sample-level tail analysis and does not adjust for batch, treatment, age, or sex.

## Next Best Step

Validate the tail endpoint in a larger independent baseline cohort, ideally with improved protein/amyloid annotation coverage and an orthogonal amyloid predictor panel for the top contributors.
