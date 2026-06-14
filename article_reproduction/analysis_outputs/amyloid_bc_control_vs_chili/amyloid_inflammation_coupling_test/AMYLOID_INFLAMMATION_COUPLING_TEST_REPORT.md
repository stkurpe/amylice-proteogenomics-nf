# Amyloid-Inflammation Coupling Test

## What Was Done

Method 9 tested whether baseline amyloid burden is coupled to inflammatory RNA modules, rather than only asking whether amyloid burden is higher in ChILI/TOL.

Data were recalculated from `project.zarr` for `time == "BC"` only, with `control -> AC` and `chili -> TOL`. Protein variants were collapsed to `sample x transcript` before expression weighting to avoid H1/H2 double-counting. The primary amyloid score was `strict_amyloid_TPM_fraction`, where strict amyloid means `Consensus_collapsed == "Amyloid"`. The continuous score used `Amyloid_Index_max`.

Immune module scores were mean z-scored `log2(gene TPM + 1)` across predefined blood/inflammation panels:

- myeloid/neutrophil
- interferon
- antigen presentation
- proteasome/UPR

The main model was `immune_score ~ amyloid_score_z * Group`; the interaction term asks whether the amyloid/immune slope differs in TOL/chili versus AC/control. Two analysis sets were run: all BC samples and coverage-QC subset excluding `protein_coding_annotation_fraction < 0.5`. Empirical interaction p-values were computed by exhaustive group-label permutation with the same group sizes.

## Samples Used

- `SRR32060276` / `N132`: chili -> TOL, excluded from QC, strict fraction=0.2036
- `SRR32060306` / `N101`: chili -> TOL, QC, strict fraction=0.08999
- `SRR32060309` / `N87`: control -> AC, QC, strict fraction=0
- `SRR32060312` / `N106`: control -> AC, QC, strict fraction=0.1391
- `SRR32060326` / `N032`: control -> AC, QC, strict fraction=0.1552
- `SRR32060336` / `N072`: chili -> TOL, QC, strict fraction=0.1517
- `SRR32060341` / `N064`: chili -> TOL, QC, strict fraction=0.1832
- `SRR32060215` / `N015`: control -> AC, QC, strict fraction=0.07492
- `SRR32060216` / `N022`: chili -> TOL, QC, strict fraction=0.1049
- `SRR32060234` / `N159`: control -> AC, QC, strict fraction=0.1026
- `SRR32060240` / `N149`: control -> AC, QC, strict fraction=0

## QC Exclusions

- `SRR32060276` (TOL), protein_coding_annotation_fraction=0.08725

## Main Result

Primary burden endpoint in the coverage-QC subset:

- AC/control mean `strict_amyloid_TPM_fraction`: 0.07864
- TOL/chili mean `strict_amyloid_TPM_fraction`: 0.1325
- Delta TOL - AC: 0.05382
- Bootstrap 95% CI: -0.002931 to 0.1184
- Welch p-value: 0.1601

Primary coupling test for myeloid/neutrophil activation in the coverage-QC subset:

- Pearson r AC/control: 0.1622
- Pearson r TOL/chili: 0.7523
- Pearson r difference TOL - AC: 0.5901
- OLS interaction coefficient TOL - AC: 0.5417
- OLS nominal interaction p-value: 0.4602
- Exact empirical interaction p-value: 0.5024
- Bootstrap interaction 95% CI: -1.156 to 4.203

## Interpretation

The burden result remains directionally consistent with the hypothesis: baseline ChILI/TOL has higher strict amyloid TPM fraction than AC/control after excluding the low-coverage sample. The coupling result is also directionally consistent with the amyloid-inflammation hypothesis: in the coverage-QC subset, amyloid burden is more strongly associated with the myeloid/neutrophil module in TOL/chili than in AC/control.

This remains discovery/exploratory, not validation. The exact label-permutation p-value is not significant, the nominal OLS interaction p-value is weak, and the QC TOL/chili group has only 4 samples. Across modules, use the interaction table to separate prioritization signals from claims: a positive interaction means amyloid burden is more positively associated with the immune module in TOL/chili than in AC/control; a negative interaction means the opposite.

## Limitations

- Very small sample size: all BC is 6 AC vs 5 TOL; coverage-QC is 6 AC vs 4 TOL.
- `SRR32060276` is biologically and statistically influential because it has very low protein/amyloid annotation coverage and was excluded from the primary QC subset.
- Module scores are compact predefined panels, not externally validated latent factors in this cohort.
- Interaction models with 10-11 samples are unstable; effect size and direction are more informative than p-values.
- The analysis is discovery/exploratory, not validation.

## Next Best Step

Validate the coupling signal in an independent or larger irAE blood RNA-seq cohort, using the same frozen amyloid burden endpoint and predefined immune module panels. If no external cohort is available, the next internal step is leave-one-subject-out stability for the interaction coefficients plus a counts-based module approach that adjusts for batch and treatment where degrees of freedom permit.
