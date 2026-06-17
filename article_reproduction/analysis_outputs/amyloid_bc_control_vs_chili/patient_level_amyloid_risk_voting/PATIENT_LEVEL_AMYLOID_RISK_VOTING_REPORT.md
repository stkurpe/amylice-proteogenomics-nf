# Patient-Level Amyloid Risk Voting

## What Was Done

Recomputed a sample-level voting endpoint for baseline blood RNA-seq samples only (`time == "BC"`), comparing AC/control vs TOL/chili. Protein variants were collapsed to sample x transcript before strict amyloid labels were used, so H1/H2 protein variants do not double-count transcript-level expression. Strict amyloid label was `Consensus_collapsed == "Amyloid"`; continuous score input was `Amyloid_Index_max` as already summarized into `continuous_expression_weighted_amyloid_index`.

Each sample received one vote for each abnormal high feature, using a control-anchored threshold: high means the value is greater than the maximum AC/control value in the same analysis set.

Votes:

- `strict_amyloid_TPM_fraction` high
- `n_strict_amyloid_transcripts` high, using expressed strict amyloid transcripts with TPM >= 1
- nonsense burden high, using unique nonsense transcript count
- S100A9/B2M/LYZ/CLU strict amyloid contributor high
- `continuous_expression_weighted_amyloid_index` high

## Samples

All BC samples: SRR32060215(AC/control, N015), SRR32060216(TOL/chili, N022), SRR32060234(AC/control, N159), SRR32060240(AC/control, N149), SRR32060276(TOL/chili, N132), SRR32060306(TOL/chili, N101), SRR32060309(AC/control, N87), SRR32060312(AC/control, N106), SRR32060326(AC/control, N032), SRR32060336(TOL/chili, N072), SRR32060341(TOL/chili, N064)

Coverage-QC exclusion threshold: `protein_coding_annotation_fraction < 0.5`.

Excluded from QC subset: SRR32060276(TOL/chili, annotation fraction=0.087).

## Main Result

All BC samples: mean risk vote score AC=0, TOL=0.6; median AC=0, TOL=1; delta mean TOL-AC=0.6; Cliff's delta=0.6; Mann-Whitney one-sided p(TOL>AC)=0.02256; exact permutation p(TOL>AC)=0.06263; bootstrap 95% CI for mean delta=[0.2, 1].

Coverage-QC subset: mean risk vote score AC=0, TOL=0.5; median AC=0, TOL=0.5; delta mean TOL-AC=0.5; Cliff's delta=0.5; Mann-Whitney one-sided p(TOL>AC)=0.04609; exact permutation p(TOL>AC)=0.1374; bootstrap 95% CI for mean delta=[0, 1].

## Interpretation

The direction is consistent with the hypothesis: TOL/chili samples have higher risk vote scores than AC/control samples in both all-BC and coverage-QC analyses. However, the exact permutation evidence is weak after QC, so this should be treated as exploratory support, not validation of a biomarker.

The signal is driven by strict amyloid burden features: strict amyloid TPM fraction and number of expressed strict amyloid transcripts. Nonsense burden, S100A9/B2M/LYZ/CLU contributor burden, and continuous amyloid index did not exceed the control-anchored high threshold in TOL/chili samples.

## Limitations

Small n, control-derived thresholds from the same dataset, treatment imbalance, batch/sex/age confounding risk, incomplete protein/amyloid annotation for SRR32060276, and no independent validation cohort.

## Next Best Step

Run a validation-oriented model on a pre-registered feature set: either exact/penalized patient-level logistic testing of the five votes or counts-based transcript-level DE followed by locked risk-score evaluation in an independent ChILI/control dataset.
