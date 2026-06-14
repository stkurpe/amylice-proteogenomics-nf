# Rank-Based Amyloid Burden Signature: BC control/AC vs chili/TOL

## What Was Done

This analysis tests whether strict consensus amyloid transcripts sit higher in the within-sample expression rank distribution at baseline (`time == "BC"`) in ChILI/TOL versus control/AC samples. It does not compare absolute TPM as the primary endpoint.

For each sample, GENCODE v48 `transcript_type == "protein_coding"` transcripts present in the expression layer were ranked by TPM. The primary metric is `rank_based_strict_amyloid_burden_score`: the mean expression rank percentile of transcripts with collapsed strict amyloid label `Consensus_collapsed == "Amyloid"`, with score set to 0 when a sample has no strict amyloid transcripts after collapse. Rank percentile is oriented so 1.0 means highest expression and 0.0 means lowest expression within that sample.

Protein variants were collapsed to sample x transcript before scoring, avoiding H1/H2 double-counting against transcript-level expression.

GENCODE source/cache: `/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/rank_based_amyloid_burden_signature/gencode_v48_transcript_annotation.tsv`.

## Samples

- All BC samples: AC/control n = 6, TOL/chili n = 5
- Coverage-QC subset (`protein_coding_annotation_fraction >= 0.5`): AC/control n = 6, TOL/chili n = 4
- QC-excluded samples: SRR32060276

## Main Result

Primary metric: `rank_based_strict_amyloid_burden_score`.

All BC samples:

- AC/control mean = 0.6183
- TOL/chili mean = 0.9141
- Delta mean TOL - AC = 0.2958
- Exact empirical one-sided p (TOL > AC) = 0.2511
- Exact empirical two-sided p = 0.4329
- Bootstrap 95% CI for delta = [-0.0167, 0.6121]
- Cliff's delta TOL vs AC = -0.067

Coverage-QC subset:

- AC/control mean = 0.6183
- TOL/chili mean = 0.9147
- Delta mean TOL - AC = 0.2964
- Exact empirical one-sided p (TOL > AC) = 0.2952
- Exact empirical two-sided p = 0.4286
- Bootstrap 95% CI for delta = [-0.0167, 0.6133]
- Cliff's delta TOL vs AC = 0.000

## Interpretation

The all-sample rank-based burden score is higher in ChILI/TOL because two AC/control samples have no strict amyloid transcripts after sample x transcript collapse. The coverage-QC subset shows the same direction. This supports the hypothesis at the level of effect direction for this discovery endpoint, but the empirical p-values and bootstrap intervals do not justify claiming a validated biomarker.

Sensitivity analysis using only samples with at least one strict amyloid transcript gives the opposite direction:

- All strict-label-positive samples: delta mean TOL - AC = -0.0133, empirical one-sided p (TOL > AC) = 0.9206
- Coverage-QC strict-label-positive samples: delta mean TOL - AC = -0.0127, empirical one-sided p (TOL > AC) = 0.8857

Therefore, the main discovery signal depends on treating absence of strict amyloid transcripts as zero burden. Among samples where strict amyloid transcripts exist, their average expression rank is not higher in ChILI/TOL.

This should be treated as a discovery/pilot signal. Validation requires a larger independent cohort and, ideally, orthogonal amyloidogenicity predictors for top contributors.

## Limitations

- Small sample size limits statistical power.
- The ranked universe now uses independently downloaded GENCODE v48 protein-coding transcripts, but amyloid labels remain available only for transcripts with local protein/amyloid predictions.
- SRR32060276 has low protein/amyloid annotation coverage and is therefore excluded from the QC subset.
- Amyloid labels are computational predictions and need external validation.

## Next Best Step

Use this rank-based burden endpoint as the primary sample-level discovery metric in a larger baseline cohort, with SRR32060276-like samples excluded or repaired by improving annotation coverage. Follow up the top strict amyloid rank contributors with WALTZ/TANGO/AGGRESCAN or an equivalent orthogonal amyloid predictor panel.
