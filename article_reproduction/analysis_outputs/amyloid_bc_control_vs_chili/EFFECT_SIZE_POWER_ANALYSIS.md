# Effect Size and Sample Size Planning

This analysis estimates effect sizes from the current pilot baseline BC comparison:

- `AC/control`: 6 samples
- `TOL/chili`: 5 samples
- coverage-QC subset: 6 AC vs 4 TOL after excluding low-coverage `SRR32060276`

The primary hypothesis-relevant metrics are:

1. `strict_amyloid_TPM_fraction`
2. `strict_amyloid_TPM`
3. `n_strict_amyloid_expressed_TPM_ge_1`

## Pilot Effect Sizes

Coverage-QC estimates:

| Metric | Direction | Hedges g | Cliff's delta | Welch p | n/group for 80% power | n/group for 90% power |
|---|---|---:|---:|---:|---:|---:|
| strict_amyloid_TPM_fraction | higher in TOL | 0.82 | -0.50 | 0.160 | 25 | 33 |
| strict_amyloid_TPM | higher in TOL | 0.68 | -0.33 | 0.220 | 36 | 47 |
| n strict amyloid expressed transcripts | higher in TOL | 0.94 | -0.58 | 0.099 | 19 | 25 |
| strict weighted index sum | more negative in TOL | 0.79 | 0.50 | 0.242 | 27 | 35 |

All-sample estimates:

| Metric | Direction | Hedges g | Cliff's delta | Welch p | n/group for 80% power | n/group for 90% power |
|---|---|---:|---:|---:|---:|---:|
| strict_amyloid_TPM_fraction | higher in TOL | 1.04 | -0.60 | 0.084 | 16 | 21 |
| strict_amyloid_TPM | higher in TOL | 0.32 | -0.20 | 0.563 | 150 | 201 |
| n strict amyloid expressed transcripts | higher in TOL | 0.45 | -0.40 | 0.433 | 80 | 107 |

## Recommended Planning Number

The current pilot estimates are unstable because the cohort is small and one TOL sample has low annotation coverage. For a defensible study design, do not plan around the most optimistic effect size (`d ~ 1.0`). Use a conservative-to-moderate large target effect:

- If expecting `d = 0.8`: about **26 samples per group** for 80% power; **34 per group** for 90% power.
- If expecting `d = 0.7`: about **34 samples per group** for 80% power; **44 per group** for 90% power.
- If expecting `d = 0.6`: about **45 samples per group** for 80% power; **60 per group** for 90% power.

Practical recommendation:

> Aim for at least **35-45 baseline samples per group** to confirm or refute the baseline amyloidogenic burden hypothesis robustly.

Minimal follow-up:

> At least **25-30 samples per group** could test the strongest pilot endpoint (`strict_amyloid_TPM_fraction`) but would remain vulnerable to overestimated pilot effects.

## Generic Two-Sample Planning Scenarios

| Assumed absolute Cohen's d | n/group for 80% power | n/group for 90% power |
|---:|---:|---:|
| 0.4 | 100 | 133 |
| 0.5 | 64 | 86 |
| 0.6 | 45 | 60 |
| 0.7 | 34 | 44 |
| 0.8 | 26 | 34 |
| 1.0 | 17 | 23 |
| 1.2 | 12 | 16 |
| 1.5 | 9 | 11 |

## Interpretation

The pilot data support a directionally higher strict amyloidogenic expression burden in baseline ChILI/TOL blood, but the current sample size is not sufficient for confirmation.

The best-powered and most biologically interpretable endpoint is:

`strict_amyloid_TPM_fraction`

Secondary endpoints:

- `n_strict_amyloid_expressed_TPM_ge_1`
- `strict_amyloid_TPM`
- top contributor burden among inflammatory/proteostasis genes

For the manuscript, frame the current data as a pilot effect-size estimation study and use the resulting sample-size calculation to justify expansion or validation.

## Files

- `effect_size_power_analysis.tsv`
- `generic_effect_size_sample_size_scenarios.tsv`
