# Counterfactual Permutation by Amyloid Labels

## What was done

Baseline (`time == BC`) blood RNA-seq samples were compared as AC/control versus TOL/chili. Expression was kept fixed, while amyloid labels and continuous amyloid index values were shuffled within bins defined by protein length and baseline mean transcript expression. Protein variants were collapsed to `sample x transcript` before calculating burden, so H1/H2-style variant duplication does not double-count expression.

Permutation count: 10,000. One-sided empirical p-value tests whether observed `TOL - AC` burden is greater than expected under random amyloid label placement among similar proteins.

## Samples

| sample      | subject | Group | raw_Group | time | protein_coding_annotation_fraction | keep_coverage_qc |
| ----------- | ------- | ----- | --------- | ---- | ---------------------------------- | ---------------- |
| SRR32060276 | N132    | TOL   | chili     | BC   | 0.0872549                          | False            |
| SRR32060306 | N101    | TOL   | chili     | BC   | 0.937622                           | True             |
| SRR32060309 | N87     | AC    | control   | BC   | 0.769414                           | True             |
| SRR32060312 | N106    | AC    | control   | BC   | 0.917182                           | True             |
| SRR32060326 | N032    | AC    | control   | BC   | 0.956891                           | True             |
| SRR32060336 | N072    | TOL   | chili     | BC   | 0.9461                             | True             |
| SRR32060341 | N064    | TOL   | chili     | BC   | 0.757099                           | True             |
| SRR32060215 | N015    | AC    | control   | BC   | 0.954963                           | True             |
| SRR32060216 | N022    | TOL   | chili     | BC   | 0.891842                           | True             |
| SRR32060234 | N159    | AC    | control   | BC   | 0.981205                           | True             |
| SRR32060240 | N149    | AC    | control   | BC   | 0.968845                           | True             |

## QC exclusions

| sample      | Group | protein_coding_annotation_fraction |
| ----------- | ----- | ---------------------------------- |
| SRR32060276 | TOL   | 0.0872549                          |

## Primary result

| analysis_set                                          | metric                      | n_permutations | observed_delta_TOL_minus_AC | null_mean_delta | null_sd_delta | null_ci2_5 | null_ci97_5 | empirical_p_greater | empirical_p_less | empirical_p_two_sided | empirical_percentile |
| ----------------------------------------------------- | --------------------------- | -------------- | --------------------------- | --------------- | ------------- | ---------- | ----------- | ------------------- | ---------------- | --------------------- | -------------------- |
| all_BC_samples                                        | strict_amyloid_TPM_fraction | 10000          | 0.0677449                   | -0.018387       | 0.0403097     | -0.0959681 | 0.0607285   | 0.0178982           | 0.982202         | 0.0357964             | 0.9822               |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | strict_amyloid_TPM_fraction | 10000          | 0.0534392                   | -0.0164675      | 0.0415452     | -0.0964671 | 0.0669848   | 0.049595            | 0.950505         | 0.0991901             | 0.9505               |

## Observed endpoint statistics

| analysis_set                                          | metric                                       | n_AC | n_TOL | mean_AC   | mean_TOL   | median_AC | median_TOL | observed_delta_TOL_minus_AC | cliffs_delta_TOL_vs_AC | bootstrap_delta_ci2_5 | bootstrap_delta_median | bootstrap_delta_ci97_5 |
| ----------------------------------------------------- | -------------------------------------------- | ---- | ----- | --------- | ---------- | --------- | ---------- | --------------------------- | ---------------------- | --------------------- | ---------------------- | ---------------------- |
| all_BC_samples                                        | strict_amyloid_TPM_fraction                  | 6    | 5     | 0.0786091 | 0.146354   | 0.0887409 | 0.150604   | 0.0677449                   | 0.6                    | 0.00848347            | 0.0674575              | 0.129649               |
| all_BC_samples                                        | strict_amyloid_TPM                           | 6    | 5     | 74691.1   | 95301      | 86089.8   | 93586      | 20609.9                     | 0.2                    | -40900.3              | 20498.5                | 80237.3                |
| all_BC_samples                                        | n_strict_amyloid_expressed_TPM_ge_1          | 6    | 5     | 1214.83   | 1663.2     | 1594.5    | 1974       | 448.367                     | 0.4                    | -528.268              | 450.017                | 1415                   |
| all_BC_samples                                        | strict_weighted_index_sum                    | 6    | 5     | -12897.5  | -18912.6   | -12644.6  | -18088.8   | -6015.08                    | -0.2                   | -21831.2              | -6030.2                | 10476.2                |
| all_BC_samples                                        | continuous_expression_weighted_amyloid_index | 6    | 5     | -0.131948 | -0.107211  | -0.179626 | -0.103945  | 0.0247374                   | 0.266667               | -0.126301             | 0.0257323              | 0.171961               |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | strict_amyloid_TPM_fraction                  | 6    | 4     | 0.0786091 | 0.132048   | 0.0887409 | 0.12777    | 0.0534392                   | 0.5                    | -0.00540187           | 0.0531054              | 0.114945               |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | strict_amyloid_TPM                           | 6    | 4     | 74691.1   | 114685     | 86089.8   | 115944     | 39994.4                     | 0.333333               | -12632.8              | 39411.8                | 92558.3                |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | n_strict_amyloid_expressed_TPM_ge_1          | 6    | 4     | 1214.83   | 2027.75    | 1594.5    | 2070.5     | 812.917                     | 0.583333               | 123.156               | 800.458                | 1558.25                |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | strict_weighted_index_sum                    | 6    | 4     | -12897.5  | -24281.8   | -12644.6  | -25455.7   | -11384.3                    | -0.5                   | -26255.1              | -11384.3               | 4022.18                |
| coverage_QC_protein_coding_annotation_fraction_ge_0_5 | continuous_expression_weighted_amyloid_index | 6    | 4     | -0.131948 | -0.0824988 | -0.179626 | -0.0860999 | 0.0494495                   | 0.333333               | -0.107803             | 0.0507536              | 0.203638               |

## Interpretation

In the coverage-QC subset, the primary burden is higher in TOL/chili and reaches a borderline one-sided empirical permutation result. Because the two-sided permutation p-value is weaker and the bootstrap CI is broad, this is discovery-level support for amyloidogenic specificity, not validation and not a biomarker claim.

## Limitations

- Small baseline sample size limits inference; effect size and direction stability are more informative than p-value alone.
- The permutation null preserves expression and coarse protein/expression bins, but bin choices are still analytical assumptions.
- Amyloid annotations are computational predictions and remain discovery-stage features.
- The low-coverage sample SRR32060276 materially affects the all-sample analysis, so the coverage-QC subset should be emphasized.

## Next best step

Validate the direction and top contributor panel in an independent ChILI baseline cohort or, if unavailable, use a pre-registered internal resampling analysis with fixed bins/endpoints and orthogonal amyloid predictors.
