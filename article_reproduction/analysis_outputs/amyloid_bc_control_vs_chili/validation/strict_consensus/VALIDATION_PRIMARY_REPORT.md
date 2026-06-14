# Validation primary amyloidogenic burden

## Cohort

| sample      | subject   | time   | raw_Group   | Group   | treatment              |   batch |
|:------------|:----------|:-------|:------------|:--------|:-----------------------|--------:|
| SRR32060218 | N013      | BC     | control     | AC      | Pembrolizumab          |       1 |
| SRR32060222 | N010      | BC     | control     | AC      | Pembrolizumab          |       1 |
| SRR32060224 | N007      | BC     | control     | AC      | Pembrolizumab          |       1 |
| SRR32060226 | N005      | BC     | control     | AC      | Ipilimumab + Nivolumab |       1 |
| SRR32060252 | N162      | BC     | chili       | TOL     | Ipilimumab + Nivolumab |       5 |
| SRR32060255 | N043      | BC     | chili       | TOL     | Ipilimumab + Nivolumab |       1 |
| SRR32060299 | N035      | BC     | chili       | TOL     | Ipilimumab + Nivolumab |       1 |
| SRR32060304 | N111      | BC     | chili       | TOL     | Ipilimumab + Nivolumab |       4 |

## GENCODE filter

GENCODE annotation used: `/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/competitive_gene_set_test_gencode_v48_protein_coding/gencode_transcript_annotation_used.tsv`. Only `biotype == protein_coding` transcripts were retained.
Primary-positive amyloid consensus label: `Amyloid`.

## Coverage QC

| sample      | Group   | subject   |   protein_coding_annotation_fraction | keep_coverage_qc   |
|:------------|:--------|:----------|-------------------------------------:|:-------------------|
| SRR32060218 | AC      | N013      |                             0.88401  | True               |
| SRR32060222 | AC      | N010      |                             0.944431 | True               |
| SRR32060224 | AC      | N007      |                             0.835855 | True               |
| SRR32060226 | AC      | N005      |                             0.907706 | True               |
| SRR32060252 | TOL     | N162      |                             0.975603 | True               |
| SRR32060255 | TOL     | N043      |                             0.799093 | True               |
| SRR32060299 | TOL     | N035      |                             0.947415 | True               |
| SRR32060304 | TOL     | N111      |                             0.779258 | True               |

## Primary Amyloid-consensus burden statistics

| metric                              |   n_AC |   n_TOL |   mean_AC |   mean_TOL |   median_AC |   median_TOL |   observed_delta_TOL_minus_AC |   mannwhitney_one_sided_TOL_gt_AC_p |   cliffs_delta_TOL_vs_AC |   rank_biserial_TOL_vs_AC |   bootstrap_delta_ci2_5 |   bootstrap_delta_median |   bootstrap_delta_ci97_5 |   sample_label_permutation_one_sided_p |
|:------------------------------------|-------:|--------:|----------:|-----------:|------------:|-------------:|------------------------------:|------------------------------------:|-------------------------:|--------------------------:|------------------------:|-------------------------:|-------------------------:|---------------------------------------:|
| strict_amyloid_TPM_fraction         |      4 |       4 |         0 |          0 |           0 |            0 |                             0 |                                   1 |                        0 |                         0 |                       0 |                        0 |                        0 |                                      1 |
| strict_amyloid_TPM                  |      4 |       4 |         0 |          0 |           0 |            0 |                             0 |                                   1 |                        0 |                         0 |                       0 |                        0 |                        0 |                                      1 |
| n_strict_amyloid_expressed_TPM_ge_1 |      4 |       4 |         0 |          0 |           0 |            0 |                             0 |                                   1 |                        0 |                         0 |                       0 |                        0 |                        0 |                                      1 |

## Amyloid label QC

| sample      | Consensus_collapsed   |   n_transcripts |   n_primary_positive | Group   | subject   |
|:------------|:----------------------|----------------:|---------------------:|:--------|:----------|
| SRR32060218 | Discordant            |            5467 |                    0 | AC      | N013      |
| SRR32060218 | Partial               |           11714 |                    0 | AC      | N013      |
| SRR32060222 | Discordant            |            4540 |                    0 | AC      | N010      |
| SRR32060222 | Partial               |            9458 |                    0 | AC      | N010      |
| SRR32060224 | Discordant            |            4725 |                    0 | AC      | N007      |
| SRR32060224 | Partial               |           10080 |                    0 | AC      | N007      |
| SRR32060226 | Discordant            |            4551 |                    0 | AC      | N005      |
| SRR32060226 | Partial               |            9570 |                    0 | AC      | N005      |
| SRR32060252 | Discordant            |            3654 |                    0 | TOL     | N162      |
| SRR32060252 | Partial               |            7718 |                    0 | TOL     | N162      |
| SRR32060255 | Discordant            |            4530 |                    0 | TOL     | N043      |
| SRR32060255 | Partial               |            9777 |                    0 | TOL     | N043      |
| SRR32060299 | Discordant            |            4734 |                    0 | TOL     | N035      |
| SRR32060299 | Partial               |            9893 |                    0 | TOL     | N035      |
| SRR32060304 | Discordant            |            4234 |                    0 | TOL     | N111      |
| SRR32060304 | Partial               |            8975 |                    0 | TOL     | N111      |

## Counterfactual amyloid-label permutation

| metric                      |   n_permutations |   observed_delta_TOL_minus_AC |   empirical_p_TOL_minus_AC_ge_observed |   observed_percentile_in_null |   null_mean |   null_ci2_5 |   null_ci97_5 |
|:----------------------------|-----------------:|------------------------------:|---------------------------------------:|------------------------------:|------------:|-------------:|--------------:|
| strict_amyloid_TPM_fraction |            10000 |                             0 |                                      1 |                           100 |           0 |            0 |             0 |

## Leave-one-sample-out stability

| left_out_sample   | left_out_subject   | metric                              |   n_AC |   n_TOL |   delta_TOL_minus_AC | direction_TOL_gt_AC   |
|:------------------|:-------------------|:------------------------------------|-------:|--------:|---------------------:|:----------------------|
| SRR32060218       | N013               | strict_amyloid_TPM_fraction         |      3 |       4 |                    0 | False                 |
| SRR32060218       | N013               | strict_amyloid_TPM                  |      3 |       4 |                    0 | False                 |
| SRR32060218       | N013               | n_strict_amyloid_expressed_TPM_ge_1 |      3 |       4 |                    0 | False                 |
| SRR32060222       | N010               | strict_amyloid_TPM_fraction         |      3 |       4 |                    0 | False                 |
| SRR32060222       | N010               | strict_amyloid_TPM                  |      3 |       4 |                    0 | False                 |
| SRR32060222       | N010               | n_strict_amyloid_expressed_TPM_ge_1 |      3 |       4 |                    0 | False                 |
| SRR32060224       | N007               | strict_amyloid_TPM_fraction         |      3 |       4 |                    0 | False                 |
| SRR32060224       | N007               | strict_amyloid_TPM                  |      3 |       4 |                    0 | False                 |
| SRR32060224       | N007               | n_strict_amyloid_expressed_TPM_ge_1 |      3 |       4 |                    0 | False                 |
| SRR32060226       | N005               | strict_amyloid_TPM_fraction         |      3 |       4 |                    0 | False                 |
| SRR32060226       | N005               | strict_amyloid_TPM                  |      3 |       4 |                    0 | False                 |
| SRR32060226       | N005               | n_strict_amyloid_expressed_TPM_ge_1 |      3 |       4 |                    0 | False                 |
| SRR32060252       | N162               | strict_amyloid_TPM_fraction         |      4 |       3 |                    0 | False                 |
| SRR32060252       | N162               | strict_amyloid_TPM                  |      4 |       3 |                    0 | False                 |
| SRR32060252       | N162               | n_strict_amyloid_expressed_TPM_ge_1 |      4 |       3 |                    0 | False                 |
| SRR32060255       | N043               | strict_amyloid_TPM_fraction         |      4 |       3 |                    0 | False                 |
| SRR32060255       | N043               | strict_amyloid_TPM                  |      4 |       3 |                    0 | False                 |
| SRR32060255       | N043               | n_strict_amyloid_expressed_TPM_ge_1 |      4 |       3 |                    0 | False                 |
| SRR32060299       | N035               | strict_amyloid_TPM_fraction         |      4 |       3 |                    0 | False                 |
| SRR32060299       | N035               | strict_amyloid_TPM                  |      4 |       3 |                    0 | False                 |

## Primary endpoint interpretation

For strict_amyloid_TPM_fraction, observed delta TOL minus AC was 0; one-sided Mann-Whitney p=1; Cliff's delta=0.
