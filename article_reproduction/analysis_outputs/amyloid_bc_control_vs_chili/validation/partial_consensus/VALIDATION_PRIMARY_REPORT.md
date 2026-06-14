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
Primary-positive amyloid consensus label: `Partial`.

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

## Primary Partial-consensus burden statistics

| metric                                 |   n_AC |   n_TOL |       mean_AC |      mean_TOL |     median_AC |    median_TOL |   observed_delta_TOL_minus_AC |   mannwhitney_one_sided_TOL_gt_AC_p |   cliffs_delta_TOL_vs_AC |   rank_biserial_TOL_vs_AC |   bootstrap_delta_ci2_5 |   bootstrap_delta_median |   bootstrap_delta_ci97_5 |   sample_label_permutation_one_sided_p |
|:---------------------------------------|-------:|--------:|--------------:|--------------:|--------------:|--------------:|------------------------------:|------------------------------------:|-------------------------:|--------------------------:|------------------------:|-------------------------:|-------------------------:|---------------------------------------:|
| partial_consensus_TPM_fraction         |      4 |       4 |      0.744164 |      0.712364 |      0.743103 |      0.701588 |                    -0.0318001 |                            0.828571 |                   -0.375 |                    -0.375 |               -0.113369 |               -0.0322591 |                0.0504841 |                               0.736326 |
| partial_consensus_TPM                  |      4 |       4 | 632167        | 618462        | 634935        | 621669        |                -13704.7       |                            0.757143 |                   -0.25  |                    -0.25  |           -86050.2      |           -13704.7       |            58526.5       |                               0.629837 |
| n_partial_consensus_expressed_TPM_ge_1 |      4 |       4 |  10205.5      |   9090.75     |   9825        |   9376        |                 -1114.75      |                            0.9      |                   -0.5   |                    -0.5   |            -2423        |            -1086.75      |                9.025     |                               0.89811  |

## Amyloid label QC

| sample      | Consensus_collapsed   |   n_transcripts |   n_primary_positive | Group   | subject   |
|:------------|:----------------------|----------------:|---------------------:|:--------|:----------|
| SRR32060218 | Discordant            |            5467 |                    0 | AC      | N013      |
| SRR32060218 | Partial               |           11714 |                11714 | AC      | N013      |
| SRR32060222 | Discordant            |            4540 |                    0 | AC      | N010      |
| SRR32060222 | Partial               |            9458 |                 9458 | AC      | N010      |
| SRR32060224 | Discordant            |            4725 |                    0 | AC      | N007      |
| SRR32060224 | Partial               |           10080 |                10080 | AC      | N007      |
| SRR32060226 | Discordant            |            4551 |                    0 | AC      | N005      |
| SRR32060226 | Partial               |            9570 |                 9570 | AC      | N005      |
| SRR32060252 | Discordant            |            3654 |                    0 | TOL     | N162      |
| SRR32060252 | Partial               |            7718 |                 7718 | TOL     | N162      |
| SRR32060255 | Discordant            |            4530 |                    0 | TOL     | N043      |
| SRR32060255 | Partial               |            9777 |                 9777 | TOL     | N043      |
| SRR32060299 | Discordant            |            4734 |                    0 | TOL     | N035      |
| SRR32060299 | Partial               |            9893 |                 9893 | TOL     | N035      |
| SRR32060304 | Discordant            |            4234 |                    0 | TOL     | N111      |
| SRR32060304 | Partial               |            8975 |                 8975 | TOL     | N111      |

## Counterfactual amyloid-label permutation

| metric                         |   n_permutations |   observed_delta_TOL_minus_AC |   empirical_p_TOL_minus_AC_ge_observed |   observed_percentile_in_null |   null_mean |   null_ci2_5 |   null_ci97_5 |
|:-------------------------------|-----------------:|------------------------------:|---------------------------------------:|------------------------------:|------------:|-------------:|--------------:|
| partial_consensus_TPM_fraction |            10000 |                    -0.0318001 |                                0.60414 |                         39.59 |  -0.0128635 |    -0.152059 |      0.128572 |

## Leave-one-sample-out stability

| left_out_sample   | left_out_subject   | metric                                 |   n_AC |   n_TOL |   delta_TOL_minus_AC | direction_TOL_gt_AC   |
|:------------------|:-------------------|:---------------------------------------|-------:|--------:|---------------------:|:----------------------|
| SRR32060218       | N013               | partial_consensus_TPM_fraction         |      3 |       4 |          -0.0411629  | False                 |
| SRR32060218       | N013               | partial_consensus_TPM                  |      3 |       4 |      -16127          | False                 |
| SRR32060218       | N013               | n_partial_consensus_expressed_TPM_ge_1 |      3 |       4 |        -611.917      | False                 |
| SRR32060222       | N010               | partial_consensus_TPM_fraction         |      3 |       4 |          -0.00797359 | False                 |
| SRR32060222       | N010               | partial_consensus_TPM                  |      3 |       4 |        4919.65       | True                  |
| SRR32060222       | N010               | n_partial_consensus_expressed_TPM_ge_1 |      3 |       4 |       -1363.92       | False                 |
| SRR32060224       | N007               | partial_consensus_TPM_fraction         |      3 |       4 |          -0.0549194  | False                 |
| SRR32060224       | N007               | partial_consensus_TPM                  |      3 |       4 |      -34174.7        | False                 |
| SRR32060224       | N007               | n_partial_consensus_expressed_TPM_ge_1 |      3 |       4 |       -1156.58       | False                 |
| SRR32060226       | N005               | partial_consensus_TPM_fraction         |      3 |       4 |          -0.0231444  | False                 |
| SRR32060226       | N005               | partial_consensus_TPM                  |      3 |       4 |       -9436.95       | False                 |
| SRR32060226       | N005               | n_partial_consensus_expressed_TPM_ge_1 |      3 |       4 |       -1326.58       | False                 |
| SRR32060252       | N162               | partial_consensus_TPM_fraction         |      4 |       3 |          -0.0392255  | False                 |
| SRR32060252       | N162               | partial_consensus_TPM                  |      4 |       3 |      -32216.5        | False                 |
| SRR32060252       | N162               | n_partial_consensus_expressed_TPM_ge_1 |      4 |       3 |        -657.167      | False                 |
| SRR32060255       | N043               | partial_consensus_TPM_fraction         |      4 |       3 |          -0.0171904  | False                 |
| SRR32060255       | N043               | partial_consensus_TPM                  |      4 |       3 |        2669.11       | True                  |
| SRR32060255       | N043               | n_partial_consensus_expressed_TPM_ge_1 |      4 |       3 |       -1343.5        | False                 |
| SRR32060299       | N035               | partial_consensus_TPM_fraction         |      4 |       3 |          -0.0621776  | False                 |
| SRR32060299       | N035               | partial_consensus_TPM                  |      4 |       3 |      -36064.6        | False                 |

## Primary endpoint interpretation

For partial_consensus_TPM_fraction, observed delta TOL minus AC was -0.0318001; one-sided Mann-Whitney p=0.828571; Cliff's delta=-0.375.
