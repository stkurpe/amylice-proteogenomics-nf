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
| partial_consensus_TPM_fraction         |      4 |       4 |      0.893001 |      0.875342 |      0.895858 |      0.873254 |                    -0.0176582 |                            0.557143 |                      0   |                       0   |                -0.10994 |               -0.0178683 |                0.0755554 |                               0.612739 |
| partial_consensus_TPM                  |      4 |       4 | 758823        | 761151        | 765819        | 744430        |                  2328.03      |                            0.557143 |                      0   |                       0   |            -93185.1     |             2328.03      |           102001         |                               0.488151 |
| n_partial_consensus_expressed_TPM_ge_1 |      4 |       4 |  15026.2      |  13378.8      |  14463        |  13758        |                 -1647.5       |                            0.9      |                     -0.5 |                      -0.5 |             -3531.5     |            -1598.25      |               -8.46875   |                               0.89811  |

## Amyloid label QC

| sample      | Consensus_collapsed   |   n_transcripts |   n_primary_positive | Group   | subject   |
|:------------|:----------------------|----------------:|---------------------:|:--------|:----------|
| SRR32060218 | Partial               |           17181 |                17181 | AC      | N013      |
| SRR32060222 | Partial               |           13998 |                13998 | AC      | N010      |
| SRR32060224 | Partial               |           14805 |                14805 | AC      | N007      |
| SRR32060226 | Partial               |           14121 |                14121 | AC      | N005      |
| SRR32060252 | Partial               |           11372 |                11372 | TOL     | N162      |
| SRR32060255 | Partial               |           14307 |                14307 | TOL     | N043      |
| SRR32060299 | Partial               |           14627 |                14627 | TOL     | N035      |
| SRR32060304 | Partial               |           13209 |                13209 | TOL     | N111      |

## Counterfactual amyloid-label permutation

| metric                         |   n_permutations |   observed_delta_TOL_minus_AC |   empirical_p_TOL_minus_AC_ge_observed |   observed_percentile_in_null |   null_mean |   null_ci2_5 |   null_ci97_5 |
|:-------------------------------|-----------------:|------------------------------:|---------------------------------------:|------------------------------:|------------:|-------------:|--------------:|
| partial_consensus_TPM_fraction |            10000 |                    -0.0176582 |                                      1 |                             0 |  -0.0176582 |   -0.0176582 |    -0.0176582 |

## Leave-one-sample-out stability

| left_out_sample   | left_out_subject   | metric                                 |   n_AC |   n_TOL |   delta_TOL_minus_AC | direction_TOL_gt_AC   |
|:------------------|:-------------------|:---------------------------------------|-------:|--------:|---------------------:|:----------------------|
| SRR32060218       | N013               | partial_consensus_TPM_fraction         |      3 |       4 |          -0.0206549  | False                 |
| SRR32060218       | N013               | partial_consensus_TPM                  |      3 |       4 |        6537.83       | True                  |
| SRR32060218       | N013               | n_partial_consensus_expressed_TPM_ge_1 |      3 |       4 |        -929.25       | False                 |
| SRR32060222       | N010               | partial_consensus_TPM_fraction         |      3 |       4 |          -0.00051457 | False                 |
| SRR32060222       | N010               | partial_consensus_TPM                  |      3 |       4 |       14946.9        | True                  |
| SRR32060222       | N010               | n_partial_consensus_expressed_TPM_ge_1 |      3 |       4 |       -1990.25       | False                 |
| SRR32060224       | N007               | partial_consensus_TPM_fraction         |      3 |       4 |          -0.0367069  | False                 |
| SRR32060224       | N007               | partial_consensus_TPM                  |      3 |       4 |      -14955.2        | False                 |
| SRR32060224       | N007               | n_partial_consensus_expressed_TPM_ge_1 |      3 |       4 |       -1721.25       | False                 |
| SRR32060226       | N005               | partial_consensus_TPM_fraction         |      3 |       4 |          -0.0127564  | False                 |
| SRR32060226       | N005               | partial_consensus_TPM                  |      3 |       4 |        2782.58       | True                  |
| SRR32060226       | N005               | n_partial_consensus_expressed_TPM_ge_1 |      3 |       4 |       -1949.25       | False                 |
| SRR32060252       | N162               | partial_consensus_TPM_fraction         |      4 |       3 |          -0.0510783  | False                 |
| SRR32060252       | N162               | partial_consensus_TPM                  |      4 |       3 |      -42311.3        | False                 |
| SRR32060252       | N162               | n_partial_consensus_expressed_TPM_ge_1 |      4 |       3 |        -978.583      | False                 |
| SRR32060255       | N043               | partial_consensus_TPM_fraction         |      4 |       3 |           0.00775813 | True                  |
| SRR32060255       | N043               | partial_consensus_TPM                  |      4 |       3 |       29202.6        | True                  |
| SRR32060255       | N043               | n_partial_consensus_expressed_TPM_ge_1 |      4 |       3 |       -1956.92       | False                 |
| SRR32060299       | N035               | partial_consensus_TPM_fraction         |      4 |       3 |          -0.0416825  | False                 |
| SRR32060299       | N035               | partial_consensus_TPM                  |      4 |       3 |      -13399.3        | False                 |

## Primary endpoint interpretation

For partial_consensus_TPM_fraction, observed delta TOL minus AC was -0.0176582; one-sided Mann-Whitney p=0.557143; Cliff's delta=0.
