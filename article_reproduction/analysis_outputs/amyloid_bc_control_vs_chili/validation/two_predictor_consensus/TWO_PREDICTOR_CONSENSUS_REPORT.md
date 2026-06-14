# Two-predictor validation consensus

Consensus was recomputed using only `AMYPred_Pred` and `AmyloGramPy_Pred`; `AmyloGram_Pred` was ignored.

## Rule

- `Amyloid`: both predictors are `Amyloid`.
- `Non-Amyloid`: both predictors are `Non-Amyloid`.
- `Discordant`: both predictors are present but disagree.
- `Partial`: exactly one of the two predictors is present.

## Raw Protein Rows

| Consensus_2pred   |      n |   fraction |
|:------------------|-------:|-----------:|
| Amyloid           | 124539 |  0.636309  |
| Discordant        |  65267 |  0.33347   |
| Non-Amyloid       |   2149 |  0.0109799 |
| Partial           |   3766 |  0.0192417 |

## GENCODE Protein-Coding Collapsed Transcripts

| Consensus_2pred_collapsed   |     n |   fraction |
|:----------------------------|------:|-----------:|
| Amyloid                     | 75356 | 0.661284   |
| Discordant                  | 36216 | 0.317812   |
| Non-Amyloid                 |   731 | 0.00641487 |
| Partial                     |  1651 | 0.0144883  |

## AC vs TOL Statistics

| metric                          |   n_AC |   n_TOL |         mean_AC |       mean_TOL |       median_AC |     median_TOL |   delta_TOL_minus_AC |   mannwhitney_one_sided_TOL_gt_AC_p |   cliffs_delta_TOL_vs_AC |   bootstrap_delta_ci2_5 |   bootstrap_delta_median |   bootstrap_delta_ci97_5 |
|:--------------------------------|-------:|--------:|----------------:|---------------:|----------------:|---------------:|---------------------:|------------------------------------:|-------------------------:|------------------------:|-------------------------:|-------------------------:|
| Amyloid_TPM_fraction            |      4 |       4 |      0.727944   |      0.69034   |      0.723682   |      0.676107  |          -0.0376039  |                            0.828571 |                   -0.375 |             -0.116485   |               -0.0377012 |               0.0451498  |
| Amyloid_TPM                     |      4 |       4 | 618319          | 599186         | 618239          | 598908         |      -19133.2        |                            0.757143 |                   -0.25  |         -91010.8        |           -19133.2       |           51642.1        |
| n_Amyloid_expressed_TPM_ge_1    |      4 |       4 |   9963          |   8876         |   9584          |   9124.5       |       -1087          |                            0.9      |                   -0.5   |          -2335.06       |            -1041.75      |              21.7875     |
| Partial_TPM_fraction            |      4 |       4 |      0.00994112 |      0.0133278 |      0.00838591 |      0.0130077 |           0.00338671 |                            0.1      |                    0.625 |             -0.00017712 |                0.0035302 |               0.00648992 |
| Partial_TPM                     |      4 |       4 |   8487.24       |  11608         |   7071.32       |  11054.3       |        3120.76       |                            0.1      |                    0.625 |           -185.051      |             3156.37      |            6131.9        |
| n_Partial_expressed_TPM_ge_1    |      4 |       4 |    214          |    198.75      |    210.5        |    209.5       |         -15.25       |                            0.557725 |                    0     |            -92          |              -15.25      |              60          |
| Discordant_TPM_fraction         |      4 |       4 |      0.1466     |      0.160717  |      0.14884    |      0.138047  |           0.0141168  |                            0.557143 |                    0     |             -0.021962   |                0.0130117 |               0.0649057  |
| Discordant_TPM                  |      4 |       4 | 124740          | 140668         | 125322          | 117416         |       15928.3        |                            0.442857 |                    0.125 |         -19866          |            14835.3       |           67054.4        |
| n_Discordant_expressed_TPM_ge_1 |      4 |       4 |   4787          |   4267         |   4613.5        |   4368         |        -520          |                            0.942857 |                   -0.625 |          -1082.25       |             -507         |             -27          |
