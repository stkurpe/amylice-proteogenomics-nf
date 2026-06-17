# Validation q17/q83 probability-group analysis

## Cohort

| sample      | subject   | Group   | raw_Group   | time   | treatment              |
|:------------|:----------|:--------|:------------|:-------|:-----------------------|
| SRR32060218 | N013      | AC      | control     | BC     | Pembrolizumab          |
| SRR32060222 | N010      | AC      | control     | BC     | Pembrolizumab          |
| SRR32060224 | N007      | AC      | control     | BC     | Pembrolizumab          |
| SRR32060226 | N005      | AC      | control     | BC     | Ipilimumab + Nivolumab |
| SRR32060252 | N162      | TOL     | chili       | BC     | Ipilimumab + Nivolumab |
| SRR32060255 | N043      | TOL     | chili       | BC     | Ipilimumab + Nivolumab |
| SRR32060299 | N035      | TOL     | chili       | BC     | Ipilimumab + Nivolumab |
| SRR32060304 | N111      | TOL     | chili       | BC     | Ipilimumab + Nivolumab |

## Rules

- AMYPred-FRL: prob <= 0.1206 Non-Amyloidogenic; prob >= 0.8802 Amyloidogenic; otherwise Intermediate.
- AmyloGramPy: prob <= 0.6947 Non-Amyloidogenic; prob >= 0.8898 Amyloidogenic; otherwise Intermediate.
- Consensus: both above high thresholds Amyloidogenic; both below low thresholds Non-Amyloidogenic; otherwise Intermediate.
- GENCODE filter: only `biotype == protein_coding`; annotation `/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/competitive_gene_set_test_gencode_v48_protein_coding/gencode_transcript_annotation_used.tsv`.

## Main Outputs

- `all_predictors_probability_class_burden_boxplots.png/pdf`
- `target_genes_gene_level_consensus_metric_boxplots.png/pdf`
- `target_gene_<GENE>_consensus_AC_vs_TOL_three_metric_panel.png/pdf`
- `validation_signature_boxplots_TOL_vs_AC.png/pdf`

## AC vs TOL Statistics: Sample-Level Burden Metrics

| predictor   | probability_class   | metric                  |   n_AC |   n_TOL |     median_AC |    median_TOL |       mean_AC |     mean_TOL |   delta_mean_TOL_minus_AC |   mannwhitney_p_value |   FDR_BH |
|:------------|:--------------------|:------------------------|-------:|--------:|--------------:|--------------:|--------------:|-------------:|--------------------------:|----------------------:|---------:|
| AMYPred-FRL | Amyloidogenic       | absolute_burden         |      4 |       4 | 10511         | 10003         | 10862         | 9901.3       |             -961.23       |               0.48571 |        1 |
| AMYPred-FRL | Amyloidogenic       | normalized_fraction     |      4 |       4 |     0.2303    |     0.22471   |     0.2296    |    0.22572   |               -0.0038784  |               0.11429 |        1 |
| AMYPred-FRL | Amyloidogenic       | KDmax_weighted_fraction |      4 |       4 |     0.22111   |     0.21716   |     0.22101   |    0.21764   |               -0.0033638  |               0.11429 |        1 |
| AMYPred-FRL | Non-Amyloidogenic   | absolute_burden         |      4 |       4 |  1725.2       |  1629.5       |  1759.3       | 1622.4       |             -136.93       |               0.34286 |        1 |
| AMYPred-FRL | Non-Amyloidogenic   | normalized_fraction     |      4 |       4 |     0.037342  |     0.037257  |     0.037278  |    0.037101  |               -0.00017668 |               1       |        1 |
| AMYPred-FRL | Non-Amyloidogenic   | KDmax_weighted_fraction |      4 |       4 |     0.040121  |     0.039819  |     0.039992  |    0.03969   |               -0.00030202 |               1       |        1 |
| AmyloGramPy | Amyloidogenic       | absolute_burden         |      4 |       4 |  7650.5       |  7511         |  7957.1       | 7348         |             -609.15       |               0.68571 |        1 |
| AmyloGramPy | Amyloidogenic       | normalized_fraction     |      4 |       4 |     0.16695   |     0.16824   |     0.16779   |    0.16755   |               -0.00024019 |               0.88571 |        1 |
| AmyloGramPy | Amyloidogenic       | KDmax_weighted_fraction |      4 |       4 |     0.19547   |     0.1974    |     0.19592   |    0.1963    |                0.00038411 |               0.88571 |        1 |
| AmyloGramPy | Non-Amyloidogenic   | absolute_burden         |      4 |       4 |  6960.1       |  6741.1       |  7065.8       | 6635.3       |             -430.55       |               0.68571 |        1 |
| AmyloGramPy | Non-Amyloidogenic   | normalized_fraction     |      4 |       4 |     0.15002   |     0.15137   |     0.1495    |    0.15115   |                0.0016573  |               0.88571 |        1 |
| AmyloGramPy | Non-Amyloidogenic   | KDmax_weighted_fraction |      4 |       4 |     0.095886  |     0.095329  |     0.095787  |    0.095799  |                1.2176e-05 |               0.68571 |        1 |
| Consensus   | Amyloidogenic       | absolute_burden         |      4 |       4 |  1571.2       |  1554.7       |  1655.6       | 1504.7       |             -150.84       |               0.68571 |        1 |
| Consensus   | Amyloidogenic       | normalized_fraction     |      4 |       4 |     0.034681  |     0.034506  |     0.034938  |    0.034285  |               -0.00065257 |               0.88571 |        1 |
| Consensus   | Amyloidogenic       | KDmax_weighted_fraction |      4 |       4 |     0.039942  |     0.039786  |     0.040118  |    0.039607  |               -0.00051129 |               0.88571 |        1 |
| Consensus   | Non-Amyloidogenic   | absolute_burden         |      4 |       4 |   143.66      |   135.82      |   139.61      |  138.26      |               -1.3445     |               0.88571 |        1 |
| Consensus   | Non-Amyloidogenic   | normalized_fraction     |      4 |       4 |     0.0029371 |     0.0030593 |     0.0029592 |    0.0031403 |                0.0001811  |               0.34286 |        1 |
| Consensus   | Non-Amyloidogenic   | KDmax_weighted_fraction |      4 |       4 |     0.0022022 |     0.0023822 |     0.0022089 |    0.0023847 |                0.00017579 |               0.68571 |        1 |

## AC vs TOL Statistics: Target Gene-Level Contributions

| gene_name   | metric                  |   n_AC |   n_TOL |   median_AC |   median_TOL |    mean_AC |   mean_TOL |   delta_mean_TOL_minus_AC |   mannwhitney_p_value |   FDR_BH |
|:------------|:------------------------|-------:|--------:|------------:|-------------:|-----------:|-----------:|--------------------------:|----------------------:|---------:|
| OR10J1      | absolute_burden         |      4 |       2 |  1.5689     |   1.6905     | 1.5523     | 1.6905     |                0.13819    |               0.8     |        1 |
| OR10J1      | normalized_fraction     |      4 |       2 |  3.4602e-05 |   4.0535e-05 | 3.3323e-05 | 4.0535e-05 |                7.2122e-06 |               0.8     |        1 |
| OR10J1      | KDmax_weighted_fraction |      4 |       2 |  5.1523e-05 |   5.9764e-05 | 4.9568e-05 | 5.9764e-05 |                1.0195e-05 |               0.8     |        1 |
| SRPK1       | absolute_burden         |      3 |       4 |  5.2618     |   5.2723     | 6.7682     | 5.2565     |               -1.5117     |               1       |        1 |
| SRPK1       | normalized_fraction     |      3 |       4 |  0.00012173 |   0.00012383 | 0.00014074 | 0.00012039 |               -2.0348e-05 |               0.85714 |        1 |
| SRPK1       | KDmax_weighted_fraction |      3 |       4 |  0.00013519 |   0.00013699 | 0.00014159 | 0.00013356 |               -8.028e-06  |               0.85714 |        1 |
| TP53BP1     | absolute_burden         |      3 |       2 |  2.1295     |   4.1929     | 2.2044     | 4.1929     |                1.9884     |               1       |        1 |
| TP53BP1     | normalized_fraction     |      3 |       2 |  4.3662e-05 |   8.9727e-05 | 4.5985e-05 | 8.9727e-05 |                4.3743e-05 |               1       |        1 |
| TP53BP1     | KDmax_weighted_fraction |      3 |       2 |  5.0608e-05 |   9.1936e-05 | 5.171e-05  | 9.1936e-05 |                4.0227e-05 |               1       |        1 |

## AC vs TOL Statistics: Transcriptional Signatures

| metric             |   n_AC |   n_TOL |   median_AC |   median_TOL |   mean_AC |   mean_TOL |   delta_mean_TOL_minus_AC |   mannwhitney_p_value |   FDR_BH |
|:-------------------|-------:|--------:|------------:|-------------:|----------:|-----------:|--------------------------:|----------------------:|---------:|
| UPR_score          |      4 |       4 |    0.057839 |    0.10971   |  0.050752 |  -0.050752 |                  -0.1015  |               0.88571 |  0.88571 |
| Inflammatory_score |      4 |       4 |    0.12049  |   -0.11855   |  0.11763  |  -0.11763  |                  -0.23526 |               0.68571 |  0.88571 |
| IFNG_score         |      4 |       4 |    0.063685 |   -0.0046603 |  0.15942  |  -0.15942  |                  -0.31885 |               0.48571 |  0.88571 |
| IFNA_score         |      4 |       4 |   -0.11609  |   -0.069133  |  0.051309 |  -0.051309 |                  -0.10262 |               0.88571 |  0.88571 |
| Myeloid_score      |      4 |       4 |    0.082334 |    0.14581   | -0.10667  |   0.10667  |                   0.21335 |               0.88571 |  0.88571 |
