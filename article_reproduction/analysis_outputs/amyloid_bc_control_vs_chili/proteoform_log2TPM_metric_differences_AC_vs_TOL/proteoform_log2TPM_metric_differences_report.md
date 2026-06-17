# Proteoform log2(TPM+1) metric differences: AC vs TOL

- Zarr: `project.zarr`
- Working class column: `ProteinCoding_BothPredictors_Q17_Q83_Class`
- Amyloidogenic class: `Amyloidogenic`
- Denominator: all rows collapsed to transcripts with `is_gencode_protein_coding == True`.
- TPM transform: `log2(TPM + 1)` for every burden/fraction endpoint.
- Zarr class attrs: `{'column': 'ProteinCoding_BothPredictors_Q17_Q83_Class', 'filter': 'GENCODE biotype == protein_coding', 'gencode_annotation': 'analysis_outputs/amyloid_bc_control_vs_chili/rank_based_amyloid_burden_signature/gencode_v48_transcript_annotation.tsv', 'source_class': 'BothPredictors_Q17_Q83_Class'}`

## Sample-level endpoints

| endpoint                | metric                                        |   n_AC |   n_TOL |      mean_AC |     mean_TOL |    median_AC |   median_TOL |   delta_AC_minus_TOL |   mannwhitney_p_value |
|:------------------------|:----------------------------------------------|-------:|--------:|-------------:|-------------:|-------------:|-------------:|---------------------:|----------------------:|
| Burden                  | amyloidogenic_burden_log2TPM                  |      6 |       5 | 1105.59      | 1402.46      | 1430.04      | 1667.44      |        -296.87       |              0.662338 |
| Fraction                | amyloidogenic_fraction_log2TPM                |      6 |       5 |    0.0229601 |    0.0317211 |    0.0330309 |    0.0339165 |          -0.00876104 |              0.930736 |
| KDmax-weighted Fraction | kdmax_weighted_amyloidogenic_fraction_log2TPM |      6 |       5 |    0.0264373 |    0.037942  |    0.0381453 |    0.0391407 |          -0.0115047  |              0.930736 |

## Top proteoform candidates

| metric                                       | gene_name   | Transcript_ID_clean   |   n_AC |   n_TOL |   delta_AC_minus_TOL |   cohen_d |     p_value |         FDR | direction     |
|:---------------------------------------------|:------------|:----------------------|-------:|--------:|---------------------:|----------:|------------:|------------:|:--------------|
| burden_contribution_log2TPM                  | OR10J1      | ENST00000642080       |      3 |       4 |          1.02377     |   3.93466 | 2.26913e-06 | 0.000730661 | higher_in_AC  |
| kdmax_weighted_fraction_contribution_log2TPM | TP53BP1     | ENST00000382044       |      3 |       4 |         -1.60906e-05 |  -3.67878 | 0.000126152 | 0.0406209   | higher_in_TOL |
| burden_contribution_log2TPM                  | SRPK1       | ENST00000373825       |      3 |       4 |          1.50416     |   2.94633 | 0.000387821 | 0.0624392   | higher_in_AC  |
| fraction_contribution_log2TPM                | TP53BP1     | ENST00000382044       |      3 |       4 |         -1.38533e-05 |  -3.4444  | 0.000342404 | 0.110254    | higher_in_TOL |
| burden_contribution_log2TPM                  | MAGEA10     | ENST00000370323       |      4 |       3 |         -0.616931    |  -2.31054 | 0.00791199  | 0.800468    | higher_in_TOL |
| burden_contribution_log2TPM                  | PPIG        | ENST00000676508       |      4 |       3 |          0.722467    |   2.24327 | 0.0099437   | 0.800468    | higher_in_AC  |
| fraction_contribution_log2TPM                | CSGALNACT2  | ENST00000374466       |      3 |       3 |          2.42758e-05 |   2.3949  | 0.016625    | 0.842632    | higher_in_AC  |
| fraction_contribution_log2TPM                | STX3        | ENST00000337979       |      3 |       3 |          2.55111e-05 |   2.24827 | 0.0245591   | 0.842632    | higher_in_AC  |
| fraction_contribution_log2TPM                | OR6C74      | ENST00000343399       |      3 |       3 |          1.58858e-05 |   2.21282 | 0.0269098   | 0.842632    | higher_in_AC  |
| fraction_contribution_log2TPM                | PPIG        | ENST00000676508       |      4 |       3 |          1.92436e-05 |   2.16608 | 0.012327    | 0.842632    | higher_in_AC  |
| fraction_contribution_log2TPM                | SFMBT1      | ENST00000394752       |      3 |       4 |         -1.39068e-05 |  -2.11892 | 0.010113    | 0.842632    | higher_in_TOL |
| fraction_contribution_log2TPM                | SRPK1       | ENST00000373825       |      3 |       4 |          3.14239e-05 |   2.06358 | 0.0311917   | 0.842632    | higher_in_AC  |
| fraction_contribution_log2TPM                | RMC1        | ENST00000269221       |      4 |       4 |         -1.22073e-05 |  -2.01417 | 0.0136308   | 0.842632    | higher_in_TOL |
| fraction_contribution_log2TPM                | LRRK2       | ENST00000298910       |      4 |       3 |          2.10482e-05 |   1.9984  | 0.0191926   | 0.842632    | higher_in_AC  |
| fraction_contribution_log2TPM                | SOS2        | ENST00000543680       |      4 |       4 |          2.75008e-05 |   1.89153 | 0.0205233   | 0.842632    | higher_in_AC  |
| fraction_contribution_log2TPM                | RRAGA       | ENST00000380527       |      4 |       4 |          1.74357e-05 |   1.84553 | 0.0238023   | 0.842632    | higher_in_AC  |
| fraction_contribution_log2TPM                | GNAQ        | ENST00000286548       |      4 |       4 |          2.21311e-05 |   1.75705 | 0.0314024   | 0.842632    | higher_in_AC  |
| fraction_contribution_log2TPM                | INTS3       | ENST00000318967       |      4 |       3 |          1.21719e-05 |   1.76223 | 0.0461892   | 0.916543    | higher_in_AC  |
| fraction_contribution_log2TPM                | PSTPIP1     | ENST00000558870       |      3 |       3 |          2.53958e-05 |   1.73126 | 0.0834048   | 0.916543    | higher_in_AC  |
| fraction_contribution_log2TPM                | OR4K17      | ENST00000641386       |      3 |       3 |          1.25054e-05 |   1.68649 | 0.0917006   | 0.916543    | higher_in_AC  |

## QC

| sample      |   n_protein_coding_rows |   n_collapsed_protein_coding_transcripts |   n_amyloidogenic_transcripts |   all_protein_coding_log2TPM_sum |   all_protein_coding_KDmax_log2TPM_sum |   n_missing_expression_after_merge |   n_missing_KDmax |
|:------------|------------------------:|-----------------------------------------:|------------------------------:|---------------------------------:|---------------------------------------:|-----------------------------------:|------------------:|
| SRR32060215 |                   23630 |                                    14630 |                            14 |                         42329    |                               105468   |                                  0 |                 3 |
| SRR32060216 |                   26027 |                                    17138 |                           612 |                         58581.5  |                               147103   |                                  0 |                28 |
| SRR32060234 |                   24126 |                                    12212 |                             0 |                         42168.5  |                               105648   |                                  0 |                 1 |
| SRR32060240 |                   21817 |                                    11151 |                           428 |                         41807.1  |                               105423   |                                  0 |                 3 |
| SRR32060276 |                    2934 |                                     2352 |                            66 |                          8200.82 |                                15970.4 |                                  0 |               478 |
| SRR32060306 |                   22753 |                                    13736 |                           490 |                         39449.3  |                                98508.1 |                                  0 |                22 |
| SRR32060309 |                   15977 |                                    13409 |                           507 |                         47092.8  |                               115983   |                                  0 |               157 |
| SRR32060312 |                   16231 |                                    13872 |                           475 |                         44922    |                               112436   |                                  0 |                45 |
| SRR32060326 |                   30436 |                                    16919 |                           615 |                         58695.6  |                               148505   |                                  0 |                 1 |
| SRR32060336 |                   26698 |                                    16555 |                           606 |                         54478.1  |                               136781   |                                  0 |                23 |
| SRR32060341 |                   17025 |                                    13961 |                           492 |                         49162.9  |                               121148   |                                  0 |               128 |

- Sample endpoint plot: `sample_level_log2TPM_metric_boxplots_AC_vs_TOL.png`
- Proteoform forest plot: `top_proteoform_metric_forest_grid_AC_vs_TOL.png`