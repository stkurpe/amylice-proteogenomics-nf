# Seven Signature Panel

## What was done
Recomputed one primary baseline comparison: `TOL/chili` versus `AC/control`. Low-coverage samples were excluded before all statistics and plots. The expression universe was explicitly filtered to GENCODE `biotype == protein_coding` transcripts before scoring or GSEA. Protein-coding transcript-level TPM was summed once per gene symbol to avoid double-counting transcript variants; scores use `log2(TPM + 1)`, per-gene z-scores across baseline samples, then the mean z-score over detected signature genes.

The protein-coding gene universe used for scoring/ranking is saved in `protein_coding_expression_gene_universe.tsv`.

The UPR score uses only the 36-gene Core UPR signature supplied in the request, not the 267-gene UPR regulon.

Robustness methods: ssGSEA and GSVA were run through `gseapy` when available; singscore was computed as rank-normalized mean signature rank.

Method status: `{"ssGSEA": "ok", "GSVA": "ok", "singscore": "ok"}`.

GSEA status: `{"TOL_vs_AC": "ok"}`.

## Samples used
Analyzed baseline samples after coverage QC: 10 total; AC n=6, TOL n=4.

| sample      | Group   | subject   | treatment              |   batch | sex    |   AGE |   protein_coding_annotation_fraction |
|:------------|:--------|:----------|:-----------------------|--------:|:-------|------:|-------------------------------------:|
| SRR32060215 | AC      | N015      | Pembrolizumab          |       1 | male   |    73 |                             0.94525  |
| SRR32060326 | AC      | N032      | Ipilimumab + Nivolumab |       2 | female |    55 |                             0.951026 |
| SRR32060312 | AC      | N106      | Ipilimumab + Nivolumab |       3 | female |    49 |                             0.900663 |
| SRR32060240 | AC      | N149      | Ipilimumab + Nivolumab |       5 | male   |    76 |                             0.964026 |
| SRR32060234 | AC      | N159      | Ipilimumab + Nivolumab |       5 | female |    56 |                             0.977679 |
| SRR32060309 | AC      | N87       | Ipilimumab + Nivolumab |       3 | male   |    63 |                             0.758511 |
| SRR32060216 | TOL     | N022      | Ipilimumab + Nivolumab |       1 | male   |    67 |                             0.884931 |
| SRR32060341 | TOL     | N064      | Ipilimumab + Nivolumab |       2 | male   |    72 |                             0.747125 |
| SRR32060336 | TOL     | N072      | Ipilimumab + Nivolumab |       2 | male   |    71 |                             0.930256 |
| SRR32060306 | TOL     | N101      | Ipilimumab + Nivolumab |       4 | male   |    50 |                             0.92886  |

## Excluded samples
Excluded samples have `protein_coding_annotation_fraction < 0.5` and are not used in any comparison, model, GSEA, or plot.

None.

## Gene coverage
| signature          |   signature_size |   genes_detected |   genes_missing | detected_gene_symbols                                                                                                                                                                                                             |   missing_gene_symbols |
|:-------------------|-----------------:|-----------------:|----------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------:|
| Amyloid_score      |               23 |               23 |               0 | AIMP1,AP2B1,CLEC2D,ILF2,LGMN,MRPS21,NLRC5,PIK3C2A,PSMD8,RARB,S100A8,SEPTIN7,SHE,SLAIN2,SNX10,STK17B,STMP1,TACC1,TIA1,TIMM10,YLPM1,ZNF580,ZNF672                                                                                   |                    nan |
| UPR_score          |               36 |               36 |               0 | ASNS,ATF3,ATF4,BHLHA15,CTH,DDIT3,DERL2,DERL3,DNAJB11,DNAJB9,DNAJC3,EDEM1,EIF2AK3,EXOSC1,FICD,GFPT1,HERPUD1,HSPA5,HYOU1,PDIA4,PPP1R15A,PREB,SEC31A,SEC61B,SEC61G,SEL1L,SELENOS,SERP1,SYVN1,TRIB3,TRIM25,UBE2J1,VCP,WFS1,WIPI1,XBP1 |                    nan |
| Inflammatory_score |               27 |               27 |               0 | CCL5,CCR5,CD274,CD3D,CD3E,CD8A,CIITA,CTLA4,CXCL10,CXCL11,CXCL13,CXCL9,GZMA,GZMB,HLA-DRA,HLA-DRB1,HLA-E,IDO1,IL2RG,ITGAL,LAG3,NKG7,PDCD1,PRF1,PTPRC,STAT1,TAGAP                                                                    |                    nan |
| IFNG_score         |                4 |                4 |               0 | CXCL10,CXCL11,CXCL9,STAT1                                                                                                                                                                                                         |                    nan |
| IFNA_score         |                7 |                7 |               0 | CXCL10,CXCL9,EOMES,GZMA,GZMB,IFNG,TBX21                                                                                                                                                                                           |                    nan |
| Myeloid_score      |                6 |                6 |               0 | CXCL1,CXCL2,CXCL3,CXCL8,IL6,PTGS2                                                                                                                                                                                                 |                    nan |

## Visualizations
Combined publication-style panel:

![Combined publication panel](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/seven_signature_publication_panel.png)

Signature score comparison:

![Signature score comparison](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/signature_scores_group_comparison.png)

Signature score boxplots:

![Signature score boxplots](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/signature_scores_boxplots_AC_vs_TOL.png)

Compact signature score boxplots:

![Compact signature score boxplots](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/signature_scores_boxplots_AC_vs_TOL_compact.png)

Heatmap of signature scores:

![Signature score heatmap](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/signature_score_heatmap.png)

Amyloid-immune coupling scatter plots:

![Amyloid immune coupling](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/amyloid_immune_coupling_scatter.png)

Correlation heatmap:

![Correlation heatmap](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/signature_score_correlation_heatmap.png)

Effect size forest plot:

![Effect size forest plot](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/signature_effect_size_forest.png)

Preranked GSEA dotplot:

![Preranked GSEA dotplot](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/gsea_prerank_dotplot.png)

Preranked GSEA enrichment panel:

![Preranked GSEA enrichment panel](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/gsea_prerank_enrichment_panel.png)

Amyloidogenic burden correlation heatmaps:

![AC burden correlation heatmap](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/amyloidogenic_burden_correlation_heatmap_AC.png)

![TOL burden correlation heatmap](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/amyloidogenic_burden_correlation_heatmap_TOL.png)

Significant amyloidogenic burden correlation regplots:

![AC significant burden correlation regplots](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/amyloidogenic_burden_significant_regplots_AC.png)

![TOL significant burden correlation regplots](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/amyloidogenic_burden_significant_regplots_TOL.png)

![Combined significant burden correlation regplots](/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/seven_signature_panel/amyloidogenic_burden_significant_regplots_combined.png)

## Amyloidogenic burden correlations
Burden metrics were recomputed from `ProteinCoding_BothPredictors_Q17_Q83_Class == Amyloidogenic` using `log2(TPM + 1)`. Heatmaps show a rho value only for Spearman correlations with p < 0.05.

| sample      |   n_protein_coding_classified_transcripts |   n_amyloidogenic_transcripts |   Amyloidogenic_Burden |   Amyloidogenic_Fraction |   KDmax_weighted_Amyloidogenic_Fraction | Group   | subject   |   protein_coding_annotation_fraction |   UPR_score |   Inflammatory_score |   IFNG_score |   IFNA_score |   Myeloid_score |
|:------------|------------------------------------------:|------------------------------:|-----------------------:|-------------------------:|----------------------------------------:|:--------|:----------|-------------------------------------:|------------:|---------------------:|-------------:|-------------:|----------------:|
| SRR32060215 |                                     17327 |                            14 |                51.8183 |               0.00104576 |                             0.000899712 | AC      | N015      |                             0.94525  |   -0.845509 |           -0.384033  |   -0.646466  |     0.386379 |       -0.602242 |
| SRR32060216 |                                     20507 |                           636 |              2101.4    |               0.0303     |                             0.0355694   | TOL     | N022      |                             0.884931 |    0.502923 |            0.631412  |    1.18192   |     0.865758 |        0.222509 |
| SRR32060234 |                                     14324 |                             0 |                 0      |               0          |                             0           | AC      | N159      |                             0.977679 |   -0.162085 |           -0.0260262 |    0.456611  |    -0.379598 |        0.138765 |
| SRR32060240 |                                     12936 |                           428 |              1443.4    |               0.0298979  |                             0.0351042   | AC      | N149      |                             0.964026 |    0.138213 |           -0.214553  |    0.367602  |     0.098197 |        0.473334 |
| SRR32060306 |                                     16275 |                           510 |              1344.95   |               0.0291272  |                             0.0342553   | TOL     | N101      |                             0.92886  |   -0.886241 |           -1.00608   |   -0.472499  |    -0.583847 |       -0.701968 |
| SRR32060309 |                                     16137 |                           518 |              1726.12   |               0.030803   |                             0.0369037   | AC      | N87       |                             0.758511 |    0.690569 |            0.4416    |   -0.505792  |     0.398249 |        0.986803 |
| SRR32060312 |                                     16552 |                           487 |              1471.14   |               0.027717   |                             0.0332267   | AC      | N106      |                             0.900663 |   -0.352995 |           -0.0888579 |   -0.813332  |    -0.570142 |       -0.625465 |
| SRR32060326 |                                     20232 |                           626 |              2085.32   |               0.0300054  |                             0.0350873   | AC      | N032      |                             0.951026 |    0.439277 |            0.0768085 |    0.0221566 |    -0.688818 |       -0.430447 |
| SRR32060336 |                                     19747 |                           627 |              1947.59   |               0.0303031  |                             0.0355283   | TOL     | N072      |                             0.930256 |    0.177077 |            0.302413  |   -0.410271  |     0.420994 |        0.383639 |
| SRR32060341 |                                     16897 |                           504 |              1725.46   |               0.0292938  |                             0.0351534   | TOL     | N064      |                             0.747125 |    0.412453 |            0.291956  |    0.129545  |     0.167859 |        0.101142 |

| Group   | metric                                | score              |   n_samples |   spearman_rho |   spearman_p_value | significant_p_lt_0_05   |
|:--------|:--------------------------------------|:-------------------|------------:|---------------:|-------------------:|:------------------------|
| AC      | Amyloidogenic_Burden                  | UPR_score          |           6 |      0.657143  |          0.156175  | False                   |
| AC      | Amyloidogenic_Burden                  | Inflammatory_score |           6 |      0.6       |          0.208     | False                   |
| AC      | Amyloidogenic_Burden                  | IFNG_score         |           6 |     -0.314286  |          0.544093  | False                   |
| AC      | Amyloidogenic_Burden                  | IFNA_score         |           6 |     -0.257143  |          0.622787  | False                   |
| AC      | Amyloidogenic_Burden                  | Myeloid_score      |           6 |      0.0857143 |          0.871743  | False                   |
| AC      | Amyloidogenic_Fraction                | UPR_score          |           6 |      0.828571  |          0.0415627 | True                    |
| AC      | Amyloidogenic_Fraction                | Inflammatory_score |           6 |      0.6       |          0.208     | False                   |
| AC      | Amyloidogenic_Fraction                | IFNG_score         |           6 |     -0.142857  |          0.787172  | False                   |
| AC      | Amyloidogenic_Fraction                | IFNA_score         |           6 |      0.142857  |          0.787172  | False                   |
| AC      | Amyloidogenic_Fraction                | Myeloid_score      |           6 |      0.485714  |          0.328723  | False                   |
| AC      | KDmax_weighted_Amyloidogenic_Fraction | UPR_score          |           6 |      0.771429  |          0.0723965 | False                   |
| AC      | KDmax_weighted_Amyloidogenic_Fraction | Inflammatory_score |           6 |      0.428571  |          0.396501  | False                   |
| AC      | KDmax_weighted_Amyloidogenic_Fraction | IFNG_score         |           6 |     -0.0857143 |          0.871743  | False                   |
| AC      | KDmax_weighted_Amyloidogenic_Fraction | IFNA_score         |           6 |      0.314286  |          0.544093  | False                   |
| AC      | KDmax_weighted_Amyloidogenic_Fraction | Myeloid_score      |           6 |      0.6       |          0.208     | False                   |
| TOL     | Amyloidogenic_Burden                  | UPR_score          |           4 |      0.8       |          0.2       | False                   |
| TOL     | Amyloidogenic_Burden                  | Inflammatory_score |           4 |      1         |          0         | True                    |
| TOL     | Amyloidogenic_Burden                  | IFNG_score         |           4 |      0.8       |          0.2       | False                   |
| TOL     | Amyloidogenic_Burden                  | IFNA_score         |           4 |      1         |          0         | True                    |
| TOL     | Amyloidogenic_Burden                  | Myeloid_score      |           4 |      0.8       |          0.2       | False                   |
| TOL     | Amyloidogenic_Fraction                | UPR_score          |           4 |      0.4       |          0.6       | False                   |
| TOL     | Amyloidogenic_Fraction                | Inflammatory_score |           4 |      0.8       |          0.2       | False                   |
| TOL     | Amyloidogenic_Fraction                | IFNG_score         |           4 |      0.4       |          0.6       | False                   |
| TOL     | Amyloidogenic_Fraction                | IFNA_score         |           4 |      0.8       |          0.2       | False                   |
| TOL     | Amyloidogenic_Fraction                | Myeloid_score      |           4 |      1         |          0         | True                    |
| TOL     | KDmax_weighted_Amyloidogenic_Fraction | UPR_score          |           4 |      0.8       |          0.2       | False                   |
| TOL     | KDmax_weighted_Amyloidogenic_Fraction | Inflammatory_score |           4 |      1         |          0         | True                    |
| TOL     | KDmax_weighted_Amyloidogenic_Fraction | IFNG_score         |           4 |      0.8       |          0.2       | False                   |
| TOL     | KDmax_weighted_Amyloidogenic_Fraction | IFNA_score         |           4 |      1         |          0         | True                    |
| TOL     | KDmax_weighted_Amyloidogenic_Fraction | Myeloid_score      |           4 |      0.8       |          0.2       | False                   |

## Preranked GSEA
Ranking metric: protein-coding gene-level Welch t-statistic for `TOL vs AC` using baseline `log2(TPM + 1)`. Positive NES means enrichment among TOL-up genes; negative NES means enrichment among AC-up genes.

| analysis_set   | signature          |   normalized_enrichment_score |   nominal_p_value |   FDR_q_value | leading_edge_genes                                                |
|:---------------|:-------------------|------------------------------:|------------------:|--------------:|:------------------------------------------------------------------|
| TOL_vs_AC      | Inflammatory_score |                      1.05177  |          0.415158 |      1        | CXCL9;CIITA;CD3E;IDO1;GZMB                                        |
| TOL_vs_AC      | Amyloid_score      |                      1.03132  |          0.436905 |      1        | PIK3C2A;ZNF672;YLPM1;LGMN;ZNF580;AP2B1;AIMP1;RARB;TACC1;SHE;PSMD8 |
| TOL_vs_AC      | IFNG_score         |                      0.984602 |          0.523741 |      0.871528 | CXCL9                                                             |
| TOL_vs_AC      | Myeloid_score      |                     -0.867303 |          0.547703 |      0.610833 | CXCL8;PTGS2                                                       |
| TOL_vs_AC      | IFNA_score         |                      0.781277 |          0.788036 |      1        | CXCL9;GZMB                                                        |
| TOL_vs_AC      | UPR_score          |                      0.703417 |          0.910615 |      0.888542 | DERL3;SEC31A;CTH;ATF3;FICD;ASNS;VCP;SELENOS;WIPI1                 |

## Group comparison results
Positive delta/effect sizes mean TOL > AC.

| signature          |   n_AC |   n_TOL |     mean_AC |   mean_TOL |   delta_mean_TOL_minus_AC |   hedges_g_TOL_vs_AC |   welch_p_value |   mannwhitney_p_value |   bootstrap_delta_mean_ci_low |   bootstrap_delta_mean_ci_high |
|:-------------------|-------:|--------:|------------:|-----------:|--------------------------:|---------------------:|----------------:|----------------------:|------------------------------:|-------------------------------:|
| Amyloid_score      |      6 |       4 | -0.0912454  | 0.115223   |                 0.206468  |            0.328069  |        0.620443 |              0.609524 |                     -0.522443 |                       0.831911 |
| UPR_score          |      6 |       4 | -0.0154216  | 0.0515529  |                 0.0669745 |            0.102573  |        0.870279 |              0.914286 |                     -0.655265 |                       0.693431 |
| Inflammatory_score |      6 |       4 | -0.0325103  | 0.0549245  |                 0.0874347 |            0.159057  |        0.830518 |              0.47619  |                     -0.632316 |                       0.631164 |
| IFNG_score         |      6 |       4 | -0.186537   | 0.107173   |                 0.29371   |            0.417507  |        0.536018 |              0.47619  |                     -0.404701 |                       1.10044  |
| IFNA_score         |      6 |       4 | -0.125956   | 0.217691   |                 0.343647  |            0.582407  |        0.382459 |              0.352381 |                     -0.307816 |                       0.94681  |
| Myeloid_score      |      6 |       4 | -0.00987532 | 0.00133055 |                 0.0112059 |            0.0169417 |        0.975999 |              0.914286 |                     -0.64412  |                       0.602594 |

## Amyloid score results
Amyloid score delta TOL-AC = 0.206, Hedges g = 0.328, Welch p = 0.62, bootstrap CI = [-0.522, 0.832].

## UPR results
UPR score delta TOL-AC = 0.067, Hedges g = 0.103, Welch p = 0.87, bootstrap CI = [-0.655, 0.693].

## Immune signature results
Inflammatory/IFN/myeloid signatures are reported in the comparison table above.

## Amyloid-Inflammation coupling
Inflammatory interaction beta = 0.503, p = 0.0661. A positive interaction would support stronger TOL coupling.

## Amyloid-UPR coupling
UPR interaction beta = -0.0178, p = 0.959.

Amyloid-UPR correlation is higher in TOL than AC directionally, but this difference is not statistically secure because bootstrap CI crosses zero.

Amyloid-IFNG correlation is also higher in TOL directionally, but again not statistically secure because bootstrap CI crosses zero.

## Spearman coupling correlations
| analysis_set   | immune_score       | Group        |   n_samples |   spearman_rho |   spearman_p_value |   correlation_difference_TOL_minus_AC |   bootstrap_ci_low |   bootstrap_ci_high |
|:---------------|:-------------------|:-------------|------------:|---------------:|-------------------:|--------------------------------------:|-------------------:|--------------------:|
| TOL_vs_AC      | UPR_score          | AC           |           6 |      0.771429  |         0.0723965  |                           nan         |          nan       |          nan        |
| TOL_vs_AC      | UPR_score          | TOL          |           4 |      0.8       |         0.2        |                           nan         |          nan       |          nan        |
| TOL_vs_AC      | UPR_score          | TOL_minus_AC |          10 |      0.0285714 |       nan          |                             0.0285714 |           -1.82353 |            1.09091  |
| TOL_vs_AC      | Inflammatory_score | AC           |           6 |      0.942857  |         0.00480466 |                           nan         |          nan       |          nan        |
| TOL_vs_AC      | Inflammatory_score | TOL          |           4 |      0.4       |         0.6        |                           nan         |          nan       |          nan        |
| TOL_vs_AC      | Inflammatory_score | TOL_minus_AC |          10 |     -0.542857  |       nan          |                            -0.542857  |           -2       |            0.387097 |
| TOL_vs_AC      | IFNG_score         | AC           |           6 |     -0.0857143 |         0.871743   |                           nan         |          nan       |          nan        |
| TOL_vs_AC      | IFNG_score         | TOL          |           4 |      0.8       |         0.2        |                           nan         |          nan       |          nan        |
| TOL_vs_AC      | IFNG_score         | TOL_minus_AC |          10 |      0.885714  |       nan          |                             0.885714  |           -1       |            1.92     |
| TOL_vs_AC      | IFNA_score         | AC           |           6 |     -0.142857  |         0.787172   |                           nan         |          nan       |          nan        |
| TOL_vs_AC      | IFNA_score         | TOL          |           4 |      0.4       |         0.6        |                           nan         |          nan       |          nan        |
| TOL_vs_AC      | IFNA_score         | TOL_minus_AC |          10 |      0.542857  |       nan          |                             0.542857  |           -1.45455 |            2        |
| TOL_vs_AC      | Myeloid_score      | AC           |           6 |      0.314286  |         0.544093   |                           nan         |          nan       |          nan        |
| TOL_vs_AC      | Myeloid_score      | TOL          |           4 |      0.2       |         0.8        |                           nan         |          nan       |          nan        |
| TOL_vs_AC      | Myeloid_score      | TOL_minus_AC |          10 |     -0.114286  |       nan          |                            -0.114286  |           -1.93548 |            1.59596  |

## Interpretation
The analysis is exploratory and underpowered. Interpret direction, effect size, and bootstrap intervals together; do not treat nominal p-values as biomarker-level evidence.

If Amyloid_score, immune signatures, or interaction terms are not consistently positive with intervals excluding zero, the stated hypothesis is not statistically confirmed in this sample set.

## Limitations
Very small sample size after QC, blood bulk RNA-seq compositional confounding, tied ranks in preranked GSEA, and no multiple-testing-powered biomarker validation. Group assignment is AC/control versus TOL/chili at baseline only.

## Next best step
Validate the score directions in an independent cohort or with subject-level/cell-composition-adjusted models, then test whether the amyloid-expression signature adds signal beyond inflammatory cell-state markers.
