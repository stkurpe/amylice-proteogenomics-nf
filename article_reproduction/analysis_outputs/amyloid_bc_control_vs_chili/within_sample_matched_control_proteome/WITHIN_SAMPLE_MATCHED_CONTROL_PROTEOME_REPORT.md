# Within-Sample Matched Control Proteome Report

## What was done

Recomputed the baseline (`time == "BC"`) control-vs-ChILI analysis directly from `project.zarr` using a within-sample matched-control proteome design.
For each strict amyloidogenic transcript (`Consensus_collapsed == "Amyloid"`), protein variants were first collapsed to `sample × transcript`, then one non-amyloid transcript from the same sample was selected without replacement by nearest-neighbor matching on observed `log2(TPM + 1)` and log protein length. This avoids H1/H2 double-counting.

GENCODE transcript annotation was used to restrict the matched proteome to `biotype == "protein_coding"`, so the amyloid and matched non-amyloid pools are exactly matched on gene biotype at the analysis-design level. GC/content and subcellular localization were not available in the local files used here, so the executable distance covariates were expression and protein length.

GENCODE source: `/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/rank_based_amyloid_burden_signature/gencode_v48_transcript_annotation.tsv`.

## Samples used

| sample | subject | time | raw_Group | Group | treatment | batch | AGE | sex |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SRR32060276 | N132 | BC | chili | TOL | Ipilimumab + Nivolumab | 4 | 44 | female |
| SRR32060306 | N101 | BC | chili | TOL | Ipilimumab + Nivolumab | 4 | 50 | male |
| SRR32060309 | N87 | BC | control | AC | Ipilimumab + Nivolumab | 3 | 63 | male |
| SRR32060312 | N106 | BC | control | AC | Ipilimumab + Nivolumab | 3 | 49 | female |
| SRR32060326 | N032 | BC | control | AC | Ipilimumab + Nivolumab | 2 | 55 | female |
| SRR32060336 | N072 | BC | chili | TOL | Ipilimumab + Nivolumab | 2 | 71 | male |
| SRR32060341 | N064 | BC | chili | TOL | Ipilimumab + Nivolumab | 2 | 72 | male |
| SRR32060215 | N015 | BC | control | AC | Pembrolizumab | 1 | 73 | male |
| SRR32060216 | N022 | BC | chili | TOL | Ipilimumab + Nivolumab | 1 | 67 | male |
| SRR32060234 | N159 | BC | control | AC | Ipilimumab + Nivolumab | 5 | 56 | female |
| SRR32060240 | N149 | BC | control | AC | Ipilimumab + Nivolumab | 5 | 76 | male |

## Coverage QC

QC threshold: `protein_coding_annotation_fraction >= 0.5`.

Excluded from the QC subset:

| sample | Group | raw_Group | protein_coding_annotation_fraction |
| --- | --- | --- | --- |
| SRR32060276 | TOL | chili | 0.08270626399219036 |

Samples with zero strict amyloid matched pairs for the primary paired-shift endpoint:

| sample | Group | raw_Group | protein_coding_annotation_fraction | strict_amyloid_TPM_fraction |
| --- | --- | --- | --- | --- |
| SRR32060309 | AC | control | 0.6567726053248654 | 0.0 |
| SRR32060240 | AC | control | 0.8274076227076321 | 0.0 |

## Main result

Primary matched endpoint: sample-level mean paired log2 expression shift, defined as:

`mean(log2(TPM_amyloid + 1) - log2(TPM_matched_non_amyloid + 1))`

All BC samples:

- n AC/control = 4; n TOL/chili = 5
- mean AC = 3.1011; mean TOL = 3.24597
- delta mean TOL - AC = 0.144874
- Cliff's delta TOL vs AC = 0.4
- exact one-sided permutation p (TOL > AC) = 0.174603
- exact two-sided permutation p = 0.373016
- bootstrap 95% CI for delta mean = [-0.0538387, 0.376494]

Coverage-QC subset:

- n AC/control = 4; n TOL/chili = 4
- mean AC = 3.1011; mean TOL = 3.143
- delta mean TOL - AC = 0.0418961
- Cliff's delta TOL vs AC = 0.25
- exact one-sided permutation p (TOL > AC) = 0.314286
- exact two-sided permutation p = 0.628571
- bootstrap 95% CI for delta mean = [-0.0947629, 0.174771]

## Interpretation

Under this within-sample matched-control design, the QC-subset direction supports the hypothesis that baseline ChILI/TOL samples have a higher amyloidogenic expression shift than control/AC samples. Because the sample size is small and the confidence interval is wide, this should be treated as discovery/exploratory evidence, not validation and not a proven biomarker.

## Limitations

- Only 11 BC samples from the metadata were present in `project.zarr`; the QC subset contains 10 samples after excluding `SRR32060276`.
- Matching used expression and protein length, with GENCODE `protein_coding` biotype enforced before matching. GC/content and subcellular localization were not locally available.
- Matching on observed expression can deliberately attenuate expression-burden differences; this method asks a conservative residual question rather than reproducing the unadjusted burden endpoint.
- Small n makes p-values unstable; effect size, direction, and leave-one-out stability are more informative here.

## Next best step

Validate the direction in an independent cohort or rerun this matched analysis after adding transcript GC and HPA/subcellular annotations, then use a counts-based model or sample-level mixed/permutation framework for inference.
