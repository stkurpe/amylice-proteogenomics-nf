# Extreme Contributor Concentration Index

## What was done
Recomputed from `project.zarr` using only baseline (`time == "BC"`) samples. Metadata `Group=control` was mapped to AC and `Group=chili` to TOL. The transcript universe was restricted to local GENCODE `biotype == "protein_coding"` transcripts. Protein variants were collapsed to `sample x transcript` before joining transcript-level TPM, avoiding H1/H2 double counting.

Strict amyloid labels were defined as `Consensus_collapsed == "Amyloid"`. The primary concentration metrics were the Gini coefficient of strict amyloid TPM contribution and the top-10 contributor share within each sample. I also computed secondary top-5/top-1 and positive weighted-contribution concentration metrics.

## Samples and QC
- SRR32060276: TOL/chili, subject N132, protein-coding annotation fraction=0.098, excluded from coverage-QC subset
- SRR32060306: TOL/chili, subject N101, protein-coding annotation fraction=0.929, kept in coverage-QC subset
- SRR32060309: AC/control, subject N87, protein-coding annotation fraction=0.759, kept in coverage-QC subset
- SRR32060312: AC/control, subject N106, protein-coding annotation fraction=0.901, kept in coverage-QC subset
- SRR32060326: AC/control, subject N032, protein-coding annotation fraction=0.951, kept in coverage-QC subset
- SRR32060336: TOL/chili, subject N072, protein-coding annotation fraction=0.938, kept in coverage-QC subset
- SRR32060341: TOL/chili, subject N064, protein-coding annotation fraction=0.750, kept in coverage-QC subset
- SRR32060215: AC/control, subject N015, protein-coding annotation fraction=0.945, kept in coverage-QC subset
- SRR32060216: TOL/chili, subject N022, protein-coding annotation fraction=0.885, kept in coverage-QC subset
- SRR32060234: AC/control, subject N159, protein-coding annotation fraction=0.978, kept in coverage-QC subset
- SRR32060240: AC/control, subject N149, protein-coding annotation fraction=0.964, kept in coverage-QC subset

## Main results
- all_BC_samples / strict_amyloid_TPM_contribution_gini: mean AC=0.8833, mean TOL=0.8929, delta TOL-AC=0.0096 (bootstrap 95% CI -0.0337 to 0.0473); Cliff's delta=0.300; MW two-sided p=0.556; exact permutation one-sided p(TOL>AC)=0.354.
- all_BC_samples / strict_amyloid_TPM_top10_share: mean AC=0.6491, mean TOL=0.7030, delta TOL-AC=0.0539 (bootstrap 95% CI -0.0878 to 0.1935); Cliff's delta=0.300; MW two-sided p=0.556; exact permutation one-sided p(TOL>AC)=0.283.
- coverage_QC_protein_coding_annotation_fraction_ge_0_5 / strict_amyloid_TPM_contribution_gini: mean AC=0.8833, mean TOL=0.8829, delta TOL-AC=-0.0004 (bootstrap 95% CI -0.0464 to 0.0384); Cliff's delta=0.125; MW two-sided p=0.886; exact permutation one-sided p(TOL>AC)=0.521.
- coverage_QC_protein_coding_annotation_fraction_ge_0_5 / strict_amyloid_TPM_top10_share: mean AC=0.6491, mean TOL=0.6502, delta TOL-AC=0.0011 (bootstrap 95% CI -0.1275 to 0.1120); Cliff's delta=0.125; MW two-sided p=0.886; exact permutation one-sided p(TOL>AC)=0.493.

## Interpretation
In the coverage-QC subset, the primary concentration metrics are not consistently higher in TOL/chili. This does not strongly support the concentration/nucleation-style hypothesis.

## Limitations
Small n, treatment/batch imbalance, and amyloid annotation coverage limit inference. Concentration metrics can be sensitive to one highly expressed transcript, so I report both all-sample and coverage-QC subsets and provide top contributor panels for inspection.

## Next best step
Validate the concentration signal in an independent baseline cohort, or test whether the top strict amyloid contributors are reproducible under counts-based models with treatment/batch covariates.
