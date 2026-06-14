# Competitive Gene-Set Test: strict amyloid transcripts in baseline ChILI

## What was done
Recomputed from `project.zarr` using only baseline (`time == "BC"`) blood RNA-seq samples. Metadata `Group=control` was mapped to AC and `Group=chili` to TOL. Strict amyloid transcripts were defined by collapsed `sample x transcript` amyloid annotations with `Consensus_collapsed == "Amyloid"`; protein variants mapping to the same transcript were collapsed to avoid H1/H2 double counting.

Transcripts were ranked by counts-derived logCPM `log2FC_TOL_vs_AC`. The competitive test asks whether strict amyloid transcripts have systematically higher log2FC than the transcript background. I report Mann-Whitney/AUC, random-set empirical p-value, and phenotype-permutation GSEA.

## Samples
- SRR32060276: TOL, subject N132, coverage=0.087, excluded from coverage-QC
- SRR32060306: TOL, subject N101, coverage=0.938, kept
- SRR32060309: AC, subject N87, coverage=0.769, kept
- SRR32060312: AC, subject N106, coverage=0.917, kept
- SRR32060326: AC, subject N032, coverage=0.957, kept
- SRR32060336: TOL, subject N072, coverage=0.946, kept
- SRR32060341: TOL, subject N064, coverage=0.757, kept
- SRR32060215: AC, subject N015, coverage=0.955, kept
- SRR32060216: TOL, subject N022, coverage=0.892, kept
- SRR32060234: AC, subject N159, coverage=0.981, kept
- SRR32060240: AC, subject N149, coverage=0.969, kept

## Results
- all_BC_samples: n=11 (AC=6, TOL=5); strict amyloid tested=2095/21045. Mean log2FC amyloid=0.2291 (95% bootstrap CI 0.1958 to 0.2626); background=0.3038; delta=-0.0747; AUC=0.477; MW one-sided p=1; random-set empirical p=1; GSEA ES=0.2084, NES=0.971, empirical p=0.267.
- coverage_QC_protein_coding_annotation_fraction_ge_0_5: n=10 (AC=6, TOL=4); strict amyloid tested=1982/20110. Mean log2FC amyloid=0.2525 (95% bootstrap CI 0.2152 to 0.2896); background=0.3222; delta=-0.0698; AUC=0.482; MW one-sided p=0.996; random-set empirical p=1; GSEA ES=0.2139, NES=0.968, empirical p=0.289.

## Interpretation
The coverage-QC competitive gene-set result does not provide strong support that strict amyloid transcripts are systematically higher among genes upregulated in ChILI/TOL.

## Limitations
Small sample size, heterogeneous treatment/batch composition, and transcript-level amyloid annotation coverage limit inference. The permutation p-values are empirical and constrained by the number of possible label rearrangements in these small groups.

## Next best step
Run a counts-based model with design covariates where feasible, then validate the amyloid gene-set direction in an independent baseline cohort or by targeted qPCR/protein-level follow-up for the top contributing strict amyloid transcripts.
