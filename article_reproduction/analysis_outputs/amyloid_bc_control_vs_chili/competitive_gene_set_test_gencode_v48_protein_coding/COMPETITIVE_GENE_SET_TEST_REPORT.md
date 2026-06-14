# Competitive Gene-Set Test: strict amyloid transcripts in baseline ChILI

## What was done
Recomputed from `project.zarr` using only baseline (`time == "BC"`) blood RNA-seq samples. Metadata `Group=control` was mapped to AC and `Group=chili` to TOL. Transcript universe was restricted using local GENCODE annotation `/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/rank_based_amyloid_burden_signature/gencode_v48_transcript_annotation.tsv` to `biotype == "protein_coding"`. Strict amyloid transcripts were defined by collapsed `sample x transcript` amyloid annotations with `Consensus_collapsed == "Amyloid"`; protein variants mapping to the same transcript were collapsed to avoid H1/H2 double counting.

Transcripts were ranked by counts-derived logCPM `log2FC_TOL_vs_AC`. The competitive test asks whether strict amyloid transcripts have systematically higher log2FC than the transcript background. I report Mann-Whitney/AUC, random-set empirical p-value, and phenotype-permutation GSEA.

## Samples
- SRR32060276: TOL, subject N132, coverage=0.098, excluded from coverage-QC
- SRR32060306: TOL, subject N101, coverage=0.929, kept
- SRR32060309: AC, subject N87, coverage=0.759, kept
- SRR32060312: AC, subject N106, coverage=0.901, kept
- SRR32060326: AC, subject N032, coverage=0.951, kept
- SRR32060336: TOL, subject N072, coverage=0.938, kept
- SRR32060341: TOL, subject N064, coverage=0.750, kept
- SRR32060215: AC, subject N015, coverage=0.945, kept
- SRR32060216: TOL, subject N022, coverage=0.885, kept
- SRR32060234: AC, subject N159, coverage=0.978, kept
- SRR32060240: AC, subject N149, coverage=0.964, kept

## Results
- all_BC_samples: n=11 (AC=6, TOL=5); strict amyloid tested=1245/17776. Mean log2FC amyloid=0.1889 (95% bootstrap CI 0.1462 to 0.2316); background=0.2940; delta=-0.1051; AUC=0.467; MW one-sided p=1; random-set empirical p=1; GSEA ES=0.1842, NES=0.844, empirical p=0.408.
- coverage_QC_protein_coding_annotation_fraction_ge_0_5: n=10 (AC=6, TOL=4); strict amyloid tested=1182/17014. Mean log2FC amyloid=0.1970 (95% bootstrap CI 0.1492 to 0.2453); background=0.3078; delta=-0.1108; AUC=0.469; MW one-sided p=1; random-set empirical p=1; GSEA ES=0.1839, NES=0.790, empirical p=0.466.

## Interpretation
The coverage-QC competitive gene-set result does not provide strong support that strict amyloid transcripts are systematically higher among genes upregulated in ChILI/TOL.

## Limitations
Small sample size, heterogeneous treatment/batch composition, and transcript-level amyloid annotation coverage limit inference. The permutation p-values are empirical and constrained by the number of possible label rearrangements in these small groups.

## Next best step
Run a counts-based model with design covariates where feasible, then validate the amyloid gene-set direction in an independent baseline cohort or by targeted qPCR/protein-level follow-up for the top contributing strict amyloid transcripts.
