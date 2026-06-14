# Image MixedLM Feature Panel: Amyloidogenicity

## What was done
Used the gene/protein IDs visible in the supplied MixedLM feature plot. The analysis was restricted to BC samples, AC/control vs TOL/chili, and to GENCODE protein-coding transcripts. H1/H2 protein IDs from the image were retained at protein level; transcript IDs are included so H1/H2 duplicates can be inspected.

Amyloidogenicity was summarized with the existing integrative index: `z(beta_propensity) + z(KD_max) + z(aromaticity) - z(abs_charge_density) - z(pI_distance_from_7)`, divided by 5. I also report AmyloGram probability, AMYPred probability, strict consensus Amyloid fraction, and expression-weighted panel index.

## Main results
- all_BC_samples / mean_Amyloid_Integrative_Index: mean AC=-0.3171, mean TOL=-0.3263, delta TOL-AC=-0.0092, Cliff's delta=-0.200, MW p=0.662.
- all_BC_samples / mean_AmyloGram_Prob: mean AC=0.6980, mean TOL=0.7085, delta TOL-AC=0.0105, Cliff's delta=0.333, MW p=0.429.
- all_BC_samples / mean_AMYPred_Prob: mean AC=0.5631, mean TOL=0.5590, delta TOL-AC=-0.0041, Cliff's delta=0.000, MW p=1.
- all_BC_samples / strict_amyloid_fraction: mean AC=0.5031, mean TOL=0.7386, delta TOL-AC=0.2355, Cliff's delta=0.300, MW p=0.462.
- all_BC_samples / expression_weighted_amyloid_index_per_panel_TPM: mean AC=0.4645, mean TOL=0.3392, delta TOL-AC=-0.1254, Cliff's delta=-0.267, MW p=0.537.
- coverage_QC_protein_coding_annotation_fraction_ge_0_5 / mean_Amyloid_Integrative_Index: mean AC=-0.3171, mean TOL=-0.3708, delta TOL-AC=-0.0537, Cliff's delta=-0.500, MW p=0.257.
- coverage_QC_protein_coding_annotation_fraction_ge_0_5 / mean_AmyloGram_Prob: mean AC=0.6980, mean TOL=0.6997, delta TOL-AC=0.0017, Cliff's delta=0.167, MW p=0.762.
- coverage_QC_protein_coding_annotation_fraction_ge_0_5 / mean_AMYPred_Prob: mean AC=0.5631, mean TOL=0.5720, delta TOL-AC=0.0089, Cliff's delta=0.250, MW p=0.61.
- coverage_QC_protein_coding_annotation_fraction_ge_0_5 / strict_amyloid_fraction: mean AC=0.5031, mean TOL=0.7482, delta TOL-AC=0.2451, Cliff's delta=0.458, MW p=0.282.
- coverage_QC_protein_coding_annotation_fraction_ge_0_5 / expression_weighted_amyloid_index_per_panel_TPM: mean AC=0.4645, mean TOL=0.2236, delta TOL-AC=-0.2410, Cliff's delta=-0.583, MW p=0.171.

## Top protein-level shifts by integrative index
- ZNF580 | H2_ENST00000545125.1: delta index=0.384, delta AmyloGram=0.137, delta AMYPred=0.084.
- STMP1 | H2_ENST00000507606.3: delta index=0.317, delta AmyloGram=0.077, delta AMYPred=-0.175.
- FXYD5 | H2_ENST00000392217.3: delta index=0.179, delta AmyloGram=-0.005, delta AMYPred=0.028.
- SEPTIN7 | H2_ENST00000485569.1: delta index=0.097, delta AmyloGram=0.000, delta AMYPred=-0.001.
- AP2B1 | H2_ENST00000590432.5: delta index=0.073, delta AmyloGram=0.000, delta AMYPred=0.035.
- YLPM1 | H2_ENST00000554107.2: delta index=0.065, delta AmyloGram=0.000, delta AMYPred=-0.009.
- SNX10 | H2_ENST00000409838.1: delta index=0.061, delta AmyloGram=0.000, delta AMYPred=-0.056.
- TACC1 | H2_ENST00000520611.1: delta index=0.053, delta AmyloGram=0.000, delta AMYPred=-0.000.
- HBP1 | H2_ENST00000479011.1: delta index=0.047, delta AmyloGram=0.000, delta AMYPred=-0.033.
- PSMD8 | H2_ENST00000591250.1: delta index=0.037, delta AmyloGram=-0.000, delta AMYPred=-0.015.
- PPP2R5C | H2_ENST00000555237.5: delta index=0.036, delta AmyloGram=0.000, delta AMYPred=0.052.
- MRPS21 | H2_ENST00000614145.5: delta index=0.031, delta AmyloGram=0.000, delta AMYPred=-0.005.

## Interpretation
This is an exploratory, image-selected panel rather than independent validation. Direction should be interpreted mainly as whether the proteins highlighted by the feature plot become more amyloidogenic in TOL/chili after combining the amyloid-related physicochemical features and expression.
