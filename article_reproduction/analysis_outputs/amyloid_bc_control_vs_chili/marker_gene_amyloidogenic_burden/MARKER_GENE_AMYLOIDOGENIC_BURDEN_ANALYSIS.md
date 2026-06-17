# Marker-gene amyloidogenic burden

This analysis recalculates amyloidogenic burden within ALB, AFP, CRP, IL6, and SLFN11 for baseline AC/control and TOL/ChILI samples. Metrics are gene-local: the normalized and KDmax-weighted fractions use each gene's classified transcripts as the denominator, not the whole protein-coding universe.

## Outputs

- `/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/marker_gene_amyloidogenic_burden/marker_gene_absolute-amyloidogenic-burden_boxplots.png`
- `/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/marker_gene_amyloidogenic_burden/marker_gene_absolute-amyloidogenic-burden_boxplots.pdf`
- `/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/marker_gene_amyloidogenic_burden/marker_gene_normalized-amyloidogenic-fraction_boxplots.png`
- `/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/marker_gene_amyloidogenic_burden/marker_gene_normalized-amyloidogenic-fraction_boxplots.pdf`
- `/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/marker_gene_amyloidogenic_burden/marker_gene_KDmax-weighted-amyloidogenic-fraction_boxplots.png`
- `/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/marker_gene_amyloidogenic_burden/marker_gene_KDmax-weighted-amyloidogenic-fraction_boxplots.pdf`
- `/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/marker_gene_amyloidogenic_burden/marker_gene_amyloidogenic_burden_sample_values.tsv`
- `/Users/user841/Documents/Article/analysis_outputs/amyloid_bc_control_vs_chili/marker_gene_amyloidogenic_burden/marker_gene_amyloidogenic_burden_statistics.tsv`

## Statistical summary

### ALB
- Absolute burden: mean AC=0.328, mean TOL=2.715, delta TOL-AC=2.387, Cliff's delta=0.417, exact Mann-Whitney p=0.352.
- Normalized fraction: mean AC=0.5, mean TOL=1.000, delta TOL-AC=0.5, Cliff's delta=0.5, exact Mann-Whitney p=0.667.
- KDmax-weighted fraction: mean AC=0.5, mean TOL=1.000, delta TOL-AC=0.5, Cliff's delta=0.5, exact Mann-Whitney p=0.667.

### AFP
- Absolute burden: mean AC=0.531, mean TOL=0, delta TOL-AC=-0.531, Cliff's delta=-0.167, exact Mann-Whitney p=0.762.
- Normalized fraction: mean AC=1.000, mean TOL=NA, delta TOL-AC=NA, Cliff's delta=NA, exact Mann-Whitney p=NA.
- KDmax-weighted fraction: mean AC=1.000, mean TOL=NA, delta TOL-AC=NA, Cliff's delta=NA, exact Mann-Whitney p=NA.

### CRP
- Absolute burden: mean AC=0, mean TOL=0, delta TOL-AC=0, Cliff's delta=0, exact Mann-Whitney p=NA.
- Normalized fraction: mean AC=NA, mean TOL=NA, delta TOL-AC=NA, Cliff's delta=NA, exact Mann-Whitney p=NA.
- KDmax-weighted fraction: mean AC=NA, mean TOL=NA, delta TOL-AC=NA, Cliff's delta=NA, exact Mann-Whitney p=NA.

### IL6
- Absolute burden: mean AC=0, mean TOL=0, delta TOL-AC=0, Cliff's delta=0, exact Mann-Whitney p=NA.
- Normalized fraction: mean AC=NA, mean TOL=NA, delta TOL-AC=NA, Cliff's delta=NA, exact Mann-Whitney p=NA.
- KDmax-weighted fraction: mean AC=NA, mean TOL=NA, delta TOL-AC=NA, Cliff's delta=NA, exact Mann-Whitney p=NA.

### SLFN11
- Absolute burden: mean AC=0, mean TOL=0, delta TOL-AC=0, Cliff's delta=0, exact Mann-Whitney p=NA.
- Normalized fraction: mean AC=0, mean TOL=0, delta TOL-AC=0, Cliff's delta=0, exact Mann-Whitney p=NA.
- KDmax-weighted fraction: mean AC=0, mean TOL=0, delta TOL-AC=0, Cliff's delta=0, exact Mann-Whitney p=NA.

## Interpretation

CRP and IL6 have no measurable baseline RNA-seq amyloidogenic burden in these samples, so their gene-local metrics are not informative here. ALB shows higher absolute amyloidogenic burden in TOL/ChILI, driven by the few samples where ALB has classified amyloidogenic expression. ALB within-gene fractions are high among samples with nonzero denominators, but the denominator is available in only two AC and two TOL samples. AFP has a single AC sample with amyloidogenic burden and no TOL signal in this matched baseline set. SLFN11 is expressed in several samples, but its classified transcripts are non-amyloidogenic in this q17/q83 consensus framework, so SLFN11 gene-local amyloidogenic burden is zero. This differs from the previous analysis where SLFN11 expression was evaluated as a marker associated with global protein-coding amyloidogenic burden.
