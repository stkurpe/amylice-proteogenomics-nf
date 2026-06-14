# Updated validation comparison

Validation Zarr:
`/Users/user841/Projects/BIoinfServer/Results_analysis/bioinfo_zarr/validation_project.zarr`

GENCODE filter:
Only `biotype == protein_coding` transcripts were retained.

Coverage QC:
All 8 samples passed `protein_coding_annotation_fraction >= 0.5`.

## Strict consensus: `Consensus_collapsed == "Amyloid"`

Result: not confirmed.

There were no strict `Amyloid` transcript labels after sample x transcript collapse.

Primary endpoint:
`strict_amyloid_TPM_fraction`

- AC mean: 0
- TOL mean: 0
- delta TOL minus AC: 0
- one-sided Mann-Whitney p: 1.0
- Cliff's delta: 0
- 10,000 label permutations empirical p: 1.0

## Partial consensus sensitivity: `Consensus_collapsed == "Partial"`

Result: not confirmed.

Primary endpoint:
`partial_consensus_TPM_fraction`

- AC mean: 0.744164
- TOL mean: 0.712364
- delta TOL minus AC: -0.0318001
- one-sided Mann-Whitney p for TOL > AC: 0.828571
- Cliff's delta: -0.375
- bootstrap 95% CI: [-0.113369, 0.0504841]
- 10,000 label permutations empirical p: 0.736326

## Label QC

The updated validation file now has `Partial` and `Discordant` consensus labels, but still no strict `Amyloid` labels. This means the original strict amyloid burden finding cannot be confirmed in this validation cohort. The Partial sensitivity also points opposite to the original TOL > AC hypothesis.
