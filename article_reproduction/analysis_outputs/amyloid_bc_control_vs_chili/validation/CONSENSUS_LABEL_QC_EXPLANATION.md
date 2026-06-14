# Consensus label QC explanation

Validation Zarr:
`/Users/user841/Projects/BIoinfServer/Results_analysis/bioinfo_zarr/validation_project.zarr`

## Why there are only Partial and Discordant labels

The updated validation amyloid table has two populated predictors:

- `AMYPred_Pred`
- `AmyloGramPy_Pred`

The legacy/original `AmyloGram_Pred` and `AmyloGram_Prob` columns are present but empty for all validation protein rows.

Because the strict consensus requires the full predictor panel to agree, validation rows cannot become strict `Amyloid`. The observed pattern is:

- `Partial`: mostly `AMYPred_Pred = Amyloid` and `AmyloGramPy_Pred = Amyloid`, with missing `AmyloGram_Pred`.
- `Discordant`: mostly `AMYPred_Pred = Non-Amyloid` and `AmyloGramPy_Pred = Amyloid`, with missing `AmyloGram_Pred`.

This is a technical annotation mismatch relative to discovery, not evidence that the validation proteome lacks amyloidogenic proteins.

## Fractions at raw protein-row level

Across all validation protein rows:

- `Partial`: 130,454 / 195,721 = 66.65%
- `Discordant`: 65,267 / 195,721 = 33.35%

Per sample, `Partial` is approximately 66.3-67.2% and `Discordant` is approximately 32.8-33.7%.

## Fractions after GENCODE protein-coding and sample-transcript collapse

Across all validation collapsed protein-coding transcripts:

- `Partial`: 77,185 / 113,620 = 67.93%
- `Discordant`: 36,435 / 113,620 = 32.07%

Per sample, `Partial` is approximately 67.6-68.3% and `Discordant` is approximately 31.7-32.4%.

## Discovery contrast

In the discovery Zarr, samples with strict `Amyloid` have all three predictor columns populated. Example discovery rows include non-empty `AMYPred_Pred`, `AmyloGram_Pred`, and `AmyloGramPy_Pred`; strict `Amyloid` appears when the populated predictors agree.

Therefore the validation/discovery mismatch is most likely caused by incomplete validation predictor coverage, specifically missing legacy/original `AmyloGram` results.

## Practical implication

The original strict consensus endpoint is not directly testable on this validation Zarr until the missing `AmyloGram` predictor is generated or the consensus rule is redefined prospectively for a two-predictor panel.
