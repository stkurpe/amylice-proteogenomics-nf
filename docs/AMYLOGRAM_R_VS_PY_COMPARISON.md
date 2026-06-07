# R AmyloGram vs AmyloGram-Py Comparison

Generated UTC: `2026-05-18T18:42:23Z`

## Executive Summary

Across the compared samples, R AmyloGram and AmyloGram-Py produced the same
number of valid predictions, the same skipped counts, and the same binary
amyloid/non-amyloid counts. No prediction-label discordance was observed in
the shared unique sequence IDs. Probability differences are tiny and consistent
with numeric/export precision differences between the R path and the Python
lookup-table path.

## Per-Sample Counts

| Sample | Input records | R predicted | Py predicted | R skipped | Py skipped | R amyloid | Py amyloid | R non-amyloid | Py non-amyloid | Discordant labels |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SRR32060215 | 27980 | 27958 | 27958 | 22 | 22 | 26772 | 26772 | 1186 | 1186 | 0 |
| SRR32060233 | 25196 | 25178 | 25178 | 18 | 18 | 24359 | 24359 | 819 | 819 | 0 |
| SRR32060234 | 28327 | 28312 | 28312 | 15 | 15 | 27248 | 27248 | 1064 | 1064 | 0 |

## Runtime And Speed

| Sample | R runtime | Py runtime | Speedup | Py records/sec |
|---|---:|---:|---:|---:|
| SRR32060215 | 67.62 min | 3.79 sec | 1069.5x | 7370.1 |
| SRR32060233 | 51.11 min | 4.02 sec | 763.1x | 6265.7 |
| SRR32060234 | 51.77 min | 3.37 sec | 922.9x | 8411.6 |

## Prediction Agreement

| Sample | Shared unique IDs | R-only unique IDs | Py-only unique IDs | R duplicate ID rows | Py duplicate ID rows | Discordant labels | Mean abs prob delta | Median abs prob delta | Max abs prob delta | Max delta sequence | Near-threshold records |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| SRR32060215 | 27956 | 0 | 0 | 2 | 2 | 0 | 0.00000059 | 0.00000048 | 0.00000506 | `H1_ENST00000511392.5` | 189 |
| SRR32060233 | 25178 | 0 | 0 | 0 | 0 | 0 | 0.00000059 | 0.00000048 | 0.00000411 | `H1_ENST00000233331.12` | 155 |
| SRR32060234 | 28312 | 0 | 0 | 0 | 0 | 0 | 0.00000059 | 0.00000047 | 0.00000411 | `H1_ENST00000233331.12` | 166 |

## Where The Algorithms Stumble

Both implementations stumble at the same biological/input boundary: records
that are empty after amino-acid cleaning or shorter than six amino acids.
R reports these as `invalid_or_short_records`; AmyloGram-Py emits a
machine-readable skipped TSV with per-record reasons. R sequence-level
prediction errors were zero for all compared samples.

| Sample | R invalid/short | Py skipped | R sequence errors | Py skipped reasons |
|---|---:|---:|---:|---|
| SRR32060215 | 22 | 22 | 0 | shorter_than_6: 22 |
| SRR32060233 | 18 | 18 | 0 | shorter_than_6: 18 |
| SRR32060234 | 15 | 15 | 0 | shorter_than_6: 15 |

## Implementation Notes

- R AmyloGram used chunked/resumable execution and deduplication; its summary
  reports `unique_sequences` and `duplicates_saved`.
- AmyloGram-Py used a precomputed `6^6 = 46656` lookup table and streaming
  rolling 6-mer encoding, so it avoids the per-window R forest traversal.
- The compared Py outputs for `SRR32060234` are intentionally stored in a
  separate smoke-test prefix: `results_amyloid_py_test/`.

## R Execution Details

| Sample | R chunk size | R chunk count | Unique sequences | Duplicates saved |
|---|---:|---:|---:|---:|
| SRR32060215 | 200 | 88 | 17496 | 10462 |
| SRR32060233 | 200 | 62 | 12314 | 12864 |
| SRR32060234 | 50 | 284 | 14199 | 14113 |

## Py Top Signal

| Sample | Py mean probability | Py max probability | Py max sequence ID |
|---|---:|---:|---|
| SRR32060215 | 0.786940 | 0.955704 | `H1_ENST00000376629.8` |
| SRR32060233 | 0.795225 | 0.955704 | `H1_ENST00000376629.8` |
| SRR32060234 | 0.790740 | 0.955704 | `H1_ENST00000376629.8` |
