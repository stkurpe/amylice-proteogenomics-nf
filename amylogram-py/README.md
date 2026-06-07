# AmyloGram-Py

Fast, testable Python reimplementation of the AmyloGram prediction path.

The goal is exact compatibility with the original R `AmyloGram` package:

```text
protein sequence
  -> overlapping 6-mers
  -> AmyloGram amino-acid group encoding
  -> binary multigram features
  -> exported ranger random forest
  -> max 6-mer amyloid probability per protein
```

This repository is intentionally separate from the main mutant proteome
pipeline. The R implementation remains the reference oracle while this package
is developed and validated.

## Development Stages

1. Export AmyloGram reference fixtures from R.
2. Match amino-acid cleaning and group encoding.
3. Match 6-mer feature vectors against R/biogram fixtures.
4. Match ranger forest probabilities against R.
5. Precompute all `6^6 = 46656` degenerate 6-mer probabilities.
6. Predict full FASTA files with lookup + max aggregation.

## Quick Tests

```bash
python3 -m unittest discover -s tests -v
```

Install the CLI locally in editable mode:

```bash
make install-dev
.venv/bin/amylogram-py --help
```

The package requires Python `>=3.10`. If your system `python3` is older, pass
an explicit interpreter and venv path:

```bash
make install-dev PYTHON=/path/to/python3.12 VENV=.venv312
.venv312/bin/amylogram-py --help
```

Common development commands:

```bash
make test
make parity
make benchmark
make demo
```

R parity fixtures can be regenerated on a machine/container with the R
`AmyloGram`, `biogram`, `seqinr`, and `jsonlite` packages:

```bash
Rscript scripts/export_reference_fixtures.R tests/fixtures/amylogram_reference.json
```

The same export can be reproduced through Docker:

```bash
bash scripts/run_reference_export_docker.sh tests/fixtures/amylogram_reference.json
```

## Build The Lookup Table

After the R fixture is available, export all `46656` degenerate 6-mer
probabilities:

```bash
PYTHONPATH=src python3 -m amylogram_py.cli \
  build-table \
  tests/fixtures/amylogram_reference.json \
  tests/fixtures/amylogram_sixmer_probabilities.json
```

For production runs, prefer the binary table:

```bash
PYTHONPATH=src python3 -m amylogram_py.cli \
  build-table \
  tests/fixtures/amylogram_reference.json \
  tests/fixtures/amylogram_sixmer_probabilities.bin \
  --format bin
```

The lookup table is the production path: each protein is encoded into
overlapping 6-mers, each 6-mer is resolved by array lookup, and the final
protein probability is the maximum window probability.

The CLI uses a streaming rolling-code implementation, so it does not allocate
the full list of 6-mer windows for each protein. FASTA records are cleaned once
by the parser and then predicted through a clean-sequence fast path.

## Predict A FASTA

```bash
PYTHONPATH=src python3 -m amylogram_py.cli \
  input.fasta \
  amylogram_predictions.csv \
  --sixmer-table tests/fixtures/amylogram_sixmer_probabilities.bin \
  --report-md amylogram_prediction_report.md \
  --report-json amylogram_prediction_report.json \
  --skipped-tsv amylogram_skipped.tsv \
  --top-k 100 \
  --top-tsv amylogram_top_hits.tsv
```

Output columns:

- `Sequence_ID`
- `AmyloGram_Prob`
- `AmyloGram_Pred`

The Markdown and JSON reports record run-level metrics such as total FASTA
records, predicted records, skipped records, elapsed seconds, throughput,
threshold, and SHA256 hashes for the input FASTA and lookup table. The skipped
TSV lists records that were too short or empty after amino-acid cleaning. The
top-hit TSV keeps the highest-probability proteins without storing all
predictions in memory.

## Evidence And Benchmarks

Parity report:

```bash
python3 evidence/run_parity_checks.py --report evidence/amylogram_parity_report.md
```

Lookup benchmark:

```bash
python3 evidence/benchmark_lookup_prediction.py --count 2000 --length 300 --report evidence/lookup_benchmark.md
```

## Docker

Build the lightweight production image:

```bash
make docker-build
```

Run the container with mounted FASTA and lookup table:

```bash
docker run --rm \
  -v "$PWD/evidence:/data/evidence:rw" \
  -v "$PWD/tests/fixtures:/data/tests/fixtures:ro" \
  amylogram-py:local \
  evidence/demo_input.fasta \
  evidence/docker_smoke_predictions.csv \
  --sixmer-table tests/fixtures/amylogram_sixmer_probabilities.bin \
  --report-md evidence/docker_smoke_prediction_report.md \
  --report-json evidence/docker_smoke_prediction_report.json \
  --skipped-tsv evidence/docker_smoke_skipped.tsv \
  --top-k 2 \
  --top-tsv evidence/docker_smoke_top_hits.tsv
```

The same smoke test is available as:

```bash
make docker-smoke
```
