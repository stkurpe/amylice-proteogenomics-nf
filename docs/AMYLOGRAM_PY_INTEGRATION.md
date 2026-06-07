# AmyloGram-Py Pipeline Integration

AmyloGram-Py is integrated as an optional third amyloidogenicity predictor in
`run_amyloid_predictors.sh`.

It is disabled by default so the existing long-running batch keeps its current
behavior. Enable it explicitly with:

```bash
RUN_AMYLOGRAM_PY=1 /home/codex/run_amyloid_predictors.sh --run-amylogram-py ...
```

## Required Server Artifacts

- Docker image: `amylogram-py:local`
- Lookup table:
  `/home/codex/amylogram-py/tests/fixtures/amylogram_sixmer_probabilities.bin`
- Updated pipeline script:
  `/home/codex/run_amyloid_predictors.sh`

## One-Sample Test Command

Run only the fast Python AmyloGram component against an already prepared
proteome:

```bash
RUN_AMYPRED=0 \
RUN_AMYLOGRAM=0 \
RUN_AMYLOGRAM_PY=1 \
AMYLOGRAM_PY_IMAGE=amylogram-py:local \
AMYLOGRAM_PY_LOOKUP=/home/codex/amylogram-py/tests/fixtures/amylogram_sixmer_probabilities.bin \
AMYLOGRAM_PY_TOP_K=1000 \
/home/codex/run_amyloid_predictors.sh \
  --sample SRR32060234 \
  --input-s3 s3://codex-test-ngsdata-calculations/prepared-bioinfo-proteome/SRR32060234/results_proteins/combine_proteome.fasta \
  --output-s3 s3://codex-test-ngsdata-calculations/prepared-bioinfo-proteome/SRR32060234/results_amyloid_py_test \
  --run-amylogram-py \
  --force
```

## Expected Local Outputs

Inside `/home/codex/amyloid_runs/SRR32060234/results/`:

- `amylogram_py_prediction.csv`
- `amylogram_py_report.md`
- `amylogram_py_report.json`
- `amylogram_py_skipped.tsv`
- `amylogram_py_top_hits.tsv`
- `amyloid_combined_predictions.csv`
- `amyloid_predictors_summary.tsv`
- `AMYLOID_PREDICTORS_SUCCESS`

## Expected S3 Outputs

Under the selected output prefix:

- `amylogram_py_prediction.csv`
- `amylogram_py_report.md`
- `amylogram_py_report.json`
- `amylogram_py_skipped.tsv`
- `amylogram_py_top_hits.tsv`
- `amyloid_combined_predictions.csv`
- `amyloid_predictors_status.tsv`
- `amyloid_predictors_summary.tsv`
- `logs/amyloid_predictors_<RUN_ID>.log`

## Safety Gate

Run the test only after the current amyloid/HLA batch is finished. Do not start
the test while any of these are active:

- `run_amyloid_after_proteomes_then_hla.sh`
- `run_amyloid_predictors.sh`
- `05_hla_consensus.sh`
- active HLA caller containers

If the batch is still active, only report status and wait for the next check.
