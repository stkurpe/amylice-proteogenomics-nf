#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$PROJECT_DIR"

if [[ -d /opt/homebrew/opt/openjdk/bin ]]; then
  export PATH="/opt/homebrew/opt/openjdk/bin:${PATH}"
fi
export NXF_HOME="${NXF_HOME:-${PROJECT_DIR}/.nextflow-home}"
mkdir -p "$NXF_HOME"

if ! command -v nextflow >/dev/null 2>&1; then
  echo "Missing required command: nextflow" >&2
  exit 127
fi

OUTDIR="${OUTDIR:-${PROJECT_DIR}/results_nextflow_amyloid_smoke}"
INPUT_FASTA="${PROJECT_DIR}/tests_nextflow/fixtures/amyloid_minimal/combine_proteome.fasta"
SIXMER_TABLE="${PROJECT_DIR}/amylogram-py/tests/fixtures/amylogram_sixmer_probabilities.bin"
RESULTS_DIR="${OUTDIR}/fixture_minimal/results_amyloid"

docker build -t amylogram-py-nextflow:local -f docker/amylogram-py-nextflow/Dockerfile amylogram-py

nextflow run nextflow/main.nf \
  -profile test,docker \
  --mode amyloid \
  --input_fasta "$INPUT_FASTA" \
  --sixmer_table "$SIXMER_TABLE" \
  --outdir "$OUTDIR"

required=(
  amylogram_py_prediction.csv
  amyloid_combined_predictions.csv
  amyloid_predictors_summary.tsv
  amyloid_predictors_status.tsv
)

for filename in "${required[@]}"; do
  test -s "${RESULTS_DIR}/${filename}"
done

python3 tests_nextflow/validate_amyloid_outputs.py "$RESULTS_DIR"

combined_rows="$(awk -F '\t' '$1 == "combined_rows" {print $2}' "${RESULTS_DIR}/amyloid_predictors_summary.tsv")"
test "${combined_rows:-0}" -gt 0
! grep -q $'\tFAIL\t' "${RESULTS_DIR}/amyloid_predictors_status.tsv"
