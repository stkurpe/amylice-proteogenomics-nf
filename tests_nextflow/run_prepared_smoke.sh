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

OUTDIR="${OUTDIR:-${PROJECT_DIR}/results_nextflow_smoke}"
SAMPLES="${PROJECT_DIR}/tests_nextflow/fixtures/prepared_minimal/samples.prepared.csv"
REF_GENOME="${PROJECT_DIR}/tests_nextflow/fixtures/prepared_minimal/genome.fa"
REF_GTF="${PROJECT_DIR}/tests_nextflow/fixtures/prepared_minimal/annotation.gtf"
RESULTS_DIR="${OUTDIR}/fixture_minimal/results_proteins"
AMYLOID_RESULTS_DIR="${OUTDIR}/fixture_minimal/results_amyloid"
FEATURE_RESULTS_DIR="${OUTDIR}/fixture_minimal/results_protein_features"
EXPECTED_DIR="${PROJECT_DIR}/tests_nextflow/fixtures/prepared_minimal/expected/results_proteins"

docker build -t amyloid-proteome-nextflow:local docker/proteome-nextflow
docker build -t amylogram-py-nextflow:local -f docker/amylogram-py-nextflow/Dockerfile amylogram-py

nextflow run nextflow/main.nf \
  -profile test,docker \
  --mode prepared \
  --samples "$SAMPLES" \
  --ref_genome "$REF_GENOME" \
  --ref_gtf "$REF_GTF" \
  --outdir "$OUTDIR"

required=(
  clean_ids.txt
  protein.fasta
  frameshift_unique.fasta
  combine_proteome.fasta
  nonsense_candidates.txt
  verification_report.tsv
  manifest.txt
)

for filename in "${required[@]}"; do
  test -s "${RESULTS_DIR}/${filename}"
done

grep -q OK "${RESULTS_DIR}/verification_report.tsv"
records="$(grep -c '^>' "${RESULTS_DIR}/combine_proteome.fasta")"
test "$records" -gt 0

python3 tests_nextflow/compare_proteome_outputs.py \
  --expected-output "$EXPECTED_DIR" \
  --nf-output "$RESULTS_DIR"

amyloid_required=(
  amylogram_py_prediction.csv
  amyloid_combined_predictions.csv
  amyloid_predictors_summary.tsv
  amyloid_predictors_status.tsv
)

for filename in "${amyloid_required[@]}"; do
  test -s "${AMYLOID_RESULTS_DIR}/${filename}"
done

python3 tests_nextflow/validate_amyloid_outputs.py "$AMYLOID_RESULTS_DIR"

feature_required=(
  protein_features_light.csv
  protein_features_light_summary.tsv
  protein_features_light_status.tsv
)

for filename in "${feature_required[@]}"; do
  test -s "${FEATURE_RESULTS_DIR}/${filename}"
done

grep -q protein_id "${FEATURE_RESULTS_DIR}/protein_features_light.csv"
grep -q OK "${FEATURE_RESULTS_DIR}/protein_features_light_status.tsv"
feature_rows="$(tail -n +2 "${FEATURE_RESULTS_DIR}/protein_features_light.csv" | grep -c .)"
test "$feature_rows" -gt 0
