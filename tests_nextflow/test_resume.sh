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

RUN_ROOT="$(mktemp -d "${TMPDIR:-/private/tmp}/nf-prepared-resume.XXXXXX")"
OUTDIR="${RUN_ROOT}/out"
WORKDIR="${RUN_ROOT}/work"
SAMPLES="${PROJECT_DIR}/tests_nextflow/fixtures/prepared_minimal/samples.prepared.csv"
REF_GENOME="${PROJECT_DIR}/tests_nextflow/fixtures/prepared_minimal/genome.fa"
REF_GTF="${PROJECT_DIR}/tests_nextflow/fixtures/prepared_minimal/annotation.gtf"
RESULTS_DIR="${OUTDIR}/fixture_minimal/results_proteins"
AMYLOID_RESULTS_DIR="${OUTDIR}/fixture_minimal/results_amyloid"
FEATURE_RESULTS_DIR="${OUTDIR}/fixture_minimal/results_protein_features"
EXPECTED_DIR="${PROJECT_DIR}/tests_nextflow/fixtures/prepared_minimal/expected/results_proteins"
TRACE_FIRST="${RUN_ROOT}/trace.first.txt"
TRACE_SECOND="${RUN_ROOT}/trace.second.txt"
REPORT_FIRST="${RUN_ROOT}/report.first.html"
REPORT_SECOND="${RUN_ROOT}/report.second.html"
TIMELINE_FIRST="${RUN_ROOT}/timeline.first.html"
TIMELINE_SECOND="${RUN_ROOT}/timeline.second.html"

docker build -t amyloid-proteome-nextflow:local docker/proteome-nextflow
docker build -t amylogram-py-nextflow:local -f docker/amylogram-py-nextflow/Dockerfile amylogram-py

nextflow run nextflow/main.nf \
  -profile test,docker \
  -work-dir "$WORKDIR" \
  -with-trace "$TRACE_FIRST" \
  -with-report "$REPORT_FIRST" \
  -with-timeline "$TIMELINE_FIRST" \
  --mode prepared \
  --samples "$SAMPLES" \
  --ref_genome "$REF_GENOME" \
  --ref_gtf "$REF_GTF" \
  --outdir "$OUTDIR"

nextflow run nextflow/main.nf \
  -profile test,docker \
  -resume \
  -work-dir "$WORKDIR" \
  -with-trace "$TRACE_SECOND" \
  -with-report "$REPORT_SECOND" \
  -with-timeline "$TIMELINE_SECOND" \
  --mode prepared \
  --samples "$SAMPLES" \
  --ref_genome "$REF_GENOME" \
  --ref_gtf "$REF_GTF" \
  --outdir "$OUTDIR"

for report_file in "$TRACE_FIRST" "$TRACE_SECOND" "$REPORT_FIRST" "$REPORT_SECOND" "$TIMELINE_FIRST" "$TIMELINE_SECOND"; do
  test -s "$report_file"
done

required_proteome=(
  clean_ids.txt
  protein.fasta
  frameshift_unique.fasta
  combine_proteome.fasta
  nonsense_candidates.txt
  verification_report.tsv
  manifest.txt
)

for filename in "${required_proteome[@]}"; do
  test -s "${RESULTS_DIR}/${filename}"
done

required_amyloid=(
  amylogram_py_prediction.csv
  amyloid_combined_predictions.csv
  amyloid_predictors_summary.tsv
  amyloid_predictors_status.tsv
)

for filename in "${required_amyloid[@]}"; do
  test -s "${AMYLOID_RESULTS_DIR}/${filename}"
done

required_features=(
  protein_features_light.csv
  protein_features_light_summary.tsv
  protein_features_light_status.tsv
)

for filename in "${required_features[@]}"; do
  test -s "${FEATURE_RESULTS_DIR}/${filename}"
done

grep -q OK "${RESULTS_DIR}/verification_report.tsv"
grep -q CACHED "$TRACE_SECOND"
grep -Eq 'CACHED|COMPLETED' "$TRACE_SECOND"

python3 tests_nextflow/compare_proteome_outputs.py \
  --expected-output "$EXPECTED_DIR" \
  --nf-output "$RESULTS_DIR"
python3 tests_nextflow/validate_amyloid_outputs.py "$AMYLOID_RESULTS_DIR"
grep -q OK "${FEATURE_RESULTS_DIR}/protein_features_light_status.tsv"

echo "Prepared resume test passed: ${RUN_ROOT}"
