#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$PROJECT_DIR"

usage() {
  cat >&2 <<'USAGE'
Usage:
  AWS_PREPARED_ABUNDANCE=/path/abundance.tsv \
  AWS_PREPARED_VCF=/path/variants_filtered.vcf \
  AWS_REFERENCE_FASTA=/path/reference/combine_proteome.fasta \
  REF_GENOME=/path/GRCh38.primary_assembly.genome.fa \
  REF_GTF=/path/gencode.v46.primary_assembly.annotation.gtf \
  bash tests_nextflow/run_aws_prepared_reference.sh

Optional:
  SAMPLE_ID=SRR32060234
  OUTDIR=/path/results_nextflow_aws_reference
  EXPECTED_RECORDS=28343
  EXPECTED_SEQUENCE_SET_MD5=de4d381238adb9c1e72e714e1e72abb1
  RUN_AMYPRED=1
  RUN_AMYLOGRAM_R=1
  RUN_PROTEIN_FEATURES=1
USAGE
}

require_file() {
  local name="$1"
  local path="$2"
  if [[ -z "$path" ]]; then
    echo "Missing required environment variable: ${name}" >&2
    usage
    exit 2
  fi
  if [[ ! -s "$path" ]]; then
    echo "${name} does not point to a non-empty file: ${path}" >&2
    exit 2
  fi
}

if [[ -d /opt/homebrew/opt/openjdk/bin ]]; then
  export PATH="/opt/homebrew/opt/openjdk/bin:${PATH}"
fi
export NXF_HOME="${NXF_HOME:-${PROJECT_DIR}/.nextflow-home}"
mkdir -p "$NXF_HOME"

if ! command -v nextflow >/dev/null 2>&1; then
  echo "Missing required command: nextflow" >&2
  exit 127
fi

SAMPLE_ID="${SAMPLE_ID:-SRR32060234}"
OUTDIR="${OUTDIR:-${PROJECT_DIR}/results_nextflow_aws_reference}"
EXPECTED_RECORDS="${EXPECTED_RECORDS:-28343}"
EXPECTED_SEQUENCE_SET_MD5="${EXPECTED_SEQUENCE_SET_MD5:-de4d381238adb9c1e72e714e1e72abb1}"
SIXMER_TABLE="${SIXMER_TABLE:-${PROJECT_DIR}/amylogram-py/tests/fixtures/amylogram_sixmer_probabilities.bin}"
RUN_AMYPRED="${RUN_AMYPRED:-0}"
RUN_AMYLOGRAM_R="${RUN_AMYLOGRAM_R:-0}"
RUN_PROTEIN_FEATURES="${RUN_PROTEIN_FEATURES:-0}"

require_file AWS_PREPARED_ABUNDANCE "${AWS_PREPARED_ABUNDANCE:-}"
require_file AWS_PREPARED_VCF "${AWS_PREPARED_VCF:-}"
require_file AWS_REFERENCE_FASTA "${AWS_REFERENCE_FASTA:-}"
require_file REF_GENOME "${REF_GENOME:-}"
require_file REF_GTF "${REF_GTF:-}"
require_file SIXMER_TABLE "$SIXMER_TABLE"

RUN_ROOT="$(mktemp -d "${TMPDIR:-/private/tmp}/nf-aws-reference.XXXXXX")"
SAMPLES="${RUN_ROOT}/samples.aws_reference.csv"
RESULTS_DIR="${OUTDIR}/${SAMPLE_ID}/results_proteins"
AMYLOID_RESULTS_DIR="${OUTDIR}/${SAMPLE_ID}/results_amyloid"
FEATURE_RESULTS_DIR="${OUTDIR}/${SAMPLE_ID}/results_protein_features"

printf 'sample_id,kind,abundance,vcf,expected_proteins_dir,notes\n' > "$SAMPLES"
printf '%s,aws_reference,%s,%s,,AWS reference regression input.\n' \
  "$SAMPLE_ID" "$AWS_PREPARED_ABUNDANCE" "$AWS_PREPARED_VCF" >> "$SAMPLES"

docker build -t amyloid-proteome-nextflow:local docker/proteome-nextflow
docker build -t amylogram-py-nextflow:local -f docker/amylogram-py-nextflow/Dockerfile amylogram-py
if [[ "$RUN_AMYPRED" == "1" || "$RUN_AMYPRED" == "true" ]]; then
  docker build -t amypred-frl-nextflow:local docker/amypred-frl-nextflow
fi
if [[ "$RUN_AMYLOGRAM_R" == "1" || "$RUN_AMYLOGRAM_R" == "true" ]]; then
  docker build -t amylogram-r-nextflow:local docker/amylogram
fi
if [[ "$RUN_PROTEIN_FEATURES" == "1" || "$RUN_PROTEIN_FEATURES" == "true" ]]; then
  docker build -t protein-features-nextflow:local docker/protein-features-nextflow
fi

run_args=()
if [[ "$RUN_AMYPRED" == "1" || "$RUN_AMYPRED" == "true" ]]; then
  run_args+=(--run_amypred true)
fi
if [[ "$RUN_AMYLOGRAM_R" == "1" || "$RUN_AMYLOGRAM_R" == "true" ]]; then
  run_args+=(--run_amylogram_r true)
fi
if [[ "$RUN_PROTEIN_FEATURES" == "1" || "$RUN_PROTEIN_FEATURES" == "true" ]]; then
  run_args+=(--run_protein_features true)
fi

nextflow run nextflow/main.nf \
  -profile docker,aws_reference \
  --mode prepared \
  --samples "$SAMPLES" \
  --ref_genome "$REF_GENOME" \
  --ref_gtf "$REF_GTF" \
  --sixmer_table "$SIXMER_TABLE" \
  --outdir "$OUTDIR" \
  "${run_args[@]}"

test -s "${RESULTS_DIR}/combine_proteome.fasta"
grep -q OK "${RESULTS_DIR}/verification_report.tsv"

python3 tests_nextflow/compare_fasta_exact.py \
  --expected "$AWS_REFERENCE_FASTA" \
  --actual "${RESULTS_DIR}/combine_proteome.fasta" \
  --expected-count "$EXPECTED_RECORDS" \
  --expected-md5 "$EXPECTED_SEQUENCE_SET_MD5"

python3 tests_nextflow/validate_amyloid_outputs.py "$AMYLOID_RESULTS_DIR"

if [[ "$RUN_AMYPRED" == "1" || "$RUN_AMYPRED" == "true" ]]; then
  test -s "${AMYLOID_RESULTS_DIR}/amypred_frl_prediction.csv"
  grep -q $'AMYPred-FRL\tOK' "${AMYLOID_RESULTS_DIR}/amyloid_predictors_status.tsv"
fi
if [[ "$RUN_AMYLOGRAM_R" == "1" || "$RUN_AMYLOGRAM_R" == "true" ]]; then
  test -s "${AMYLOID_RESULTS_DIR}/amylogram_prediction_fast.csv"
  test -s "${AMYLOID_RESULTS_DIR}/amylogram_fast_summary.tsv"
  grep -q $'AmyloGram\tOK' "${AMYLOID_RESULTS_DIR}/amyloid_predictors_status.tsv"
fi

test -s "${FEATURE_RESULTS_DIR}/protein_features_light.csv"
grep -q OK "${FEATURE_RESULTS_DIR}/protein_features_light_status.tsv"
if [[ "$RUN_PROTEIN_FEATURES" == "1" || "$RUN_PROTEIN_FEATURES" == "true" ]]; then
  test -s "${FEATURE_RESULTS_DIR}/protein_features.csv"
  test -s "${FEATURE_RESULTS_DIR}/protein_features_status.tsv"
  grep -q OK "${FEATURE_RESULTS_DIR}/protein_features_status.tsv"
fi
