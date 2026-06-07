#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SAMPLES_FILE="${SAMPLES_FILE:-${SCRIPT_DIR}/samples.txt}"
S3_BUCKET="${S3_BUCKET:-s3://codex-test-ngsdata-calculations/public-sra-test}"
AWS_PROFILE="${AWS_PROFILE:-codex-sandbox}"
AWS_REGION="${AWS_REGION:-us-east-1}"
PROJECT_DIR="${PROJECT_DIR:-/home/codex}"
RUN_AMYLOID="${RUN_AMYLOID:-true}"
RUN_PROTEIN_FEATURES="${RUN_PROTEIN_FEATURES:-true}"
FORCE_AMYLOID="${FORCE_AMYLOID:-0}"

usage() {
  cat <<EOF
Usage:
  ./run_full_pipeline.sh --samples samples.txt --s3-bucket s3://bucket/prefix

Options:
  --samples PATH          File with one SRA accession per line.
  --s3-bucket URI         Output S3 prefix.
  --aws-profile NAME      AWS CLI profile name.
  --aws-region REGION     AWS region.
  --project-dir PATH      Working directory on the compute server.
  --skip-amyloid          Stop after mutant proteome generation.
  --skip-protein-features Do not calculate R protein physicochemical features.
  --force-amyloid         Recompute amyloid predictors even if outputs exist.
  -h, --help              Show this help.

Environment:
  THREADS=8 CLEANUP_AFTER_UPLOAD=true PREPARE_REFS=true
  AMYLOGRAM_CORES=8 AMYLOGRAM_CHUNK_SIZE=50
  PROTEIN_FEATURES_IMAGE=codex/amylogram:1.1-r4.3.3 PROTEIN_FEATURES_R_LIBS=/path/to/R/library
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples)
      SAMPLES_FILE="$2"
      shift 2
      ;;
    --s3-bucket)
      S3_BUCKET="$2"
      shift 2
      ;;
    --aws-profile)
      AWS_PROFILE="$2"
      shift 2
      ;;
    --aws-region)
      AWS_REGION="$2"
      shift 2
      ;;
    --project-dir)
      PROJECT_DIR="$2"
      shift 2
      ;;
    --skip-amyloid)
      RUN_AMYLOID="false"
      shift
      ;;
    --skip-protein-features)
      RUN_PROTEIN_FEATURES="false"
      shift
      ;;
    --force-amyloid)
      FORCE_AMYLOID=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

test -s "$SAMPLES_FILE"

export SAMPLES_FILE
export S3_BUCKET
export AWS_PROFILE
export AWS_REGION
export PROJECT_DIR

bash "${SCRIPT_DIR}/pipeline_scripts/run_public_sra_to_proteome.sh" "$S3_BUCKET"

if [[ "$RUN_AMYLOID" != "true" ]]; then
  exit 0
fi

while IFS= read -r raw_sample_id || [[ -n "$raw_sample_id" ]]; do
  sample_id="$(printf '%s' "$raw_sample_id" | tr -d '[:space:]')"
  [[ -z "$sample_id" || "$sample_id" == \#* ]] && continue

  amyloid_args=(
    --sample "$sample_id"
    --output-s3 "${S3_BUCKET}/${sample_id}/results_amyloid"
    --protein-features-output-s3 "${S3_BUCKET}/${sample_id}/results_protein_features"
  )
  if [[ "$FORCE_AMYLOID" == "1" ]]; then
    amyloid_args+=(--force)
  fi
  if [[ "$RUN_PROTEIN_FEATURES" != "true" ]]; then
    amyloid_args+=(--skip-protein-features)
  fi

  bash "${SCRIPT_DIR}/run_amyloid_predictors.sh" "${amyloid_args[@]}"
done < "$SAMPLES_FILE"
