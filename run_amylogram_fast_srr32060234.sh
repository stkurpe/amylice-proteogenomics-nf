#!/usr/bin/env bash
set -Eeuo pipefail

SAMPLE_ID="${SAMPLE_ID:-SRR32060234}"
WORK_DIR="${WORK_DIR:-/home/codex/amyloid_runs/${SAMPLE_ID}}"
IMAGE="${IMAGE:-codex/amylogram:1.1-r4.3.3}"
AWS_PROFILE="${AWS_PROFILE:-codex-sandbox}"
AWS_REGION="${AWS_REGION:-us-east-1}"
S3_PREFIX="${S3_PREFIX:-s3://codex-test-ngsdata-calculations/public-sra-test/${SAMPLE_ID}/results_amyloid}"
AMYLOGRAM_CORES="${AMYLOGRAM_CORES:-7}"
AMYLOGRAM_CHUNK_SIZE="${AMYLOGRAM_CHUNK_SIZE:-250}"

INPUT_FASTA="${INPUT_FASTA:-${WORK_DIR}/input/combine_proteome.fasta}"
CODE_DIR="${CODE_DIR:-${WORK_DIR}/code/AmyloGram}"
RESULTS_DIR="${RESULTS_DIR:-${WORK_DIR}/results}"
LOG_DIR="${LOG_DIR:-${WORK_DIR}/logs}"
CHUNK_DIR="${CHUNK_DIR:-${RESULTS_DIR}/amylogram_fast_chunks}"
OUTPUT_CSV="${OUTPUT_CSV:-${RESULTS_DIR}/amylogram_prediction_fast.csv}"
SUMMARY_TSV="${SUMMARY_TSV:-${RESULTS_DIR}/amylogram_fast_summary.tsv}"
LOG_FILE="${LOG_FILE:-${LOG_DIR}/amylogram_fast_$(date -u +%Y%m%dT%H%M%SZ).log}"

mkdir -p "$RESULTS_DIR" "$LOG_DIR" "$CHUNK_DIR"

{
  echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] Starting fast AmyloGram run"
  echo "sample=${SAMPLE_ID}"
  echo "image=${IMAGE}"
  echo "input=${INPUT_FASTA}"
  echo "output=${OUTPUT_CSV}"
  echo "summary=${SUMMARY_TSV}"
  echo "chunks=${CHUNK_DIR}"
  echo "cores=${AMYLOGRAM_CORES}"
  echo "chunk_size=${AMYLOGRAM_CHUNK_SIZE}"

  test -s "$INPUT_FASTA"
  test -s "${CODE_DIR}/predict_amylogram_fast.R"

  docker run --rm \
    --name "amylogram-fast-${SAMPLE_ID}" \
    -e AMYLOGRAM_CORES="${AMYLOGRAM_CORES}" \
    -e AMYLOGRAM_CHUNK_SIZE="${AMYLOGRAM_CHUNK_SIZE}" \
    -v "${WORK_DIR}:${WORK_DIR}" \
    -w "${WORK_DIR}" \
    "${IMAGE}" \
    Rscript "${CODE_DIR}/predict_amylogram_fast.R" \
      "${INPUT_FASTA}" \
      "${OUTPUT_CSV}" \
      "${CHUNK_DIR}" \
      "${SUMMARY_TSV}"

  echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] Docker run finished"
  wc -l "$OUTPUT_CSV" "$SUMMARY_TSV"

  aws s3 cp "$OUTPUT_CSV" "${S3_PREFIX}/amylogram_prediction_fast.csv" \
    --profile "$AWS_PROFILE" \
    --region "$AWS_REGION" \
    --quiet

  aws s3 cp "$SUMMARY_TSV" "${S3_PREFIX}/amylogram_fast_summary.tsv" \
    --profile "$AWS_PROFILE" \
    --region "$AWS_REGION" \
    --quiet

  aws s3 cp "${CHUNK_DIR}/chunk_stats.tsv" "${S3_PREFIX}/amylogram_fast_chunk_stats.tsv" \
    --profile "$AWS_PROFILE" \
    --region "$AWS_REGION" \
    --quiet

  aws s3 cp "$LOG_FILE" "${S3_PREFIX}/amylogram_fast.log" \
    --profile "$AWS_PROFILE" \
    --region "$AWS_REGION" \
    --quiet

  echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] Uploaded fast AmyloGram outputs to ${S3_PREFIX}"
} 2>&1 | tee "$LOG_FILE"
