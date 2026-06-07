#!/usr/bin/env bash
set -Eeuo pipefail

SOURCE_S3="${SOURCE_S3:-s3://bioinfo-data-amylice-2026}"
DEST_S3="${DEST_S3:-s3://codex-test-ngsdata-calculations/prepared-bioinfo-proteome}"
AWS_PROFILE="${AWS_PROFILE:-codex-sandbox}"
AWS_REGION="${AWS_REGION:-us-east-1}"
PROJECT_DIR="${PROJECT_DIR:-/home/codex}"
PIPELINE_DIR="${PIPELINE_DIR:-${PROJECT_DIR}/pipeline_scripts}"
WORK_BASE="${WORK_BASE:-${PROJECT_DIR}/prepared_bioinfo_work}"
SAMPLES_FILE="${SAMPLES_FILE:-${WORK_BASE}/samples.prepared.txt}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
LOG_DIR="${LOG_DIR:-${PROJECT_DIR}/run_logs/prepared_bioinfo_proteome_${RUN_ID}}"
SKIP_EXISTING="${SKIP_EXISTING:-true}"

mkdir -p "$WORK_BASE" "$LOG_DIR"

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "${LOG_DIR}/batch.log"
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    log "Missing command: $1"
    exit 1
  }
}

fetch_samples() {
  if [[ -s "$SAMPLES_FILE" ]]; then
    log "Using existing samples file: ${SAMPLES_FILE}"
    return 0
  fi
  log "Downloading samples list from ${SOURCE_S3}/samples.txt"
  aws s3 cp "${SOURCE_S3}/samples.txt" "$SAMPLES_FILE" \
    --profile "$AWS_PROFILE" \
    --region "$AWS_REGION" \
    --quiet
}

sample_done_in_s3() {
  local sample_id="$1"
  aws s3 ls "${DEST_S3}/${sample_id}/results_proteins/UPLOAD_SUCCESS" \
    --profile "$AWS_PROFILE" \
    --region "$AWS_REGION" >/dev/null 2>&1
}

fetch_sample_inputs() {
  local sample_id="$1"
  local sample_dir="${WORK_BASE}/${sample_id}"
  mkdir -p "${sample_dir}/results_expression" "${sample_dir}/results_gatk"

  log "Fetching prepared inputs for ${sample_id}"
  aws s3 cp "${SOURCE_S3}/${sample_id}/results_expression/abundance.tsv" \
    "${sample_dir}/results_expression/abundance.tsv" \
    --profile "$AWS_PROFILE" \
    --region "$AWS_REGION" \
    --quiet

  if aws s3 ls "${SOURCE_S3}/${sample_id}/results_gatk/variants_filtered.vcf.gz" \
      --profile "$AWS_PROFILE" \
      --region "$AWS_REGION" >/dev/null 2>&1; then
    aws s3 cp "${SOURCE_S3}/${sample_id}/results_gatk/variants_filtered.vcf.gz" \
      "${sample_dir}/results_gatk/variants_filtered.vcf.gz" \
      --profile "$AWS_PROFILE" \
      --region "$AWS_REGION" \
      --quiet
  else
    aws s3 cp "${SOURCE_S3}/${sample_id}/results_gatk/variants_filtered.vcf" \
      "${sample_dir}/results_gatk/variants_filtered.vcf" \
      --profile "$AWS_PROFILE" \
      --region "$AWS_REGION" \
      --quiet
  fi
}

run_one_sample() {
  local sample_id="$1"
  local single_file="${WORK_BASE}/${sample_id}.samples.txt"
  printf '%s\n' "$sample_id" > "$single_file"

  fetch_sample_inputs "$sample_id"

  log "Starting proteome generation for ${sample_id}"
  SAMPLES_FILE="$single_file" \
  PROJECT_DIR="$PROJECT_DIR" \
  BASE_DIR="$WORK_BASE" \
  REF_DIR="${PROJECT_DIR}/references" \
  AWS_PROFILE="$AWS_PROFILE" \
  CLEANUP_AFTER_UPLOAD=true \
  PREPARE_REFS=false \
  PYTHONUNBUFFERED=1 \
    bash "${PIPELINE_DIR}/run_local_proteome.sh" "$DEST_S3"

  log "Finished proteome generation for ${sample_id}"
}

main() {
  require_cmd aws
  require_cmd python3
  require_cmd bcftools
  require_cmd samtools
  require_cmd gffread
  test -s "${PIPELINE_DIR}/run_local_proteome.sh"
  test -s "${PROJECT_DIR}/references/GRCh38.primary_assembly.genome.fa"
  test -s "${PROJECT_DIR}/references/gencode.v46.primary_assembly.annotation.gtf"

  fetch_samples

  while IFS= read -r raw_sample_id || [[ -n "$raw_sample_id" ]]; do
    sample_id="$(printf '%s' "$raw_sample_id" | tr -d '[:space:]')"
    [[ -z "$sample_id" || "$sample_id" == \#* ]] && continue

    if [[ "$SKIP_EXISTING" == "true" ]] && sample_done_in_s3 "$sample_id"; then
      log "SKIP ${sample_id}: existing UPLOAD_SUCCESS in ${DEST_S3}"
      continue
    fi

    if run_one_sample "$sample_id"; then
      log "OK ${sample_id}"
    else
      log "FAIL ${sample_id}; continuing with next sample"
    fi
  done < "$SAMPLES_FILE"

  log "Batch completed"
}

main "$@"
