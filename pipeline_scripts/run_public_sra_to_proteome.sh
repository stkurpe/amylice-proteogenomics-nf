#!/bin/bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$PROJECT_DIR"

source "${SCRIPT_DIR}/config_s3.sh"
source "${SCRIPT_DIR}/common.sh"
source "${SCRIPT_DIR}/pipeline_logging.sh"

SAMPLES_FILE="${SAMPLES_FILE:-samples.txt}"
S3_BUCKET="${1:-${S3_BUCKET:-s3://codex-test-ngsdata-calculations/results}}"
AWS_PROFILE="${AWS_PROFILE:-codex-sandbox}"
PREPARE_REFS="${PREPARE_REFS:-true}"
CLEANUP_AFTER_UPLOAD="${CLEANUP_AFTER_UPLOAD:-true}"
DRY_RUN="${DRY_RUN:-false}"

require_file "$SAMPLES_FILE" "samples list"

run_workflow_step() {
    local sample_id="$1"
    local label="$2"
    shift 2
    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "DRY_RUN: $*"
        return 0
    fi
    run_logged "$label" "$@"
}

if [[ "$PREPARE_REFS" == "true" ]]; then
    echo "Preparing references if missing..."
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "DRY_RUN: bash ${SCRIPT_DIR}/00_prepare_refs.sh"
    else
        bash "${SCRIPT_DIR}/00_prepare_refs.sh"
    fi
fi

while IFS= read -r raw_sample_id || [[ -n "$raw_sample_id" ]]; do
    SAMPLE_ID="$(printf '%s' "$raw_sample_id" | tr -d '[:space:]')"
    [[ -z "$SAMPLE_ID" || "$SAMPLE_ID" == \#* ]] && continue
    SAMPLE_ID="$(require_sample_id "$SAMPLE_ID")"

    SAMPLE_DIR="${BASE_DIR}/${SAMPLE_ID}"
    LOG_DIR="${SAMPLE_DIR}/logs"
    init_pipeline_logging "$SAMPLE_ID" "$LOG_DIR" "public_sra_to_proteome"

    echo "=================================================="
    echo "PUBLIC SRA TO MUTANT PROTEOME: ${SAMPLE_ID}"
    echo "=================================================="
    log_info "sample=${SAMPLE_ID}"
    log_info "output_s3=${S3_BUCKET}"
    log_info "aws_profile=${AWS_PROFILE}"
    log_info "cleanup_after_upload=${CLEANUP_AFTER_UPLOAD}"
    log_info "dry_run=${DRY_RUN}"

    run_workflow_step "$SAMPLE_ID" "download FASTQ from SRA" bash "${SCRIPT_DIR}/01_download.sh" "$SAMPLE_ID"
    run_workflow_step "$SAMPLE_ID" "FastQC" bash "${SCRIPT_DIR}/02_qc.sh" "$SAMPLE_ID"
    run_workflow_step "$SAMPLE_ID" "Kallisto expression" bash "${SCRIPT_DIR}/03_expression.sh" "$SAMPLE_ID"
    run_workflow_step "$SAMPLE_ID" "STAR and GATK variants" bash "${SCRIPT_DIR}/04_variants.sh" "$SAMPLE_ID"

    SAMPLES_FILE_TMP="${SAMPLE_DIR}/single_sample.txt"
    mkdir -p "$SAMPLE_DIR"
    printf "%s\n" "$SAMPLE_ID" > "$SAMPLES_FILE_TMP"
    log_info "single-sample file: ${SAMPLES_FILE_TMP}"

    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "DRY_RUN: SAMPLES_FILE=${SAMPLES_FILE_TMP} S3_BUCKET=${S3_BUCKET} bash ${SCRIPT_DIR}/run_local_proteome.sh"
    else
        SAMPLES_FILE="$SAMPLES_FILE_TMP" \
        S3_BUCKET="$S3_BUCKET" \
        AWS_PROFILE="$AWS_PROFILE" \
        CLEANUP_AFTER_UPLOAD="$CLEANUP_AFTER_UPLOAD" \
            bash "${SCRIPT_DIR}/run_local_proteome.sh" "$S3_BUCKET"
    fi
done < "$SAMPLES_FILE"
