#!/bin/bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$PROJECT_DIR"

source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/common.sh"
source "${SCRIPT_DIR}/pipeline_logging.sh"

SAMPLES_FILE="${SAMPLES_FILE:-samples.txt}"
DRY_RUN="${DRY_RUN:-false}"

require_file "$SAMPLES_FILE" "samples list"

run_or_echo() {
    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "DRY_RUN: $*"
        return 0
    fi
    run_logged "$*" "$@"
}

while IFS= read -r raw_sample_id || [[ -n "$raw_sample_id" ]]; do
    SAMPLE_ID="$(printf '%s' "$raw_sample_id" | tr -d '[:space:]')"
    [[ -z "$SAMPLE_ID" || "$SAMPLE_ID" == \#* ]] && continue
    SAMPLE_ID="$(require_sample_id "$SAMPLE_ID")"

    SAMPLE_DIR="${BASE_DIR}/${SAMPLE_ID}"
    LOG_DIR="${SAMPLE_DIR}/logs"
    init_pipeline_logging "$SAMPLE_ID" "$LOG_DIR" "download_public_inputs"

    log_info "public source: SRA accession ${SAMPLE_ID}"
    log_info "output dir: ${SAMPLE_DIR}/downloaded"
    run_or_echo bash "${SCRIPT_DIR}/01_download.sh" "$SAMPLE_ID"
done < "$SAMPLES_FILE"
