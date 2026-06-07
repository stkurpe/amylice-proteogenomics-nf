#!/bin/bash

# Shared logging helpers for pipeline scripts.
# Source this file from a pipeline step, then call:
#   init_pipeline_logging "$SAMPLE_ID" "$LOG_DIR" "$STEP_NAME"
#   log_info "message"
#   run_logged "label" command arg1 arg2

pipeline_log_ts() {
    date -u +"%Y-%m-%dT%H:%M:%SZ"
}

init_pipeline_logging() {
    PIPELINE_SAMPLE_ID="${1:-unknown}"
    PIPELINE_LOG_DIR="${2:-logs}"
    PIPELINE_STEP_NAME="${3:-pipeline}"

    mkdir -p "$PIPELINE_LOG_DIR"
    PIPELINE_LOG_FILE="${PIPELINE_LOG_DIR}/${PIPELINE_STEP_NAME}.log"
    PIPELINE_EVENT_FILE="${PIPELINE_LOG_DIR}/${PIPELINE_STEP_NAME}.events.tsv"

    if [[ ! -f "$PIPELINE_EVENT_FILE" ]]; then
        printf "timestamp\tsample\tstep\tlevel\tmessage\n" > "$PIPELINE_EVENT_FILE"
    fi

    log_info "logging initialized: log=${PIPELINE_LOG_FILE}"
}

log_event() {
    local level="${1:-INFO}"
    shift || true
    local message="$*"
    local ts
    ts="$(pipeline_log_ts)"

    printf "[%s] [%s] [%s] %s\n" "$ts" "$level" "${PIPELINE_STEP_NAME:-pipeline}" "$message" | tee -a "${PIPELINE_LOG_FILE:-/dev/stderr}" >&2
    if [[ -n "${PIPELINE_EVENT_FILE:-}" ]]; then
        printf "%s\t%s\t%s\t%s\t%s\n" "$ts" "${PIPELINE_SAMPLE_ID:-unknown}" "${PIPELINE_STEP_NAME:-pipeline}" "$level" "$message" >> "$PIPELINE_EVENT_FILE"
    fi
}

log_info() {
    log_event "INFO" "$@"
}

log_warn() {
    log_event "WARN" "$@"
}

log_error() {
    log_event "ERROR" "$@"
}

log_step() {
    log_event "STEP" "$@"
}

run_logged() {
    local label="$1"
    shift
    local start_ts
    local end_ts
    local status

    start_ts="$(pipeline_log_ts)"
    log_step "START ${label}"

    if "$@" >> "${PIPELINE_LOG_FILE:-/dev/stderr}" 2>&1; then
        status=0
    else
        status=$?
    fi

    end_ts="$(pipeline_log_ts)"
    if [[ "$status" -eq 0 ]]; then
        log_step "OK ${label} start=${start_ts} end=${end_ts}"
    else
        log_error "FAIL ${label} status=${status} start=${start_ts} end=${end_ts}"
    fi

    return "$status"
}

write_pipeline_manifest() {
    local manifest_path="$1"
    shift

    mkdir -p "$(dirname "$manifest_path")"
    {
        printf "timestamp=%s\n" "$(pipeline_log_ts)"
        printf "sample=%s\n" "${PIPELINE_SAMPLE_ID:-unknown}"
        printf "step=%s\n" "${PIPELINE_STEP_NAME:-pipeline}"
        printf "log_file=%s\n" "${PIPELINE_LOG_FILE:-}"
        printf "events_file=%s\n" "${PIPELINE_EVENT_FILE:-}"
        while [[ "$#" -gt 0 ]]; do
            printf "%s\n" "$1"
            shift
        done
    } > "$manifest_path"
    log_info "manifest written: ${manifest_path}"
}
