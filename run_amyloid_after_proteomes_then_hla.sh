#!/usr/bin/env bash
set -Eeuo pipefail

AWS_PROFILE="${AWS_PROFILE:-codex-sandbox}"
AWS_REGION="${AWS_REGION:-us-east-1}"
PROTEOME_PREFIX="${PROTEOME_PREFIX:-s3://codex-test-ngsdata-calculations/prepared-bioinfo-proteome}"
SAMPLES_FILE="${SAMPLES_FILE:-/home/codex/prepared_bioinfo_work/samples.prepared.txt}"
LOG_ROOT="${LOG_ROOT:-/home/codex/run_logs/amyloid_then_hla_$(date -u +%Y%m%dT%H%M%SZ)}"
AMYLOID_SCRIPT="${AMYLOID_SCRIPT:-/home/codex/run_amyloid_predictors.sh}"
HLA_SCRIPT="${HLA_SCRIPT:-/home/codex/amyloid-mutant-proteome-pipeline/pipeline_scripts/05_hla_consensus.sh}"
HLA_DEBUG_SAMPLE="${HLA_DEBUG_SAMPLE:-SRR32060234}"
AMYLOGRAM_CORES="${AMYLOGRAM_CORES:-8}"
AMYLOGRAM_CHUNK_SIZE="${AMYLOGRAM_CHUNK_SIZE:-50}"
WAIT_SECONDS="${WAIT_SECONDS:-300}"

mkdir -p "$LOG_ROOT"
MAIN_LOG="$LOG_ROOT/main.log"
STATUS_TSV="$LOG_ROOT/status.tsv"
printf 'time_utc\tsample\tstage\tstatus\tdetail\n' > "$STATUS_TSV"

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$MAIN_LOG"
}

status() {
  printf '%s\t%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" "$2" "$3" "$4" >> "$STATUS_TSV"
}

s3_exists() {
  aws s3 ls "$1" --profile "$AWS_PROFILE" --region "$AWS_REGION" >/dev/null 2>&1
}

wait_for_proteome_batch() {
  log "Waiting for prepared proteome batch to finish"
  while pgrep -f "bash /home/codex/run_prepared_proteomes_from_bioinfo.sh" >/dev/null 2>&1; do
    pid_line="$(pgrep -af "bash /home/codex/run_prepared_proteomes_from_bioinfo.sh" | tr '\n' '; ')"
    log "Proteome batch still running: ${pid_line}"
    sleep "$WAIT_SECONDS"
  done
  log "Prepared proteome batch is not running; starting amyloidogenicity queue"
}

run_amyloid_sample() {
  local sample="$1"
  local fasta_s3="${PROTEOME_PREFIX}/${sample}/results_proteins/combine_proteome.fasta"
  local marker_s3="${PROTEOME_PREFIX}/${sample}/results_proteins/UPLOAD_SUCCESS"
  local output_s3="${PROTEOME_PREFIX}/${sample}/results_amyloid"
  local combined_s3="${output_s3}/amyloid_combined_predictions.csv"
  local sample_log="$LOG_ROOT/${sample}.amyloid.log"

  if ! s3_exists "$marker_s3"; then
    log "SKIP ${sample}: proteome UPLOAD_SUCCESS is missing"
    status "$sample" "amyloid" "SKIP" "missing proteome marker"
    return 0
  fi
  if ! s3_exists "$fasta_s3"; then
    log "SKIP ${sample}: combine_proteome.fasta is missing"
    status "$sample" "amyloid" "SKIP" "missing combine_proteome.fasta"
    return 0
  fi
  if s3_exists "$combined_s3"; then
    log "SKIP ${sample}: amyloid combined output already exists"
    status "$sample" "amyloid" "SKIP" "already exists"
    return 0
  fi

  log "START amyloidogenicity ${sample}"
  status "$sample" "amyloid" "START" "$fasta_s3"
  if AMYLOGRAM_CORES="$AMYLOGRAM_CORES" AMYLOGRAM_CHUNK_SIZE="$AMYLOGRAM_CHUNK_SIZE" \
    AWS_PROFILE="$AWS_PROFILE" AWS_REGION="$AWS_REGION" \
    "$AMYLOID_SCRIPT" \
      --sample "$sample" \
      --input-s3 "$fasta_s3" \
      --output-s3 "$output_s3" \
      > "$sample_log" 2>&1; then
    log "OK amyloidogenicity ${sample}"
    status "$sample" "amyloid" "OK" "$output_s3"
  else
    local rc=$?
    log "FAIL amyloidogenicity ${sample}; rc=${rc}; see ${sample_log}"
    status "$sample" "amyloid" "FAIL" "rc=${rc}; log=${sample_log}"
  fi
}

run_hla_debug() {
  local sample="$HLA_DEBUG_SAMPLE"
  local sample_log="$LOG_ROOT/${sample}.hla_debug.log"

  log "START HLA debug ${sample}"
  status "$sample" "hla_debug" "START" "$HLA_SCRIPT"
  if SOURCE_S3=s3://bioinfo-data-amylice-2026 \
    DEST_S3="$PROTEOME_PREFIX" \
    AWS_PROFILE="$AWS_PROFILE" AWS_REGION="$AWS_REGION" \
    bash "$HLA_SCRIPT" "$sample" > "$sample_log" 2>&1; then
    log "OK HLA debug ${sample}"
    status "$sample" "hla_debug" "OK" "${PROTEOME_PREFIX}/${sample}/results_hla/"
  else
    local rc=$?
    log "FAIL HLA debug ${sample}; rc=${rc}; see ${sample_log}"
    status "$sample" "hla_debug" "FAIL" "rc=${rc}; log=${sample_log}"
  fi
}

main() {
  log "Workflow log root: ${LOG_ROOT}"
  wait_for_proteome_batch

  while IFS= read -r sample; do
    [[ -z "$sample" || "$sample" == \#* ]] && continue
    run_amyloid_sample "$sample"
  done < "$SAMPLES_FILE"

  log "Amyloidogenicity queue finished; starting one-sample HLA debug"
  run_hla_debug
  log "Workflow finished"
}

main "$@"
