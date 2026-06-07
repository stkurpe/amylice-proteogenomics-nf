#!/usr/bin/env bash
set -Eeuo pipefail

AWS_PROFILE="${AWS_PROFILE:-codex-sandbox}"
AWS_REGION="${AWS_REGION:-us-east-1}"
PROTEOME_PREFIX="${PROTEOME_PREFIX:-s3://codex-test-ngsdata-calculations/prepared-bioinfo-proteome}"
SAMPLES_FILE="${SAMPLES_FILE:-/home/codex/prepared_bioinfo_work/samples.prepared.txt}"
LOG_ROOT="${LOG_ROOT:-/home/codex/run_logs/amyloid_retry_$(date -u +%Y%m%dT%H%M%SZ)}"
AMYLOID_SCRIPT="${AMYLOID_SCRIPT:-/home/codex/run_amyloid_predictors.sh}"
AMYLOGRAM_CORES="${AMYLOGRAM_CORES:-8}"
AMYLOGRAM_CHUNK_SIZE="${AMYLOGRAM_CHUNK_SIZE:-50}"

mkdir -p "$LOG_ROOT"
MAIN_LOG="$LOG_ROOT/main.log"
STATUS_TSV="$LOG_ROOT/status.tsv"
printf 'time_utc\tsample\tstatus\tdetail\n' > "$STATUS_TSV"

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$MAIN_LOG"
}

status() {
  printf '%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" "$2" "$3" >> "$STATUS_TSV"
}

s3_exists() {
  aws s3 ls "$1" --profile "$AWS_PROFILE" --region "$AWS_REGION" >/dev/null 2>&1
}

while pgrep -f "bash /home/codex/run_amyloid_predictors.sh" >/dev/null 2>&1; do
  log "Waiting for active amyloid predictor process to finish"
  sleep 300
done

while IFS= read -r sample; do
  [[ -z "$sample" || "$sample" == \#* ]] && continue
  marker="${PROTEOME_PREFIX}/${sample}/results_proteins/UPLOAD_SUCCESS"
  fasta="${PROTEOME_PREFIX}/${sample}/results_proteins/combine_proteome.fasta"
  output="${PROTEOME_PREFIX}/${sample}/results_amyloid"
  combined="${output}/amyloid_combined_predictions.csv"

  if ! s3_exists "$marker"; then
    log "SKIP ${sample}: missing proteome marker"
    status "$sample" "SKIP" "missing proteome marker"
    continue
  fi
  if ! s3_exists "$fasta"; then
    log "SKIP ${sample}: missing combine FASTA"
    status "$sample" "SKIP" "missing fasta"
    continue
  fi
  if s3_exists "$combined"; then
    log "SKIP ${sample}: combined exists"
    status "$sample" "SKIP" "combined exists"
    continue
  fi

  log "START ${sample}"
  status "$sample" "START" "$fasta"
  if AMYLOGRAM_CORES="$AMYLOGRAM_CORES" AMYLOGRAM_CHUNK_SIZE="$AMYLOGRAM_CHUNK_SIZE" \
    AWS_PROFILE="$AWS_PROFILE" AWS_REGION="$AWS_REGION" \
    "$AMYLOID_SCRIPT" --sample "$sample" --input-s3 "$fasta" --output-s3 "$output" > "$LOG_ROOT/${sample}.log" 2>&1; then
    log "OK ${sample}"
    status "$sample" "OK" "$output"
  else
    rc=$?
    log "FAIL ${sample} rc=${rc}"
    status "$sample" "FAIL" "rc=${rc} log=${LOG_ROOT}/${sample}.log"
  fi
done < "$SAMPLES_FILE"

log "Retry queue finished"
