#!/bin/bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/common.sh"

SAMPLE_ID="$(require_sample_id "${1:-}")"
FASTQ="${BASE_DIR}/${SAMPLE_ID}/downloaded/${SAMPLE_ID}.fastq.gz"

require_file "$FASTQ" "FASTQ"

echo ">>> [STEP 2] QC: ${SAMPLE_ID}"
mkdir -p "${BASE_DIR}/${SAMPLE_ID}/results_qc"

docker run $DOCKER_OPTS \
    quay.io/biocontainers/fastqc:0.11.9--0 \
    fastqc -o "/data/${SAMPLE_ID}/results_qc" -t 4 "/data/${SAMPLE_ID}/downloaded/${SAMPLE_ID}.fastq.gz"
