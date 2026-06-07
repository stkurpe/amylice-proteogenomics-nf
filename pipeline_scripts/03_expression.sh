#!/bin/bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/common.sh"

SAMPLE_ID="$(require_sample_id "${1:-}")"
FASTQ="${BASE_DIR}/${SAMPLE_ID}/downloaded/${SAMPLE_ID}.fastq.gz"
INDEX="${REF_DIR}/transcriptome.idx"

require_file "$FASTQ" "FASTQ"
require_file "$INDEX" "kallisto index"

echo ">>> [STEP 3] KALLISTO: ${SAMPLE_ID}"
mkdir -p "${BASE_DIR}/${SAMPLE_ID}/results_expression"

docker run $DOCKER_OPTS \
    quay.io/biocontainers/kallisto:0.46.1--h4f7b962_0 \
    kallisto quant \
    -i /data/references/transcriptome.idx \
    -o "/data/${SAMPLE_ID}/results_expression" \
    --single -l 200 -s 20 \
    "/data/${SAMPLE_ID}/downloaded/${SAMPLE_ID}.fastq.gz"
