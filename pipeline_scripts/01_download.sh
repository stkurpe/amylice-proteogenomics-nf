#!/bin/bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/common.sh"

SAMPLE_ID="$(require_sample_id "${1:-}")"
SAMPLE_DIR="${BASE_DIR}/${SAMPLE_ID}"
FASTQ_GZ="${SAMPLE_DIR}/downloaded/${SAMPLE_ID}.fastq.gz"

echo ">>> [STEP 1] DOWNLOAD: ${SAMPLE_ID}"
mkdir -p "${SAMPLE_DIR}/downloaded"

if [[ -f "$FASTQ_GZ" ]]; then
    echo "FASTQ already exists: ${FASTQ_GZ}"
    exit 0
fi

echo "--- Running fasterq-dump..."
docker run $DOCKER_OPTS \
    staphb/sratoolkit:3.3.0 \
    fasterq-dump --split-files --outdir "/data/${SAMPLE_ID}/downloaded" --threads "$THREADS" --progress "$SAMPLE_ID"

cd "${SAMPLE_DIR}/downloaded"

if [[ -f "${SAMPLE_ID}_1.fastq" ]]; then
    mv "${SAMPLE_ID}_1.fastq" "${SAMPLE_ID}.fastq"
fi

if [[ -f "${SAMPLE_ID}_2.fastq" ]]; then
    echo "WARNING: paired-end mate 2 detected and preserved as ${SAMPLE_ID}_2.fastq.gz; downstream steps currently use mate 1 only."
    gzip -f "${SAMPLE_ID}_2.fastq"
fi

if [[ -f "${SAMPLE_ID}.fastq" ]]; then
    gzip -f "${SAMPLE_ID}.fastq"
else
    pipeline_die "FASTQ file was not produced for ${SAMPLE_ID}"
fi
