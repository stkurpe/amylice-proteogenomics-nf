#!/usr/bin/env bash
set -Eeuo pipefail

AWS_PROFILE="${AWS_PROFILE:-codex-sandbox}"
AWS_REGION="${AWS_REGION:-us-east-1}"
HISAT_INDEX_DIR="${HISAT_INDEX_DIR:-/home/codex/tools/src/hisat-genotype/indicies}"
HISAT_INDEX_URL="${HISAT_INDEX_URL:-ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat-genotype/data/genotype_genome_20180128.tar.gz}"
S3_CACHE_URI="${S3_CACHE_URI:-s3://codex-test-ngsdata-calculations/tool-cache/hisat-genotype/genotype_genome_20180128/}"
THREADS="${THREADS:-4}"

log() {
    printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"
}

index_ready() {
    [[ -s "${HISAT_INDEX_DIR}/genotype_genome.fa" ]] || return 1
    find "${HISAT_INDEX_DIR}" -maxdepth 1 -type f \( -name 'genotype_genome*.ht2' -o -name 'genotype_genome*.ht2l' \) | grep -q .
}

mkdir -p "${HISAT_INDEX_DIR}"

if index_ready; then
    log "HISAT genotype_genome index already exists in ${HISAT_INDEX_DIR}"
else
    staging="${HISAT_INDEX_DIR}/.genotype_genome_staging_$(date -u +%Y%m%dT%H%M%SZ)"
    mkdir -p "${staging}"
    log "Downloading and extracting ${HISAT_INDEX_URL} into ${staging}"
    curl -L --fail --show-error "${HISAT_INDEX_URL}" | tar -xz -C "${staging}"

    log "Moving genotype_genome files into ${HISAT_INDEX_DIR}"
    find "${staging}" -maxdepth 1 -type f -name 'genotype_genome*' -exec mv {} "${HISAT_INDEX_DIR}/" \;
fi

if index_ready; then
    log "Index is ready; listing files"
    find "${HISAT_INDEX_DIR}" -maxdepth 1 -type f -name 'genotype_genome*' -printf '%f\t%s\n' | sort
else
    log "Index is still incomplete"
    exit 1
fi

log "Uploading HISAT index cache to ${S3_CACHE_URI}"
aws s3 cp "${HISAT_INDEX_DIR}/" "${S3_CACHE_URI}" \
    --recursive \
    --exclude '*' \
    --include 'genotype_genome*' \
    --profile "${AWS_PROFILE}" \
    --region "${AWS_REGION}" \
    --no-progress

log "DONE"
