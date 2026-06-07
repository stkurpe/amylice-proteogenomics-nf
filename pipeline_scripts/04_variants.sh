#!/bin/bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/common.sh"

SAMPLE_ID="$(require_sample_id "${1:-}")"
S_DIR="${BASE_DIR}/${SAMPLE_ID}"
FASTQ="${S_DIR}/downloaded/${SAMPLE_ID}.fastq.gz"
GENOME="${REF_DIR}/${REF_GENOME}"

require_file "$FASTQ" "FASTQ"
require_file "$GENOME" "reference genome"

echo ">>> [STEP 4] STAR & GATK: ${SAMPLE_ID}"
mkdir -p "$S_DIR/results_star" "$S_DIR/results_gatk"

if [[ ! -f "$S_DIR/results_star/Aligned.sortedByCoord.out.bam" ]]; then
    echo "--- STAR Mapping..."
    docker run $DOCKER_OPTS \
        quay.io/biocontainers/star:2.7.11b--h5ca1c30_8 \
        STAR --runThreadN "$THREADS" \
        --genomeDir /data/references/star_index_full \
        --readFilesIn "/data/${SAMPLE_ID}/downloaded/${SAMPLE_ID}.fastq.gz" \
        --readFilesCommand zcat \
        --outFileNamePrefix "/data/${SAMPLE_ID}/results_star/" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattrRGline "ID:${SAMPLE_ID}" "SM:${SAMPLE_ID}" PL:ILLUMINA
fi

if [[ ! -f "$S_DIR/results_star/dedup.bam" ]]; then
    echo "--- MarkDuplicates..."
    docker run $DOCKER_OPTS \
        quay.io/biocontainers/picard:2.27.5--hdfd78af_0 \
        picard MarkDuplicates \
        "I=/data/${SAMPLE_ID}/results_star/Aligned.sortedByCoord.out.bam" \
        "O=/data/${SAMPLE_ID}/results_star/dedup.bam" \
        "M=/data/${SAMPLE_ID}/results_star/metrics.txt"
    docker run $DOCKER_OPTS quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1 \
        samtools index "/data/${SAMPLE_ID}/results_star/dedup.bam"
fi

if [[ ! -f "$S_DIR/results_gatk/split.bam" ]]; then
    echo "--- SplitNCigarReads..."
    docker run $DOCKER_OPTS \
        broadinstitute/gatk:4.5.0.0 \
        gatk SplitNCigarReads \
        -R "/data/references/${REF_GENOME}" \
        -I "/data/${SAMPLE_ID}/results_star/dedup.bam" \
        -O "/data/${SAMPLE_ID}/results_gatk/split.bam"
fi

if [[ ! -f "$S_DIR/results_gatk/variants.vcf.gz" ]]; then
    echo "--- HaplotypeCaller..."
    docker run $DOCKER_OPTS \
        broadinstitute/gatk:4.5.0.0 \
        gatk HaplotypeCaller \
        -R "/data/references/${REF_GENOME}" \
        -I "/data/${SAMPLE_ID}/results_gatk/split.bam" \
        -O "/data/${SAMPLE_ID}/results_gatk/variants.vcf.gz" \
        --dont-use-soft-clipped-bases \
        -stand-call-conf 20.0
fi

if [[ ! -f "$S_DIR/results_gatk/variants_filtered.vcf.gz" ]]; then
    echo "--- VariantFiltration..."
    docker run $DOCKER_OPTS \
        broadinstitute/gatk:4.5.0.0 \
        gatk VariantFiltration \
        -R "/data/references/${REF_GENOME}" \
        -V "/data/${SAMPLE_ID}/results_gatk/variants.vcf.gz" \
        -O "/data/${SAMPLE_ID}/results_gatk/variants_filtered.vcf.gz" \
        --window 35 --cluster 3 \
        --filter-name "FS" --filter-expression "FS > 30.0" \
        --filter-name "QD" --filter-expression "QD < 2.0"
fi
