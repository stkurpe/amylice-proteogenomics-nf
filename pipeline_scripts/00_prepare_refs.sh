#!/bin/bash
set -e
source $(dirname "$0")/config.sh

echo ">>> [SETUP] ПРОВЕРКА И ПОДГОТОВКА РЕФЕРЕНСОВ..."
mkdir -p "$REF_DIR"
mkdir -p "$REF_DIR/star_index_full"

if [ ! -f "$REF_DIR/$REF_GENOME" ]; then
    echo "--- Скачивание генома..."
    wget -O "$REF_DIR/$REF_GENOME.gz" https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
    gzip -d "$REF_DIR/$REF_GENOME.gz"
fi

if [ ! -f "$REF_DIR/$REF_GTF" ]; then
    echo "--- Скачивание GTF..."
    wget -O "$REF_DIR/$REF_GTF.gz" https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz
    gzip -d "$REF_DIR/$REF_GTF.gz"
fi

if [ ! -f "$REF_DIR/$REF_VCF" ]; then
    echo "--- Скачивание dbSNP..."
    wget -O "$REF_DIR/$REF_VCF.gz" https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
    gzip -d "$REF_DIR/$REF_VCF.gz"
    wget -nc https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz.tbi -O "$REF_DIR/$REF_VCF.gz.tbi"
fi

if [ ! -f "$REF_DIR/$REF_TRANSCRIPT" ]; then
    echo "--- Скачивание транскриптома..."
    wget -O "$REF_DIR/$REF_TRANSCRIPT.gz" https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.pc_transcripts.fa.gz
    gzip -d "$REF_DIR/$REF_TRANSCRIPT.gz"
fi

if [ ! -f "$REF_DIR/transcriptome.idx" ]; then
    echo "--- Генерация индекса Kallisto..."
    docker run $DOCKER_OPTS \
        quay.io/biocontainers/kallisto:0.46.1--h4f7b962_0 \
        kallisto index -i /data/references/transcriptome.idx /data/references/$REF_TRANSCRIPT
fi

if [ -z "$(ls -A $REF_DIR/star_index_full)" ]; then
    echo "--- Генерация индекса STAR (занимает время)..."
    docker run $DOCKER_OPTS \
        quay.io/biocontainers/star:2.7.11b--h5ca1c30_8 \
        STAR --runMode genomeGenerate \
        --runThreadN $THREADS \
        --genomeDir /data/references/star_index_full \
        --genomeFastaFiles /data/references/$REF_GENOME \
        --sjdbGTFfile /data/references/$REF_GTF \
        --sjdbOverhang 99 \
        --genomeSAsparseD 2
fi

if [ ! -f "$REF_DIR/${REF_GENOME%.fa}.dict" ]; then
    docker run $DOCKER_OPTS quay.io/biocontainers/picard:2.27.5--hdfd78af_0 \
        picard CreateSequenceDictionary R=/data/references/$REF_GENOME O=/data/references/${REF_GENOME%.fa}.dict
fi
if [ ! -f "$REF_DIR/$REF_GENOME.fai" ]; then
    docker run $DOCKER_OPTS quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1 \
        samtools faidx /data/references/$REF_GENOME
fi
