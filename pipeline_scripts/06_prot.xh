cat << 'EOF' > pipeline_scripts/05_proteo.sh
#!/bin/bash
set -e

# === НАСТРОЙКИ ===
SAMPLE_ID=$1
BASE_DIR="/home/ubuntu"
SAMPLE_DIR="${BASE_DIR}/${SAMPLE_ID}"
REF_DIR="${BASE_DIR}/references"
PROTEO_DIR="${SAMPLE_DIR}/results_proteins"
EXPRESSION_FILE="${SAMPLE_DIR}/results_expression/abundance.tsv"
VCF_FILE="${SAMPLE_DIR}/results_gatk/variants_filtered.vcf.gz"

# Создаем папку заново
mkdir -p "$PROTEO_DIR"
# Пишем лог
exec > "${SAMPLE_DIR}/logs/05_proteo.log" 2>&1

echo "🚀 [PROTEO] Начинаем создание белков для $SAMPLE_ID"

# 1. Фильтруем экспрессирующиеся гены (TPM > 1)
echo "--- Фильтрация генов (TPM > 1)..."
awk -F "\t" 'NR>1 && $5 > 1 {print $1}' "$EXPRESSION_FILE" > "${PROTEO_DIR}/clean_ids.txt"

COUNT=$(wc -l < "${PROTEO_DIR}/clean_ids.txt")
echo "✅ Найдено активных транскриптов: $COUNT"

# 2. Создаем консенсусный геном (включаем мутации из VCF)
echo "--- Создание консенсусного генома (bcftools)..."
# Создаем индекс для референса, если его нет (нужен для bcftools)
if [ ! -f "${REF_DIR}/GRCh38.primary_assembly.genome.fa.fai" ]; then
    samtools faidx "${REF_DIR}/GRCh38.primary_assembly.genome.fa"
fi

bcftools consensus -f "${REF_DIR}/GRCh38.primary_assembly.genome.fa" "$VCF_FILE" > "${PROTEO_DIR}/consensus_genome.fa"

# 3. Трансляция в белки (gffread)
echo "--- Трансляция в белки..."
# gffread берет GTF (разметку генов), накладывает её на НАШ новый геном с мутациями
# и транслирует последовательности в аминокислоты (-y)
gffread "${REF_DIR}/gencode.v46.primary_assembly.annotation.gtf" \
    -g "${PROTEO_DIR}/consensus_genome.fa" \
    -y "${PROTEO_DIR}/personal_proteins_full.faa"

# 4. Финальная фильтрация (оставляем только то, что экспрессируется)
echo "--- Очистка по списку экспрессии..."
# Используем seqtk (он есть в контейнере) или простой grep, если seqtk нет.
# В вашем контейнере seqtk должен быть.
seqtk subseq "${PROTEO_DIR}/personal_proteins_full.faa" "${PROTEO_DIR}/clean_ids.txt" > "${PROTEO_DIR}/personal_proteins.faa"

echo "✅ Готово! Файл создан: ${PROTEO_DIR}/personal_proteins.faa"
EOF
chmod +x pipeline_scripts/05_proteo.sh

