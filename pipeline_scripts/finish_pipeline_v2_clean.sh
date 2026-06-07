#!/bin/bash
set -e

# ==============================================================================
# 🚀 ЗАПУСК V2.4: Двойной протеом + СИСТЕМНЫЙ ФИЛЬТР КОРОТКИХ БЕЛКОВ
# ==============================================================================

# 1. ЗАГРУЗКА КОНФИГА
CONFIG_FILE="pipeline_scripts/config_s3.sh"
if [ -f "$CONFIG_FILE" ]; then
    source "$CONFIG_FILE"
else
    echo "❌ Ошибка: Конфиг не найден ($CONFIG_FILE)"
    exit 1
fi

SAMPLES_FILE="samples.txt"
GENOME_PATH="${REF_DIR}/${REF_GENOME}"
GTF_PATH="${REF_DIR}/${REF_GTF}"

# Проверки
if [ -z "$S3_BUCKET" ]; then echo "❌ S3_BUCKET не задан."; exit 1; fi
if [ ! -f "$GENOME_PATH" ]; then echo "❌ Нет генома."; exit 1; fi
if [ ! -f "$GTF_PATH" ]; then echo "❌ Нет GTF."; exit 1; fi

# ------------------------------------------------------------------------------
# 🐍 PYTHON 1: ФИЛЬТР КОРОТКИХ ТРАНСКРИПТОВ (< 3 кодонов)
# ------------------------------------------------------------------------------
cat << 'EOF' > filter_short_transcripts.py
import sys
import re

# Аргументы: 1 = GTF файл, 2 = Входной список ID, 3 = Выходной список ID
gtf_file = sys.argv[1]
input_ids = sys.argv[2]
output_ids = sys.argv[3]
MIN_BP = 10 # Минимальная длина (3 кодона * 3 = 9. Мы требуем > 9, то есть 10+)

print(f"   📏 Фильтрация: оставляем только CDS длиннее {MIN_BP} bp...")

# 1. Читаем список желаемых ID
target_ids = set()
with open(input_ids, 'r') as f:
    for line in f:
        target_ids.add(line.strip())

# 2. Считаем длину CDS для каждого транскрипта в GTF
cds_lengths = {} # {transcript_id: total_bp}

with open(gtf_file, 'r') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        
        # Нас интересуют только CDS (кодирующие регионы)
        if len(parts) < 9 or parts[2] != 'CDS': 
            continue
            
        # Парсим ID
        attr = parts[8]
        tid_match = re.search(r'transcript_id "([^"]+)"', attr)
        if not tid_match: continue
        tid = tid_match.group(1)
        
        # Если этого ID нет в нашем списке - пропускаем, экономим время
        if tid not in target_ids:
            continue
            
        length = int(parts[4]) - int(parts[3]) + 1
        cds_lengths[tid] = cds_lengths.get(tid, 0) + length

# 3. Фильтруем и записываем
kept_count = 0
removed_count = 0

with open(output_ids, 'w') as out:
    for tid in target_ids:
        # Если в GTF вообще не нашли CDS для этого ID (например, это ncRNA),
        # то длина будет 0, и он удалится. Это правильно для протеомики.
        length = cds_lengths.get(tid, 0)
        
        if length >= MIN_BP:
            out.write(f"{tid}\n")
            kept_count += 1
        else:
            removed_count += 1
            # Можно раскомментировать для отладки:
            # print(f"Removing {tid}: length {length} bp")

print(f"      ✅ Оставлено: {kept_count}")
print(f"      🗑  Удалено (слишком короткие): {removed_count}")

EOF

# ------------------------------------------------------------------------------
# 🐍 PYTHON 2: ОЧИСТКА ПОСЛЕДОВАТЕЛЬНОСТЕЙ
# ------------------------------------------------------------------------------
cat << 'EOF' > clean_proteins_embedded.py
import sys
infile = sys.argv[1]
outfile = sys.argv[2]
kept_count = 0
removed_count = 0
with open(infile, 'r') as fin, open(outfile, 'w') as fout:
    header = None
    seq_parts = []
    def process_seq(h, s):
        global kept_count, removed_count
        if not h: return
        dots = s.count('.')
        if dots > 1:
            removed_count += 1
            return 
        if '.' in s: s = s.split('.')[0]
        if len(s) == 0: return
        kept_count += 1
        fout.write(f"{h}\n{s}\n")
    for line in fin:
        line = line.strip()
        if line.startswith('>'):
            if header: process_seq(header, "".join(seq_parts))
            header = line
            seq_parts = []
        else: seq_parts.append(line)
    if header: process_seq(header, "".join(seq_parts))
EOF

# ------------------------------------------------------------------------------
# 🐍 PYTHON 3: АНАЛИЗ NONSENSE
# ------------------------------------------------------------------------------
cat << 'EOF' > analyze_nonsense_embedded.py
import sys
protein_file = sys.argv[1]
output_file = sys.argv[2]
data = {} 
count = 0
with open(protein_file, 'r') as f, open(output_file, 'w') as out:
    out.write("TRANSCRIPT_ID\tLEN_H1\tLEN_H2\tLOSS\tSTATUS\n")
    header = None
    seq_len = 0
    def save(h, length):
        clean = h.strip().split()[0][1:] 
        if '_' not in clean: return
        haplo = clean.split('_')[0]
        tid = "_".join(clean.split('_')[1:])
        if tid not in data: data[tid] = {}
        data[tid][haplo] = length
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if header: save(header, seq_len)
            header = line
            seq_len = 0
        else: seq_len += len(line)
    if header: save(header, seq_len)
    for tid, lens in data.items():
        if 'H1' not in lens or 'H2' not in lens: continue
        l1, l2 = lens['H1'], lens['H2']
        if l1 == l2: continue
        diff = l1 - l2
        maxlen = max(l1, l2)
        if maxlen == 0: continue
        if abs(diff) > 10 and (abs(diff)/maxlen > 0.05):
            status = "H2_TRUNCATED" if l1 > l2 else "H1_TRUNCATED"
            out.write(f"{tid}\t{l1}\t{l2}\t-{abs(diff)}\t{status}\n")
            count += 1   
EOF

# ==============================================================================
# ОСНОВНОЙ ЦИКЛ
# ==============================================================================

while read SAMPLE_ID; do
    [ -z "$SAMPLE_ID" ] && continue
    echo "--------------------------------------------------"
    echo "🧬 ОБРАБОТКА: $SAMPLE_ID"
    
    WORK_DIR="${BASE_DIR:-/home/ubuntu}/${SAMPLE_ID}"
    rm -rf "$WORK_DIR"
    mkdir -p "$WORK_DIR/results_gatk" "$WORK_DIR/results_proteins" "$WORK_DIR/results_expression"
    
    # 1. СКАЧИВАНИЕ
    echo "⬇️  Скачивание..."
    aws s3 cp "$S3_BUCKET/${SAMPLE_ID}/results_gatk/variants_filtered.vcf.gz" "$WORK_DIR/results_gatk/" --quiet || true
    if [ ! -f "$WORK_DIR/results_gatk/variants_filtered.vcf.gz" ]; then
         aws s3 cp "$S3_BUCKET/${SAMPLE_ID}/results_gatk/variants_filtered.vcf" "$WORK_DIR/results_gatk/" --quiet || true
         [ -f "$WORK_DIR/results_gatk/variants_filtered.vcf" ] && bcftools view "$WORK_DIR/results_gatk/variants_filtered.vcf" -Oz -o "$WORK_DIR/results_gatk/variants_filtered.vcf.gz"
    fi
    aws s3 cp "$S3_BUCKET/${SAMPLE_ID}/results_expression/abundance.tsv" "$WORK_DIR/results_expression/" --quiet || true
    
    VCF_IN="$WORK_DIR/results_gatk/variants_filtered.vcf.gz"
    ABUNDANCE_FILE="$WORK_DIR/results_expression/abundance.tsv"
    
    if [ ! -f "$VCF_IN" ]; then echo "⚠️ Нет VCF. Пропуск."; continue; fi
    bcftools index -f "$VCF_IN"

    # 2. ПОДГОТОВКА VCF (Только SNP)
    echo "📦 Подготовка VCF..."
    VCF_READY="$WORK_DIR/results_gatk/ready.vcf.gz"
    bcftools view -v snps "$VCF_IN" | bcftools norm -m -any -f "$GENOME_PATH" -Oz -o "$VCF_READY"
    bcftools index "$VCF_READY"

    # 3. ФИЛЬТР ID (БЕЗ МИТОХОНДРИЙ + ЭКСПРЕССИЯ)
    echo "🛡  Первичный отбор ID..."
    mkdir -p "$WORK_DIR/results_proteins"
    RAW_IDS="$WORK_DIR/raw_ids.tmp"
    CLEAN_IDS="$WORK_DIR/results_proteins/clean_ids.txt"

    grep -vE "^chrM|^M|^MT" "$GTF_PATH" | grep -o 'transcript_id "[^"]*"' | cut -d'"' -f2 | sort -u > "$WORK_DIR/nuclear_ids.tmp"
    
    if [ -f "$ABUNDANCE_FILE" ]; then
        awk -F "\t" 'NR>1 && $5 > 1 {split($1, a, "|"); print a[1]}' "$ABUNDANCE_FILE" | sort -u > "$WORK_DIR/expressed_ids.tmp"
        comm -12 "$WORK_DIR/nuclear_ids.tmp" "$WORK_DIR/expressed_ids.tmp" > "$RAW_IDS"
    else
        cat "$WORK_DIR/nuclear_ids.tmp" > "$RAW_IDS"
    fi

    # === [НОВЫЙ ЭТАП] ФИЛЬТР ПО ДЛИНЕ (СИСТЕМНЫЙ) ===
    # Удаляем всё, что короче 3 кодонов (>9 bp)
    echo "📏 Фильтрация коротких транскриптов..."
    python3 filter_short_transcripts.py "$GTF_PATH" "$RAW_IDS" "$CLEAN_IDS"
    # ================================================
    
    PROTEO_DIR="$WORK_DIR/results_proteins"
    
    # 4. ТРАНСЛЯЦИЯ
    echo "🔬 Генерация протеомов..."
    
    # H1
    bcftools consensus -H 1 -f "$GENOME_PATH" "$VCF_READY" -o "$PROTEO_DIR/genome_h1.fa" 2>/dev/null
    samtools faidx "$PROTEO_DIR/genome_h1.fa"
    gffread "$GTF_PATH" -g "$PROTEO_DIR/genome_h1.fa" -y "$PROTEO_DIR/prot_h1.fa" --ids "$CLEAN_IDS"
    sed -i 's/>/>H1_/' "$PROTEO_DIR/prot_h1.fa"
    
    # H2
    bcftools consensus -H 2 -f "$GENOME_PATH" "$VCF_READY" -o "$PROTEO_DIR/genome_h2.fa" 2>/dev/null
    samtools faidx "$PROTEO_DIR/genome_h2.fa"
    gffread "$GTF_PATH" -g "$PROTEO_DIR/genome_h2.fa" -y "$PROTEO_DIR/prot_h2.fa" --ids "$CLEAN_IDS"
    sed -i 's/>/>H2_/' "$PROTEO_DIR/prot_h2.fa"
    
    # Объединение
    cat "$PROTEO_DIR/prot_h1.fa" "$PROTEO_DIR/prot_h2.fa" > "$PROTEO_DIR/protein_raw.fasta"
    
    # 4.5 ЧИСТКА
    python3 clean_proteins_embedded.py "$PROTEO_DIR/protein_raw.fasta" "$PROTEO_DIR/protein.fasta"
    
    # 5. АНАЛИЗ
    echo "🕵️  Поиск Nonsense-мутаций..."
    REPORT_FILE="$PROTEO_DIR/nonsense_candidates.txt"
    python3 analyze_nonsense_embedded.py "$PROTEO_DIR/protein.fasta" "$REPORT_FILE"
    
    # 6. ВЫГРУЗКА
    echo "☁️  Загрузка в S3..."
    aws s3 cp "$PROTEO_DIR/protein.fasta" "$S3_BUCKET/${SAMPLE_ID}/results_proteins/protein.fasta" --quiet
    aws s3 cp "$PROTEO_DIR/clean_ids.txt" "$S3_BUCKET/${SAMPLE_ID}/results_proteins/clean_ids.txt" --quiet
    aws s3 cp "$REPORT_FILE" "$S3_BUCKET/${SAMPLE_ID}/results_proteins/nonsense_candidates.txt" --quiet
    
    rm -rf "$WORK_DIR"
    echo "✨ Готово: $SAMPLE_ID"

done < "$SAMPLES_FILE"

# Удаляем временные скрипты
rm -f analyze_nonsense_embedded.py clean_proteins_embedded.py filter_short_transcripts.py
echo "🎉 ПАЙПЛАЙН V2.4 ЗАВЕРШЕН!"
