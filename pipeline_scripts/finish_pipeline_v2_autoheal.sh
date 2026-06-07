#!/bin/bash
set -e

# ==============================================================================
# 🚀 ЗАПУСК V2.6: Двойной протеом + AUTO-HEAL (АВТО-ИСПРАВЛЕНИЕ ОШИБОК)
# ==============================================================================
# Если gffread падает (SegFault), скрипт сам находит виновника, удаляет его 
# и перезапускается.

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
# 🐍 PYTHON 1: ДЕТЕКТИВ (Поиск транскрипта, убивающего gffread)
# ------------------------------------------------------------------------------
cat << 'EOF' > find_killer_embedded.py
import sys, subprocess, os

gtf_path = sys.argv[1]
genome_path = sys.argv[2]
ids_file = sys.argv[3]
temp_out = "debug_temp.fa"
temp_ids = "debug_ids.txt"

# Читаем все ID
with open(ids_file, 'r') as f:
    all_ids = [line.strip() for line in f if line.strip()]

def test_batch(batch_ids):
    with open(temp_ids, 'w') as f:
        for tid in batch_ids: f.write(f"{tid}\n")
    cmd = ["gffread", gtf_path, "-g", genome_path, "-y", temp_out, "--ids", temp_ids]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return False # Успех
    except subprocess.CalledProcessError:
        return True # Упал

# Бинарный поиск
candidates = all_ids
while len(candidates) > 1:
    mid = len(candidates) // 2
    left = candidates[:mid]
    if test_batch(left): candidates = left
    else: candidates = candidates[mid:]

if candidates:
    killer = candidates[0]
    print(killer) # Выводим только ID для захвата в bash
EOF

# ------------------------------------------------------------------------------
# 🐍 PYTHON 2: ФИЛЬТР КОРОТКИХ ТРАНСКРИПТОВ
# ------------------------------------------------------------------------------
cat << 'EOF' > filter_short_transcripts.py
import sys, re
gtf_file = sys.argv[1]
input_ids = sys.argv[2]
output_ids = sys.argv[3]
MIN_BP = 10 
target_ids = set()
with open(input_ids, 'r') as f:
    for line in f: target_ids.add(line.strip())
cds_lengths = {} 
with open(gtf_file, 'r') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'CDS': continue
        attr = parts[8]
        tid_match = re.search(r'transcript_id "([^"]+)"', attr)
        if not tid_match: continue
        tid = tid_match.group(1)
        if tid not in target_ids: continue
        length = int(parts[4]) - int(parts[3]) + 1
        cds_lengths[tid] = cds_lengths.get(tid, 0) + length
with open(output_ids, 'w') as out:
    for tid in target_ids:
        if cds_lengths.get(tid, 0) >= MIN_BP: out.write(f"{tid}\n")
EOF

# ------------------------------------------------------------------------------
# 🐍 PYTHON 3: ОЧИСТКА ПОСЛЕДОВАТЕЛЬНОСТЕЙ
# ------------------------------------------------------------------------------
cat << 'EOF' > clean_proteins_embedded.py
import sys
infile = sys.argv[1]
outfile = sys.argv[2]
with open(infile, 'r') as fin, open(outfile, 'w') as fout:
    header = None
    seq_parts = []
    def process_seq(h, s):
        if not h: return
        if s.count('.') > 1: return 
        if '.' in s: s = s.split('.')[0]
        if len(s) == 0: return
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
# 🐍 PYTHON 4: АНАЛИЗ NONSENSE
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
# ФУНКЦИЯ: БЕЗОПАСНЫЙ ЗАПУСК GFFREAD (С АВТО-ЛЕЧЕНИЕМ)
# ==============================================================================
run_safe_gffread() {
    local G_GTF="$1"
    local G_GENOME="$2"
    local G_OUT="$3"
    local G_IDS="$4"
    local HAPLO="$5"

    local MAX_RETRIES=5
    local ATTEMPT=1

    while [ $ATTEMPT -le $MAX_RETRIES ]; do
        echo "   ▶️  gffread $HAPLO (Попытка $ATTEMPT/$MAX_RETRIES)..."
        
        # Запуск gffread. Если успешно (код 0) - выходим из цикла
        if gffread "$G_GTF" -g "$G_GENOME" -y "$G_OUT" --ids "$G_IDS" 2>/dev/null; then
            return 0
        fi

        echo "   💥 ОШИБКА: gffread упал! Запускаем детектива..."
        
        # Запускаем Python-скрипт для поиска "убийцы"
        # Перенаправляем вывод скрипта в переменную KILLER_ID
        KILLER_ID=$(python3 find_killer_embedded.py "$G_GTF" "$G_GENOME" "$G_IDS")
        
        if [ -n "$KILLER_ID" ]; then
            echo "   🚑 НАЙДЕН ВИНОВНИК: $KILLER_ID. Удаляем из списка..."
            # Удаляем ID из файла clean_ids.txt
            sed -i "/$KILLER_ID/d" "$G_IDS"
        else
            echo "   ❌ Не удалось найти виновника автоматически."
            return 1
        fi
        
        ATTEMPT=$((ATTEMPT + 1))
    done

    echo "   ❌ Превышено количество попыток восстановления."
    return 1
}

# ==============================================================================
# ОСНОВНОЙ ЦИКЛ
# ==============================================================================

while read SAMPLE_ID; do
    [ -z "$SAMPLE_ID" ] && continue
    echo "--------------------------------------------------"
    echo "🧬 ОБРАБОТКА: $SAMPLE_ID"
    
    WORK_DIR="${BASE_DIR:-/home/ubuntu}/${SAMPLE_ID}"
    # НЕ удаляем всю папку, если это повторный запуск, чтобы не качать заново
    # rm -rf "$WORK_DIR" 
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

    # 3. ФИЛЬТР ID
    echo "🛡  Создание списка ID..."
    PROTEO_DIR="$WORK_DIR/results_proteins"
    RAW_IDS="$WORK_DIR/raw_ids.tmp"
    CLEAN_IDS="$WORK_DIR/results_proteins/clean_ids.txt"
    mkdir -p "$PROTEO_DIR"

    grep -vE "^chrM|^M|^MT" "$GTF_PATH" | grep -o 'transcript_id "[^"]*"' | cut -d'"' -f2 | sort -u > "$WORK_DIR/nuclear_ids.tmp"
    
    if [ -f "$ABUNDANCE_FILE" ]; then
        awk -F "\t" 'NR>1 && $5 > 1 {split($1, a, "|"); print a[1]}' "$ABUNDANCE_FILE" | sort -u > "$WORK_DIR/expressed_ids.tmp"
        comm -12 "$WORK_DIR/nuclear_ids.tmp" "$WORK_DIR/expressed_ids.tmp" > "$RAW_IDS"
    else
        cat "$WORK_DIR/nuclear_ids.tmp" > "$RAW_IDS"
    fi

    # Фильтр длины
    echo "📏 Фильтрация коротких транскриптов..."
    python3 filter_short_transcripts.py "$GTF_PATH" "$RAW_IDS" "$CLEAN_IDS"

    # Старые известные "убийцы" (для скорости сразу удаляем)
    sed -i '/ENST00000467006/d' "$CLEAN_IDS"
    sed -i '/ENST00000262126/d' "$CLEAN_IDS"

    # 4. ТРАНСЛЯЦИЯ С АВТО-ЛЕЧЕНИЕМ
    echo "🔬 Генерация протеомов (с защитой от сбоев)..."
    
    # H1
    bcftools consensus -H 1 -f "$GENOME_PATH" "$VCF_READY" -o "$PROTEO_DIR/genome_h1.fa" 2>/dev/null
    samtools faidx "$PROTEO_DIR/genome_h1.fa"
    
    # ЗАПУСК ФУНКЦИИ БЕЗОПАСНОГО GFFREAD
    run_safe_gffread "$GTF_PATH" "$PROTEO_DIR/genome_h1.fa" "$PROTEO_DIR/prot_h1.fa" "$CLEAN_IDS" "H1"
    
    sed -i 's/>/>H1_/' "$PROTEO_DIR/prot_h1.fa"
    
    # H2
    bcftools consensus -H 2 -f "$GENOME_PATH" "$VCF_READY" -o "$PROTEO_DIR/genome_h2.fa" 2>/dev/null
    samtools faidx "$PROTEO_DIR/genome_h2.fa"
    
    # ЗАПУСК ФУНКЦИИ БЕЗОПАСНОГО GFFREAD
    run_safe_gffread "$GTF_PATH" "$PROTEO_DIR/genome_h2.fa" "$PROTEO_DIR/prot_h2.fa" "$CLEAN_IDS" "H2"
    
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
rm -f analyze_nonsense_embedded.py clean_proteins_embedded.py filter_short_transcripts.py find_killer_embedded.py debug_temp.fa debug_ids.txt
echo "🎉 ПАЙПЛАЙН V2.6 (AUTO-HEAL) ЗАВЕРШЕН!"
