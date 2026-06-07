#!/bin/bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/common.sh"

SAMPLE_ID="$(require_sample_id "${1:-}")"
SAMPLE_DIR="${BASE_DIR}/${SAMPLE_ID}"
REF_GENOME_PATH="${REF_DIR}/${REF_GENOME}"
REF_GTF_PATH="${REF_DIR}/${REF_GTF}"
PROTEO_DIR="${SAMPLE_DIR}/results_proteins"
EXPRESSION_FILE="${SAMPLE_DIR}/results_expression/abundance.tsv"
VCF_GZ="${SAMPLE_DIR}/results_gatk/variants_filtered.vcf.gz"
VCF_PLAIN="${SAMPLE_DIR}/results_gatk/variants_filtered.vcf"

require_file "$REF_GENOME_PATH" "reference genome"
require_file "$REF_GTF_PATH" "reference GTF"
require_file "$EXPRESSION_FILE" "expression table"

mkdir -p "$PROTEO_DIR"

echo ">>> [STEP 5] PROTEO: ${SAMPLE_ID}"

echo "--- Filtering expressed transcripts..."
awk -F "\t" 'NR>1 && $5 > 1 {split($1, a, "|"); print a[1]}' "$EXPRESSION_FILE" | sort -u > "${PROTEO_DIR}/clean_ids.txt"
require_file "${PROTEO_DIR}/clean_ids.txt" "clean transcript id list"

if [[ ! -f "${REF_GENOME_PATH}.fai" ]]; then
    samtools faidx "$REF_GENOME_PATH"
fi

if [[ -f "$VCF_GZ" ]]; then
    VCF_FILE="$VCF_GZ"
elif [[ -f "$VCF_PLAIN" ]]; then
    VCF_FILE="$VCF_PLAIN"
else
    pipeline_die "Missing filtered VCF for ${SAMPLE_ID}"
fi

echo "--- Generating consensus genome..."
bcftools consensus -f "$REF_GENOME_PATH" "$VCF_FILE" > "${PROTEO_DIR}/consensus_genome.fa"

echo "--- Translating proteins..."
gffread "$REF_GTF_PATH" \
    -g "${PROTEO_DIR}/consensus_genome.fa" \
    -y "${PROTEO_DIR}/personal_proteins_full.faa"

echo "--- Selecting expressed proteins..."
python3 - "${PROTEO_DIR}/personal_proteins_full.faa" "${PROTEO_DIR}/clean_ids.txt" "${PROTEO_DIR}/personal_proteins.faa" <<'PY'
import sys
from pathlib import Path

fasta_path = Path(sys.argv[1])
ids_path = Path(sys.argv[2])
out_path = Path(sys.argv[3])
wanted = {line.strip().split('|')[0] for line in ids_path.read_text().splitlines() if line.strip()}

written = 0
with fasta_path.open() as src, out_path.open("w") as out:
    header = None
    seq = []

    def flush():
        global written
        if not header:
            return
        transcript_id = header[1:].split()[0].split('|')[0]
        if transcript_id in wanted:
            out.write(header + "\n")
            out.write("".join(seq) + "\n")
            written += 1

    for raw in src:
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header = line
            seq = []
        else:
            seq.append(line)
    flush()

if written == 0:
    raise SystemExit("No expressed proteins were selected")
PY

echo "Protein FASTA created: ${PROTEO_DIR}/personal_proteins.faa"
