#!/bin/bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$PROJECT_DIR"

source "${SCRIPT_DIR}/config_s3.sh"
source "${SCRIPT_DIR}/common.sh"
source "${SCRIPT_DIR}/pipeline_logging.sh"

SAMPLES_FILE="${SAMPLES_FILE:-samples.txt}"
GENOME_PATH="${REF_DIR}/${REF_GENOME}"
GTF_PATH="${REF_DIR}/${REF_GTF}"
S3_BUCKET="${1:-${S3_BUCKET:-s3://codex-test-ngsdata-calculations/results}}"
AWS_PROFILE="${AWS_PROFILE:-codex-sandbox}"
CLEANUP_AFTER_UPLOAD="${CLEANUP_AFTER_UPLOAD:-true}"
DRY_RUN="${DRY_RUN:-false}"

require_file "$SAMPLES_FILE" "samples list"
require_file "$GENOME_PATH" "reference genome"
require_file "$GTF_PATH" "reference GTF"

run_py() {
    PYTHONPATH="$PROJECT_DIR" python3 -m proteome_pipeline.cli "$@"
}

run_or_echo() {
    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "DRY_RUN: $*"
        return 0
    fi
    run_logged "$*" "$@"
}

ensure_vcf_gz() {
    local gatk_dir="$1"
    if [[ -f "${gatk_dir}/variants_filtered.vcf.gz" ]]; then
        return 0
    fi
    if [[ -f "${gatk_dir}/variants_filtered.vcf" ]]; then
        run_or_echo bcftools view "${gatk_dir}/variants_filtered.vcf" -Oz -o "${gatk_dir}/variants_filtered.vcf.gz"
        run_or_echo bcftools index -f "${gatk_dir}/variants_filtered.vcf.gz"
    fi
}

generate_snp_proteome() {
    local sample_id="$1"
    local work_dir="${BASE_DIR}/${sample_id}"
    local proteo_dir="${work_dir}/results_proteins"
    local vcf_in="${work_dir}/results_gatk/variants_filtered.vcf.gz"
    local ready_vcf="${work_dir}/results_gatk/ready_snps.vcf.gz"

    require_file "$vcf_in" "filtered VCF"
    mkdir -p "$proteo_dir"
    run_or_echo bcftools index -f "$vcf_in"

    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "DRY_RUN: bcftools view -v snps ${vcf_in} | bcftools norm -m -any -f ${GENOME_PATH} -Oz -o ${ready_vcf}"
    else
        log_step "START prepare SNP VCF"
        bcftools view -v snps "$vcf_in" | bcftools norm -m -any -f "$GENOME_PATH" -Oz -o "$ready_vcf" >> "$PIPELINE_LOG_FILE" 2>&1
        log_step "OK prepare SNP VCF"
    fi
    run_or_echo bcftools index -f "$ready_vcf"

    run_or_echo bcftools consensus -H 1 -f "$GENOME_PATH" "$ready_vcf" -o "${proteo_dir}/genome_h1.fa"
    run_or_echo samtools faidx "${proteo_dir}/genome_h1.fa"
    run_or_echo gffread "$GTF_PATH" -g "${proteo_dir}/genome_h1.fa" -y "${proteo_dir}/prot_h1.raw.fa" --ids "${proteo_dir}/clean_ids.txt"
    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "DRY_RUN: prefix H1 FASTA headers"
    else
        sed 's/>/>H1_/' "${proteo_dir}/prot_h1.raw.fa" > "${proteo_dir}/prot_h1.fa"
    fi

    run_or_echo bcftools consensus -H 2 -f "$GENOME_PATH" "$ready_vcf" -o "${proteo_dir}/genome_h2.fa"
    run_or_echo samtools faidx "${proteo_dir}/genome_h2.fa"
    run_or_echo gffread "$GTF_PATH" -g "${proteo_dir}/genome_h2.fa" -y "${proteo_dir}/prot_h2.raw.fa" --ids "${proteo_dir}/clean_ids.txt"
    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "DRY_RUN: prefix H2 FASTA headers"
    else
        sed 's/>/>H2_/' "${proteo_dir}/prot_h2.raw.fa" > "${proteo_dir}/prot_h2.fa"
        cat "${proteo_dir}/prot_h1.fa" "${proteo_dir}/prot_h2.fa" > "${proteo_dir}/protein_raw.fasta"
    fi

    run_or_echo run_py clean-proteins --input "${proteo_dir}/protein_raw.fasta" --out "${proteo_dir}/protein.fasta"
    run_or_echo run_py nonsense-report --input "${proteo_dir}/protein.fasta" --out "${proteo_dir}/nonsense_candidates.txt"
}

upload_file_to_prefix() {
    local sample_id="$1"
    local file="$2"
    local remote_prefix="$3"
    if [[ -f "$file" ]]; then
        run_or_echo aws s3 cp "$file" "${remote_prefix}/$(basename "$file")" --profile "$AWS_PROFILE" --region us-east-1 --quiet
    else
        log_warn "upload skipped, missing file: ${file}"
    fi
}

upload_proteome_if_exists() {
    local sample_id="$1"
    local file="$2"
    upload_file_to_prefix "$sample_id" "$file" "${S3_BUCKET}/${sample_id}/results_proteins"
}

upload_dir_if_exists() {
    local sample_id="$1"
    local dir="$2"
    local remote_prefix="$3"
    if [[ -d "$dir" ]]; then
        run_or_echo aws s3 cp "$dir" "${remote_prefix}/" --recursive --profile "$AWS_PROFILE" --region us-east-1 --quiet
    else
        log_warn "directory upload skipped, missing directory: ${dir}"
    fi
}

upload_analysis_outputs() {
    local sample_id="$1"
    local work_dir="${BASE_DIR}/${sample_id}"

    upload_dir_if_exists "$sample_id" "${work_dir}/results_qc" "${S3_BUCKET}/${sample_id}/results_qc"
    upload_dir_if_exists "$sample_id" "${work_dir}/results_expression" "${S3_BUCKET}/${sample_id}/results_expression"
    upload_dir_if_exists "$sample_id" "${work_dir}/results_star" "${S3_BUCKET}/${sample_id}/results_star"
    upload_dir_if_exists "$sample_id" "${work_dir}/results_gatk" "${S3_BUCKET}/${sample_id}/results_gatk"
    upload_dir_if_exists "$sample_id" "${work_dir}/logs" "${S3_BUCKET}/${sample_id}/logs"
}

upload_results() {
    local sample_id="$1"
    local proteo_dir="$2"

    upload_analysis_outputs "$sample_id"
    upload_proteome_if_exists "$sample_id" "${proteo_dir}/protein.fasta"
    upload_proteome_if_exists "$sample_id" "${proteo_dir}/clean_ids.txt"
    upload_proteome_if_exists "$sample_id" "${proteo_dir}/nonsense_candidates.txt"
    upload_proteome_if_exists "$sample_id" "${proteo_dir}/frameshift_unique.fasta"
    upload_proteome_if_exists "$sample_id" "${proteo_dir}/combine_proteome.fasta"
    upload_proteome_if_exists "$sample_id" "${proteo_dir}/verification_report.tsv"
    upload_proteome_if_exists "$sample_id" "${proteo_dir}/manifest.txt"
    upload_proteome_if_exists "$sample_id" "$PIPELINE_EVENT_FILE"
    upload_proteome_if_exists "$sample_id" "$PIPELINE_LOG_FILE"
    upload_proteome_if_exists "$sample_id" "${proteo_dir}/UPLOAD_SUCCESS"
}

cleanup_after_success() {
    local work_dir="$1"
    local marker="$2"
    if [[ "$CLEANUP_AFTER_UPLOAD" != "true" ]]; then
        log_info "cleanup disabled: CLEANUP_AFTER_UPLOAD=${CLEANUP_AFTER_UPLOAD}"
        return 0
    fi
    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "DRY_RUN: would remove ${work_dir}"
        return 0
    fi
    require_file "$marker" "upload success marker"
    case "$work_dir" in
        "$BASE_DIR"/*) ;;
        *) pipeline_die "Refusing cleanup outside BASE_DIR: ${work_dir}" ;;
    esac
    log_warn "removing local work directory after verified upload: ${work_dir}"
    rm -rf "$work_dir"
}

while IFS= read -r raw_sample_id || [[ -n "$raw_sample_id" ]]; do
    SAMPLE_ID="$(printf '%s' "$raw_sample_id" | tr -d '[:space:]')"
    [[ -z "$SAMPLE_ID" || "$SAMPLE_ID" == \#* ]] && continue
    SAMPLE_ID="$(require_sample_id "$SAMPLE_ID")"

    WORK_DIR="${BASE_DIR}/${SAMPLE_ID}"
    PROTEO_DIR="${WORK_DIR}/results_proteins"
    LOG_DIR="${WORK_DIR}/logs"
    init_pipeline_logging "$SAMPLE_ID" "$LOG_DIR" "local_proteome"

    log_info "sample=${SAMPLE_ID}"
    log_info "output_s3=${S3_BUCKET}"
    log_info "cleanup_after_upload=${CLEANUP_AFTER_UPLOAD}"
    log_info "dry_run=${DRY_RUN}"

    mkdir -p "$PROTEO_DIR"
    run_logged "ensure filtered VCF gzip" ensure_vcf_gz "${WORK_DIR}/results_gatk"
    require_file "${WORK_DIR}/results_expression/abundance.tsv" "abundance table"

    run_logged "clean transcript ids" run_py clean-ids \
        --gtf "$GTF_PATH" \
        --abundance "${WORK_DIR}/results_expression/abundance.tsv" \
        --out "${PROTEO_DIR}/clean_ids.txt" \
        --min-tpm "${MIN_TPM:-1}" \
        --min-cds-bp "${MIN_CDS_BP:-10}"

    run_logged "generate SNP proteome" generate_snp_proteome "$SAMPLE_ID"
    run_logged "generate frameshift proteome" run_py frameshifts \
        --vcf "${WORK_DIR}/results_gatk/variants_filtered.vcf.gz" \
        --gtf "$GTF_PATH" \
        --genome "$GENOME_PATH" \
        --ids "${PROTEO_DIR}/clean_ids.txt" \
        --out "${PROTEO_DIR}/frameshift_unique.fasta" \
        --min-aa "${MIN_FRAMESHIFT_AA:-6}"

    run_logged "combine proteomes" run_py combine \
        --out "${PROTEO_DIR}/combine_proteome.fasta" \
        "${PROTEO_DIR}/protein.fasta" \
        "${PROTEO_DIR}/frameshift_unique.fasta"

    write_pipeline_manifest "${PROTEO_DIR}/manifest.txt" \
        "input_source=local_public_sra_workflow" \
        "output_s3=${S3_BUCKET}" \
        "sample=${SAMPLE_ID}" \
        "base_dir=${BASE_DIR}" \
        "genome=${GENOME_PATH}" \
        "gtf=${GTF_PATH}" \
        "min_tpm=${MIN_TPM:-1}" \
        "min_cds_bp=${MIN_CDS_BP:-10}" \
        "min_frameshift_aa=${MIN_FRAMESHIFT_AA:-6}"

    if run_py verify --proteo-dir "$PROTEO_DIR" --out "${PROTEO_DIR}/verification_report.tsv"; then
        log_info "verification OK"
    else
        log_error "verification requires attention; uploading diagnostics and skipping cleanup"
        upload_results "$SAMPLE_ID" "$PROTEO_DIR"
        continue
    fi

    touch "${PROTEO_DIR}/UPLOAD_SUCCESS"
    upload_results "$SAMPLE_ID" "$PROTEO_DIR"
    cleanup_after_success "$WORK_DIR" "${PROTEO_DIR}/UPLOAD_SUCCESS"
done < "$SAMPLES_FILE"
