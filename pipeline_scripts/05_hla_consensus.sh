#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${SCRIPT_DIR}/common.sh"

export PATH="/home/codex/tools/bin:${PATH}"

AWS_PROFILE="${AWS_PROFILE:-codex-sandbox}"
AWS_REGION="${AWS_REGION:-us-east-1}"
SOURCE_S3="${SOURCE_S3:-s3://bioinfo-data-amylice-2026}"
DEST_S3="${DEST_S3:-s3://codex-test-ngsdata-calculations/prepared-bioinfo-proteome}"
WORK_ROOT="${WORK_ROOT:-/home/codex/hla_runs}"
DRY_RUN="${DRY_RUN:-false}"
RUN_ARCAS="${RUN_ARCAS:-true}"
RUN_OPTITYPE="${RUN_OPTITYPE:-true}"
RUN_HISAT="${RUN_HISAT:-true}"
HLA_THREADS="${HLA_THREADS:-4}"
HLA_SINGLE_END="${HLA_SINGLE_END:-true}"

ARCASHLA_CMD="${ARCASHLA_CMD:-arcasHLA}"
OPTITYPE_CMD="${OPTITYPE_CMD:-OptiTypePipeline.py}"
HISAT_CMD="${HISAT_CMD:-hisatgenotype}"
HISAT_INDEX_DIR="${HISAT_INDEX_DIR:-/home/codex/tools/src/hisat-genotype/indicies}"
HISAT_GENOTYPE_GENOME="${HISAT_GENOTYPE_GENOME:-genotype_genome}"
HLA_PY="${HLA_PY:-${PROJECT_DIR}/proteome_pipeline/hla_consensus.py}"
HLA_EXTRACTED_FASTQ_MANIFEST="${HLA_EXTRACTED_FASTQ_MANIFEST:-hla_extracted_fastq.txt}"

usage() {
    cat <<EOF
Usage:
  SAMPLES="SRR32060234 SRR32060239" ./pipeline_scripts/05_hla_consensus.sh
  ./pipeline_scripts/05_hla_consensus.sh SRR32060234 SRR32060239

Environment:
  SOURCE_S3      Input bucket with STAR/GATK BAMs. Default: ${SOURCE_S3}
  DEST_S3        Output bucket/prefix. Default: ${DEST_S3}
  DRY_RUN        true prints planned commands and validates inputs.
  RUN_ARCAS      true/false for arcasHLA.
  RUN_OPTITYPE   true/false for OptiType.
  RUN_HISAT      true/false for HISAT-genotype.
  HLA_SINGLE_END  true adds arcasHLA --single for single-end RNA-seq BAMs.
EOF
}

log() {
    printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"
}

run_or_print() {
    if [[ "${DRY_RUN}" == "true" ]]; then
        printf 'DRY_RUN:'
        printf ' %q' "$@"
        printf '\n'
    else
        "$@"
    fi
}

tool_available() {
    command -v "$1" >/dev/null 2>&1
}

aws_ls_key() {
    local uri="$1"
    aws s3 ls "${uri}" --profile "${AWS_PROFILE}" --region "${AWS_REGION}" 2>/dev/null | awk '{print $4}' | head -n 1
}

resolve_bam_uri() {
    local sample_id="$1"
    local candidates=(
        "${SOURCE_S3}/${sample_id}/results_star/dedup.bam"
        "${SOURCE_S3}/${sample_id}/results_star/Aligned.sortedByCoord.out.bam"
        "${SOURCE_S3}/${sample_id}/results_gatk/split.bam"
    )
    local uri
    for uri in "${candidates[@]}"; do
        if [[ -n "$(aws_ls_key "${uri}")" ]]; then
            printf '%s\n' "${uri}"
            return 0
        fi
    done
    return 1
}

write_status() {
    local path="$1"
    local sample_id="$2"
    local caller="$3"
    local status="$4"
    local detail="$5"
    printf '%s\t%s\t%s\t%s\n' "${sample_id}" "${caller}" "${status}" "${detail}" >> "${path}"
}

upload_results() {
    local out_dir="$1"
    local dest_uri="$2"
    run_or_print aws s3 cp "${out_dir}/" "${dest_uri}" \
        --recursive \
        --exclude 'bam/*' \
        --exclude 'hisat_bam_smoke/*' \
        --exclude 'hisat_debug_*/*' \
        --profile "${AWS_PROFILE}" \
        --region "${AWS_REGION}" \
        --no-progress
}

resolve_arcas_fastqs() {
    local sample_id="$1"
    local out_dir="$2"
    [[ -d "${out_dir}/arcas/extract" ]] || return 0
    if [[ "${HLA_SINGLE_END}" == "true" ]]; then
        find "${out_dir}/arcas/extract" -type f \( -name "${sample_id}.extracted.fq.gz" -o -name "${sample_id}.extracted.fastq.gz" \) | sort
    else
        find "${out_dir}/arcas/extract" -type f \( -name "${sample_id}.extracted.1.fq.gz" -o -name "${sample_id}.extracted.2.fq.gz" -o -name "${sample_id}.extracted.1.fastq.gz" -o -name "${sample_id}.extracted.2.fastq.gz" \) | sort
    fi
}

resolve_arcas_genotype_json() {
    local sample_id="$1"
    local out_dir="$2"
    local preferred="${out_dir}/arcas/genotype/${sample_id}.genotype.json"
    if [[ -f "${preferred}" ]]; then
        printf '%s\n' "${preferred}"
        return 0
    fi
    [[ -d "${out_dir}/arcas/genotype" ]] || return 0
    find "${out_dir}/arcas/genotype" -type f -name '*.genotype.json' | sort | head -n 1
}

hisat_index_ready() {
    [[ -s "${HISAT_INDEX_DIR}/${HISAT_GENOTYPE_GENOME}.fa" ]] || return 1
    find "${HISAT_INDEX_DIR}" -maxdepth 1 -type f \( -name "${HISAT_GENOTYPE_GENOME}*.ht2" -o -name "${HISAT_GENOTYPE_GENOME}*.ht2l" \) | grep -q .
}

tsv_data_rows() {
    local path="$1"
    [[ -s "${path}" ]] || {
        printf '0\n'
        return 0
    }
    awk 'NR > 1 && $0 !~ /^[[:space:]]*$/ { count++ } END { print count + 0 }' "${path}"
}

has_available_hla_caller() {
    [[ "${RUN_ARCAS}" == "true" && -n "$(command -v "${ARCASHLA_CMD}" 2>/dev/null || true)" ]] && return 0
    [[ "${RUN_OPTITYPE}" == "true" && -n "$(command -v "${OPTITYPE_CMD}" 2>/dev/null || true)" ]] && return 0
    [[ "${RUN_HISAT}" == "true" && -n "$(command -v "${HISAT_CMD}" 2>/dev/null || true)" ]] && return 0
    return 1
}

run_arcas() {
    local sample_id="$1"
    local bam="$2"
    local out_dir="$3"
    local status_file="$4"
    if [[ "${RUN_ARCAS}" != "true" ]]; then
        write_status "${status_file}" "${sample_id}" "arcasHLA" "SKIPPED" "disabled"
        return 0
    fi
    if ! tool_available "${ARCASHLA_CMD}"; then
        write_status "${status_file}" "${sample_id}" "arcasHLA" "MISSING_TOOL" "${ARCASHLA_CMD}"
        return 0
    fi
    if [[ "${DRY_RUN}" == "true" ]]; then
        local arcas_single_args=()
        if [[ "${HLA_SINGLE_END}" == "true" ]]; then
            arcas_single_args=(--single)
        fi
        run_or_print "${ARCASHLA_CMD}" extract "${bam}" "${arcas_single_args[@]}" -o "${out_dir}/arcas/extract" -t "${HLA_THREADS}"
        run_or_print "${ARCASHLA_CMD}" genotype "${out_dir}/arcas/extract/${sample_id}.extracted.fq.gz" "${arcas_single_args[@]}" -o "${out_dir}/arcas/genotype" -t "${HLA_THREADS}"
        printf '%s\n' "${out_dir}/arcas/extract/${sample_id}.extracted.fq.gz" > "${out_dir}/${HLA_EXTRACTED_FASTQ_MANIFEST}"
        write_status "${status_file}" "${sample_id}" "arcasHLA" "DRY_RUN" "${ARCASHLA_CMD}"
        return 0
    fi
    mkdir -p "${out_dir}/arcas/extract" "${out_dir}/arcas/genotype"
    local arcas_single_args=()
    if [[ "${HLA_SINGLE_END}" == "true" ]]; then
        arcas_single_args=(--single)
    fi
    if ! "${ARCASHLA_CMD}" extract "${bam}" "${arcas_single_args[@]}" -o "${out_dir}/arcas/extract" -t "${HLA_THREADS}" > "${out_dir}/logs/arcas_extract.log" 2>&1; then
        write_status "${status_file}" "${sample_id}" "arcasHLA" "FAIL" "extract failed; see logs/arcas_extract.log"
        return 0
    fi
    local extracted
    if [[ "${HLA_SINGLE_END}" == "true" ]]; then
        mapfile -t extracted < <(find "${out_dir}/arcas/extract" -type f \( -name "${sample_id}.extracted.fq.gz" -o -name "${sample_id}.extracted.fastq.gz" -o -name '*.alignment.p' \) | sort)
    else
        mapfile -t extracted < <(find "${out_dir}/arcas/extract" -type f \( -name "${sample_id}.extracted.1.fq.gz" -o -name "${sample_id}.extracted.2.fq.gz" -o -name "${sample_id}.extracted.1.fastq.gz" -o -name "${sample_id}.extracted.2.fastq.gz" -o -name '*.alignment.p' \) | sort)
    fi
    if [[ ${#extracted[@]} -eq 0 ]]; then
        write_status "${status_file}" "${sample_id}" "arcasHLA" "NO_OUTPUT" "extract produced no FASTQ/alignment.p"
        return 0
    fi
    if ! "${ARCASHLA_CMD}" genotype "${extracted[@]}" "${arcas_single_args[@]}" -o "${out_dir}/arcas/genotype" -t "${HLA_THREADS}" > "${out_dir}/logs/arcas_genotype.log" 2>&1; then
        write_status "${status_file}" "${sample_id}" "arcasHLA" "FAIL" "genotype failed; see logs/arcas_genotype.log"
        return 0
    fi
    local json
    local extracted_fastqs=()
    json="$(resolve_arcas_genotype_json "${sample_id}" "${out_dir}" || true)"
    mapfile -t extracted_fastqs < <(resolve_arcas_fastqs "${sample_id}" "${out_dir}")
    if [[ ${#extracted_fastqs[@]} -gt 0 ]]; then
        printf '%s\n' "${extracted_fastqs[@]}" > "${out_dir}/${HLA_EXTRACTED_FASTQ_MANIFEST}"
    fi
    if [[ -n "${json}" ]]; then
        python3 "${HLA_PY}" normalize --caller arcas --sample "${sample_id}" --input "${json}" --output "${out_dir}/normalized/arcasHLA.tsv"
        write_status "${status_file}" "${sample_id}" "arcasHLA" "OK" "${json}"
    else
        write_status "${status_file}" "${sample_id}" "arcasHLA" "NO_OUTPUT" "no json output found"
    fi
}

run_optitype() {
    local sample_id="$1"
    local bam="$2"
    local out_dir="$3"
    local status_file="$4"
    if [[ "${RUN_OPTITYPE}" != "true" ]]; then
        write_status "${status_file}" "${sample_id}" "OptiType" "SKIPPED" "disabled"
        return 0
    fi
    if ! tool_available "${OPTITYPE_CMD}"; then
        write_status "${status_file}" "${sample_id}" "OptiType" "MISSING_TOOL" "${OPTITYPE_CMD}"
        return 0
    fi
    local optitype_inputs=()
    if [[ -f "${out_dir}/${HLA_EXTRACTED_FASTQ_MANIFEST}" ]]; then
        mapfile -t optitype_inputs < "${out_dir}/${HLA_EXTRACTED_FASTQ_MANIFEST}"
    fi
    if [[ ${#optitype_inputs[@]} -eq 0 ]]; then
        mapfile -t optitype_inputs < <(resolve_arcas_fastqs "${sample_id}" "${out_dir}")
    fi
    if [[ ${#optitype_inputs[@]} -eq 0 ]]; then
        write_status "${status_file}" "${sample_id}" "OptiType" "MISSING_INPUT" "No arcas extracted FASTQ found; not using genome BAM"
        return 0
    fi
    if [[ "${DRY_RUN}" == "true" ]]; then
        run_or_print "${OPTITYPE_CMD}" -i "${optitype_inputs[@]}" --rna -o "${out_dir}/optitype" -p "${sample_id}"
        write_status "${status_file}" "${sample_id}" "OptiType" "DRY_RUN" "${OPTITYPE_CMD}"
        return 0
    fi
    mkdir -p "${out_dir}/optitype"
    if ! "${OPTITYPE_CMD}" -i "${optitype_inputs[@]}" --rna -o "${out_dir}/optitype" -p "${sample_id}" > "${out_dir}/logs/optitype.log" 2>&1; then
        write_status "${status_file}" "${sample_id}" "OptiType" "FAIL" "OptiType failed; see logs/optitype.log"
        return 0
    fi
    local result
    result="$(find "${out_dir}/optitype" -type f \( -name '*result.tsv' -o -name '*.tsv' \) | head -n 1 || true)"
    if [[ -n "${result}" ]]; then
        python3 "${HLA_PY}" normalize --caller optitype --sample "${sample_id}" --input "${result}" --output "${out_dir}/normalized/OptiType.tsv"
        write_status "${status_file}" "${sample_id}" "OptiType" "OK" "${result}"
    else
        write_status "${status_file}" "${sample_id}" "OptiType" "NO_OUTPUT" "no tsv output found"
    fi
}

run_hisat() {
    local sample_id="$1"
    local bam="$2"
    local out_dir="$3"
    local status_file="$4"
    if [[ "${RUN_HISAT}" != "true" ]]; then
        write_status "${status_file}" "${sample_id}" "HISAT-genotype" "SKIPPED" "disabled"
        return 0
    fi
    if ! tool_available "${HISAT_CMD}"; then
        write_status "${status_file}" "${sample_id}" "HISAT-genotype" "MISSING_TOOL" "${HISAT_CMD}"
        return 0
    fi
    if ! hisat_index_ready; then
        write_status "${status_file}" "${sample_id}" "HISAT-genotype" "MISSING_INDEX" "Missing ${HISAT_GENOTYPE_GENOME}.fa or ${HISAT_GENOTYPE_GENOME}*.ht2 in ${HISAT_INDEX_DIR}"
        return 0
    fi
    local hisat_inputs=()
    if [[ -f "${out_dir}/${HLA_EXTRACTED_FASTQ_MANIFEST}" ]]; then
        mapfile -t hisat_inputs < "${out_dir}/${HLA_EXTRACTED_FASTQ_MANIFEST}"
    fi
    if [[ ${#hisat_inputs[@]} -eq 0 ]]; then
        mapfile -t hisat_inputs < <(resolve_arcas_fastqs "${sample_id}" "${out_dir}")
    fi
    if [[ ${#hisat_inputs[@]} -eq 0 ]]; then
        write_status "${status_file}" "${sample_id}" "HISAT-genotype" "MISSING_INPUT" "No arcas extracted FASTQ found; not using genome BAM"
        return 0
    fi
    if [[ "${DRY_RUN}" == "true" ]]; then
        if [[ "${HLA_SINGLE_END}" == "true" ]]; then
            run_or_print "${HISAT_CMD}" --base hla -x "${HISAT_GENOTYPE_GENOME}" --in-dir "" -U "${hisat_inputs[0]}" --single-end --out-dir "${out_dir}/hisat" -z "${HISAT_INDEX_DIR}" -p "${HLA_THREADS}"
        else
            run_or_print "${HISAT_CMD}" --base hla -x "${HISAT_GENOTYPE_GENOME}" --in-dir "" -1 "${hisat_inputs[0]}" -2 "${hisat_inputs[1]}" --out-dir "${out_dir}/hisat" -z "${HISAT_INDEX_DIR}" -p "${HLA_THREADS}"
        fi
        write_status "${status_file}" "${sample_id}" "HISAT-genotype" "DRY_RUN" "${HISAT_CMD}"
        return 0
    else
        local hisat_dir="${out_dir}/hisat/run_$(date -u +%Y%m%dT%H%M%SZ)"
        mkdir -p "${hisat_dir}"
        local hisat_args=(--base hla -x "${HISAT_GENOTYPE_GENOME}" --in-dir "" --out-dir "${hisat_dir}" -z "${HISAT_INDEX_DIR}" -p "${HLA_THREADS}")
        if [[ "${HLA_SINGLE_END}" == "true" ]]; then
            hisat_args+=(-U "${hisat_inputs[0]}" --single-end)
        else
            hisat_args+=(-1 "${hisat_inputs[0]}" -2 "${hisat_inputs[1]}")
        fi
        if ! "${HISAT_CMD}" "${hisat_args[@]}" > "${hisat_dir}/hisat_genotype.log" 2>&1; then
            write_status "${status_file}" "${sample_id}" "HISAT-genotype" "FAIL" "HISAT-genotype failed; see ${hisat_dir}/hisat_genotype.log"
            return 0
        fi
    fi
    local hisat_log
    hisat_log="$(find "${out_dir}/hisat" -type f -name 'hisat_genotype.log' -printf '%T@ %p\n' | sort -nr | awk 'NR == 1 { print $2 }')"
    python3 "${HLA_PY}" normalize --caller hisat --sample "${sample_id}" --input "${hisat_log}" --output "${out_dir}/normalized/HISAT-genotype.tsv"
    if [[ "$(tsv_data_rows "${out_dir}/normalized/HISAT-genotype.tsv")" -gt 0 ]]; then
        write_status "${status_file}" "${sample_id}" "HISAT-genotype" "OK" "${hisat_log}"
    else
        write_status "${status_file}" "${sample_id}" "HISAT-genotype" "NO_CALLS" "HISAT finished but no HLA alleles were parsed; see ${hisat_log}"
    fi
}

process_sample() {
    local sample_id
    sample_id="$(require_sample_id "$1")"
    local work_dir="${WORK_ROOT}/${sample_id}"
    local out_dir="${work_dir}/results_hla"
    local status_file="${out_dir}/hla_status.tsv"
    local bam_uri

    log "HLA start: ${sample_id}"
    bam_uri="$(resolve_bam_uri "${sample_id}")" || {
        mkdir -p "${out_dir}"
        printf 'sample\tcaller\tstatus\tdetail\n' > "${status_file}"
        write_status "${status_file}" "${sample_id}" "input" "MISSING_INPUT" "No BAM found in ${SOURCE_S3}/${sample_id}"
        upload_results "${out_dir}" "${DEST_S3}/${sample_id}/results_hla/"
        return 0
    }

    mkdir -p "${out_dir}/"{bam,arcas,optitype,hisat,normalized,logs}
    printf 'sample\tcaller\tstatus\tdetail\n' > "${status_file}"
    printf 'sample_id=%s\nsource_bam=%s\ncreated_utc=%s\n' "${sample_id}" "${bam_uri}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "${out_dir}/manifest.txt"

    if ! has_available_hla_caller; then
        write_status "${status_file}" "${sample_id}" "arcasHLA" "MISSING_TOOL" "${ARCASHLA_CMD}"
        write_status "${status_file}" "${sample_id}" "OptiType" "MISSING_TOOL" "${OPTITYPE_CMD}"
        write_status "${status_file}" "${sample_id}" "HISAT-genotype" "MISSING_TOOL" "${HISAT_CMD}"
        printf 'sample\tlocus\tallele\tconsensus_status\tsupport_count\tsupport_callers\tcaller_count\tsource_files\n' > "${out_dir}/hla_consensus.tsv"
        log "No enabled HLA caller is available on PATH; uploading status without downloading BAM"
        upload_results "${out_dir}" "${DEST_S3}/${sample_id}/results_hla/"
        log "HLA done: ${sample_id}"
        return 0
    fi

    local bam="${out_dir}/bam/${sample_id}.bam"
    if [[ "${DRY_RUN}" == "true" ]]; then
        log "DRY_RUN input BAM: ${bam_uri}"
    elif [[ -s "${bam}" ]]; then
        log "Using existing local BAM: ${bam}"
    else
        run_or_print aws s3 cp "${bam_uri}" "${bam}" --profile "${AWS_PROFILE}" --region "${AWS_REGION}" --no-progress
    fi

    run_arcas "${sample_id}" "${bam}" "${out_dir}" "${status_file}"
    run_optitype "${sample_id}" "${bam}" "${out_dir}" "${status_file}"
    run_hisat "${sample_id}" "${bam}" "${out_dir}" "${status_file}"

    local normalized_files=("${out_dir}"/normalized/*.tsv)
    if [[ -e "${normalized_files[0]}" ]]; then
        python3 "${HLA_PY}" consensus --inputs "${normalized_files[@]}" --output "${out_dir}/hla_consensus.tsv"
    else
        printf 'sample\tlocus\tallele\tconsensus_status\tsupport_count\tsupport_callers\tcaller_count\tsource_files\n' > "${out_dir}/hla_consensus.tsv"
    fi

    upload_results "${out_dir}" "${DEST_S3}/${sample_id}/results_hla/"
    log "HLA done: ${sample_id}"
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

samples=("$@")
if [[ ${#samples[@]} -eq 0 && -n "${SAMPLES:-}" ]]; then
    # shellcheck disable=SC2206
    samples=(${SAMPLES})
fi
if [[ ${#samples[@]} -eq 0 ]]; then
    pipeline_die "No samples provided. Pass sample ids as arguments or SAMPLES='SRR... SRR...'."
fi

for sample in "${samples[@]}"; do
    process_sample "${sample}"
done
