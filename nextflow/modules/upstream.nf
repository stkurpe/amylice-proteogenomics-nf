process DOWNLOAD_FASTQ {
    tag "$sample_id"
    container 'staphb/sratoolkit:3.3.0'

    input:
    val sample_id

    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz")

    script:
    """
    set -Eeuo pipefail
    fasterq-dump --split-files --outdir . --threads "${params.threads}" --progress "${sample_id}"
    if [ -f "${sample_id}_1.fastq" ]; then
      mv "${sample_id}_1.fastq" "${sample_id}.fastq"
    fi
    if [ -f "${sample_id}_2.fastq" ]; then
      gzip -f "${sample_id}_2.fastq"
    fi
    test -s "${sample_id}.fastq"
    gzip -f "${sample_id}.fastq"
    test -s "${sample_id}.fastq.gz"
    """
}

process FASTQC_READS {
    tag "$sample_id"
    container 'quay.io/biocontainers/fastqc:0.11.9--0'
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${sample_id}/${filename}" }

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path(fastq), path('results_qc')

    script:
    """
    set -Eeuo pipefail
    mkdir -p results_qc
    fastqc -o results_qc -t "${params.threads}" "${fastq}"
    test -d results_qc
    """
}

process KALLISTO_QUANT {
    tag "$sample_id"
    container 'quay.io/biocontainers/kallisto:0.46.1--h4f7b962_0'
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${sample_id}/${filename}" }

    input:
    tuple val(sample_id), path(fastq), path(qc_dir)
    path kallisto_index

    output:
    tuple val(sample_id), path('results_expression/abundance.tsv'), path(fastq)

    script:
    """
    set -Eeuo pipefail
    mkdir -p results_expression
    kallisto quant \
      -i "${kallisto_index}" \
      -o results_expression \
      --single -l "${params.kallisto_fragment_length}" -s "${params.kallisto_fragment_sd}" \
      "${fastq}"
    test -s results_expression/abundance.tsv
    """
}

process STAR_ALIGN {
    tag "$sample_id"
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${sample_id}/${filename}" }

    input:
    tuple val(sample_id), path(abundance), path(fastq)
    path star_index

    output:
    tuple val(sample_id), path(abundance), path('results_star/Aligned.sortedByCoord.out.bam')

    script:
    """
    set -Eeuo pipefail
    mkdir -p results_star
    STAR --runThreadN "${params.threads}" \
      --genomeDir "${star_index}" \
      --readFilesIn "${fastq}" \
      --readFilesCommand zcat \
      --outFileNamePrefix "results_star/" \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattrRGline "ID:${sample_id}" "SM:${sample_id}" PL:ILLUMINA
    test -s results_star/Aligned.sortedByCoord.out.bam
    """
}

process MARK_DUPLICATES {
    tag "$sample_id"
    container 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${sample_id}/${filename}" }

    input:
    tuple val(sample_id), path(abundance), path(aligned_bam)

    output:
    tuple val(sample_id), path(abundance), path('results_star/dedup.bam')

    script:
    """
    set -Eeuo pipefail
    mkdir -p results_star
    picard MarkDuplicates \
      "I=${aligned_bam}" \
      "O=results_star/dedup.bam" \
      "M=results_star/metrics.txt"
    test -s results_star/dedup.bam
    """
}

process INDEX_DEDUP_BAM {
    tag "$sample_id"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(sample_id), path(abundance), path(dedup_bam)

    output:
    tuple val(sample_id), path(abundance), path(dedup_bam), path('results_star/dedup.bam.bai')

    script:
    """
    set -Eeuo pipefail
    samtools index "${dedup_bam}"
    test -s results_star/dedup.bam.bai
    """
}

process GATK_SPLIT_N_CIGAR {
    tag "$sample_id"
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${sample_id}/${filename}" }

    input:
    tuple val(sample_id), path(abundance), path(dedup_bam), path(dedup_bai)
    path ref_dir

    output:
    tuple val(sample_id), path(abundance), path('results_gatk/split.bam')

    script:
    """
    set -Eeuo pipefail
    mkdir -p results_gatk
    gatk SplitNCigarReads \
      -R "${ref_dir}/${params.ref_genome_name}" \
      -I "${dedup_bam}" \
      -O results_gatk/split.bam
    test -s results_gatk/split.bam
    """
}

process INDEX_SPLIT_BAM {
    tag "$sample_id"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(sample_id), path(abundance), path(split_bam)

    output:
    tuple val(sample_id), path(abundance), path(split_bam), path('results_gatk/split.bam.bai')

    script:
    """
    set -Eeuo pipefail
    samtools index "${split_bam}"
    test -s results_gatk/split.bam.bai
    """
}

process GATK_HAPLOTYPE_CALLER {
    tag "$sample_id"
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(sample_id), path(abundance), path(split_bam), path(split_bai)
    path ref_dir

    output:
    tuple val(sample_id), path(abundance), path('variants.vcf.gz')

    script:
    """
    set -Eeuo pipefail
    gatk HaplotypeCaller \
      -R "${ref_dir}/${params.ref_genome_name}" \
      -I "${split_bam}" \
      -O variants.vcf.gz \
      --dont-use-soft-clipped-bases \
      -stand-call-conf 20.0
    test -s variants.vcf.gz
    """
}

process GATK_VARIANT_FILTRATION {
    tag "$sample_id"
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${sample_id}/${filename}/results_gatk/variants_filtered.vcf.gz" }

    input:
    tuple val(sample_id), path(abundance), path(variants_vcf)
    path ref_dir

    output:
    tuple val(sample_id), path(abundance), path('variants_filtered.vcf.gz')

    script:
    """
    set -Eeuo pipefail
    gatk VariantFiltration \
      -R "${ref_dir}/${params.ref_genome_name}" \
      -V "${variants_vcf}" \
      -O variants_filtered.vcf.gz \
      --window 35 --cluster 3 \
      --filter-name "FS" --filter-expression "FS > 30.0" \
      --filter-name "QD" --filter-expression "QD < 2.0"
    test -s variants_filtered.vcf.gz
    """
}
