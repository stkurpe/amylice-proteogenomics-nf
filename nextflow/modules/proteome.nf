process CLEAN_IDS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src)

    output:
    tuple val(sample_id), path('clean_ids.txt'), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src)

    script:
    """
    set -Eeuo pipefail
    PYTHONPATH="\$PWD" python3 -m proteome_pipeline.cli clean-ids \
      --gtf "${gtf}" \
      --abundance "${abundance}" \
      --out clean_ids.txt \
      --min-tpm "${params.min_tpm}" \
      --min-cds-bp "${params.min_cds_bp}"
    test -s clean_ids.txt
    """
}

process PREPARE_SNP_VCF {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(clean_ids), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src)

    output:
    tuple val(sample_id), path(clean_ids), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src), path('ready_snps.vcf.gz'), path('ready_snps.vcf.gz.tbi')

    script:
    """
    set -Eeuo pipefail
    samtools faidx "${genome}"
    bcftools view -v snps "${vcf}" | bcftools norm -m -any -f "${genome}" -Oz -o ready_snps.vcf.gz
    bcftools index -t -f ready_snps.vcf.gz
    test -s ready_snps.vcf.gz
    test -s ready_snps.vcf.gz.tbi
    """
}

process CONSENSUS_H1 {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(clean_ids), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src), path(ready_snps), path(ready_snps_tbi)

    output:
    tuple val(sample_id), path(clean_ids), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src), path('genome_h1.fa'), path('genome_h1.fa.fai')

    script:
    """
    set -Eeuo pipefail
    bcftools consensus -H 1 -f "${genome}" "${ready_snps}" -o genome_h1.fa
    samtools faidx genome_h1.fa
    test -s genome_h1.fa
    test -s genome_h1.fa.fai
    """
}

process CONSENSUS_H2 {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(clean_ids), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src), path(ready_snps), path(ready_snps_tbi)

    output:
    tuple val(sample_id), path(clean_ids), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src), path('genome_h2.fa'), path('genome_h2.fa.fai')

    script:
    """
    set -Eeuo pipefail
    bcftools consensus -H 2 -f "${genome}" "${ready_snps}" -o genome_h2.fa
    samtools faidx genome_h2.fa
    test -s genome_h2.fa
    test -s genome_h2.fa.fai
    """
}

process TRANSLATE_H1 {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(clean_ids), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src), path(genome_h1), path(genome_h1_fai)

    output:
    tuple val(sample_id), path(clean_ids), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src), path('prot_h1.fa')

    script:
    """
    set -Eeuo pipefail
    gffread "${gtf}" -g "${genome_h1}" -y prot_h1.raw.fa --ids "${clean_ids}"
    sed 's/>/>H1_/' prot_h1.raw.fa > prot_h1.fa
    test -s prot_h1.fa
    """
}

process TRANSLATE_H2 {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(clean_ids), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src), path(genome_h2), path(genome_h2_fai)

    output:
    tuple val(sample_id), path('prot_h2.fa')

    script:
    """
    set -Eeuo pipefail
    gffread "${gtf}" -g "${genome_h2}" -y prot_h2.raw.fa --ids "${clean_ids}"
    sed 's/>/>H2_/' prot_h2.raw.fa > prot_h2.fa
    test -s prot_h2.fa
    """
}

process CLEAN_SNP_PROTEINS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(clean_ids), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src), path(prot_h1), path(prot_h2)

    output:
    tuple val(sample_id), path(clean_ids), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src), path('protein.fasta'), path('nonsense_candidates.txt')

    script:
    """
    set -Eeuo pipefail
    cat "${prot_h1}" "${prot_h2}" > protein_raw.fasta
    PYTHONPATH="\$PWD" python3 -m proteome_pipeline.cli clean-proteins \
      --input protein_raw.fasta \
      --out protein.fasta
    PYTHONPATH="\$PWD" python3 -m proteome_pipeline.cli nonsense-report \
      --input protein.fasta \
      --out nonsense_candidates.txt
    test -s protein.fasta
    test -s nonsense_candidates.txt
    """
}

process FRAMESHIFT_PROTEOME {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(clean_ids), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src)

    output:
    tuple val(sample_id), path('frameshift_unique.fasta')

    script:
    """
    set -Eeuo pipefail
    samtools faidx "${genome}"
    PYTHONPATH="\$PWD" python3 -m proteome_pipeline.cli frameshifts \
      --vcf "${vcf}" \
      --gtf "${gtf}" \
      --genome "${genome}" \
      --ids "${clean_ids}" \
      --out frameshift_unique.fasta \
      --min-aa "${params.min_frameshift_aa}"
    test -f frameshift_unique.fasta
    """
}

process COMBINE_PROTEOME {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(clean_ids), path(abundance), path(vcf), path(genome), path(gtf), path(pipeline_src), path(protein_fasta), path(nonsense_candidates), path(frameshift_fasta)

    output:
    tuple val(sample_id), path(clean_ids), path(protein_fasta), path(frameshift_fasta), path(nonsense_candidates), path(pipeline_src), path('combine_proteome.fasta')

    script:
    """
    set -Eeuo pipefail
    PYTHONPATH="\$PWD" python3 -m proteome_pipeline.cli combine \
      --out combine_proteome.fasta \
      "${protein_fasta}" \
      "${frameshift_fasta}"
    test -s combine_proteome.fasta
    """
}

process VERIFY_PROTEOME {
    tag "$sample_id"
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${sample_id}/${filename}" }

    input:
    tuple val(sample_id), path(clean_ids), path(protein_fasta), path(frameshift_fasta), path(nonsense_candidates), path(pipeline_src), path(combine_fasta)

    output:
    tuple val(sample_id), path('results_proteins')

    script:
    """
    set -Eeuo pipefail
    mkdir -p results_proteins
    cp "${clean_ids}" results_proteins/clean_ids.txt
    cp "${protein_fasta}" results_proteins/protein.fasta
    cp "${frameshift_fasta}" results_proteins/frameshift_unique.fasta
    cp "${nonsense_candidates}" results_proteins/nonsense_candidates.txt
    cp "${combine_fasta}" results_proteins/combine_proteome.fasta
    cat > results_proteins/manifest.txt <<EOF
input_source=nextflow_prepared_fixture_or_manifest
sample=${sample_id}
min_tpm=${params.min_tpm}
min_cds_bp=${params.min_cds_bp}
min_frameshift_aa=${params.min_frameshift_aa}
EOF
    PYTHONPATH="\$PWD" python3 -m proteome_pipeline.cli verify \
      --proteo-dir results_proteins \
      --out results_proteins/verification_report.tsv
    test -s results_proteins/verification_report.tsv
    ! grep -q '^ATTENTION' results_proteins/verification_report.tsv
    """
}
