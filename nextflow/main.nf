#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {
    CLEAN_IDS
    PREPARE_SNP_VCF
    CONSENSUS_H1
    CONSENSUS_H2
    TRANSLATE_H1
    TRANSLATE_H2
    CLEAN_SNP_PROTEINS
    FRAMESHIFT_PROTEOME
    COMBINE_PROTEOME
    VERIFY_PROTEOME
} from './modules/proteome.nf'

include {
    DOWNLOAD_FASTQ
    FASTQC_READS
    KALLISTO_QUANT
    STAR_ALIGN
    MARK_DUPLICATES
    INDEX_DEDUP_BAM
    GATK_SPLIT_N_CIGAR
    INDEX_SPLIT_BAM
    GATK_HAPLOTYPE_CALLER
    GATK_VARIANT_FILTRATION
} from './modules/upstream.nf'

include {
    AMYLOGRAM_PY
    AMYPRED_FRL
    AMYLOGRAM_R
    MERGE_AMYLOID_PREDICTIONS
    MERGE_AMYLOID_PREDICTIONS_AMYPRED_PY
    MERGE_AMYLOID_PREDICTIONS_FULL
    PROTEIN_FEATURES_LIGHT
    PROTEIN_FEATURES_FULL
} from './modules/amyloid.nf'

workflow {
    def usage = '''
    Amyloid mutant proteome Nextflow migration

    Usage:
      nextflow run nextflow/main.nf -profile docker --mode full --samples samples.txt
      nextflow run nextflow/main.nf -profile test,docker --mode prepared
      nextflow run nextflow/main.nf -profile test,docker --mode amyloid --input_fasta PATH

    Required for full mode:
      --samples PATH          Text file with one SRA accession per line
      --ref_dir PATH          Reference directory containing genome FASTA sidecars for GATK
      --ref_genome_name NAME  Reference genome FASTA filename inside --ref_dir
      --ref_gtf PATH          Reference GTF annotation
      --star_index PATH       STAR genome index directory
      --kallisto_index PATH   Kallisto transcriptome index
      --sixmer_table PATH     AmyloGram-Py six-mer lookup table
      --outdir PATH           Output directory

    Required for prepared mode:
      --samples PATH       CSV with sample_id, abundance, vcf columns
      --ref_genome PATH    Reference genome FASTA
      --ref_gtf PATH       Reference GTF annotation
      --outdir PATH        Output directory

    Required for amyloid mode:
      --input_fasta PATH    combine_proteome.fasta input
      --sixmer_table PATH   AmyloGram-Py six-mer lookup table

    Tunable prepared-mode params:
      --min_tpm FLOAT
      --min_cds_bp INT
      --min_frameshift_aa INT
      --run_amypred true|false
      --run_amylogram_r true|false
      --run_protein_features true|false
    '''

    if (params.help) {
        log.info usage.strip()
        return
    }

    if (!params.outdir) {
        error "Workflow requires --outdir."
    }

    if (params.mode == 'full') {
        if (!params.samples) {
            error "Full mode requires --samples."
        }
        if (!params.ref_dir || !params.ref_genome_name || !params.ref_gtf || !params.star_index || !params.kallisto_index) {
            error "Full mode requires --ref_dir, --ref_genome_name, --ref_gtf, --star_index, and --kallisto_index."
        }
        if (!params.sixmer_table) {
            error "Full mode requires --sixmer_table for AmyloGram-Py."
        }

        full_samples = Channel
            .fromPath(params.samples)
            .splitText()
            .map { line -> line.trim() }
            .filter { line -> line && !line.startsWith('#') }

        fastq = DOWNLOAD_FASTQ(full_samples)
        qc = FASTQC_READS(fastq)
        expression = KALLISTO_QUANT(qc, file(params.kallisto_index))
        aligned = STAR_ALIGN(expression, file(params.star_index))
        dedup = MARK_DUPLICATES(aligned)
        indexed = INDEX_DEDUP_BAM(dedup)
        split_bam = GATK_SPLIT_N_CIGAR(indexed, file(params.ref_dir))
        indexed_split = INDEX_SPLIT_BAM(split_bam)
        variants = GATK_HAPLOTYPE_CALLER(indexed_split, file(params.ref_dir))
        filtered = GATK_VARIANT_FILTRATION(variants, file(params.ref_dir))

        full_prepared_samples = filtered.map { sample_id, abundance, vcf ->
            tuple(
                sample_id as String,
                abundance,
                vcf,
                file("${params.ref_dir}/${params.ref_genome_name}"),
                file(params.ref_gtf),
                file("${params.repo_dir}/proteome_pipeline")
            )
        }

        clean_ids = CLEAN_IDS(full_prepared_samples)
        ready_snps = PREPARE_SNP_VCF(clean_ids)
        consensus_h1 = CONSENSUS_H1(ready_snps)
        consensus_h2 = CONSENSUS_H2(ready_snps)
        translated_h1 = TRANSLATE_H1(consensus_h1)
        translated_h2 = TRANSLATE_H2(consensus_h2)
        translated = translated_h1.join(translated_h2)
        snp_proteins = CLEAN_SNP_PROTEINS(translated)
        frameshifts = FRAMESHIFT_PROTEOME(clean_ids)
        combined = COMBINE_PROTEOME(snp_proteins.join(frameshifts))
        verified = VERIFY_PROTEOME(combined)

        full_amyloid_input = verified.map { sample_id, results_dir ->
            tuple(sample_id, file("${results_dir}/combine_proteome.fasta"), file(params.sixmer_table))
        }
        full_amypred_input = verified.map { sample_id, results_dir ->
            tuple(sample_id, file("${results_dir}/combine_proteome.fasta"), file("${params.repo_dir}/amyloid_predictors/AMYPred-FRL"))
        }
        full_amylogram_r_input = verified.map { sample_id, results_dir ->
            tuple(sample_id, file("${results_dir}/combine_proteome.fasta"), file("${params.repo_dir}/amyloid_predictors/AmyloGram"))
        }
        full_features_input = verified.map { sample_id, results_dir ->
            tuple(sample_id, file("${results_dir}/combine_proteome.fasta"), file("${params.repo_dir}/amyloid_predictors/ProteinFeatures"))
        }
        full_light_feature_input = verified.map { sample_id, results_dir ->
            tuple(sample_id, file("${results_dir}/combine_proteome.fasta"))
        }
        amypred = AMYPRED_FRL(full_amypred_input)
        amylogram_r = AMYLOGRAM_R(full_amylogram_r_input)
        amylogram_py = AMYLOGRAM_PY(full_amyloid_input)
        full_merge_input = amypred
            .join(amylogram_r)
            .join(amylogram_py)
            .map { sample_id, input_fasta, amypred_csv, amylogram_r_csv, amylogram_r_summary, py_input_fasta, amylogram_py_csv, amylogram_report_json, amylogram_report_md, amylogram_skipped_tsv, amylogram_top_tsv ->
                tuple(sample_id, input_fasta, amypred_csv, amylogram_r_csv, amylogram_r_summary, amylogram_py_csv, amylogram_report_json, amylogram_report_md, amylogram_skipped_tsv, amylogram_top_tsv)
            }
        MERGE_AMYLOID_PREDICTIONS_FULL(full_merge_input)
        PROTEIN_FEATURES_FULL(full_features_input)
        PROTEIN_FEATURES_LIGHT(full_light_feature_input)
    } else if (params.mode == 'prepared') {
        if (!params.ref_genome || !params.ref_gtf) {
            error "Prepared mode requires --ref_genome and --ref_gtf."
        }
        if (!params.samples) {
            error "Prepared mode requires --samples."
        }
        if (!params.sixmer_table) {
            error "Prepared mode requires --sixmer_table for the prepared end-to-end amyloid branch."
        }
        if (params.min_tpm as BigDecimal < 0) {
            error "Parameter --min_tpm must be >= 0."
        }
        if (params.min_cds_bp as Integer < 1) {
            error "Parameter --min_cds_bp must be >= 1."
        }
        if (params.min_frameshift_aa as Integer < 1) {
            error "Parameter --min_frameshift_aa must be >= 1."
        }
        if (params.run_amylogram_r && !params.run_amypred) {
            error "Prepared mode currently requires --run_amypred true when --run_amylogram_r true."
        }

        prepared_samples = Channel
            .fromPath(params.samples)
            .splitCsv(header: true)
            .filter { row -> row.sample_id && row.abundance && row.vcf }
            .map { row ->
                tuple(
                    row.sample_id as String,
                    file(row.abundance),
                    file(row.vcf),
                    file(params.ref_genome),
                    file(params.ref_gtf),
                    file("${params.repo_dir}/proteome_pipeline")
                )
            }

        clean_ids = CLEAN_IDS(prepared_samples)
        ready_snps = PREPARE_SNP_VCF(clean_ids)
        consensus_h1 = CONSENSUS_H1(ready_snps)
        consensus_h2 = CONSENSUS_H2(ready_snps)
        translated_h1 = TRANSLATE_H1(consensus_h1)
        translated_h2 = TRANSLATE_H2(consensus_h2)
        translated = translated_h1.join(translated_h2)
        snp_proteins = CLEAN_SNP_PROTEINS(translated)
        frameshifts = FRAMESHIFT_PROTEOME(clean_ids)
        combined = COMBINE_PROTEOME(snp_proteins.join(frameshifts))
        verified = VERIFY_PROTEOME(combined)

        prepared_amyloid_input = verified.map { sample_id, results_dir ->
            tuple(sample_id, file("${results_dir}/combine_proteome.fasta"), file(params.sixmer_table))
        }
        prepared_amypred_input = verified.map { sample_id, results_dir ->
            tuple(sample_id, file("${results_dir}/combine_proteome.fasta"), file("${params.repo_dir}/amyloid_predictors/AMYPred-FRL"))
        }
        prepared_amylogram_r_input = verified.map { sample_id, results_dir ->
            tuple(sample_id, file("${results_dir}/combine_proteome.fasta"), file("${params.repo_dir}/amyloid_predictors/AmyloGram"))
        }
        prepared_feature_input = verified.map { sample_id, results_dir ->
            tuple(sample_id, file("${results_dir}/combine_proteome.fasta"))
        }
        prepared_full_feature_input = verified.map { sample_id, results_dir ->
            tuple(sample_id, file("${results_dir}/combine_proteome.fasta"), file("${params.repo_dir}/amyloid_predictors/ProteinFeatures"))
        }
        amylogram_py = AMYLOGRAM_PY(prepared_amyloid_input)
        if (params.run_amypred && params.run_amylogram_r) {
            amypred = AMYPRED_FRL(prepared_amypred_input)
            amylogram_r = AMYLOGRAM_R(prepared_amylogram_r_input)
            prepared_full_merge_input = amypred
                .join(amylogram_r)
                .join(amylogram_py)
                .map { sample_id, input_fasta, amypred_csv, amylogram_r_csv, amylogram_r_summary, py_input_fasta, amylogram_py_csv, amylogram_report_json, amylogram_report_md, amylogram_skipped_tsv, amylogram_top_tsv ->
                    tuple(sample_id, input_fasta, amypred_csv, amylogram_r_csv, amylogram_r_summary, amylogram_py_csv, amylogram_report_json, amylogram_report_md, amylogram_skipped_tsv, amylogram_top_tsv)
                }
            MERGE_AMYLOID_PREDICTIONS_FULL(prepared_full_merge_input)
        } else if (params.run_amypred) {
            amypred = AMYPRED_FRL(prepared_amypred_input)
            prepared_merge_input = amypred
                .join(amylogram_py)
                .map { sample_id, input_fasta, amypred_csv, py_input_fasta, amylogram_py_csv, amylogram_report_json, amylogram_report_md, amylogram_skipped_tsv, amylogram_top_tsv ->
                    tuple(sample_id, input_fasta, amypred_csv, amylogram_py_csv, amylogram_report_json, amylogram_report_md, amylogram_skipped_tsv, amylogram_top_tsv)
                }
            MERGE_AMYLOID_PREDICTIONS_AMYPRED_PY(prepared_merge_input)
        } else {
            MERGE_AMYLOID_PREDICTIONS(amylogram_py)
        }
        if (params.run_protein_features) {
            PROTEIN_FEATURES_FULL(prepared_full_feature_input)
        }
        PROTEIN_FEATURES_LIGHT(prepared_feature_input)
    } else if (params.mode == 'amyloid') {
        if (!params.input_fasta) {
            error "Amyloid mode requires --input_fasta."
        }
        if (!params.sixmer_table) {
            error "Amyloid mode requires --sixmer_table."
        }
        if (!params.run_amylogram_py) {
            error "Phase 4 amyloid smoke requires --run_amylogram_py true."
        }
        if (params.amylogram_py_top_k as Integer < 0) {
            error "Parameter --amylogram_py_top_k must be >= 0."
        }

        amyloid_input = Channel.of(tuple('fixture_minimal', file(params.input_fasta), file(params.sixmer_table)))
        amylogram_py = AMYLOGRAM_PY(amyloid_input)
        MERGE_AMYLOID_PREDICTIONS(amylogram_py)
    } else {
        error "Unsupported mode '${params.mode}'. Supported modes: full, prepared, amyloid."
    }
}
