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
    AMYLOGRAM_PY
    MERGE_AMYLOID_PREDICTIONS
    PROTEIN_FEATURES
} from './modules/amyloid.nf'

workflow {
    def usage = '''
    Amyloid mutant proteome Nextflow migration

    Usage:
      nextflow run nextflow/main.nf -profile test,docker --mode prepared
      nextflow run nextflow/main.nf -profile test,docker --mode amyloid --input_fasta PATH

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
    '''

    if (params.help) {
        log.info usage.strip()
        return
    }

    if (!params.outdir) {
        error "Workflow requires --outdir."
    }

    if (params.mode == 'prepared') {
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
        prepared_feature_input = verified.map { sample_id, results_dir ->
            tuple(sample_id, file("${results_dir}/combine_proteome.fasta"))
        }
        amylogram_py = AMYLOGRAM_PY(prepared_amyloid_input)
        MERGE_AMYLOID_PREDICTIONS(amylogram_py)
        PROTEIN_FEATURES(prepared_feature_input)
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
        error "Unsupported mode '${params.mode}'. Supported modes: prepared, amyloid."
    }
}
