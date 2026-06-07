# Product Pipeline Defaults

Updated: 2026-05-19

## Amyloidogenicity predictors

Default predictors:

- AMYPred-FRL: enabled by default (`RUN_AMYPRED=1`).
- AmyloGram-Py: primary AmyloGram implementation, enabled by default (`RUN_AMYLOGRAM_PY=1`).
- Legacy R AmyloGram: disabled by default (`RUN_AMYLOGRAM=0`). Enable explicitly with `--run-amylogram` or `RUN_AMYLOGRAM=1`.

Typical run:

```bash
/home/codex/run_amyloid_predictors.sh --sample SRR32060234
```

Run with legacy R AmyloGram as an extra predictor:

```bash
/home/codex/run_amyloid_predictors.sh --sample SRR32060234 --run-amylogram
```

Skip AmyloGram-Py if needed:

```bash
/home/codex/run_amyloid_predictors.sh --sample SRR32060234 --skip-amylogram-py
```

## HLA consensus

Default HLA callers are all enabled:

- arcasHLA (`RUN_ARCAS=true`)
- OptiType (`RUN_OPTITYPE=true`)
- HISAT-genotype (`RUN_HISAT=true`)

The workflow uses the STAR/GATK BAM as the source for arcasHLA extraction. OptiType and HISAT-genotype then consume the arcasHLA extracted FASTQ rather than the generic genome BAM.

Important HISAT-genotype defaults:

- `HISAT_GENOTYPE_GENOME=genotype_genome`
- `HISAT_INDEX_DIR=/home/codex/tools/src/hisat-genotype/indicies`
- `--in-dir ""` is passed so absolute FASTQ paths are not corrupted by HISAT-genotype path joining.

Dry-run example:

```bash
PATH="/home/codex/micromamba/envs/hla-arcas/bin:/home/codex/micromamba/envs/hla-optitype/bin:/home/codex/micromamba/envs/hla-hisat/bin:/home/codex/tools/bin:$PATH" DRY_RUN=true SAMPLES="SRR32060234" /home/codex/pipeline_scripts/05_hla_consensus.sh
```

Outputs are uploaded to:

```text
s3://codex-test-ngsdata-calculations/prepared-bioinfo-proteome/<SAMPLE_ID>/results_hla/
```

The upload excludes local input BAMs and debug scratch directories.
