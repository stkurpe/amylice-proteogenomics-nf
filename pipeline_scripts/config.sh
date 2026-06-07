#!/bin/bash

# Core pipeline configuration. Defaults keep project data under /home/codex.
export PROJECT_DIR="${PROJECT_DIR:-/home/codex}"
export BASE_DIR="${BASE_DIR:-${PROJECT_DIR}}"
export REF_DIR="${REF_DIR:-${PROJECT_DIR}/references}"

export THREADS="${THREADS:-8}"
export JAVA_OPTS="${JAVA_OPTS:--Xmx12g}"
export CLEANUP_AFTER_UPLOAD="${CLEANUP_AFTER_UPLOAD:-false}"
export PREPARE_REFS="${PREPARE_REFS:-true}"

export DOCKER_OPTS="${DOCKER_OPTS:---rm -v ${BASE_DIR}:/data -w /data --user $(id -u):$(id -g)}"

export REF_GENOME="${REF_GENOME:-GRCh38.primary_assembly.genome.fa}"
export REF_GTF="${REF_GTF:-gencode.v46.primary_assembly.annotation.gtf}"
export REF_VCF="${REF_VCF:-common_all_20180418.vcf}"
export REF_TRANSCRIPT="${REF_TRANSCRIPT:-gencode.v46.pc_transcripts.fa}"
