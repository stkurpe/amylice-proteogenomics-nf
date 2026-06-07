#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<EOF
Usage:
  ./run_complete_processing_pipeline.sh --samples samples.txt --s3-bucket s3://bucket/prefix

Runs the full workflow for each sample:
  download -> QC -> expression -> alignment/variants -> mutant proteome ->
  amyloidogenicity predictors -> protein physicochemical features -> S3 archive.

This is a named wrapper around run_full_pipeline.sh with protein features enabled
by default. Pass through any run_full_pipeline.sh option after the arguments above.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

RUN_PROTEIN_FEATURES=true exec "${SCRIPT_DIR}/run_full_pipeline.sh" "$@"
