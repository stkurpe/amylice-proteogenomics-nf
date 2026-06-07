#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLE_ID="${SAMPLE_ID:-SRR32060234}"
BASE_DIR="${BASE_DIR:-/home/codex}"
WORK_DIR="${WORK_DIR:-${BASE_DIR}/amyloid_runs/${SAMPLE_ID}}"
INPUT_FASTA="${INPUT_FASTA:-${WORK_DIR}/input/combine_proteome.fasta}"
INPUT_S3="${INPUT_S3:-s3://codex-test-ngsdata-calculations/public-sra-test/${SAMPLE_ID}/results_proteins/combine_proteome.fasta}"
OUTPUT_S3="${OUTPUT_S3:-s3://codex-test-ngsdata-calculations/public-sra-test/${SAMPLE_ID}/results_amyloid}"
PROTEIN_FEATURES_OUTPUT_S3="${PROTEIN_FEATURES_OUTPUT_S3:-}"
AWS_PROFILE="${AWS_PROFILE:-codex-sandbox}"
AWS_REGION="${AWS_REGION:-us-east-1}"

AMYPRED_IMAGE="${AMYPRED_IMAGE:-python:3.8-slim}"
AMYLOGRAM_IMAGE="${AMYLOGRAM_IMAGE:-codex/amylogram:1.1-r4.3.3}"
AMYLOGRAM_PY_IMAGE="${AMYLOGRAM_PY_IMAGE:-amylogram-py:local}"
AMYLOGRAM_PY_IMAGE_TAR="${AMYLOGRAM_PY_IMAGE_TAR:-/home/codex/docker_images/amylogram-py-local-linux-amd64.tar.gz}"
PROTEIN_FEATURES_IMAGE="${PROTEIN_FEATURES_IMAGE:-${AMYLOGRAM_IMAGE}}"
AMYLOGRAM_CORES="${AMYLOGRAM_CORES:-8}"
AMYLOGRAM_CHUNK_SIZE="${AMYLOGRAM_CHUNK_SIZE:-200}"
AMYLOGRAM_RETRY_SPLIT_DEPTH="${AMYLOGRAM_RETRY_SPLIT_DEPTH:-10}"
AMYLOGRAM_ALLOW_ERRORS="${AMYLOGRAM_ALLOW_ERRORS:-0}"
AMYLOGRAM_PY_TOP_K="${AMYLOGRAM_PY_TOP_K:-1000}"
AMYLOGRAM_PY_LOOKUP="${AMYLOGRAM_PY_LOOKUP:-}"
if [[ -z "$AMYLOGRAM_PY_LOOKUP" ]]; then
  if [[ -s "/home/codex/amylogram-py/tests/fixtures/amylogram_sixmer_probabilities.bin" ]]; then
    AMYLOGRAM_PY_LOOKUP="/home/codex/amylogram-py/tests/fixtures/amylogram_sixmer_probabilities.bin"
  else
    AMYLOGRAM_PY_LOOKUP="${BASE_DIR}/amylogram-py/tests/fixtures/amylogram_sixmer_probabilities.bin"
  fi
fi
PROTEIN_FEATURES_WINDOW="${PROTEIN_FEATURES_WINDOW:-7}"
PROTEIN_FEATURES_PH="${PROTEIN_FEATURES_PH:-7.0}"
PROTEIN_FEATURES_R_LIBS="${PROTEIN_FEATURES_R_LIBS:-}"

FORCE="${FORCE:-0}"
RUN_AMYPRED="${RUN_AMYPRED:-1}"
RUN_AMYLOGRAM="${RUN_AMYLOGRAM:-0}"
RUN_AMYLOGRAM_PY="${RUN_AMYLOGRAM_PY:-1}"
RUN_PROTEIN_FEATURES="${RUN_PROTEIN_FEATURES:-1}"
UPLOAD="${UPLOAD:-1}"

CODE_DIR="${WORK_DIR}/code"
AMYPRED_DIR="${AMYPRED_DIR:-${CODE_DIR}/AMYPred-FRL}"
AMYLOGRAM_DIR="${AMYLOGRAM_DIR:-${CODE_DIR}/AmyloGram}"
PROTEIN_FEATURES_DIR="${PROTEIN_FEATURES_DIR:-${CODE_DIR}/ProteinFeatures}"
AMYPRED_FALLBACK_DIR="${AMYPRED_FALLBACK_DIR:-${BASE_DIR}/amyloid_predictors/AMYPred-FRL}"
AMYLOGRAM_FALLBACK_DIR="${AMYLOGRAM_FALLBACK_DIR:-${BASE_DIR}/amyloid_predictors/AmyloGram}"
PROTEIN_FEATURES_FALLBACK_DIR="${PROTEIN_FEATURES_FALLBACK_DIR:-${BASE_DIR}/amyloid_predictors/ProteinFeatures}"
if [[ ! -s "${AMYPRED_FALLBACK_DIR}/predict.py" && -s "${SCRIPT_DIR}/amyloid_predictors/AMYPred-FRL/predict.py" ]]; then
  AMYPRED_FALLBACK_DIR="${SCRIPT_DIR}/amyloid_predictors/AMYPred-FRL"
fi
if [[ ! -s "${AMYLOGRAM_FALLBACK_DIR}/predict_amylogram_fast.R" && -s "${SCRIPT_DIR}/amyloid_predictors/AmyloGram/predict_amylogram_fast.R" ]]; then
  AMYLOGRAM_FALLBACK_DIR="${SCRIPT_DIR}/amyloid_predictors/AmyloGram"
fi
if [[ ! -s "${PROTEIN_FEATURES_FALLBACK_DIR}/calc_protein_features.R" && -s "${SCRIPT_DIR}/amyloid_predictors/ProteinFeatures/calc_protein_features.R" ]]; then
  PROTEIN_FEATURES_FALLBACK_DIR="${SCRIPT_DIR}/amyloid_predictors/ProteinFeatures"
fi
if [[ ! -s "${AMYLOGRAM_FALLBACK_DIR}/predict_amylogram_fast.R" ]]; then
  AMYLOGRAM_FALLBACK_DIR="${BASE_DIR}/amyloid_predictors_one_button_20260516T171835Z"
fi
INPUT_DIR="${WORK_DIR}/input"
RESULTS_DIR="${WORK_DIR}/results"
LOG_DIR="${WORK_DIR}/logs"
DEPS_DIR="${DEPS_DIR:-${BASE_DIR}/amyloid_runs/shared/python_deps_py38}"

RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
STATUS_TSV="${RESULTS_DIR}/amyloid_predictors_status.tsv"
SUMMARY_TSV="${RESULTS_DIR}/amyloid_predictors_summary.tsv"
COMBINED_CSV="${RESULTS_DIR}/amyloid_combined_predictions.csv"
MAIN_LOG="${LOG_DIR}/amyloid_predictors_${RUN_ID}.log"

usage() {
  cat <<EOF
Usage:
  /home/codex/run_amyloid_predictors.sh [options]

Options:
  --sample ID             Sample ID. Default: ${SAMPLE_ID}
  --input-fasta PATH      Local FASTA. Default: ${INPUT_FASTA}
  --input-s3 URI          Download FASTA if local input is missing.
  --output-s3 URI         Upload output prefix. Default: ${OUTPUT_S3}
  --protein-features-output-s3 URI
                          Upload protein feature outputs to this prefix.
  --force                 Recompute outputs even if they already exist.
  --no-upload             Do not upload results to S3.
  --skip-amypred          Skip AMYPred-FRL.
  --run-amylogram         Run legacy R AmyloGram predictor. Disabled by default.
  --skip-amylogram        Skip legacy R AmyloGram.
  --run-amylogram-py      Run primary fast AmyloGram-Py lookup predictor. Enabled by default.
  --skip-amylogram-py     Skip primary AmyloGram-Py.
  --run-protein-features  Run R protein physicochemical feature extraction. Enabled by default.
  --skip-protein-features Skip R protein physicochemical feature extraction.
  -h, --help              Show this help.

Environment overrides:
  AMYLOGRAM_CORES=8 AMYLOGRAM_CHUNK_SIZE=200 AMYLOGRAM_RETRY_SPLIT_DEPTH=10 AMYLOGRAM_ALLOW_ERRORS=0 AMYLOGRAM_PY_IMAGE=amylogram-py:local AMYLOGRAM_PY_LOOKUP=/home/codex/amylogram-py/tests/fixtures/amylogram_sixmer_probabilities.bin AMYLOGRAM_PY_TOP_K=1000 PROTEIN_FEATURES_IMAGE=codex/amylogram:1.1-r4.3.3 PROTEIN_FEATURES_WINDOW=7 PROTEIN_FEATURES_PH=7.0 PROTEIN_FEATURES_R_LIBS=/path/to/R/library AWS_PROFILE=codex-sandbox AWS_REGION=us-east-1
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample)
      SAMPLE_ID="$2"
      WORK_DIR="${BASE_DIR}/amyloid_runs/${SAMPLE_ID}"
      INPUT_FASTA="${WORK_DIR}/input/combine_proteome.fasta"
      INPUT_S3="s3://codex-test-ngsdata-calculations/public-sra-test/${SAMPLE_ID}/results_proteins/combine_proteome.fasta"
      OUTPUT_S3="s3://codex-test-ngsdata-calculations/public-sra-test/${SAMPLE_ID}/results_amyloid"
      shift 2
      ;;
    --input-fasta)
      INPUT_FASTA="$2"
      shift 2
      ;;
    --input-s3)
      INPUT_S3="$2"
      shift 2
      ;;
    --output-s3)
      OUTPUT_S3="$2"
      shift 2
      ;;
    --protein-features-output-s3)
      PROTEIN_FEATURES_OUTPUT_S3="$2"
      shift 2
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --no-upload)
      UPLOAD=0
      shift
      ;;
    --skip-amypred)
      RUN_AMYPRED=0
      shift
      ;;
    --skip-amylogram)
      RUN_AMYLOGRAM=0
      shift
      ;;
    --run-amylogram)
      RUN_AMYLOGRAM=1
      shift
      ;;
    --run-amylogram-py)
      RUN_AMYLOGRAM_PY=1
      shift
      ;;
    --skip-amylogram-py)
      RUN_AMYLOGRAM_PY=0
      shift
      ;;
    --run-protein-features)
      RUN_PROTEIN_FEATURES=1
      shift
      ;;
    --skip-protein-features)
      RUN_PROTEIN_FEATURES=0
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

CODE_DIR="${WORK_DIR}/code"
AMYPRED_DIR="${CODE_DIR}/AMYPred-FRL"
AMYLOGRAM_DIR="${CODE_DIR}/AmyloGram"
PROTEIN_FEATURES_DIR="${CODE_DIR}/ProteinFeatures"
AMYPRED_FALLBACK_DIR="${AMYPRED_FALLBACK_DIR:-${BASE_DIR}/amyloid_predictors/AMYPred-FRL}"
AMYLOGRAM_FALLBACK_DIR="${AMYLOGRAM_FALLBACK_DIR:-${BASE_DIR}/amyloid_predictors/AmyloGram}"
PROTEIN_FEATURES_FALLBACK_DIR="${PROTEIN_FEATURES_FALLBACK_DIR:-${BASE_DIR}/amyloid_predictors/ProteinFeatures}"
if [[ ! -s "${AMYPRED_FALLBACK_DIR}/predict.py" && -s "${SCRIPT_DIR}/amyloid_predictors/AMYPred-FRL/predict.py" ]]; then
  AMYPRED_FALLBACK_DIR="${SCRIPT_DIR}/amyloid_predictors/AMYPred-FRL"
fi
if [[ ! -s "${AMYLOGRAM_FALLBACK_DIR}/predict_amylogram_fast.R" && -s "${SCRIPT_DIR}/amyloid_predictors/AmyloGram/predict_amylogram_fast.R" ]]; then
  AMYLOGRAM_FALLBACK_DIR="${SCRIPT_DIR}/amyloid_predictors/AmyloGram"
fi
if [[ ! -s "${PROTEIN_FEATURES_FALLBACK_DIR}/calc_protein_features.R" && -s "${SCRIPT_DIR}/amyloid_predictors/ProteinFeatures/calc_protein_features.R" ]]; then
  PROTEIN_FEATURES_FALLBACK_DIR="${SCRIPT_DIR}/amyloid_predictors/ProteinFeatures"
fi
if [[ ! -s "${AMYLOGRAM_FALLBACK_DIR}/predict_amylogram_fast.R" ]]; then
  AMYLOGRAM_FALLBACK_DIR="${BASE_DIR}/amyloid_predictors_one_button_20260516T171835Z"
fi
INPUT_DIR="${WORK_DIR}/input"
RESULTS_DIR="${WORK_DIR}/results"
LOG_DIR="${WORK_DIR}/logs"
DEPS_DIR="${DEPS_DIR:-${BASE_DIR}/amyloid_runs/shared/python_deps_py38}"
STATUS_TSV="${RESULTS_DIR}/amyloid_predictors_status.tsv"
SUMMARY_TSV="${RESULTS_DIR}/amyloid_predictors_summary.tsv"
COMBINED_CSV="${RESULTS_DIR}/amyloid_combined_predictions.csv"
MAIN_LOG="${LOG_DIR}/amyloid_predictors_${RUN_ID}.log"

if [[ -z "$PROTEIN_FEATURES_OUTPUT_S3" && "$OUTPUT_S3" == */results_amyloid ]]; then
  PROTEIN_FEATURES_OUTPUT_S3="${OUTPUT_S3%/results_amyloid}/results_protein_features"
fi

mkdir -p "$INPUT_DIR" "$RESULTS_DIR" "$LOG_DIR"

if [[ ! -s "${AMYPRED_DIR}/predict.py" && -s "${AMYPRED_FALLBACK_DIR}/predict.py" ]]; then
  AMYPRED_DIR="$AMYPRED_FALLBACK_DIR"
fi
if [[ ! -s "${AMYLOGRAM_DIR}/predict_amylogram_fast.R" && -s "${AMYLOGRAM_FALLBACK_DIR}/predict_amylogram_fast.R" ]]; then
  AMYLOGRAM_DIR="$AMYLOGRAM_FALLBACK_DIR"
fi
if [[ ! -s "${PROTEIN_FEATURES_DIR}/calc_protein_features.R" && -s "${PROTEIN_FEATURES_FALLBACK_DIR}/calc_protein_features.R" ]]; then
  PROTEIN_FEATURES_DIR="$PROTEIN_FEATURES_FALLBACK_DIR"
fi

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$MAIN_LOG"
}

record_status() {
  local predictor="$1"
  local status="$2"
  local message="$3"
  printf '%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$predictor" "$status" "$message" >> "$STATUS_TSV"
}

run_step() {
  local predictor="$1"
  local label="$2"
  shift 2
  log "START ${predictor}: ${label}"
  set +e
  "$@" >> "${LOG_DIR}/${predictor}_${RUN_ID}.log" 2>&1
  local rc=$?
  set -e
  if [[ "$rc" -eq 0 ]]; then
    log "OK ${predictor}: ${label}"
    record_status "$predictor" "OK" "$label"
  else
    log "FAIL ${predictor}: ${label} rc=${rc}"
    record_status "$predictor" "FAIL" "${label} rc=${rc}"
  fi
  return "$rc"
}

line_count_without_header() {
  local path="$1"
  if [[ -s "$path" ]]; then
    local lines
    lines="$(wc -l < "$path")"
    if [[ "$lines" -gt 0 ]]; then
      echo $((lines - 1))
    else
      echo 0
    fi
  else
    echo 0
  fi
}

ensure_input() {
  if [[ -s "$INPUT_FASTA" ]]; then
    log "Input FASTA exists: ${INPUT_FASTA}"
    return 0
  fi
  log "Input FASTA missing; downloading from ${INPUT_S3}"
  aws s3 cp "$INPUT_S3" "$INPUT_FASTA" \
    --profile "$AWS_PROFILE" \
    --region "$AWS_REGION" \
    --quiet
}

run_amypred() {
  local output_csv="${RESULTS_DIR}/amypred_frl_prediction.csv"
  if [[ "$FORCE" != "1" && -s "$output_csv" ]]; then
    log "SKIP AMYPred-FRL: existing ${output_csv}"
    record_status "AMYPred-FRL" "SKIP" "existing output"
    return 0
  fi

  test -s "${AMYPRED_DIR}/predict.py"
  test -s "$INPUT_FASTA"

  run_step "AMYPred-FRL" "predict" \
    docker run --rm \
      --name "amypred-frl-${SAMPLE_ID}-${RUN_ID}" \
      --user "$(id -u):$(id -g)" \
      -e PYTHONPATH="/deps" \
      -e PIP_CACHE_DIR="/tmp/pip-cache" \
      -v "${WORK_DIR}:${WORK_DIR}" \
      -v "${BASE_DIR}/amyloid_predictors:${BASE_DIR}/amyloid_predictors:ro" \
      -v "${BASE_DIR}/amyloid_runs/shared:${BASE_DIR}/amyloid_runs/shared" \
      -w "${AMYPRED_DIR}" \
      "$AMYPRED_IMAGE" \
      bash -lc "set -Eeuo pipefail; if [ ! -f '${DEPS_DIR}/pandas/__init__.py' ]; then python -m pip install --target '${DEPS_DIR}' numpy scipy scikit-learn pandas xgboost tabulate; fi; PYTHONPATH='${DEPS_DIR}' python predict.py '${INPUT_FASTA}' '${output_csv}'"
}

run_amylogram() {
  local output_csv="${RESULTS_DIR}/amylogram_prediction_fast.csv"
  local summary_tsv="${RESULTS_DIR}/amylogram_fast_summary.tsv"
  local chunk_dir="${RESULTS_DIR}/amylogram_fast_chunks_${RUN_ID}"
  if [[ "$FORCE" != "1" && -s "$output_csv" && -s "$summary_tsv" ]]; then
    log "SKIP AmyloGram: existing ${output_csv}"
    record_status "AmyloGram" "SKIP" "existing output"
    return 0
  fi

  test -s "${AMYLOGRAM_DIR}/predict_amylogram_fast.R"
  test -s "$INPUT_FASTA"

  mkdir -p "$chunk_dir"
  run_step "AmyloGram" "predict fast" \
    docker run --rm \
      --name "amylogram-fast-${SAMPLE_ID}-${RUN_ID}" \
      -e AMYLOGRAM_CORES="${AMYLOGRAM_CORES}" \
      -e AMYLOGRAM_CHUNK_SIZE="${AMYLOGRAM_CHUNK_SIZE}" \
      -e AMYLOGRAM_RETRY_SPLIT_DEPTH="${AMYLOGRAM_RETRY_SPLIT_DEPTH}" \
      -e AMYLOGRAM_ALLOW_ERRORS="${AMYLOGRAM_ALLOW_ERRORS}" \
      -v "${WORK_DIR}:${WORK_DIR}" \
      -v "${BASE_DIR}/amyloid_predictors:${BASE_DIR}/amyloid_predictors:ro" \
      -v "${BASE_DIR}/amyloid_predictors_one_button_20260516T171835Z:${BASE_DIR}/amyloid_predictors_one_button_20260516T171835Z:ro" \
      -w "${WORK_DIR}" \
      "$AMYLOGRAM_IMAGE" \
      Rscript "${AMYLOGRAM_DIR}/predict_amylogram_fast.R" \
        "$INPUT_FASTA" \
        "$output_csv" \
        "$chunk_dir" \
        "$summary_tsv"
}

run_amylogram_py() {
  local output_csv="${RESULTS_DIR}/amylogram_py_prediction.csv"
  local report_md="${RESULTS_DIR}/amylogram_py_report.md"
  local report_json="${RESULTS_DIR}/amylogram_py_report.json"
  local skipped_tsv="${RESULTS_DIR}/amylogram_py_skipped.tsv"
  local top_tsv="${RESULTS_DIR}/amylogram_py_top_hits.tsv"
  local lookup_dir
  lookup_dir="$(dirname "$AMYLOGRAM_PY_LOOKUP")"

  if [[ "$FORCE" != "1" && -s "$output_csv" && -s "$report_json" ]]; then
    log "SKIP AmyloGram-Py: existing ${output_csv}"
    record_status "AmyloGram-Py" "SKIP" "existing output"
    return 0
  fi

  test -s "$INPUT_FASTA"
  if [[ ! -s "$AMYLOGRAM_PY_LOOKUP" ]]; then
    log "AmyloGram-Py lookup missing: ${AMYLOGRAM_PY_LOOKUP}"
    return 1
  fi
  if ! docker image inspect "$AMYLOGRAM_PY_IMAGE" >/dev/null 2>&1; then
    if [[ -s "$AMYLOGRAM_PY_IMAGE_TAR" ]]; then
      log "Loading AmyloGram-Py Docker image from ${AMYLOGRAM_PY_IMAGE_TAR}"
      if ! gzip -cd "$AMYLOGRAM_PY_IMAGE_TAR" | docker load >> "${LOG_DIR}/AmyloGram-Py-image-load_${RUN_ID}.log" 2>&1; then
        log "Failed to load AmyloGram-Py Docker image from ${AMYLOGRAM_PY_IMAGE_TAR}"
        return 1
      fi
    fi
  fi
  if ! docker image inspect "$AMYLOGRAM_PY_IMAGE" >/dev/null 2>&1; then
    log "AmyloGram-Py Docker image unavailable: ${AMYLOGRAM_PY_IMAGE}"
    return 1
  fi

  run_step "AmyloGram-Py" "predict lookup" \
    docker run --rm \
      --name "amylogram-py-${SAMPLE_ID}-${RUN_ID}" \
      --user "$(id -u):$(id -g)" \
      -v "${WORK_DIR}:${WORK_DIR}" \
      -v "${lookup_dir}:${lookup_dir}:ro" \
      "$AMYLOGRAM_PY_IMAGE" \
      "$INPUT_FASTA" \
      "$output_csv" \
      --sixmer-table "$AMYLOGRAM_PY_LOOKUP" \
      --report-md "$report_md" \
      --report-json "$report_json" \
      --skipped-tsv "$skipped_tsv" \
      --top-k "$AMYLOGRAM_PY_TOP_K" \
      --top-tsv "$top_tsv"
}

run_protein_features() {
  local output_csv="${RESULTS_DIR}/protein_features.csv"
  local docker_args=()
  if [[ "$FORCE" != "1" && -s "$output_csv" ]]; then
    log "SKIP ProteinFeatures: existing ${output_csv}"
    record_status "ProteinFeatures" "SKIP" "existing output"
    return 0
  fi

  test -s "${PROTEIN_FEATURES_DIR}/calc_protein_features.R"
  test -s "$INPUT_FASTA"

  if [[ -n "$PROTEIN_FEATURES_R_LIBS" ]]; then
    docker_args+=(-e "R_LIBS_USER=${PROTEIN_FEATURES_R_LIBS}" -v "${PROTEIN_FEATURES_R_LIBS}:${PROTEIN_FEATURES_R_LIBS}:ro")
  fi

  run_step "ProteinFeatures" "calculate R physicochemical features" \
    docker run --rm \
      --name "protein-features-${SAMPLE_ID}-${RUN_ID}" \
      --user "$(id -u):$(id -g)" \
      "${docker_args[@]}" \
      -v "${WORK_DIR}:${WORK_DIR}" \
      -v "${BASE_DIR}/amyloid_predictors:${BASE_DIR}/amyloid_predictors:ro" \
      -v "${PROTEIN_FEATURES_DIR}:${PROTEIN_FEATURES_DIR}:ro" \
      -w "${WORK_DIR}" \
      "$PROTEIN_FEATURES_IMAGE" \
      Rscript "${PROTEIN_FEATURES_DIR}/calc_protein_features.R" \
        "$INPUT_FASTA" \
        "$output_csv" \
        "$PROTEIN_FEATURES_WINDOW" \
        "$PROTEIN_FEATURES_PH"
}

merge_predictions() {
  local amypred_csv="${RESULTS_DIR}/amypred_frl_prediction.csv"
  local amylogram_csv="${RESULTS_DIR}/amylogram_prediction_fast.csv"
  local amylogram_py_csv="${RESULTS_DIR}/amylogram_py_prediction.csv"
  if [[ ! -s "$amypred_csv" && ! -s "$amylogram_csv" && ! -s "$amylogram_py_csv" ]]; then
    log "Combined table skipped: no predictor outputs are available"
    record_status "combined" "SKIP" "missing predictor output"
    return 0
  fi

  run_step "combined" "merge predictor outputs" \
    python3 - "$amypred_csv" "$amylogram_csv" "$amylogram_py_csv" "$COMBINED_CSV" <<'PY'
import csv
import sys
from pathlib import Path

amypred_csv, amylogram_csv, amylogram_py_csv, output_csv = sys.argv[1:5]

def norm_pred(value):
    value = (value or "").strip().lower().replace("_", "-")
    if value in {"amyloid", "amyloidogenic", "amyl"}:
        return "Amyloid"
    if value in {"non-amyloid", "non amyloid", "nonamyloid", "not-amyloid"}:
        return "Non-Amyloid"
    return value or "NA"

amypred = {}
if Path(amypred_csv).is_file() and Path(amypred_csv).stat().st_size:
    with open(amypred_csv, newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            sid = row.get("Sequence ID") or row.get("Sequence_ID")
            if not sid:
                continue
            amypred[sid] = {
                "AMYPred_Prob": row.get("Amyloid Prob", ""),
                "AMYPred_Pred": norm_pred(row.get("Prediction", "")),
            }

amylogram = {}
if Path(amylogram_csv).is_file() and Path(amylogram_csv).stat().st_size:
    with open(amylogram_csv, newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            sid = row.get("Sequence_ID") or row.get("Sequence ID")
            if not sid:
                continue
            amylogram[sid] = {
                "AmyloGram_Prob": row.get("AmyloGram_Prob", ""),
                "AmyloGram_Pred": norm_pred(row.get("AmyloGram_Pred", "")),
            }

amylogram_py = {}
if Path(amylogram_py_csv).is_file() and Path(amylogram_py_csv).stat().st_size:
    with open(amylogram_py_csv, newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            sid = row.get("Sequence_ID") or row.get("Sequence ID")
            if not sid:
                continue
            amylogram_py[sid] = {
                "AmyloGramPy_Prob": row.get("AmyloGram_Prob", ""),
                "AmyloGramPy_Pred": norm_pred(row.get("AmyloGram_Pred", "")),
            }

ids = sorted(set(amypred) | set(amylogram) | set(amylogram_py))
fields = [
    "Sequence_ID",
    "AMYPred_Prob",
    "AMYPred_Pred",
    "AmyloGram_Prob",
    "AmyloGram_Pred",
    "AmyloGramPy_Prob",
    "AmyloGramPy_Pred",
    "Consensus",
]

with open(output_csv, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    for sid in ids:
        a = amypred.get(sid, {"AMYPred_Prob": "", "AMYPred_Pred": "NA"})
        g = amylogram.get(sid, {"AmyloGram_Prob": "", "AmyloGram_Pred": "NA"})
        gp = amylogram_py.get(sid, {"AmyloGramPy_Prob": "", "AmyloGramPy_Pred": "NA"})
        preds = [a["AMYPred_Pred"], g["AmyloGram_Pred"], gp["AmyloGramPy_Pred"]]
        available = [pred for pred in preds if pred != "NA"]
        if not available:
            consensus = "NA"
        elif len(set(available)) == 1 and len(available) == len(preds):
            consensus = available[0]
        elif len(set(available)) == 1:
            consensus = "Partial"
        else:
            consensus = "Discordant"
        writer.writerow({
            "Sequence_ID": sid,
            "AMYPred_Prob": a["AMYPred_Prob"],
            "AMYPred_Pred": a["AMYPred_Pred"],
            "AmyloGram_Prob": g["AmyloGram_Prob"],
            "AmyloGram_Pred": g["AmyloGram_Pred"],
            "AmyloGramPy_Prob": gp["AmyloGramPy_Prob"],
            "AmyloGramPy_Pred": gp["AmyloGramPy_Pred"],
            "Consensus": consensus,
        })
PY
}

write_summary() {
  local input_records=0
  if [[ -s "$INPUT_FASTA" ]]; then
    input_records="$(grep -c '^>' "$INPUT_FASTA" || true)"
  fi

  {
    printf 'metric\tvalue\n'
    printf 'sample\t%s\n' "$SAMPLE_ID"
    printf 'run_id\t%s\n' "$RUN_ID"
    printf 'input_fasta\t%s\n' "$INPUT_FASTA"
    printf 'input_records\t%s\n' "$input_records"
    printf 'amypred_rows\t%s\n' "$(line_count_without_header "${RESULTS_DIR}/amypred_frl_prediction.csv")"
    printf 'amylogram_rows\t%s\n' "$(line_count_without_header "${RESULTS_DIR}/amylogram_prediction_fast.csv")"
    printf 'amylogram_py_rows\t%s\n' "$(line_count_without_header "${RESULTS_DIR}/amylogram_py_prediction.csv")"
    printf 'protein_features_rows\t%s\n' "$(line_count_without_header "${RESULTS_DIR}/protein_features.csv")"
    printf 'combined_rows\t%s\n' "$(line_count_without_header "$COMBINED_CSV")"
    printf 'status_file\t%s\n' "$STATUS_TSV"
    printf 'main_log\t%s\n' "$MAIN_LOG"
    printf 'output_s3\t%s\n' "$OUTPUT_S3"
  } > "$SUMMARY_TSV"
}

upload_if_exists() {
  local src="$1"
  local dst="$2"
  if [[ -s "$src" ]]; then
    aws s3 cp "$src" "$dst" --profile "$AWS_PROFILE" --region "$AWS_REGION" --quiet
  else
    log "Optional upload skipped, missing file: ${src}"
  fi
}

upload_outputs() {
  if [[ "$UPLOAD" != "1" ]]; then
    log "Upload disabled"
    return 0
  fi
  log "Uploading predictor outputs to ${OUTPUT_S3}"
  upload_if_exists "${RESULTS_DIR}/amypred_frl_prediction.csv" "${OUTPUT_S3}/amypred_frl_prediction.csv"
  upload_if_exists "${RESULTS_DIR}/amylogram_prediction_fast.csv" "${OUTPUT_S3}/amylogram_prediction_fast.csv"
  upload_if_exists "${RESULTS_DIR}/amylogram_fast_summary.tsv" "${OUTPUT_S3}/amylogram_fast_summary.tsv"
  upload_if_exists "${RESULTS_DIR}/amylogram_py_prediction.csv" "${OUTPUT_S3}/amylogram_py_prediction.csv"
  upload_if_exists "${RESULTS_DIR}/amylogram_py_report.md" "${OUTPUT_S3}/amylogram_py_report.md"
  upload_if_exists "${RESULTS_DIR}/amylogram_py_report.json" "${OUTPUT_S3}/amylogram_py_report.json"
  upload_if_exists "${RESULTS_DIR}/amylogram_py_skipped.tsv" "${OUTPUT_S3}/amylogram_py_skipped.tsv"
  upload_if_exists "${RESULTS_DIR}/amylogram_py_top_hits.tsv" "${OUTPUT_S3}/amylogram_py_top_hits.tsv"
  upload_if_exists "${RESULTS_DIR}/protein_features.csv" "${OUTPUT_S3}/protein_features.csv"
  if [[ -n "$PROTEIN_FEATURES_OUTPUT_S3" ]]; then
    upload_if_exists "${RESULTS_DIR}/protein_features.csv" "${PROTEIN_FEATURES_OUTPUT_S3}/protein_features.csv"
    upload_if_exists "$STATUS_TSV" "${PROTEIN_FEATURES_OUTPUT_S3}/protein_features_status.tsv"
    upload_if_exists "$SUMMARY_TSV" "${PROTEIN_FEATURES_OUTPUT_S3}/protein_features_summary.tsv"
    upload_if_exists "$MAIN_LOG" "${PROTEIN_FEATURES_OUTPUT_S3}/logs/amyloid_predictors_${RUN_ID}.log"
    upload_if_exists "${LOG_DIR}/ProteinFeatures_${RUN_ID}.log" "${PROTEIN_FEATURES_OUTPUT_S3}/logs/ProteinFeatures_${RUN_ID}.log"
  fi
  upload_if_exists "${RESULTS_DIR}/amyloid_combined_predictions.csv" "${OUTPUT_S3}/amyloid_combined_predictions.csv"
  aws s3 cp "$STATUS_TSV" "${OUTPUT_S3}/amyloid_predictors_status.tsv" --profile "$AWS_PROFILE" --region "$AWS_REGION" --quiet
  aws s3 cp "$SUMMARY_TSV" "${OUTPUT_S3}/amyloid_predictors_summary.tsv" --profile "$AWS_PROFILE" --region "$AWS_REGION" --quiet
  aws s3 cp "$MAIN_LOG" "${OUTPUT_S3}/logs/amyloid_predictors_${RUN_ID}.log" --profile "$AWS_PROFILE" --region "$AWS_REGION" --quiet
}

main() {
  printf 'timestamp\tpredictor\tstatus\tmessage\n' > "$STATUS_TSV"
  log "sample=${SAMPLE_ID}"
  log "work_dir=${WORK_DIR}"
  log "input_fasta=${INPUT_FASTA}"
  log "output_s3=${OUTPUT_S3}"

  ensure_input

  if [[ "$RUN_AMYPRED" == "1" || "$RUN_AMYPRED" == "true" ]]; then
    run_amypred || true
  else
    log "AMYPred-FRL skipped by option"
    record_status "AMYPred-FRL" "SKIP" "disabled"
  fi

  if [[ "$RUN_AMYLOGRAM" == "1" || "$RUN_AMYLOGRAM" == "true" ]]; then
    run_amylogram || true
  else
    log "Legacy R AmyloGram skipped by option/default"
    record_status "AmyloGram" "SKIP" "disabled; enable with --run-amylogram or RUN_AMYLOGRAM=1"
  fi

  if [[ "$RUN_AMYLOGRAM_PY" == "1" || "$RUN_AMYLOGRAM_PY" == "true" ]]; then
    run_amylogram_py || true
  else
    log "Primary AmyloGram-Py skipped by option"
    record_status "AmyloGram-Py" "SKIP" "disabled"
  fi

  if [[ "$RUN_PROTEIN_FEATURES" == "1" || "$RUN_PROTEIN_FEATURES" == "true" ]]; then
    run_protein_features || true
  else
    log "ProteinFeatures skipped by option"
    record_status "ProteinFeatures" "SKIP" "disabled"
  fi

  merge_predictions || true
  write_summary

  if grep -q $'\tFAIL\t' "$STATUS_TSV"; then
    log "One or more steps failed; uploading diagnostics"
    failed=1
  else
    touch "${RESULTS_DIR}/AMYLOID_PREDICTORS_SUCCESS"
    log "All enabled amyloid predictor steps completed or were safely skipped"
    failed=0
  fi

  upload_outputs
  log "Done"
  return "$failed"
}

main "$@"
