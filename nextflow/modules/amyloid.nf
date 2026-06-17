process AMYLOGRAM_PY {
    tag "$sample_id"
    container 'amylogram-py-nextflow:local'

    input:
    tuple val(sample_id), path(input_fasta), path(sixmer_table)

    output:
    tuple val(sample_id), path(input_fasta), path('amylogram_py_prediction.csv'), path('amylogram_py_report.json'), path('amylogram_py_report.md'), path('amylogram_py_skipped.tsv'), path('amylogram_py_top_hits.tsv')

    script:
    """
    set -Eeuo pipefail
    amylogram-py "${input_fasta}" amylogram_py_prediction.csv \
      --sixmer-table "${sixmer_table}" \
      --report-json amylogram_py_report.json \
      --report-md amylogram_py_report.md \
      --skipped-tsv amylogram_py_skipped.tsv \
      --top-k "${params.amylogram_py_top_k}" \
      --top-tsv amylogram_py_top_hits.tsv

    python3 - <<'PY'
import csv
with open("amylogram_py_prediction.csv", newline="", encoding="utf-8") as handle:
    fields = next(csv.reader(handle))
expected = ["Sequence_ID", "AmyloGram_Prob", "AmyloGram_Pred"]
if fields != expected:
    raise SystemExit(f"Unexpected AmyloGram-Py columns: {fields}")
PY
    """
}

process AMYPRED_FRL {
    tag "$sample_id"
    container "${params.amypred_container}"

    input:
    tuple val(sample_id), path(input_fasta), path(amypred_dir)

    output:
    tuple val(sample_id), path(input_fasta), path('amypred_frl_prediction.csv')

    script:
    """
    set -Eeuo pipefail
    test -s "${amypred_dir}/predict.py"
    test -s "${amypred_dir}/data/PAAC.txt"
    test -s "${amypred_dir}/data/TR_P_132.fasta"
    test -s "${amypred_dir}/data/TR_N_305.fasta"
    test -s "${amypred_dir}/model/pima.pickle_model_svm_PF.dat"
    cp "${input_fasta}" input.fasta
    work_input="\$PWD/input.fasta"
    work_output="\$PWD/amypred_frl_prediction.csv"
    cd "${amypred_dir}"
    python predict.py "\${work_input}" "\${work_output}"
    cd "\$OLDPWD"
    test -s amypred_frl_prediction.csv
    """
}

process AMYLOGRAM_R {
    tag "$sample_id"
    container "${params.amylogram_r_container}"

    input:
    tuple val(sample_id), path(input_fasta), path(amylogram_dir)

    output:
    tuple val(sample_id), path('amylogram_prediction_fast.csv'), path('amylogram_fast_summary.tsv')

    script:
    """
    set -Eeuo pipefail
    mkdir -p amylogram_fast_chunks
    Rscript "${amylogram_dir}/predict_amylogram_fast.R" \
      "${input_fasta}" \
      amylogram_prediction_fast.csv \
      amylogram_fast_chunks \
      amylogram_fast_summary.tsv
    test -s amylogram_prediction_fast.csv
    test -s amylogram_fast_summary.tsv
    """
}

process MERGE_AMYLOID_PREDICTIONS {
    tag "$sample_id"
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${sample_id}/${filename}" }

    input:
    tuple val(sample_id), path(input_fasta), path(amylogram_py_csv), path(amylogram_report_json), path(amylogram_report_md), path(amylogram_skipped_tsv), path(amylogram_top_tsv)

    output:
    tuple val(sample_id), path('results_amyloid')

    script:
    """
    set -Eeuo pipefail
    mkdir -p results_amyloid
    cp "${amylogram_py_csv}" results_amyloid/amylogram_py_prediction.csv
    cp "${amylogram_report_json}" results_amyloid/amylogram_py_report.json
    cp "${amylogram_report_md}" results_amyloid/amylogram_py_report.md
    cp "${amylogram_skipped_tsv}" results_amyloid/amylogram_py_skipped.tsv
    cp "${amylogram_top_tsv}" results_amyloid/amylogram_py_top_hits.tsv

    python3 - "${input_fasta}" results_amyloid/amylogram_py_prediction.csv results_amyloid/amylogram_py_report.json <<'PY'
import csv
import json
import sys
from pathlib import Path

input_fasta, amylogram_py_csv, report_json = sys.argv[1:4]

def norm_pred(value):
    value = (value or "").strip().lower().replace("_", "-")
    if value in {"amyloid", "amyloidogenic", "amyl"}:
        return "Amyloid"
    if value in {"non-amyloid", "non amyloid", "nonamyloid", "not-amyloid"}:
        return "Non-Amyloid"
    return value or "NA"

rows = []
with open(amylogram_py_csv, newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    if reader.fieldnames != ["Sequence_ID", "AmyloGram_Prob", "AmyloGram_Pred"]:
        raise SystemExit(f"Unexpected AmyloGram-Py columns: {reader.fieldnames}")
    for row in reader:
        pred = norm_pred(row.get("AmyloGram_Pred"))
        rows.append({
            "Sequence_ID": row.get("Sequence_ID", ""),
            "AMYPred_Prob": "",
            "AMYPred_Pred": "NA",
            "AmyloGramPy_Prob": row.get("AmyloGram_Prob", ""),
            "AmyloGramPy_Pred": pred,
            "Consensus": "Partial" if pred != "NA" else "NA",
        })

fields = ["Sequence_ID", "AMYPred_Prob", "AMYPred_Pred", "AmyloGramPy_Prob", "AmyloGramPy_Pred", "Consensus"]
with open("results_amyloid/amyloid_combined_predictions.csv", "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

metrics = json.loads(Path(report_json).read_text(encoding="utf-8"))
with open("results_amyloid/amyloid_predictors_summary.tsv", "w", encoding="utf-8") as handle:
    handle.write("metric\\tvalue\\n")
    handle.write(f"input_fasta\\t{input_fasta}\\n")
    handle.write(f"input_records\\t{metrics.get('total_records', 0)}\\n")
    handle.write(f"amylogram_py_rows\\t{len(rows)}\\n")
    handle.write(f"combined_rows\\t{len(rows)}\\n")

with open("results_amyloid/amyloid_predictors_status.tsv", "w", encoding="utf-8") as handle:
    handle.write("timestamp\\tpredictor\\tstatus\\tmessage\\n")
    handle.write("NA\\tAMYPred-FRL\\tSKIP\\tpending Nextflow process integration\\n")
    handle.write("NA\\tAmyloGram-Py\\tOK\\tpredicted records\\n")
    handle.write("NA\\tProteinFeaturesFull\\tSKIP\\tpending Nextflow process integration\\n")
    handle.write("NA\\tcombined\\tOK\\tmerged available predictor outputs\\n")

if not rows:
    raise SystemExit("No combined rows were produced")
PY

    test -s results_amyloid/amylogram_py_prediction.csv
    test -s results_amyloid/amyloid_combined_predictions.csv
    test -s results_amyloid/amyloid_predictors_summary.tsv
    test -s results_amyloid/amyloid_predictors_status.tsv
    ! grep -q \$'\\tFAIL\\t' results_amyloid/amyloid_predictors_status.tsv
    """
}

process MERGE_AMYLOID_PREDICTIONS_AMYPRED_PY {
    tag "$sample_id"
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${sample_id}/${filename}" }

    input:
    tuple val(sample_id), path(input_fasta), path(amypred_csv), path(amylogram_py_csv), path(amylogram_report_json), path(amylogram_report_md), path(amylogram_skipped_tsv), path(amylogram_top_tsv)

    output:
    tuple val(sample_id), path('results_amyloid')

    script:
    """
    set -Eeuo pipefail
    mkdir -p results_amyloid
    cp "${amypred_csv}" results_amyloid/amypred_frl_prediction.csv
    cp "${amylogram_py_csv}" results_amyloid/amylogram_py_prediction.csv
    cp "${amylogram_report_json}" results_amyloid/amylogram_py_report.json
    cp "${amylogram_report_md}" results_amyloid/amylogram_py_report.md
    cp "${amylogram_skipped_tsv}" results_amyloid/amylogram_py_skipped.tsv
    cp "${amylogram_top_tsv}" results_amyloid/amylogram_py_top_hits.tsv

    python3 - "${input_fasta}" results_amyloid/amypred_frl_prediction.csv results_amyloid/amylogram_py_prediction.csv results_amyloid/amylogram_py_report.json <<'PY'
import csv
import json
import sys
from pathlib import Path

input_fasta, amypred_csv, amylogram_py_csv, report_json = sys.argv[1:5]

def norm_pred(value):
    value = (value or "").strip().lower().replace("_", "-")
    if value in {"amyloid", "amyloidogenic", "amyl"}:
        return "Amyloid"
    if value in {"non-amyloid", "non amyloid", "nonamyloid", "not-amyloid"}:
        return "Non-Amyloid"
    return value or "NA"

def count_rows(path):
    if not Path(path).is_file() or not Path(path).stat().st_size:
        return 0
    with open(path, encoding="utf-8-sig") as handle:
        return max(sum(1 for _ in handle) - 1, 0)

amypred = {}
with open(amypred_csv, newline="", encoding="utf-8-sig") as handle:
    reader = csv.DictReader(handle)
    if reader.fieldnames != ["Sequence ID", "Amyloid Prob", "Prediction"]:
        raise SystemExit(f"Unexpected AMYPred-FRL columns: {reader.fieldnames}")
    for row in reader:
        sid = row.get("Sequence ID", "")
        if sid:
            amypred[sid] = {
                "AMYPred_Prob": row.get("Amyloid Prob", ""),
                "AMYPred_Pred": norm_pred(row.get("Prediction")),
            }

amylogram_py = {}
with open(amylogram_py_csv, newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    if reader.fieldnames != ["Sequence_ID", "AmyloGram_Prob", "AmyloGram_Pred"]:
        raise SystemExit(f"Unexpected AmyloGram-Py columns: {reader.fieldnames}")
    for row in reader:
        sid = row.get("Sequence_ID", "")
        if sid:
            amylogram_py[sid] = {
                "AmyloGramPy_Prob": row.get("AmyloGram_Prob", ""),
                "AmyloGramPy_Pred": norm_pred(row.get("AmyloGram_Pred")),
            }

ids = sorted(set(amypred) | set(amylogram_py))
fields = ["Sequence_ID", "AMYPred_Prob", "AMYPred_Pred", "AmyloGramPy_Prob", "AmyloGramPy_Pred", "Consensus"]
with open("results_amyloid/amyloid_combined_predictions.csv", "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    for sid in ids:
        a = amypred.get(sid, {"AMYPred_Prob": "", "AMYPred_Pred": "NA"})
        gp = amylogram_py.get(sid, {"AmyloGramPy_Prob": "", "AmyloGramPy_Pred": "NA"})
        preds = [a["AMYPred_Pred"], gp["AmyloGramPy_Pred"]]
        available = [pred for pred in preds if pred != "NA"]
        if not available:
            consensus = "NA"
        elif len(available) == 2 and len(set(available)) == 1:
            consensus = available[0]
        elif len(available) == 2:
            consensus = "Discordant"
        else:
            consensus = "Partial"
        writer.writerow({
            "Sequence_ID": sid,
            "AMYPred_Prob": a["AMYPred_Prob"],
            "AMYPred_Pred": a["AMYPred_Pred"],
            "AmyloGramPy_Prob": gp["AmyloGramPy_Prob"],
            "AmyloGramPy_Pred": gp["AmyloGramPy_Pred"],
            "Consensus": consensus,
        })

metrics = json.loads(Path(report_json).read_text(encoding="utf-8"))
with open("results_amyloid/amyloid_predictors_summary.tsv", "w", encoding="utf-8") as handle:
    handle.write("metric\\tvalue\\n")
    handle.write(f"input_fasta\\t{input_fasta}\\n")
    handle.write(f"input_records\\t{metrics.get('total_records', 0)}\\n")
    handle.write(f"amypred_rows\\t{count_rows(amypred_csv)}\\n")
    handle.write(f"amylogram_py_rows\\t{count_rows(amylogram_py_csv)}\\n")
    handle.write(f"combined_rows\\t{len(ids)}\\n")

with open("results_amyloid/amyloid_predictors_status.tsv", "w", encoding="utf-8") as handle:
    handle.write("timestamp\\tpredictor\\tstatus\\tmessage\\n")
    handle.write("NA\\tAMYPred-FRL\\tOK\\tpredicted records\\n")
    handle.write("NA\\tAmyloGram-Py\\tOK\\tpredicted records\\n")
    handle.write("NA\\tProteinFeaturesFull\\tSKIP\\tpending Nextflow process integration\\n")
    handle.write("NA\\tcombined\\tOK\\tmerged predictor outputs\\n")

if not ids:
    raise SystemExit("No combined rows were produced")
PY

    test -s results_amyloid/amypred_frl_prediction.csv
    test -s results_amyloid/amylogram_py_prediction.csv
    test -s results_amyloid/amyloid_combined_predictions.csv
    test -s results_amyloid/amyloid_predictors_summary.tsv
    test -s results_amyloid/amyloid_predictors_status.tsv
    ! grep -q \$'\\tFAIL\\t' results_amyloid/amyloid_predictors_status.tsv
    """
}

process MERGE_AMYLOID_PREDICTIONS_FULL {
    tag "$sample_id"
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${sample_id}/${filename}" }

    input:
    tuple val(sample_id), path(input_fasta), path(amypred_csv), path(amylogram_r_csv), path(amylogram_r_summary), path(amylogram_py_csv), path(amylogram_report_json), path(amylogram_report_md), path(amylogram_skipped_tsv), path(amylogram_top_tsv)

    output:
    tuple val(sample_id), path('results_amyloid')

    script:
    """
    set -Eeuo pipefail
    mkdir -p results_amyloid
    cp "${amypred_csv}" results_amyloid/amypred_frl_prediction.csv
    cp "${amylogram_r_csv}" results_amyloid/amylogram_prediction_fast.csv
    cp "${amylogram_r_summary}" results_amyloid/amylogram_fast_summary.tsv
    cp "${amylogram_py_csv}" results_amyloid/amylogram_py_prediction.csv
    cp "${amylogram_report_json}" results_amyloid/amylogram_py_report.json
    cp "${amylogram_report_md}" results_amyloid/amylogram_py_report.md
    cp "${amylogram_skipped_tsv}" results_amyloid/amylogram_py_skipped.tsv
    cp "${amylogram_top_tsv}" results_amyloid/amylogram_py_top_hits.tsv

    python3 - "${input_fasta}" "${amypred_csv}" "${amylogram_r_csv}" "${amylogram_py_csv}" <<'PY'
import csv
import sys
from pathlib import Path

input_fasta, amypred_csv, amylogram_csv, amylogram_py_csv = sys.argv[1:5]

def norm_pred(value):
    value = (value or "").strip().lower().replace("_", "-")
    if value in {"amyloid", "amyloidogenic", "amyl"}:
        return "Amyloid"
    if value in {"non-amyloid", "non amyloid", "nonamyloid", "not-amyloid"}:
        return "Non-Amyloid"
    return value or "NA"

def read_amypred(path):
    rows = {}
    if Path(path).is_file() and Path(path).stat().st_size:
        with open(path, newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                sid = row.get("Sequence ID") or row.get("Sequence_ID")
                if sid:
                    rows[sid] = {
                        "AMYPred_Prob": row.get("Amyloid Prob", ""),
                        "AMYPred_Pred": norm_pred(row.get("Prediction", "")),
                    }
    return rows

def read_amylogram(path, prob_col, pred_col, out_prob, out_pred):
    rows = {}
    if Path(path).is_file() and Path(path).stat().st_size:
        with open(path, newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                sid = row.get("Sequence_ID") or row.get("Sequence ID")
                if sid:
                    rows[sid] = {
                        out_prob: row.get(prob_col, ""),
                        out_pred: norm_pred(row.get(pred_col, "")),
                    }
    return rows

amypred = read_amypred(amypred_csv)
amylogram = read_amylogram(amylogram_csv, "AmyloGram_Prob", "AmyloGram_Pred", "AmyloGram_Prob", "AmyloGram_Pred")
amylogram_py = read_amylogram(amylogram_py_csv, "AmyloGram_Prob", "AmyloGram_Pred", "AmyloGramPy_Prob", "AmyloGramPy_Pred")

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
with open("results_amyloid/amyloid_combined_predictions.csv", "w", newline="", encoding="utf-8") as handle:
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

def count_rows(path):
    if not Path(path).is_file() or not Path(path).stat().st_size:
        return 0
    with open(path, encoding="utf-8-sig") as handle:
        return max(sum(1 for _ in handle) - 1, 0)

input_records = 0
with open(input_fasta, encoding="utf-8") as handle:
    input_records = sum(1 for line in handle if line.startswith(">"))

with open("results_amyloid/amyloid_predictors_summary.tsv", "w", encoding="utf-8") as handle:
    handle.write("metric\\tvalue\\n")
    handle.write(f"input_fasta\\t{input_fasta}\\n")
    handle.write(f"input_records\\t{input_records}\\n")
    handle.write(f"amypred_rows\\t{count_rows(amypred_csv)}\\n")
    handle.write(f"amylogram_rows\\t{count_rows(amylogram_csv)}\\n")
    handle.write(f"amylogram_py_rows\\t{count_rows(amylogram_py_csv)}\\n")
    handle.write(f"combined_rows\\t{len(ids)}\\n")
PY

    cat > results_amyloid/amyloid_predictors_status.tsv <<'EOF'
timestamp	predictor	status	message
NA	AMYPred-FRL	OK	predicted records
NA	AmyloGram	OK	predicted records
NA	AmyloGram-Py	OK	predicted records
NA	combined	OK	merged predictor outputs
EOF
    test -s results_amyloid/amyloid_combined_predictions.csv
    test -s results_amyloid/amyloid_predictors_summary.tsv
    test -s results_amyloid/amyloid_predictors_status.tsv
    """
}

process PROTEIN_FEATURES_LIGHT {
    tag "$sample_id"
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${sample_id}/${filename}" }

    input:
    tuple val(sample_id), path(input_fasta)

    output:
    tuple val(sample_id), path('results_protein_features')

    script:
    """
    set -Eeuo pipefail
    mkdir -p results_protein_features

    python3 - "${input_fasta}" <<'PY'
import csv
import math
import sys

fasta = sys.argv[1]

kd = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}

charged = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1}
aromatic = set("FWY")
alpha = set("ALEMQKRH")
beta = set("VIYFWTC")

def read_fasta(path):
    name = None
    parts = []
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts)

fields = [
    "protein_id",
    "protein_length",
    "KD_mean",
    "KD_max",
    "KD_sd",
    "charge",
    "charge_density",
    "aromaticity",
    "alpha_propensity",
    "beta_propensity",
    "chameleon_score",
    "pI",
]
rows = []
for protein_id, seq in read_fasta(fasta):
    clean = "".join(aa for aa in seq.upper() if aa in kd)
    length = len(clean)
    if not length:
        continue
    kd_values = [kd[aa] for aa in clean]
    kd_mean = sum(kd_values) / length
    kd_sd = math.sqrt(sum((value - kd_mean) ** 2 for value in kd_values) / length)
    charge = sum(charged.get(aa, 0.0) for aa in clean)
    rows.append({
        "protein_id": protein_id,
        "protein_length": str(length),
        "KD_mean": f"{kd_mean:.6f}",
        "KD_max": f"{max(kd_values):.6f}",
        "KD_sd": f"{kd_sd:.6f}",
        "charge": f"{charge:.6f}",
        "charge_density": f"{charge / length:.6f}",
        "aromaticity": f"{sum(aa in aromatic for aa in clean) / length:.6f}",
        "alpha_propensity": f"{sum(aa in alpha for aa in clean) / length:.6f}",
        "beta_propensity": f"{sum(aa in beta for aa in clean) / length:.6f}",
        "chameleon_score": f"{abs(kd_mean) + kd_sd:.6f}",
        "pI": "",
    })

with open("results_protein_features/protein_features_light.csv", "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

with open("results_protein_features/protein_features_light_summary.tsv", "w", encoding="utf-8") as handle:
    handle.write("metric\\tvalue\\n")
    handle.write(f"input_fasta\\t{fasta}\\n")
    handle.write(f"feature_rows\\t{len(rows)}\\n")

with open("results_protein_features/protein_features_light_status.tsv", "w", encoding="utf-8") as handle:
    handle.write("timestamp\\tcomponent\\tstatus\\tmessage\\n")
    handle.write("NA\\tProteinFeaturesLight\\tOK\\tcomputed lightweight sequence features\\n")

if not rows:
    raise SystemExit("No protein feature rows were produced")
PY

    test -s results_protein_features/protein_features_light.csv
    test -s results_protein_features/protein_features_light_summary.tsv
    test -s results_protein_features/protein_features_light_status.tsv
    ! grep -q \$'\\tFAIL\\t' results_protein_features/protein_features_light_status.tsv
    """
}

process PROTEIN_FEATURES_FULL {
    tag "$sample_id"
    container "${params.protein_features_container}"
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${sample_id}/${filename}" }

    input:
    tuple val(sample_id), path(input_fasta), path(protein_features_dir)

    output:
    tuple val(sample_id), path('results_protein_features')

    script:
    """
    set -Eeuo pipefail
    mkdir -p results_protein_features
    Rscript "${protein_features_dir}/calc_protein_features.R" \
      "${input_fasta}" \
      results_protein_features/protein_features.csv \
      "${params.protein_features_window}" \
      "${params.protein_features_ph}"
    cat > results_protein_features/protein_features_status.tsv <<'EOF'
timestamp	component	status	message
NA	ProteinFeatures	OK	computed full protein features
EOF
    {
      printf 'metric\\tvalue\\n'
      printf 'input_fasta\\t%s\\n' "${input_fasta}"
      printf 'feature_rows\\t%s\\n' "\$(tail -n +2 results_protein_features/protein_features.csv | grep -c . || true)"
    } > results_protein_features/protein_features_summary.tsv
    test -s results_protein_features/protein_features.csv
    test -s results_protein_features/protein_features_status.tsv
    test -s results_protein_features/protein_features_summary.tsv
    """
}
