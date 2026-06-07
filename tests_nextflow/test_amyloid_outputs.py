from __future__ import annotations

import csv
from pathlib import Path

import pytest

from tests_nextflow.validate_amyloid_outputs import (
    AMYPRED_FRL_COLUMNS,
    COMBINED_COLUMNS,
    PROTEIN_FEATURES_COLUMNS,
    validate_outputs,
)


def write_csv(path: Path, columns: list[str], row: dict[str, str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerow(row)


def test_combined_schema_contract() -> None:
    assert COMBINED_COLUMNS == [
        "Sequence_ID",
        "AMYPred_Prob",
        "AMYPred_Pred",
        "AmyloGramPy_Prob",
        "AmyloGramPy_Pred",
        "Consensus",
    ]


def test_optional_predictor_schema_contracts() -> None:
    assert AMYPRED_FRL_COLUMNS == ["Sequence ID", "Amyloid Prob", "Prediction"]
    assert PROTEIN_FEATURES_COLUMNS == [
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


def test_validate_amyloid_outputs_accepts_expected_schema(tmp_path: Path) -> None:
    write_csv(
        tmp_path / "amylogram_py_prediction.csv",
        ["Sequence_ID", "AmyloGram_Prob", "AmyloGram_Pred"],
        {"Sequence_ID": "p1", "AmyloGram_Prob": "0.7", "AmyloGram_Pred": "AMYLOID"},
    )
    write_csv(
        tmp_path / "amyloid_combined_predictions.csv",
        COMBINED_COLUMNS,
        {
            "Sequence_ID": "p1",
            "AMYPred_Prob": "",
            "AMYPred_Pred": "NA",
            "AmyloGramPy_Prob": "0.7",
            "AmyloGramPy_Pred": "Amyloid",
            "Consensus": "Partial",
        },
    )
    (tmp_path / "summary.tsv").write_text("metric\tvalue\ncombined_rows\t1\n", encoding="utf-8")
    (tmp_path / "status.tsv").write_text("timestamp\tpredictor\tstatus\tmessage\nNA\tcombined\tOK\tok\n", encoding="utf-8")

    validate_outputs(tmp_path)


def test_validate_amyloid_outputs_rejects_bad_combined_schema(tmp_path: Path) -> None:
    write_csv(
        tmp_path / "amylogram_py_prediction.csv",
        ["Sequence_ID", "AmyloGram_Prob", "AmyloGram_Pred"],
        {"Sequence_ID": "p1", "AmyloGram_Prob": "0.7", "AmyloGram_Pred": "AMYLOID"},
    )
    write_csv(
        tmp_path / "amyloid_combined_predictions.csv",
        ["Sequence_ID", "AmyloGramPy_Prob"],
        {"Sequence_ID": "p1", "AmyloGramPy_Prob": "0.7"},
    )
    (tmp_path / "summary.tsv").write_text("metric\tvalue\ncombined_rows\t1\n", encoding="utf-8")
    (tmp_path / "status.tsv").write_text("timestamp\tpredictor\tstatus\tmessage\nNA\tcombined\tOK\tok\n", encoding="utf-8")

    with pytest.raises(AssertionError):
        validate_outputs(tmp_path)
