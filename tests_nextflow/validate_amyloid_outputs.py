#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path
import sys


AMYLOGRAM_PY_COLUMNS = ["Sequence_ID", "AmyloGram_Prob", "AmyloGram_Pred"]
AMYPRED_FRL_COLUMNS = ["Sequence ID", "Amyloid Prob", "Prediction"]
PROTEIN_FEATURES_LIGHT_COLUMNS = [
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
COMBINED_COLUMNS = [
    "Sequence_ID",
    "AMYPred_Prob",
    "AMYPred_Pred",
    "AmyloGramPy_Prob",
    "AmyloGramPy_Pred",
    "Consensus",
]
COMBINED_FULL_COLUMNS = [
    "Sequence_ID",
    "AMYPred_Prob",
    "AMYPred_Pred",
    "AmyloGram_Prob",
    "AmyloGram_Pred",
    "AmyloGramPy_Prob",
    "AmyloGramPy_Pred",
    "Consensus",
]


def read_columns(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle).fieldnames or [])


def require_columns(path: Path, expected: list[str]) -> None:
    actual = read_columns(path)
    if actual != expected:
        raise AssertionError(f"{path.name}: expected columns {expected}, got {actual}")


def require_any_columns(path: Path, expected_options: list[list[str]]) -> None:
    actual = read_columns(path)
    if actual not in expected_options:
        raise AssertionError(f"{path.name}: expected one of {expected_options}, got {actual}")


def validate_outputs(root: Path) -> None:
    amypred = root / "amypred_frl_prediction.csv"
    if amypred.exists():
        require_columns(amypred, AMYPRED_FRL_COLUMNS)

    require_columns(root / "amylogram_py_prediction.csv", AMYLOGRAM_PY_COLUMNS)

    protein_features = root / "protein_features_light.csv"
    if protein_features.exists():
        require_columns(protein_features, PROTEIN_FEATURES_LIGHT_COLUMNS)

    require_any_columns(root / "amyloid_combined_predictions.csv", [COMBINED_COLUMNS, COMBINED_FULL_COLUMNS])

    summary = {}
    with (root / "amyloid_predictors_summary.tsv").open(encoding="utf-8") as handle:
        next(handle)
        for raw in handle:
            key, value = raw.rstrip("\n").split("\t", 1)
            summary[key] = value
    if int(summary.get("combined_rows", "0")) <= 0:
        raise AssertionError("amyloid_predictors_summary.tsv combined_rows must be > 0")

    with (root / "amyloid_predictors_status.tsv").open(encoding="utf-8") as handle:
        if any("\tFAIL\t" in line for line in handle):
            raise AssertionError("amyloid_predictors_status.tsv contains FAIL")


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("Usage: validate_amyloid_outputs.py RESULTS_AMYLOID_DIR", file=sys.stderr)
        return 2
    try:
        validate_outputs(Path(argv[0]))
    except AssertionError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print("OK: amyloid outputs match expected schemas")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
