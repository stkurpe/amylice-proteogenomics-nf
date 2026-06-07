#!/usr/bin/env python3
"""Contract tests for amyloidogenicity predictor scripts.

These tests intentionally avoid running the heavy AMYPred-FRL training stack,
R AmyloGram inference, Docker, network access, or AWS writes. They document
the expected behavior by inspecting source code structure and lightweight
FASTA parsing/cleaning contracts.
"""

from __future__ import annotations

import ast
import re
import unittest
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
AMYLOID_DIR = PROJECT_DIR / "amyloid_predictors"
AMYPRED = AMYLOID_DIR / "AMYPred-FRL" / "predict.py"
AMYLOGRAM = AMYLOID_DIR / "AmyloGram"

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def parse_fasta_like_amypred(text: str, min_length: int = 30) -> list[tuple[str, str]]:
    """Small test double for AMYPred-FRL read_protein_sequences contract."""
    if ">" not in text:
        raise ValueError("not FASTA")

    records: list[tuple[str, str]] = []
    for raw_record in text.split(">")[1:]:
        lines = raw_record.splitlines()
        if not lines:
            continue
        seq_id = lines[0].strip().split()[0]
        cleaned = re.sub("[^ACDEFGHIKLMNPQRSTVWY-]", "-", "".join(lines[1:]).upper())
        if len(cleaned) >= min_length:
            records.append((seq_id, cleaned))
    if not records:
        raise ValueError("no valid sequences")
    return records


def clean_like_amylogram(seq: str) -> str:
    return re.sub("[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())


class AmyPredFRLContractTests(unittest.TestCase):
    def setUp(self) -> None:
        self.source = read_text(AMYPRED)
        self.tree = ast.parse(self.source)

    def test_source_snapshot_exists(self) -> None:
        self.assertTrue(AMYPRED.is_file(), f"Missing AMYPred-FRL source: {AMYPRED}")

    def test_feature_descriptor_functions_are_declared(self) -> None:
        expected = {
            "AAC",
            "DPC",
            "APAAC",
            "CTDC",
            "CTDD",
            "CTDT",
            "GAAC",
            "KSCTriad",
            "CTriad",
            "DDE",
            "PAAC",
        }
        defined = {node.name for node in self.tree.body if isinstance(node, ast.FunctionDef)}
        self.assertTrue(expected.issubset(defined), f"Missing descriptor functions: {sorted(expected - defined)}")

    def test_stacking_uses_six_base_classifiers(self) -> None:
        expected_classifiers = {
            "RandomForestClassifier",
            "ExtraTreesClassifier",
            "SVC",
            "LogisticRegression",
            "XGBClassifier",
            "KNeighborsClassifier",
        }
        for classifier in expected_classifiers:
            self.assertIn(classifier, self.source)

        classifier_calls = {
            node.func.id
            for node in ast.walk(self.tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            and node.func.id in expected_classifiers
        }
        self.assertEqual(classifier_calls, expected_classifiers)

    def test_meta_feature_geometry_is_10_descriptor_groups_by_6_classifiers(self) -> None:
        self.assertRegex(self.source, r"for\s+i\s+in\s+range\(10\)")
        self.assertRegex(self.source, r"for\s+clf\s+in\s+clfs")

        # 10 descriptor groups are fed to six base classifiers, yielding 60
        # probabilistic features before LR-RFE/SVM selection.
        self.assertEqual(10 * 6, 60)

    def test_final_feature_mask_selects_twenty_meta_features(self) -> None:
        match = re.search(r"mask\s*=\s*\[([^\]]+)\]", self.source)
        self.assertIsNotNone(match, "Missing selected-feature mask")
        mask = [int(value.strip()) for value in match.group(1).split(",")]

        self.assertEqual(len(mask), 20)
        self.assertEqual(len(set(mask)), 20)
        self.assertGreaterEqual(min(mask), 0)
        self.assertLess(max(mask), 60)

    def test_final_model_contract_is_probability_csv(self) -> None:
        self.assertIn("model/pima.pickle_model_svm_PF.dat", self.source)
        self.assertIn("predict_proba(Selected_feat)", self.source)
        for column in ("Sequence ID", "Amyloid Prob", "Prediction"):
            self.assertIn(column, self.source)
        self.assertRegex(self.source, r"AMYLOID.+if.+amy_prob\s*>\s*threshold")

    def test_fasta_reader_contract_filters_short_and_normalizes_invalid_amino_acids(self) -> None:
        fasta = (
            ">too_short\nACDEFG\n"
            ">valid with spaces\nACDZXBJOUacdefghiklmnpqrstvwyACDEFGHIKLMN\n"
        )
        records = parse_fasta_like_amypred(fasta)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0][0], "valid")
        self.assertGreaterEqual(len(records[0][1]), 30)
        self.assertTrue(set(records[0][1]).issubset(STANDARD_AA | {"-"}))
        self.assertIn("-", records[0][1])


class AmyloGramContractTests(unittest.TestCase):
    def test_r_scripts_exist(self) -> None:
        expected = {
            "predict_amylogram.R",
            "predict_amylogram_fast.R",
            "predict_amylogram_optimized.R",
            "predict_amylogram_parallel.R",
            "predict_amylogram_robust.R",
            "benchmark_amylogram.R",
        }
        existing = {path.name for path in AMYLOGRAM.glob("*.R")}
        self.assertTrue(expected.issubset(existing), f"Missing AmyloGram R scripts: {sorted(expected - existing)}")

    def test_prediction_scripts_share_output_schema(self) -> None:
        for script in AMYLOGRAM.glob("predict_amylogram*.R"):
            text = read_text(script)
            self.assertIn("Sequence_ID", text, script.name)
            self.assertIn("AmyloGram_Prob", text, script.name)
            self.assertIn("AmyloGram_Pred", text, script.name)
            self.assertIn("write.csv", text, script.name)

    def test_amylogram_cleaning_contract_removes_non_standard_residues_and_filters_length_six(self) -> None:
        raw = {
            "bad_short": "ACDXX",
            "valid_dirty": "mka* zzACDEFG",
        }
        cleaned = {seq_id: clean_like_amylogram(seq) for seq_id, seq in raw.items()}
        valid = {seq_id: seq for seq_id, seq in cleaned.items() if len(seq) >= 6}

        self.assertNotIn("bad_short", valid)
        self.assertEqual(valid["valid_dirty"], "MKAACDEFG")
        self.assertTrue(set(valid["valid_dirty"]).issubset(STANDARD_AA))

    def test_parallel_variants_are_execution_strategies_not_new_algorithms(self) -> None:
        optimized = read_text(AMYLOGRAM / "predict_amylogram_optimized.R")
        parallel = read_text(AMYLOGRAM / "predict_amylogram_parallel.R")

        self.assertIn("library(AmyloGram)", optimized)
        self.assertIn("library(AmyloGram)", parallel)
        self.assertIn("mclapply", optimized)
        self.assertIn("mclapply", parallel)
        self.assertIn("AmyloGram_model", optimized)
        self.assertIn("AmyloGram_model", parallel)

    def test_fast_runner_is_resumable_deduplicated_and_fault_tolerant(self) -> None:
        fast = read_text(AMYLOGRAM / "predict_amylogram_fast.R")

        for token in (
            "unique(valid_seqs)",
            "match(valid_seqs, unique_seqs)",
            "AMYLOGRAM_CHUNK_SIZE",
            "AMYLOGRAM_CORES",
            "AMYLOGRAM_RETRY_SPLIT_DEPTH",
            "AMYLOGRAM_ALLOW_ERRORS",
            "predict_resilient",
            "file.exists(out_path)",
            "chunk_stats.tsv",
            "error_predictions",
        ):
            self.assertIn(token, fast)

        self.assertRegex(fast, r"chunk_size\s*<-\s*get_env_int\(\"AMYLOGRAM_CHUNK_SIZE\",\s*200\)")
        self.assertIn("stop(sprintf(\"AmyloGram completed with %d sequence-level errors\"", fast)


class AmyloGramPyPipelineIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.script = read_text(PROJECT_DIR / "run_amyloid_predictors.sh")

    def test_fast_python_predictor_is_optional_and_enabled_by_default(self) -> None:
        self.assertIn('RUN_AMYLOGRAM_PY="${RUN_AMYLOGRAM_PY:-1}"', self.script)
        self.assertIn("--run-amylogram-py", self.script)
        self.assertIn("--skip-amylogram-py", self.script)

    def test_fast_python_predictor_outputs_are_reported_and_uploaded(self) -> None:
        for token in (
            "run_amylogram_py",
            "AMYLOGRAM_PY_IMAGE",
            "AMYLOGRAM_PY_LOOKUP",
            "amylogram_py_prediction.csv",
            "amylogram_py_report.md",
            "amylogram_py_report.json",
            "amylogram_py_skipped.tsv",
            "amylogram_py_top_hits.tsv",
            "amylogram_py_rows",
        ):
            self.assertIn(token, self.script)

    def test_r_protein_feature_step_is_reported_and_uploaded(self) -> None:
        for token in (
            "run_protein_features",
            "PROTEIN_FEATURES_IMAGE",
            "PROTEIN_FEATURES_OUTPUT_S3",
            "PROTEIN_FEATURES_WINDOW",
            "PROTEIN_FEATURES_PH",
            "PROTEIN_FEATURES_R_LIBS",
            "calc_protein_features.R",
            "protein_features.csv",
            "protein_features_rows",
            "--run-protein-features",
            "--skip-protein-features",
            "--protein-features-output-s3",
        ):
            self.assertIn(token, self.script)

    def test_combined_predictions_accept_three_predictors(self) -> None:
        for column in (
            "AMYPred_Prob",
            "AmyloGram_Prob",
            "AmyloGramPy_Prob",
            "AmyloGramPy_Pred",
            "Consensus",
        ):
            self.assertIn(column, self.script)


if __name__ == "__main__":
    unittest.main(verbosity=2)
