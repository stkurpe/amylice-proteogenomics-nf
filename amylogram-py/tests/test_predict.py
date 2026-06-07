from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from amylogram_py.cli import main
from amylogram_py.encoding import encode_groups, rolling_sixmer_codes
from amylogram_py.lookup import write_probability_table_binary
from amylogram_py.predict import predict_clean_sequence_probability, predict_sequence_probability, probability_to_label


class PredictTests(unittest.TestCase):
    def test_sequence_prediction_uses_max_over_sixmers(self) -> None:
        table = [0.0] * (6**6)
        codes = rolling_sixmer_codes(encode_groups("VQIVYKACD"))
        table[codes[0]] = 0.25
        table[codes[1]] = 0.91
        table[codes[2]] = 0.10
        table[codes[3]] = 0.50
        self.assertEqual(predict_sequence_probability("VQIVYKACD", table), 0.91)
        self.assertEqual(predict_sequence_probability("xxVQIVYKACD***", table), 0.91)
        self.assertEqual(predict_clean_sequence_probability("VQIVYKACD", table), 0.91)

    def test_probability_label_uses_strict_greater_than_threshold(self) -> None:
        self.assertEqual(probability_to_label(0.50001), "AMYLOID")
        self.assertEqual(probability_to_label(0.5), "Non-Amyloid")

    def test_cli_writes_expected_schema(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            table_path = tmp_path / "table.json"
            table_bin_path = tmp_path / "table.bin"
            fasta_path = tmp_path / "input.fasta"
            output_path = tmp_path / "out.csv"
            output_bin_path = tmp_path / "out_bin.csv"
            report_path = tmp_path / "report.md"
            report_json_path = tmp_path / "report.json"
            skipped_path = tmp_path / "skipped.tsv"
            top_path = tmp_path / "top.tsv"
            table = [0.1] * (6**6)
            table[rolling_sixmer_codes(encode_groups("VQIVYK"))[0]] = 0.9
            table[rolling_sixmer_codes(encode_groups("ACDEFG"))[0]] = 0.7
            table_path.write_text("[" + ",".join(str(value) for value in table) + "]", encoding="utf-8")
            write_probability_table_binary(table, table_bin_path)
            fasta_path.write_text(">seq1\nVQIVYK\n>seq2\nACDEFG\n>short\nACD\n>dirty\n***\n", encoding="utf-8")

            rc = main(
                [
                    str(fasta_path),
                    str(output_path),
                    "--sixmer-table",
                    str(table_path),
                    "--report-md",
                    str(report_path),
                    "--report-json",
                    str(report_json_path),
                    "--skipped-tsv",
                    str(skipped_path),
                    "--top-k",
                    "1",
                    "--top-tsv",
                    str(top_path),
                ]
            )
            self.assertEqual(rc, 0)

            with output_path.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["Sequence_ID"], "seq1")
            self.assertEqual(set(rows[0]), {"Sequence_ID", "AmyloGram_Prob", "AmyloGram_Pred"})
            self.assertIn("| Total FASTA records | 4 |", report_path.read_text(encoding="utf-8"))
            self.assertIn("| Skipped records | 2 |", report_path.read_text(encoding="utf-8"))
            self.assertIn("| Max probability | 0.900000 |", report_path.read_text(encoding="utf-8"))
            self.assertIn("Input FASTA SHA256", report_path.read_text(encoding="utf-8"))
            report_json = json.loads(report_json_path.read_text(encoding="utf-8"))
            self.assertEqual(report_json["total_records"], 4)
            self.assertEqual(report_json["skipped_records"], 2)
            self.assertEqual(report_json["max_sequence_id"], "seq1")
            self.assertAlmostEqual(report_json["mean_probability"], 0.8)
            self.assertEqual(len(report_json["input_fasta_sha256"]), 64)
            self.assertEqual(len(report_json["sixmer_table_sha256"]), 64)
            self.assertIn("shorter_than_6", skipped_path.read_text(encoding="utf-8"))
            self.assertIn("empty_after_cleaning", skipped_path.read_text(encoding="utf-8"))
            with top_path.open(encoding="utf-8") as handle:
                top_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(top_rows), 1)
            self.assertEqual(top_rows[0]["Sequence_ID"], "seq1")

            rc = main([str(fasta_path), str(output_bin_path), "--sixmer-table", str(table_bin_path)])
            self.assertEqual(rc, 0)


if __name__ == "__main__":
    unittest.main()
