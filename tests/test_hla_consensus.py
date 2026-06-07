#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from proteome_pipeline.hla_consensus import (
    build_consensus,
    parse_arcas,
    parse_hisat,
    parse_optitype,
    read_normalized,
    write_normalized,
)


class HlaParserTests(unittest.TestCase):
    def test_parse_arcas_json_with_locus_keys(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "arcas.json"
            path.write_text(json.dumps({"A": ["A*02:01", "HLA-A*24:02"], "B": ["B*07:02"]}), encoding="utf-8")
            calls = parse_arcas("SRR1", path)

        self.assertEqual([call.allele for call in calls], ["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02"])
        self.assertEqual([call.rank for call in calls if call.locus == "A"], [1, 2])

    def test_parse_optitype_result_tsv(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "result.tsv"
            path.write_text("A1\tA2\tB1\tC1\nA*01:01\tA*02:01\tB*08:01\tC*07:01\n", encoding="utf-8")
            calls = parse_optitype("SRR1", path)

        self.assertEqual({call.allele for call in calls}, {"HLA-A*01:01", "HLA-A*02:01", "HLA-B*08:01", "HLA-C*07:01"})

    def test_parse_hisat_report_text(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "hisat.log"
            path.write_text("ranked allele HLA-DRB1*15:01 score\nsecond DQB1*06:02\n", encoding="utf-8")
            calls = parse_hisat("SRR1", path)

        self.assertEqual([call.allele for call in calls], ["HLA-DRB1*15:01", "HLA-DQB1*06:02"])


class HlaConsensusTests(unittest.TestCase):
    def test_consensus_marks_two_caller_support(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            normalized = Path(tmp) / "calls.tsv"
            calls = parse_arcas_from_values(
                [
                    ("SRR1", "arcasHLA", "A", "HLA-A*02:01", 1),
                    ("SRR1", "OptiType", "A", "HLA-A*02:01", 1),
                    ("SRR1", "HISAT-genotype", "A", "HLA-A*24:02", 1),
                ],
                normalized,
            )
            rows = build_consensus(calls)

        self.assertEqual(rows[0]["allele"], "HLA-A*02:01")
        self.assertEqual(rows[0]["consensus_status"], "CONSENSUS")
        self.assertEqual(rows[0]["support_count"], 2)
        self.assertEqual(rows[1]["consensus_status"], "LOW_SUPPORT")

    def test_conflict_when_callers_disagree(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            normalized = Path(tmp) / "calls.tsv"
            calls = parse_arcas_from_values(
                [
                    ("SRR1", "arcasHLA", "B", "HLA-B*07:02", 1),
                    ("SRR1", "OptiType", "B", "HLA-B*08:01", 1),
                    ("SRR1", "HISAT-genotype", "B", "HLA-B*15:01", 1),
                ],
                normalized,
            )
            rows = build_consensus(calls)

        self.assertEqual({row["consensus_status"] for row in rows}, {"CONFLICT"})


def parse_arcas_from_values(values: list[tuple[str, str, str, str, int]], output: Path):
    from proteome_pipeline.hla_consensus import HlaCall

    calls = [
        HlaCall(sample=sample, caller=caller, locus=locus, allele=allele, rank=rank, source_file="fixture")
        for sample, caller, locus, allele, rank in values
    ]
    write_normalized(calls, output)
    round_tripped = read_normalized([output])
    with output.open(newline="", encoding="utf-8") as handle:
        self_check = list(csv.DictReader(handle, delimiter="\t"))
    assert self_check
    return round_tripped


if __name__ == "__main__":
    unittest.main()
