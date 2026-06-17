#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import sys


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
REQUIRED_FILES = [
    "clean_ids.txt",
    "protein.fasta",
    "frameshift_unique.fasta",
    "combine_proteome.fasta",
    "verification_report.tsv",
]
FASTA_FILES = [
    "protein.fasta",
    "frameshift_unique.fasta",
    "combine_proteome.fasta",
]


@dataclass(frozen=True)
class FastaComparison:
    filename: str
    record_count: int


def fail(message: str) -> None:
    raise AssertionError(message)


def read_clean_ids(path: Path) -> list[str]:
    return sorted(line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip())


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    seq_parts: list[str] = []

    def flush() -> None:
        nonlocal header, seq_parts
        if header is None:
            return
        seq_id = header.split()[0].lstrip(">")
        sequence = "".join(seq_parts).strip().upper()
        if not seq_id:
            fail(f"{path}: FASTA record has empty ID")
        if seq_id in records:
            fail(f"{path}: duplicate FASTA ID: {seq_id}")
        if not sequence:
            fail(f"{path}: FASTA record has empty sequence: {seq_id}")
        records[seq_id] = sequence

    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
    flush()
    return records


def invalid_residue_summary(records: dict[str, str]) -> dict[str, str]:
    invalid: dict[str, str] = {}
    for seq_id, sequence in records.items():
        bad = sorted(set(sequence) - STANDARD_AA)
        if bad:
            invalid[seq_id] = "".join(bad)
    return invalid


def assert_no_invalid_aa_explosion(label: str, records: dict[str, str]) -> None:
    invalid = invalid_residue_summary(records)
    if not invalid:
        return
    total_records = len(records)
    invalid_records = len(invalid)
    if invalid_records / max(total_records, 1) > 0.05:
        examples = ", ".join(f"{seq_id}:{chars}" for seq_id, chars in list(invalid.items())[:5])
        fail(f"{label}: invalid amino acid explosion: {invalid_records}/{total_records} records; {examples}")


def assert_verification_all_ok(path: Path) -> None:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        statuses = [row.get("status", "") for row in reader]
    if not statuses:
        fail(f"{path}: verification report has no check rows")
    bad = [status for status in statuses if status != "OK"]
    if bad:
        fail(f"{path}: verification report contains non-OK statuses: {bad[:5]}")


def require_files(root: Path) -> None:
    for filename in REQUIRED_FILES:
        path = root / filename
        if not path.is_file():
            fail(f"Missing required file: {path}")
        if path.stat().st_size == 0:
            fail(f"Required file is empty: {path}")


def compare_fasta_file(expected_root: Path, nf_root: Path, filename: str) -> FastaComparison:
    expected_records = read_fasta(expected_root / filename)
    nf_records = read_fasta(nf_root / filename)

    assert_no_invalid_aa_explosion(f"expected {filename}", expected_records)
    assert_no_invalid_aa_explosion(f"nextflow {filename}", nf_records)

    expected_ids = set(expected_records)
    nf_ids = set(nf_records)
    if expected_ids != nf_ids:
        missing = sorted(expected_ids - nf_ids)
        extra = sorted(nf_ids - expected_ids)
        fail(f"{filename}: sequence IDs differ; missing={missing[:10]} extra={extra[:10]}")

    mismatched = [
        seq_id
        for seq_id in sorted(expected_ids)
        if expected_records[seq_id] != nf_records[seq_id]
    ]
    if mismatched:
        fail(f"{filename}: sequences differ for IDs: {mismatched[:10]}")

    if len(expected_records) != len(nf_records):
        fail(f"{filename}: record count differs: expected={len(expected_records)} nextflow={len(nf_records)}")

    return FastaComparison(filename=filename, record_count=len(nf_records))


def compare_outputs(expected_root: Path, nf_root: Path) -> list[FastaComparison]:
    require_files(expected_root)
    require_files(nf_root)

    expected_clean = read_clean_ids(expected_root / "clean_ids.txt")
    nf_clean = read_clean_ids(nf_root / "clean_ids.txt")
    if expected_clean != nf_clean:
        fail("clean_ids.txt differs after sorting")

    assert_verification_all_ok(expected_root / "verification_report.tsv")
    assert_verification_all_ok(nf_root / "verification_report.tsv")

    return [compare_fasta_file(expected_root, nf_root, filename) for filename in FASTA_FILES]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare expected fixture and Nextflow proteome outputs.")
    parser.add_argument("--expected-output", required=True, type=Path, help="Expected results_proteins fixture directory.")
    parser.add_argument("--nf-output", required=True, type=Path, help="Nextflow results_proteins directory.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        comparisons = compare_outputs(args.expected_output, args.nf_output)
    except AssertionError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1

    print("OK: proteome outputs match")
    for comparison in comparisons:
        print(f"{comparison.filename}\trecords={comparison.record_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
