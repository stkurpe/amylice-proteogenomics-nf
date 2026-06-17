#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys


def fail(message: str) -> None:
    raise AssertionError(message)


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


def sequence_set_md5(records: dict[str, str]) -> str:
    payload = "\n".join(f"{seq_id}\t{records[seq_id]}" for seq_id in sorted(records))
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


def compare_exact(expected: Path, actual: Path, expected_count: int | None, expected_md5: str | None) -> str:
    expected_records = read_fasta(expected)
    actual_records = read_fasta(actual)

    missing = sorted(set(expected_records) - set(actual_records))
    extra = sorted(set(actual_records) - set(expected_records))
    common = sorted(set(expected_records) & set(actual_records))
    mismatched = [seq_id for seq_id in common if expected_records[seq_id] != actual_records[seq_id]]

    expected_digest = sequence_set_md5(expected_records)
    actual_digest = sequence_set_md5(actual_records)

    lines = [
        f"expected_records\t{len(expected_records)}",
        f"actual_records\t{len(actual_records)}",
        f"expected_unique_ids\t{len(expected_records)}",
        f"actual_unique_ids\t{len(actual_records)}",
        f"missing_in_actual\t{len(missing)}",
        f"extra_in_actual\t{len(extra)}",
        f"sequence_mismatches_common_ids\t{len(mismatched)}",
        f"expected_sequence_set_md5\t{expected_digest}",
        f"actual_sequence_set_md5\t{actual_digest}",
    ]

    if expected_count is not None and len(actual_records) != expected_count:
        fail(f"actual record count differs from expected count: actual={len(actual_records)} expected={expected_count}")
    if expected_md5 is not None and actual_digest != expected_md5:
        fail(f"actual sequence-set md5 differs: actual={actual_digest} expected={expected_md5}")
    if missing or extra:
        fail(f"FASTA sequence IDs differ; missing={missing[:10]} extra={extra[:10]}")
    if mismatched:
        fail(f"FASTA sequences differ for IDs: {mismatched[:10]}")
    if expected_digest != actual_digest:
        fail(f"FASTA sequence-set md5 differs: expected={expected_digest} actual={actual_digest}")

    return "\n".join(["OK: FASTA files match exactly", *lines])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare two FASTA files by ID and sequence.")
    parser.add_argument("--expected", required=True, type=Path, help="Reference FASTA.")
    parser.add_argument("--actual", required=True, type=Path, help="Generated FASTA.")
    parser.add_argument("--expected-count", type=int, help="Expected record count for the generated FASTA.")
    parser.add_argument("--expected-md5", help="Expected sorted ID/sequence MD5 for the generated FASTA.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        print(compare_exact(args.expected, args.actual, args.expected_count, args.expected_md5))
    except AssertionError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
