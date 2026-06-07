#!/usr/bin/env python3
"""Normalize HLA caller outputs and build a caller-level consensus genotype."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence


HLA_RE = re.compile(r"(?:HLA[-_])?([A-Z0-9]+)\*([0-9]{2,3}(?::[0-9A-Z]{2,3}){1,3})")
NORMALIZED_FIELDS = ["sample", "caller", "locus", "allele", "rank", "source_file"]
CONSENSUS_FIELDS = [
    "sample",
    "locus",
    "allele",
    "consensus_status",
    "support_count",
    "support_callers",
    "caller_count",
    "source_files",
]


@dataclass(frozen=True)
class HlaCall:
    sample: str
    caller: str
    locus: str
    allele: str
    rank: int
    source_file: str


def normalize_locus(value: str) -> str:
    value = value.strip().upper().replace("_", "-")
    if value.startswith("HLA-"):
        value = value[4:]
    if "*" in value:
        value = value.split("*", 1)[0]
    return value


def normalize_allele(value: str, fallback_locus: str | None = None) -> tuple[str, str] | None:
    value = value.strip().upper().replace("_", "-")
    match = HLA_RE.search(value)
    if not match:
        return None
    locus = normalize_locus(match.group(1) or fallback_locus or "")
    if not locus:
        return None
    allele = f"HLA-{locus}*{match.group(2)}"
    return locus, allele


def iter_strings(value: object) -> Iterator[str]:
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, child in value.items():
            yield str(key)
            yield from iter_strings(child)
    elif isinstance(value, (list, tuple, set)):
        for child in value:
            yield from iter_strings(child)


def unique_calls(sample: str, caller: str, source_file: Path, values: Iterable[tuple[str, str]]) -> list[HlaCall]:
    seen: set[tuple[str, str]] = set()
    calls: list[HlaCall] = []
    for locus, allele in values:
        key = (locus, allele)
        if key in seen:
            continue
        seen.add(key)
        rank = 1 + sum(1 for existing in calls if existing.locus == locus)
        calls.append(HlaCall(sample, caller, locus, allele, rank, str(source_file)))
    return calls


def parse_arcas(sample: str, path: Path) -> list[HlaCall]:
    data = json.loads(path.read_text(encoding="utf-8"))
    values: list[tuple[str, str]] = []
    if isinstance(data, dict):
        for key, value in data.items():
            key_locus = normalize_locus(str(key))
            for text in iter_strings(value):
                call = normalize_allele(text, key_locus)
                if call:
                    values.append(call)
    return unique_calls(sample, "arcasHLA", path, values)


def parse_optitype(sample: str, path: Path) -> list[HlaCall]:
    values: list[tuple[str, str]] = []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        dialect = csv.Sniffer().sniff(handle.read(2048), delimiters="\t,")
        handle.seek(0)
        rows = list(csv.DictReader(handle, dialect=dialect))
    if not rows:
        return []
    for key, value in rows[0].items():
        if not value:
            continue
        key_locus = normalize_locus(key.rstrip("12"))
        call = normalize_allele(value, key_locus)
        if call:
            values.append(call)
    return unique_calls(sample, "OptiType", path, values)


def parse_hisat(sample: str, path: Path) -> list[HlaCall]:
    values: list[tuple[str, str]] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        for match in HLA_RE.finditer(raw):
            call = normalize_allele(match.group(0))
            if call:
                values.append(call)
    return unique_calls(sample, "HISAT-genotype", path, values)


def write_normalized(calls: Sequence[HlaCall], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=NORMALIZED_FIELDS, delimiter="\t")
        writer.writeheader()
        for call in calls:
            writer.writerow(
                {
                    "sample": call.sample,
                    "caller": call.caller,
                    "locus": call.locus,
                    "allele": call.allele,
                    "rank": call.rank,
                    "source_file": call.source_file,
                }
            )


def read_normalized(paths: Sequence[Path]) -> list[HlaCall]:
    calls: list[HlaCall] = []
    for path in paths:
        if not path.exists() or path.stat().st_size == 0:
            continue
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                try:
                    calls.append(
                        HlaCall(
                            sample=row["sample"],
                            caller=row["caller"],
                            locus=normalize_locus(row["locus"]),
                            allele=row["allele"],
                            rank=int(row.get("rank") or 999),
                            source_file=row.get("source_file") or str(path),
                        )
                    )
                except (KeyError, ValueError):
                    continue
    return calls


def build_consensus(calls: Sequence[HlaCall], min_support: int = 2, ploidy: int = 2) -> list[dict[str, str | int]]:
    grouped: dict[tuple[str, str], list[HlaCall]] = defaultdict(list)
    for call in calls:
        grouped[(call.sample, call.locus)].append(call)

    rows: list[dict[str, str | int]] = []
    for (sample, locus), locus_calls in sorted(grouped.items()):
        caller_count = len({call.caller for call in locus_calls})
        by_allele: dict[str, list[HlaCall]] = defaultdict(list)
        for call in locus_calls:
            by_allele[call.allele].append(call)

        ranked = sorted(
            by_allele.items(),
            key=lambda item: (
                -len({call.caller for call in item[1]}),
                min(call.rank for call in item[1]),
                item[0],
            ),
        )
        selected = ranked[:ploidy]
        has_supported = any(len({call.caller for call in allele_calls}) >= min_support for _, allele_calls in selected)

        for allele, allele_calls in selected:
            callers = sorted({call.caller for call in allele_calls})
            support_count = len(callers)
            if support_count >= min_support:
                status = "CONSENSUS"
            elif caller_count <= 1:
                status = "SINGLE_CALLER"
            elif has_supported:
                status = "LOW_SUPPORT"
            else:
                status = "CONFLICT"

            rows.append(
                {
                    "sample": sample,
                    "locus": locus,
                    "allele": allele,
                    "consensus_status": status,
                    "support_count": support_count,
                    "support_callers": ",".join(callers),
                    "caller_count": caller_count,
                    "source_files": ",".join(sorted({call.source_file for call in allele_calls})),
                }
            )
    return rows


def write_consensus(rows: Sequence[dict[str, str | int]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CONSENSUS_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def normalize_command(args: argparse.Namespace) -> int:
    parser = {
        "arcas": parse_arcas,
        "optitype": parse_optitype,
        "hisat": parse_hisat,
    }[args.caller]
    calls = parser(args.sample, args.input)
    write_normalized(calls, args.output)
    return 0


def consensus_command(args: argparse.Namespace) -> int:
    calls = read_normalized(args.inputs)
    rows = build_consensus(calls, min_support=args.min_support, ploidy=args.ploidy)
    write_consensus(rows, args.output)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(required=True)

    normalize = subparsers.add_parser("normalize", help="Normalize a caller-specific HLA output")
    normalize.add_argument("--caller", choices=["arcas", "optitype", "hisat"], required=True)
    normalize.add_argument("--sample", required=True)
    normalize.add_argument("--input", type=Path, required=True)
    normalize.add_argument("--output", type=Path, required=True)
    normalize.set_defaults(func=normalize_command)

    consensus = subparsers.add_parser("consensus", help="Build consensus from normalized TSV files")
    consensus.add_argument("--inputs", type=Path, nargs="+", required=True)
    consensus.add_argument("--output", type=Path, required=True)
    consensus.add_argument("--min-support", type=int, default=2)
    consensus.add_argument("--ploidy", type=int, default=2)
    consensus.set_defaults(func=consensus_command)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
