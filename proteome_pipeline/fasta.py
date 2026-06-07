from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str

    @property
    def id(self) -> str:
        return self.header.split()[0].lstrip(">")


def read_fasta(path: str | Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    seq_parts: list[str] = []

    with Path(path).open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(FastaRecord(header, "".join(seq_parts)))
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        records.append(FastaRecord(header, "".join(seq_parts)))

    return records


def write_fasta(records: Iterable[FastaRecord], path: str | Path, width: int = 80) -> int:
    count = 0
    with Path(path).open("w") as handle:
        for record in records:
            header = record.header if record.header.startswith(">") else f">{record.header}"
            handle.write(f"{header}\n")
            seq = record.sequence
            for i in range(0, len(seq), width):
                handle.write(f"{seq[i:i + width]}\n")
            count += 1
    return count


def clean_stop_codons(records: Iterable[FastaRecord]) -> list[FastaRecord]:
    cleaned: list[FastaRecord] = []
    for record in records:
        seq = record.sequence.strip()
        if seq.count(".") > 1:
            continue
        if "." in seq:
            seq = seq.split(".", 1)[0]
        if seq:
            cleaned.append(FastaRecord(record.header, seq))
    return cleaned


def combine_fastas(paths: Iterable[str | Path]) -> list[FastaRecord]:
    combined: list[FastaRecord] = []
    seen: set[tuple[str, str]] = set()
    for path in paths:
        p = Path(path)
        if not p.exists() or p.stat().st_size == 0:
            continue
        for record in read_fasta(p):
            key = (record.id, record.sequence)
            if key in seen:
                continue
            seen.add(key)
            combined.append(record)
    return combined
