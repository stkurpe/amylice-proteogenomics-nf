"""Small FASTA reader used by the AmyloGram-Py CLI."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from .encoding import clean_sequence


@dataclass(frozen=True)
class FastaRecord:
    sequence_id: str
    sequence: str
    raw_length: int
    clean_length: int
    status: str
    reason: str


def iter_fasta_records(path: str | Path, min_length: int = 6) -> Iterator[FastaRecord]:
    """Yield every FASTA record with cleaning status and skip reason."""
    current_id: str | None = None
    chunks: list[str] = []

    def emit() -> FastaRecord | None:
        if current_id is None:
            return None
        raw_sequence = "".join(chunks)
        sequence = clean_sequence(raw_sequence)
        if not sequence:
            status = "skipped"
            reason = "empty_after_cleaning"
        elif len(sequence) < min_length:
            status = "skipped"
            reason = f"shorter_than_{min_length}"
        else:
            status = "ok"
            reason = ""
        return FastaRecord(
            sequence_id=current_id,
            sequence=sequence,
            raw_length=len(raw_sequence),
            clean_length=len(sequence),
            status=status,
            reason=reason,
        )

    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                record = emit()
                if record is not None:
                    yield record
                current_id = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line)

    record = emit()
    if record is not None:
        yield record


def read_fasta(path: str | Path, min_length: int = 6) -> Iterator[tuple[str, str]]:
    """Yield cleaned FASTA records with sequence length >= min_length."""
    for record in iter_fasta_records(path, min_length=min_length):
        if record.status == "ok":
            yield record.sequence_id, record.sequence
