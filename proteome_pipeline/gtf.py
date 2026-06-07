from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re


MITO_CHROMS = {"chrM", "M", "MT"}


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)
    cds_bp: int = 0

    @property
    def clean_id(self) -> str:
        return self.transcript_id.split(".", 1)[0]


def parse_attributes(attr: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for key, value in re.findall(r'([A-Za-z0-9_]+) "([^"]+)"', attr):
        parsed[key] = value
    return parsed


def load_transcripts(gtf_path: str | Path, allowed_ids: set[str] | None = None) -> dict[str, Transcript]:
    transcripts: dict[str, Transcript] = {}
    with Path(gtf_path).open() as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, feature, start, end, strand, attr = (
                parts[0], parts[2], int(parts[3]), int(parts[4]), parts[6], parts[8]
            )
            if chrom in MITO_CHROMS or feature not in {"exon", "CDS"}:
                continue
            attrs = parse_attributes(attr)
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            clean_tid = tid.split(".", 1)[0]
            if allowed_ids is not None and tid not in allowed_ids and clean_tid not in allowed_ids:
                continue
            gene = attrs.get("gene_name") or attrs.get("gene_id") or tid
            tx = transcripts.setdefault(tid, Transcript(tid, gene, chrom, strand))
            if feature == "exon":
                tx.exons.append((start, end))
            elif feature == "CDS":
                tx.cds_bp += end - start + 1
    return transcripts


def load_id_set(path: str | Path, strip_version: bool = False) -> set[str]:
    ids: set[str] = set()
    with Path(path).open() as handle:
        for raw in handle:
            value = raw.strip()
            if not value:
                continue
            ids.add(value.split(".", 1)[0] if strip_version else value)
    return ids


def expression_ids(abundance_path: str | Path, min_tpm: float = 1.0) -> set[str]:
    ids: set[str] = set()
    with Path(abundance_path).open() as handle:
        for line_no, raw in enumerate(handle, start=1):
            if line_no == 1:
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            try:
                tpm = float(parts[4])
            except ValueError:
                continue
            if tpm > min_tpm:
                ids.add(parts[0].split("|", 1)[0])
    return ids


def clean_transcript_ids(
    gtf_path: str | Path,
    abundance_path: str | Path | None,
    output_path: str | Path,
    min_tpm: float = 1.0,
    min_cds_bp: int = 10,
) -> list[str]:
    expressed = expression_ids(abundance_path, min_tpm) if abundance_path else None
    transcripts = load_transcripts(gtf_path, expressed)
    clean = sorted(
        tx.transcript_id
        for tx in transcripts.values()
        if tx.cds_bp >= min_cds_bp and tx.exons
    )
    Path(output_path).write_text("\n".join(clean) + ("\n" if clean else ""))
    return clean
