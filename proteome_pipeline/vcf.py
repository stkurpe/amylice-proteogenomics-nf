from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import gzip


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def is_frameshift(self) -> bool:
        diff = len(self.alt) - len(self.ref)
        return diff != 0 and abs(diff) % 3 != 0


def _open_text(path: str | Path):
    p = Path(path)
    return gzip.open(p, "rt") if p.suffix == ".gz" else p.open()


def iter_variants(path: str | Path):
    with _open_text(path) as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            yield Variant(parts[0], int(parts[1]), parts[3].upper(), parts[4].split(",", 1)[0].upper())


def iter_frameshift_variants(path: str | Path):
    for variant in iter_variants(path):
        if variant.is_frameshift:
            yield variant
