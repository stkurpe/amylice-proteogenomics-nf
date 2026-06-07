"""Build and use the 6^6 AmyloGram probability lookup table."""

from __future__ import annotations

import json
from pathlib import Path
import struct
from typing import Sequence

from .features import feature_vector_for_code
from .forest import RangerForest

SIXMER_SPACE = 6**6
LOOKUP_MAGIC = b"AMYLOGRAM_PY_LOOKUP_V1\n"
LOOKUP_STRUCT = struct.Struct(f"<{SIXMER_SPACE}d")


def build_probability_table(forest: RangerForest, selected_features: list[str]) -> list[float]:
    """Precompute probabilities for every degenerate 6-mer code."""
    return [
        forest.predict_probability(feature_vector_for_code(code, selected_features))
        for code in range(SIXMER_SPACE)
    ]


def write_probability_table(table: Sequence[float], path: str | Path) -> None:
    if len(table) != SIXMER_SPACE:
        raise ValueError(f"Expected {SIXMER_SPACE} probabilities, got {len(table)}")
    Path(path).write_text(json.dumps([float(value) for value in table]), encoding="utf-8")


def read_probability_table(path: str | Path) -> list[float]:
    values = json.loads(Path(path).read_text(encoding="utf-8"))
    if len(values) != SIXMER_SPACE:
        raise ValueError(f"Expected {SIXMER_SPACE} probabilities, got {len(values)}")
    return [float(value) for value in values]


def write_probability_table_binary(table: Sequence[float], path: str | Path) -> None:
    """Write a compact float64 lookup table with a small magic header."""
    if len(table) != SIXMER_SPACE:
        raise ValueError(f"Expected {SIXMER_SPACE} probabilities, got {len(table)}")
    with Path(path).open("wb") as handle:
        handle.write(LOOKUP_MAGIC)
        handle.write(LOOKUP_STRUCT.pack(*[float(value) for value in table]))


def read_probability_table_binary(path: str | Path) -> list[float]:
    """Read a compact float64 lookup table."""
    with Path(path).open("rb") as handle:
        magic = handle.read(len(LOOKUP_MAGIC))
        if magic != LOOKUP_MAGIC:
            raise ValueError("Not an AmyloGram-Py binary lookup table")
        payload = handle.read()
    expected_size = LOOKUP_STRUCT.size
    if len(payload) != expected_size:
        raise ValueError(f"Expected {expected_size} lookup bytes, got {len(payload)}")
    return list(LOOKUP_STRUCT.unpack(payload))


def read_probability_table_auto(path: str | Path) -> list[float]:
    """Read JSON or binary lookup tables by inspecting the file header."""
    path = Path(path)
    with path.open("rb") as handle:
        prefix = handle.read(len(LOOKUP_MAGIC))
    if prefix == LOOKUP_MAGIC:
        return read_probability_table_binary(path)
    return read_probability_table(path)
