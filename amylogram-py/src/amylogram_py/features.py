"""AmyloGram 6-mer feature encoder.

The original R package builds binary multigram features through `biogram` with:

```r
ns = c(1, rep(2, 4), rep(3, 3))
ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0))
```

This module intentionally exposes each small block so tests can compare Python
features against R fixtures one layer at a time.
"""

from __future__ import annotations

from itertools import product

from .encoding import decode_sixmer_code

FeatureSpec = tuple[int, tuple[int, ...]]

FEATURE_SPECS: tuple[FeatureSpec, ...] = (
    (1, (0,)),
    (2, (0,)),
    (2, (1,)),
    (2, (2,)),
    (2, (3,)),
    (3, (0, 0)),
    (3, (0, 1)),
    (3, (1, 0)),
)


def feature_name(pattern: tuple[int, ...], distances: tuple[int, ...]) -> str:
    """Return the R/biogram-like feature name for a degenerate multigram."""
    pattern_part = ".".join(str(value + 1) for value in pattern)
    distance_part = ".".join(str(value) for value in distances)
    return f"{pattern_part}_{distance_part}"


def _patterns_in_biogram_order(ngram_size: int):
    """Yield group patterns in the column order used by R biogram."""
    for reversed_pattern in product(range(6), repeat=ngram_size):
        yield tuple(reversed(reversed_pattern))


def all_feature_names() -> list[str]:
    """Generate feature names in the expected biogram order."""
    names: list[str] = []
    for ngram_size, distances in FEATURE_SPECS:
        for pattern in _patterns_in_biogram_order(ngram_size):
            names.append(feature_name(pattern, distances))
    return names


ALL_FEATURE_NAMES = all_feature_names()
FEATURE_INDEX = {name: idx for idx, name in enumerate(ALL_FEATURE_NAMES)}


def _positions_for_spec(start: int, distances: tuple[int, ...]) -> list[int]:
    positions = [start]
    for distance in distances:
        positions.append(positions[-1] + distance + 1)
    return positions


def active_feature_names_for_groups(groups: tuple[int, int, int, int, int, int]) -> set[str]:
    """Return active binary multigram feature names for one 6-mer."""
    active: set[str] = set()
    length = len(groups)
    for ngram_size, distances in FEATURE_SPECS:
        if ngram_size == 1:
            for value in groups:
                active.add(feature_name((value,), distances))
            continue

        max_span = sum(distances) + ngram_size
        for start in range(0, length - max_span + 1):
            positions = _positions_for_spec(start, distances)
            pattern = tuple(groups[pos] for pos in positions)
            active.add(feature_name(pattern, distances))
    return active


def feature_vector_for_code(code: int, selected_features: list[str] | None = None) -> list[int]:
    """Return a binary feature vector for a degenerate 6-mer code."""
    groups = decode_sixmer_code(code)
    active = active_feature_names_for_groups(groups)
    features = selected_features if selected_features is not None else ALL_FEATURE_NAMES
    return [1 if name in active else 0 for name in features]
