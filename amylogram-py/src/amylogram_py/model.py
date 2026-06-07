"""Model-loading primitives for exported AmyloGram ranger forests."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

from .forest import RangerForest, forest_from_export


@dataclass(frozen=True)
class AmyloGramModelMetadata:
    """Minimal metadata needed before implementing tree traversal."""

    imp_features: list[str]
    independent_variable_names: list[str]
    enc: dict[str, list[str]]

    @property
    def feature_count(self) -> int:
        return len(self.imp_features)


@dataclass(frozen=True)
class AmyloGramModel(AmyloGramModelMetadata):
    """Full AmyloGram model exported from R, including the ranger forest."""

    forest: RangerForest


def _strip_r_dataframe_prefix(name: str) -> str:
    return name[1:] if name.startswith("X") else name


def load_model_metadata(path: str | Path) -> AmyloGramModelMetadata:
    """Load and validate the JSON exported from R AmyloGram."""
    payload: dict[str, Any] = json.loads(Path(path).read_text(encoding="utf-8"))
    return model_metadata_from_payload(payload)


def model_metadata_from_payload(payload: dict[str, Any]) -> AmyloGramModelMetadata:
    """Validate metadata shared by the full-model and metadata-only loaders."""
    imp_features = [str(value) for value in payload["imp_features"]]
    independent = [str(value) for value in payload["independent_variable_names"]]
    enc = {str(key): [str(item) for item in value] for key, value in payload["enc"].items()}

    if len(imp_features) != 262:
        raise ValueError(f"Expected 262 AmyloGram selected features, got {len(imp_features)}")
    if len(independent) != 262:
        raise ValueError(f"Expected 262 ranger independent variables, got {len(independent)}")
    if set(enc) != {"1", "2", "3", "4", "5", "6"}:
        raise ValueError(f"Unexpected AmyloGram encoding groups: {sorted(enc)}")
    if [_strip_r_dataframe_prefix(value) for value in independent] != imp_features:
        raise ValueError("Ranger independent variables do not match AmyloGram selected features")

    return AmyloGramModelMetadata(
        imp_features=imp_features,
        independent_variable_names=independent,
        enc=enc,
    )


def load_amylogram_model(path: str | Path) -> AmyloGramModel:
    """Load the full exported AmyloGram model, including the ranger forest."""
    payload: dict[str, Any] = json.loads(Path(path).read_text(encoding="utf-8"))
    metadata = model_metadata_from_payload(payload)
    forest = forest_from_export(payload["forest"], feature_count=metadata.feature_count)
    return AmyloGramModel(
        imp_features=metadata.imp_features,
        independent_variable_names=metadata.independent_variable_names,
        enc=metadata.enc,
        forest=forest,
    )
