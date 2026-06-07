"""Prediction helpers for the lookup-table AmyloGram-Py path."""

from __future__ import annotations

from .encoding import (
    encode_groups,
    iter_rolling_sixmer_codes_from_clean_sequence,
    iter_rolling_sixmer_codes_from_sequence,
    rolling_sixmer_codes,
)
from .features import feature_vector_for_code
from .model import AmyloGramModel


def predict_sequence_probability(sequence: str, sixmer_probabilities: list[float]) -> float:
    """Predict a protein by max probability across its overlapping 6-mers."""
    return predict_clean_sequence_probability_from_codes(
        iter_rolling_sixmer_codes_from_sequence(sequence),
        sixmer_probabilities,
    )


def predict_clean_sequence_probability(sequence: str, sixmer_probabilities: list[float]) -> float:
    """Predict an already cleaned uppercase protein sequence."""
    return predict_clean_sequence_probability_from_codes(
        iter_rolling_sixmer_codes_from_clean_sequence(sequence),
        sixmer_probabilities,
    )


def predict_clean_sequence_probability_from_codes(codes, sixmer_probabilities: list[float]) -> float:
    """Predict from a stream of six-mer codes."""
    max_probability: float | None = None
    for code in codes:
        probability = sixmer_probabilities[code]
        if max_probability is None or probability > max_probability:
            max_probability = probability

    if max_probability is None:
        raise ValueError("AmyloGram requires sequences of at least six amino acids")
    return max_probability


def predict_sixmer_probability_from_model(sequence: str, model: AmyloGramModel) -> float:
    """Predict one six-amino-acid window directly through the ranger forest."""
    groups = encode_groups(sequence)
    codes = rolling_sixmer_codes(groups)
    if len(codes) != 1:
        raise ValueError("Expected exactly one six-amino-acid window")
    features = feature_vector_for_code(codes[0], model.imp_features)
    return model.forest.predict_probability(features)


def predict_sequence_probability_from_model(sequence: str, model: AmyloGramModel) -> float:
    """Predict a protein directly through the ranger forest without a lookup table."""
    groups = encode_groups(sequence)
    codes = rolling_sixmer_codes(groups)
    if not codes:
        raise ValueError("AmyloGram requires sequences of at least six amino acids")
    return max(
        model.forest.predict_probability(feature_vector_for_code(code, model.imp_features))
        for code in codes
    )


def probability_to_label(probability: float, threshold: float = 0.5) -> str:
    return "AMYLOID" if probability > threshold else "Non-Amyloid"
