"""Python implementation scaffold for AmyloGram-compatible prediction."""

from .encoding import AA_TO_GROUP, clean_sequence, encode_groups, rolling_sixmer_codes

__all__ = [
    "AA_TO_GROUP",
    "clean_sequence",
    "encode_groups",
    "rolling_sixmer_codes",
]

