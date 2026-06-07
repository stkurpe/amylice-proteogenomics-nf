"""AmyloGram amino-acid cleaning and degenerate group encoding."""

from __future__ import annotations

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")

GROUPS: tuple[str, ...] = (
    "G",
    "KPR",
    "ILV",
    "FWY",
    "ACHM",
    "DENQST",
)

AA_TO_GROUP: dict[str, int] = {
    aa: group_index
    for group_index, group in enumerate(GROUPS)
    for aa in group
}


def clean_sequence(sequence: str) -> str:
    """Uppercase and remove residues outside the 20 standard amino acids."""
    return "".join(aa for aa in sequence.upper() if aa in STANDARD_AA)


def encode_groups(sequence: str) -> list[int]:
    """Map a cleaned amino-acid sequence to AmyloGram groups 0..5."""
    return [AA_TO_GROUP[aa] for aa in sequence]


def sixmer_code(groups: list[int] | tuple[int, ...]) -> int:
    """Encode six group IDs as a base-6 integer."""
    if len(groups) != 6:
        raise ValueError("sixmer_code expects exactly six groups")
    code = 0
    for value in groups:
        if value < 0 or value > 5:
            raise ValueError(f"group value outside 0..5: {value}")
        code = code * 6 + value
    return code


def decode_sixmer_code(code: int) -> tuple[int, int, int, int, int, int]:
    """Decode a base-6 integer into six group IDs."""
    if code < 0 or code >= 6**6:
        raise ValueError(f"sixmer code outside 0..{6**6 - 1}: {code}")
    out = [0] * 6
    for index in range(5, -1, -1):
        out[index] = code % 6
        code //= 6
    return tuple(out)  # type: ignore[return-value]


def rolling_sixmer_codes(groups: list[int]) -> list[int]:
    """Return base-6 codes for all overlapping 6-mers in a group sequence."""
    if len(groups) < 6:
        return []
    return [sixmer_code(groups[i : i + 6]) for i in range(len(groups) - 5)]


def iter_clean_group_ids(sequence: str):
    """Yield AmyloGram group IDs while applying the standard amino-acid filter."""
    for aa in sequence.upper():
        group = AA_TO_GROUP.get(aa)
        if group is not None:
            yield group


def iter_rolling_sixmer_codes_from_sequence(sequence: str):
    """Yield overlapping 6-mer codes without allocating the full group list."""
    code = 0
    length = 0
    modulus = 6**5
    for group in iter_clean_group_ids(sequence):
        if length < 6:
            code = code * 6 + group
            length += 1
            if length == 6:
                yield code
        else:
            code = (code % modulus) * 6 + group
            yield code


def iter_rolling_sixmer_codes_from_clean_sequence(sequence: str):
    """Yield overlapping 6-mer codes for an already cleaned uppercase sequence."""
    if len(sequence) < 6:
        return
    code = 0
    modulus = 6**5
    for index, aa in enumerate(sequence):
        group = AA_TO_GROUP[aa]
        if index < 6:
            code = code * 6 + group
            if index == 5:
                yield code
        else:
            code = (code % modulus) * 6 + group
            yield code
