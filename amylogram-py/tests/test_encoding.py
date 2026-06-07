from __future__ import annotations

import unittest

from amylogram_py.encoding import (
    AA_TO_GROUP,
    clean_sequence,
    decode_sixmer_code,
    encode_groups,
    iter_rolling_sixmer_codes_from_sequence,
    iter_rolling_sixmer_codes_from_clean_sequence,
    rolling_sixmer_codes,
    sixmer_code,
)


class EncodingTests(unittest.TestCase):
    def test_clean_sequence_matches_current_pipeline_contract(self) -> None:
        self.assertEqual(clean_sequence("mka* zzACDEFG"), "MKAACDEFG")

    def test_amylogram_groups_are_expected(self) -> None:
        expected = {
            "G": 0,
            "K": 1,
            "P": 1,
            "R": 1,
            "I": 2,
            "L": 2,
            "V": 2,
            "F": 3,
            "W": 3,
            "Y": 3,
            "A": 4,
            "C": 4,
            "H": 4,
            "M": 4,
            "D": 5,
            "E": 5,
            "N": 5,
            "Q": 5,
            "S": 5,
            "T": 5,
        }
        self.assertEqual(AA_TO_GROUP, expected)

    def test_sixmer_code_round_trip(self) -> None:
        groups = [2, 5, 2, 2, 3, 1]
        code = sixmer_code(groups)
        self.assertEqual(decode_sixmer_code(code), tuple(groups))

    def test_rolling_codes(self) -> None:
        groups = encode_groups("VQIVYK")
        self.assertEqual(len(rolling_sixmer_codes(groups)), 1)

    def test_streaming_rolling_codes_match_list_path(self) -> None:
        sequence = "xxVQIVYKACD***"
        expected = rolling_sixmer_codes(encode_groups(clean_sequence(sequence)))
        self.assertEqual(list(iter_rolling_sixmer_codes_from_sequence(sequence)), expected)
        self.assertEqual(list(iter_rolling_sixmer_codes_from_clean_sequence(clean_sequence(sequence))), expected)


if __name__ == "__main__":
    unittest.main()
