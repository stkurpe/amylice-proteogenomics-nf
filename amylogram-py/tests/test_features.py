from __future__ import annotations

import json
import unittest
from pathlib import Path

from amylogram_py.encoding import encode_groups, sixmer_code
from amylogram_py.features import ALL_FEATURE_NAMES, active_feature_names_for_groups, feature_vector_for_code


FIXTURE = Path(__file__).parent / "fixtures" / "amylogram_reference.json"


class FeatureTests(unittest.TestCase):
    def test_all_feature_geometry_matches_amylogram_biogram_setup(self) -> None:
        # 6 monograms + 4 * 36 bigrams + 3 * 216 trigrams
        self.assertEqual(len(ALL_FEATURE_NAMES), 798)
        self.assertEqual(ALL_FEATURE_NAMES[:6], ["1_0", "2_0", "3_0", "4_0", "5_0", "6_0"])
        self.assertIn("1.1_0", ALL_FEATURE_NAMES)
        self.assertIn("6.6.6_1.0", ALL_FEATURE_NAMES)

    def test_active_feature_names_for_single_sixmer(self) -> None:
        groups = tuple(encode_groups("VQIVYK"))
        active = active_feature_names_for_groups(groups)  # type: ignore[arg-type]
        self.assertIn("3_0", active)
        self.assertIn("6_0", active)
        self.assertIn("3.6_0", active)
        self.assertIn("3.3_1", active)
        self.assertIn("3.6.3_0.0", active)

    @unittest.skipUnless(FIXTURE.exists(), "R AmyloGram fixture not generated")
    def test_feature_vectors_match_r_biogram_fixture(self) -> None:
        payload = json.loads(FIXTURE.read_text(encoding="utf-8"))
        for example in payload["examples"]:
            code = sixmer_code(encode_groups(example["sequence"][:6]))
            observed = feature_vector_for_code(code)
            self.assertEqual(observed, example["all_feature_values"], example["id"])


if __name__ == "__main__":
    unittest.main()

