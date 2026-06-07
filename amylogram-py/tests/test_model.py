from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from amylogram_py.model import load_amylogram_model, load_model_metadata
from amylogram_py.predict import predict_sequence_probability_from_model


FIXTURE = Path(__file__).parent / "fixtures" / "amylogram_reference.json"


class ModelMetadataTests(unittest.TestCase):
    def test_load_model_metadata_validates_feature_geometry(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "model.json"
            payload = {
                "enc": {
                    "1": ["g"],
                    "2": ["k", "p", "r"],
                    "3": ["i", "l", "v"],
                    "4": ["f", "w", "y"],
                    "5": ["a", "c", "h", "m"],
                    "6": ["d", "e", "n", "q", "s", "t"],
                },
                "imp_features": [f"f{i}" for i in range(262)],
                "independent_variable_names": [f"Xf{i}" for i in range(262)],
            }
            path.write_text(json.dumps(payload), encoding="utf-8")

            metadata = load_model_metadata(path)
            self.assertEqual(metadata.feature_count, 262)
            self.assertEqual(metadata.enc["1"], ["g"])

    def test_load_model_metadata_rejects_wrong_feature_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.json"
            path.write_text(
                json.dumps(
                    {
                        "enc": {str(i): [] for i in range(1, 7)},
                        "imp_features": ["x"],
                        "independent_variable_names": ["Xx"],
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "262"):
                load_model_metadata(path)

    @unittest.skipUnless(FIXTURE.exists(), "R AmyloGram fixture not generated")
    def test_full_model_matches_r_amylogram_example_probabilities(self) -> None:
        payload = json.loads(FIXTURE.read_text(encoding="utf-8"))
        model = load_amylogram_model(FIXTURE)

        for example in payload["examples"]:
            observed = predict_sequence_probability_from_model(example["sequence"], model)
            self.assertAlmostEqual(observed, example["probability"], places=4, msg=example["id"])


if __name__ == "__main__":
    unittest.main()
