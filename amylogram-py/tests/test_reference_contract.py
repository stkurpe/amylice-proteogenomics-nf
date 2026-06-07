from __future__ import annotations

import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]


class ReferenceContractTests(unittest.TestCase):
    def test_r_fixture_exporter_exists(self) -> None:
        script = PROJECT / "scripts" / "export_reference_fixtures.R"
        text = script.read_text(encoding="utf-8")
        self.assertIn("AmyloGram_model", text)
        self.assertIn("count_multigrams", text)
        self.assertIn("predict(model", text)
        self.assertIn("imp_features", text)

    def test_project_documents_parity_blocks(self) -> None:
        readme = (PROJECT / "README.md").read_text(encoding="utf-8")
        for token in (
            "Export AmyloGram reference fixtures",
            "Match amino-acid cleaning and group encoding",
            "Match 6-mer feature vectors",
            "Match ranger forest probabilities",
            "Precompute all `6^6 = 46656`",
        ):
            self.assertIn(token, readme)


if __name__ == "__main__":
    unittest.main()

