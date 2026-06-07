from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from amylogram_py.provenance import sha256_file


class ProvenanceTests(unittest.TestCase):
    def test_sha256_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "example.txt"
            path.write_text("amylogram\n", encoding="utf-8")
            self.assertEqual(
                sha256_file(path),
                "911100c78648526b947f527f3a4c6e9a9eea57126c2f91e32b69b377a6b0f797",
            )


if __name__ == "__main__":
    unittest.main()
