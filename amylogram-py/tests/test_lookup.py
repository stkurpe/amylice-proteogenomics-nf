from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from amylogram_py.forest import RangerForest, RangerTree
from amylogram_py.lookup import (
    SIXMER_SPACE,
    build_probability_table,
    read_probability_table,
    read_probability_table_auto,
    read_probability_table_binary,
    write_probability_table,
    write_probability_table_binary,
)


class LookupTests(unittest.TestCase):
    def test_build_probability_table_covers_full_degenerate_sixmer_space(self) -> None:
        forest = RangerForest([RangerTree([0], [0], [0], [0.0], [[3, 7]])])
        table = build_probability_table(forest, selected_features=[])
        self.assertEqual(len(table), SIXMER_SPACE)
        self.assertTrue(all(value == 0.7 for value in table))

    def test_probability_table_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.json"
            table = [0.25] * SIXMER_SPACE
            write_probability_table(table, path)
            self.assertEqual(read_probability_table(path), table)

    def test_binary_probability_table_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.bin"
            table = [float(code) / SIXMER_SPACE for code in range(SIXMER_SPACE)]
            write_probability_table_binary(table, path)
            self.assertEqual(read_probability_table_binary(path), table)
            self.assertEqual(read_probability_table_auto(path), table)

    def test_auto_probability_table_reader_accepts_json(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.json"
            table = [0.125] * SIXMER_SPACE
            write_probability_table(table, path)
            self.assertEqual(read_probability_table_auto(path), table)


if __name__ == "__main__":
    unittest.main()
