from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from amylogram_py.fasta import iter_fasta_records, read_fasta


class FastaTests(unittest.TestCase):
    def test_iter_fasta_records_reports_skip_reasons(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.fasta"
            path.write_text(">ok description\nVQIVYK\n>short\nACD\n>dirty\n***zzz\n", encoding="utf-8")

            records = list(iter_fasta_records(path))
            self.assertEqual([record.sequence_id for record in records], ["ok", "short", "dirty"])
            self.assertEqual(records[0].status, "ok")
            self.assertEqual(records[1].reason, "shorter_than_6")
            self.assertEqual(records[2].reason, "empty_after_cleaning")

    def test_read_fasta_preserves_ok_only_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.fasta"
            path.write_text(">ok\nVQIVYK\n>short\nACD\n", encoding="utf-8")
            self.assertEqual(list(read_fasta(path)), [("ok", "VQIVYK")])


if __name__ == "__main__":
    unittest.main()
