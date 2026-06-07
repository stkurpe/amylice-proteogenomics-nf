#!/usr/bin/env python3
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from proteome_pipeline.fasta import FastaRecord, clean_stop_codons, combine_fastas, read_fasta, write_fasta
from proteome_pipeline.frameshift import generate_frameshift_records
from proteome_pipeline.genetics import reverse_complement, translate_until_stop
from proteome_pipeline.gtf import clean_transcript_ids
from proteome_pipeline.verify import write_verification_report


class FastaCleaningTests(unittest.TestCase):
    def test_trims_single_stop_and_drops_multiple_stops(self) -> None:
        cleaned = clean_stop_codons([
            FastaRecord(">ok", "MAAA.BBB"),
            FastaRecord(">bad", "MAA.BB.CC"),
            FastaRecord(">plain", "MPEPTIDE"),
        ])
        self.assertEqual([r.id for r in cleaned], ["ok", "plain"])
        self.assertEqual(cleaned[0].sequence, "MAAA")

    def test_combine_deduplicates_by_id_and_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            a = Path(tmp) / "a.fa"
            b = Path(tmp) / "b.fa"
            write_fasta([FastaRecord(">p1", "AAA")], a)
            write_fasta([FastaRecord(">p1", "AAA"), FastaRecord(">p2", "BBB")], b)
            self.assertEqual([r.id for r in combine_fastas([a, b])], ["p1", "p2"])


class GeneticsTests(unittest.TestCase):
    def test_translate_until_stop(self) -> None:
        self.assertEqual(translate_until_stop("ATGAAATAAGGG"), "MK")

    def test_reverse_complement(self) -> None:
        self.assertEqual(reverse_complement("ATGCCN"), "NGGCAT")


class CleanIdTests(unittest.TestCase):
    def test_clean_ids_remove_mito_short_cds_and_unexpressed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gtf = root / "test.gtf"
            abundance = root / "abundance.tsv"
            out = root / "clean_ids.txt"
            gtf.write_text(
                "\n".join([
                    'chr1\tx\ttranscript\t1\t30\t.\t+\t.\tgene_id "G1"; transcript_id "TX1.1"; gene_name "G1";',
                    'chr1\tx\texon\t1\t30\t.\t+\t.\tgene_id "G1"; transcript_id "TX1.1"; gene_name "G1";',
                    'chr1\tx\tCDS\t1\t30\t.\t+\t0\tgene_id "G1"; transcript_id "TX1.1"; gene_name "G1";',
                    'chr1\tx\texon\t40\t45\t.\t+\t.\tgene_id "G2"; transcript_id "SHORT.1"; gene_name "G2";',
                    'chr1\tx\tCDS\t40\t45\t.\t+\t0\tgene_id "G2"; transcript_id "SHORT.1"; gene_name "G2";',
                    'chrM\tx\texon\t1\t30\t.\t+\t.\tgene_id "MT"; transcript_id "MTTX.1"; gene_name "MT";',
                    'chrM\tx\tCDS\t1\t30\t.\t+\t0\tgene_id "MT"; transcript_id "MTTX.1"; gene_name "MT";',
                    "",
                ])
            )
            abundance.write_text("target_id\tlength\teff_length\test_counts\ttpm\nTX1.1|x\t1\t1\t1\t2\nSHORT.1|x\t1\t1\t1\t9\nMTTX.1|x\t1\t1\t1\t9\n")
            ids = clean_transcript_ids(gtf, abundance, out, min_tpm=1, min_cds_bp=10)
            self.assertEqual(ids, ["TX1.1"])
            self.assertEqual(out.read_text().strip(), "TX1.1")


class FrameshiftTests(unittest.TestCase):
    def test_frameshift_uses_clean_ids_and_deduplicates_per_gene(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gtf = root / "test.gtf"
            ids = root / "clean_ids.txt"
            vcf = root / "test.vcf"
            ids.write_text("TX1.1\nTX2.1\n")
            gtf.write_text(
                "\n".join([
                    'chr1\tx\texon\t1\t24\t.\t+\t.\tgene_id "G1"; transcript_id "TX1.1"; gene_name "GENE1";',
                    'chr1\tx\tCDS\t1\t24\t.\t+\t0\tgene_id "G1"; transcript_id "TX1.1"; gene_name "GENE1";',
                    'chr1\tx\texon\t1\t24\t.\t+\t.\tgene_id "G1"; transcript_id "TX2.1"; gene_name "GENE1";',
                    'chr1\tx\tCDS\t1\t24\t.\t+\t0\tgene_id "G1"; transcript_id "TX2.1"; gene_name "GENE1";',
                    'chr1\tx\texon\t1\t24\t.\t+\t.\tgene_id "G2"; transcript_id "BLOCKED.1"; gene_name "G2";',
                    'chr1\tx\tCDS\t1\t24\t.\t+\t0\tgene_id "G2"; transcript_id "BLOCKED.1"; gene_name "G2";',
                    "",
                ])
            )
            vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t4\t.\tA\tAT\t.\tPASS\t.\n")

            genome = {"chr1": "ATGAAAAAAAAAAAAAAAAAAAAA"}

            def fetch(chrom: str, start: int, end: int) -> str:
                return genome[chrom][start - 1:end]

            records = generate_frameshift_records(vcf, gtf, ids, fetch, min_aa=6)
            self.assertEqual(len(records), 1)
            self.assertTrue(records[0].header.startswith(">FRAMESHIFT_GENE1_TX1.1"))
            self.assertGreaterEqual(len(records[0].sequence), 6)


class VerificationTests(unittest.TestCase):
    def test_verification_reports_attention_for_missing_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "clean_ids.txt").write_text("TX1\n")
            write_fasta([FastaRecord(">p1", "AAA")], root / "protein.fasta")
            report = root / "verification.tsv"
            ok = write_verification_report(root, report)
            self.assertFalse(ok)
            text = report.read_text()
            self.assertIn("ATTENTION\texists:frameshift", text)
            self.assertIn("ATTENTION\texists:combined", text)

    def test_verification_passes_complete_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "clean_ids.txt").write_text("TX1\n")
            (root / "nonsense_candidates.txt").write_text("TRANSCRIPT_ID\tLEN_H1\tLEN_H2\tLOSS\tSTATUS\n")
            write_fasta([FastaRecord(">p1", "AAA")], root / "protein.fasta")
            write_fasta([FastaRecord(">fs1", "BBBBBB")], root / "frameshift_unique.fasta")
            write_fasta([FastaRecord(">p1", "AAA"), FastaRecord(">fs1", "BBBBBB")], root / "combine_proteome.fasta")
            report = root / "verification.tsv"
            self.assertTrue(write_verification_report(root, report))
            self.assertNotIn("ATTENTION", report.read_text())


if __name__ == "__main__":
    unittest.main(verbosity=2)
