from __future__ import annotations

import csv
import importlib
import unittest
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
NEXTFLOW_DIR = PROJECT_DIR / "nextflow"
SAMPLES_CSV = NEXTFLOW_DIR / "samples.prepared.example.csv"
FIXTURE_DIR = PROJECT_DIR / "tests_nextflow" / "fixtures" / "prepared_minimal"


class Phase0ContractTests(unittest.TestCase):
    def test_phase0_nextflow_scaffold_files_exist(self) -> None:
        expected = [
            NEXTFLOW_DIR / "main.nf",
            PROJECT_DIR / "nextflow.config",
            NEXTFLOW_DIR / "modules" / "proteome.nf",
            SAMPLES_CSV,
            PROJECT_DIR / "docker" / "proteome-nextflow" / "Dockerfile",
            PROJECT_DIR / "tests_nextflow" / "run_prepared_smoke.sh",
            PROJECT_DIR / "tests_nextflow" / "test_resume.sh",
        ]
        missing = [str(path.relative_to(PROJECT_DIR)) for path in expected if not path.is_file()]
        self.assertEqual(missing, [])

    def test_phase1_prepared_processes_are_declared(self) -> None:
        module_text = (NEXTFLOW_DIR / "modules" / "proteome.nf").read_text(encoding="utf-8")
        for process_name in [
            "CLEAN_IDS",
            "PREPARE_SNP_VCF",
            "CONSENSUS_H1",
            "CONSENSUS_H2",
            "TRANSLATE_H1",
            "TRANSLATE_H2",
            "CLEAN_SNP_PROTEINS",
            "FRAMESHIFT_PROTEOME",
            "COMBINE_PROTEOME",
            "VERIFY_PROTEOME",
        ]:
            self.assertIn(f"process {process_name}", module_text)

    def test_prepared_smoke_script_runs_nextflow_and_checks_outputs(self) -> None:
        smoke = (PROJECT_DIR / "tests_nextflow" / "run_prepared_smoke.sh").read_text(encoding="utf-8")
        self.assertIn("docker build -t amyloid-proteome-nextflow:local", smoke)
        self.assertIn("docker build -t amylogram-py-nextflow:local", smoke)
        self.assertIn("nextflow run nextflow/main.nf", smoke)
        self.assertIn("-profile test,docker", smoke)
        for filename in [
            "clean_ids.txt",
            "protein.fasta",
            "frameshift_unique.fasta",
            "combine_proteome.fasta",
            "nonsense_candidates.txt",
            "verification_report.tsv",
            "manifest.txt",
        ]:
            self.assertIn(filename, smoke)
        self.assertIn("grep -q OK", smoke)
        self.assertIn("grep -c '^>'", smoke)
        self.assertIn("results_amyloid", smoke)
        self.assertIn("results_protein_features", smoke)

    def test_phase5_resume_script_checks_cache_and_reports(self) -> None:
        resume = (PROJECT_DIR / "tests_nextflow" / "test_resume.sh").read_text(encoding="utf-8")
        self.assertIn("-resume", resume)
        self.assertIn("-with-trace", resume)
        self.assertIn("-with-report", resume)
        self.assertIn("-with-timeline", resume)
        self.assertIn("grep -q CACHED", resume)
        self.assertIn("results_amyloid", resume)
        self.assertIn("results_protein_features", resume)

    def test_samples_manifest_is_valid_and_selects_phase0_samples(self) -> None:
        with SAMPLES_CSV.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))

        self.assertTrue(rows)
        self.assertEqual(
            set(rows[0]),
            {"sample_id", "kind", "abundance", "vcf", "expected_proteins_dir", "notes"},
        )

        by_id = {row["sample_id"]: row for row in rows}
        self.assertIn("fixture_minimal", by_id)
        self.assertIn("SRR32060234", by_id)
        self.assertEqual(by_id["fixture_minimal"]["kind"], "fixture")
        self.assertEqual(by_id["SRR32060234"]["kind"], "validated_sra")

    def test_local_fixture_inputs_exist_and_have_expected_shapes(self) -> None:
        abundance = FIXTURE_DIR / "abundance.tsv"
        vcf = FIXTURE_DIR / "variants_filtered.vcf"
        genome = FIXTURE_DIR / "genome.fa"
        gtf = FIXTURE_DIR / "annotation.gtf"

        for path in [abundance, vcf, genome, gtf]:
            self.assertTrue(path.is_file(), path)
            self.assertGreater(path.stat().st_size, 0, path)

        header = abundance.read_text(encoding="utf-8").splitlines()[0].split("\t")
        self.assertTrue({"target_id", "length", "eff_length", "est_counts", "tpm"}.issubset(header))

        vcf_text = vcf.read_text(encoding="utf-8")
        self.assertTrue(vcf_text.startswith("##fileformat=VCF"))
        self.assertIn(
            "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfixture_minimal\n",
            vcf_text,
        )

        gtf_lines = [line for line in gtf.read_text(encoding="utf-8").splitlines() if line.strip()]
        self.assertTrue(any("\ttranscript\t" in line for line in gtf_lines))
        self.assertTrue(any("\tCDS\t" in line for line in gtf_lines))

    def test_proteome_pipeline_cli_imports(self) -> None:
        module = importlib.import_module("proteome_pipeline.cli")
        self.assertTrue(hasattr(module, "build_parser"))
        self.assertTrue(hasattr(module, "main"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
