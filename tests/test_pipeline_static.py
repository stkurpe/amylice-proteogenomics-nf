#!/usr/bin/env python3
"""Static tests for the /home/codex shell pipeline.

These tests intentionally do not run Docker, AWS, downloads, sudo, or the
pipeline itself. They validate script syntax and configuration consistency,
then report known risky patterns as explicit findings.
"""

from __future__ import annotations

import re
import subprocess
import unittest
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
PIPELINE_DIR = PROJECT_DIR / "pipeline_scripts"

PROJECT_FILES = [
    PROJECT_DIR / "run_full_pipeline.sh",
    PROJECT_DIR / "run_complete_processing_pipeline.sh",
    PROJECT_DIR / "run_prepared_proteomes_from_bioinfo.sh",
    PROJECT_DIR / "run_amyloid_predictors.sh",
    PROJECT_DIR / "samples.txt.example",
    PIPELINE_DIR / "00_prepare_refs.sh",
    PIPELINE_DIR / "proteom_generation_analysis.sh",
    PIPELINE_DIR / "pipeline_logging.sh",
    PIPELINE_DIR / "download_public_inputs.sh",
    PIPELINE_DIR / "run_local_proteome.sh",
    PIPELINE_DIR / "run_public_sra_to_proteome.sh",
    PIPELINE_DIR / "common.sh",
    PIPELINE_DIR / "01_download.sh",
    PIPELINE_DIR / "02_qc.sh",
    PIPELINE_DIR / "03_expression.sh",
    PIPELINE_DIR / "04_variants.sh",
    PIPELINE_DIR / "05_hla_consensus.sh",
    PIPELINE_DIR / "05_proteo.sh",
    PIPELINE_DIR / "06_prot.xh",
    PIPELINE_DIR / "config.sh",
    PIPELINE_DIR / "config_s3.sh",
    PIPELINE_DIR / "finish_pipeline_combined.sh",
    PIPELINE_DIR / "finish_pipeline_v2_autoheal.sh",
    PIPELINE_DIR / "finish_pipeline_v2_clean.sh",
    PIPELINE_DIR / "run_gatk_step.sh",
]

SHELL_FILES = [path for path in PROJECT_FILES if path.suffix in {".sh", ".xh"}]
CONFIG_FILES = [PIPELINE_DIR / "config.sh", PIPELINE_DIR / "config_s3.sh"]
CORE_STEP_FILES = [
    PIPELINE_DIR / "01_download.sh",
    PIPELINE_DIR / "02_qc.sh",
    PIPELINE_DIR / "03_expression.sh",
    PIPELINE_DIR / "04_variants.sh",
    PIPELINE_DIR / "05_proteo.sh",
]

EXPECTED_RISK_PATTERNS = {
    "destructive_delete": re.compile(r"\brm\s+-(?:[A-Za-z]*r[A-Za-z]*f|[A-Za-z]*f[A-Za-z]*r)\b|\brm\s+-f\b"),
    "sudo_or_system_change": re.compile(r"\bsudo\b|\bapt-get\b|\busermod\b|\bsystemctl\b"),
    "network_or_remote_io": re.compile(r"\bwget\b|\bcurl\b|\baws\s+s3\s+cp\b|\bdocker\s+run\b"),
    "inline_rewrite": re.compile(r"\bcat\s+<<|\bsed\s+-i\b|\bchmod\b"),
}


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


class FilePresenceTests(unittest.TestCase):
    def test_project_root_exists(self) -> None:
        self.assertTrue(PROJECT_DIR.is_dir(), f"Missing project dir: {PROJECT_DIR}")

    def test_expected_files_exist(self) -> None:
        missing = [str(path) for path in PROJECT_FILES if not path.exists()]
        self.assertEqual(missing, [], f"Missing expected project files: {missing}")

    def test_tests_are_outside_pipeline_dir(self) -> None:
        tests_dir = PROJECT_DIR / "tests"
        self.assertTrue(tests_dir.is_dir(), f"Missing tests dir: {tests_dir}")
        self.assertFalse(str(tests_dir).startswith(str(PIPELINE_DIR)))


class ShellSyntaxTests(unittest.TestCase):
    def test_shell_scripts_parse_with_bash_n(self) -> None:
        failures: list[str] = []
        for path in SHELL_FILES:
            result = subprocess.run(
                ["bash", "-n", str(path)],
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            if result.returncode != 0:
                failures.append(f"{path}: {result.stderr.strip()}")
        self.assertEqual(failures, [], "\n".join(failures))

    def test_scripts_have_shebang_or_are_intentional_fragments(self) -> None:
        failures: list[str] = []
        for path in SHELL_FILES:
            if path.name in {"06_prot.xh", "run_gatk_step.sh"}:
                continue
            first_nonempty = next(
                (line for line in read_text(path).splitlines() if line.strip()),
                "",
            )
            if first_nonempty not in {"#!/bin/bash", "#!/usr/bin/env bash"}:
                failures.append(str(path))
        self.assertEqual(failures, [], f"Scripts without bash shebang: {failures}")


class ConfigConsistencyTests(unittest.TestCase):
    def test_configs_export_expected_variables(self) -> None:
        expected = {
            "BASE_DIR",
            "REF_DIR",
            "THREADS",
            "DOCKER_OPTS",
            "REF_GENOME",
            "REF_GTF",
            "REF_VCF",
            "REF_TRANSCRIPT",
        }
        config_text = "\n".join(read_text(path) for path in CONFIG_FILES if path.exists())
        missing = [name for name in sorted(expected) if f"export {name}=" not in config_text]
        self.assertEqual(missing, [], f"Missing config exports: {missing}")

    def test_s3_bucket_is_configured_without_credentials(self) -> None:
        config_s3 = read_text(PIPELINE_DIR / "config_s3.sh")
        self.assertIn("S3_BUCKET", config_s3)
        forbidden = re.compile(r"AWS_SECRET|AWS_ACCESS_KEY|SECRET_ACCESS_KEY|TOKEN|PASSWORD", re.I)
        self.assertIsNone(forbidden.search(config_s3), "config_s3.sh appears to contain credential material")

    def test_core_steps_source_shared_helpers(self) -> None:
        failures: list[str] = []
        for path in CORE_STEP_FILES:
            text = read_text(path)
            if "config.sh" not in text:
                failures.append(f"{path}: missing config.sh source")
            if "common.sh" not in text:
                failures.append(f"{path}: missing common.sh source")
            if "require_sample_id" not in text:
                failures.append(f"{path}: missing sample id validation")
        self.assertEqual(failures, [], "\n".join(failures))


class UploadPolicyTests(unittest.TestCase):
    def test_local_proteome_uploads_intermediate_analysis_outputs_before_cleanup(self) -> None:
        text = read_text(PIPELINE_DIR / "run_local_proteome.sh")

        self.assertIn("upload_analysis_outputs()", text)
        for directory in ("results_qc", "results_expression", "results_star", "results_gatk", "logs"):
            self.assertIn(f"${{work_dir}}/{directory}", text)
            self.assertIn(f"${{S3_BUCKET}}/${{sample_id}}/{directory}", text)

        upload_pos = text.index("upload_results \"$SAMPLE_ID\" \"$PROTEO_DIR\"")
        cleanup_pos = text.index("cleanup_after_success \"$WORK_DIR\"")
        self.assertLess(upload_pos, cleanup_pos, "analysis outputs must upload before local cleanup")

    def test_upload_success_marker_is_uploaded_after_all_results(self) -> None:
        text = read_text(PIPELINE_DIR / "run_local_proteome.sh")
        upload_results = text[text.index("upload_results()") : text.index("cleanup_after_success()")]

        marker_pos = upload_results.rindex("UPLOAD_SUCCESS")
        for required in (
            "upload_analysis_outputs",
            "protein.fasta",
            "combine_proteome.fasta",
            "verification_report.tsv",
            "PIPELINE_LOG_FILE",
        ):
            self.assertLess(
                upload_results.index(required),
                marker_pos,
                f"{required} should be uploaded before UPLOAD_SUCCESS marker",
            )


class SampleFileTests(unittest.TestCase):
    def test_samples_are_nonempty_comments_or_safe_ids(self) -> None:
        samples = PROJECT_DIR / "samples.txt.example"
        unsafe: list[str] = []
        usable = 0
        sample_re = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")

        for number, raw in enumerate(read_text(samples).splitlines(), start=1):
            value = raw.strip()
            if not value or value.startswith("#"):
                continue
            usable += 1
            if not sample_re.fullmatch(value) or "/" in value or ".." in value:
                unsafe.append(f"{number}:{value}")

        self.assertGreater(usable, 0, "samples.txt has no usable sample ids")
        self.assertEqual(unsafe, [], f"Unsafe sample ids: {unsafe}")


class StaticRiskInventoryTests(unittest.TestCase):
    def test_risky_patterns_are_inventory_not_surprises(self) -> None:
        findings: dict[str, list[str]] = {name: [] for name in EXPECTED_RISK_PATTERNS}
        for path in PROJECT_FILES:
            if not path.exists() or path.is_dir():
                continue
            text = read_text(path)
            for line_no, line in enumerate(text.splitlines(), start=1):
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                for name, pattern in EXPECTED_RISK_PATTERNS.items():
                    if pattern.search(stripped):
                        findings[name].append(f"{path}:{line_no}: {stripped}")

        for name, entries in findings.items():
            print(f"\n[{name}] {len(entries)} finding(s)")
            for entry in entries:
                print(entry)

        self.assertTrue(
            any(findings.values()),
            "No risky patterns found; update this inventory test after safety refactors.",
        )

    def test_no_obvious_secret_files_are_part_of_project_files(self) -> None:
        secret_names = {".env", "id_rsa", "id_ed25519", "credentials", "authorized_keys"}
        included = [str(path) for path in PROJECT_FILES if path.name in secret_names]
        self.assertEqual(included, [], f"Secret-like files should not be test fixtures: {included}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
