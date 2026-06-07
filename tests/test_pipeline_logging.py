#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
LOGGING_SCRIPT = PROJECT_DIR / "pipeline_scripts" / "pipeline_logging.sh"


class PipelineLoggingTests(unittest.TestCase):
    def test_logging_script_exists_and_parses(self) -> None:
        self.assertTrue(LOGGING_SCRIPT.is_file(), f"Missing logging script: {LOGGING_SCRIPT}")
        result = subprocess.run(
            ["bash", "-n", str(LOGGING_SCRIPT)],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_logging_script_writes_log_events_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            driver = tmp_path / "driver.sh"
            log_dir = tmp_path / "logs"
            manifest = tmp_path / "manifest.txt"
            driver.write_text(
                textwrap.dedent(
                    f"""\
                    #!/bin/bash
                    set -Eeuo pipefail
                    source "{LOGGING_SCRIPT}"
                    init_pipeline_logging "SRR_TEST" "{log_dir}" "unit"
                    log_info "hello world"
                    run_logged "true command" bash -c 'printf "payload\\\\n"'
                    write_pipeline_manifest "{manifest}" "custom_key=custom_value"
                    """
                )
            )
            result = subprocess.run(["bash", str(driver)], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(result.returncode, 0, result.stderr)

            log_file = log_dir / "unit.log"
            event_file = log_dir / "unit.events.tsv"
            self.assertTrue(log_file.is_file())
            self.assertTrue(event_file.is_file())
            self.assertTrue(manifest.is_file())

            log_text = log_file.read_text()
            events_text = event_file.read_text()
            manifest_text = manifest.read_text()

            self.assertIn("hello world", log_text)
            self.assertIn("payload", log_text)
            self.assertIn("SRR_TEST\tunit\tINFO\thello world", events_text)
            self.assertIn("custom_key=custom_value", manifest_text)


if __name__ == "__main__":
    unittest.main(verbosity=2)
