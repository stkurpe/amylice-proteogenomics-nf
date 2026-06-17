from __future__ import annotations

import re
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
NEXTFLOW_FILES = [
    PROJECT_DIR / "nextflow.config",
    PROJECT_DIR / "nextflow" / "main.nf",
    PROJECT_DIR / "nextflow" / "modules" / "proteome.nf",
]


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_config_profiles_are_declared() -> None:
    config = read(PROJECT_DIR / "nextflow.config")
    for profile in ["local", "test", "aws_reference", "docker"]:
        assert re.search(rf"\b{profile}\s*\{{", config), f"missing {profile} profile"


def test_process_container_is_set_for_docker_profile() -> None:
    config = read(PROJECT_DIR / "nextflow.config")
    assert "process.container" in config
    assert "amyloid-proteome-nextflow:local" in config


def test_required_params_have_defaults_or_validation() -> None:
    config = read(PROJECT_DIR / "nextflow.config")
    main = read(PROJECT_DIR / "nextflow" / "main.nf")

    assert "outdir" in config
    for param in ["min_tpm", "min_cds_bp", "min_frameshift_aa"]:
        assert re.search(rf"\b{param}\s*=", config), f"missing default for {param}"
        assert f"params.{param}" in main, f"missing validation/use for {param}"


def test_help_mode_is_supported_without_running_processes() -> None:
    main = read(PROJECT_DIR / "nextflow" / "main.nf")
    assert "params.help" in main
    assert "return" in main
    assert "--samples PATH" in main


def test_no_home_codex_hardcoding_in_new_nextflow_files() -> None:
    offenders: list[str] = []
    for path in NEXTFLOW_FILES:
        for line_no, line in enumerate(read(path).splitlines(), start=1):
            if "/home/codex" in line:
                offenders.append(f"{path.relative_to(PROJECT_DIR)}:{line_no}:{line}")
    assert offenders == []


def test_proteome_container_declares_required_tools() -> None:
    dockerfile = read(PROJECT_DIR / "docker" / "proteome-nextflow" / "Dockerfile")
    for token in ["python:3.11-slim", "bcftools", "samtools", "gffread"]:
        assert token in dockerfile
