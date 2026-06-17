from __future__ import annotations

import subprocess
from pathlib import Path

from tests_nextflow.compare_proteome_outputs import compare_outputs


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_DIR / "tests_nextflow" / "compare_proteome_outputs.py"
BASELINE = PROJECT_DIR / "tests_nextflow" / "fixtures" / "prepared_minimal" / "expected" / "results_proteins"


def write_fixture(root: Path, protein_seq: str = "MAIVMGR") -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / "clean_ids.txt").write_text("ENST000001\n", encoding="utf-8")
    (root / "protein.fasta").write_text(f">H1_ENST000001\n{protein_seq}\n", encoding="utf-8")
    (root / "frameshift_unique.fasta").write_text(">FRAMESHIFT_GENE_FIXTURE_ENST000001\nMAIVNGPLKGCPI\n", encoding="utf-8")
    (root / "combine_proteome.fasta").write_text(
        f">H1_ENST000001\n{protein_seq}\n>FRAMESHIFT_GENE_FIXTURE_ENST000001\nMAIVNGPLKGCPI\n",
        encoding="utf-8",
    )
    (root / "verification_report.tsv").write_text(
        "status\tcheck\tpath\tmessage\nOK\texists:combined\tcombine_proteome.fasta\trequired output exists\n",
        encoding="utf-8",
    )


def test_compare_outputs_accepts_matching_normalized_fastas(tmp_path: Path) -> None:
    expected_root = tmp_path / "expected"
    nf_root = tmp_path / "nf"
    write_fixture(expected_root)
    write_fixture(nf_root)

    comparisons = compare_outputs(expected_root, nf_root)
    assert {item.filename for item in comparisons} == {
        "protein.fasta",
        "frameshift_unique.fasta",
        "combine_proteome.fasta",
    }


def test_compare_outputs_rejects_sequence_mismatch(tmp_path: Path) -> None:
    expected_root = tmp_path / "expected"
    nf_root = tmp_path / "nf"
    write_fixture(expected_root)
    write_fixture(nf_root, protein_seq="MTIVMGR")

    result = subprocess.run(
        [
            str(SCRIPT),
            "--expected-output",
            str(expected_root),
            "--nf-output",
            str(nf_root),
        ],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert result.returncode == 1
    assert "sequences differ" in result.stderr


def test_baseline_fixture_is_comparable_to_itself() -> None:
    comparisons = compare_outputs(BASELINE, BASELINE)
    assert sum(item.record_count for item in comparisons) > 0
