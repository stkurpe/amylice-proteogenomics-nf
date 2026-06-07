from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .fasta import read_fasta


@dataclass(frozen=True)
class Check:
    name: str
    ok: bool
    message: str
    path: str = ""


REQUIRED_OUTPUTS = {
    "clean_ids": "clean_ids.txt",
    "protein": "protein.fasta",
    "frameshift": "frameshift_unique.fasta",
    "combined": "combine_proteome.fasta",
    "nonsense": "nonsense_candidates.txt",
}


def _line_count(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for line in handle if line.strip())


def _fasta_count(path: Path) -> int:
    return len(read_fasta(path))


def verify_proteome_outputs(proteo_dir: str | Path) -> list[Check]:
    root = Path(proteo_dir)
    checks: list[Check] = []

    for name, filename in REQUIRED_OUTPUTS.items():
        path = root / filename
        checks.append(Check(f"exists:{name}", path.is_file(), "required output exists" if path.is_file() else "missing required output", str(path)))
        if path.is_file():
            checks.append(Check(f"nonempty:{name}", path.stat().st_size > 0, "required output is non-empty" if path.stat().st_size > 0 else "required output is empty", str(path)))

    clean_ids = root / REQUIRED_OUTPUTS["clean_ids"]
    if clean_ids.is_file():
        count = _line_count(clean_ids)
        checks.append(Check("content:clean_ids", count > 0, f"clean transcript ids: {count}", str(clean_ids)))

    fasta_files = {
        "protein": root / REQUIRED_OUTPUTS["protein"],
        "frameshift": root / REQUIRED_OUTPUTS["frameshift"],
        "combined": root / REQUIRED_OUTPUTS["combined"],
    }
    fasta_counts: dict[str, int] = {}
    for name, path in fasta_files.items():
        if path.is_file() and path.stat().st_size > 0:
            count = _fasta_count(path)
            fasta_counts[name] = count
            checks.append(Check(f"content:{name}", count > 0, f"FASTA records: {count}", str(path)))

    if "protein" in fasta_counts and "combined" in fasta_counts:
        ok = fasta_counts["combined"] >= fasta_counts["protein"]
        checks.append(Check("consistency:combined_contains_protein", ok, f"combined={fasta_counts['combined']} protein={fasta_counts['protein']}", str(fasta_files["combined"])))

    if "frameshift" in fasta_counts and "combined" in fasta_counts:
        ok = fasta_counts["combined"] >= fasta_counts["frameshift"]
        checks.append(Check("consistency:combined_contains_frameshift", ok, f"combined={fasta_counts['combined']} frameshift={fasta_counts['frameshift']}", str(fasta_files["combined"])))

    return checks


def write_verification_report(proteo_dir: str | Path, output_path: str | Path) -> bool:
    checks = verify_proteome_outputs(proteo_dir)
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        handle.write("status\tcheck\tpath\tmessage\n")
        for check in checks:
            handle.write(f"{'OK' if check.ok else 'ATTENTION'}\t{check.name}\t{check.path}\t{check.message}\n")
    return all(check.ok for check in checks)
