from __future__ import annotations

from pathlib import Path

from .fasta import read_fasta


def write_nonsense_report(protein_fasta: str | Path, output_path: str | Path) -> int:
    data: dict[str, dict[str, int]] = {}
    for record in read_fasta(protein_fasta):
        clean = record.id
        if "_" not in clean:
            continue
        haplo, transcript = clean.split("_", 1)
        data.setdefault(transcript, {})[haplo] = len(record.sequence)

    count = 0
    with Path(output_path).open("w") as out:
        out.write("TRANSCRIPT_ID\tLEN_H1\tLEN_H2\tLOSS\tSTATUS\n")
        for transcript, lens in sorted(data.items()):
            if "H1" not in lens or "H2" not in lens:
                continue
            l1, l2 = lens["H1"], lens["H2"]
            if l1 == l2:
                continue
            max_len = max(l1, l2)
            if max_len and abs(l1 - l2) > 10 and abs(l1 - l2) / max_len > 0.05:
                status = "H2_TRUNCATED" if l1 > l2 else "H1_TRUNCATED"
                out.write(f"{transcript}\t{l1}\t{l2}\t-{abs(l1 - l2)}\t{status}\n")
                count += 1
    return count
