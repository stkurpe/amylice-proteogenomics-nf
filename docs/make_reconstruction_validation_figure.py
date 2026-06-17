#!/usr/bin/env python3
"""Build cohort-level reconstruction validation graphics from Amylice outputs.

The script intentionally uses only the Python standard library so the figure can
be regenerated on lightweight systems without plotting dependencies.
"""

from __future__ import annotations

import csv
import html
import math
import statistics
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
WORK_ROOT = REPO_ROOT / "ValCalculation" / "work"
OUT_DIR = Path(__file__).resolve().parent
METRICS_CSV = OUT_DIR / "reconstruction_validation_metrics.csv"
FIGURE_SVG = OUT_DIR / "figure_reconstruction_validation.svg"


def read_summary(path: Path) -> dict[str, str]:
    rows: dict[str, str] = {}
    if not path.exists():
        return rows
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows[row["metric"]] = row["value"]
    return rows


def read_status(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts)))
                name = line[1:]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        records.append((name, "".join(parts)))
    return records


def sequence_class(header: str) -> str:
    if header.startswith("FRAMESHIFT_"):
        return "Frameshift"
    if header.startswith("H1_"):
        return "Haplotype 1"
    if header.startswith("H2_"):
        return "Haplotype 2"
    return "Other"


def int_value(summary: dict[str, str], key: str) -> int:
    value = summary.get(key, "0")
    try:
        return int(float(value))
    except ValueError:
        return 0


def collect_metrics() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    samples: list[dict[str, object]] = []
    lengths: list[dict[str, object]] = []
    for sample_dir in sorted(WORK_ROOT.glob("SRR*")):
        fasta = sample_dir / "input" / "combine_proteome.fasta"
        if not fasta.exists():
            continue
        sample_id = sample_dir.name
        records = fasta_records(fasta)
        class_counts = {"Haplotype 1": 0, "Haplotype 2": 0, "Frameshift": 0, "Other": 0}
        sample_lengths: list[int] = []
        for header, seq in records:
            cls = sequence_class(header)
            class_counts[cls] += 1
            length = len(seq)
            sample_lengths.append(length)
            lengths.append({"sample": sample_id, "class": cls, "length": length})

        summary = read_summary(sample_dir / "results" / "amyloid_predictors_summary.tsv")
        status_rows = read_status(sample_dir / "results" / "amyloid_predictors_status.tsv")
        status_values = {row.get("status", "") for row in status_rows}
        completed = "FAIL" not in status_values and bool(status_rows)
        input_records = int_value(summary, "input_records") or len(records)
        combined_rows = int_value(summary, "combined_rows")
        amypred_rows = int_value(summary, "amypred_rows")
        retained_pct = (combined_rows / input_records * 100.0) if input_records else 0.0

        samples.append(
            {
                "sample": sample_id,
                "total_sequences": len(records),
                "haplotype_1": class_counts["Haplotype 1"],
                "haplotype_2": class_counts["Haplotype 2"],
                "frameshift": class_counts["Frameshift"],
                "other": class_counts["Other"],
                "altered_sequences": class_counts["Haplotype 1"]
                + class_counts["Haplotype 2"]
                + class_counts["Frameshift"],
                "median_length": statistics.median(sample_lengths) if sample_lengths else 0,
                "mean_length": statistics.mean(sample_lengths) if sample_lengths else 0,
                "input_records": input_records,
                "combined_rows": combined_rows,
                "amypred_rows": amypred_rows,
                "retained_pct": retained_pct,
                "completed": completed,
            }
        )
    return samples, lengths


def write_metrics(samples: list[dict[str, object]]) -> None:
    fields = [
        "sample",
        "total_sequences",
        "haplotype_1",
        "haplotype_2",
        "frameshift",
        "other",
        "altered_sequences",
        "median_length",
        "mean_length",
        "input_records",
        "combined_rows",
        "amypred_rows",
        "retained_pct",
        "completed",
    ]
    with METRICS_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(samples)


def fmt_int(value: float | int) -> str:
    return f"{int(round(float(value))):,}"


def svg_text(x: float, y: float, text: str, size: int = 14, weight: int = 400, fill: str = "#1f2933", anchor: str = "start") -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" font-weight="{weight}" '
        f'fill="{fill}" text-anchor="{anchor}">{html.escape(text)}</text>'
    )


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str = "none", rx: float = 8, opacity: float = 1.0) -> str:
    return (
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx:.1f}" '
        f'fill="{fill}" stroke="{stroke}" opacity="{opacity:.3f}"/>'
    )


def line(x1: float, y1: float, x2: float, y2: float, stroke: str = "#c7d0da", width: float = 1.0) -> str:
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="{width:.1f}"/>'


def panel_label(x: float, y: float, label: str, title: str) -> list[str]:
    return [
        svg_text(x, y, label, size=13, weight=800, fill="#64748b"),
        svg_text(x + 28, y, title, size=18, weight=750, fill="#111827"),
    ]


def nice_max(value: float) -> float:
    if value <= 0:
        return 1
    exp = math.floor(math.log10(value))
    base = 10**exp
    for mult in (1, 2, 5, 10):
        if value <= mult * base:
            return mult * base
    return 10 * base


def draw_summary_cards(samples: list[dict[str, object]], x: int, y: int, w: int, h: int) -> list[str]:
    total_samples = len(samples)
    completed = sum(1 for row in samples if row["completed"])
    failed = total_samples - completed
    seqs = [int(row["total_sequences"]) for row in samples]
    frameshifts = [int(row["frameshift"]) for row in samples]
    retained = [float(row["retained_pct"]) for row in samples]
    cards = [
        ("Samples processed", fmt_int(total_samples), f"{completed} complete, {failed} incomplete"),
        ("Median sequences", fmt_int(statistics.median(seqs)), "reconstructed per sample"),
        ("Median frameshifts", fmt_int(statistics.median(frameshifts)), "altered protein products"),
        ("Median retained", f"{statistics.median(retained):.1f}%", "after prediction merge"),
    ]
    gap = 14
    card_w = (w - gap * 3) / 4
    out: list[str] = []
    for i, (title, value, sub) in enumerate(cards):
        cx = x + i * (card_w + gap)
        out.append(rect(cx, y, card_w, h, "#eef6ff", "#d7e4f2", rx=10))
        out.append(svg_text(cx + 16, y + 28, title, size=13, weight=700, fill="#475569"))
        out.append(svg_text(cx + 16, y + 68, value, size=30, weight=800, fill="#0f172a"))
        out.append(svg_text(cx + 16, y + 96, sub, size=13, fill="#64748b"))
    return out


def draw_bar_panel(samples: list[dict[str, object]], x: int, y: int, w: int, h: int) -> list[str]:
    out = panel_label(x, y, "A", "Reconstructed protein sequences per sample")
    plot_x, plot_y = x + 48, y + 32
    plot_w, plot_h = w - 72, h - 70
    max_v = nice_max(max(int(row["total_sequences"]) for row in samples))
    out.append(line(plot_x, plot_y + plot_h, plot_x + plot_w, plot_y + plot_h, "#94a3b8", 1.4))
    out.append(line(plot_x, plot_y, plot_x, plot_y + plot_h, "#94a3b8", 1.4))
    bar_gap = 8
    bar_w = (plot_w - bar_gap * (len(samples) - 1)) / len(samples)
    for i, row in enumerate(samples):
        value = int(row["total_sequences"])
        bh = value / max_v * plot_h
        bx = plot_x + i * (bar_w + bar_gap)
        by = plot_y + plot_h - bh
        out.append(rect(bx, by, bar_w, bh, "#5b8def", rx=3))
        out.append(svg_text(bx + bar_w / 2, plot_y + plot_h + 18, str(i + 1), size=10, fill="#475569", anchor="middle"))
    for frac in (0, 0.5, 1.0):
        gy = plot_y + plot_h - plot_h * frac
        out.append(line(plot_x, gy, plot_x + plot_w, gy, "#e2e8f0", 0.8))
        out.append(svg_text(plot_x - 8, gy + 4, fmt_int(max_v * frac), size=10, fill="#64748b", anchor="end"))
    out.append(svg_text(plot_x, y + h - 8, "Samples are ordered by SRR ID; x-axis labels use compact indices.", size=11, fill="#64748b"))
    return out


def draw_stacked_panel(samples: list[dict[str, object]], x: int, y: int, w: int, h: int) -> list[str]:
    out = panel_label(x, y, "B", "Proteoform classes in each reconstructed proteome")
    colors = {"Haplotype 1": "#7dd3fc", "Haplotype 2": "#38bdf8", "Frameshift": "#f59e0b", "Other": "#cbd5e1"}
    keys = [("haplotype_1", "Haplotype 1"), ("haplotype_2", "Haplotype 2"), ("frameshift", "Frameshift"), ("other", "Other")]
    plot_x, plot_y = x + 48, y + 32
    plot_w, plot_h = w - 72, h - 76
    max_v = nice_max(max(int(row["total_sequences"]) for row in samples))
    bar_gap = 8
    bar_w = (plot_w - bar_gap * (len(samples) - 1)) / len(samples)
    out.append(line(plot_x, plot_y + plot_h, plot_x + plot_w, plot_y + plot_h, "#94a3b8", 1.4))
    out.append(line(plot_x, plot_y, plot_x, plot_y + plot_h, "#94a3b8", 1.4))
    for i, row in enumerate(samples):
        bx = plot_x + i * (bar_w + bar_gap)
        current = plot_y + plot_h
        for key, label in keys:
            value = int(row[key])
            seg_h = value / max_v * plot_h
            current -= seg_h
            if seg_h > 0:
                out.append(rect(bx, current, bar_w, seg_h, colors[label], rx=2))
    lx = plot_x
    ly = y + h - 34
    for i, (_, label) in enumerate(keys[:-1]):
        out.append(rect(lx + i * 120, ly, 12, 12, colors[label], rx=2))
        out.append(svg_text(lx + 18 + i * 120, ly + 11, label, size=11, fill="#475569"))
    return out


def draw_length_panel(lengths: list[dict[str, object]], x: int, y: int, w: int, h: int) -> list[str]:
    out = panel_label(x, y, "C", "Protein length distribution")
    values = [int(row["length"]) for row in lengths if int(row["length"]) > 0]
    bins = [0, 25, 50, 100, 200, 400, 800, 1600, 3200, 6400]
    counts = [0 for _ in range(len(bins) - 1)]
    for value in values:
        for i in range(len(bins) - 1):
            if bins[i] <= value < bins[i + 1]:
                counts[i] += 1
                break
    plot_x, plot_y = x + 54, y + 34
    plot_w, plot_h = w - 86, h - 78
    max_c = nice_max(max(counts))
    bar_gap = 6
    bar_w = (plot_w - bar_gap * (len(counts) - 1)) / len(counts)
    out.append(line(plot_x, plot_y + plot_h, plot_x + plot_w, plot_y + plot_h, "#94a3b8", 1.4))
    out.append(line(plot_x, plot_y, plot_x, plot_y + plot_h, "#94a3b8", 1.4))
    for i, count in enumerate(counts):
        bh = count / max_c * plot_h
        bx = plot_x + i * (bar_w + bar_gap)
        by = plot_y + plot_h - bh
        out.append(rect(bx, by, bar_w, bh, "#86efac", rx=3))
        label = f"{bins[i]}-{bins[i + 1]}"
        out.append(svg_text(bx + bar_w / 2, plot_y + plot_h + 18, label, size=9, fill="#475569", anchor="middle"))
    out.append(svg_text(plot_x - 8, plot_y + 4, fmt_int(max_c), size=10, fill="#64748b", anchor="end"))
    out.append(svg_text(plot_x, y + h - 8, f"{fmt_int(len(values))} protein sequences across all samples", size=11, fill="#64748b"))
    return out


def draw_retained_panel(samples: list[dict[str, object]], x: int, y: int, w: int, h: int) -> list[str]:
    out = panel_label(x, y, "D", "Prediction-ready retained fraction")
    plot_x, plot_y = x + 48, y + 32
    plot_w, plot_h = w - 72, h - 70
    out.append(line(plot_x, plot_y + plot_h, plot_x + plot_w, plot_y + plot_h, "#94a3b8", 1.4))
    out.append(line(plot_x, plot_y, plot_x, plot_y + plot_h, "#94a3b8", 1.4))
    bar_gap = 8
    bar_w = (plot_w - bar_gap * (len(samples) - 1)) / len(samples)
    for i, row in enumerate(samples):
        value = float(row["retained_pct"])
        bh = value / 100.0 * plot_h
        bx = plot_x + i * (bar_w + bar_gap)
        by = plot_y + plot_h - bh
        out.append(rect(bx, by, bar_w, bh, "#fb7185", rx=3))
        out.append(svg_text(bx + bar_w / 2, plot_y + plot_h + 18, str(i + 1), size=10, fill="#475569", anchor="middle"))
    for pct in (0, 50, 100):
        gy = plot_y + plot_h - plot_h * pct / 100
        out.append(line(plot_x, gy, plot_x + plot_w, gy, "#e2e8f0", 0.8))
        out.append(svg_text(plot_x - 8, gy + 4, f"{pct}%", size=10, fill="#64748b", anchor="end"))
    return out


def draw_status_panel(samples: list[dict[str, object]], x: int, y: int, w: int, h: int) -> list[str]:
    out = panel_label(x, y, "E", "Sample completion")
    completed = sum(1 for row in samples if row["completed"])
    failed = len(samples) - completed
    total = max(len(samples), 1)
    cx, cy = x + 140, y + 116
    radius = 72
    complete_angle = 2 * math.pi * completed / total
    out.append(f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="none" stroke="#e2e8f0" stroke-width="26"/>')
    dash = complete_angle * radius
    circumference = 2 * math.pi * radius
    out.append(
        f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="none" stroke="#22c55e" '
        f'stroke-width="26" stroke-dasharray="{dash:.1f} {circumference - dash:.1f}" '
        f'transform="rotate(-90 {cx} {cy})"/>'
    )
    out.append(svg_text(cx, cy - 6, f"{completed}/{total}", size=30, weight=800, fill="#0f172a", anchor="middle"))
    out.append(svg_text(cx, cy + 22, "complete", size=13, fill="#64748b", anchor="middle"))
    out.append(svg_text(x + 250, y + 98, f"Completed samples: {completed}", size=15, weight=700, fill="#166534"))
    out.append(svg_text(x + 250, y + 126, f"Incomplete / failed: {failed}", size=15, weight=700, fill="#9f1239"))
    out.append(svg_text(x + 250, y + 162, "Status is read from amyloid_predictors_status.tsv.", size=12, fill="#64748b"))
    return out


def draw_sample_map(samples: list[dict[str, object]], x: int, y: int, w: int, h: int) -> list[str]:
    out = panel_label(x, y, "F", "Sample index")
    row_h = 23
    for i, row in enumerate(samples):
        yy = y + 34 + i * row_h
        out.append(svg_text(x, yy, str(i + 1), size=11, weight=800, fill="#475569"))
        out.append(svg_text(x + 24, yy, str(row["sample"]), size=11, fill="#334155"))
        out.append(svg_text(x + 132, yy, fmt_int(int(row["total_sequences"])), size=11, fill="#334155"))
    out.append(svg_text(x + 132, y + 20, "sequences", size=10, weight=700, fill="#64748b"))
    return out


def write_figure(samples: list[dict[str, object]], lengths: list[dict[str, object]]) -> None:
    width, height = 1500, 1040
    out: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfcfe"/>',
        svg_text(54, 58, "Technical validation of patient-specific expressed proteome reconstruction", size=27, weight=800),
        svg_text(
            54,
            86,
            "Cohort-level metrics derived from Amylice reconstructed proteomes and prediction-ready outputs.",
            size=15,
            fill="#64748b",
        ),
    ]
    out.extend(draw_summary_cards(samples, 54, 112, 1392, 118))
    out.extend(draw_bar_panel(samples, 54, 292, 650, 270))
    out.extend(draw_stacked_panel(samples, 796, 292, 650, 270))
    out.extend(draw_length_panel(lengths, 54, 632, 650, 280))
    out.extend(draw_retained_panel(samples, 796, 632, 410, 280))
    out.extend(draw_status_panel(samples, 54, 930, 520, 92))
    out.extend(draw_sample_map(samples, 1238, 632, 208, 280))
    out.append(svg_text(54, 1018, f"Source: {WORK_ROOT}", size=11, fill="#94a3b8"))
    out.append("</svg>")
    FIGURE_SVG.write_text("\n".join(out), encoding="utf-8")


def main() -> None:
    samples, lengths = collect_metrics()
    if not samples:
        raise SystemExit(f"No sample FASTA files found under {WORK_ROOT}")
    write_metrics(samples)
    write_figure(samples, lengths)
    total_lengths = len(lengths)
    print(f"Wrote {METRICS_CSV}")
    print(f"Wrote {FIGURE_SVG}")
    print(f"Samples: {len(samples)}; protein sequences: {total_lengths:,}")


if __name__ == "__main__":
    main()
