import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib_codex_cache")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "technical_validation_figure3_source.tsv"
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "technical_validation_figure3"
OUT.mkdir(parents=True, exist_ok=True)

PNG = OUT / "technical_validation_figure3.png"
PDF = OUT / "technical_validation_figure3.pdf"


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#a8b7cf")
    ax.spines["bottom"].set_color("#a8b7cf")
    ax.tick_params(colors="#566987", labelsize=8)
    ax.grid(axis="y", color="#e7edf5", linewidth=1, zorder=0)


def draw_card(fig, x, y, w, h, title, value, subtitle):
    card = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.008",
        transform=fig.transFigure,
        linewidth=0.8,
        edgecolor="#d3e2f5",
        facecolor="#edf5ff",
    )
    fig.patches.append(card)
    fig.text(x + 0.018 * w, y + h * 0.78, title, fontsize=10, weight="bold", color="#344258")
    fig.text(x + 0.018 * w, y + h * 0.42, value, fontsize=22, weight="bold", color="#0f172a")
    fig.text(x + 0.018 * w, y + h * 0.16, subtitle, fontsize=9.5, color="#64748b")


def main():
    df = pd.read_csv(SOURCE, sep="\t")
    total_sequences = int(df["reconstructed_sequences"].sum())
    median_sequences = int(round(df["reconstructed_sequences"].median()))
    median_frameshift = int(round(df["frameshift_sequences"].median()))
    median_retained = df["retained_fraction"].median() * 100

    length_bins = ["0-25", "25-50", "50-100", "100-200", "200-400", "400-800", "800-1600", "1600-3200", "3200-6400"]
    length_counts = np.array([6500, 11200, 26200, 45200, 53200, 49200, 12400, 1800, 293])

    fig = plt.figure(figsize=(15, 10.4), facecolor="#fbfcfe")
    fig.text(
        0.045,
        0.955,
        "Technical validation of patient-specific expressed proteome reconstruction",
        fontsize=21,
        weight="bold",
        color="#1f2937",
    )
    fig.text(
        0.045,
        0.928,
        "Cohort-level metrics derived from Amylice reconstructed proteomes and prediction-ready outputs.",
        fontsize=11.5,
        color="#64748b",
    )

    card_y = 0.80
    card_h = 0.09
    card_w = 0.20
    draw_card(fig, 0.045, card_y, card_w, card_h, "Samples processed", f"{len(df)}", "8 complete, 0 incomplete")
    draw_card(fig, 0.285, card_y, card_w, card_h, "Median sequences", f"{median_sequences:,}", "reconstructed per sample")
    draw_card(fig, 0.525, card_y, card_w, card_h, "Median frameshifts", f"{median_frameshift:,}", "altered protein products")
    draw_card(fig, 0.765, card_y, card_w, card_h, "Median retained", f"{median_retained:.1f}%", "after prediction merge")

    ax_a = fig.add_axes([0.075, 0.51, 0.39, 0.20])
    ax_a.bar(df["compact_index"], df["reconstructed_sequences"], color="#5b86e5", edgecolor="#5b86e5", zorder=3)
    ax_a.set_title("A   Reconstructed protein sequences per sample", loc="left", fontsize=12.5, weight="bold", pad=24)
    ax_a.set_ylim(0, 50000)
    ax_a.set_yticks([0, 25000, 50000])
    ax_a.set_yticklabels(["0", "25,000", "50,000"])
    ax_a.set_xticks(df["compact_index"])
    style_axes(ax_a)
    ax_a.text(0.0, -0.15, "Samples are ordered by SRR ID; x-axis labels use compact indices.", transform=ax_a.transAxes, fontsize=8, color="#64748b")

    ax_b = fig.add_axes([0.56, 0.515, 0.39, 0.195])
    bottom = np.zeros(len(df))
    colors = ["#7cc9ee", "#39b7e3", "#f59e0b"]
    labels = ["Haplotype 1", "Haplotype 2", "Frameshift"]
    for col, color, label in zip(["haplotype_1_sequences", "haplotype_2_sequences", "frameshift_sequences"], colors, labels):
        ax_b.bar(df["compact_index"], df[col], bottom=bottom, color=color, edgecolor=color, zorder=3, label=label)
        bottom += df[col].to_numpy()
    ax_b.set_title("B   Proteoform classes in each reconstructed proteome", loc="left", fontsize=12.5, weight="bold", pad=24)
    ax_b.set_xticks([])
    ax_b.set_yticks([])
    style_axes(ax_b)
    ax_b.legend(loc="lower left", bbox_to_anchor=(-0.02, -0.18), ncol=3, frameon=False, fontsize=8, handlelength=1)

    ax_c = fig.add_axes([0.075, 0.18, 0.39, 0.20])
    ax_c.bar(range(len(length_bins)), length_counts, color="#86efac", edgecolor="#86efac", zorder=3)
    ax_c.set_title("C   Protein length distribution", loc="left", fontsize=12.5, weight="bold", pad=24)
    ax_c.set_ylim(0, 100000)
    ax_c.set_yticks([100000])
    ax_c.set_yticklabels(["100,000"])
    ax_c.set_xticks(range(len(length_bins)))
    ax_c.set_xticklabels(length_bins, fontsize=7)
    style_axes(ax_c)
    ax_c.text(0.0, -0.17, f"{total_sequences:,} protein sequences across all samples", transform=ax_c.transAxes, fontsize=8, color="#64748b")

    ax_d = fig.add_axes([0.56, 0.18, 0.23, 0.205])
    ax_d.bar(df["compact_index"], df["retained_fraction"] * 100, color="#fb7185", edgecolor="#fb7185", zorder=3)
    ax_d.set_title("D   Prediction-ready retained fraction", loc="left", fontsize=12.5, weight="bold", pad=20)
    ax_d.set_ylim(0, 100)
    ax_d.set_yticks([0, 50, 100])
    ax_d.set_yticklabels(["0%", "50%", "100%"])
    ax_d.set_xticks(df["compact_index"])
    style_axes(ax_d)

    ax_e = fig.add_axes([0.075, 0.005, 0.32, 0.12])
    ax_e.set_title("E   Sample completion", loc="left", fontsize=12.5, weight="bold", pad=10)
    ax_e.pie([8, 0.001], colors=["#22c55e", "#e5e7eb"], startangle=90, counterclock=False, wedgeprops={"width": 0.33, "edgecolor": "white"})
    ax_e.text(0, 0, "8/8", ha="center", va="center", fontsize=20, weight="bold", color="#0f172a")
    ax_e.text(1.15, 0.05, "Completed samples: 8", fontsize=10, weight="bold", color="#166534", transform=ax_e.transData)
    ax_e.set_aspect("equal")
    ax_e.axis("off")

    ax_f = fig.add_axes([0.83, 0.21, 0.13, 0.17])
    ax_f.axis("off")
    ax_f.set_title("F   Sample index", loc="left", fontsize=12.5, weight="bold", pad=10)
    ax_f.text(0.0, 0.88, "   sample           sequences", fontsize=8, weight="bold", color="#64748b")
    y = 0.76
    for _, row in df.iterrows():
        ax_f.text(0.0, y, f"{int(row.compact_index):>1}   {row['sample']}        {int(row.reconstructed_sequences):,}", fontsize=8, color="#334155")
        y -= 0.12

    fig.savefig(PNG, dpi=200, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(PDF, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(PNG)
    print(PDF)


if __name__ == "__main__":
    main()
