from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyBboxPatch


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "validation" / "two_predictor_consensus"
COMBO_PATH = OUT / "two_predictor_prediction_combination_counts.tsv"
COLLAPSED_PATH = OUT / "two_predictor_collapsed_gencode_pc_annotations.tsv"


COLORS = {
    "navy": "#1B3557",
    "blue": "#2C7FB8",
    "teal": "#22A7A0",
    "amber": "#E6A23C",
    "coral": "#E66A5B",
    "lavender": "#6F6BAA",
    "gray": "#596579",
    "light": "#F6F8FB",
    "line": "#D7DFEA",
}


def fmt_int(x):
    return f"{int(x):,}".replace(",", " ")


def fmt_pct(x):
    return f"{100 * x:.1f}%"


def predictor_counts_from_combo(combo):
    amypred = int(combo.loc[combo["AMYPred_Pred_norm"].eq("Amyloid"), "n_protein_rows"].sum())
    amylogram = int(combo.loc[combo["AmyloGramPy_Pred_norm"].eq("Amyloid"), "n_protein_rows"].sum())
    both = int(
        combo.loc[
            combo["AMYPred_Pred_norm"].eq("Amyloid")
            & combo["AmyloGramPy_Pred_norm"].eq("Amyloid"),
            "n_protein_rows",
        ].sum()
    )
    partial_amylogram = int(
        combo.loc[
            combo["AMYPred_Pred_norm"].isna()
            & combo["AmyloGramPy_Pred_norm"].eq("Amyloid"),
            "n_protein_rows",
        ].sum()
    )
    total = int(combo["n_protein_rows"].sum())
    return {
        "amypred": amypred,
        "amylogram": amylogram,
        "both": both,
        "amypred_only": amypred - both,
        "amylogram_only": amylogram - both,
        "partial_amylogram": partial_amylogram,
        "total": total,
    }


def collapsed_counts(collapsed):
    amypred = int(collapsed["AMYPred_Pred_any"].eq("Amyloid").sum())
    amylogram = int(collapsed["AmyloGramPy_Pred_any"].eq("Amyloid").sum())
    both = int(
        (
            collapsed["AMYPred_Pred_any"].eq("Amyloid")
            & collapsed["AmyloGramPy_Pred_any"].eq("Amyloid")
        ).sum()
    )
    total = int(len(collapsed))
    return {
        "amypred": amypred,
        "amylogram": amylogram,
        "both": both,
        "amypred_only": amypred - both,
        "amylogram_only": amylogram - both,
        "total": total,
    }


def add_card(ax, x, y, w, h, title, value, subtitle, color):
    card = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.035",
        facecolor="white",
        edgecolor=COLORS["line"],
        linewidth=1.1,
        transform=ax.transAxes,
        clip_on=False,
    )
    card.set_path_effects([pe.SimplePatchShadow(offset=(1, -1), alpha=0.08), pe.Normal()])
    ax.add_patch(card)
    ax.text(x + 0.04, y + h - 0.06, title, transform=ax.transAxes, fontsize=8.8, color=COLORS["gray"], weight="bold", va="top")
    ax.text(x + 0.04, y + 0.15, value, transform=ax.transAxes, fontsize=18.5, color=color, weight="bold", va="bottom")
    ax.text(x + 0.04, y + 0.045, subtitle, transform=ax.transAxes, fontsize=8.5, color=COLORS["gray"])


def draw_venn(ax, counts, title):
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.9, 1.9)
    ax.set_ylim(-1.25, 1.35)

    left = Circle((-0.56, 0.04), 0.98, facecolor=COLORS["blue"], alpha=0.58, edgecolor=COLORS["navy"], lw=1.8)
    right = Circle((0.56, 0.04), 0.98, facecolor=COLORS["teal"], alpha=0.58, edgecolor="#11645F", lw=1.8)
    ax.add_patch(left)
    ax.add_patch(right)

    ax.text(-1.1, 1.16, "AMYPred-FRL", ha="center", va="center", fontsize=12, color=COLORS["navy"], weight="bold")
    ax.text(1.1, 1.16, "AmyloGram-Py", ha="center", va="center", fontsize=12, color="#11645F", weight="bold")
    ax.text(-1.18, 0.08, fmt_int(counts["amypred_only"]), ha="center", va="center", fontsize=13.5, color="white", weight="bold")
    ax.text(0, 0.08, fmt_int(counts["both"]), ha="center", va="center", fontsize=15, color="white", weight="bold")
    ax.text(1.18, 0.08, fmt_int(counts["amylogram_only"]), ha="center", va="center", fontsize=13.5, color="white", weight="bold")
    ax.text(-1.18, -0.18, "AMYPred only", ha="center", va="center", fontsize=8.3, color="white")
    ax.text(0, -0.18, "shared", ha="center", va="center", fontsize=8.8, color="white")
    ax.text(1.18, -0.18, "AmyloGram only", ha="center", va="center", fontsize=8.3, color="white")
    ax.text(0, -1.03, title, ha="center", va="center", fontsize=8.8, color=COLORS["gray"], linespacing=1.25)


def main():
    combo = pd.read_csv(COMBO_PATH, sep="\t")
    collapsed = pd.read_csv(COLLAPSED_PATH, sep="\t")
    counts = predictor_counts_from_combo(combo)
    collapsed_set = collapsed_counts(collapsed)

    consensus = (
        combo.groupby("Consensus_2pred", dropna=False)["n_protein_rows"]
        .sum()
        .reindex(["Amyloid", "Discordant", "Partial", "Non-Amyloid"])
        .fillna(0)
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "figure.dpi": 160,
        }
    )

    fig = plt.figure(figsize=(13.5, 10.2), facecolor=COLORS["light"])
    gs = fig.add_gridspec(3, 4, height_ratios=[1.45, 2.15, 2.2], hspace=0.5, wspace=0.55)

    ax_title = fig.add_subplot(gs[0, :])
    ax_title.axis("off")
    ax_title.text(
        0.01,
        0.92,
        "Two-predictor amyloidogenicity map",
        transform=ax_title.transAxes,
        fontsize=24,
        weight="bold",
        color=COLORS["navy"],
    )
    ax_title.text(
        0.01,
        0.74,
        "AMYPred-FRL and AmyloGram-Py predictions across validation protein rows",
        transform=ax_title.transAxes,
        fontsize=12.5,
        color=COLORS["gray"],
    )

    add_card(
        ax_title,
        0.01,
        0.16,
        0.21,
        0.42,
        "Total rows",
        fmt_int(counts["total"]),
        "raw protein-level records",
        COLORS["navy"],
    )
    add_card(
        ax_title,
        0.255,
        0.16,
        0.21,
        0.42,
        "AMYPred-FRL amyloid",
        fmt_int(counts["amypred"]),
        fmt_pct(counts["amypred"] / counts["total"]),
        COLORS["blue"],
    )
    add_card(
        ax_title,
        0.50,
        0.16,
        0.21,
        0.42,
        "AmyloGram-Py amyloid",
        fmt_int(counts["amylogram"]),
        fmt_pct(counts["amylogram"] / counts["total"]),
        COLORS["teal"],
    )
    add_card(
        ax_title,
        0.745,
        0.16,
        0.21,
        0.42,
        "Shared amyloid calls",
        fmt_int(counts["both"]),
        fmt_pct(counts["both"] / counts["total"]),
        COLORS["lavender"],
    )

    ax_bar = fig.add_subplot(gs[1, :2])
    labels = ["AMYPred-FRL", "AmyloGram-Py", "Both amyloid"]
    values = [counts["amypred"], counts["amylogram"], counts["both"]]
    colors = [COLORS["blue"], COLORS["teal"], COLORS["lavender"]]
    bars = ax_bar.bar(labels, values, color=colors, width=0.64, edgecolor="white", linewidth=1.5)
    ax_bar.set_title("Amyloid predictions by tool", loc="left", color=COLORS["navy"], weight="bold", pad=12)
    ax_bar.set_ylabel("Number of protein rows")
    ax_bar.grid(axis="y", color=COLORS["line"], linewidth=0.9)
    ax_bar.set_axisbelow(True)
    ax_bar.spines[["top", "right", "left"]].set_visible(False)
    ax_bar.tick_params(axis="y", length=0)
    ax_bar.set_ylim(0, max(values) * 1.22)
    for bar, value in zip(bars, values):
        ax_bar.text(
            bar.get_x() + bar.get_width() / 2,
            value + max(values) * 0.035,
            fmt_int(value),
            ha="center",
            va="bottom",
            fontsize=12,
            weight="bold",
            color=COLORS["navy"],
        )

    ax_venn = fig.add_subplot(gs[1, 2:])
    draw_venn(
        ax_venn,
        counts,
        f"Set overlap among amyloid calls\n{fmt_int(counts['partial_amylogram'])} AmyloGram-Py amyloid rows had missing AMYPred-FRL calls.",
    )
    ax_venn.set_title("Overlap of amyloid predictions", loc="left", color=COLORS["navy"], weight="bold", pad=10)

    ax_consensus = fig.add_subplot(gs[2, :2])
    consensus_colors = [COLORS["lavender"], COLORS["coral"], COLORS["amber"], COLORS["gray"]]
    y = np.arange(len(consensus))
    ax_consensus.barh(y, consensus.values, color=consensus_colors, edgecolor="white", linewidth=1.4)
    ax_consensus.set_yticks(y, consensus.index)
    ax_consensus.invert_yaxis()
    ax_consensus.set_xlabel("Number of protein rows")
    ax_consensus.set_title("Consensus classes", loc="left", color=COLORS["navy"], weight="bold", pad=12)
    ax_consensus.grid(axis="x", color=COLORS["line"], linewidth=0.9)
    ax_consensus.set_axisbelow(True)
    ax_consensus.spines[["top", "right", "left"]].set_visible(False)
    ax_consensus.tick_params(axis="y", length=0)
    for yi, value in enumerate(consensus.values):
        ax_consensus.text(
            value + consensus.max() * 0.02,
            yi,
            f"{fmt_int(value)} ({fmt_pct(value / counts['total'])})",
            va="center",
            fontsize=10.5,
            color=COLORS["navy"],
            weight="bold",
        )
    ax_consensus.set_xlim(0, consensus.max() * 1.3)

    ax_note = fig.add_subplot(gs[2, 2:])
    ax_note.axis("off")
    note_box = FancyBboxPatch(
        (0.02, 0.08),
        0.96,
        0.82,
        boxstyle="round,pad=0.025,rounding_size=0.04",
        facecolor="white",
        edgecolor=COLORS["line"],
        linewidth=1.1,
        transform=ax_note.transAxes,
    )
    note_box.set_path_effects([pe.SimplePatchShadow(offset=(1, -1), alpha=0.08), pe.Normal()])
    ax_note.add_patch(note_box)
    ax_note.text(
        0.08,
        0.77,
        "Collapsed protein-coding transcript check",
        transform=ax_note.transAxes,
        fontsize=14,
        color=COLORS["navy"],
        weight="bold",
    )
    collapsed_lines = [
        ("Total collapsed transcripts", collapsed_set["total"]),
        ("AMYPred-FRL amyloid", collapsed_set["amypred"]),
        ("AmyloGram-Py amyloid", collapsed_set["amylogram"]),
        ("Shared amyloid calls", collapsed_set["both"]),
    ]
    for i, (label, value) in enumerate(collapsed_lines):
        y0 = 0.64 - i * 0.14
        ax_note.text(0.09, y0, label, transform=ax_note.transAxes, fontsize=10.8, color=COLORS["gray"])
        ax_note.text(0.88, y0, fmt_int(value), transform=ax_note.transAxes, fontsize=11.5, color=COLORS["navy"], weight="bold", ha="right")

    ax_note.text(
        0.08,
        0.09,
        "Interpretation: AmyloGram-Py gives broader amyloid coverage, while the shared set defines the strict two-predictor amyloid-positive core.",
        transform=ax_note.transAxes,
        fontsize=9.5,
        color=COLORS["gray"],
        wrap=True,
    )

    fig.text(
        0.012,
        0.012,
        "Data source: two_predictor_consensus validation tables. Consensus: Amyloid = both tools amyloid; Discordant = both present but disagree; Partial = one tool missing.",
        fontsize=8.5,
        color=COLORS["gray"],
    )

    png = OUT / "two_predictor_prediction_overlap_panel.png"
    pdf = OUT / "two_predictor_prediction_overlap_panel.pdf"
    summary = OUT / "two_predictor_prediction_overlap_panel_summary.tsv"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(
        [
            {"level": "raw_protein_rows", **counts},
            {"level": "collapsed_gencode_protein_coding_transcripts", **collapsed_set, "partial_amylogram": np.nan},
        ]
    ).to_csv(summary, sep="\t", index=False)
    print(png)
    print(pdf)
    print(summary)


if __name__ == "__main__":
    main()
