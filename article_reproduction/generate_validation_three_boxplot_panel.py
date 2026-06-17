import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib_codex_cache")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili"
FIG3 = BASE / "figure3_graphical_evidence"
OUT = BASE / "validation" / "q17q83_probability_groups"

BURDEN_PATH = FIG3 / "sample_level_pc_q17q83_amyloidogenic_burden.tsv"
TARGET_GENE_PATH = OUT / "target_gene_consensus_metric_contributions.tsv"
SIGNATURE_PATH = OUT / "validation_signature_scores_selected.tsv"

GROUP_ORDER = ["AC", "TOL"]
PALETTE = {"AC": "#3B6EA8", "TOL": "#B24A3B"}

METRICS = [
    (
        "pc_q17q83_amyloidogenic_log2TPM_burden",
        "absolute_burden",
        "Burden",
        "Amyloidogenic burden\nΣ log2(TPM+1)",
    ),
    (
        "pc_q17q83_amyloidogenic_log2TPM_fraction",
        "normalized_fraction",
        "Fraction",
        "Amyloidogenic fraction\nΣ log2(TPM+1) / all protein-coding",
    ),
    (
        "pc_q17q83_kdmax_weighted_amyloidogenic_log2TPM_fraction",
        "KDmax_weighted_fraction",
        "KDmax-weighted Fraction",
        "KDmax-weighted fraction\nΣ(KDmax×log2(TPM+1)) / all protein-coding",
    ),
]

GENES = ["OR10J1", "SRPK1", "TP53BP1"]
SIGNATURES = [
    ("UPR_score", "UPR"),
    ("Inflammatory_score", "Inflammatory"),
    ("IFNG_score", "Interferon-gamma"),
    ("IFNA_score", "Interferon-alpha"),
    ("Myeloid_score", "Myeloid"),
]


def cliffs_delta(x, y):
    x = np.asarray(pd.Series(x).dropna().astype(float))
    y = np.asarray(pd.Series(y).dropna().astype(float))
    if len(x) == 0 or len(y) == 0:
        return np.nan
    gt = sum(np.sum(xi > y) for xi in x)
    lt = sum(np.sum(xi < y) for xi in x)
    return (gt - lt) / (len(x) * len(y))


def fmt(value):
    if pd.isna(value):
        return "NA"
    if abs(value) >= 100:
        return f"{value:.1f}"
    if abs(value) >= 1:
        return f"{value:.2f}"
    return f"{value:.3g}"


def metric_stats(df, metric):
    ac = df.loc[df["Group"].eq("AC"), metric].dropna().astype(float)
    tol = df.loc[df["Group"].eq("TOL"), metric].dropna().astype(float)
    if len(ac) and len(tol):
        p = mannwhitneyu(tol, ac, alternative="two-sided").pvalue
        delta = tol.mean() - ac.mean()
        cliff = cliffs_delta(tol, ac)
    else:
        p = delta = cliff = np.nan
    return delta, p, cliff


def box_strip(ax, df, y, ylabel, title=None, annotate=False):
    plot_df = df.dropna(subset=[y]).copy()
    sns.boxplot(
        data=plot_df,
        x="Group",
        y=y,
        order=GROUP_ORDER,
        hue="Group",
        hue_order=GROUP_ORDER,
        palette=PALETTE,
        width=0.48,
        fliersize=0,
        linewidth=1.1,
        legend=False,
        ax=ax,
    )
    sns.stripplot(
        data=plot_df,
        x="Group",
        y=y,
        order=GROUP_ORDER,
        hue="Group",
        hue_order=GROUP_ORDER,
        palette=PALETTE,
        jitter=0.12,
        size=4.2,
        edgecolor="white",
        linewidth=0.55,
        alpha=0.95,
        legend=False,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel(ylabel, fontsize=8.2)
    if title:
        ax.set_title(title, fontsize=10.5, pad=7)
    ax.tick_params(axis="x", labelsize=8.5)
    ax.tick_params(axis="y", labelsize=7.5)
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.7)
    ax.set_box_aspect(1)
    sns.despine(ax=ax)
    if annotate:
        delta, p, cliff = metric_stats(plot_df, y)
        ax.text(
            0.02,
            0.98,
            f"Δ TOL-AC = {fmt(delta)}\nMW p = {fmt(p)}\nCliff's δ = {fmt(cliff)}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.4,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#CBD5DF", "alpha": 0.95},
        )


def main():
    for path in [BURDEN_PATH, TARGET_GENE_PATH, SIGNATURE_PATH]:
        if not path.exists():
            raise FileNotFoundError(path)

    burden = pd.read_csv(BURDEN_PATH, sep="\t")
    target_genes = pd.read_csv(TARGET_GENE_PATH, sep="\t")
    signatures = pd.read_csv(SIGNATURE_PATH, sep="\t")

    sns.set_theme(style="white", context="paper")
    fig = plt.figure(figsize=(15.5, 22.0), constrained_layout=False)
    outer = fig.add_gridspec(3, 1, height_ratios=[1, 3, 0.65], hspace=0.28)

    gs_a = outer[0].subgridspec(1, 3, wspace=0.28)
    first_ax = None
    for col, (burden_metric, _, title, ylabel) in enumerate(METRICS):
        ax = fig.add_subplot(gs_a[0, col])
        if first_ax is None:
            first_ax = ax
        box_strip(ax, burden, burden_metric, ylabel, title=title, annotate=True)
    first_ax.text(-0.16, 1.08, "A", transform=first_ax.transAxes, fontsize=18, fontweight="bold", va="top")

    gs_b = outer[1].subgridspec(3, 3, wspace=0.28, hspace=0.34)
    first_ax = None
    for row, gene in enumerate(GENES):
        gene_df = target_genes[target_genes["gene_name"].eq(gene)].copy()
        for col, (_, target_metric, title, ylabel) in enumerate(METRICS):
            ax = fig.add_subplot(gs_b[row, col])
            if first_ax is None:
                first_ax = ax
            subplot_title = f"{gene} | {title}" if row == 0 else gene
            box_strip(ax, gene_df, target_metric, ylabel, title=subplot_title, annotate=False)
            if row < len(GENES) - 1:
                ax.set_xlabel("")
                ax.set_xticklabels([])
    first_ax.text(-0.16, 1.08, "B", transform=first_ax.transAxes, fontsize=18, fontweight="bold", va="top")

    gs_c = outer[2].subgridspec(1, 5, wspace=0.28)
    first_ax = None
    for col, (signature, title) in enumerate(SIGNATURES):
        ax = fig.add_subplot(gs_c[0, col])
        if first_ax is None:
            first_ax = ax
        sig_df = signatures[["sample", "Group", "subject", signature]].rename(columns={signature: "score"})
        box_strip(ax, sig_df, "score", "Signature z-score", title=title, annotate=False)
    first_ax.text(-0.16, 1.08, "C", transform=first_ax.transAxes, fontsize=18, fontweight="bold", va="top")

    fig.subplots_adjust(left=0.065, right=0.985, top=0.985, bottom=0.04)
    fig.savefig(OUT / "validation_three_boxplot_panel.png", dpi=300)
    fig.savefig(OUT / "validation_three_boxplot_panel.pdf")
    print(OUT / "validation_three_boxplot_panel.png")
    print(OUT / "validation_three_boxplot_panel.pdf")


if __name__ == "__main__":
    main()
