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
INPUT = (
    ROOT
    / "analysis_outputs"
    / "amyloid_bc_control_vs_chili"
    / "seven_signature_panel"
    / "amyloidogenic_burden_scores.tsv"
)
OUT = (
    ROOT
    / "analysis_outputs"
    / "amyloid_bc_control_vs_chili"
    / "baseline_expression_weighted_amyloidogenic_burden"
)

GROUP_ORDER = ["AC", "TOL"]
GROUP_LABELS = {"AC": "AC/control", "TOL": "TOL/ChILI"}
PALETTE = {"AC": "#3B6EA8", "TOL": "#B24A3B"}

METRICS = [
    {
        "column": "Amyloidogenic_Burden",
        "label": "Absolute burden",
        "ylabel": "Σ log2(TPM + 1)\namyloidogenic transcripts",
        "description": "absolute amyloidogenic burden",
    },
    {
        "column": "Amyloidogenic_Fraction",
        "label": "Normalized fraction",
        "ylabel": "Amyloidogenic fraction",
        "description": "normalized amyloidogenic fraction",
    },
    {
        "column": "KDmax_weighted_Amyloidogenic_Fraction",
        "label": "KDmax-weighted fraction",
        "ylabel": "KDmax-weighted\namyloidogenic fraction",
        "description": "KDmax-weighted amyloidogenic fraction",
    },
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
    if abs(value) >= 10:
        return f"{value:.2f}"
    if abs(value) >= 1:
        return f"{value:.3f}"
    return f"{value:.3g}"


def summarize_metric(df, metric):
    ac = df.loc[df["Group"].eq("AC"), metric].dropna().astype(float)
    tol = df.loc[df["Group"].eq("TOL"), metric].dropna().astype(float)
    row = {
        "metric": metric,
        "n_AC": len(ac),
        "mean_AC": ac.mean(),
        "median_AC": ac.median(),
        "sd_AC": ac.std(ddof=1),
        "n_TOL": len(tol),
        "mean_TOL": tol.mean(),
        "median_TOL": tol.median(),
        "sd_TOL": tol.std(ddof=1),
        "delta_mean_TOL_minus_AC": tol.mean() - ac.mean(),
        "delta_median_TOL_minus_AC": tol.median() - ac.median(),
        "cliffs_delta_TOL_vs_AC": cliffs_delta(tol, ac),
    }
    if len(ac) and len(tol):
        row["mannwhitney_exact_p_two_sided"] = mannwhitneyu(
            tol,
            ac,
            alternative="two-sided",
            method="exact",
        ).pvalue
    else:
        row["mannwhitney_exact_p_two_sided"] = np.nan
    return row


def draw_panel(df, stats):
    sns.set_theme(style="white", context="paper")
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.4), constrained_layout=True)

    for ax, metric in zip(axes, METRICS):
        col = metric["column"]
        plot_df = df[["sample", "Group", "subject", col]].dropna().copy()
        sns.boxplot(
            data=plot_df,
            x="Group",
            y=col,
            order=GROUP_ORDER,
            hue="Group",
            hue_order=GROUP_ORDER,
            palette=PALETTE,
            width=0.5,
            fliersize=0,
            linewidth=1.1,
            legend=False,
            ax=ax,
        )
        sns.stripplot(
            data=plot_df,
            x="Group",
            y=col,
            order=GROUP_ORDER,
            hue="Group",
            hue_order=GROUP_ORDER,
            palette=PALETTE,
            jitter=0.11,
            size=4.8,
            edgecolor="white",
            linewidth=0.6,
            alpha=0.96,
            legend=False,
            ax=ax,
        )
        stat = stats.loc[stats["metric"].eq(col)].iloc[0]
        ax.text(
            0.03,
            0.97,
            (
                f"Δ={fmt(stat['delta_mean_TOL_minus_AC'])}\n"
                f"p={fmt(stat['mannwhitney_exact_p_two_sided'])}\n"
                f"δ={fmt(stat['cliffs_delta_TOL_vs_AC'])}"
            ),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.4,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "#CBD5E1",
                "alpha": 0.94,
            },
        )
        ax.set_title(metric["label"], fontsize=10.5, pad=8)
        ax.set_xlabel("")
        ax.set_xticks(range(len(GROUP_ORDER)))
        ax.set_xticklabels([GROUP_LABELS[g] for g in GROUP_ORDER], fontsize=8.5)
        ax.set_ylabel(metric["ylabel"], fontsize=8.5)
        ax.tick_params(axis="y", labelsize=7.8)
        ax.grid(axis="y", color="#E5E7EB", linewidth=0.7)
        sns.despine(ax=ax)

    fig.suptitle(
        "Baseline expression-weighted amyloidogenic burden",
        fontsize=12.5,
        fontweight="bold",
    )
    png = OUT / "baseline_expression_weighted_amyloidogenic_burden_boxplots.png"
    pdf = OUT / "baseline_expression_weighted_amyloidogenic_burden_boxplots.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


def write_report(stats, png, pdf):
    lines = [
        "# Baseline expression-weighted amyloidogenic burden",
        "",
        "Boxplots compare baseline AC/control and TOL/ChILI blood RNA-seq samples using three sample-level amyloidogenic burden metrics: absolute amyloidogenic burden, normalized amyloidogenic fraction, and KDmax-weighted amyloidogenic fraction.",
        "",
        f"- Figure PNG: `{png}`",
        f"- Figure PDF: `{pdf}`",
        f"- Statistics table: `{OUT / 'baseline_expression_weighted_amyloidogenic_burden_statistics.tsv'}`",
        "",
        "## Statistical summary",
        "",
    ]
    for metric in METRICS:
        row = stats.loc[stats["metric"].eq(metric["column"])].iloc[0]
        lines.append(
            "- "
            f"{metric['description']}: "
            f"mean AC={fmt(row['mean_AC'])}, mean TOL={fmt(row['mean_TOL'])}, "
            f"delta TOL-AC={fmt(row['delta_mean_TOL_minus_AC'])}, "
            f"Cliff's delta={fmt(row['cliffs_delta_TOL_vs_AC'])}, "
            f"exact Mann-Whitney p={fmt(row['mannwhitney_exact_p_two_sided'])}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "TOL/ChILI samples show consistently higher baseline expression-weighted amyloidogenic burden across the normalized and KDmax-weighted fractions, with the same direction for absolute burden. Because the cohort is small, these results should be treated as effect-size-oriented discovery evidence rather than a definitive inferential test.",
        ]
    )
    report = OUT / "BASELINE_EXPRESSION_WEIGHTED_AMYLOIDOGENIC_BURDEN_ANALYSIS.md"
    report.write_text("\n".join(lines) + "\n")
    return report


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT, sep="\t")
    df = df[df["Group"].isin(GROUP_ORDER)].copy()
    stats = pd.DataFrame([summarize_metric(df, metric["column"]) for metric in METRICS])
    stats.to_csv(
        OUT / "baseline_expression_weighted_amyloidogenic_burden_statistics.tsv",
        sep="\t",
        index=False,
    )
    df.to_csv(
        OUT / "baseline_expression_weighted_amyloidogenic_burden_sample_values.tsv",
        sep="\t",
        index=False,
    )
    png, pdf = draw_panel(df, stats)
    report = write_report(stats, png, pdf)
    print(png)
    print(pdf)
    print(report)


if __name__ == "__main__":
    main()
