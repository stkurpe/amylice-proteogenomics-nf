import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib_codex_cache")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import zarr
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parent
ZARR_PATH = ROOT / "project.zarr"
META_PATH = ROOT / "GSE287540_SraRunTable.csv"
BASELINE_BURDEN_PATH = (
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
    / "marker_gene_amyloidogenic_burden"
)

TARGET_GENES = ["ALB", "AFP", "CRP", "IL6", "SLFN11"]
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
GROUP_ORDER = ["AC", "TOL"]
GROUP_LABELS = {"AC": "AC/control", "TOL": "TOL/ChILI"}
PALETTE = {"AC": "#3B6EA8", "TOL": "#B24A3B"}

METRICS = [
    (
        "absolute_amyloidogenic_burden",
        "Absolute burden",
        "Σ log2(TPM + 1)\namyloidogenic transcripts",
    ),
    (
        "normalized_amyloidogenic_fraction",
        "Normalized fraction",
        "Amyloidogenic fraction\nwithin gene",
    ),
    (
        "KDmax_weighted_amyloidogenic_fraction",
        "KDmax-weighted fraction",
        "KDmax-weighted\namyloidogenic fraction",
    ),
]


def clean_enst(value):
    if pd.isna(value):
        return np.nan
    match = re.search(r"ENST\d+(?:\.\d+)?", str(value))
    return match.group(0).split(".")[0] if match else np.nan


def zarr_group_to_df(group):
    return pd.DataFrame({col: group[col][:] for col in group.array_keys()})


def collapse_pc_class(values):
    vals = set(pd.Series(values).dropna().astype(str))
    if "Amyloidogenic" in vals:
        return "Amyloidogenic"
    if "Intermediate" in vals:
        return "Intermediate"
    if "Non-Amyloidogenic" in vals:
        return "Non-Amyloidogenic"
    return np.nan


def load_baseline_metadata(root):
    meta = pd.read_csv(META_PATH)
    meta["raw_Group"] = meta["Group"].astype(str)
    meta["Group"] = meta["raw_Group"].map(GROUP_MAP)
    meta["Run"] = meta["Run"].astype(str)
    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    meta = meta[
        meta["Run"].isin(zarr_samples)
        & meta["time"].eq("BC")
        & meta["Group"].isin(GROUP_ORDER)
    ].copy()
    selected = [sample for sample in zarr_samples if sample in set(meta["Run"])]
    return meta.rename(columns={"Run": "sample"}).set_index("sample").loc[selected].reset_index(), selected


def gene_sample_metrics(root, meta, samples):
    rows = []
    meta_small = meta[["sample", "Group", "subject"]].copy()

    for sample in samples:
        expr = zarr_group_to_df(root[f"samples/{sample}/expression"])[
            ["transcript_id", "tpm"]
        ].copy()
        expr["Transcript_ID_clean"] = expr["transcript_id"].map(clean_enst)
        expr["log2_TPM_plus1"] = np.log2(
            pd.to_numeric(expr["tpm"], errors="coerce").fillna(0) + 1
        )
        expr = (
            expr[["Transcript_ID_clean", "log2_TPM_plus1"]]
            .dropna()
            .drop_duplicates("Transcript_ID_clean")
        )

        amy = zarr_group_to_df(root[f"samples/{sample}/amyloid"]).copy()
        protein_features = zarr_group_to_df(root[f"samples/{sample}/protein_features"])[
            ["protein_id", "KD_max"]
        ].copy()
        amy = amy.merge(protein_features, left_on="Sequence_ID", right_on="protein_id", how="left")
        amy["Transcript_ID_clean"] = amy["Transcript_ID_clean"].map(clean_enst)
        amy["gene_name"] = amy["gencode_gene_name"].astype(str)
        amy["KD_max"] = pd.to_numeric(amy["KD_max"], errors="coerce")
        amy = amy[
            amy["gene_name"].isin(TARGET_GENES)
            & amy["ProteinCoding_BothPredictors_Q17_Q83_Class"].notna()
        ].copy()

        collapsed = (
            amy.dropna(subset=["Transcript_ID_clean"])
            .groupby(["gene_name", "Transcript_ID_clean"], as_index=False)
            .agg(
                ProteinCoding_BothPredictors_Q17_Q83_Class=(
                    "ProteinCoding_BothPredictors_Q17_Q83_Class",
                    collapse_pc_class,
                ),
                KD_max=("KD_max", "max"),
            )
        )
        tx = collapsed.merge(expr, on="Transcript_ID_clean", how="left")
        tx["log2_TPM_plus1"] = tx["log2_TPM_plus1"].fillna(0)
        tx["is_amyloidogenic"] = tx[
            "ProteinCoding_BothPredictors_Q17_Q83_Class"
        ].eq("Amyloidogenic")
        tx["kd_weighted_log2_TPM"] = tx["KD_max"] * tx["log2_TPM_plus1"]

        for gene in TARGET_GENES:
            g = tx[tx["gene_name"].eq(gene)].copy()
            all_log = g["log2_TPM_plus1"].sum()
            amy_log = g.loc[g["is_amyloidogenic"], "log2_TPM_plus1"].sum()
            kd_all = g["kd_weighted_log2_TPM"].sum(skipna=True)
            kd_amy = g.loc[g["is_amyloidogenic"], "kd_weighted_log2_TPM"].sum(skipna=True)
            rows.append(
                {
                    "sample": sample,
                    "gene_name": gene,
                    "n_classified_transcripts": int(g.shape[0]),
                    "n_amyloidogenic_transcripts": int(g["is_amyloidogenic"].sum()),
                    "gene_total_log2_TPM_plus1": all_log,
                    "absolute_amyloidogenic_burden": amy_log,
                    "normalized_amyloidogenic_fraction": amy_log / all_log if all_log else np.nan,
                    "KDmax_weighted_total": kd_all,
                    "KDmax_weighted_amyloidogenic": kd_amy,
                    "KDmax_weighted_amyloidogenic_fraction": kd_amy / kd_all if kd_all else np.nan,
                }
            )

    return pd.DataFrame(rows).merge(meta_small, on="sample", how="left")


def cliffs_delta(x, y):
    x = np.asarray(pd.Series(x).dropna().astype(float))
    y = np.asarray(pd.Series(y).dropna().astype(float))
    if len(x) == 0 or len(y) == 0:
        return np.nan
    gt = sum(np.sum(xi > y) for xi in x)
    lt = sum(np.sum(xi < y) for xi in x)
    return (gt - lt) / (len(x) * len(y))


def summarize(df):
    rows = []
    for gene in TARGET_GENES:
        gene_df = df[df["gene_name"].eq(gene)]
        for metric, _, _ in METRICS:
            ac = gene_df.loc[gene_df["Group"].eq("AC"), metric].dropna().astype(float)
            tol = gene_df.loc[gene_df["Group"].eq("TOL"), metric].dropna().astype(float)
            row = {
                "gene_name": gene,
                "metric": metric,
                "n_AC": len(ac),
                "mean_AC": ac.mean(),
                "median_AC": ac.median(),
                "n_TOL": len(tol),
                "mean_TOL": tol.mean(),
                "median_TOL": tol.median(),
                "delta_mean_TOL_minus_AC": tol.mean() - ac.mean(),
                "cliffs_delta_TOL_vs_AC": cliffs_delta(tol, ac),
            }
            if len(ac) and len(tol) and (ac.nunique() > 1 or tol.nunique() > 1):
                row["mannwhitney_exact_p_two_sided"] = mannwhitneyu(
                    tol,
                    ac,
                    alternative="two-sided",
                    method="exact",
                ).pvalue
            else:
                row["mannwhitney_exact_p_two_sided"] = np.nan
            rows.append(row)
    return pd.DataFrame(rows)


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


def draw_metric_panels(df, stats):
    sns.set_theme(style="white", context="paper")
    out_paths = []
    for metric, title, ylabel in METRICS:
        fig, ax = plt.subplots(figsize=(9.8, 4.2), constrained_layout=True)
        plot_df = df.dropna(subset=[metric]).copy()
        sns.boxplot(
            data=plot_df,
            x="gene_name",
            y=metric,
            hue="Group",
            order=TARGET_GENES,
            hue_order=GROUP_ORDER,
            palette=PALETTE,
            width=0.68,
            fliersize=0,
            linewidth=1.1,
            ax=ax,
        )
        sns.stripplot(
            data=plot_df,
            x="gene_name",
            y=metric,
            hue="Group",
            order=TARGET_GENES,
            hue_order=GROUP_ORDER,
            palette=PALETTE,
            dodge=True,
            jitter=0.09,
            size=4.2,
            edgecolor="white",
            linewidth=0.55,
            alpha=0.95,
            ax=ax,
        )
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            handles[:2],
            [GROUP_LABELS[g] for g in GROUP_ORDER],
            frameon=False,
            loc="upper right",
            title="",
        )
        ax.set_title(f"{title} by marker gene", fontsize=12, fontweight="bold", pad=9)
        ax.set_xlabel("")
        ax.set_ylabel(ylabel, fontsize=9.2)
        ax.tick_params(axis="x", labelsize=9)
        ax.tick_params(axis="y", labelsize=8)
        ax.grid(axis="y", color="#E5E7EB", linewidth=0.7)
        sns.despine(ax=ax)

        stem = metric.replace("_", "-")
        png = OUT / f"marker_gene_{stem}_boxplots.png"
        pdf = OUT / f"marker_gene_{stem}_boxplots.pdf"
        fig.savefig(png, dpi=300)
        fig.savefig(pdf)
        plt.close(fig)
        out_paths.extend([png, pdf])
    return out_paths


def write_report(stats, figures):
    lines = [
        "# Marker-gene amyloidogenic burden",
        "",
        "This analysis recalculates amyloidogenic burden within ALB, AFP, CRP, IL6, and SLFN11 for baseline AC/control and TOL/ChILI samples. Metrics are gene-local: the normalized and KDmax-weighted fractions use each gene's classified transcripts as the denominator, not the whole protein-coding universe.",
        "",
        "## Outputs",
        "",
    ]
    for path in figures:
        lines.append(f"- `{path}`")
    lines.extend(
        [
            f"- `{OUT / 'marker_gene_amyloidogenic_burden_sample_values.tsv'}`",
            f"- `{OUT / 'marker_gene_amyloidogenic_burden_statistics.tsv'}`",
            "",
            "## Statistical summary",
            "",
        ]
    )
    for gene in TARGET_GENES:
        lines.append(f"### {gene}")
        for metric, label, _ in METRICS:
            row = stats[(stats["gene_name"].eq(gene)) & (stats["metric"].eq(metric))].iloc[0]
            lines.append(
                "- "
                f"{label}: mean AC={fmt(row['mean_AC'])}, mean TOL={fmt(row['mean_TOL'])}, "
                f"delta TOL-AC={fmt(row['delta_mean_TOL_minus_AC'])}, "
                f"Cliff's delta={fmt(row['cliffs_delta_TOL_vs_AC'])}, "
                f"exact Mann-Whitney p={fmt(row['mannwhitney_exact_p_two_sided'])}."
            )
        lines.append("")
    lines.extend(
        [
            "## Interpretation",
            "",
            "CRP and IL6 have no measurable baseline RNA-seq amyloidogenic burden in these samples, so their gene-local metrics are not informative here. ALB shows higher absolute amyloidogenic burden in TOL/ChILI, driven by the few samples where ALB has classified amyloidogenic expression. ALB within-gene fractions are high among samples with nonzero denominators, but the denominator is available in only two AC and two TOL samples. AFP has a single AC sample with amyloidogenic burden and no TOL signal in this matched baseline set. SLFN11 is expressed in several samples, but its classified transcripts are non-amyloidogenic in this q17/q83 consensus framework, so SLFN11 gene-local amyloidogenic burden is zero. This differs from the previous analysis where SLFN11 expression was evaluated as a marker associated with global protein-coding amyloidogenic burden.",
        ]
    )
    report = OUT / "MARKER_GENE_AMYLOIDOGENIC_BURDEN_ANALYSIS.md"
    report.write_text("\n".join(lines) + "\n")
    return report


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    meta, samples = load_baseline_metadata(root)
    if BASELINE_BURDEN_PATH.exists():
        baseline_samples = pd.read_csv(BASELINE_BURDEN_PATH, sep="\t")["sample"].astype(str).tolist()
        samples = [sample for sample in baseline_samples if sample in set(samples)]
        meta = meta.set_index("sample").loc[samples].reset_index()
    df = gene_sample_metrics(root, meta, samples)
    stats = summarize(df)
    df.to_csv(OUT / "marker_gene_amyloidogenic_burden_sample_values.tsv", sep="\t", index=False)
    stats.to_csv(OUT / "marker_gene_amyloidogenic_burden_statistics.tsv", sep="\t", index=False)
    figures = draw_metric_panels(df, stats)
    report = write_report(stats, figures)
    for path in figures:
        print(path)
    print(report)


if __name__ == "__main__":
    main()
