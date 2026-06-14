import json
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib_codex_cache")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf
import zarr
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parent
ZARR_PATH = ROOT / "project.zarr"
BASE = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili"
META_PATH = BASE / "sample_metadata_BC.tsv"
OUT = BASE / "proteoform_log2TPM_metric_differences_AC_vs_TOL"
OUT.mkdir(parents=True, exist_ok=True)

CLASS_COL = "ProteinCoding_BothPredictors_Q17_Q83_Class"
AMY_CLASS = "Amyloidogenic"
GROUP_ORDER = ["TOL", "AC"]
PALETTE = {"AC": "#2B6CB0", "TOL": "#C2410C"}
DIRECTION_PALETTE = {"higher_in_AC": "#2B6CB0", "higher_in_TOL": "#C2410C"}

METRICS = [
    "burden_contribution_log2TPM",
    "fraction_contribution_log2TPM",
    "kdmax_weighted_fraction_contribution_log2TPM",
]
METRIC_LABELS = {
    "burden_contribution_log2TPM": "Burden",
    "fraction_contribution_log2TPM": "Fraction",
    "kdmax_weighted_fraction_contribution_log2TPM": "KDmax-weighted Fraction",
}


def zarr_df(group):
    return pd.DataFrame({col: group[col][:] for col in group.attrs["columns"]})


def clean_transcript(value):
    match = re.search(r"(ENST\d+)(?:\.\d+)?", str(value))
    return match.group(1) if match else pd.NA


def as_bool(value):
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def collapse_first(values):
    clean = pd.Series(values, dtype="string").replace("", pd.NA).dropna()
    return clean.iloc[0] if len(clean) else pd.NA


def load_long_table():
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    frames = []
    sample_qc_rows = []

    for sample in root["samples"].group_keys():
        amy = zarr_df(root[f"samples/{sample}/amyloid"])
        expr = zarr_df(root[f"samples/{sample}/expression"])
        features = zarr_df(root[f"samples/{sample}/protein_features"])

        amy[CLASS_COL] = amy[CLASS_COL].astype("string").replace("", pd.NA)
        amy["Transcript_ID_clean"] = amy["Transcript_ID_clean"].astype("string").replace("", pd.NA)
        amy["gencode_gene_id"] = amy["gencode_gene_id"].astype("string").replace("", pd.NA)
        amy["gencode_gene_name"] = amy["gencode_gene_name"].astype("string").replace("", pd.NA)
        amy["is_gencode_protein_coding_bool"] = amy["is_gencode_protein_coding"].map(as_bool)

        features["Sequence_ID"] = features["protein_id"].astype(str)
        feature_cols = [
            "sample",
            "Sequence_ID",
            "protein_length",
            "KD_max",
            "charge",
            "charge_density",
            "beta_propensity",
            "chameleon_score",
            "aromaticity",
            "pI",
        ]
        features = features[[c for c in feature_cols if c in features.columns]].copy()

        amy_features = amy.merge(features, on=["sample", "Sequence_ID"], how="left")
        pc = amy_features[amy_features["is_gencode_protein_coding_bool"]].copy()

        collapsed = (
            pc.groupby(["sample", "Transcript_ID_clean"], dropna=False)
            .agg(
                gencode_gene_id=("gencode_gene_id", collapse_first),
                gencode_gene_name=("gencode_gene_name", collapse_first),
                protein_coding_class=(CLASS_COL, collapse_first),
                KDmax=("KD_max", "max"),
                protein_length=("protein_length", "max"),
                n_protein_rows=("Sequence_ID", "nunique"),
            )
            .reset_index()
        )

        expr["Transcript_ID_clean"] = expr["transcript_id"].map(clean_transcript).astype("string")
        expr = expr[["sample", "Transcript_ID_clean", "transcript_id", "gene_id", "gene_name", "tpm"]].copy()
        d = collapsed.merge(expr, on=["sample", "Transcript_ID_clean"], how="left")
        d["tpm"] = pd.to_numeric(d["tpm"], errors="coerce").fillna(0.0)
        d["KDmax"] = pd.to_numeric(d["KDmax"], errors="coerce")
        d["log2_TPM_plus_1"] = np.log2(d["tpm"].clip(lower=0) + 1.0)
        d["kdmax_x_log2_TPM_plus_1"] = d["KDmax"].fillna(0.0) * d["log2_TPM_plus_1"]
        d["is_amyloidogenic"] = d["protein_coding_class"].eq(AMY_CLASS).fillna(False).astype(bool)

        total_log = d["log2_TPM_plus_1"].sum()
        total_weighted_log = d["kdmax_x_log2_TPM_plus_1"].sum()
        d["burden_contribution_log2TPM"] = np.where(d["is_amyloidogenic"], d["log2_TPM_plus_1"], 0.0)
        d["fraction_contribution_log2TPM"] = np.where(
            d["is_amyloidogenic"] & (total_log > 0),
            d["log2_TPM_plus_1"] / total_log,
            0.0,
        )
        d["kdmax_weighted_fraction_contribution_log2TPM"] = np.where(
            d["is_amyloidogenic"] & (total_weighted_log > 0),
            d["kdmax_x_log2_TPM_plus_1"] / total_weighted_log,
            0.0,
        )

        frames.append(d)
        sample_qc_rows.append(
            {
                "sample": sample,
                "n_protein_coding_rows": int(len(pc)),
                "n_collapsed_protein_coding_transcripts": int(collapsed["Transcript_ID_clean"].nunique()),
                "n_amyloidogenic_transcripts": int(d["is_amyloidogenic"].sum()),
                "all_protein_coding_log2TPM_sum": float(total_log),
                "all_protein_coding_KDmax_log2TPM_sum": float(total_weighted_log),
                "n_missing_expression_after_merge": int(d["transcript_id"].isna().sum()),
                "n_missing_KDmax": int(d["KDmax"].isna().sum()),
            }
        )

    return root, pd.concat(frames, ignore_index=True), pd.DataFrame(sample_qc_rows)


def attach_metadata(long_df):
    meta = pd.read_csv(META_PATH, sep="\t")
    keep_cols = ["sample", "Group", "raw_Group", "subject", "time", "batch", "treatment", "sex", "AGE"]
    keep_cols = [c for c in keep_cols if c in meta.columns]
    out = long_df.merge(meta[keep_cols], on="sample", how="left")
    out = out[out["time"].eq("BC") & out["Group"].isin(["AC", "TOL"])].copy()
    out["Group"] = pd.Categorical(out["Group"], categories=GROUP_ORDER, ordered=True)
    return out


def build_sample_metrics(long_df):
    rows = []
    for sample, g in long_df.groupby("sample", observed=True):
        amy = g[g["is_amyloidogenic"]]
        log_total = g["log2_TPM_plus_1"].sum()
        weighted_total = g["kdmax_x_log2_TPM_plus_1"].sum()
        rows.append(
            {
                "sample": sample,
                "Group": g["Group"].iloc[0],
                "subject": g["subject"].iloc[0] if "subject" in g.columns else pd.NA,
                "time": g["time"].iloc[0] if "time" in g.columns else pd.NA,
                "amyloidogenic_burden_log2TPM": amy["log2_TPM_plus_1"].sum(),
                "amyloidogenic_fraction_log2TPM": amy["log2_TPM_plus_1"].sum() / log_total if log_total else np.nan,
                "kdmax_weighted_amyloidogenic_fraction_log2TPM": amy["kdmax_x_log2_TPM_plus_1"].sum() / weighted_total
                if weighted_total
                else np.nan,
                "all_protein_coding_log2TPM_sum": log_total,
                "all_protein_coding_KDmax_log2TPM_sum": weighted_total,
                "n_amyloidogenic_transcripts": int(amy["Transcript_ID_clean"].nunique()),
                "n_protein_coding_transcripts": int(g["Transcript_ID_clean"].nunique()),
            }
        )
    return pd.DataFrame(rows)


def cohen_d(ac, tol):
    ac = np.asarray(ac, dtype=float)
    tol = np.asarray(tol, dtype=float)
    if len(ac) < 2 or len(tol) < 2:
        return np.nan
    sd_ac = np.std(ac, ddof=1)
    sd_tol = np.std(tol, ddof=1)
    pooled = np.sqrt((sd_ac**2 + sd_tol**2) / 2)
    return np.nan if pooled == 0 or np.isnan(pooled) else float((np.mean(ac) - np.mean(tol)) / pooled)


def analyze_one(data, transcript_id, metric):
    d = data.loc[data["Transcript_ID_clean"].eq(transcript_id)].copy()
    d = d.dropna(subset=[metric, "Group", "sample"])
    if d.empty:
        return None

    ac = d.loc[d["Group"].eq("AC"), metric].to_numpy(dtype=float)
    tol = d.loc[d["Group"].eq("TOL"), metric].to_numpy(dtype=float)
    if len(ac) < 3 or len(tol) < 3:
        return None

    mean_ac = float(np.mean(ac))
    mean_tol = float(np.mean(tol))
    delta = mean_ac - mean_tol

    if np.all(ac == ac[0]) and np.all(tol == tol[0]) and ac[0] == tol[0]:
        t_p = np.nan
        mw_p = np.nan
    else:
        t_p = float(ttest_ind(ac, tol, equal_var=False, nan_policy="omit").pvalue)
        mw_p = float(mannwhitneyu(ac, tol, alternative="two-sided").pvalue)

    try:
        model = smf.ols(f"{metric} ~ C(Group, Treatment(reference='TOL'))", data=d).fit(cov_type="HC3")
        term = "C(Group, Treatment(reference='TOL'))[T.AC]"
        estimate = float(model.params.get(term, np.nan))
        pvalue = float(model.pvalues.get(term, np.nan))
        ci_low, ci_high = model.conf_int().loc[term].astype(float).tolist()
    except Exception:
        estimate, pvalue, ci_low, ci_high = [np.nan] * 4

    expressed_ac = int((d.loc[d["Group"].eq("AC"), "tpm"] >= 1).sum())
    expressed_tol = int((d.loc[d["Group"].eq("TOL"), "tpm"] >= 1).sum())
    return {
        "Transcript_ID_clean": transcript_id,
        "metric": metric,
        "gene_id": collapse_first(d["gencode_gene_id"]),
        "gene_name": collapse_first(d["gencode_gene_name"]),
        "protein_coding_class": collapse_first(d["protein_coding_class"]),
        "n": int(len(d)),
        "n_AC": int(len(ac)),
        "n_TOL": int(len(tol)),
        "n_AC_TPM_ge_1": expressed_ac,
        "n_TOL_TPM_ge_1": expressed_tol,
        "mean_AC": mean_ac,
        "mean_TOL": mean_tol,
        "median_AC": float(np.median(ac)),
        "median_TOL": float(np.median(tol)),
        "delta_AC_minus_TOL": delta,
        "cohen_d": cohen_d(ac, tol),
        "welch_p_value": t_p,
        "mannwhitney_p_value": mw_p,
        "estimate_adjusted_AC_minus_TOL": estimate,
        "ci_low": float(ci_low),
        "ci_high": float(ci_high),
        "p_value": pvalue,
        "mean_TPM_AC": float(d.loc[d["Group"].eq("AC"), "tpm"].mean()),
        "mean_TPM_TOL": float(d.loc[d["Group"].eq("TOL"), "tpm"].mean()),
        "mean_log2TPM_AC": float(d.loc[d["Group"].eq("AC"), "log2_TPM_plus_1"].mean()),
        "mean_log2TPM_TOL": float(d.loc[d["Group"].eq("TOL"), "log2_TPM_plus_1"].mean()),
        "KDmax": float(pd.to_numeric(d["KDmax"], errors="coerce").max()),
        "n_protein_rows_max_per_sample": int(d["n_protein_rows"].max()),
    }


def build_differential_stats(long_df):
    amy = long_df[long_df["is_amyloidogenic"]].copy()
    results = []
    for transcript_id in sorted(amy["Transcript_ID_clean"].dropna().unique()):
        for metric in METRICS:
            row = analyze_one(amy, transcript_id, metric)
            if row is not None:
                results.append(row)

    stats = pd.DataFrame(results)
    if stats.empty:
        return stats

    tiny = np.nextafter(0, 1)
    stats["p_value"] = stats["p_value"].replace(0, tiny)
    stats["FDR"] = np.nan
    for metric in METRICS:
        mask = stats["metric"].eq(metric)
        valid = mask & stats["p_value"].notna()
        if valid.any():
            _, fdr, _, _ = multipletests(stats.loc[valid, "p_value"], method="fdr_bh")
            stats.loc[valid, "FDR"] = fdr
    stats["FDR"] = stats["FDR"].replace(0, tiny)
    stats["direction"] = np.where(stats["delta_AC_minus_TOL"] > 0, "higher_in_AC", "higher_in_TOL")
    stats["abs_effect"] = stats["cohen_d"].abs()
    stats["minus_log10_p"] = -np.log10(stats["p_value"])
    stats["minus_log10_FDR"] = -np.log10(stats["FDR"])
    return stats.sort_values(["FDR", "abs_effect"], ascending=[True, False])


def metric_group_stats(sample_metrics):
    metric_map = {
        "amyloidogenic_burden_log2TPM": "Burden",
        "amyloidogenic_fraction_log2TPM": "Fraction",
        "kdmax_weighted_amyloidogenic_fraction_log2TPM": "KDmax-weighted Fraction",
    }
    rows = []
    for metric, label in metric_map.items():
        ac = sample_metrics.loc[sample_metrics["Group"].eq("AC"), metric].dropna().to_numpy(dtype=float)
        tol = sample_metrics.loc[sample_metrics["Group"].eq("TOL"), metric].dropna().to_numpy(dtype=float)
        p = mannwhitneyu(ac, tol, alternative="two-sided").pvalue if len(ac) and len(tol) else np.nan
        rows.append(
            {
                "endpoint": label,
                "metric": metric,
                "n_AC": int(len(ac)),
                "n_TOL": int(len(tol)),
                "mean_AC": float(np.mean(ac)) if len(ac) else np.nan,
                "mean_TOL": float(np.mean(tol)) if len(tol) else np.nan,
                "median_AC": float(np.median(ac)) if len(ac) else np.nan,
                "median_TOL": float(np.median(tol)) if len(tol) else np.nan,
                "delta_AC_minus_TOL": float(np.mean(ac) - np.mean(tol)) if len(ac) and len(tol) else np.nan,
                "mannwhitney_p_value": float(p) if np.isfinite(p) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def plot_sample_boxplots(sample_metrics):
    sns.set_theme(style="white", context="talk", font_scale=0.78)
    plot_metrics = [
        ("amyloidogenic_burden_log2TPM", "Burden\nΣ log2(TPM+1)"),
        ("amyloidogenic_fraction_log2TPM", "Fraction\nΣ log2(TPM+1) amyloidogenic / all protein-coding"),
        (
            "kdmax_weighted_amyloidogenic_fraction_log2TPM",
            "KDmax-weighted Fraction\nΣ(KDmax x log2(TPM+1)) amyloidogenic / all protein-coding",
        ),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 5.2), constrained_layout=True)
    for ax, (metric, title) in zip(axes, plot_metrics):
        d = sample_metrics.dropna(subset=[metric])
        sns.boxplot(
            data=d,
            x="Group",
            y=metric,
            order=GROUP_ORDER,
            palette=PALETTE,
            fliersize=0,
            width=0.48,
            linewidth=1.2,
            ax=ax,
        )
        sns.stripplot(
            data=d,
            x="Group",
            y=metric,
            order=GROUP_ORDER,
            palette=PALETTE,
            size=7,
            jitter=0.12,
            edgecolor="white",
            linewidth=0.8,
            ax=ax,
        )
        ax.set_title(title, fontsize=12, pad=12)
        ax.set_xlabel("")
        ax.set_ylabel("")
        sns.despine(ax=ax)
    png = OUT / "sample_level_log2TPM_metric_boxplots_AC_vs_TOL.png"
    pdf = OUT / "sample_level_log2TPM_metric_boxplots_AC_vs_TOL.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def plot_top_forest_grid(stats, top_n=20):
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)
    fig, axes = plt.subplots(1, len(METRICS), figsize=(23.0, max(7.2, top_n * 0.38)), sharey=False)
    for ax, metric in zip(axes, METRICS):
        d = stats.loc[stats["metric"].eq(metric)].dropna(subset=["estimate_adjusted_AC_minus_TOL", "ci_low", "ci_high"]).copy()
        d = d.sort_values(["FDR", "abs_effect"], ascending=[True, False]).head(top_n)
        d = d.sort_values("estimate_adjusted_AC_minus_TOL")
        labels = d["gene_name"].fillna("").astype(str) + " | " + d["Transcript_ID_clean"].astype(str)
        y = np.arange(len(d))
        colors = d["direction"].map(DIRECTION_PALETTE).fillna("#A0AEC0").to_numpy()
        alphas = np.where(d["FDR"].lt(0.10).fillna(False), 1.0, 0.48)
        for pos, (_, row) in enumerate(d.iterrows()):
            ax.errorbar(
                row["estimate_adjusted_AC_minus_TOL"],
                pos,
                xerr=[
                    [row["estimate_adjusted_AC_minus_TOL"] - row["ci_low"]],
                    [row["ci_high"] - row["estimate_adjusted_AC_minus_TOL"]],
                ],
                fmt="o",
                color=colors[pos],
                ecolor=colors[pos],
                alpha=alphas[pos],
                capsize=3,
                markersize=4.2,
            )
        ax.axvline(0, linestyle="--", color="#111827", linewidth=1)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=6.8)
        ax.set_xlabel("Adjusted delta AC - TOL")
        ax.set_title(METRIC_LABELS[metric], fontsize=10.5, pad=10)
        ax.grid(axis="x", alpha=0.28)
    axes[0].set_ylabel("Amyloidogenic protein-coding transcript")
    fig.suptitle("Top amyloidogenic proteoform metric differences: AC vs TOL", fontsize=14, y=1.02)
    fig.text(
        0.55,
        0.035,
        "Burden = log2(TPM+1); Fraction = log2(TPM+1) / all protein-coding log2 sum; "
        "KDmax-weighted Fraction = KDmax x log2(TPM+1) / all protein-coding weighted sum",
        ha="center",
        va="bottom",
        fontsize=9,
        color="#334155",
    )
    fig.subplots_adjust(left=0.11, right=0.985, bottom=0.16, top=0.88, wspace=0.88)
    png = OUT / "top_proteoform_metric_forest_grid_AC_vs_TOL.png"
    pdf = OUT / "top_proteoform_metric_forest_grid_AC_vs_TOL.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


def write_report(root, sample_group_stats, stats, sample_box_png, forest_png, sample_qc):
    top = stats.sort_values(["FDR", "abs_effect"], ascending=[True, False]).head(20)
    lines = [
        "# Proteoform log2(TPM+1) metric differences: AC vs TOL",
        "",
        f"- Zarr: `{ZARR_PATH.name}`",
        f"- Working class column: `{CLASS_COL}`",
        f"- Amyloidogenic class: `{AMY_CLASS}`",
        "- Denominator: all rows collapsed to transcripts with `is_gencode_protein_coding == True`.",
        "- TPM transform: `log2(TPM + 1)` for every burden/fraction endpoint.",
        f"- Zarr class attrs: `{root.attrs.get('protein_coding_consensus_class', {})}`",
        "",
        "## Sample-level endpoints",
        "",
        sample_group_stats.to_markdown(index=False),
        "",
        "## Top proteoform candidates",
        "",
        top[
            [
                "metric",
                "gene_name",
                "Transcript_ID_clean",
                "n_AC",
                "n_TOL",
                "delta_AC_minus_TOL",
                "cohen_d",
                "p_value",
                "FDR",
                "direction",
            ]
        ].to_markdown(index=False),
        "",
        "## QC",
        "",
        sample_qc.to_markdown(index=False),
        "",
        f"- Sample endpoint plot: `{sample_box_png.name}`",
        f"- Proteoform forest plot: `{forest_png.name}`",
    ]
    report = OUT / "proteoform_log2TPM_metric_differences_report.md"
    report.write_text("\n".join(lines), encoding="utf-8")
    return report


def main():
    root, long_df, sample_qc = load_long_table()
    long_df = attach_metadata(long_df)
    sample_metrics = build_sample_metrics(long_df)
    sample_group_stats = metric_group_stats(sample_metrics)
    stats = build_differential_stats(long_df)

    sample_metrics.to_csv(OUT / "sample_level_log2TPM_metrics.tsv", sep="\t", index=False)
    sample_group_stats.to_csv(OUT / "sample_level_log2TPM_metric_group_statistics.tsv", sep="\t", index=False)
    sample_qc.to_csv(OUT / "sample_level_zarr_merge_qc.tsv", sep="\t", index=False)
    long_df[
        [
            "sample",
            "Group",
            "Transcript_ID_clean",
            "gencode_gene_id",
            "gencode_gene_name",
            "protein_coding_class",
            "is_amyloidogenic",
            "tpm",
            "log2_TPM_plus_1",
            "KDmax",
            "kdmax_x_log2_TPM_plus_1",
            "burden_contribution_log2TPM",
            "fraction_contribution_log2TPM",
            "kdmax_weighted_fraction_contribution_log2TPM",
        ]
    ].to_csv(OUT / "proteoform_metric_contributions_long.tsv", sep="\t", index=False)
    stats.to_csv(OUT / "proteoform_metric_differential_statistics.tsv", sep="\t", index=False)
    stats.sort_values(["FDR", "abs_effect"], ascending=[True, False]).head(100).to_csv(
        OUT / "top100_candidate_proteoform_metric_differences.tsv", sep="\t", index=False
    )

    sample_box_png, sample_box_pdf = plot_sample_boxplots(sample_metrics)
    forest_png, forest_pdf = plot_top_forest_grid(stats)
    report = write_report(root, sample_group_stats, stats, sample_box_png, forest_png, sample_qc)

    summary = {
        "out_dir": str(OUT),
        "class_column": CLASS_COL,
        "n_samples": int(sample_metrics["sample"].nunique()),
        "n_AC": int(sample_metrics["Group"].eq("AC").sum()),
        "n_TOL": int(sample_metrics["Group"].eq("TOL").sum()),
        "n_amyloidogenic_transcripts_tested": int(stats["Transcript_ID_clean"].nunique()) if not stats.empty else 0,
        "n_tests": int(len(stats)),
        "n_FDR_lt_0_10": int(stats["FDR"].lt(0.10).sum()) if not stats.empty else 0,
        "sample_boxplot_png": str(sample_box_png),
        "sample_boxplot_pdf": str(sample_box_pdf),
        "forest_png": str(forest_png),
        "forest_pdf": str(forest_pdf),
        "report": str(report),
    }
    with open(OUT / "proteoform_log2TPM_metric_differences_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
