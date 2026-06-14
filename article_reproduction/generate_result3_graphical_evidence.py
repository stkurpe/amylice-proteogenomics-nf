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
import zarr
from scipy.stats import mannwhitneyu, spearmanr


ROOT = Path(__file__).resolve().parent
ZARR_PATH = ROOT / "project.zarr"
BASE = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili"
META_PATH = BASE / "sample_metadata_BC.tsv"
QC_PATH = BASE / "sample_level_amyloid_burden.tsv"
SCORES_PATH = BASE / "seven_signature_panel" / "signature_scores.tsv"
OUT = BASE / "figure3_graphical_evidence"
OUT.mkdir(parents=True, exist_ok=True)

CLASS_COL = "ProteinCoding_BothPredictors_Q17_Q83_Class"
AMY_CLASS = "Amyloidogenic"
NON_CLASS = "Non-Amyloidogenic"
MID_CLASS = "Intermediate"
CLASS_PRIORITY = {AMY_CLASS: 3, MID_CLASS: 2, NON_CLASS: 1}
PALETTE = {"AC": "#2B6CB0", "TOL": "#C2410C"}
LOW_SCORE_AC_EXCLUDE = {"SRR32060215": "low AC amyloidogenic score outlier"}

BURDEN = [
    "pc_q17q83_amyloidogenic_log2TPM_fraction",
    "pc_q17q83_amyloidogenic_log2TPM_burden",
    "pc_q17q83_kdmax_weighted_amyloidogenic_log2TPM_fraction",
    "n_pc_q17q83_amyloidogenic_TPM_ge_1",
    "pc_q17q83_mean_two_predictor_probability",
]
BURDEN_LABELS = {
    "pc_q17q83_amyloidogenic_log2TPM_fraction": "Amyloidogenic\nlog2(TPM+1) fraction",
    "pc_q17q83_amyloidogenic_log2TPM_burden": "Amyloidogenic burden\nsum log2(TPM+1)",
    "pc_q17q83_kdmax_weighted_amyloidogenic_log2TPM_fraction": "KDmax-weighted\namyloidogenic log2(TPM+1) fraction",
    "n_pc_q17q83_amyloidogenic_TPM_ge_1": "N amyloidogenic\nTPM >= 1",
    "pc_q17q83_mean_two_predictor_probability": "Mean two-predictor\nprobability",
}
SIGS = ["UPR_score", "Inflammatory_score", "IFNG_score", "IFNA_score", "Myeloid_score"]
SIG_LABELS = {
    "UPR_score": "UPR",
    "Inflammatory_score": "Inflammatory",
    "IFNG_score": "IFN-gamma",
    "IFNA_score": "IFN-alpha",
    "Myeloid_score": "Myeloid",
}


def zarr_df(group):
    return pd.DataFrame({col: group[col][:] for col in group.attrs["columns"]})


def clean_transcript(value):
    match = re.search(r"(ENST\d+)(?:\.\d+)?", str(value))
    return match.group(1) if match else pd.NA


def collapse_class(values):
    clean = pd.Series(values, dtype="string").replace("", pd.NA).dropna()
    if clean.empty:
        return pd.NA
    return max(clean, key=lambda x: CLASS_PRIORITY.get(str(x), 0))


def cliffs_delta(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan
    stat, _ = mannwhitneyu(x, y, alternative="two-sided")
    return float(2 * stat / (len(x) * len(y)) - 1)


def as_bool(value):
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def bootstrap_delta_ci(a, b, n=10000, seed=20260613):
    rng = np.random.default_rng(seed)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) == 0 or len(b) == 0:
        return np.nan, np.nan
    deltas = []
    for _ in range(n):
        aa = rng.choice(a, size=len(a), replace=True)
        bb = rng.choice(b, size=len(b), replace=True)
        deltas.append(np.mean(bb) - np.mean(aa))
    return tuple(np.percentile(deltas, [2.5, 97.5]))


def load_q17q83_protein_coding_data():
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    all_frames = []
    frames = []
    expr_frames = []
    feature_frames = []
    for sample in root["samples"].group_keys():
        amy = zarr_df(root[f"samples/{sample}/amyloid"])
        amy[CLASS_COL] = amy[CLASS_COL].astype("string").replace("", pd.NA)
        amy["Transcript_ID_clean"] = amy["Transcript_ID_clean"].astype("string").replace("", pd.NA)
        amy["gencode_gene_name"] = amy["gencode_gene_name"].astype("string").replace("", pd.NA)
        all_frames.append(amy.copy())
        frames.append(amy[amy[CLASS_COL].notna()].copy())

        expr = zarr_df(root[f"samples/{sample}/expression"])
        expr["Transcript_ID_clean"] = expr["transcript_id"].map(clean_transcript)
        expr_frames.append(expr[["sample", "Transcript_ID_clean", "tpm"]])

        features = zarr_df(root[f"samples/{sample}/protein_features"])
        features["Sequence_ID"] = features["protein_id"].astype(str)
        features["Transcript_ID_clean"] = features["Sequence_ID"].map(clean_transcript)
        feature_frames.append(features[["sample", "Sequence_ID", "Transcript_ID_clean", "KD_max"]])

    amyloid_rows = pd.concat(all_frames, ignore_index=True)
    classified_rows = pd.concat(frames, ignore_index=True)
    expression = pd.concat(expr_frames, ignore_index=True)
    features = pd.concat(feature_frames, ignore_index=True)

    class_features = classified_rows[["sample", "Sequence_ID", "Transcript_ID_clean", CLASS_COL]].merge(
        features,
        on=["sample", "Sequence_ID", "Transcript_ID_clean"],
        how="left",
    )
    collapsed_features = (
        class_features.groupby(["sample", "Transcript_ID_clean"], dropna=False)
        .agg(KDmax_max=("KD_max", "max"), KDmax_mean=("KD_max", "mean"))
        .reset_index()
    )

    collapsed = (
        classified_rows.groupby(["sample", "Transcript_ID_clean"], dropna=False)
        .agg(
            ProteinCoding_BothPredictors_Q17_Q83_Class=(CLASS_COL, collapse_class),
            gencode_gene_name=("gencode_gene_name", lambda x: x.dropna().iloc[0] if x.dropna().size else pd.NA),
            gencode_gene_id=("gencode_gene_id", lambda x: x.dropna().iloc[0] if x.dropna().size else pd.NA),
            mean_AMYPred_Prob=("AMYPred_Prob", "mean"),
            mean_AmyloGramPy_Prob=("AmyloGramPy_Prob", "mean"),
            max_AMYPred_Prob=("AMYPred_Prob", "max"),
            max_AmyloGramPy_Prob=("AmyloGramPy_Prob", "max"),
            n_protein_variants=("Sequence_ID", "nunique"),
        )
        .reset_index()
    )
    collapsed["mean_two_predictor_probability"] = collapsed[
        ["mean_AMYPred_Prob", "mean_AmyloGramPy_Prob"]
    ].mean(axis=1)

    expr_amy = expression.merge(collapsed, on=["sample", "Transcript_ID_clean"], how="inner")
    expr_amy = expr_amy.merge(collapsed_features, on=["sample", "Transcript_ID_clean"], how="left")
    expr_amy["tpm"] = pd.to_numeric(expr_amy["tpm"], errors="coerce").fillna(0.0)
    expr_amy["KDmax_max"] = pd.to_numeric(expr_amy["KDmax_max"], errors="coerce").fillna(0.0)
    expr_amy["log2_tpm_plus_1"] = np.log2(expr_amy["tpm"] + 1.0)
    expr_amy["kdmax_x_tpm"] = expr_amy["KDmax_max"] * expr_amy["tpm"]
    expr_amy["kdmax_x_log2_tpm_plus_1"] = expr_amy["KDmax_max"] * expr_amy["log2_tpm_plus_1"]
    return root, amyloid_rows, expr_amy


def sample_burden(expr_amy, meta):
    rows = []
    for sample, g in expr_amy.groupby("sample", observed=True):
        pc_total = g["tpm"].sum()
        pc_log_total = g["log2_tpm_plus_1"].sum()
        pc_kdmax_tpm_total = g["kdmax_x_tpm"].sum()
        pc_kdmax_log_total = g["kdmax_x_log2_tpm_plus_1"].sum()
        amy = g[g[CLASS_COL].eq(AMY_CLASS)].copy()
        amy_tpm = amy["tpm"].sum()
        amy_log_burden = amy["log2_tpm_plus_1"].sum()
        amy_kdmax_tpm = amy["kdmax_x_tpm"].sum()
        amy_kdmax_log = amy["kdmax_x_log2_tpm_plus_1"].sum()
        expressed = amy[amy["tpm"].ge(1)]
        rows.append(
            {
                "sample": sample,
                "pc_all_protein_coding_TPM": pc_total,
                "pc_all_protein_coding_log2TPM_sum": pc_log_total,
                "pc_all_protein_coding_KDmax_TPM_sum": pc_kdmax_tpm_total,
                "pc_all_protein_coding_KDmax_log2TPM_sum": pc_kdmax_log_total,
                "pc_q17q83_amyloidogenic_TPM": amy_tpm,
                "pc_q17q83_amyloidogenic_TPM_fraction": amy_tpm / pc_total if pc_total else np.nan,
                "pc_q17q83_amyloidogenic_log2TPM_burden": amy_log_burden,
                "pc_q17q83_amyloidogenic_log2TPM_fraction": amy_log_burden / pc_log_total if pc_log_total else np.nan,
                "pc_q17q83_kdmax_weighted_amyloidogenic_TPM": amy_kdmax_tpm,
                "pc_q17q83_kdmax_weighted_amyloidogenic_TPM_fraction": amy_kdmax_tpm / pc_kdmax_tpm_total
                if pc_kdmax_tpm_total
                else np.nan,
                "pc_q17q83_kdmax_weighted_amyloidogenic_log2TPM": amy_kdmax_log,
                "pc_q17q83_kdmax_weighted_amyloidogenic_log2TPM_fraction": amy_kdmax_log / pc_kdmax_log_total
                if pc_kdmax_log_total
                else np.nan,
                "n_pc_q17q83_amyloidogenic_TPM_ge_1": int(len(expressed)),
                "pc_q17q83_mean_two_predictor_probability": float(
                    np.average(amy["mean_two_predictor_probability"], weights=amy["tpm"])
                )
                if amy_tpm > 0
                else np.nan,
            }
        )
    burden = pd.DataFrame(rows)
    out = meta.merge(burden, on="sample", how="inner")
    out = out[out["time"].eq("BC") & out["Group"].isin(["AC", "TOL"])].copy()
    if "keep_coverage_qc" in out.columns:
        out = out[out["keep_coverage_qc"].map(as_bool)].copy()
    out = out[~out["sample"].astype(str).isin(LOW_SCORE_AC_EXCLUDE)].copy()
    out["Group"] = pd.Categorical(out["Group"], categories=["AC", "TOL"], ordered=True)
    return out


def group_stats(df, metric):
    a = df.loc[df["Group"].eq("AC"), metric].dropna().to_numpy(dtype=float)
    b = df.loc[df["Group"].eq("TOL"), metric].dropna().to_numpy(dtype=float)
    if len(a) and len(b):
        _, p = mannwhitneyu(a, b, alternative="two-sided")
        delta = float(np.mean(b) - np.mean(a))
        ci_low, ci_high = bootstrap_delta_ci(a, b)
        cliff = cliffs_delta(b, a)
    else:
        p, delta, ci_low, ci_high, cliff = [np.nan] * 5
    return {
        "n_AC": int(len(a)),
        "n_TOL": int(len(b)),
        "mean_AC": float(np.mean(a)) if len(a) else np.nan,
        "mean_TOL": float(np.mean(b)) if len(b) else np.nan,
        "median_AC": float(np.median(a)) if len(a) else np.nan,
        "median_TOL": float(np.median(b)) if len(b) else np.nan,
        "delta_mean_TOL_minus_AC": delta,
        "bootstrap_ci_low": float(ci_low),
        "bootstrap_ci_high": float(ci_high),
        "mannwhitney_p": float(p),
        "cliffs_delta_TOL_vs_AC": float(cliff),
    }


def draw_figure(root, amyloid_rows, expr_amy, burden, merged, contrib, primary_stats):
    sns.set_theme(style="white", context="talk", font_scale=0.72)
    fig = plt.figure(figsize=(14, 11), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[0.75, 1.3, 1.15], width_ratios=[1, 1])

    ax_a = fig.add_subplot(gs[0, :])
    ax_a.axis("off")
    ax_a.text(-0.13, 1.12, "A", transform=ax_a.transAxes, fontsize=15, fontweight="bold", va="top")
    steps = [
        "Baseline\nblood RNA-seq",
        "Amylice\nproteoform reconstruction",
        "GENCODE protein-coding\nQ17/Q83 consensus class",
        "Expression-weighted\namyloidogenic burden",
        "AC/control vs\nChILI/TOL",
        "Inflammatory and\nproteostasis context",
    ]
    xs = np.linspace(0.06, 0.94, len(steps))
    for i, (x, text) in enumerate(zip(xs, steps)):
        ax_a.text(
            x,
            0.52,
            text,
            ha="center",
            va="center",
            fontsize=10,
            color="#111827",
            bbox=dict(boxstyle="round,pad=0.38", fc="#F8FAFC", ec="#64748B", lw=1.2),
        )
        if i < len(steps) - 1:
            ax_a.annotate(
                "",
                xy=(xs[i + 1] - 0.07, 0.52),
                xytext=(x + 0.07, 0.52),
                arrowprops=dict(arrowstyle="->", lw=1.3, color="#334155"),
            )
    class_series = pd.Series(amyloid_rows[CLASS_COL].astype("string").replace("", pd.NA)).fillna("NA")
    class_counts = (
        class_series.value_counts()
        .reindex([AMY_CLASS, MID_CLASS, NON_CLASS, "NA"])
        .fillna(0)
        .astype(int)
    )
    attr = root.attrs.get("protein_coding_consensus_class", {})
    ax_a.text(
        0.5,
        0.09,
        f"Working class: {attr.get('column', CLASS_COL)}; "
        f"{AMY_CLASS}={class_counts[AMY_CLASS]:,}, {MID_CLASS}={class_counts[MID_CLASS]:,}, "
        f"{NON_CLASS}={class_counts[NON_CLASS]:,}, NA={class_counts['NA']:,} rows.".replace(",", " "),
        ha="center",
        va="center",
        fontsize=9.5,
        color="#475569",
    )

    ax_b = fig.add_subplot(gs[1, 0])
    metric = "pc_q17q83_amyloidogenic_log2TPM_fraction"
    plot_df = burden.dropna(subset=[metric])
    sns.boxplot(
        data=plot_df,
        x="Group",
        y=metric,
        order=["AC", "TOL"],
        palette=PALETTE,
        width=0.46,
        fliersize=0,
        linewidth=1.2,
        ax=ax_b,
    )
    sns.stripplot(
        data=plot_df,
        x="Group",
        y=metric,
        order=["AC", "TOL"],
        palette=PALETTE,
        size=7,
        jitter=0.12,
        edgecolor="white",
        linewidth=0.8,
        ax=ax_b,
    )
    ax_b.set_title("B. Protein-coding Q17/Q83 amyloidogenic burden", fontsize=12, pad=14)
    ax_b.set_ylabel("Amyloidogenic fraction, sum log2(TPM+1)")
    ax_b.set_xlabel("")
    ax_b.text(
        0.02,
        0.98,
        f"TOL-AC delta = {primary_stats['delta_mean_TOL_minus_AC']:.3f}\n"
        f"95% bootstrap CI [{primary_stats['bootstrap_ci_low']:.3f}, {primary_stats['bootstrap_ci_high']:.3f}]\n"
        f"MW p = {primary_stats['mannwhitney_p']:.3f}; Cliff's delta = {primary_stats['cliffs_delta_TOL_vs_AC']:.2f}",
        transform=ax_b.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#CBD5E1"),
    )

    ax_c = fig.add_subplot(gs[1, 1])
    top = contrib.sort_values("mean_log2TPM_contribution", ascending=True).tail(12)
    sizes = 35 + 30 * top["n_samples_expressed_TPM_ge_1"]
    sc = ax_c.scatter(
        top["mean_log2TPM_contribution"],
        top["label"],
        s=sizes,
        c=top["mean_two_predictor_probability"],
        cmap="viridis",
        vmin=0,
        vmax=1,
        edgecolor="#1F2937",
        linewidth=0.5,
    )
    ax_c.set_title("C. Top amyloidogenic burden-driving proteoforms", fontsize=12, pad=14)
    ax_c.set_xlabel("Mean log2(TPM+1) contribution")
    ax_c.set_ylabel("")
    ax_c.grid(axis="x", color="#E2E8F0", linewidth=0.8)
    cbar = fig.colorbar(sc, ax=ax_c, shrink=0.85, pad=0.02)
    cbar.set_label("Mean two-predictor probability")

    ax_d = fig.add_subplot(gs[2, 0])
    corr = pd.DataFrame(index=[BURDEN_LABELS[x] for x in BURDEN], columns=[SIG_LABELS[x] for x in SIGS], dtype=float)
    annot = corr.copy().astype(str)
    for burden_metric in BURDEN:
        for sig in SIGS:
            d = merged[[burden_metric, sig]].dropna()
            rho, _ = spearmanr(d[burden_metric], d[sig]) if len(d) >= 4 else (np.nan, np.nan)
            corr.loc[BURDEN_LABELS[burden_metric], SIG_LABELS[sig]] = rho
            annot.loc[BURDEN_LABELS[burden_metric], SIG_LABELS[sig]] = "" if np.isnan(rho) else f"{rho:.2f}"
    sns.heatmap(
        corr,
        vmin=-1,
        vmax=1,
        cmap="vlag",
        center=0,
        annot=annot,
        fmt="",
        linewidths=0.7,
        linecolor="white",
        cbar_kws={"label": "Spearman rho"},
        ax=ax_d,
    )
    ax_d.set_title("D. Q17/Q83 burden coupling with transcriptional signatures", fontsize=12, pad=14)
    ax_d.set_xlabel("")
    ax_d.set_ylabel("")
    ax_d.tick_params(axis="x", rotation=35)
    ax_d.tick_params(axis="y", rotation=0)

    ax_e = fig.add_subplot(gs[2, 1])
    xmetric = "pc_q17q83_amyloidogenic_log2TPM_fraction"
    ymetric = "Inflammatory_score"
    d = merged.dropna(subset=[xmetric, ymetric])
    for group, sub in d.groupby("Group", observed=True):
        ax_e.scatter(
            sub[xmetric],
            sub[ymetric],
            s=70,
            color=PALETTE[str(group)],
            edgecolor="white",
            linewidth=0.9,
            label=f"{group} (n={len(sub)})",
            zorder=3,
        )
    if len(d) >= 4:
        x = d[xmetric].to_numpy()
        y = d[ymetric].to_numpy()
        fit = np.polyfit(x, y, deg=1)
        xs_line = np.linspace(x.min(), x.max(), 100)
        ax_e.plot(xs_line, fit[0] * xs_line + fit[1], color="#111827", lw=1.4)
        rho, p = spearmanr(x, y)
    else:
        rho, p = (np.nan, np.nan)
    ax_e.axhline(0, color="#CBD5E1", lw=0.9, zorder=1)
    ax_e.set_title("E. Representative burden-signature association", fontsize=12, pad=14)
    ax_e.set_xlabel("Amyloidogenic fraction, sum log2(TPM+1)")
    ax_e.set_ylabel("Inflammatory signature score")
    ax_e.legend(frameon=False, loc="best")
    ax_e.text(
        0.02,
        0.98,
        f"QC baseline samples: Spearman rho = {rho:.2f}, p = {p:.3f}",
        transform=ax_e.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#CBD5E1"),
    )

    for ax in [ax_b, ax_c, ax_e]:
        sns.despine(ax=ax)

    fig.suptitle(
        "Figure 3. Baseline protein-coding amyloidogenic burden and inflammatory coupling in ChILI/TOL samples",
        fontsize=15,
        fontweight="bold",
    )
    png = OUT / "figure3_baseline_pc_q17q83_amyloidogenic_burden_evidence.png"
    pdf = OUT / "figure3_baseline_pc_q17q83_amyloidogenic_burden_evidence.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    return png, pdf, rho, p


def draw_boxplot_figure(burden):
    sns.set_theme(style="white", context="talk", font_scale=0.82)
    metrics = [
        (
            "pc_q17q83_amyloidogenic_log2TPM_burden",
            "Amyloidogenic burden\nΣ log2(TPM+1)",
            "Burden",
        ),
        (
            "pc_q17q83_amyloidogenic_log2TPM_fraction",
            "Amyloidogenic fraction\nΣ log2(TPM+1) amyloidogenic / Σ log2(TPM+1) all protein-coding",
            "Fraction",
        ),
        (
            "pc_q17q83_kdmax_weighted_amyloidogenic_log2TPM_fraction",
            "KDmax-weighted amyloidogenic fraction\nΣ(KDmax×log2(TPM+1)) amyloidogenic / Σ(KDmax×log2(TPM+1)) all protein-coding",
            "KDmax-weighted Fraction",
        ),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15.8, 5.2), constrained_layout=True)
    stats = {}
    for ax, (metric, ylabel, title) in zip(axes, metrics):
        plot_df = burden.dropna(subset=[metric])
        sns.boxplot(
            data=plot_df,
            x="Group",
            y=metric,
            order=["AC", "TOL"],
            palette=PALETTE,
            width=0.48,
            fliersize=0,
            linewidth=1.25,
            ax=ax,
        )
        sns.stripplot(
            data=plot_df,
            x="Group",
            y=metric,
            order=["AC", "TOL"],
            palette=PALETTE,
            size=7,
            jitter=0.13,
            edgecolor="white",
            linewidth=0.8,
            ax=ax,
        )
        metric_stats = group_stats(burden, metric)
        stats[metric] = metric_stats
        ax.set_title(title, fontsize=13, pad=12)
        ax.set_xlabel("")
        ax.set_ylabel(ylabel, fontsize=11)
        ax.text(
            0.02,
            0.98,
            f"Δ TOL-AC = {metric_stats['delta_mean_TOL_minus_AC']:.3g}\n"
            f"MW p = {metric_stats['mannwhitney_p']:.3f}\n"
            f"Cliff's δ = {metric_stats['cliffs_delta_TOL_vs_AC']:.2f}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9.5,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#CBD5E1"),
        )
        sns.despine(ax=ax)
    png = OUT / "pc_q17q83_log2TPM_kdmax_burden_boxplots.png"
    pdf = OUT / "pc_q17q83_log2TPM_kdmax_burden_boxplots.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    return png, pdf, stats


def main():
    root, amyloid_rows, expr_amy = load_q17q83_protein_coding_data()
    meta = pd.read_csv(META_PATH, sep="\t")
    if QC_PATH.exists():
        qc = pd.read_csv(
            QC_PATH,
            sep="\t",
            usecols=["sample", "protein_coding_annotation_fraction", "keep_coverage_qc"],
        )
        meta = meta.drop(columns=["protein_coding_annotation_fraction", "keep_coverage_qc"], errors="ignore").merge(
            qc, on="sample", how="left"
        )
    burden = sample_burden(expr_amy, meta)
    scores = pd.read_csv(SCORES_PATH, sep="\t")
    merged = scores.merge(burden[["sample"] + BURDEN], on="sample", how="inner")
    merged = merged[merged["time"].eq("BC") & merged["Group"].isin(["AC", "TOL"])].copy()
    merged["Group"] = pd.Categorical(merged["Group"], categories=["AC", "TOL"], ordered=True)

    contrib_source = expr_amy[expr_amy[CLASS_COL].eq(AMY_CLASS)].copy()
    contrib_source = contrib_source.merge(meta[["sample", "Group", "time"]], on="sample", how="left")
    contrib_source = contrib_source[contrib_source["time"].eq("BC") & contrib_source["Group"].isin(["AC", "TOL"])].copy()
    contrib = (
        contrib_source.groupby(["Transcript_ID_clean", "gencode_gene_name"], dropna=False)
        .agg(
            mean_TPM_contribution=("tpm", "mean"),
            mean_log2TPM_contribution=("log2_tpm_plus_1", "mean"),
            median_TPM_contribution=("tpm", "median"),
            n_samples_expressed_TPM_ge_1=("tpm", lambda x: int((x >= 1).sum())),
            mean_two_predictor_probability=("mean_two_predictor_probability", "mean"),
        )
        .reset_index()
    )
    contrib["gene_label"] = contrib["gencode_gene_name"].fillna(contrib["Transcript_ID_clean"])
    contrib["label"] = contrib["gene_label"] + " | " + contrib["Transcript_ID_clean"].str.replace("ENST", "ENST ", regex=False)

    primary_metric = "pc_q17q83_amyloidogenic_log2TPM_fraction"
    primary_stats = group_stats(burden, primary_metric)
    png, pdf, rho, p = draw_figure(root, amyloid_rows, expr_amy, burden, merged, contrib, primary_stats)
    boxplot_png, boxplot_pdf, boxplot_stats = draw_boxplot_figure(burden)

    class_counts = pd.Series(amyloid_rows[CLASS_COL].astype("string").replace("", pd.NA)).fillna("NA").value_counts().to_dict()
    baseline_meta = meta[meta["time"].eq("BC") & meta["Group"].isin(["AC", "TOL"])].copy()
    classified_samples = set(expr_amy["sample"].astype(str))
    burden_samples = set(burden["sample"].astype(str))
    excluded = []
    for _, row in baseline_meta.iterrows():
        sample = str(row["sample"])
        if sample in burden_samples:
            continue
        reason = []
        if sample in LOW_SCORE_AC_EXCLUDE:
            reason.append(LOW_SCORE_AC_EXCLUDE[sample])
        if sample not in classified_samples:
            reason.append("no non-empty ProteinCoding_BothPredictors_Q17_Q83_Class rows")
        if "keep_coverage_qc" in row and not as_bool(row["keep_coverage_qc"]):
            reason.append("protein_coding_annotation_fraction < 0.5")
        excluded.append({"sample": sample, "Group": row["Group"], "reason": "; ".join(reason) or "excluded"})
    summary = {
        "figure_png": str(png),
        "figure_pdf": str(pdf),
        "boxplot_png": str(boxplot_png),
        "boxplot_pdf": str(boxplot_pdf),
        "class_column": CLASS_COL,
        "class_counts_rows": {str(k): int(v) for k, v in class_counts.items()},
        "n_baseline_samples": int(len(burden)),
        "n_AC": int((burden["Group"] == "AC").sum()),
        "n_TOL": int((burden["Group"] == "TOL").sum()),
        "excluded_baseline_samples": excluded,
        "primary_metric": primary_metric,
        **primary_stats,
        "boxplot_stats": boxplot_stats,
        "scatter_spearman_rho": float(rho),
        "scatter_spearman_p": float(p),
    }
    burden.to_csv(OUT / "sample_level_pc_q17q83_amyloidogenic_burden.tsv", sep="\t", index=False)
    contrib.sort_values("mean_TPM_contribution", ascending=False).to_csv(
        OUT / "top_pc_q17q83_amyloidogenic_contributors.tsv", sep="\t", index=False
    )
    with open(OUT / "figure3_pc_q17q83_graphical_evidence_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
