import json
import math
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
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parent
ZARR_PATH = Path("/Users/user841/Projects/BIoinfServer/Results_analysis/bioinfo_zarr/validation_project.zarr")
META_PATH = ROOT / "GSE287540_SraRunTable.csv"
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "validation" / "q17q83_probability_groups"

GENCODE_CANDIDATES = [
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "competitive_gene_set_test_gencode_v48_protein_coding" / "gencode_transcript_annotation_used.tsv",
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "rank_based_amyloid_burden_signature" / "gencode_v48_transcript_annotation.tsv",
]

VALIDATION_SAMPLES = [
    "SRR32060218",
    "SRR32060222",
    "SRR32060224",
    "SRR32060226",
    "SRR32060252",
    "SRR32060255",
    "SRR32060299",
    "SRR32060304",
]

GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
TIMEPOINT = "BC"
CLASS_ORDER = ["Non-Amyloidogenic", "Intermediate", "Amyloidogenic"]
EXTREME_CLASS_ORDER = ["Non-Amyloidogenic", "Amyloidogenic"]
CLASS_PALETTE = {
    "Non-Amyloidogenic": "#4C78A8",
    "Intermediate": "#7A7A7A",
    "Amyloidogenic": "#C85245",
}
GROUP_PALETTE = {"AC": "#3B6EA8", "TOL": "#B24A3B"}
PREDICTOR_LABELS = {
    "AMYPred-FRL": "AMYPred-FRL",
    "AmyloGramPy": "AmyloGramPy",
    "Consensus": "Two-predictor consensus",
}
METRIC_LABELS = {
    "absolute_burden": "Absolute burden\nsum log2(TPM+1)",
    "normalized_fraction": "Normalized fraction",
    "KDmax_weighted_fraction": "KDmax-weighted fraction",
}
SIGNATURES = {
    "UPR_score": "UPR",
    "Inflammatory_score": "Inflammatory",
    "IFNG_score": "Interferon-gamma",
    "IFNA_score": "Interferon-alpha",
    "Myeloid_score": "Myeloid",
}
TARGET_TRANSCRIPTS = {
    "OR10J1 | ENST00000642080": "ENST00000642080",
    "SRPK1 | ENST00000373825": "ENST00000373825",
    "TP53BP1 | ENST00000382044": "ENST00000382044",
}
TARGET_GENES = ["OR10J1", "SRPK1", "TP53BP1"]

THRESHOLDS = {
    "AMYPred-FRL": {"low": 0.1206, "high": 0.8802, "prob_col": "AMYPred_Prob"},
    "AmyloGramPy": {"low": 0.6947, "high": 0.8898, "prob_col": "AmyloGramPy_Prob"},
}


def clean_enst(x):
    if pd.isna(x):
        return np.nan
    m = re.search(r"ENST\d+(?:\.\d+)?", str(x))
    return m.group(0).split(".")[0] if m else np.nan


def zarr_group_to_df(group):
    return pd.DataFrame({col: group[col][:] for col in group.array_keys()})


def classify_probability(prob, low, high):
    if pd.isna(prob):
        return np.nan
    if prob <= low:
        return "Non-Amyloidogenic"
    if prob >= high:
        return "Amyloidogenic"
    return "Intermediate"


def consensus_class(amypred_prob, amylogrampy_prob):
    if pd.isna(amypred_prob) or pd.isna(amylogrampy_prob):
        return np.nan
    if amypred_prob >= THRESHOLDS["AMYPred-FRL"]["high"] and amylogrampy_prob >= THRESHOLDS["AmyloGramPy"]["high"]:
        return "Amyloidogenic"
    if amypred_prob <= THRESHOLDS["AMYPred-FRL"]["low"] and amylogrampy_prob <= THRESHOLDS["AmyloGramPy"]["low"]:
        return "Non-Amyloidogenic"
    return "Intermediate"


def collapse_class(values):
    vals = set(pd.Series(values).dropna().astype(str))
    for label in ["Amyloidogenic", "Intermediate", "Non-Amyloidogenic"]:
        if label in vals:
            return label
    return np.nan


def load_gencode_pc():
    needed = {"Transcript_ID_clean", "gene_id", "gene_name", "biotype"}
    for path in GENCODE_CANDIDATES:
        if path.exists():
            ann = pd.read_csv(path, sep="\t")
            if needed.issubset(ann.columns):
                ann = ann[list(needed)].drop_duplicates("Transcript_ID_clean")
                return ann[ann["biotype"].eq("protein_coding")].copy(), path
    raise FileNotFoundError("No GENCODE protein-coding annotation table was found.")


def load_sample_metadata():
    meta = pd.read_csv(META_PATH)
    meta["Run"] = meta["Run"].astype(str)
    meta["raw_Group"] = meta["Group"].astype(str)
    meta["Group"] = meta["raw_Group"].map(GROUP_MAP)
    meta = meta[
        meta["Run"].isin(VALIDATION_SAMPLES)
        & meta["time"].eq(TIMEPOINT)
        & meta["Group"].isin(["AC", "TOL"])
    ].rename(columns={"Run": "sample"})
    return meta.set_index("sample").loc[VALIDATION_SAMPLES].reset_index()


def load_expression_matrix(root, samples):
    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    sample_idx = [zarr_samples.index(s) for s in samples]
    tx = pd.Series(root["layers/expression/transcript_ids"][:]).astype(str).map(clean_enst)
    tpm = pd.DataFrame(
        root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=tx,
        columns=samples,
    )
    return tpm[~tpm.index.isna()].groupby(level=0).sum()


def build_collapsed_probability_classes(root, samples, gencode_pc):
    frames = []
    for sample in samples:
        amy = zarr_group_to_df(root[f"samples/{sample}/amyloid"])
        pf = zarr_group_to_df(root[f"samples/{sample}/protein_features"])[["protein_id", "KD_max"]]
        amy = amy.merge(pf, left_on="Sequence_ID", right_on="protein_id", how="left")
        amy["sample"] = sample
        amy["Transcript_ID_clean"] = amy["Sequence_ID"].map(clean_enst)
        amy["AMYPred-FRL_class"] = amy["AMYPred_Prob"].map(
            lambda x: classify_probability(x, THRESHOLDS["AMYPred-FRL"]["low"], THRESHOLDS["AMYPred-FRL"]["high"])
        )
        amy["AmyloGramPy_class"] = amy["AmyloGramPy_Prob"].map(
            lambda x: classify_probability(x, THRESHOLDS["AmyloGramPy"]["low"], THRESHOLDS["AmyloGramPy"]["high"])
        )
        amy["Consensus_class"] = [
            consensus_class(a, b) for a, b in zip(amy["AMYPred_Prob"], amy["AmyloGramPy_Prob"])
        ]
        frames.append(amy)

    raw = pd.concat(frames, ignore_index=True)
    raw = raw.dropna(subset=["Transcript_ID_clean"])
    raw = raw.merge(gencode_pc, on="Transcript_ID_clean", how="inner")
    raw.to_csv(OUT / "protein_coding_raw_probability_classes.tsv", sep="\t", index=False)

    collapsed = (
        raw.groupby(["sample", "Transcript_ID_clean"], as_index=False)
        .agg(
            AMYPred_FRL_class=("AMYPred-FRL_class", collapse_class),
            AmyloGramPy_class=("AmyloGramPy_class", collapse_class),
            Consensus_class=("Consensus_class", collapse_class),
            AMYPred_Prob_max=("AMYPred_Prob", "max"),
            AmyloGramPy_Prob_max=("AmyloGramPy_Prob", "max"),
            KDmax=("KD_max", "max"),
            n_protein_rows=("Sequence_ID", "nunique"),
            gene_id=("gene_id", "first"),
            gene_name=("gene_name", "first"),
            biotype=("biotype", "first"),
        )
    )
    collapsed.to_csv(OUT / "protein_coding_collapsed_transcript_probability_classes.tsv", sep="\t", index=False)
    return collapsed


def compute_sample_class_metrics(tpm, collapsed, sample_meta, protein_coding_transcripts):
    class_cols = {
        "AMYPred-FRL": "AMYPred_FRL_class",
        "AmyloGramPy": "AmyloGramPy_class",
        "Consensus": "Consensus_class",
    }
    rows = []
    contribution_rows = []

    for sample in sample_meta["sample"]:
        pc_tpm = tpm.loc[tpm.index.intersection(protein_coding_transcripts), sample]
        ann = collapsed[collapsed["sample"].eq(sample)].set_index("Transcript_ID_clean")
        ids = pc_tpm.index.intersection(ann.index)
        ann = ann.loc[ids].copy()
        ann.index.name = "Transcript_ID_clean"
        ann["tpm"] = pc_tpm.loc[ids].astype(float)
        ann["log2_TPM_plus_1"] = np.log2(ann["tpm"].fillna(0) + 1)
        ann["KDmax"] = pd.to_numeric(ann["KDmax"], errors="coerce")
        denominator = float(ann["log2_TPM_plus_1"].sum())
        kd_denominator = float((ann["KDmax"] * ann["log2_TPM_plus_1"]).sum(skipna=True))

        for predictor, col in class_cols.items():
            for cls in CLASS_ORDER:
                mask = ann[col].eq(cls)
                absolute = float(ann.loc[mask, "log2_TPM_plus_1"].sum())
                kd_num = float((ann.loc[mask, "KDmax"] * ann.loc[mask, "log2_TPM_plus_1"]).sum(skipna=True))
                rows.append(
                    {
                        "sample": sample,
                        "predictor": predictor,
                        "probability_class": cls,
                        "absolute_burden": absolute,
                        "normalized_fraction": absolute / denominator if denominator else np.nan,
                        "KDmax_weighted_fraction": kd_num / kd_denominator if kd_denominator else np.nan,
                        "n_transcripts": int(mask.sum()),
                        "denominator_log2_TPM_plus_1": denominator,
                        "kdmax_weighted_denominator": kd_denominator,
                    }
                )

            tx_contrib = ann.reset_index()[
                [
                    "Transcript_ID_clean",
                    "gene_id",
                    "gene_name",
                    "biotype",
                    "KDmax",
                    "tpm",
                    "log2_TPM_plus_1",
                    col,
                ]
            ].copy()
            tx_contrib = tx_contrib.rename(columns={col: "probability_class"})
            tx_contrib["sample"] = sample
            tx_contrib["predictor"] = predictor
            tx_contrib["absolute_burden"] = tx_contrib["log2_TPM_plus_1"]
            tx_contrib["normalized_fraction"] = tx_contrib["log2_TPM_plus_1"] / denominator if denominator else np.nan
            tx_contrib["KDmax_weighted_fraction"] = (
                tx_contrib["KDmax"] * tx_contrib["log2_TPM_plus_1"] / kd_denominator if kd_denominator else np.nan
            )
            contribution_rows.append(tx_contrib)

    metrics = pd.DataFrame(rows).merge(sample_meta[["sample", "Group", "subject", "treatment"]], on="sample", how="left")
    contributions = pd.concat(contribution_rows, ignore_index=True).merge(
        sample_meta[["sample", "Group", "subject", "treatment"]], on="sample", how="left"
    )
    metrics.to_csv(OUT / "sample_level_probability_class_burden_metrics.tsv", sep="\t", index=False)
    contributions.to_csv(OUT / "transcript_level_probability_class_metric_contributions.tsv", sep="\t", index=False)
    return metrics, contributions


def group_statistics(df, group_cols, value_cols):
    rows = []
    grouped = [((), df)] if not group_cols else df.groupby(group_cols, dropna=False)
    for keys, sub in grouped:
        if not isinstance(keys, tuple):
            keys = (keys,)
        base = dict(zip(group_cols, keys))
        for value in value_cols:
            ac = sub.loc[sub["Group"].eq("AC"), value].dropna().astype(float)
            tol = sub.loc[sub["Group"].eq("TOL"), value].dropna().astype(float)
            p = np.nan
            if len(ac) > 0 and len(tol) > 0:
                p = mannwhitneyu(tol, ac, alternative="two-sided").pvalue
            rows.append(
                {
                    **base,
                    "metric": value,
                    "n_AC": len(ac),
                    "n_TOL": len(tol),
                    "median_AC": ac.median() if len(ac) else np.nan,
                    "median_TOL": tol.median() if len(tol) else np.nan,
                    "mean_AC": ac.mean() if len(ac) else np.nan,
                    "mean_TOL": tol.mean() if len(tol) else np.nan,
                    "delta_mean_TOL_minus_AC": (tol.mean() - ac.mean()) if len(ac) and len(tol) else np.nan,
                    "mannwhitney_p_value": p,
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["FDR_BH"] = multipletests(out["mannwhitney_p_value"].fillna(1.0), method="fdr_bh")[1]
    return out


def cliffs_delta(x, y):
    x = np.asarray(pd.Series(x).dropna().astype(float))
    y = np.asarray(pd.Series(y).dropna().astype(float))
    if len(x) == 0 or len(y) == 0:
        return np.nan
    gt = sum(np.sum(xi > y) for xi in x)
    lt = sum(np.sum(xi < y) for xi in x)
    return (gt - lt) / (len(x) * len(y))


def format_stat_number(value, digits=3):
    if pd.isna(value):
        return "NA"
    if abs(value) >= 100:
        return f"{value:.1f}"
    if abs(value) >= 1:
        return f"{value:.2f}"
    return f"{value:.3g}"


def metric_panel_stats(df):
    rows = []
    for metric in METRIC_LABELS:
        ac = df.loc[df["Group"].eq("AC"), metric].dropna().astype(float)
        tol = df.loc[df["Group"].eq("TOL"), metric].dropna().astype(float)
        if len(ac) and len(tol):
            p = mannwhitneyu(tol, ac, alternative="two-sided").pvalue
            delta = tol.mean() - ac.mean()
            cliff = cliffs_delta(tol, ac)
        else:
            p = delta = cliff = np.nan
        rows.append(
            {
                "metric": metric,
                "n_AC": len(ac),
                "n_TOL": len(tol),
                "mean_AC": ac.mean() if len(ac) else np.nan,
                "mean_TOL": tol.mean() if len(tol) else np.nan,
                "delta_mean_TOL_minus_AC": delta,
                "mannwhitney_p_value": p,
                "cliffs_delta_TOL_vs_AC": cliff,
            }
        )
    return pd.DataFrame(rows)


def draw_three_metric_ac_tol_panel(df, title, out_stem, y_prefix="Amyloidogenic"):
    stats = metric_panel_stats(df)
    titles = {
        "absolute_burden": "Burden",
        "normalized_fraction": "Fraction",
        "KDmax_weighted_fraction": "KDmax-weighted Fraction",
    }
    ylabels = {
        "absolute_burden": f"{y_prefix} burden\nΣ log2(TPM+1)",
        "normalized_fraction": f"{y_prefix} fraction\nΣ log2(TPM+1) / Σ log2(TPM+1) all protein-coding",
        "KDmax_weighted_fraction": (
            f"KDmax-weighted {y_prefix.lower()} fraction\n"
            "Σ(KDmax×log2(TPM+1)) / Σ(KDmax×log2(TPM+1)) all protein-coding"
        ),
    }
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), sharex=False)
    for ax, metric in zip(axes, METRIC_LABELS):
        d = df[["sample", "Group", metric]].dropna().copy()
        sns.boxplot(
            data=d,
            x="Group",
            y=metric,
            order=["AC", "TOL"],
            palette=GROUP_PALETTE,
            width=0.48,
            fliersize=0,
            ax=ax,
        )
        sns.stripplot(
            data=d,
            x="Group",
            y=metric,
            order=["AC", "TOL"],
            palette=GROUP_PALETTE,
            edgecolor="white",
            linewidth=0.7,
            size=5,
            alpha=0.95,
            ax=ax,
        )
        stat = stats[stats["metric"].eq(metric)].iloc[0]
        annotation = (
            f"Δ TOL-AC = {format_stat_number(stat['delta_mean_TOL_minus_AC'])}\n"
            f"MW p = {format_stat_number(stat['mannwhitney_p_value'])}\n"
            f"Cliff's δ = {format_stat_number(stat['cliffs_delta_TOL_vs_AC'])}"
        )
        ax.text(
            0.02,
            0.96,
            annotation,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.5,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#CBD5DF", "alpha": 0.95},
        )
        ax.set_title(titles[metric], fontsize=12, pad=8)
        ax.set_xlabel("")
        ax.set_ylabel(ylabels[metric], fontsize=9)
        ax.grid(False)
        sns.despine(ax=ax)
    fig.suptitle(title, fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / f"{out_stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{out_stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    stats.to_csv(OUT / f"{out_stem}_statistics.tsv", sep="\t", index=False)
    return stats


def draw_sample_metric_boxplots(metrics):
    metrics = metrics[metrics["probability_class"].isin(EXTREME_CLASS_ORDER)].copy()
    long = metrics.melt(
        id_vars=["sample", "Group", "subject", "predictor", "probability_class"],
        value_vars=list(METRIC_LABELS),
        var_name="metric",
        value_name="value",
    )
    long["probability_class"] = pd.Categorical(long["probability_class"], EXTREME_CLASS_ORDER, ordered=True)
    long["metric_label"] = long["metric"].map(METRIC_LABELS)

    for predictor, sub in long.groupby("predictor", sort=False):
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), sharex=False)
        for ax, metric in zip(axes, METRIC_LABELS):
            d = sub[sub["metric"].eq(metric)].copy()
            sns.boxplot(
                data=d,
                x="probability_class",
                y="value",
                hue="Group",
                hue_order=["AC", "TOL"],
                palette=GROUP_PALETTE,
                width=0.68,
                fliersize=0,
                ax=ax,
            )
            sns.stripplot(
                data=d,
                x="probability_class",
                y="value",
                hue="Group",
                hue_order=["AC", "TOL"],
                dodge=True,
                palette=GROUP_PALETTE,
                edgecolor="#222222",
                linewidth=0.45,
                size=4.5,
                alpha=0.9,
                ax=ax,
            )
            ax.set_title(METRIC_LABELS[metric], fontsize=11)
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.tick_params(axis="x", rotation=18, labelsize=9)
            ax.grid(axis="y", color="#E7E7E7", linewidth=0.8)
            if ax.legend_:
                ax.legend_.remove()
        handles, labels = axes[0].get_legend_handles_labels()
        handles, labels = handles[:2], labels[:2]
        fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
        fig.suptitle(f"Validation protein-coding transcripts: {PREDICTOR_LABELS[predictor]}", fontsize=14, y=0.98)
        fig.tight_layout(rect=[0, 0.08, 1, 0.94])
        safe = predictor.replace("/", "_").replace(" ", "_")
        fig.savefig(OUT / f"{safe}_probability_class_burden_boxplots.png", dpi=300)
        fig.savefig(OUT / f"{safe}_probability_class_burden_boxplots.pdf")
        plt.close(fig)

    g = sns.catplot(
        data=long,
        x="probability_class",
        y="value",
        hue="Group",
        col="metric_label",
        row="predictor",
        kind="box",
        hue_order=["AC", "TOL"],
        palette=GROUP_PALETTE,
        sharey=False,
        height=3.0,
        aspect=1.35,
        fliersize=0,
    )
    g.set_axis_labels("", "")
    g.set_titles(row_template="{row_name}", col_template="{col_name}")
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=18, labelsize=8)
        ax.grid(axis="y", color="#E7E7E7", linewidth=0.8)
    g.fig.suptitle("Validation: burden metrics by q17/q83 probability class", fontsize=15, y=1.01)
    g.fig.tight_layout()
    g.fig.savefig(OUT / "all_predictors_probability_class_burden_boxplots.png", dpi=300, bbox_inches="tight")
    g.fig.savefig(OUT / "all_predictors_probability_class_burden_boxplots.pdf", bbox_inches="tight")
    plt.close(g.fig)


def draw_extreme_class_ac_tol_panels(metrics):
    rows = []
    for (predictor, cls), sub in metrics[metrics["probability_class"].isin(EXTREME_CLASS_ORDER)].groupby(
        ["predictor", "probability_class"], sort=False
    ):
        safe_predictor = predictor.replace("/", "_").replace(" ", "_")
        safe_class = cls.replace("-", "_")
        title = f"Validation {PREDICTOR_LABELS[predictor]}: {cls}"
        prefix = cls.replace("-Amyloidogenic", "-amyloidogenic")
        stats = draw_three_metric_ac_tol_panel(
            sub,
            title=title,
            out_stem=f"{safe_predictor}_{safe_class}_AC_vs_TOL_three_metric_panel",
            y_prefix=prefix,
        )
        stats.insert(0, "probability_class", cls)
        stats.insert(0, "predictor", predictor)
        rows.append(stats)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(OUT / "extreme_probability_class_three_metric_panel_statistics.tsv", sep="\t", index=False)
    return out


def draw_target_transcript_boxplots(contributions):
    target = contributions[
        contributions["predictor"].eq("Consensus")
        & contributions["Transcript_ID_clean"].isin(TARGET_TRANSCRIPTS.values())
    ].copy()
    target["target"] = target["Transcript_ID_clean"].map({v: k for k, v in TARGET_TRANSCRIPTS.items()})
    target.to_csv(OUT / "target_transcript_consensus_metric_contributions.tsv", sep="\t", index=False)

    long = target.melt(
        id_vars=["sample", "Group", "subject", "target", "probability_class"],
        value_vars=list(METRIC_LABELS),
        var_name="metric",
        value_name="value",
    )
    long["metric_label"] = long["metric"].map(METRIC_LABELS)

    g = sns.catplot(
        data=long,
        x="Group",
        y="value",
        col="metric_label",
        row="target",
        order=["AC", "TOL"],
        kind="box",
        palette=GROUP_PALETTE,
        sharey=False,
        height=3.0,
        aspect=1.1,
        fliersize=0,
    )
    g.map_dataframe(
        sns.stripplot,
        x="Group",
        y="value",
        order=["AC", "TOL"],
        palette=GROUP_PALETTE,
        edgecolor="#222222",
        linewidth=0.45,
        size=4.5,
        alpha=0.9,
    )
    for ax in g.axes.flat:
        ax.grid(axis="y", color="#E7E7E7", linewidth=0.8)
    g.set_axis_labels("", "")
    g.set_titles(row_template="{row_name}", col_template="{col_name}")
    g.fig.suptitle("Validation target transcripts: consensus metric contributions", fontsize=15, y=1.02)
    g.fig.tight_layout()
    g.fig.savefig(OUT / "target_transcripts_consensus_metric_boxplots.png", dpi=300, bbox_inches="tight")
    g.fig.savefig(OUT / "target_transcripts_consensus_metric_boxplots.pdf", bbox_inches="tight")
    plt.close(g.fig)
    return target


def draw_target_gene_boxplots(contributions):
    target = contributions[
        contributions["predictor"].eq("Consensus")
        & contributions["gene_name"].isin(TARGET_GENES)
    ].copy()
    gene_level = (
        target.groupby(["sample", "Group", "subject", "gene_name"], as_index=False)
        .agg(
            absolute_burden=("absolute_burden", "sum"),
            normalized_fraction=("normalized_fraction", "sum"),
            KDmax_weighted_fraction=("KDmax_weighted_fraction", "sum"),
            n_protein_coding_transcripts=("Transcript_ID_clean", "nunique"),
        )
    )
    gene_level.to_csv(OUT / "target_gene_consensus_metric_contributions.tsv", sep="\t", index=False)

    stats = group_statistics(gene_level, ["gene_name"], list(METRIC_LABELS))
    stats.to_csv(OUT / "target_gene_consensus_group_statistics.tsv", sep="\t", index=False)

    panel_stats = []
    for gene, sub in gene_level.groupby("gene_name", sort=True):
        s = draw_three_metric_ac_tol_panel(
            sub,
            title=f"Validation gene-level consensus metrics: {gene}",
            out_stem=f"target_gene_{gene}_consensus_AC_vs_TOL_three_metric_panel",
            y_prefix=f"{gene}",
        )
        s.insert(0, "gene_name", gene)
        panel_stats.append(s)
    pd.concat(panel_stats, ignore_index=True).to_csv(
        OUT / "target_gene_three_metric_panel_statistics.tsv", sep="\t", index=False
    )

    long = gene_level.melt(
        id_vars=["sample", "Group", "subject", "gene_name", "n_protein_coding_transcripts"],
        value_vars=list(METRIC_LABELS),
        var_name="metric",
        value_name="value",
    )
    long["metric_label"] = long["metric"].map(
        {
            "absolute_burden": "Burden",
            "normalized_fraction": "Fraction",
            "KDmax_weighted_fraction": "KDmax-weighted Fraction",
        }
    )
    g = sns.catplot(
        data=long,
        x="Group",
        y="value",
        col="metric_label",
        row="gene_name",
        order=["AC", "TOL"],
        kind="box",
        palette=GROUP_PALETTE,
        sharey=False,
        height=3.1,
        aspect=1.1,
        fliersize=0,
    )
    g.map_dataframe(
        sns.stripplot,
        x="Group",
        y="value",
        order=["AC", "TOL"],
        palette=GROUP_PALETTE,
        edgecolor="white",
        linewidth=0.7,
        size=4.8,
        alpha=0.95,
    )
    for ax in g.axes.flat:
        ax.grid(False)
        sns.despine(ax=ax)
    g.set_axis_labels("", "")
    g.set_titles(row_template="{row_name}", col_template="{col_name}")
    g.fig.suptitle("Validation target genes: summed gene-level consensus metrics", fontsize=15, y=1.02)
    g.fig.tight_layout()
    g.fig.savefig(OUT / "target_genes_gene_level_consensus_metric_boxplots.png", dpi=300, bbox_inches="tight")
    g.fig.savefig(OUT / "target_genes_gene_level_consensus_metric_boxplots.pdf", bbox_inches="tight")
    plt.close(g.fig)
    return gene_level, stats


def draw_signature_boxplots(sample_meta, tpm, gencode_pc):
    from run_seven_signature_panel import SIGNATURES as SIGNATURE_GENE_SETS

    tx_ann = gencode_pc[["Transcript_ID_clean", "gene_name"]].drop_duplicates("Transcript_ID_clean")
    expr = tpm.loc[tpm.index.intersection(tx_ann["Transcript_ID_clean"])].copy()
    expr = expr.merge(tx_ann, left_index=True, right_on="Transcript_ID_clean", how="left")
    gene_tpm = expr.dropna(subset=["gene_name"]).groupby("gene_name")[sample_meta["sample"].tolist()].sum()
    log_gene = np.log2(gene_tpm + 1)
    gene_mean = log_gene.mean(axis=1)
    gene_sd = log_gene.std(axis=1, ddof=0).replace(0, np.nan)
    z_gene = log_gene.sub(gene_mean, axis=0).div(gene_sd, axis=0)

    scores = sample_meta[["sample", "Group", "subject"]].copy()
    coverage_rows = []
    for sig in SIGNATURES:
        genes = list(dict.fromkeys(SIGNATURE_GENE_SETS[sig]))
        present = [g for g in genes if g in z_gene.index]
        missing = [g for g in genes if g not in z_gene.index]
        scores[sig] = z_gene.loc[present].mean(axis=0).reindex(scores["sample"]).values if present else np.nan
        coverage_rows.append(
            {
                "signature": sig,
                "signature_label": SIGNATURES[sig],
                "signature_size": len(genes),
                "genes_detected": len(present),
                "genes_missing": len(missing),
                "detected_gene_symbols": ",".join(present),
                "missing_gene_symbols": ",".join(missing),
            }
        )
    scores.to_csv(OUT / "validation_signature_scores_selected.tsv", sep="\t", index=False)
    pd.DataFrame(coverage_rows).to_csv(OUT / "validation_signature_gene_coverage.tsv", sep="\t", index=False)

    stats = group_statistics(scores, [], list(SIGNATURES))
    stats.to_csv(OUT / "validation_signature_group_statistics.tsv", sep="\t", index=False)

    long = scores.melt(
        id_vars=["sample", "Group", "subject"],
        value_vars=list(SIGNATURES),
        var_name="signature",
        value_name="score",
    )
    long["signature_label"] = long["signature"].map(SIGNATURES)
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    sns.boxplot(
        data=long,
        x="signature_label",
        y="score",
        hue="Group",
        hue_order=["AC", "TOL"],
        palette=GROUP_PALETTE,
        width=0.68,
        fliersize=0,
        ax=ax,
    )
    sns.stripplot(
        data=long,
        x="signature_label",
        y="score",
        hue="Group",
        hue_order=["AC", "TOL"],
        dodge=True,
        palette=GROUP_PALETTE,
        edgecolor="#222222",
        linewidth=0.45,
        size=4.8,
        alpha=0.9,
        ax=ax,
    )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], loc="best", frameon=False)
    ax.set_xlabel("")
    ax.set_ylabel("Signature z-score")
    ax.set_title("Validation transcriptional signatures: TOL vs AC", fontsize=14)
    ax.tick_params(axis="x", rotation=20)
    ax.grid(axis="y", color="#E7E7E7", linewidth=0.8)
    fig.tight_layout()
    fig.savefig(OUT / "validation_signature_boxplots_TOL_vs_AC.png", dpi=300)
    fig.savefig(OUT / "validation_signature_boxplots_TOL_vs_AC.pdf")
    plt.close(fig)
    return scores, stats


def write_report(sample_meta, gencode_path, metrics_stats, target_stats, signature_stats):
    def table(df):
        if df.empty:
            return "_No rows._"
        d = df.copy()
        for col in d.columns:
            if pd.api.types.is_float_dtype(d[col]):
                d[col] = d[col].map(lambda x: "" if pd.isna(x) else f"{x:.5g}")
        return d.to_markdown(index=False)

    report = [
        "# Validation q17/q83 probability-group analysis",
        "",
        "## Cohort",
        "",
        table(sample_meta[["sample", "subject", "Group", "raw_Group", "time", "treatment"]]),
        "",
        "## Rules",
        "",
        "- AMYPred-FRL: prob <= 0.1206 Non-Amyloidogenic; prob >= 0.8802 Amyloidogenic; otherwise Intermediate.",
        "- AmyloGramPy: prob <= 0.6947 Non-Amyloidogenic; prob >= 0.8898 Amyloidogenic; otherwise Intermediate.",
        "- Consensus: both above high thresholds Amyloidogenic; both below low thresholds Non-Amyloidogenic; otherwise Intermediate.",
        f"- GENCODE filter: only `biotype == protein_coding`; annotation `{gencode_path}`.",
        "",
        "## Main Outputs",
        "",
        "- `all_predictors_probability_class_burden_boxplots.png/pdf`",
        "- `target_genes_gene_level_consensus_metric_boxplots.png/pdf`",
        "- `target_gene_<GENE>_consensus_AC_vs_TOL_three_metric_panel.png/pdf`",
        "- `validation_signature_boxplots_TOL_vs_AC.png/pdf`",
        "",
        "## AC vs TOL Statistics: Sample-Level Burden Metrics",
        "",
        table(metrics_stats),
        "",
        "## AC vs TOL Statistics: Target Gene-Level Contributions",
        "",
        table(target_stats),
        "",
        "## AC vs TOL Statistics: Transcriptional Signatures",
        "",
        table(signature_stats),
    ]
    (OUT / "VALIDATION_Q17Q83_PROBABILITY_GROUP_REPORT.md").write_text("\n".join(report) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    sample_meta = load_sample_metadata()
    gencode_pc, gencode_path = load_gencode_pc()
    tpm = load_expression_matrix(root, sample_meta["sample"].tolist())
    collapsed = build_collapsed_probability_classes(root, sample_meta["sample"].tolist(), gencode_pc)
    metrics, contributions = compute_sample_class_metrics(
        tpm=tpm,
        collapsed=collapsed,
        sample_meta=sample_meta,
        protein_coding_transcripts=set(gencode_pc["Transcript_ID_clean"]),
    )

    metrics_stats_all = group_statistics(metrics, ["predictor", "probability_class"], list(METRIC_LABELS))
    metrics_stats_all.to_csv(OUT / "sample_level_probability_class_group_statistics_all_classes.tsv", sep="\t", index=False)
    metrics_stats = metrics_stats_all[metrics_stats_all["probability_class"].isin(EXTREME_CLASS_ORDER)].copy()
    metrics_stats.to_csv(OUT / "sample_level_probability_class_group_statistics.tsv", sep="\t", index=False)
    draw_sample_metric_boxplots(metrics)
    draw_extreme_class_ac_tol_panels(metrics)

    draw_target_transcript_boxplots(contributions)
    _, target_stats = draw_target_gene_boxplots(contributions)

    _, signature_stats = draw_signature_boxplots(sample_meta, tpm, gencode_pc)

    summary = {
        "zarr_path": str(ZARR_PATH),
        "output_dir": str(OUT),
        "validation_samples": VALIDATION_SAMPLES,
        "thresholds": THRESHOLDS,
        "n_protein_coding_collapsed_rows": int(collapsed.shape[0]),
        "n_samples": int(sample_meta.shape[0]),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_report(sample_meta, gencode_path, metrics_stats, target_stats, signature_stats)


if __name__ == "__main__":
    main()
