from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
import zarr
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "protein_coding_q17q83_consensus_feature_plots"
ZARR_PATH = ROOT / "project.zarr"
SAMPLE_META_PATH = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "sample_metadata_BC.tsv"
GENCODE_PATH = (
    ROOT
    / "analysis_outputs"
    / "amyloid_bc_control_vs_chili"
    / "rank_based_amyloid_burden_signature"
    / "gencode_v48_transcript_annotation.tsv"
)

FEATURES = [
    "protein_length",
    "beta_propensity",
    "alpha_propensity",
    "chameleon_score",
    "KD_max",
    "KD_mean",
    "KD_sd",
    "charge",
    "aromaticity",
    "pI",
]

CLASS_COL = "BothPredictors_Q17_Q83_Class"
ORDER = ["Amyloidogenic", "Non-Amyloidogenic"]
ALL_CLASSES = ["Amyloidogenic", "Intermediate", "Non-Amyloidogenic"]
PALETTE = {"Amyloidogenic": "#D9824B", "Non-Amyloidogenic": "#2C7A7B", "Intermediate": "#9AA3AD"}
LOG_FEATURES = {"protein_length", "charge"}


def fmt_int(x):
    return f"{int(x):,}".replace(",", " ")


def zarr_df(group):
    return pd.DataFrame({col: group[col][:] for col in group.attrs["columns"]})


def transcript_id_clean(value):
    match = re.search(r"(ENST\d+(?:\.\d+)?)", str(value))
    if not match:
        return pd.NA
    return match.group(1).split(".")[0]


def cliffs_delta(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan
    stat, _ = mannwhitneyu(x, y, alternative="two-sided")
    return float(2 * stat / (len(x) * len(y)) - 1)


def bh_fdr(pvalues):
    p = np.asarray(pvalues, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return q
    idx = np.where(valid)[0]
    order = idx[np.argsort(p[valid])]
    ranked = p[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    q[order] = np.clip(adjusted, 0, 1)
    return q


def load_protein_coding_consensus():
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    frames = []
    for sample in root["samples"].group_keys():
        pf = zarr_df(root[f"samples/{sample}/protein_features"])
        amy = zarr_df(root[f"samples/{sample}/amyloid"])
        pf["ID"] = pf["sample"].astype(str) + "|" + pf["protein_id"].astype(str)
        amy["ID"] = amy["sample"].astype(str) + "|" + amy["Sequence_ID"].astype(str)
        merged = pf.merge(amy.drop(columns=["sample"], errors="ignore"), on="ID", how="inner")
        frames.append(merged)

    raw = pd.concat(frames, ignore_index=True)
    raw[CLASS_COL] = raw[CLASS_COL].astype("string").replace("", pd.NA)
    raw["Transcript_ID_clean"] = raw["Sequence_ID"].map(transcript_id_clean)

    ann = pd.read_csv(GENCODE_PATH, sep="\t")
    ann = ann[["Transcript_ID_clean", "gene_id", "gene_name", "biotype"]].drop_duplicates("Transcript_ID_clean")
    raw = raw.merge(ann, on="Transcript_ID_clean", how="left")
    raw = raw[raw["biotype"].eq("protein_coding")].copy()

    sample_meta = pd.read_csv(SAMPLE_META_PATH, sep="\t")
    zarr_samples = set(raw["sample"].astype(str).unique())
    sample_meta = sample_meta[sample_meta["sample"].astype(str).isin(zarr_samples)].copy()
    return raw, sample_meta, root.attrs.get("q17_q83_probability_thresholds", {})


def build_stats(df):
    rows = []
    for feature in FEATURES:
        plot_df = df[[feature, CLASS_COL]].dropna()
        amy = plot_df.loc[plot_df[CLASS_COL].eq("Amyloidogenic"), feature].astype(float).to_numpy()
        non = plot_df.loc[plot_df[CLASS_COL].eq("Non-Amyloidogenic"), feature].astype(float).to_numpy()
        if len(amy) and len(non):
            stat, p = mannwhitneyu(amy, non, alternative="two-sided")
            delta = cliffs_delta(amy, non)
        else:
            stat, p, delta = np.nan, np.nan, np.nan
        rows.append(
            {
                "feature": feature,
                "n_amyloidogenic": int(len(amy)),
                "n_non_amyloidogenic": int(len(non)),
                "median_amyloidogenic": float(np.median(amy)) if len(amy) else np.nan,
                "median_non_amyloidogenic": float(np.median(non)) if len(non) else np.nan,
                "mannwhitney_u": float(stat) if np.isfinite(stat) else np.nan,
                "p_value": float(p) if np.isfinite(p) else np.nan,
                "cliffs_delta": float(delta) if np.isfinite(delta) else np.nan,
            }
        )
    stats = pd.DataFrame(rows)
    stats["q_value"] = bh_fdr(stats["p_value"].to_numpy())
    return stats


def plot_value(series, feature):
    values = pd.to_numeric(series, errors="coerce")
    if feature == "charge":
        return np.sign(values) * np.log10(np.abs(values) + 1)
    if feature == "protein_length":
        return np.log10(values.clip(lower=1))
    return values


def plot_ylabel(feature):
    if feature == "charge":
        return "charge, signed log10(|x|+1)"
    if feature == "protein_length":
        return "protein_length, log10"
    return feature


def downsample_points(plot_df, feature, max_per_group=900):
    pieces = []
    for cls, sub in plot_df.groupby(CLASS_COL, observed=True):
        if len(sub) > max_per_group:
            pieces.append(sub.sample(max_per_group, random_state=17))
        else:
            pieces.append(sub)
    return pd.concat(pieces, ignore_index=True) if pieces else plot_df.iloc[0:0].copy()


def draw_consensus_distribution_panel(df, stats, sample_summary, threshold_text):
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.22)
    plot_source = df[df[CLASS_COL].isin(ORDER)].copy()
    plot_source[CLASS_COL] = pd.Categorical(plot_source[CLASS_COL], categories=ORDER, ordered=True)

    fig, axes = plt.subplots(3, 4, figsize=(12.6, 12.6))
    axes = axes.flatten()

    for i, feature in enumerate(FEATURES):
        ax = axes[i]
        plot_df = plot_source[[feature, CLASS_COL]].dropna().copy()
        plot_df["plot_value"] = plot_value(plot_df[feature], feature)
        plot_df = plot_df.replace([np.inf, -np.inf], np.nan).dropna(subset=["plot_value"])

        sns.boxplot(
            data=plot_df,
            x=CLASS_COL,
            y="plot_value",
            order=ORDER,
            hue=CLASS_COL,
            hue_order=ORDER,
            palette=PALETTE,
            showfliers=False,
            width=0.84,
            linewidth=1.15,
            boxprops={"alpha": 0.72, "edgecolor": "#263238"},
            medianprops={"color": "#263238", "linewidth": 1.35},
            whiskerprops={"color": "#263238", "linewidth": 1.0},
            capprops={"color": "#263238", "linewidth": 1.0},
            ax=ax,
            legend=False,
        )
        point_df = downsample_points(plot_df, feature)
        sns.stripplot(
            data=point_df,
            x=CLASS_COL,
            y="plot_value",
            order=ORDER,
            ax=ax,
            color="#16324F",
            alpha=0.12,
            size=1.6,
            jitter=0.26,
        )

        row = stats.loc[stats["feature"].eq(feature)].iloc[0]
        ax.set_title(feature, fontsize=13, weight="bold", pad=12)
        ax.text(
            0.02,
            0.97,
            f"q={row['q_value']:.2e}\ndelta={row['cliffs_delta']:.2f}",
            transform=ax.transAxes,
            fontsize=9.5,
            color="#263238",
            va="top",
            ha="left",
            bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "none", "alpha": 0.78},
        )
        ax.set_xlabel("")
        ax.set_ylabel(plot_ylabel(feature), fontsize=10.5)
        ax.tick_params(axis="x", rotation=0, labelsize=10)
        ax.tick_params(axis="y", labelsize=10)
        for label in ax.get_xticklabels():
            label.set_horizontalalignment("center")
        ax.grid(axis="y", color="#D8E1EA", linewidth=0.8)
        ax.spines[["top", "right", "left"]].set_visible(False)

    counts = plot_source[CLASS_COL].value_counts().reindex(ORDER).fillna(0).astype(int)
    for j in range(len(FEATURES), len(axes)):
        axes[j].axis("off")

    fig.tight_layout(rect=(0, 0, 1, 1), w_pad=0.9, h_pad=1.05)
    out_prefix = OUT / "protein_coding_q17q83_consensus"
    fig.savefig(f"{out_prefix}_boxplots.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_prefix}_boxplots.pdf", bbox_inches="tight")
    plt.close(fig)


def repel_texts(ax, texts, xpad=0.012, ypad_frac=0.035, iterations=140):
    fig = ax.figure
    for _ in range(iterations):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        moved = False
        boxes = [text.get_window_extent(renderer=renderer).expanded(1.05, 1.15) for text in texts]
        y0, y1 = ax.get_ylim()
        dy = (y1 - y0) * ypad_frac
        x0, x1 = ax.get_xlim()
        dx = (x1 - x0) * xpad
        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                if boxes[i].overlaps(boxes[j]):
                    xi, yi = texts[i].get_position()
                    xj, yj = texts[j].get_position()
                    if yi <= yj:
                        texts[i].set_position((xi - dx * 0.35, yi - dy))
                        texts[j].set_position((xj + dx * 0.35, yj + dy))
                    else:
                        texts[i].set_position((xi + dx * 0.35, yi + dy))
                        texts[j].set_position((xj - dx * 0.35, yj - dy))
                    moved = True
        if not moved:
            break


def draw_consensus_volcano(stats, sample_summary):
    plot_df = stats.copy()
    plot_df["q_value_adj"] = plot_df["q_value"].replace(0, 1e-300)
    plot_df["neg_log10_q"] = -np.log10(plot_df["q_value_adj"])
    effect_threshold = 0.1
    sig_threshold = 0.05

    conditions = [
        (plot_df["q_value"] < sig_threshold) & (plot_df["cliffs_delta"] >= effect_threshold),
        (plot_df["q_value"] < sig_threshold) & (plot_df["cliffs_delta"] <= -effect_threshold),
    ]
    choices = ["Significant Up", "Significant Down"]
    plot_df["status"] = np.select(conditions, choices, default="Not significant")
    colors = {"Significant Up": "#D95F02", "Significant Down": "#1F78B4", "Not significant": "#A9B2BD"}

    fig, ax = plt.subplots(figsize=(9.2, 6.8))
    for status in ["Not significant", "Significant Down", "Significant Up"]:
        sub = plot_df[plot_df["status"].eq(status)]
        ax.scatter(
            sub["cliffs_delta"],
            sub["neg_log10_q"],
            s=98 if status != "Not significant" else 74,
            color=colors[status],
            edgecolor="white",
            linewidth=0.9,
            label=status,
            zorder=3 if status != "Not significant" else 2,
        )

    for x, label, ha in [
        (-effect_threshold, "Cliff's Delta -0.1", "right"),
        (effect_threshold, "Cliff's Delta +0.1", "left"),
    ]:
        ax.axvline(x, color="#596579", linestyle=":", linewidth=1.3)
        ax.text(x, 0.985, label, transform=ax.get_xaxis_transform(), fontsize=9.5, color="#596579", ha=ha, va="top")
    ax.axvline(0, color="#263238", linestyle="--", linewidth=1.0)

    for q, label, xpos in [(0.05, "q=0.05", 0.985), (0.01, "q=0.01", 0.90)]:
        y = -np.log10(q)
        ax.axhline(y, color="#7C8794", linestyle="--", linewidth=0.9, alpha=0.85)
        ax.text(xpos, y + 0.8, label, transform=ax.get_yaxis_transform(), fontsize=9.5, color="#596579", ha="right")

    texts = []
    for _, row in plot_df.iterrows():
        x = row["cliffs_delta"]
        y = row["neg_log10_q"]
        xoff = 0.018 if x >= 0 else -0.018
        yoff = 0.12
        text = ax.annotate(
            row["feature"],
            xy=(x, y),
            xytext=(x + xoff, y + yoff),
            fontsize=11,
            color="#1B3557",
            arrowprops={"arrowstyle": "-", "color": "#7C8794", "lw": 0.7, "shrinkA": 2, "shrinkB": 4},
        )
        texts.append(text)

    ax.set_xlabel("Cliff's Delta (Amyloidogenic vs Non-Amyloidogenic)", fontsize=12)
    ax.set_ylabel("-log10(q-value)", fontsize=12)
    ax.set_title("")
    ax.tick_params(axis="both", labelsize=11)
    ax.legend(frameon=False, loc="upper left", title="Statistical class", fontsize=10, title_fontsize=11)
    ax.grid(color="#C9D4DF", linewidth=0.8, alpha=0.72)
    ax.spines[["top", "right"]].set_visible(False)
    ax.margins(x=0.18, y=0.18)
    ax.set_ylim(bottom=0)
    repel_texts(ax, texts)
    fig.tight_layout()

    out_prefix = OUT / "protein_coding_q17q83_consensus"
    fig.savefig(f"{out_prefix}_volcano.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_prefix}_volcano.pdf", bbox_inches="tight")
    plt.close(fig)


def write_summary_tables(df, stats, thresholds):
    OUT.mkdir(parents=True, exist_ok=True)
    counts = (
        df[CLASS_COL]
        .astype("string")
        .fillna("Missing")
        .value_counts()
        .reindex([*ALL_CLASSES, "Missing"])
        .fillna(0)
        .astype(int)
        .rename_axis("class")
        .reset_index(name="count")
    )
    counts["filter"] = "GENCODE protein_coding"
    counts["classifier"] = CLASS_COL
    counts.to_csv(OUT / "protein_coding_q17q83_consensus_class_counts.tsv", sep="\t", index=False)
    stats.to_csv(OUT / "protein_coding_q17q83_consensus_feature_statistics.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"predictor": col, "q17": vals["q17"], "q83": vals["q83"]}
            for col, vals in thresholds.items()
        ]
    ).to_csv(OUT / "protein_coding_q17q83_thresholds_used.tsv", sep="\t", index=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df, sample_meta, thresholds = load_protein_coding_consensus()
    n_samples = int(sample_meta["sample"].nunique())
    n_ac = int(sample_meta.loc[sample_meta["Group"].eq("AC"), "sample"].nunique())
    n_tol = int(sample_meta.loc[sample_meta["Group"].eq("TOL"), "sample"].nunique())
    n_subjects = int(sample_meta["subject"].nunique()) if "subject" in sample_meta.columns else n_samples
    sample_summary = f"n={n_samples} BC study samples ({n_ac} AC, {n_tol} TOL), {n_subjects} subjects; protein-coding rows={fmt_int(len(df))}"

    thresholds = dict(thresholds)
    threshold_text = (
        f"AMYPred q17/q83: {thresholds['AMYPred_Prob']['q17']:.3f}/{thresholds['AMYPred_Prob']['q83']:.3f}\n"
        f"AmyloGramPy q17/q83: {thresholds['AmyloGramPy_Prob']['q17']:.3f}/{thresholds['AmyloGramPy_Prob']['q83']:.3f}"
    )

    compare_df = df[df[CLASS_COL].isin(ORDER)].copy()
    stats = build_stats(compare_df)
    write_summary_tables(df, stats, thresholds)
    draw_consensus_distribution_panel(compare_df, stats, sample_summary, threshold_text)
    draw_consensus_volcano(stats, sample_summary)

    for path in sorted(OUT.glob("protein_coding_q17q83_consensus*")):
        print(path)


if __name__ == "__main__":
    main()
