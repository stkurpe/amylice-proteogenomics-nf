from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import zarr
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "study_cohort_predictor_feature_plots"
ZARR_PATH = ROOT / "project.zarr"
SAMPLE_META_PATH = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "sample_metadata_BC.tsv"

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

PREDICTORS = ["AMYPred_Pred", "AmyloGramPy_Pred"]
ORDER = ["Amyloid", "Non-Amyloid"]
PALETTE = {"Amyloid": "#6F6BAA", "Non-Amyloid": "#596579"}
GENCODE_PATH = (
    ROOT
    / "analysis_outputs"
    / "amyloid_bc_control_vs_chili"
    / "rank_based_amyloid_burden_signature"
    / "gencode_v48_transcript_annotation.tsv"
)


def fmt_int(x):
    return f"{int(x):,}".replace(",", " ")


def cliffs_delta(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan
    gt = sum(np.sum(xi > y) for xi in x)
    lt = sum(np.sum(xi < y) for xi in x)
    return float((gt - lt) / (len(x) * len(y)))


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


def zarr_df(group):
    cols = list(group.attrs["columns"])
    return pd.DataFrame({col: group[col][:] for col in cols})


def clean_prediction_column(s):
    return s.astype("string").replace("", pd.NA)


def transcript_id_clean(s):
    match = re.search(r"(ENST\d+(?:\.\d+)?)", str(s))
    if not match:
        return pd.NA
    return match.group(1).split(".")[0]


def add_canonical_columns(raw):
    raw = raw.copy()
    if "AmyloGramPy_Pred" in raw.columns:
        amylogram = clean_prediction_column(raw["AmyloGramPy_Pred"])
    else:
        amylogram = pd.Series(pd.NA, index=raw.index, dtype="string")
    if "AmyloGram_Pred" in raw.columns:
        amylogram = amylogram.fillna(clean_prediction_column(raw["AmyloGram_Pred"]))
    raw["AmyloGram_canonical_Pred"] = amylogram

    amypred = clean_prediction_column(raw["AMYPred_Pred"])
    consensus = pd.Series("Partial", index=raw.index, dtype="string")
    both_present = amypred.notna() & amylogram.notna()
    consensus.loc[both_present & amypred.eq(amylogram) & amypred.eq("Amyloid")] = "Amyloid"
    consensus.loc[both_present & amypred.eq(amylogram) & amypred.eq("Non-Amyloid")] = "Non-Amyloid"
    consensus.loc[both_present & ~amypred.eq(amylogram)] = "Discordant"
    raw["Consensus"] = consensus
    return raw


def load_raw_from_zarr(zarr_path):
    root = zarr.open_group(str(zarr_path), mode="r")
    frames = []
    for sample in root["samples"].group_keys():
        pf = zarr_df(root[f"samples/{sample}/protein_features"])
        amy = zarr_df(root[f"samples/{sample}/amyloid"])
        pf["ID"] = pf["sample"].astype(str) + "|" + pf["protein_id"].astype(str)
        amy["ID"] = amy["sample"].astype(str) + "|" + amy["Sequence_ID"].astype(str)
        merged = pf.merge(amy.drop(columns=["sample"], errors="ignore"), on="ID", how="inner")
        frames.append(merged)
    raw = pd.concat(frames, ignore_index=True)
    for predictor in PREDICTORS:
        raw[predictor] = clean_prediction_column(raw[predictor])
    raw["Consensus"] = clean_prediction_column(raw["Consensus"])
    raw = add_canonical_columns(raw)
    return raw


def write_prediction_count_tables(raw):
    predictors = ["AMYPred_Pred", "AmyloGramPy_Pred", "AmyloGram_canonical_Pred", "Consensus"]
    rows = []
    for predictor in predictors:
        s = clean_prediction_column(raw[predictor])
        counts = s.value_counts(dropna=False)
        present = int(s.notna().sum())
        rows.append(
            {
                "predictor": predictor,
                "n_rows_total": int(len(raw)),
                "n_present": present,
                "n_missing": int(s.isna().sum()),
                "n_amyloid": int(counts.get("Amyloid", 0)),
                "n_non_amyloid": int(counts.get("Non-Amyloid", 0)),
                "n_discordant": int(counts.get("Discordant", 0)),
                "n_partial": int(counts.get("Partial", 0)),
                "amyloid_fraction_of_present": float(counts.get("Amyloid", 0) / present) if present else np.nan,
            }
        )
    pd.DataFrame(rows).to_csv(OUT / "study_cohort_amyloidogenic_counts_by_predictor.tsv", sep="\t", index=False)

    amypred_amyloid = clean_prediction_column(raw["AMYPred_Pred"]).eq("Amyloid")
    amylogram_py_amyloid = clean_prediction_column(raw["AmyloGramPy_Pred"]).eq("Amyloid")
    amylogram_canonical_amyloid = clean_prediction_column(raw["AmyloGram_canonical_Pred"]).eq("Amyloid")
    overlap_rows = [
        ("AMYPred amyloid", int(amypred_amyloid.sum())),
        ("AmyloGramPy amyloid", int(amylogram_py_amyloid.sum())),
        ("AmyloGram canonical amyloid", int(amylogram_canonical_amyloid.sum())),
        ("AMYPred & AmyloGramPy amyloid", int((amypred_amyloid & amylogram_py_amyloid).sum())),
        ("AMYPred & AmyloGram canonical amyloid", int((amypred_amyloid & amylogram_canonical_amyloid).sum())),
    ]
    pd.DataFrame(overlap_rows, columns=["set", "count"]).to_csv(
        OUT / "study_cohort_amyloidogenic_predictor_overlaps.tsv", sep="\t", index=False
    )

    if GENCODE_PATH.exists():
        ann = pd.read_csv(GENCODE_PATH, sep="\t")
        ann = ann[["Transcript_ID_clean", "biotype"]].drop_duplicates("Transcript_ID_clean")
        biotype_df = raw[["Sequence_ID"]].copy()
        biotype_df["Transcript_ID_clean"] = biotype_df["Sequence_ID"].map(transcript_id_clean)
        biotype_df = biotype_df.merge(ann, on="Transcript_ID_clean", how="left")
        biotype_df["biotype"] = biotype_df["biotype"].fillna("unannotated")
        biotype_df["biotype"].value_counts().rename_axis("biotype").reset_index(name="count").to_csv(
            OUT / "study_cohort_prediction_rows_by_biotype.tsv", sep="\t", index=False
        )

        long_rows = []
        for predictor in predictors:
            pred = clean_prediction_column(raw[predictor])
            tmp = biotype_df[["biotype"]].copy()
            tmp["predictor"] = predictor
            tmp["prediction"] = pred.fillna("Missing").astype(str)
            long_rows.append(tmp)
        (
            pd.concat(long_rows, ignore_index=True)
            .groupby(["predictor", "biotype", "prediction"], dropna=False)
            .size()
            .reset_index(name="count")
            .to_csv(OUT / "study_cohort_prediction_counts_by_biotype_and_predictor.tsv", sep="\t", index=False)
        )


def build_stats(df, predictor):
    rows = []
    for feature in FEATURES:
        plot_df = df[[feature, predictor]].dropna()
        amy = plot_df.loc[plot_df[predictor].eq("Amyloid"), feature].astype(float).to_numpy()
        non = plot_df.loc[plot_df[predictor].eq("Non-Amyloid"), feature].astype(float).to_numpy()
        if len(amy) and len(non):
            stat, p = mannwhitneyu(amy, non, alternative="two-sided")
            delta = cliffs_delta(amy, non)
        else:
            stat, p, delta = np.nan, np.nan, np.nan
        rows.append(
            {
                "predictor": predictor,
                "feature": feature,
                "n_amyloid": int(len(amy)),
                "n_non_amyloid": int(len(non)),
                "median_amyloid": float(np.median(amy)) if len(amy) else np.nan,
                "median_non_amyloid": float(np.median(non)) if len(non) else np.nan,
                "mannwhitney_u": float(stat) if np.isfinite(stat) else np.nan,
                "p_value": float(p) if np.isfinite(p) else np.nan,
                "cliffs_delta": float(delta) if np.isfinite(delta) else np.nan,
            }
        )
    res = pd.DataFrame(rows)
    res["q_value"] = bh_fdr(res["p_value"].to_numpy())
    return res


def filter_for_predictor(raw, predictor):
    df = raw.dropna(subset=[predictor]).copy()
    df = df[df["Consensus"].isin(ORDER)].copy()
    df[predictor] = pd.Categorical(df[predictor], categories=ORDER, ordered=True)
    return df


def draw_boxplots(df, stats, predictor, sample_summary, out_prefix):
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.05)
    ncols = 3
    nrows = int(np.ceil(len(FEATURES) / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 4.2 * nrows))
    axes = axes.flatten()

    for i, feature in enumerate(FEATURES):
        ax = axes[i]
        plot_df = df[[feature, predictor]].dropna()
        sns.boxplot(
            data=plot_df,
            x=predictor,
            y=feature,
            hue=predictor,
            order=ORDER,
            hue_order=ORDER,
            ax=ax,
            showfliers=False,
            palette=PALETTE,
            width=0.58,
            linewidth=1.1,
            legend=False,
        )
        sns.stripplot(
            data=plot_df,
            x=predictor,
            y=feature,
            order=ORDER,
            ax=ax,
            color="#1B3557",
            alpha=0.16,
            size=2,
            jitter=0.22,
        )
        row = stats.loc[stats["feature"].eq(feature)].iloc[0]
        ax.set_title(f"{feature}\np={row['p_value']:.2e}, q={row['q_value']:.2e}", fontsize=11, weight="bold")
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=18)
        ax.grid(axis="y", color="#D7DFEA", linewidth=0.8)
        ax.spines[["top", "right", "left"]].set_visible(False)

    for j in range(len(FEATURES), len(axes)):
        axes[j].set_visible(False)

    counts = df[predictor].value_counts().reindex(ORDER).fillna(0).astype(int)
    fig.suptitle(
        f"Study cohort feature distributions by {predictor}\n"
        f"{sample_summary}; protein rows: Amyloid={fmt_int(counts['Amyloid'])}, Non-Amyloid={fmt_int(counts['Non-Amyloid'])}",
        fontsize=18,
        weight="bold",
        color="#1B3557",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(f"{out_prefix}_boxplots.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_prefix}_boxplots.pdf", bbox_inches="tight")
    plt.close(fig)


def draw_volcano(stats, predictor, sample_summary, out_prefix):
    plot_df = stats.copy()
    plot_df["q_value_adj"] = plot_df["q_value"].replace(0, 1e-300)
    plot_df["neg_log10_q"] = -np.log10(plot_df["q_value_adj"])

    fig, ax = plt.subplots(figsize=(8.5, 6.3))
    colors = np.where(plot_df["cliffs_delta"] >= 0, "#6F6BAA", "#2C7FB8")
    ax.scatter(plot_df["cliffs_delta"], plot_df["neg_log10_q"], s=80, color=colors, edgecolor="white", linewidth=0.8)

    label_offsets = {
        "alpha_propensity": (0.012, 0.10),
        "charge": (0.012, -0.03),
        "KD_mean": (0.012, 0.04),
        "KD_max": (0.012, 0.02),
    }
    for _, row in plot_df.iterrows():
        dx, dy = label_offsets.get(row["feature"], (0.008, 0.04))
        ax.text(row["cliffs_delta"] + dx, row["neg_log10_q"] + dy, row["feature"], fontsize=9, color="#1B3557")

    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    ax.axvline(0.1, color="#596579", linestyle=":", linewidth=1.2)
    ax.axvline(-0.1, color="#596579", linestyle=":", linewidth=1.2)
    ax.set_xlabel("Cliff's Delta (Amyloid vs Non-Amyloid)")
    ax.set_ylabel("-log10(q-value)")
    ax.set_title(f"Volcano plot: {predictor} amyloid vs non-amyloid features\n{sample_summary}", fontsize=14, weight="bold", color="#1B3557")
    ax.grid(color="#D7DFEA", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_volcano.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_prefix}_volcano.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    raw = load_raw_from_zarr(ZARR_PATH)
    write_prediction_count_tables(raw)
    sample_meta = pd.read_csv(SAMPLE_META_PATH, sep="\t")
    zarr_samples = set(raw["sample"].astype(str).unique())
    sample_meta = sample_meta[sample_meta["sample"].astype(str).isin(zarr_samples)].copy()

    n_samples = int(sample_meta["sample"].nunique())
    n_ac = int(sample_meta.loc[sample_meta["Group"].eq("AC"), "sample"].nunique())
    n_tol = int(sample_meta.loc[sample_meta["Group"].eq("TOL"), "sample"].nunique())
    n_subjects = int(sample_meta["subject"].nunique()) if "subject" in sample_meta.columns else np.nan
    sample_summary = f"n={n_samples} BC study samples ({n_ac} AC, {n_tol} TOL), {n_subjects} subjects"

    all_stats = []
    summary_rows = []
    for predictor in PREDICTORS:
        df = filter_for_predictor(raw, predictor)
        stats = build_stats(df, predictor)
        all_stats.append(stats)
        counts = df[predictor].value_counts(dropna=False).reindex(ORDER).fillna(0).astype(int)
        summary_rows.append(
            {
                "predictor": predictor,
                "n_samples": n_samples,
                "n_subjects": n_subjects,
                "n_AC_samples": n_ac,
                "n_TOL_samples": n_tol,
                "n_protein_rows_after_filter": int(len(df)),
                "n_amyloid_rows": int(counts["Amyloid"]),
                "n_non_amyloid_rows": int(counts["Non-Amyloid"]),
                "source_zarr": str(ZARR_PATH.relative_to(ROOT)),
            }
        )
        out_prefix = OUT / f"study_cohort_{predictor}"
        draw_boxplots(df, stats, predictor, sample_summary, out_prefix)
        draw_volcano(stats, predictor, sample_summary, out_prefix)

    pd.concat(all_stats, ignore_index=True).to_csv(OUT / "study_cohort_feature_volcano_statistics.tsv", sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(OUT / "study_cohort_predictor_feature_summary.tsv", sep="\t", index=False)

    for path in sorted(OUT.glob("study_cohort_*")):
        print(path)
    print(OUT / "study_cohort_feature_volcano_statistics.tsv")
    print(OUT / "study_cohort_predictor_feature_summary.tsv")


if __name__ == "__main__":
    main()
