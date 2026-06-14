import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parent
ZARR_PATH = Path("/Users/user841/Projects/BIoinfServer/Results_analysis/bioinfo_zarr/validation_project.zarr")
ANNOTATION_PATH = ROOT / "GSE287540_SraRunTable.csv"
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "validation" / "two_predictor_consensus"
GENCODE_CANDIDATES = [
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "competitive_gene_set_test_gencode_v48_protein_coding" / "gencode_transcript_annotation_used.tsv",
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "rank_based_amyloid_burden_signature" / "gencode_v48_transcript_annotation.tsv",
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "gencode_v48_transcript_annotation.tsv",
    ROOT / "gencode_v48_transcript_annotation.tsv",
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
COVERAGE_MIN = 0.5
N_BOOTSTRAPS = 10_000
RANDOM_SEED = 20260606


def clean_enst(x):
    if pd.isna(x):
        return np.nan
    match = re.search(r"ENST\d+(?:\.\d+)?", str(x))
    return match.group(0).split(".")[0] if match else np.nan


def zarr_group_to_df(group):
    return pd.DataFrame({col: group[col][:] for col in group.array_keys()})


def load_gencode():
    needed = {"Transcript_ID_clean", "gene_id", "gene_name", "biotype"}
    for path in GENCODE_CANDIDATES:
        if path.exists():
            ann = pd.read_csv(path, sep="\t")
            if needed.issubset(ann.columns):
                return ann[list(needed)].drop_duplicates("Transcript_ID_clean"), path
    raise FileNotFoundError("No local GENCODE transcript annotation found.")


def normalize_pred(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip()
    if not s:
        return np.nan
    if s in {"Amyloid", "Non-Amyloid"}:
        return s
    return s


def two_predictor_consensus(amypred, amylogrampy):
    a = normalize_pred(amypred)
    b = normalize_pred(amylogrampy)
    if pd.isna(a) and pd.isna(b):
        return np.nan
    if pd.isna(a) or pd.isna(b):
        return "Partial"
    if a == b:
        return a
    return "Discordant"


def consensus_priority(values):
    vals = set(pd.Series(values).dropna().astype(str))
    if "Amyloid" in vals:
        return "Amyloid"
    if "Partial" in vals:
        return "Partial"
    if "Discordant" in vals:
        return "Discordant"
    if "Non-Amyloid" in vals:
        return "Non-Amyloid"
    return np.nan


def cliffs_delta(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan
    gt = sum(np.sum(xi > y) for xi in x)
    lt = sum(np.sum(xi < y) for xi in x)
    return (gt - lt) / (len(x) * len(y))


def bootstrap_delta_ci(ac, tol, rng):
    ac = np.asarray(ac, dtype=float)
    tol = np.asarray(tol, dtype=float)
    deltas = np.empty(N_BOOTSTRAPS, dtype=float)
    for i in range(N_BOOTSTRAPS):
        deltas[i] = rng.choice(tol, len(tol), replace=True).mean() - rng.choice(ac, len(ac), replace=True).mean()
    return np.percentile(deltas, [2.5, 50, 97.5])


def group_stats(sample_df, metrics, rng):
    rows = []
    kept = sample_df[sample_df["keep_coverage_qc"]].copy()
    for metric in metrics:
        ac = kept.loc[kept["Group"].eq("AC"), metric].astype(float).to_numpy()
        tol = kept.loc[kept["Group"].eq("TOL"), metric].astype(float).to_numpy()
        ci = bootstrap_delta_ci(ac, tol, rng)
        rows.append(
            {
                "metric": metric,
                "n_AC": len(ac),
                "n_TOL": len(tol),
                "mean_AC": float(ac.mean()),
                "mean_TOL": float(tol.mean()),
                "median_AC": float(np.median(ac)),
                "median_TOL": float(np.median(tol)),
                "delta_TOL_minus_AC": float(tol.mean() - ac.mean()),
                "mannwhitney_one_sided_TOL_gt_AC_p": float(mannwhitneyu(tol, ac, alternative="greater").pvalue),
                "cliffs_delta_TOL_vs_AC": float(cliffs_delta(tol, ac)),
                "bootstrap_delta_ci2_5": float(ci[0]),
                "bootstrap_delta_median": float(ci[1]),
                "bootstrap_delta_ci97_5": float(ci[2]),
            }
        )
    return pd.DataFrame(rows)


def markdown_table(df):
    if df.empty:
        return "_No rows._"
    display = df.copy()
    for col in display.columns:
        if pd.api.types.is_float_dtype(display[col]):
            display[col] = display[col].map(lambda x: "" if pd.isna(x) else f"{x:.6g}")
    return display.to_markdown(index=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    gencode, gencode_path = load_gencode()
    gencode_pc = gencode[gencode["biotype"].eq("protein_coding")].copy()
    protein_coding = set(gencode_pc["Transcript_ID_clean"])

    meta = pd.read_csv(ANNOTATION_PATH)
    meta["Run"] = meta["Run"].astype(str)
    meta["raw_Group"] = meta["Group"].astype(str)
    meta["Group"] = meta["raw_Group"].map(GROUP_MAP)
    sample_meta = meta[
        meta["Run"].isin(VALIDATION_SAMPLES)
        & meta["time"].eq(TIMEPOINT)
        & meta["Group"].isin(["AC", "TOL"])
    ].rename(columns={"Run": "sample"})
    sample_meta = sample_meta.set_index("sample").loc[VALIDATION_SAMPLES].reset_index()
    sample_meta.to_csv(OUT / "validation_sample_metadata.tsv", sep="\t", index=False)

    raw_frames = []
    for sample in VALIDATION_SAMPLES:
        amy = zarr_group_to_df(root[f"samples/{sample}/amyloid"])
        amy["sample"] = sample
        amy["AMYPred_Pred_norm"] = amy["AMYPred_Pred"].map(normalize_pred)
        amy["AmyloGramPy_Pred_norm"] = amy["AmyloGramPy_Pred"].map(normalize_pred)
        amy["Consensus_2pred"] = [
            two_predictor_consensus(a, b) for a, b in zip(amy["AMYPred_Pred_norm"], amy["AmyloGramPy_Pred_norm"])
        ]
        amy["Transcript_ID_clean"] = amy["Sequence_ID"].map(clean_enst)
        raw_frames.append(amy)
    raw = pd.concat(raw_frames, ignore_index=True)

    raw_counts = (
        raw.groupby(["sample", "Consensus_2pred"], dropna=False)
        .size()
        .reset_index(name="n_protein_rows")
        .merge(raw.groupby("sample").size().rename("sample_total_rows"), on="sample")
    )
    raw_counts["fraction"] = raw_counts["n_protein_rows"] / raw_counts["sample_total_rows"]
    raw_counts = raw_counts.merge(sample_meta[["sample", "Group", "subject"]], on="sample", how="left")
    raw_counts.to_csv(OUT / "two_predictor_raw_protein_row_consensus_fractions.tsv", sep="\t", index=False)

    combo = (
        raw.groupby(["Consensus_2pred", "AMYPred_Pred_norm", "AmyloGramPy_Pred_norm"], dropna=False)
        .size()
        .reset_index(name="n_protein_rows")
    )
    combo["fraction"] = combo["n_protein_rows"] / combo["n_protein_rows"].sum()
    combo.to_csv(OUT / "two_predictor_prediction_combination_counts.tsv", sep="\t", index=False)

    collapsed = (
        raw.dropna(subset=["sample", "Transcript_ID_clean"])
        .groupby(["sample", "Transcript_ID_clean"], as_index=False)
        .agg(
            Consensus_2pred_collapsed=("Consensus_2pred", consensus_priority),
            AMYPred_Pred_any=("AMYPred_Pred_norm", consensus_priority),
            AmyloGramPy_Pred_any=("AmyloGramPy_Pred_norm", consensus_priority),
        )
    )
    collapsed = collapsed.merge(gencode_pc, on="Transcript_ID_clean", how="inner")
    collapsed.to_csv(OUT / "two_predictor_collapsed_gencode_pc_annotations.tsv", sep="\t", index=False)

    collapsed_counts = (
        collapsed.groupby(["sample", "Consensus_2pred_collapsed"], dropna=False)
        .size()
        .reset_index(name="n_transcripts")
        .merge(collapsed.groupby("sample").size().rename("sample_total_transcripts"), on="sample")
    )
    collapsed_counts["fraction"] = collapsed_counts["n_transcripts"] / collapsed_counts["sample_total_transcripts"]
    collapsed_counts = collapsed_counts.merge(sample_meta[["sample", "Group", "subject"]], on="sample", how="left")
    collapsed_counts.to_csv(OUT / "two_predictor_collapsed_gencode_pc_consensus_fractions.tsv", sep="\t", index=False)

    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    sample_idx = [zarr_samples.index(s) for s in VALIDATION_SAMPLES]
    tx_clean = pd.Series(root["layers/expression/transcript_ids"][:]).astype(str).map(clean_enst)
    tpm = pd.DataFrame(
        root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=tx_clean,
        columns=VALIDATION_SAMPLES,
    )
    tpm = tpm[~tpm.index.isna()].groupby(level=0).sum()

    sample_rows = []
    labels = ["Amyloid", "Partial", "Discordant", "Non-Amyloid"]
    for sample in VALIDATION_SAMPLES:
        pc_tpm = tpm.loc[tpm.index.intersection(protein_coding), sample]
        ann = collapsed[collapsed["sample"].eq(sample)].set_index("Transcript_ID_clean")
        annotated_ids = pc_tpm.index.intersection(ann.index)
        ann = ann.loc[annotated_ids]
        annotated_tpm = pc_tpm.loc[annotated_ids]
        row = {
            "sample": sample,
            "protein_coding_TPM_gencode": float(pc_tpm.sum()),
            "protein_coding_annotated_TPM": float(annotated_tpm.sum()),
            "protein_coding_annotation_fraction": float(annotated_tpm.sum() / pc_tpm.sum()),
        }
        for label in labels:
            mask = ann["Consensus_2pred_collapsed"].eq(label)
            tpm_sum = float(annotated_tpm.loc[mask].sum())
            row[f"{label}_TPM"] = tpm_sum
            row[f"{label}_TPM_fraction"] = tpm_sum / row["protein_coding_TPM_gencode"]
            row[f"n_{label}_expressed_TPM_ge_1"] = int(((annotated_tpm >= 1.0) & mask).sum())
        sample_rows.append(row)
    sample_level = pd.DataFrame(sample_rows).merge(sample_meta, on="sample", how="left")
    sample_level["keep_coverage_qc"] = sample_level["protein_coding_annotation_fraction"] >= COVERAGE_MIN
    sample_level.to_csv(OUT / "two_predictor_sample_level_consensus_burden.tsv", sep="\t", index=False)

    metrics = [
        "Amyloid_TPM_fraction",
        "Amyloid_TPM",
        "n_Amyloid_expressed_TPM_ge_1",
        "Partial_TPM_fraction",
        "Partial_TPM",
        "n_Partial_expressed_TPM_ge_1",
        "Discordant_TPM_fraction",
        "Discordant_TPM",
        "n_Discordant_expressed_TPM_ge_1",
    ]
    stats = group_stats(sample_level, metrics, rng)
    stats.to_csv(OUT / "two_predictor_group_statistics.tsv", sep="\t", index=False)

    summary = {
        "zarr_path": str(ZARR_PATH),
        "out_dir": str(OUT),
        "gencode_annotation": str(gencode_path),
        "consensus_rule": {
            "predictors_used": ["AMYPred_Pred", "AmyloGramPy_Pred"],
            "ignored": ["AmyloGram_Pred"],
            "Amyloid": "both predictors are Amyloid",
            "Non-Amyloid": "both predictors are Non-Amyloid",
            "Discordant": "both predictors present but disagree",
            "Partial": "exactly one of the two predictors is present",
        },
        "n_samples": len(VALIDATION_SAMPLES),
        "n_samples_passing_coverage_qc": int(sample_level["keep_coverage_qc"].sum()),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    overall_raw = raw_counts.groupby("Consensus_2pred")["n_protein_rows"].sum()
    overall_raw = pd.DataFrame({"n": overall_raw, "fraction": overall_raw / overall_raw.sum()}).reset_index()
    overall_collapsed = collapsed_counts.groupby("Consensus_2pred_collapsed")["n_transcripts"].sum()
    overall_collapsed = pd.DataFrame(
        {"n": overall_collapsed, "fraction": overall_collapsed / overall_collapsed.sum()}
    ).reset_index()
    report = [
        "# Two-predictor validation consensus",
        "",
        "Consensus was recomputed using only `AMYPred_Pred` and `AmyloGramPy_Pred`; `AmyloGram_Pred` was ignored.",
        "",
        "## Rule",
        "",
        "- `Amyloid`: both predictors are `Amyloid`.",
        "- `Non-Amyloid`: both predictors are `Non-Amyloid`.",
        "- `Discordant`: both predictors are present but disagree.",
        "- `Partial`: exactly one of the two predictors is present.",
        "",
        "## Raw Protein Rows",
        "",
        markdown_table(overall_raw),
        "",
        "## GENCODE Protein-Coding Collapsed Transcripts",
        "",
        markdown_table(overall_collapsed),
        "",
        "## AC vs TOL Statistics",
        "",
        markdown_table(stats),
    ]
    (OUT / "TWO_PREDICTOR_CONSENSUS_REPORT.md").write_text("\n".join(report) + "\n")


if __name__ == "__main__":
    main()
