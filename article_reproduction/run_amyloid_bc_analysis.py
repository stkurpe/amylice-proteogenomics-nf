import json
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")


ROOT = Path(__file__).resolve().parent
ZARR_PATH = ROOT / "project.zarr"
ANNOTATION_CANDIDATES = [
    ROOT / "GSE287540_SraRunTable.csv",
    ROOT / "SraRunTable.csv",
    Path("/Users/user841/Downloads/SraRunTable (1).csv"),
    Path("/Users/user841/Downloads/SraRunTable.csv"),
]
HPA_CANDIDATES = [
    ROOT / "subcellular_location.tsv",
    Path("/Users/user841/Downloads/subcellular_location.tsv"),
]
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili"
OUT.mkdir(parents=True, exist_ok=True)

TIMEPOINT = "BC"
GROUPS = ["TOL", "AC"]
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
GROUP_LABELS = {"AC": "control", "TOL": "chili"}
COVERAGE_MIN = 0.5
MIN_COUNT_FOR_DE = 10
MIN_SAMPLES_FOR_DE = 3
PSEUDOCOUNT = 1.0


def clean_enst(x):
    if pd.isna(x):
        return np.nan
    m = re.search(r"ENST\d+(?:\.\d+)?", str(x))
    return m.group(0).split(".")[0] if m else np.nan


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


def zscore_series(s):
    s = pd.Series(s, dtype=float)
    sd = s.std(skipna=True, ddof=0)
    if sd == 0 or pd.isna(sd):
        return pd.Series(np.nan, index=s.index)
    return (s - s.mean(skipna=True)) / sd


def zarr_group_to_df(g):
    return pd.DataFrame({col: g[col][:] for col in g.array_keys()})


def load_annotation():
    for p in ANNOTATION_CANDIDATES:
        if p.exists():
            df = pd.read_csv(p)
            if {"Run", "subject", "time", "Group"}.issubset(df.columns):
                return df, p
    raise FileNotFoundError("Could not find SraRunTable with Run/subject/time/Group columns.")


def load_gencode_annotation():
    cache = OUT / "gencode_v48_transcript_annotation.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t")

    gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz"
    try:
        gtf = pd.read_csv(
            gtf_url,
            sep="\t",
            comment="#",
            header=None,
            names=[
                "seqname",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute",
            ],
            compression="gzip",
        )
    except Exception as e:
        print(f"WARNING: could not load GENCODE from network: {e}")
        return pd.DataFrame(columns=["Transcript_ID_clean", "gene_id", "gene_name", "biotype"])

    transcripts = gtf[gtf["feature"] == "transcript"].copy()

    def attr(s, key):
        m = re.search(fr'{key} "([^"]+)"', str(s))
        return m.group(1) if m else np.nan

    transcripts["Transcript_ID_clean"] = transcripts["attribute"].map(lambda x: clean_enst(attr(x, "transcript_id")))
    transcripts["gene_id"] = transcripts["attribute"].map(lambda x: attr(x, "gene_id"))
    transcripts["gene_name"] = transcripts["attribute"].map(lambda x: attr(x, "gene_name"))
    transcripts["biotype"] = transcripts["attribute"].map(lambda x: attr(x, "transcript_type"))
    out = transcripts[["Transcript_ID_clean", "gene_id", "gene_name", "biotype"]].drop_duplicates("Transcript_ID_clean")
    out.to_csv(cache, sep="\t", index=False)
    return out


def load_hpa():
    for p in HPA_CANDIDATES:
        if p.exists():
            hpa = pd.read_csv(p, sep="\t")
            all_loc = (
                hpa.get("Main location", pd.Series("", index=hpa.index)).fillna("")
                + ";"
                + hpa.get("Additional location", pd.Series("", index=hpa.index)).fillna("")
                + ";"
                + hpa.get("Extracellular location", pd.Series("", index=hpa.index)).fillna("")
            ).str.lower()
            hpa["all_locations"] = all_loc
            return hpa, p
    return pd.DataFrame(), None


def compare_sample_metric(df, metrics, label):
    rows = []
    for metric in metrics:
        d = df.dropna(subset=[metric, "Group"])
        ac = d.loc[d["Group"] == "AC", metric].astype(float)
        tol = d.loc[d["Group"] == "TOL", metric].astype(float)
        if len(ac) < 2 or len(tol) < 2:
            continue
        rows.append(
            {
                "analysis_set": label,
                "metric": metric,
                "n_AC": len(ac),
                "n_TOL": len(tol),
                "mean_AC": ac.mean(),
                "mean_TOL": tol.mean(),
                "median_AC": ac.median(),
                "median_TOL": tol.median(),
                "delta_mean_AC_minus_TOL": ac.mean() - tol.mean(),
                "mannwhitney_p": mannwhitneyu(ac, tol, alternative="two-sided").pvalue,
                "welch_t_p": ttest_ind(ac, tol, equal_var=False, nan_policy="omit").pvalue,
                "cliffs_delta_AC_vs_TOL": cliffs_delta(ac, tol),
            }
        )
    res = pd.DataFrame(rows)
    if not res.empty:
        res["FDR_mannwhitney"] = multipletests(res["mannwhitney_p"], method="fdr_bh")[1]
        res["FDR_welch"] = multipletests(res["welch_t_p"], method="fdr_bh")[1]
    return res


def main():
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    annotation, annotation_path = load_annotation()
    annotation = annotation.copy()
    annotation["raw_Group"] = annotation["Group"].astype(str)
    annotation["Group"] = annotation["raw_Group"].map(GROUP_MAP)
    annotation["Run"] = annotation["Run"].astype(str)

    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    sample_meta = annotation[
        annotation["Run"].isin(zarr_samples)
        & annotation["time"].eq(TIMEPOINT)
        & annotation["Group"].isin(GROUPS)
    ].copy()
    sample_meta = sample_meta.rename(columns={"Run": "sample"})
    selected_samples = [s for s in zarr_samples if s in set(sample_meta["sample"])]
    sample_meta = sample_meta.set_index("sample").loc[selected_samples].reset_index()
    sample_meta["Group"] = pd.Categorical(sample_meta["Group"], categories=GROUPS, ordered=True)
    sample_meta.to_csv(OUT / "sample_metadata_BC.tsv", sep="\t", index=False)

    print("Annotation:", annotation_path)
    print("Selected BC samples:", len(selected_samples))
    print(sample_meta[["sample", "subject", "time", "raw_Group", "Group"]])

    sample_idx = [zarr_samples.index(s) for s in selected_samples]
    transcript_ids = pd.Series(root["layers/expression/transcript_ids"][:], name="Transcript_ID").astype(str)
    transcript_clean = transcript_ids.map(clean_enst)
    expr_ann = zarr_group_to_df(root[f"samples/{selected_samples[0]}/expression"])[
        ["transcript_id", "gene_id", "gene_name"]
    ].copy()
    expr_ann["Transcript_ID_clean"] = expr_ann["transcript_id"].map(clean_enst)
    expr_ann = expr_ann[["Transcript_ID_clean", "gene_id", "gene_name"]].drop_duplicates("Transcript_ID_clean")
    tpm = pd.DataFrame(root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
                       index=transcript_ids, columns=selected_samples)
    counts = pd.DataFrame(root["layers/expression/est_counts"].get_orthogonal_selection((sample_idx, slice(None))).T,
                          index=transcript_ids, columns=selected_samples)
    tpm.index.name = "Transcript_ID"
    counts.index.name = "Transcript_ID"

    # Build sample-level amyloid/protein feature table from per-sample groups to preserve string predictions.
    amy_dfs = []
    pf_dfs = []
    hla_dfs = []
    nonsense_dfs = []
    for sample in selected_samples:
        sg = root[f"samples/{sample}"]
        amy = zarr_group_to_df(sg["amyloid"])
        pf = zarr_group_to_df(sg["protein_features"])
        amy_dfs.append(amy)
        pf_dfs.append(pf)
        if "hla" in sg:
            hla_dfs.append(zarr_group_to_df(sg["hla"]))
        if "nonsense" in sg:
            nd = zarr_group_to_df(sg["nonsense"])
            if len(nd):
                nonsense_dfs.append(nd)

    amyloid_all = pd.concat(amy_dfs, ignore_index=True)
    protein_features = pd.concat(pf_dfs, ignore_index=True)
    hla = pd.concat(hla_dfs, ignore_index=True) if hla_dfs else pd.DataFrame()
    nonsense = pd.concat(nonsense_dfs, ignore_index=True) if nonsense_dfs else pd.DataFrame()

    protein_features["ID"] = protein_features["sample"].astype(str) + "|" + protein_features["protein_id"].astype(str)
    amyloid_all["ID"] = amyloid_all["sample"].astype(str) + "|" + amyloid_all["Sequence_ID"].astype(str)
    combined = protein_features.merge(amyloid_all.drop(columns=["sample"], errors="ignore"), on="ID", how="inner")
    combined["Transcript_ID_clean"] = combined["Sequence_ID"].map(clean_enst)

    # Add continuous index.
    combined["abs_charge_density"] = combined["charge_density"].abs()
    combined["pI_distance_from_7"] = (combined["pI"] - 7).abs()
    for col in ["beta_propensity", "KD_max", "aromaticity", "abs_charge_density", "pI_distance_from_7"]:
        combined[f"{col}_z"] = zscore_series(combined[col])
    combined["Amyloid_Integrative_Index"] = (
        combined["beta_propensity_z"]
        + combined["KD_max_z"]
        + combined["aromaticity_z"]
        - combined["abs_charge_density_z"]
        - combined["pI_distance_from_7_z"]
    ) / 5

    transcript_annotation = load_gencode_annotation()
    gencode_loaded = not transcript_annotation.empty
    if gencode_loaded:
        transcript_annotation = transcript_annotation.merge(
            expr_ann.rename(columns={"gene_id": "gene_id_expr", "gene_name": "gene_name_expr"}),
            on="Transcript_ID_clean",
            how="outer",
        )
        transcript_annotation["gene_id"] = transcript_annotation["gene_id"].fillna(transcript_annotation["gene_id_expr"])
        transcript_annotation["gene_name"] = transcript_annotation["gene_name"].fillna(transcript_annotation["gene_name_expr"])
        transcript_annotation = transcript_annotation[["Transcript_ID_clean", "gene_id", "gene_name", "biotype"]]
    else:
        transcript_annotation = expr_ann.copy()
        transcript_annotation["biotype"] = np.nan

    amy_tx = (
        combined.dropna(subset=["sample", "Transcript_ID_clean"])
        .groupby(["sample", "Transcript_ID_clean"], as_index=False)
        .agg(
            Amyloid_Index_max=("Amyloid_Integrative_Index", "max"),
            Amyloid_Index_mean=("Amyloid_Integrative_Index", "mean"),
            AmyloGram_Prob_max=("AmyloGram_Prob", "max"),
            AMYPred_Prob_max=("AMYPred_Prob", "max"),
            beta_propensity_mean=("beta_propensity", "mean"),
            KD_max_mean=("KD_max", "mean"),
            aromaticity_mean=("aromaticity", "mean"),
            charge_density_mean=("charge_density", "mean"),
            pI_mean=("pI", "mean"),
            protein_length_mean=("protein_length", "mean"),
            n_protein_variants=("protein_id", "nunique"),
            Consensus_collapsed=("Consensus", consensus_priority),
        )
    )
    amy_tx = amy_tx.merge(transcript_annotation, on="Transcript_ID_clean", how="left")
    if gencode_loaded:
        amy_tx = amy_tx[amy_tx["biotype"].eq("protein_coding")].copy()
    amy_tx["is_amyloid"] = amy_tx["Consensus_collapsed"].eq("Amyloid")
    amy_tx["is_non_amyloid"] = amy_tx["Consensus_collapsed"].eq("Non-Amyloid")
    amy_tx.to_csv(OUT / "collapsed_sample_transcript_amyloid_annotations.tsv", sep="\t", index=False)

    expr_long = (
        tpm.reset_index()
        .melt(id_vars="Transcript_ID", var_name="sample", value_name="TPM")
        .merge(
            counts.reset_index().melt(id_vars="Transcript_ID", var_name="sample", value_name="est_counts"),
            on=["Transcript_ID", "sample"],
            how="left",
        )
    )
    expr_long["Transcript_ID_clean"] = expr_long["Transcript_ID"].map(clean_enst)
    expr_long = expr_long.merge(sample_meta[["sample", "subject", "time", "raw_Group", "Group"]], on="sample", how="left")
    expr_long["log2_TPM"] = np.log2(expr_long["TPM"].fillna(0) + PSEUDOCOUNT)

    expr_amy = expr_long.merge(amy_tx, on=["sample", "Transcript_ID_clean"], how="left")
    expr_amy["has_amyloid_annotation"] = expr_amy["Amyloid_Index_max"].notna()
    expr_amy["weighted_amyloid_index"] = expr_amy["TPM"].fillna(0) * expr_amy["Amyloid_Index_max"]
    expr_amy["strict_amyloid_TPM"] = np.where(expr_amy["Consensus_collapsed"].eq("Amyloid"), expr_amy["TPM"].fillna(0), 0.0)
    expr_amy["strict_amyloid_weighted_index"] = np.where(
        expr_amy["Consensus_collapsed"].eq("Amyloid"),
        expr_amy["TPM"].fillna(0) * expr_amy["Amyloid_Index_max"],
        0.0,
    )
    expr_amy["amyloid_TPM"] = np.where(expr_amy["is_amyloid"].fillna(False), expr_amy["TPM"].fillna(0), 0.0)

    def sample_summary(g):
        total_tpm = g["TPM"].sum()
        ann = g[g["has_amyloid_annotation"]]
        ann_tpm = ann["TPM"].sum()
        strict_tpm = ann["strict_amyloid_TPM"].sum()
        strict_weighted = ann["strict_amyloid_weighted_index"].sum()
        weighted = ann["weighted_amyloid_index"].sum()
        return pd.Series(
            {
                "total_TPM": total_tpm,
                "protein_coding_annotated_TPM": ann_tpm,
                "protein_coding_annotation_fraction": ann_tpm / total_tpm if total_tpm else np.nan,
                "strict_amyloid_TPM": strict_tpm,
                "strict_amyloid_TPM_fraction": strict_tpm / ann_tpm if ann_tpm else np.nan,
                "strict_weighted_index_sum": strict_weighted,
                "strict_expression_weighted_amyloid_index": strict_weighted / strict_tpm if strict_tpm else np.nan,
                "continuous_weighted_index_sum": weighted,
                "continuous_expression_weighted_amyloid_index": weighted / ann_tpm if ann_tpm else np.nan,
                "n_strict_amyloid_expressed_TPM_ge_1": ((ann["Consensus_collapsed"].eq("Amyloid")) & (ann["TPM"] >= 1)).sum(),
                "n_annotated_expressed_TPM_ge_1": (ann["TPM"] >= 1).sum(),
            }
        )

    sample_burden = expr_amy.groupby("sample").apply(sample_summary).reset_index()
    sample_burden = sample_burden.merge(sample_meta, on="sample", how="left")
    sample_burden["keep_coverage_qc"] = sample_burden["protein_coding_annotation_fraction"] >= COVERAGE_MIN
    sample_burden.to_csv(OUT / "sample_level_amyloid_burden.tsv", sep="\t", index=False)

    burden_metrics = [
        "strict_amyloid_TPM",
        "strict_amyloid_TPM_fraction",
        "strict_weighted_index_sum",
        "strict_expression_weighted_amyloid_index",
        "continuous_expression_weighted_amyloid_index",
        "n_strict_amyloid_expressed_TPM_ge_1",
    ]
    burden_stats_all = compare_sample_metric(sample_burden, burden_metrics, "all_BC_samples")
    burden_stats_qc = compare_sample_metric(sample_burden[sample_burden["keep_coverage_qc"]], burden_metrics, "coverage_qc")
    burden_stats = pd.concat([burden_stats_all, burden_stats_qc], ignore_index=True)
    burden_stats.to_csv(OUT / "sample_level_burden_statistics.tsv", sep="\t", index=False)

    # Counts-based DE fallback: library-size normalized logCPM Welch.
    lib_size = counts.sum(axis=0)
    count_filter = (counts >= MIN_COUNT_FOR_DE).sum(axis=1) >= MIN_SAMPLES_FOR_DE
    counts_f = counts[count_filter].copy()
    cpm = counts_f.div(lib_size, axis=1) * 1e6
    logcpm = np.log2(cpm + PSEUDOCOUNT)
    ac_samples = sample_meta.loc[sample_meta["Group"].eq("AC"), "sample"].tolist()
    tol_samples = sample_meta.loc[sample_meta["Group"].eq("TOL"), "sample"].tolist()

    de_rows = []
    for tid, row in logcpm.iterrows():
        ac = row[ac_samples].astype(float)
        tol = row[tol_samples].astype(float)
        de_rows.append(
            {
                "Transcript_ID": tid,
                "Transcript_ID_clean": clean_enst(tid),
                "baseMean_counts": counts_f.loc[tid].mean(),
                "mean_logCPM_AC": ac.mean(),
                "mean_logCPM_TOL": tol.mean(),
                "log2FC_AC_minus_TOL": ac.mean() - tol.mean(),
                "p_value": ttest_ind(ac, tol, equal_var=False, nan_policy="omit").pvalue,
                "mean_TPM_AC": tpm.loc[tid, ac_samples].mean() if tid in tpm.index else np.nan,
                "mean_TPM_TOL": tpm.loc[tid, tol_samples].mean() if tid in tpm.index else np.nan,
            }
        )
    de = pd.DataFrame(de_rows)
    de["FDR"] = multipletests(de["p_value"].fillna(1), method="fdr_bh")[1]

    tx_global = (
        amy_tx.groupby("Transcript_ID_clean", as_index=False)
        .agg(
            Amyloid_Index_max=("Amyloid_Index_max", "max"),
            Amyloid_Index_mean=("Amyloid_Index_mean", "mean"),
            Consensus_any=("Consensus_collapsed", consensus_priority),
            AmyloGram_Prob_max=("AmyloGram_Prob_max", "max"),
            AMYPred_Prob_max=("AMYPred_Prob_max", "max"),
            gene_name=("gene_name", "first"),
            gene_id=("gene_id", "first"),
            biotype=("biotype", "first"),
        )
    )
    de = de.merge(tx_global, on="Transcript_ID_clean", how="left")
    de["strict_consensus_amyloid"] = de["Consensus_any"].eq("Amyloid")
    de["amyloid_high_top10pct"] = de["Amyloid_Index_max"] >= de["Amyloid_Index_max"].quantile(0.9)
    de["minus_log10_FDR"] = -np.log10(de["FDR"].clip(lower=np.nextafter(0, 1)))
    de.to_csv(OUT / "counts_logCPM_DE_with_amyloid_annotation.tsv", sep="\t", index=False)
    counts_f.to_csv(OUT / "counts_for_external_DESeq2.tsv", sep="\t")
    sample_meta.to_csv(OUT / "metadata_for_external_DESeq2.tsv", sep="\t", index=False)

    # Strict top contributors.
    contrib = expr_amy[expr_amy["Consensus_collapsed"].eq("Amyloid")].copy()
    contrib["strict_contribution"] = contrib["TPM"].fillna(0) * contrib["Amyloid_Index_max"]
    contrib_summary = (
        contrib.groupby("Transcript_ID_clean", as_index=False)
        .agg(
            mean_strict_contribution_AC=("strict_contribution", lambda x: x[contrib.loc[x.index, "Group"].eq("AC")].mean()),
            mean_strict_contribution_TOL=("strict_contribution", lambda x: x[contrib.loc[x.index, "Group"].eq("TOL")].mean()),
            mean_TPM_AC=("TPM", lambda x: x[contrib.loc[x.index, "Group"].eq("AC")].mean()),
            mean_TPM_TOL=("TPM", lambda x: x[contrib.loc[x.index, "Group"].eq("TOL")].mean()),
            Amyloid_Index_max=("Amyloid_Index_max", "max"),
            gene_name=("gene_name", "first"),
            gene_id=("gene_id", "first"),
            biotype=("biotype", "first"),
        )
    )
    contrib_summary["delta_strict_contribution_AC_minus_TOL"] = (
        contrib_summary["mean_strict_contribution_AC"] - contrib_summary["mean_strict_contribution_TOL"]
    )
    contrib_summary["delta_TPM_AC_minus_TOL"] = contrib_summary["mean_TPM_AC"] - contrib_summary["mean_TPM_TOL"]
    contrib_summary = contrib_summary.merge(
        de[["Transcript_ID_clean", "log2FC_AC_minus_TOL", "FDR", "baseMean_counts"]],
        on="Transcript_ID_clean",
        how="left",
    )
    contrib_summary.to_csv(OUT / "strict_amyloid_top_contributors.tsv", sep="\t", index=False)

    # HLA and nonsense summaries.
    if not hla.empty:
        hla = hla.merge(sample_meta[["sample", "Group", "raw_Group", "subject"]], on="sample", how="left")
        hla.to_csv(OUT / "hla_calls_BC.tsv", sep="\t", index=False)
        hla_summary = hla.groupby(["Group", "locus", "allele"], as_index=False).agg(n_samples=("sample", "nunique"))
        hla_summary.to_csv(OUT / "hla_allele_group_summary.tsv", sep="\t", index=False)
    if not nonsense.empty:
        nonsense["Transcript_ID_clean"] = nonsense["TRANSCRIPT_ID"].map(clean_enst)
        nonsense = nonsense.merge(sample_meta[["sample", "Group", "raw_Group", "subject"]], on="sample", how="left")
        nonsense.to_csv(OUT / "nonsense_candidates_BC.tsv", sep="\t", index=False)
        nonsense_sample = nonsense.groupby(["sample", "Group"], as_index=False).agg(n_nonsense=("TRANSCRIPT_ID", "size"))
        nonsense_stats = compare_sample_metric(nonsense_sample, ["n_nonsense"], "nonsense_candidates")
        nonsense_sample.to_csv(OUT / "nonsense_count_per_sample.tsv", sep="\t", index=False)
        nonsense_stats.to_csv(OUT / "nonsense_count_statistics.tsv", sep="\t", index=False)

    # HPA compartment annotation for top contributors if available.
    hpa, hpa_path = load_hpa()
    if not hpa.empty and "gene_name" in contrib_summary:
        def has(text, kws):
            return any(k in str(text).lower() for k in kws)

        loc = hpa[["Gene name", "Reliability", "all_locations"]].drop_duplicates("Gene name").copy()
        loc["loc_membrane"] = loc["all_locations"].map(lambda x: has(x, ["plasma membrane", "cell junction"]))
        loc["loc_cytosol"] = loc["all_locations"].map(lambda x: has(x, ["cytosol"]))
        loc["loc_nucleus"] = loc["all_locations"].map(lambda x: has(x, ["nucleoplasm", "nucleoli", "nuclear", "nuclear speckle"]))
        loc["loc_mito"] = loc["all_locations"].map(lambda x: has(x, ["mitochond"]))
        loc["loc_extracellular"] = loc["all_locations"].map(lambda x: has(x, ["extracellular", "secreted"]))
        contrib_hpa = contrib_summary.merge(loc, left_on="gene_name", right_on="Gene name", how="left")
        contrib_hpa.to_csv(OUT / "strict_amyloid_top_contributors_with_hpa.tsv", sep="\t", index=False)

    # Compact report JSON.
    report = {
        "annotation_path": str(annotation_path),
        "gencode_loaded": bool(gencode_loaded),
        "hpa_path": str(hpa_path) if hpa_path else None,
        "n_selected_samples": len(selected_samples),
        "group_counts": sample_meta["Group"].astype(str).value_counts().to_dict(),
        "low_coverage_samples": sample_burden.loc[~sample_burden["keep_coverage_qc"], ["sample", "Group", "protein_coding_annotation_fraction"]].to_dict("records"),
        "n_transcripts_de_tested": int(de.shape[0]),
        "n_strict_amyloid_de_tested": int(de["strict_consensus_amyloid"].sum()),
        "burden_stats": burden_stats.to_dict("records"),
        "top_contributors_AC": contrib_summary.sort_values("delta_strict_contribution_AC_minus_TOL", ascending=False).head(15).to_dict("records"),
        "top_contributors_TOL": contrib_summary.sort_values("delta_strict_contribution_AC_minus_TOL", ascending=True).head(15).to_dict("records"),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")

    print("DONE")
    print("Outputs:", OUT)
    print("Burden stats")
    print(burden_stats)
    print("Top AC contributors")
    print(contrib_summary.sort_values("delta_strict_contribution_AC_minus_TOL", ascending=False).head(10))
    print("Top TOL contributors")
    print(contrib_summary.sort_values("delta_strict_contribution_AC_minus_TOL", ascending=True).head(10))


if __name__ == "__main__":
    main()
