import json
import math
import os
import re
import warnings
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib_codex_cache")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf
import zarr
from scipy.stats import mannwhitneyu, spearmanr, ttest_ind
from statsmodels.stats.multitest import multipletests

try:
    import gseapy as gp
except Exception:  # pragma: no cover - optional robustness dependency
    gp = None

warnings.filterwarnings("ignore")


ROOT = Path(__file__).resolve().parent
ZARR_PATH = ROOT / "project.zarr"
META_PATH = ROOT / "GSE287540_SraRunTable.csv"
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "seven_signature_panel"
OUT.mkdir(parents=True, exist_ok=True)
GENCODE_CANDIDATES = [
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "rank_based_amyloid_burden_signature" / "gencode_v48_transcript_annotation.tsv",
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "competitive_gene_set_test_gencode_v48_protein_coding" / "gencode_transcript_annotation_used.tsv",
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "extreme_contributor_concentration_index" / "gencode_transcript_annotation_used.tsv",
]

TIMEPOINT = "BC"
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
COVERAGE_MIN = 0.5
BOOT_N = 10000
SEED = 20260606
PRIMARY_LABEL = "TOL_vs_AC"
GROUP_PALETTE = {"AC": "#3B6EA8", "TOL": "#B24A3B"}
GROUP_PALETTE_LIST = [GROUP_PALETTE["AC"], GROUP_PALETTE["TOL"]]

SIGNATURES = {
    "Amyloid_score": [
        "AIMP1", "AP2B1", "CLEC2D", "ILF2", "LGMN", "MRPS21",
        "NLRC5", "PIK3C2A", "PSMD8", "RARB", "S100A8",
        "SEPTIN7", "SHE", "SLAIN2", "SNX10", "STK17B",
        "STMP1", "TACC1", "TIA1", "TIMM10", "YLPM1",
        "ZNF580", "ZNF672",
    ],
    "UPR_score": [
        "ASNS", "ATF3", "ATF4", "BHLHA15", "CTH", "DDIT3", "DERL2",
        "DERL3", "DNAJB11", "DNAJB9", "DNAJC3", "EDEM1", "EIF2AK3",
        "EXOSC1", "FICD", "GFPT1", "HERPUD1", "HSPA5", "HYOU1",
        "PDIA4", "PPP1R15A", "PREB", "SEC31A", "SEC61B", "SEC61G",
        "SEL1L", "SELENOS", "SERP1", "SYVN1", "TRIB3", "TRIM25",
        "UBE2J1", "VCP", "WFS1", "WIPI1", "XBP1",
    ],
    "Inflammatory_score": [
        "CCL5", "CCR5", "CD274", "CD3D", "CD3E", "CD8A", "CIITA",
        "CTLA4", "CXCL10", "CXCL11", "CXCL13", "CXCL9", "GZMA",
        "GZMB", "HLA-DRA", "HLA-DRB1", "HLA-E", "IDO1", "IL2RG",
        "ITGAL", "LAG3", "NKG7", "PDCD1", "PRF1", "PTPRC", "STAT1",
        "TAGAP",
    ],
    "IFNG_score": ["CXCL10", "CXCL11", "CXCL9", "STAT1"],
    "IFNA_score": ["CXCL10", "CXCL9", "EOMES", "GZMA", "GZMB", "IFNG", "TBX21"],
    "Myeloid_score": ["CXCL1", "CXCL2", "CXCL3", "CXCL8", "IL6", "PTGS2"],
}

SIG_ORDER = list(SIGNATURES)
BOXPLOT_SIG_ORDER = [sig for sig in SIG_ORDER if sig != "Amyloid_score"]
SIGNATURE_LABELS = {
    "Amyloid_score": "Amyloid",
    "UPR_score": "UPR",
    "Inflammatory_score": "Inflammatory",
    "IFNG_score": "IFNG",
    "IFNA_score": "IFNA",
    "Myeloid_score": "Myeloid",
}
IMMUNE_SIGS = [
    "UPR_score",
    "Inflammatory_score",
    "IFNG_score",
    "IFNA_score",
    "Myeloid_score",
]
BURDEN_METRICS = [
    "Amyloidogenic_Burden",
    "Amyloidogenic_Fraction",
    "KDmax_weighted_Amyloidogenic_Fraction",
]
BURDEN_METRIC_LABELS = {
    "Amyloidogenic_Burden": "Burden",
    "Amyloidogenic_Fraction": "Fraction",
    "KDmax_weighted_Amyloidogenic_Fraction": "KDmax-weighted\nFraction",
}
CORRELATION_SCORE_LABELS = {
    "UPR_score": "UPR",
    "Inflammatory_score": "Inflammatory",
    "IFNG_score": "IFNG",
    "IFNA_score": "IFNA",
    "Myeloid_score": "Myeloid",
}


def clean_enst(x):
    if pd.isna(x):
        return np.nan
    m = re.search(r"ENST\d+(?:\.\d+)?", str(x))
    return m.group(0).split(".")[0] if m else np.nan


def zarr_group_to_df(g):
    return pd.DataFrame({col: g[col][:] for col in g.array_keys()})


def load_gencode_annotation():
    for path in GENCODE_CANDIDATES:
        if path.exists():
            ann = pd.read_csv(path, sep="\t")
            required = {"Transcript_ID_clean", "gene_name", "biotype"}
            if required.issubset(ann.columns):
                ann = ann.copy()
                ann["Transcript_ID_clean"] = ann["Transcript_ID_clean"].astype(str)
                ann["gene_name"] = ann["gene_name"].astype(str)
                ann["biotype"] = ann["biotype"].astype(str)
                return ann[["Transcript_ID_clean", "gene_name", "biotype"]].drop_duplicates("Transcript_ID_clean"), path
    raise FileNotFoundError("No local GENCODE transcript annotation with Transcript_ID_clean/gene_name/biotype was found.")


def consensus_priority(values):
    vals = set(pd.Series(values).dropna().astype(str))
    for v in ["Amyloid", "Partial", "Discordant", "Non-Amyloid"]:
        if v in vals:
            return v
    return np.nan


def zscore_frame(df):
    mu = df.mean(axis=1)
    sd = df.std(axis=1, ddof=0)
    return df.sub(mu, axis=0).div(sd.replace(0, np.nan), axis=0)


def cohen_d_tolved(ac, tol):
    ac = pd.Series(ac).dropna().astype(float)
    tol = pd.Series(tol).dropna().astype(float)
    if len(ac) < 2 or len(tol) < 2:
        return np.nan
    pooled = math.sqrt(((len(ac) - 1) * ac.var(ddof=1) + (len(tol) - 1) * tol.var(ddof=1)) / (len(ac) + len(tol) - 2))
    if pooled == 0:
        return np.nan
    return (tol.mean() - ac.mean()) / pooled


def hedges_g(d, n1, n2):
    if not np.isfinite(d):
        return np.nan
    df = n1 + n2 - 2
    if df <= 1:
        return np.nan
    return d * (1 - 3 / (4 * df - 1))


def format_p_value(p):
    if not np.isfinite(p):
        return "p=NA"
    if p < 0.001:
        return "p<0.001"
    return f"p={p:.2g}"


def bootstrap_delta_ci(ac, tol, rng, n=BOOT_N):
    ac = np.asarray(pd.Series(ac).dropna().astype(float))
    tol = np.asarray(pd.Series(tol).dropna().astype(float))
    if len(ac) == 0 or len(tol) == 0:
        return np.nan, np.nan
    vals = np.empty(n)
    for i in range(n):
        vals[i] = rng.choice(tol, len(tol), replace=True).mean() - rng.choice(ac, len(ac), replace=True).mean()
    return np.percentile(vals, [2.5, 97.5])


def bootstrap_hedges_ci(ac, tol, rng, n=BOOT_N):
    ac = np.asarray(pd.Series(ac).dropna().astype(float))
    tol = np.asarray(pd.Series(tol).dropna().astype(float))
    if len(ac) < 2 or len(tol) < 2:
        return np.nan, np.nan
    vals = np.empty(n)
    for i in range(n):
        ac_b = rng.choice(ac, len(ac), replace=True)
        tol_b = rng.choice(tol, len(tol), replace=True)
        d = cohen_d_tolved(ac_b, tol_b)
        vals[i] = hedges_g(d, len(ac_b), len(tol_b))
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return np.nan, np.nan
    return np.percentile(vals, [2.5, 97.5])


def bootstrap_corr_diff_ci(df, immune, rng, n=BOOT_N):
    ac = df[df["Group"].eq("AC")][["Amyloid_score", immune]].dropna()
    tol = df[df["Group"].eq("TOL")][["Amyloid_score", immune]].dropna()
    if len(ac) < 3 or len(tol) < 3:
        return np.nan, np.nan
    vals = np.empty(n)
    for i in range(n):
        ac_b = ac.iloc[rng.integers(0, len(ac), len(ac))]
        tol_b = tol.iloc[rng.integers(0, len(tol), len(tol))]
        r_ac = spearmanr(ac_b["Amyloid_score"], ac_b[immune]).statistic
        r_tol = spearmanr(tol_b["Amyloid_score"], tol_b[immune]).statistic
        vals[i] = r_tol - r_ac
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return np.nan, np.nan
    return np.percentile(vals, [2.5, 97.5])


def run_gseapy_method(method, expr, signatures):
    if gp is None:
        return pd.DataFrame(index=expr.columns), "gseapy_not_available"
    try:
        fn = gp.ssgsea if method == "ssGSEA" else gp.gsva
        res = fn(
            data=expr,
            gene_sets=signatures,
            outdir=None,
            min_size=1,
            max_size=1000,
            threads=1,
            seed=SEED,
            no_plot=True,
            verbose=False,
        )
        out = res.res2d.copy()
        term_col = "Term" if "Term" in out.columns else "Name"
        sample_col = "Name" if term_col == "Term" else "Sample"
        score_col = "NES" if "NES" in out.columns else "ES"
        if sample_col not in out.columns:
            sample_col = "Sample"
        wide = out.pivot(index=sample_col, columns=term_col, values=score_col)
        wide = wide.reindex(index=expr.columns, columns=SIG_ORDER)
        wide = wide.apply(pd.to_numeric, errors="coerce")
        return wide, "ok"
    except Exception as e:
        return pd.DataFrame(index=expr.columns), f"failed: {e}"


def singscore(expr, signatures):
    ranks = expr.rank(axis=0, method="average", ascending=True)
    n = ranks.shape[0]
    out = pd.DataFrame(index=expr.columns)
    for sig, genes in signatures.items():
        present = [g for g in genes if g in ranks.index]
        if not present:
            out[sig] = np.nan
            continue
        raw = ranks.loc[present].mean(axis=0)
        out[sig] = (raw - 1) / (n - 1)
    return out


def differential_gene_rank(log_gene, sample_info):
    ac_samples = sample_info.loc[sample_info["Group"].eq("AC"), "sample"].tolist()
    tol_samples = sample_info.loc[sample_info["Group"].eq("TOL"), "sample"].tolist()
    rows = []
    for gene, row in log_gene.iterrows():
        ac = row[ac_samples].astype(float)
        tol = row[tol_samples].astype(float)
        test = ttest_ind(tol, ac, equal_var=False, nan_policy="omit")
        t_stat = test.statistic
        p_val = test.pvalue
        logfc = tol.mean() - ac.mean()
        if not np.isfinite(t_stat):
            t_stat = 0.0 if not np.isfinite(logfc) else logfc
        rows.append(
            {
                "gene": gene,
                "mean_log2_TPM_plus1_AC": ac.mean(),
                "mean_log2_TPM_plus1_TOL": tol.mean(),
                "delta_log2_TPM_plus1_TOL_minus_AC": logfc,
                "welch_t_TOL_vs_AC": t_stat,
                "welch_p_value": p_val if np.isfinite(p_val) else 1.0,
                "rank_metric": t_stat,
            }
        )
    rank_df = pd.DataFrame(rows)
    rank_df["rank_metric"] = pd.to_numeric(rank_df["rank_metric"], errors="coerce").fillna(0)
    rank_df["FDR"] = multipletests(rank_df["welch_p_value"].fillna(1), method="fdr_bh")[1]
    rank_df = rank_df.sort_values("rank_metric", ascending=False)
    return rank_df


def run_preranked_gsea(log_gene, scores):
    gsea_dir = OUT / "gsea_prerank"
    gsea_dir.mkdir(exist_ok=True)
    summary = []
    statuses = {}
    for label, sample_info in [(PRIMARY_LABEL, scores)]:
        rank_df = differential_gene_rank(log_gene, sample_info)
        rank_df.to_csv(gsea_dir / f"{label}_ranked_genes.tsv", sep="\t", index=False)
        rnk = rank_df[["gene", "rank_metric"]]
        if gp is None:
            statuses[label] = "gseapy_not_available"
            continue
        try:
            res = gp.prerank(
                rnk=rnk,
                gene_sets=SIGNATURES,
                outdir=None,
                min_size=1,
                max_size=1000,
                permutation_num=1000,
                weight=1.0,
                ascending=False,
                threads=1,
                seed=SEED,
                no_plot=True,
                verbose=False,
            )
            table = res.res2d.copy()
            table["analysis_set"] = label
            summary.append(table)
            statuses[label] = "ok"
        except Exception as e:
            statuses[label] = f"failed: {e}"
    if summary:
        gsea = pd.concat(summary, ignore_index=True)
        rename = {
            "Term": "signature",
            "Name": "comparison",
            "ES": "enrichment_score",
            "NES": "normalized_enrichment_score",
            "NOM p-val": "nominal_p_value",
            "FDR q-val": "FDR_q_value",
            "FWER p-val": "FWER_p_value",
            "Tag %": "tag_percent",
            "Gene %": "gene_percent",
            "Lead_genes": "leading_edge_genes",
        }
        gsea = gsea.rename(columns={k: v for k, v in rename.items() if k in gsea.columns})
        first_cols = [
            "analysis_set",
            "signature",
            "enrichment_score",
            "normalized_enrichment_score",
            "nominal_p_value",
            "FDR_q_value",
            "FWER_p_value",
            "tag_percent",
            "gene_percent",
            "leading_edge_genes",
        ]
        ordered = [c for c in first_cols if c in gsea.columns] + [c for c in gsea.columns if c not in first_cols]
        gsea = gsea[ordered]
    else:
        gsea = pd.DataFrame()
    return gsea, statuses


def load_expression_and_metadata():
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    meta = pd.read_csv(META_PATH)
    meta["raw_Group"] = meta["Group"].astype(str)
    meta["Group"] = meta["raw_Group"].map(GROUP_MAP)
    meta["Run"] = meta["Run"].astype(str)

    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    meta = meta[
        meta["Run"].isin(zarr_samples)
        & meta["time"].eq(TIMEPOINT)
        & meta["Group"].isin(["AC", "TOL"])
    ].copy()
    selected = [s for s in zarr_samples if s in set(meta["Run"])]
    meta = meta.rename(columns={"Run": "sample"}).set_index("sample").loc[selected].reset_index()

    sample_idx = [zarr_samples.index(s) for s in selected]
    tpm = pd.DataFrame(
        root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=pd.Series(root["layers/expression/transcript_ids"][:]).astype(str),
        columns=selected,
    )
    gencode, gencode_path = load_gencode_annotation()
    tx_ann = pd.DataFrame({"Transcript_ID": tpm.index.astype(str)})
    tx_ann["Transcript_ID_clean"] = tx_ann["Transcript_ID"].map(clean_enst)
    tx_ann = tx_ann.merge(gencode, on="Transcript_ID_clean", how="left")
    protein_coding_tx = tx_ann["biotype"].eq("protein_coding")
    tpm = tpm.loc[protein_coding_tx.values].copy()
    tx_ann = tx_ann.loc[protein_coding_tx.values].copy()
    tpm["gene_name"] = tx_ann["gene_name"].values
    gene_tpm = tpm.dropna(subset=["gene_name"]).groupby("gene_name")[selected].sum()
    log_gene = np.log2(gene_tpm + 1)
    universe = pd.DataFrame(
        {
            "gene_name": log_gene.index,
            "n_protein_coding_transcripts_in_expression": tpm.dropna(subset=["gene_name"]).groupby("gene_name").size().reindex(log_gene.index).values,
            "gencode_annotation_path": str(gencode_path),
        }
    )
    universe.to_csv(OUT / "protein_coding_expression_gene_universe.tsv", sep="\t", index=False)

    return root, meta, selected, log_gene


def recompute_coverage_qc(root, meta, selected, log_gene):
    # This rebuilds the coverage flag from project.zarr amyloid/protein predictions,
    # using Consensus_collapsed == Amyloid and Amyloid_Index_max for burden summaries.
    transcript_ids = pd.Series(root["layers/expression/transcript_ids"][:]).astype(str)
    sample_idx = [list(root["layers/expression/sample_ids"][:]).index(s) for s in selected]
    tpm = pd.DataFrame(
        root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=transcript_ids,
        columns=selected,
    )
    first_expr = zarr_group_to_df(root[f"samples/{selected[0]}/expression"])[["transcript_id", "gene_id", "gene_name"]].copy()
    first_expr["Transcript_ID_clean"] = first_expr["transcript_id"].map(clean_enst)
    expr_ann = first_expr[["Transcript_ID_clean", "gene_id", "gene_name"]].drop_duplicates("Transcript_ID_clean")
    gencode, _ = load_gencode_annotation()
    expr_ann = expr_ann.merge(
        gencode[["Transcript_ID_clean", "gene_name", "biotype"]].rename(columns={"gene_name": "gene_name_gencode"}),
        on="Transcript_ID_clean",
        how="left",
    )
    expr_ann["gene_name"] = expr_ann["gene_name_gencode"].fillna(expr_ann["gene_name"])
    expr_ann = expr_ann[["Transcript_ID_clean", "gene_id", "gene_name", "biotype"]]

    amy_dfs = []
    pf_dfs = []
    for sample in selected:
        sg = root[f"samples/{sample}"]
        amy_dfs.append(zarr_group_to_df(sg["amyloid"]))
        pf_dfs.append(zarr_group_to_df(sg["protein_features"]))
    amyloid_all = pd.concat(amy_dfs, ignore_index=True)
    protein_features = pd.concat(pf_dfs, ignore_index=True)
    protein_features["ID"] = protein_features["sample"].astype(str) + "|" + protein_features["protein_id"].astype(str)
    amyloid_all["ID"] = amyloid_all["sample"].astype(str) + "|" + amyloid_all["Sequence_ID"].astype(str)
    combined = protein_features.merge(amyloid_all.drop(columns=["sample"], errors="ignore"), on="ID", how="inner")
    combined["Transcript_ID_clean"] = combined["Sequence_ID"].map(clean_enst)

    combined["abs_charge_density"] = combined["charge_density"].abs()
    combined["pI_distance_from_7"] = (combined["pI"] - 7).abs()
    for col in ["beta_propensity", "KD_max", "aromaticity", "abs_charge_density", "pI_distance_from_7"]:
        sd = combined[col].std(ddof=0)
        combined[f"{col}_z"] = (combined[col] - combined[col].mean()) / sd if sd else np.nan
    combined["Amyloid_Integrative_Index"] = (
        combined["beta_propensity_z"]
        + combined["KD_max_z"]
        + combined["aromaticity_z"]
        - combined["abs_charge_density_z"]
        - combined["pI_distance_from_7_z"]
    ) / 5

    amy_tx = (
        combined.dropna(subset=["sample", "Transcript_ID_clean"])
        .groupby(["sample", "Transcript_ID_clean"], as_index=False)
        .agg(
            Amyloid_Index_max=("Amyloid_Integrative_Index", "max"),
            Consensus_collapsed=("Consensus", consensus_priority),
        )
        .merge(expr_ann, on="Transcript_ID_clean", how="left")
    )
    amy_tx = amy_tx[amy_tx["biotype"].eq("protein_coding")].copy()

    expr_long = tpm.reset_index(names="Transcript_ID").melt(id_vars="Transcript_ID", var_name="sample", value_name="TPM")
    expr_long["Transcript_ID_clean"] = expr_long["Transcript_ID"].map(clean_enst)
    expr_long = expr_long.merge(expr_ann[["Transcript_ID_clean", "biotype"]], on="Transcript_ID_clean", how="left")
    expr_long = expr_long[expr_long["biotype"].eq("protein_coding")].copy()
    expr_amy = expr_long.merge(amy_tx, on=["sample", "Transcript_ID_clean"], how="left")
    expr_amy["has_amyloid_annotation"] = expr_amy["Amyloid_Index_max"].notna()
    expr_amy["strict_amyloid_TPM"] = np.where(expr_amy["Consensus_collapsed"].eq("Amyloid"), expr_amy["TPM"], 0.0)
    expr_amy["strict_weighted_index"] = expr_amy["strict_amyloid_TPM"] * expr_amy["Amyloid_Index_max"].fillna(0)

    rows = []
    for sample, g in expr_amy.groupby("sample"):
        total_tpm = g["TPM"].sum()
        ann_tpm = g.loc[g["has_amyloid_annotation"], "TPM"].sum()
        strict_tpm = g["strict_amyloid_TPM"].sum()
        rows.append(
            {
                "sample": sample,
                "total_TPM": total_tpm,
                "protein_coding_annotated_TPM": ann_tpm,
                "protein_coding_annotation_fraction": ann_tpm / total_tpm if total_tpm else np.nan,
                "strict_amyloid_TPM": strict_tpm,
                "strict_amyloid_TPM_fraction": strict_tpm / ann_tpm if ann_tpm else np.nan,
                "strict_expression_weighted_amyloid_index": (
                    g["strict_weighted_index"].sum() / strict_tpm if strict_tpm else np.nan
                ),
            }
        )
    qc = pd.DataFrame(rows).merge(meta, on="sample", how="left")
    qc["keep_coverage_qc"] = qc["protein_coding_annotation_fraction"] >= COVERAGE_MIN
    return qc


def score_signatures(log_gene, meta, qc):
    z = zscore_frame(log_gene)
    scores = meta.copy()
    coverage_rows = []
    for sig, genes in SIGNATURES.items():
        unique_genes = list(dict.fromkeys(genes))
        present = [g for g in unique_genes if g in z.index]
        missing = [g for g in unique_genes if g not in z.index]
        scores[sig] = z.loc[present].mean(axis=0).reindex(scores["sample"]).values if present else np.nan
        coverage_rows.append(
            {
                "signature": sig,
                "signature_size": len(unique_genes),
                "genes_detected": len(present),
                "genes_missing": len(missing),
                "detected_gene_symbols": ",".join(present),
                "missing_gene_symbols": ",".join(missing),
            }
        )

    ssgsea_scores, ssgsea_status = run_gseapy_method("ssGSEA", log_gene, SIGNATURES)
    gsva_scores, gsva_status = run_gseapy_method("GSVA", log_gene, SIGNATURES)
    sing_scores = singscore(log_gene, SIGNATURES)

    for prefix, df in [("ssGSEA", ssgsea_scores), ("GSVA", gsva_scores), ("singscore", sing_scores)]:
        for sig in SIG_ORDER:
            scores[f"{prefix}_{sig}"] = df[sig].reindex(scores["sample"]).values if sig in df else np.nan

    scores = scores.merge(
        qc[
            [
                "sample",
                "protein_coding_annotation_fraction",
                "keep_coverage_qc",
                "strict_amyloid_TPM",
                "strict_amyloid_TPM_fraction",
                "strict_expression_weighted_amyloid_index",
            ]
        ],
        on="sample",
        how="left",
    )
    return scores, pd.DataFrame(coverage_rows), {"ssGSEA": ssgsea_status, "GSVA": gsva_status, "singscore": "ok"}


def collapse_pc_class(values):
    vals = set(pd.Series(values).dropna().astype(str))
    if "Amyloidogenic" in vals:
        return "Amyloidogenic"
    if "Intermediate" in vals:
        return "Intermediate"
    if "Non-Amyloidogenic" in vals:
        return "Non-Amyloidogenic"
    return np.nan


def compute_amyloidogenic_burden_metrics(root, scores):
    rows = []
    for sample in scores["sample"].astype(str):
        expr = zarr_group_to_df(root[f"samples/{sample}/expression"])[["transcript_id", "tpm"]].copy()
        expr["Transcript_ID_clean"] = expr["transcript_id"].map(clean_enst)
        expr["log2_TPM_plus1"] = np.log2(pd.to_numeric(expr["tpm"], errors="coerce").fillna(0) + 1)
        expr = expr[["Transcript_ID_clean", "log2_TPM_plus1"]].dropna().drop_duplicates("Transcript_ID_clean")

        amy = zarr_group_to_df(root[f"samples/{sample}/amyloid"])
        pf = zarr_group_to_df(root[f"samples/{sample}/protein_features"])[["protein_id", "KD_max"]].copy()
        amy = amy.merge(pf, left_on="Sequence_ID", right_on="protein_id", how="left")
        amy = amy[amy["ProteinCoding_BothPredictors_Q17_Q83_Class"].notna()].copy()
        amy["Transcript_ID_clean"] = amy["Transcript_ID_clean"].map(clean_enst)
        amy["KD_max"] = pd.to_numeric(amy["KD_max"], errors="coerce")

        collapsed = (
            amy.dropna(subset=["Transcript_ID_clean"])
            .groupby("Transcript_ID_clean", as_index=False)
            .agg(
                ProteinCoding_BothPredictors_Q17_Q83_Class=("ProteinCoding_BothPredictors_Q17_Q83_Class", collapse_pc_class),
                KD_max=("KD_max", "max"),
            )
        )
        tx = collapsed.merge(expr, on="Transcript_ID_clean", how="left")
        tx["log2_TPM_plus1"] = tx["log2_TPM_plus1"].fillna(0)
        tx["is_amyloidogenic"] = tx["ProteinCoding_BothPredictors_Q17_Q83_Class"].eq("Amyloidogenic")
        all_log_sum = tx["log2_TPM_plus1"].sum()
        amy_log_sum = tx.loc[tx["is_amyloidogenic"], "log2_TPM_plus1"].sum()
        kd_weight = tx["KD_max"] * tx["log2_TPM_plus1"]
        kd_all = kd_weight.sum(skipna=True)
        kd_amy = kd_weight[tx["is_amyloidogenic"]].sum(skipna=True)

        rows.append(
            {
                "sample": sample,
                "n_protein_coding_classified_transcripts": int(tx.shape[0]),
                "n_amyloidogenic_transcripts": int(tx["is_amyloidogenic"].sum()),
                "Amyloidogenic_Burden": amy_log_sum,
                "Amyloidogenic_Fraction": amy_log_sum / all_log_sum if all_log_sum else np.nan,
                "KDmax_weighted_Amyloidogenic_Fraction": kd_amy / kd_all if kd_all else np.nan,
            }
        )
    return pd.DataFrame(rows).merge(
        scores[["sample", "Group", "subject", "protein_coding_annotation_fraction"] + IMMUNE_SIGS],
        on="sample",
        how="left",
    )


def burden_score_correlations(burden_scores):
    rows = []
    for group, df in burden_scores.groupby("Group", observed=False):
        for metric in BURDEN_METRICS:
            for score in IMMUNE_SIGS:
                d = df[[metric, score]].dropna()
                rho = p = np.nan
                if len(d) >= 3 and d[metric].nunique() > 1 and d[score].nunique() > 1:
                    sp = spearmanr(d[metric], d[score])
                    rho, p = sp.statistic, sp.pvalue
                rows.append(
                    {
                        "Group": group,
                        "metric": metric,
                        "score": score,
                        "n_samples": len(d),
                        "spearman_rho": rho,
                        "spearman_p_value": p,
                        "significant_p_lt_0_05": bool(np.isfinite(p) and p < 0.05),
                    }
                )
    return pd.DataFrame(rows)


def save_burden_correlation_heatmaps(corr_df):
    for group in ["AC", "TOL"]:
        sub = corr_df[corr_df["Group"].eq(group)].copy()
        rho = sub.pivot(index="metric", columns="score", values="spearman_rho").reindex(index=BURDEN_METRICS, columns=IMMUNE_SIGS)
        pval = sub.pivot(index="metric", columns="score", values="spearman_p_value").reindex(index=BURDEN_METRICS, columns=IMMUNE_SIGS)
        sig = pval < 0.05
        labels = rho.copy().astype(object)
        for idx in labels.index:
            for col in labels.columns:
                labels.loc[idx, col] = f"{rho.loc[idx, col]:.2f}" if bool(sig.loc[idx, col]) and pd.notna(rho.loc[idx, col]) else ""
        masked = rho.where(sig)
        fig, ax = plt.subplots(figsize=(8.2, 4.2))
        sns.heatmap(
            masked.rename(index=BURDEN_METRIC_LABELS, columns=CORRELATION_SCORE_LABELS),
            cmap="vlag",
            center=0,
            vmin=-1,
            vmax=1,
            annot=labels.rename(index=BURDEN_METRIC_LABELS, columns=CORRELATION_SCORE_LABELS),
            fmt="",
            linewidths=1,
            linecolor="white",
            cbar_kws={"label": "Spearman rho\n(significant cells only)", "shrink": 0.8},
            annot_kws={"fontsize": 13, "fontweight": "bold"},
            ax=ax,
        )
        ax.set_title(f"{group}: Amyloidogenic burden correlations (p < 0.05 shown)", fontsize=14, pad=10)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=30, labelsize=11)
        ax.tick_params(axis="y", rotation=0, labelsize=11)
        fig.tight_layout()
        fig.savefig(OUT / f"amyloidogenic_burden_correlation_heatmap_{group}.png", dpi=300)
        fig.savefig(OUT / f"amyloidogenic_burden_correlation_heatmap_{group}.pdf")
        plt.close(fig)


def save_significant_burden_regplots(burden_scores, corr_df):
    sig_corr = corr_df[corr_df["significant_p_lt_0_05"].astype(bool)].copy()

    def draw_panel(ax, row):
        group = row["Group"]
        metric = row["metric"]
        score = row["score"]
        all_data = burden_scores[[metric, score, "Group", "sample"]].dropna().copy()
        line_data = all_data[all_data["Group"].eq(group)].copy()
        sns.scatterplot(
            data=all_data,
            x=metric,
            y=score,
            hue="Group",
            hue_order=["AC", "TOL"],
            palette=GROUP_PALETTE,
            s=58,
            edgecolor="#2B2B2B",
            linewidth=0.6,
            alpha=0.92,
            ax=ax,
        )
        sns.regplot(
            data=line_data,
            x=metric,
            y=score,
            scatter=False,
            line_kws={"color": GROUP_PALETTE[group], "linewidth": 2.0},
            ci=None,
            ax=ax,
        )
        ax.set_title(
            f"{BURDEN_METRIC_LABELS[metric].replace(chr(10), ' ')} vs {CORRELATION_SCORE_LABELS[score]}\n"
            f"{group}: rho={row['spearman_rho']:.2f}, {format_p_value(row['spearman_p_value'])}, n={int(row['n_samples'])}",
            fontsize=10.5,
        )
        ax.set_xlabel(BURDEN_METRIC_LABELS[metric].replace("\n", " "))
        ax.set_ylabel(CORRELATION_SCORE_LABELS[score])
        ax.grid(True, axis="both", color="#E5E5E5", linewidth=0.8)
        if ax.legend_:
            ax.legend_.remove()

    for group in ["AC", "TOL"]:
        group_corr = sig_corr[sig_corr["Group"].eq(group)].copy()
        out_png = OUT / f"amyloidogenic_burden_significant_regplots_{group}.png"
        out_pdf = OUT / f"amyloidogenic_burden_significant_regplots_{group}.pdf"
        if group_corr.empty:
            fig, ax = plt.subplots(figsize=(5.5, 3.8))
            ax.axis("off")
            ax.text(
                0.5,
                0.5,
                f"{group}: no significant Spearman correlations (p < 0.05)",
                ha="center",
                va="center",
                fontsize=12,
            )
            fig.tight_layout()
            fig.savefig(out_png, dpi=300)
            fig.savefig(out_pdf)
            plt.close(fig)
            continue

        n_panels = len(group_corr)
        n_cols = min(3, n_panels)
        n_rows = math.ceil(n_panels / n_cols)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.7 * n_cols, 3.8 * n_rows), squeeze=False)
        for ax, (_, row) in zip(axes.flat, group_corr.iterrows()):
            draw_panel(ax, row)
        for ax in axes.flat[n_panels:]:
            ax.axis("off")
        handles = [
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=GROUP_PALETTE["AC"], markeredgecolor="#2B2B2B", markersize=8, label="AC"),
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=GROUP_PALETTE["TOL"], markeredgecolor="#2B2B2B", markersize=8, label="TOL"),
        ]
        fig.legend(handles=handles, loc="lower center", frameon=False, bbox_to_anchor=(0.5, -0.005), ncol=2)
        fig.suptitle(f"{group}: significant amyloidogenic burden correlations", fontsize=14, y=0.995)
        fig.tight_layout(rect=[0, 0.04, 1, 0.97])
        fig.savefig(out_png, dpi=300)
        fig.savefig(out_pdf)
        plt.close(fig)

    if not sig_corr.empty:
        sig_corr = sig_corr.sort_values(["Group", "metric", "score"]).reset_index(drop=True)
        n_panels = len(sig_corr)
        n_cols = 3
        n_rows = math.ceil(n_panels / n_cols)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.7 * n_cols, 3.8 * n_rows), squeeze=False)
        for ax, (_, row) in zip(axes.flat, sig_corr.iterrows()):
            draw_panel(ax, row)
        for ax in axes.flat[n_panels:]:
            ax.axis("off")
        handles = [
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=GROUP_PALETTE["AC"], markeredgecolor="#2B2B2B", markersize=8, label="AC"),
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=GROUP_PALETTE["TOL"], markeredgecolor="#2B2B2B", markersize=8, label="TOL"),
            plt.Line2D([0], [0], color=GROUP_PALETTE["AC"], linewidth=2.0, label="AC significant fit"),
            plt.Line2D([0], [0], color=GROUP_PALETTE["TOL"], linewidth=2.0, label="TOL significant fit"),
        ]
        fig.legend(handles=handles, loc="lower center", frameon=False, bbox_to_anchor=(0.5, -0.005), ncol=4)
        fig.suptitle("Significant amyloidogenic burden correlations by group", fontsize=14, y=0.995)
        fig.tight_layout(rect=[0, 0.04, 1, 0.97])
        fig.savefig(OUT / "amyloidogenic_burden_significant_regplots_combined.png", dpi=300)
        fig.savefig(OUT / "amyloidogenic_burden_significant_regplots_combined.pdf")
        plt.close(fig)


def save_publication_summary_panel():
    panels = [
        OUT / "signature_scores_boxplots_AC_vs_TOL.png",
        OUT / "amyloidogenic_burden_significant_regplots_combined.png",
    ]
    missing = [str(path) for path in panels if not path.exists()]
    if missing:
        raise FileNotFoundError("Cannot build summary panel; missing files: " + ", ".join(missing))

    fig = plt.figure(figsize=(16, 22), constrained_layout=False)
    gs = fig.add_gridspec(3, 2, height_ratios=[0.95, 0.72, 1.45], hspace=0.23, wspace=0.18)

    ax_box = fig.add_subplot(gs[0, :])
    ax_box.imshow(plt.imread(panels[0]))
    ax_box.axis("off")

    corr_df = pd.read_csv(OUT / "amyloidogenic_burden_correlations.tsv", sep="\t")
    heat_axes = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    cbar_ax = fig.add_axes([0.925, 0.535, 0.012, 0.12])
    for i, (group, ax) in enumerate(zip(["AC", "TOL"], heat_axes)):
        sub = corr_df[corr_df["Group"].eq(group)].copy()
        rho = sub.pivot(index="metric", columns="score", values="spearman_rho").reindex(index=BURDEN_METRICS, columns=IMMUNE_SIGS)
        pval = sub.pivot(index="metric", columns="score", values="spearman_p_value").reindex(index=BURDEN_METRICS, columns=IMMUNE_SIGS)
        sig = pval < 0.05
        labels = rho.copy().astype(object)
        for idx in labels.index:
            for col in labels.columns:
                labels.loc[idx, col] = f"{rho.loc[idx, col]:.2f}" if bool(sig.loc[idx, col]) and pd.notna(rho.loc[idx, col]) else ""
        sns.heatmap(
            rho.where(sig).rename(index=BURDEN_METRIC_LABELS, columns=CORRELATION_SCORE_LABELS),
            cmap="vlag",
            center=0,
            vmin=-1,
            vmax=1,
            annot=labels.rename(index=BURDEN_METRIC_LABELS, columns=CORRELATION_SCORE_LABELS),
            fmt="",
            linewidths=1,
            linecolor="white",
            cbar=i == 1,
            cbar_ax=cbar_ax if i == 1 else None,
            cbar_kws={"label": "Spearman rho\n(significant cells only)"},
            annot_kws={"fontsize": 13, "fontweight": "bold"},
            ax=ax,
        )
        ax.set_title(f"{group}: amyloidogenic burden correlations (p < 0.05 shown)", fontsize=12, pad=8)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=30, labelsize=10)
        ax.tick_params(axis="y", rotation=0, labelsize=10)

    ax_reg = fig.add_subplot(gs[2, :])
    ax_reg.imshow(plt.imread(panels[1]))
    ax_reg.axis("off")

    fig.suptitle(
        "Baseline TOL vs AC: protein-coding signature scores and amyloidogenic burden coupling",
        fontsize=18,
        y=0.992,
    )
    fig.subplots_adjust(top=0.975, bottom=0.01, left=0.045, right=0.91)
    fig.savefig(OUT / "seven_signature_publication_panel.png", dpi=300)
    fig.savefig(OUT / "seven_signature_publication_panel.pdf")
    plt.close(fig)


def comparison_stats(scores):
    rng = np.random.default_rng(SEED)
    rows = []
    for label, df in [(PRIMARY_LABEL, scores)]:
        for sig in SIG_ORDER:
            ac = df[df["Group"].eq("AC")][sig].dropna()
            tol = df[df["Group"].eq("TOL")][sig].dropna()
            if len(ac) == 0 or len(tol) == 0:
                continue
            d = cohen_d_tolved(ac, tol)
            ci_low, ci_high = bootstrap_delta_ci(ac, tol, rng)
            g_ci_low, g_ci_high = bootstrap_hedges_ci(ac, tol, rng)
            rows.append(
                {
                    "analysis_set": label,
                    "signature": sig,
                    "n_AC": len(ac),
                    "n_TOL": len(tol),
                    "mean_AC": ac.mean(),
                    "mean_TOL": tol.mean(),
                    "median_AC": ac.median(),
                    "median_TOL": tol.median(),
                    "delta_mean_TOL_minus_AC": tol.mean() - ac.mean(),
                    "direction": "TOL>AC" if tol.mean() > ac.mean() else "TOL<AC",
                    "cohens_d_TOL_vs_AC": d,
                    "hedges_g_TOL_vs_AC": hedges_g(d, len(ac), len(tol)),
                    "welch_p_value": ttest_ind(tol, ac, equal_var=False, nan_policy="omit").pvalue,
                    "mannwhitney_p_value": mannwhitneyu(tol, ac, alternative="two-sided").pvalue,
                    "bootstrap_delta_mean_ci_low": ci_low,
                    "bootstrap_delta_mean_ci_high": ci_high,
                    "bootstrap_hedges_g_ci_low": g_ci_low,
                    "bootstrap_hedges_g_ci_high": g_ci_high,
                    "signature_size": len(dict.fromkeys(SIGNATURES[sig])),
                    "genes_detected": int(scores.attrs["coverage"].set_index("signature").loc[sig, "genes_detected"]),
                }
            )
    return pd.DataFrame(rows)


def coupling_models(scores):
    rng = np.random.default_rng(SEED + 1)
    model_rows = []
    corr_rows = []
    for label, df in [(PRIMARY_LABEL, scores)]:
        d = df.copy()
        d["Group"] = pd.Categorical(d["Group"], categories=["AC", "TOL"])
        for immune in IMMUNE_SIGS:
            model_df = d[["sample", "Group", "Amyloid_score", immune]].dropna()
            if model_df["Group"].nunique() == 2 and len(model_df) >= 6:
                fit = smf.ols(f"{immune} ~ Amyloid_score * Group", data=model_df).fit()
                term = "Amyloid_score:Group[T.TOL]"
                model_rows.append(
                    {
                        "analysis_set": label,
                        "immune_score": immune,
                        "n_samples": len(model_df),
                        "interaction_term": term,
                        "interaction_beta": fit.params.get(term, np.nan),
                        "interaction_p_value": fit.pvalues.get(term, np.nan),
                        "amyloid_main_beta_AC": fit.params.get("Amyloid_score", np.nan),
                        "group_TOL_beta": fit.params.get("Group[T.TOL]", np.nan),
                        "model_r_squared": fit.rsquared,
                        "positive_interaction_supports_stronger_TOL_coupling": bool(fit.params.get(term, np.nan) > 0),
                    }
                )
            for group in ["AC", "TOL"]:
                sub = model_df[model_df["Group"].eq(group)]
                rho, p = (np.nan, np.nan)
                if len(sub) >= 3:
                    sp = spearmanr(sub["Amyloid_score"], sub[immune])
                    rho, p = sp.statistic, sp.pvalue
                corr_rows.append(
                    {
                        "analysis_set": label,
                        "immune_score": immune,
                        "Group": group,
                        "n_samples": len(sub),
                        "spearman_rho": rho,
                        "spearman_p_value": p,
                    }
                )
            ac = model_df[model_df["Group"].eq("AC")]
            tol = model_df[model_df["Group"].eq("TOL")]
            r_ac = spearmanr(ac["Amyloid_score"], ac[immune]).statistic if len(ac) >= 3 else np.nan
            r_tol = spearmanr(tol["Amyloid_score"], tol[immune]).statistic if len(tol) >= 3 else np.nan
            ci_low, ci_high = bootstrap_corr_diff_ci(model_df, immune, rng)
            corr_rows.append(
                {
                    "analysis_set": label,
                    "immune_score": immune,
                    "Group": "TOL_minus_AC",
                    "n_samples": len(model_df),
                    "spearman_rho": r_tol - r_ac if np.isfinite(r_tol) and np.isfinite(r_ac) else np.nan,
                    "spearman_p_value": np.nan,
                    "correlation_difference_TOL_minus_AC": r_tol - r_ac if np.isfinite(r_tol) and np.isfinite(r_ac) else np.nan,
                    "bootstrap_ci_low": ci_low,
                    "bootstrap_ci_high": ci_high,
                }
            )
    return pd.DataFrame(model_rows), pd.DataFrame(corr_rows)


def correlation_matrices(scores):
    rows = []
    sets = [(PRIMARY_LABEL, scores)]
    for label, df in sets:
        corr = df[SIG_ORDER].corr(method="spearman")
        for s1 in SIG_ORDER:
            for s2 in SIG_ORDER:
                rows.append({"analysis_set": label, "signature_1": s1, "signature_2": s2, "spearman_rho": corr.loc[s1, s2]})
    return pd.DataFrame(rows)


def ordered_scores(scores):
    ordered = scores.copy()
    ordered["Group"] = pd.Categorical(ordered["Group"], categories=["AC", "TOL"], ordered=True)
    sort_cols = ["Group"]
    for col in ["subject", "batch", "sample"]:
        if col in ordered.columns:
            sort_cols.append(col)
    return ordered.sort_values(sort_cols).reset_index(drop=True)


def running_enrichment(rank_df, genes):
    ranked_genes = rank_df["gene"].astype(str).tolist()
    stats = rank_df["rank_metric"].astype(float).to_numpy()
    gene_set = set(genes)
    hits = np.array([g in gene_set for g in ranked_genes], dtype=bool)
    n_hits = hits.sum()
    n_miss = len(hits) - n_hits
    if n_hits == 0 or n_miss == 0:
        return np.arange(len(hits)), np.zeros(len(hits)), hits
    hit_weights = np.abs(stats) * hits
    hit_norm = hit_weights.sum()
    hit_step = hit_weights / hit_norm if hit_norm else hits / n_hits
    miss_step = (~hits) / n_miss
    es = np.cumsum(hit_step - miss_step)
    return np.arange(len(hits)), es, hits


def save_gsea_figures(gsea_df):
    if gsea_df.empty:
        return
    plot_df = gsea_df.copy()
    plot_df["normalized_enrichment_score"] = pd.to_numeric(plot_df["normalized_enrichment_score"], errors="coerce")
    plot_df["FDR_q_value"] = pd.to_numeric(plot_df["FDR_q_value"], errors="coerce")
    plot_df["signature"] = pd.Categorical(plot_df["signature"], SIG_ORDER, ordered=True)
    plot_df = plot_df.sort_values(["analysis_set", "signature"])
    pretty_sig = {s: s.replace("_score", "") for s in SIG_ORDER}
    pretty_set = {PRIMARY_LABEL: "TOL vs AC"}
    nes = plot_df.pivot(index="analysis_set", columns="signature", values="normalized_enrichment_score").reindex(
        index=list(pretty_set), columns=SIG_ORDER
    )
    fdr = plot_df.pivot(index="analysis_set", columns="signature", values="FDR_q_value").reindex(
        index=list(pretty_set), columns=SIG_ORDER
    )
    nes = nes.rename(index=pretty_set, columns=pretty_sig)
    fdr = fdr.rename(index=pretty_set, columns=pretty_sig)
    annot = nes.copy().astype(object)
    for row in nes.index:
        for col in nes.columns:
            n = nes.loc[row, col]
            q = fdr.loc[row, col]
            annot.loc[row, col] = f"NES {n:.2f}\nFDR {q:.2g}" if pd.notna(n) else ""

    fig, ax = plt.subplots(figsize=(13.5, 4.2))
    sns.heatmap(
        nes,
        cmap="vlag",
        center=0,
        vmin=-2,
        vmax=2,
        annot=annot,
        fmt="",
        annot_kws={"fontsize": 11},
        linewidths=1,
        linecolor="white",
        cbar_kws={"label": "NES (TOL-up positive)", "shrink": 0.8},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Preranked GSEA: positive NES means TOL-up enrichment", fontsize=15, pad=12)
    ax.tick_params(axis="x", rotation=25, labelsize=12)
    ax.tick_params(axis="y", rotation=0, labelsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "gsea_prerank_dotplot.png", dpi=300)
    fig.savefig(OUT / "gsea_prerank_dotplot.pdf")
    plt.close(fig)

    rank_path = OUT / "gsea_prerank" / f"{PRIMARY_LABEL}_ranked_genes.tsv"
    if not rank_path.exists():
        return
    rank_df = pd.read_csv(rank_path, sep="\t")
    top = (
        plot_df[plot_df["analysis_set"].eq(PRIMARY_LABEL)]
        .assign(abs_nes=lambda d: d["normalized_enrichment_score"].abs())
        .sort_values("abs_nes", ascending=False)
    )
    curve_sigs = [s for s in ["IFNG_score", "Inflammatory_score", "Amyloid_score", "UPR_score"] if s in set(top["signature"])]
    if len(curve_sigs) < 4:
        curve_sigs = top["signature"].head(4).tolist()

    fig = plt.figure(figsize=(18, 11))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1], height_ratios=[0.9, 1.1], wspace=0.32, hspace=0.35)
    ax_heat = fig.add_subplot(gs[0, 0])
    ax_table = fig.add_subplot(gs[1, 0])
    sub_gs = gs[:, 1].subgridspec(len(curve_sigs), 1, hspace=0.38)
    curve_axes = [fig.add_subplot(sub_gs[i, 0]) for i in range(len(curve_sigs))]

    all_bc = plot_df[plot_df["analysis_set"].eq(PRIMARY_LABEL)].copy()
    all_bc["signature"] = pd.Categorical(all_bc["signature"], SIG_ORDER, ordered=True)
    all_bc = all_bc.sort_values("signature")
    nes_one = all_bc.set_index("signature")["normalized_enrichment_score"].reindex(SIG_ORDER).to_frame("TOL vs AC").T
    fdr_one = all_bc.set_index("signature")["FDR_q_value"].reindex(SIG_ORDER).to_frame("TOL vs AC").T
    nes_one = nes_one.rename(columns={s: s.replace("_score", "") for s in SIG_ORDER})
    fdr_one = fdr_one.rename(columns={s: s.replace("_score", "") for s in SIG_ORDER})
    annot = nes_one.copy().astype(object)
    for col in nes_one.columns:
        annot.loc["TOL vs AC", col] = f"NES {nes_one.loc['TOL vs AC', col]:.2f}\nFDR {fdr_one.loc['TOL vs AC', col]:.2g}"
    sns.heatmap(
        nes_one,
        cmap="vlag",
        center=0,
        vmin=-2,
        vmax=2,
        annot=annot,
        fmt="",
        annot_kws={"fontsize": 8},
        linewidths=1,
        linecolor="white",
        cbar_kws={"label": "NES", "shrink": 0.82},
        ax=ax_heat,
    )
    ax_heat.set_title("A  GSEA summary, TOL vs AC", loc="left", fontweight="bold")
    ax_heat.set_xlabel("")
    ax_heat.set_ylabel("")
    ax_heat.tick_params(axis="x", rotation=28, labelsize=9)
    ax_heat.tick_params(axis="y", rotation=0)

    table_df = all_bc[["signature", "normalized_enrichment_score", "nominal_p_value", "FDR_q_value", "leading_edge_genes"]].copy()
    table_df["signature"] = table_df["signature"].astype(str).str.replace("_score", "", regex=False)
    table_df = table_df.sort_values("normalized_enrichment_score", ascending=False)
    ax_table.axis("off")
    ax_table.set_title("C  Preranked GSEA table and leading edge", loc="left", fontweight="bold", pad=8)
    cell_text = []
    for _, row in table_df.iterrows():
        lead = str(row["leading_edge_genes"])
        if len(lead) > 38:
            lead = lead[:35] + "..."
        cell_text.append(
            [
                row["signature"],
                f"{row['normalized_enrichment_score']:.2f}",
                f"{row['nominal_p_value']:.2g}",
                f"{row['FDR_q_value']:.2g}",
                lead,
            ]
        )
    tbl = ax_table.table(
        cellText=cell_text,
        colLabels=["Pathway", "NES", "pval", "padj", "leading edge"],
        loc="center",
        cellLoc="left",
        colLoc="left",
        colWidths=[0.22, 0.11, 0.11, 0.11, 0.45],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.scale(1, 1.3)

    for ax, sig in zip(curve_axes, curve_sigs):
        x, es, hits = running_enrichment(rank_df, SIGNATURES[sig])
        row = all_bc[all_bc["signature"].astype(str).eq(sig)].iloc[0]
        ax.plot(x, es, color="#39D300", lw=2)
        ax.axhline(0, color="0.55", lw=0.8)
        ax.axhline(es.max(), color="#ff4da0", lw=1, ls=":")
        ax.axhline(es.min(), color="#ff4da0", lw=1, ls=":")
        hit_pos = x[hits]
        ymin, ymax = ax.get_ylim()
        ax.vlines(hit_pos, ymin=ymin, ymax=ymin + (ymax - ymin) * 0.18, color="0.15", lw=0.35)
        ax.set_xlim(0, len(x))
        ax.set_ylabel("ES", fontsize=8)
        ax.set_title(
            f"B  {sig.replace('_score', '')}: NES={row['normalized_enrichment_score']:.2f}, FDR={row['FDR_q_value']:.2g}",
            loc="left",
            fontsize=9,
        )
        ax.tick_params(axis="both", labelsize=8)
    curve_axes[-1].set_xlabel("rank in TOL vs AC preranked protein-coding genes", fontsize=9)

    fig.suptitle("Preranked GSEA panel: protein-coding gene universe", fontsize=16, y=0.98)
    fig.savefig(OUT / "gsea_prerank_enrichment_panel.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "gsea_prerank_enrichment_panel.pdf", bbox_inches="tight")
    plt.close(fig)


def save_figures(scores, stats_df, corr_df, gsea_df=None):
    scores = ordered_scores(scores)
    scores.attrs.clear()
    sns.set_theme(style="whitegrid", context="talk")
    long = scores.melt(
        id_vars=["sample", "Group", "subject", "treatment", "batch", "sex", "AGE", "keep_coverage_qc"],
        value_vars=SIG_ORDER,
        var_name="signature",
        value_name="score",
    )
    long["signature"] = pd.Categorical(long["signature"], SIG_ORDER, ordered=True)
    long_box = long[long["signature"].isin(BOXPLOT_SIG_ORDER)].copy()
    long_box["signature"] = pd.Categorical(long_box["signature"], BOXPLOT_SIG_ORDER, ordered=True)

    fig, ax = plt.subplots(figsize=(15, 7))
    sns.violinplot(data=long_box, x="signature", y="score", hue="Group", split=True, inner=None, cut=0, palette=GROUP_PALETTE_LIST, ax=ax)
    sns.boxplot(data=long_box, x="signature", y="score", hue="Group", dodge=True, width=0.28, showcaps=True, fliersize=0, palette=GROUP_PALETTE_LIST, ax=ax)
    sns.stripplot(data=long_box, x="signature", y="score", hue="Group", dodge=True, color="#202020", alpha=0.68, size=4, ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], title="Group", frameon=False, loc="upper right")
    ax.set_xlabel("")
    ax.set_ylabel("Mean z-score")
    ax.set_xticks(np.arange(len(BOXPLOT_SIG_ORDER)))
    ax.set_xticklabels([SIGNATURE_LABELS[sig] for sig in BOXPLOT_SIG_ORDER], rotation=28, ha="right", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "signature_scores_group_comparison.png", dpi=300)
    fig.savefig(OUT / "signature_scores_group_comparison.pdf")
    plt.close(fig)

    stats_primary = stats_df[stats_df["analysis_set"].eq(PRIMARY_LABEL)].set_index("signature")
    fig, ax = plt.subplots(figsize=(14, 6.8))
    sns.boxplot(
        data=long_box,
        x="signature",
        y="score",
        hue="Group",
        hue_order=["AC", "TOL"],
        palette=GROUP_PALETTE,
        width=0.52,
        fliersize=0,
        linewidth=1.15,
        saturation=0.95,
        ax=ax,
    )
    sns.stripplot(
        data=long_box,
        x="signature",
        y="score",
        hue="Group",
        hue_order=["AC", "TOL"],
        dodge=True,
        color="#1A1A1A",
        alpha=0.72,
        size=4.2,
        jitter=0.12,
        ax=ax,
    )
    y_min, y_max = long_box["score"].min(), long_box["score"].max()
    y_pad = max((y_max - y_min) * 0.08, 0.18)
    ax.set_ylim(y_min - y_pad, y_max + y_pad * 2.15)
    for i, sig in enumerate(BOXPLOT_SIG_ORDER):
        if sig in stats_primary.index:
            row = stats_primary.loc[sig]
            ax.text(
                i,
                y_max + y_pad * 0.62,
                f"p={row['welch_p_value']:.2g}\ng={row['hedges_g_TOL_vs_AC']:.2f}",
                ha="center",
                va="bottom",
                fontsize=8.5,
                color="#2B2B2B",
            )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], title="", frameon=False, loc="upper left", bbox_to_anchor=(1.005, 1.02), borderaxespad=0)
    ax.set_xlabel("")
    ax.set_ylabel("Mean z-score")
    ax.set_title("Signature scores in baseline AC and TOL samples", fontsize=15, pad=12)
    ax.set_xticks(np.arange(len(BOXPLOT_SIG_ORDER)))
    ax.set_xticklabels([SIGNATURE_LABELS[sig] for sig in BOXPLOT_SIG_ORDER], rotation=28, ha="right", fontsize=10)
    ax.grid(axis="x", visible=False)
    fig.tight_layout(rect=[0, 0, 0.94, 1])
    fig.savefig(OUT / "signature_scores_boxplots_AC_vs_TOL.png", dpi=300)
    fig.savefig(OUT / "signature_scores_boxplots_AC_vs_TOL.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, len(BOXPLOT_SIG_ORDER), figsize=(13.5, 4.6), sharey=False)
    flat_axes = np.ravel(axes)
    for ax, sig in zip(flat_axes, BOXPLOT_SIG_ORDER):
        sub = long_box[long_box["signature"].eq(sig)]
        sns.boxplot(
            data=sub,
            x="Group",
            y="score",
            order=["AC", "TOL"],
            palette=GROUP_PALETTE,
            width=0.46,
            fliersize=0,
            linewidth=1.1,
            saturation=0.95,
            ax=ax,
        )
        sns.stripplot(
            data=sub,
            x="Group",
            y="score",
            order=["AC", "TOL"],
            color="#1A1A1A",
            alpha=0.72,
            size=4,
            jitter=0.13,
            ax=ax,
        )
        if sig in stats_primary.index:
            row = stats_primary.loc[sig]
            ax.set_title(
                f"{SIGNATURE_LABELS[sig]}\np={row['welch_p_value']:.2g}, g={row['hedges_g_TOL_vs_AC']:.2f}",
                fontsize=10.5,
            )
        else:
            ax.set_title(SIGNATURE_LABELS[sig], fontsize=10.5)
        ax.set_xlabel("")
        ax.set_ylabel("Mean z-score" if ax is flat_axes[0] else "")
        ax.grid(axis="x", visible=False)
    fig.suptitle("Baseline signature score comparison: TOL vs AC", fontsize=15, y=0.995)
    fig.tight_layout()
    fig.savefig(OUT / "signature_scores_boxplots_AC_vs_TOL_compact.png", dpi=300)
    fig.savefig(OUT / "signature_scores_boxplots_AC_vs_TOL_compact.pdf")
    plt.close(fig)

    heat = scores.set_index("sample")[SIG_ORDER].T
    fig, ax = plt.subplots(figsize=(13, 6))
    sns.heatmap(heat, cmap="vlag", center=0, linewidths=0.5, linecolor="white", ax=ax)
    ax.set_xlabel("Samples")
    ax.set_ylabel("Signatures")
    ax.set_title("Signature scores by baseline sample")
    sample_labels = []
    for _, row in scores.iterrows():
        flag = "QC+" if row["keep_coverage_qc"] else "QC-"
        sample_labels.append(
            f"{row['sample']}\n{row['Group']} | {row['subject']} | b{row['batch']} | {row['sex']} | {row['AGE']} | {flag}"
        )
    ax.set_xticklabels(sample_labels, rotation=45, ha="right", fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "signature_score_heatmap.png", dpi=300)
    fig.savefig(OUT / "signature_score_heatmap.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True)
    for ax, immune in zip(axes.flat, IMMUNE_SIGS):
        sns.scatterplot(
            data=scores,
            x="Amyloid_score",
            y=immune,
            hue="Group",
            hue_order=["AC", "TOL"],
            style="keep_coverage_qc",
            palette=GROUP_PALETTE,
            s=80,
            ax=ax,
        )
        sns.regplot(data=scores[scores["Group"].eq("AC")], x="Amyloid_score", y=immune, scatter=False, color=GROUP_PALETTE["AC"], ax=ax)
        sns.regplot(data=scores[scores["Group"].eq("TOL")], x="Amyloid_score", y=immune, scatter=False, color=GROUP_PALETTE["TOL"], ax=ax)
        text = []
        for group in ["AC", "TOL"]:
            sub = scores[scores["Group"].eq(group)]
            rho = spearmanr(sub["Amyloid_score"], sub[immune]).statistic if len(sub) >= 3 else np.nan
            text.append(f"{group} rho={rho:.2f}" if np.isfinite(rho) else f"{group} rho=NA")
        ax.text(0.03, 0.97, "\n".join(text), transform=ax.transAxes, va="top", ha="left", fontsize=9)
        ax.set_title(immune.replace("_score", ""))
        ax.legend_.remove() if ax.legend_ else None
    for ax in axes.flat[len(IMMUNE_SIGS):]:
        ax.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(OUT / "amyloid_immune_coupling_scatter.png", dpi=300)
    fig.savefig(OUT / "amyloid_immune_coupling_scatter.pdf")
    plt.close(fig)

    pretty = {
        "Amyloid_score": "Amyloid",
        "UPR_score": "UPR",
        "Inflammatory_score": "Inflammatory",
        "IFNG_score": "IFNG",
        "IFNA_score": "IFNA",
        "Myeloid_score": "Myeloid",
    }
    mat = corr_df[corr_df["analysis_set"].eq(PRIMARY_LABEL)].pivot(index="signature_1", columns="signature_2", values="spearman_rho").loc[SIG_ORDER, SIG_ORDER]
    mat = mat.rename(index=pretty, columns=pretty)
    fig, ax = plt.subplots(figsize=(8.5, 7.5))
    sns.heatmap(
        mat,
        cmap="vlag",
        center=0,
        vmin=-1,
        vmax=1,
        annot=True,
        fmt=".2f",
        square=True,
        cbar_kws={"label": "Spearman rho", "shrink": 0.72},
        annot_kws={"fontsize": 13},
        linewidths=1,
        linecolor="white",
        ax=ax,
    )
    ax.set_title("Signature correlation matrix, TOL vs AC analysis set", fontsize=15, pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=35, labelsize=12)
    ax.tick_params(axis="y", rotation=0, labelsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "signature_score_correlation_heatmap.png", dpi=300)
    fig.savefig(OUT / "signature_score_correlation_heatmap.pdf")
    plt.close(fig)

    forest = stats_df[stats_df["analysis_set"].eq(PRIMARY_LABEL)].copy()
    forest["signature"] = pd.Categorical(forest["signature"], SIG_ORDER, ordered=True)
    forest = forest.sort_values("signature")
    fig, ax = plt.subplots(figsize=(8, 5.5))
    y = np.arange(len(forest))
    ax.axvline(0, color="black", lw=1)
    ax.errorbar(
        forest["hedges_g_TOL_vs_AC"],
        y,
        xerr=[
            np.abs(forest["hedges_g_TOL_vs_AC"] - forest["bootstrap_hedges_g_ci_low"]),
            np.abs(forest["bootstrap_hedges_g_ci_high"] - forest["hedges_g_TOL_vs_AC"]),
        ],
        fmt="o",
        color="#2F4B7C",
        ecolor="#8CB1D9",
        capsize=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(forest["signature"].str.replace("_score", "", regex=False))
    ax.set_xlabel("Hedges g (TOL vs AC), bootstrap 95% CI")
    ax.set_title("TOL vs AC effect sizes")
    fig.tight_layout()
    fig.savefig(OUT / "signature_effect_size_forest.png", dpi=300)
    fig.savefig(OUT / "signature_effect_size_forest.pdf")
    plt.close(fig)

    if gsea_df is not None:
        save_gsea_figures(gsea_df)


def write_report(
    scores,
    coverage,
    stats_df,
    models,
    corrs,
    qc,
    method_status,
    gsea_df=None,
    gsea_status=None,
    burden_scores=None,
    burden_corrs=None,
):
    scores = ordered_scores(scores)
    qc = ordered_scores(qc)
    excluded = qc[~qc["keep_coverage_qc"]][["sample", "Group", "subject", "protein_coding_annotation_fraction"]]
    analyzed = scores[scores["keep_coverage_qc"]].copy()
    main_stats = stats_df[stats_df["analysis_set"].eq(PRIMARY_LABEL)].set_index("signature")
    upr_model = models[(models["analysis_set"].eq(PRIMARY_LABEL)) & (models["immune_score"].eq("UPR_score"))]
    infl_model = models[(models["analysis_set"].eq(PRIMARY_LABEL)) & (models["immune_score"].eq("Inflammatory_score"))]

    def fmt(x):
        return "NA" if pd.isna(x) else f"{x:.3g}"

    def img(name, alt):
        return f"![{alt}]({OUT / name})"

    lines = [
        "# Seven Signature Panel",
        "",
        "## What was done",
        "Recomputed one primary baseline comparison: `TOL/chili` versus `AC/control`. Low-coverage samples were excluded before all statistics and plots. The expression universe was explicitly filtered to GENCODE `biotype == protein_coding` transcripts before scoring or GSEA. Protein-coding transcript-level TPM was summed once per gene symbol to avoid double-counting transcript variants; scores use `log2(TPM + 1)`, per-gene z-scores across baseline samples, then the mean z-score over detected signature genes.",
        "",
        "The protein-coding gene universe used for scoring/ranking is saved in `protein_coding_expression_gene_universe.tsv`.",
        "",
        "The UPR score uses only the 36-gene Core UPR signature supplied in the request, not the 267-gene UPR regulon.",
        "",
        "Robustness methods: ssGSEA and GSVA were run through `gseapy` when available; singscore was computed as rank-normalized mean signature rank.",
        "",
        f"Method status: `{json.dumps(method_status)}`.",
        "",
        f"GSEA status: `{json.dumps(gsea_status or {})}`.",
        "",
        "## Samples used",
        f"Analyzed baseline samples after coverage QC: {analyzed.shape[0]} total; AC n={(analyzed['Group'] == 'AC').sum()}, TOL n={(analyzed['Group'] == 'TOL').sum()}.",
        "",
        analyzed[["sample", "Group", "subject", "treatment", "batch", "sex", "AGE", "protein_coding_annotation_fraction"]].to_markdown(index=False),
        "",
        "## Excluded samples",
        "Excluded samples have `protein_coding_annotation_fraction < 0.5` and are not used in any comparison, model, GSEA, or plot.",
        "",
        excluded.to_markdown(index=False) if not excluded.empty else "None.",
        "",
        "## Gene coverage",
        coverage.to_markdown(index=False),
        "",
        "## Visualizations",
        "Combined publication-style panel:",
        "",
        img("seven_signature_publication_panel.png", "Combined publication panel"),
        "",
        "Signature score comparison:",
        "",
        img("signature_scores_group_comparison.png", "Signature score comparison"),
        "",
        "Signature score boxplots:",
        "",
        img("signature_scores_boxplots_AC_vs_TOL.png", "Signature score boxplots"),
        "",
        "Compact signature score boxplots:",
        "",
        img("signature_scores_boxplots_AC_vs_TOL_compact.png", "Compact signature score boxplots"),
        "",
        "Heatmap of signature scores:",
        "",
        img("signature_score_heatmap.png", "Signature score heatmap"),
        "",
        "Amyloid-immune coupling scatter plots:",
        "",
        img("amyloid_immune_coupling_scatter.png", "Amyloid immune coupling"),
        "",
        "Correlation heatmap:",
        "",
        img("signature_score_correlation_heatmap.png", "Correlation heatmap"),
        "",
        "Effect size forest plot:",
        "",
        img("signature_effect_size_forest.png", "Effect size forest plot"),
        "",
        "Preranked GSEA dotplot:",
        "",
        img("gsea_prerank_dotplot.png", "Preranked GSEA dotplot"),
        "",
        "Preranked GSEA enrichment panel:",
        "",
        img("gsea_prerank_enrichment_panel.png", "Preranked GSEA enrichment panel"),
        "",
        "Amyloidogenic burden correlation heatmaps:",
        "",
        img("amyloidogenic_burden_correlation_heatmap_AC.png", "AC burden correlation heatmap"),
        "",
        img("amyloidogenic_burden_correlation_heatmap_TOL.png", "TOL burden correlation heatmap"),
        "",
        "Significant amyloidogenic burden correlation regplots:",
        "",
        img("amyloidogenic_burden_significant_regplots_AC.png", "AC significant burden correlation regplots"),
        "",
        img("amyloidogenic_burden_significant_regplots_TOL.png", "TOL significant burden correlation regplots"),
        "",
        img("amyloidogenic_burden_significant_regplots_combined.png", "Combined significant burden correlation regplots"),
        "",
        "## Amyloidogenic burden correlations",
        "Burden metrics were recomputed from `ProteinCoding_BothPredictors_Q17_Q83_Class == Amyloidogenic` using `log2(TPM + 1)`. Heatmaps show a rho value only for Spearman correlations with p < 0.05.",
        "",
        (
            burden_scores.to_markdown(index=False)
            if burden_scores is not None and not burden_scores.empty
            else "Amyloidogenic burden metrics were not available."
        ),
        "",
        (
            burden_corrs.to_markdown(index=False)
            if burden_corrs is not None and not burden_corrs.empty
            else "Amyloidogenic burden correlations were not available."
        ),
        "",
        "## Preranked GSEA",
        "Ranking metric: protein-coding gene-level Welch t-statistic for `TOL vs AC` using baseline `log2(TPM + 1)`. Positive NES means enrichment among TOL-up genes; negative NES means enrichment among AC-up genes.",
        "",
        (
            gsea_df[
                [
                    "analysis_set",
                    "signature",
                    "normalized_enrichment_score",
                    "nominal_p_value",
                    "FDR_q_value",
                    "leading_edge_genes",
                ]
            ].to_markdown(index=False)
            if gsea_df is not None and not gsea_df.empty
            else "GSEA was not available or returned no results."
        ),
        "",
        "## Group comparison results",
        "Positive delta/effect sizes mean TOL > AC.",
        "",
        main_stats[[
            "n_AC", "n_TOL", "mean_AC", "mean_TOL", "delta_mean_TOL_minus_AC",
            "hedges_g_TOL_vs_AC", "welch_p_value", "mannwhitney_p_value",
            "bootstrap_delta_mean_ci_low", "bootstrap_delta_mean_ci_high",
        ]].reset_index().to_markdown(index=False),
        "",
        "## Amyloid score results",
        f"Amyloid score delta TOL-AC = {fmt(main_stats.loc['Amyloid_score', 'delta_mean_TOL_minus_AC'])}, Hedges g = {fmt(main_stats.loc['Amyloid_score', 'hedges_g_TOL_vs_AC'])}, Welch p = {fmt(main_stats.loc['Amyloid_score', 'welch_p_value'])}, bootstrap CI = [{fmt(main_stats.loc['Amyloid_score', 'bootstrap_delta_mean_ci_low'])}, {fmt(main_stats.loc['Amyloid_score', 'bootstrap_delta_mean_ci_high'])}].",
        "",
        "## UPR results",
        f"UPR score delta TOL-AC = {fmt(main_stats.loc['UPR_score', 'delta_mean_TOL_minus_AC'])}, Hedges g = {fmt(main_stats.loc['UPR_score', 'hedges_g_TOL_vs_AC'])}, Welch p = {fmt(main_stats.loc['UPR_score', 'welch_p_value'])}, bootstrap CI = [{fmt(main_stats.loc['UPR_score', 'bootstrap_delta_mean_ci_low'])}, {fmt(main_stats.loc['UPR_score', 'bootstrap_delta_mean_ci_high'])}].",
        "",
        "## Immune signature results",
        "Inflammatory/IFN/myeloid signatures are reported in the comparison table above.",
        "",
        "## Amyloid-Inflammation coupling",
        f"Inflammatory interaction beta = {fmt(infl_model['interaction_beta'].iloc[0] if not infl_model.empty else np.nan)}, p = {fmt(infl_model['interaction_p_value'].iloc[0] if not infl_model.empty else np.nan)}. A positive interaction would support stronger TOL coupling.",
        "",
        "## Amyloid-UPR coupling",
        f"UPR interaction beta = {fmt(upr_model['interaction_beta'].iloc[0] if not upr_model.empty else np.nan)}, p = {fmt(upr_model['interaction_p_value'].iloc[0] if not upr_model.empty else np.nan)}.",
        "",
        "Amyloid-UPR correlation is higher in TOL than AC directionally, but this difference is not statistically secure because bootstrap CI crosses zero.",
        "",
        "Amyloid-IFNG correlation is also higher in TOL directionally, but again not statistically secure because bootstrap CI crosses zero.",
        "",
        "## Spearman coupling correlations",
        corrs.to_markdown(index=False),
        "",
        "## Interpretation",
        "The analysis is exploratory and underpowered. Interpret direction, effect size, and bootstrap intervals together; do not treat nominal p-values as biomarker-level evidence.",
        "",
        "If Amyloid_score, immune signatures, or interaction terms are not consistently positive with intervals excluding zero, the stated hypothesis is not statistically confirmed in this sample set.",
        "",
        "## Limitations",
        "Very small sample size after QC, blood bulk RNA-seq compositional confounding, tied ranks in preranked GSEA, and no multiple-testing-powered biomarker validation. Group assignment is AC/control versus TOL/chili at baseline only.",
        "",
        "## Next best step",
        "Validate the score directions in an independent cohort or with subject-level/cell-composition-adjusted models, then test whether the amyloid-expression signature adds signal beyond inflammatory cell-state markers.",
    ]
    (OUT / "analysis_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def cleanup_old_outputs():
    stale_patterns = [
        "signature_score_correlation_heatmap_all_BC.*",
        "signature_score_correlation_heatmap_AC_only.*",
        "signature_score_correlation_heatmap_TOL_only.*",
        "signature_score_correlation_heatmap_coverage_QC_subset.*",
        "test_*.png",
    ]
    for pattern in stale_patterns:
        for path in OUT.glob(pattern):
            path.unlink(missing_ok=True)
    gsea_dir = OUT / "gsea_prerank"
    if gsea_dir.exists():
        for path in gsea_dir.glob("*_ranked_genes.tsv"):
            if path.name != f"{PRIMARY_LABEL}_ranked_genes.tsv":
                path.unlink(missing_ok=True)


def main():
    cleanup_old_outputs()
    root, meta, selected, log_gene = load_expression_and_metadata()
    qc = recompute_coverage_qc(root, meta, selected, log_gene)
    scores_all, coverage, method_status = score_signatures(log_gene, meta, qc)
    scores = scores_all[scores_all["keep_coverage_qc"]].copy()
    scores.attrs["coverage"] = coverage
    stats_df = comparison_stats(scores)
    models, corrs = coupling_models(scores)
    corr_df = correlation_matrices(scores)
    gsea_df, gsea_status = run_preranked_gsea(log_gene, scores)
    burden_scores = compute_amyloidogenic_burden_metrics(root, scores)
    burden_corrs = burden_score_correlations(burden_scores)

    excluded = scores_all[~scores_all["keep_coverage_qc"]].copy()
    scores.to_csv(OUT / "signature_scores.tsv", sep="\t", index=False)
    stats_df.to_csv(OUT / "signature_score_statistics.tsv", sep="\t", index=False)
    coverage.to_csv(OUT / "signature_gene_coverage.tsv", sep="\t", index=False)
    models.to_csv(OUT / "amyloid_immune_coupling_models.tsv", sep="\t", index=False)
    corrs.to_csv(OUT / "amyloid_immune_correlations.tsv", sep="\t", index=False)
    corr_df.to_csv(OUT / "signature_score_correlation_matrix.tsv", sep="\t", index=False)
    gsea_df.to_csv(OUT / "gsea_prerank_results.tsv", sep="\t", index=False)
    burden_scores.to_csv(OUT / "amyloidogenic_burden_scores.tsv", sep="\t", index=False)
    burden_corrs.to_csv(OUT / "amyloidogenic_burden_correlations.tsv", sep="\t", index=False)
    excluded.to_csv(OUT / "excluded_samples_qc.tsv", sep="\t", index=False)
    (OUT / "method_status.json").write_text(json.dumps(method_status, indent=2), encoding="utf-8")
    (OUT / "gsea_status.json").write_text(json.dumps(gsea_status, indent=2), encoding="utf-8")

    save_figures(scores, stats_df, corr_df, gsea_df)
    save_burden_correlation_heatmaps(burden_corrs)
    save_significant_burden_regplots(burden_scores, burden_corrs)
    save_publication_summary_panel()
    write_report(
        scores_all,
        coverage,
        stats_df,
        models,
        corrs,
        scores_all,
        method_status,
        gsea_df,
        gsea_status,
        burden_scores,
        burden_corrs,
    )

    print("DONE")
    print(f"Output: {OUT}")
    print(stats_df[stats_df["analysis_set"].eq(PRIMARY_LABEL)][["signature", "delta_mean_TOL_minus_AC", "hedges_g_TOL_vs_AC", "welch_p_value"]])
    print("Excluded QC samples:")
    print(excluded[["sample", "Group", "protein_coding_annotation_fraction"]])


if __name__ == "__main__":
    main()
