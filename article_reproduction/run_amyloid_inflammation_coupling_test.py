import json
import itertools
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.stats import mannwhitneyu, pearsonr, spearmanr, ttest_ind
from statsmodels.formula.api import ols

warnings.filterwarnings("ignore")


ROOT = Path(__file__).resolve().parent
ZARR_PATH = ROOT / "project.zarr"
ANNOTATION_PATH = ROOT / "GSE287540_SraRunTable.csv"
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "amyloid_inflammation_coupling_test"
OUT.mkdir(parents=True, exist_ok=True)

TIMEPOINT = "BC"
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
GROUPS = ["AC", "TOL"]
COVERAGE_MIN = 0.5
PSEUDOCOUNT = 1.0
RNG_SEED = 20260604
N_BOOT = 1000
N_PERM = 3000

MODULES = {
    "myeloid_neutrophil": [
        "S100A8",
        "S100A9",
        "LYZ",
        "LST1",
        "FCGR3B",
        "FCGR1A",
        "CSF3R",
        "CXCR2",
        "MMP9",
        "ELANE",
        "MPO",
        "CEACAM8",
        "ITGAM",
        "TREM1",
        "MNDA",
        "AIF1",
    ],
    "interferon": [
        "ISG15",
        "IFI6",
        "IFI27",
        "IFI44",
        "IFI44L",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "MX1",
        "MX2",
        "OAS1",
        "OAS2",
        "OAS3",
        "RSAD2",
        "IRF7",
        "STAT1",
        "CXCL10",
    ],
    "antigen_presentation": [
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "HLA-DRA",
        "HLA-DRB1",
        "HLA-DPA1",
        "HLA-DPB1",
        "HLA-DQA1",
        "HLA-DQB1",
        "B2M",
        "TAP1",
        "TAP2",
        "CD74",
        "CIITA",
        "PSME1",
        "PSME2",
    ],
    "proteasome_UPR": [
        "PSMB8",
        "PSMB9",
        "PSMB10",
        "PSMB5",
        "PSMC1",
        "PSMC2",
        "PSMC3",
        "PSMC4",
        "PSMC5",
        "PSMC6",
        "PSMD1",
        "PSMD2",
        "PSMD3",
        "PSMD4",
        "PSMD7",
        "PSMD11",
        "HSPA5",
        "HSP90B1",
        "XBP1",
        "DDIT3",
        "ATF4",
        "CALR",
        "CANX",
        "PDIA3",
        "ERN1",
        "EIF2AK3",
    ],
}


def clean_enst(x):
    if pd.isna(x):
        return np.nan
    m = re.search(r"ENST\d+(?:\.\d+)?", str(x))
    return m.group(0).split(".")[0] if m else np.nan


def decode_array(values):
    out = []
    for v in values:
        if isinstance(v, bytes):
            out.append(v.decode("utf-8"))
        else:
            out.append(str(v))
    return out


def zarr_group_to_df(g):
    return pd.DataFrame({col: g[col][:] for col in g.array_keys()})


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


def zscore_series(s):
    s = pd.Series(s, dtype=float)
    sd = s.std(skipna=True, ddof=0)
    if sd == 0 or pd.isna(sd):
        return pd.Series(np.nan, index=s.index)
    return (s - s.mean(skipna=True)) / sd


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


def safe_corr(x, y, method):
    d = pd.DataFrame({"x": x, "y": y}).dropna()
    if d.shape[0] < 3 or d["x"].nunique() < 2 or d["y"].nunique() < 2:
        return np.nan, np.nan
    if method == "pearson":
        r, p = pearsonr(d["x"], d["y"])
    else:
        r, p = spearmanr(d["x"], d["y"])
    return float(r), float(p)


def fast_interaction_fit(df):
    d = df.dropna(subset=["immune_score", "amyloid_score_z", "Group"]).copy()
    if d.shape[0] < 4 or d["Group"].nunique() < 2:
        return np.nan, np.nan, np.nan
    g = d["Group"].astype(str).eq("TOL").astype(float).to_numpy()
    a = d["amyloid_score_z"].astype(float).to_numpy()
    y = d["immune_score"].astype(float).to_numpy()
    x = np.column_stack([np.ones(len(d)), a, g, a * g])
    try:
        beta = np.linalg.lstsq(x, y, rcond=None)[0]
    except Exception:
        return np.nan, np.nan, np.nan
    slope_ac = float(beta[1])
    interaction = float(beta[3])
    return slope_ac, slope_ac + interaction, interaction


def fisher_z(r):
    if pd.isna(r) or abs(r) >= 1:
        return np.nan
    return float(np.arctanh(r))


def mean_ci(values):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return np.nan, np.nan
    return float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))


def bootstrap_group_delta(df, metric, rng):
    ac = df.loc[df["Group"].eq("AC"), metric].dropna().to_numpy(float)
    tol = df.loc[df["Group"].eq("TOL"), metric].dropna().to_numpy(float)
    if len(ac) < 2 or len(tol) < 2:
        return np.nan, np.nan
    vals = []
    for _ in range(N_BOOT):
        vals.append(rng.choice(tol, len(tol), replace=True).mean() - rng.choice(ac, len(ac), replace=True).mean())
    return mean_ci(vals)


def build_bc_data():
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    annotation = pd.read_csv(ANNOTATION_PATH)
    annotation["raw_Group"] = annotation["Group"].astype(str)
    annotation["Group"] = annotation["raw_Group"].map(GROUP_MAP)
    annotation["Run"] = annotation["Run"].astype(str)

    zarr_samples = decode_array(root["layers/expression/sample_ids"][:])
    sample_meta = annotation[
        annotation["Run"].isin(zarr_samples)
        & annotation["time"].eq(TIMEPOINT)
        & annotation["Group"].isin(GROUPS)
    ].copy()
    sample_meta = sample_meta.rename(columns={"Run": "sample"})
    selected_samples = [s for s in zarr_samples if s in set(sample_meta["sample"])]
    sample_meta = sample_meta.set_index("sample").loc[selected_samples].reset_index()
    sample_meta["Group"] = pd.Categorical(sample_meta["Group"], categories=GROUPS, ordered=True)

    sample_idx = [zarr_samples.index(s) for s in selected_samples]
    transcript_ids = pd.Series(decode_array(root["layers/expression/transcript_ids"][:]), name="Transcript_ID")
    transcript_clean = transcript_ids.map(clean_enst)

    expr_ann = zarr_group_to_df(root[f"samples/{selected_samples[0]}/expression"])[
        ["transcript_id", "gene_id", "gene_name"]
    ].copy()
    expr_ann["Transcript_ID_clean"] = expr_ann["transcript_id"].map(clean_enst)
    expr_ann["gene_name"] = expr_ann["gene_name"].astype(str)
    expr_ann = expr_ann[["Transcript_ID_clean", "gene_id", "gene_name"]].drop_duplicates("Transcript_ID_clean")

    tpm = pd.DataFrame(
        root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=transcript_ids,
        columns=selected_samples,
    )
    tpm.index.name = "Transcript_ID"

    amy_dfs = []
    pf_dfs = []
    for sample in selected_samples:
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
        combined[f"{col}_z"] = zscore_series(combined[col])
    combined["Amyloid_Index_max_source"] = (
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
            Amyloid_Index_max=("Amyloid_Index_max_source", "max"),
            Amyloid_Index_mean=("Amyloid_Index_max_source", "mean"),
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
    amy_tx = amy_tx.merge(expr_ann, on="Transcript_ID_clean", how="left")
    amy_tx["is_amyloid"] = amy_tx["Consensus_collapsed"].eq("Amyloid")

    transcript_base = pd.DataFrame(
        {
            "Transcript_ID": transcript_ids,
            "Transcript_ID_clean": transcript_clean,
        }
    )

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

    sample_rows = []
    contrib_chunks = []
    meta_small = sample_meta[["sample", "subject", "time", "raw_Group", "Group"]].copy()
    for sample in selected_samples:
        g = transcript_base.copy()
        g["sample"] = sample
        g["TPM"] = tpm[sample].to_numpy()
        g = g.merge(meta_small, on="sample", how="left")
        g = g.merge(amy_tx[amy_tx["sample"].eq(sample)], on=["sample", "Transcript_ID_clean"], how="left")
        g["has_amyloid_annotation"] = g["Amyloid_Index_max"].notna()
        g["weighted_amyloid_index"] = g["TPM"].fillna(0) * g["Amyloid_Index_max"]
        g["strict_amyloid_TPM"] = np.where(g["Consensus_collapsed"].eq("Amyloid"), g["TPM"].fillna(0), 0.0)
        g["strict_amyloid_weighted_index"] = np.where(
            g["Consensus_collapsed"].eq("Amyloid"),
            g["TPM"].fillna(0) * g["Amyloid_Index_max"],
            0.0,
        )
        row = sample_summary(g)
        row["sample"] = sample
        sample_rows.append(row)
        strict = g[g["Consensus_collapsed"].eq("Amyloid")].copy()
        if len(strict):
            strict["strict_contribution"] = strict["TPM"].fillna(0) * strict["Amyloid_Index_max"]
            contrib_chunks.append(
                strict[
                    [
                        "sample",
                        "Group",
                        "Transcript_ID_clean",
                        "TPM",
                        "Amyloid_Index_max",
                        "gene_name",
                        "strict_contribution",
                    ]
                ]
            )

    sample_burden = pd.DataFrame(sample_rows)
    sample_burden = sample_burden.merge(sample_meta, on="sample", how="left")
    sample_burden["keep_coverage_qc"] = sample_burden["protein_coding_annotation_fraction"] >= COVERAGE_MIN

    tx_gene = pd.DataFrame({"Transcript_ID": transcript_ids, "Transcript_ID_clean": transcript_clean})
    tx_gene = tx_gene.merge(expr_ann, on="Transcript_ID_clean", how="left")
    tpm_gene = tpm.copy()
    tpm_gene["gene_name"] = tx_gene.set_index("Transcript_ID").loc[tpm.index, "gene_name"].values
    tpm_gene = tpm_gene.dropna(subset=["gene_name"]).groupby("gene_name")[selected_samples].sum()
    tpm_gene = tpm_gene[~tpm_gene.index.astype(str).str.lower().isin(["nan", "", "none"])]

    contrib = pd.concat(contrib_chunks, ignore_index=True) if contrib_chunks else pd.DataFrame()
    top_contrib = (
        contrib.groupby("Transcript_ID_clean", as_index=False)
        .agg(
            mean_strict_contribution_AC=("strict_contribution", lambda x: x[contrib.loc[x.index, "Group"].eq("AC")].mean()),
            mean_strict_contribution_TOL=("strict_contribution", lambda x: x[contrib.loc[x.index, "Group"].eq("TOL")].mean()),
            mean_TPM_AC=("TPM", lambda x: x[contrib.loc[x.index, "Group"].eq("AC")].mean()),
            mean_TPM_TOL=("TPM", lambda x: x[contrib.loc[x.index, "Group"].eq("TOL")].mean()),
            Amyloid_Index_max=("Amyloid_Index_max", "max"),
            gene_name=("gene_name", "first"),
            n_samples_detected=("sample", "nunique"),
        )
    )
    top_contrib["delta_strict_contribution_TOL_minus_AC"] = (
        top_contrib["mean_strict_contribution_TOL"] - top_contrib["mean_strict_contribution_AC"]
    )

    return sample_meta, sample_burden, amy_tx, tpm_gene, top_contrib


def score_modules(tpm_gene, samples):
    log_gene = np.log2(tpm_gene[samples].astype(float) + PSEUDOCOUNT)
    rows = []
    score_df = pd.DataFrame({"sample": samples})
    for module, genes in MODULES.items():
        present = [g for g in genes if g in log_gene.index]
        missing = [g for g in genes if g not in log_gene.index]
        if present:
            x = log_gene.loc[present].T
            z = x.apply(zscore_series, axis=0)
            score = z.mean(axis=1, skipna=True)
            raw = x.mean(axis=1, skipna=True)
            score_df[f"{module}_score"] = score_df["sample"].map(score)
            score_df[f"{module}_mean_log2TPM"] = score_df["sample"].map(raw)
        else:
            score_df[f"{module}_score"] = np.nan
            score_df[f"{module}_mean_log2TPM"] = np.nan
        rows.append(
            {
                "module": module,
                "n_genes_defined": len(genes),
                "n_genes_present": len(present),
                "genes_present": ",".join(present),
                "genes_missing": ",".join(missing),
            }
        )
    return score_df, pd.DataFrame(rows)


def fit_interactions(data, label, rng):
    rows = []
    modules = list(MODULES)
    df = data.copy()
    df["amyloid_score"] = df["strict_amyloid_TPM_fraction"].astype(float)
    df["amyloid_score_z"] = zscore_series(df["amyloid_score"])
    df["group_TOL"] = df["Group"].astype(str).eq("TOL").astype(int)

    for module in modules:
        score_col = f"{module}_score"
        d = df[["sample", "Group", "amyloid_score", "amyloid_score_z", score_col]].dropna().copy()
        if d.shape[0] < 6 or d["Group"].nunique() < 2:
            continue
        d = d.rename(columns={score_col: "immune_score"})
        try:
            model = ols("immune_score ~ amyloid_score_z * C(Group)", data=d).fit()
            term = "amyloid_score_z:C(Group)[T.TOL]"
            interaction_coef = float(model.params.get(term, np.nan))
            interaction_p = float(model.pvalues.get(term, np.nan))
            slope_ac = float(model.params.get("amyloid_score_z", np.nan))
            slope_tol = slope_ac + interaction_coef
            r2 = float(model.rsquared)
        except Exception:
            interaction_coef = interaction_p = slope_ac = slope_tol = r2 = np.nan

        ac = d[d["Group"].eq("AC")]
        tol = d[d["Group"].eq("TOL")]
        pearson_ac, pearson_p_ac = safe_corr(ac["amyloid_score"], ac["immune_score"], "pearson")
        pearson_tol, pearson_p_tol = safe_corr(tol["amyloid_score"], tol["immune_score"], "pearson")
        spearman_ac, spearman_p_ac = safe_corr(ac["amyloid_score"], ac["immune_score"], "spearman")
        spearman_tol, spearman_p_tol = safe_corr(tol["amyloid_score"], tol["immune_score"], "spearman")
        pearson_diff = pearson_tol - pearson_ac if np.isfinite(pearson_tol) and np.isfinite(pearson_ac) else np.nan
        spearman_diff = spearman_tol - spearman_ac if np.isfinite(spearman_tol) and np.isfinite(spearman_ac) else np.nan
        fisher_diff = fisher_z(pearson_tol) - fisher_z(pearson_ac)

        boot_interactions = []
        boot_pearson_diffs = []
        for _ in range(N_BOOT):
            parts = []
            for group in GROUPS:
                gd = d[d["Group"].eq(group)]
                parts.append(gd.iloc[rng.integers(0, len(gd), len(gd))])
            bd = pd.concat(parts, ignore_index=True)
            _, _, boot_interaction = fast_interaction_fit(bd)
            boot_interactions.append(boot_interaction)
            bac = bd[bd["Group"].eq("AC")]
            btol = bd[bd["Group"].eq("TOL")]
            br_ac, _ = safe_corr(bac["amyloid_score"], bac["immune_score"], "pearson")
            br_tol, _ = safe_corr(btol["amyloid_score"], btol["immune_score"], "pearson")
            boot_pearson_diffs.append(br_tol - br_ac if np.isfinite(br_tol) and np.isfinite(br_ac) else np.nan)
        interaction_ci_low, interaction_ci_high = mean_ci(boot_interactions)
        pearson_diff_ci_low, pearson_diff_ci_high = mean_ci(boot_pearson_diffs)

        obs_abs_interaction = abs(interaction_coef) if np.isfinite(interaction_coef) else np.nan
        obs_abs_pearson_diff = abs(pearson_diff) if np.isfinite(pearson_diff) else np.nan
        perm_interaction = []
        perm_pearson_diff = []
        n_tol = int(d["Group"].astype(str).eq("TOL").sum())
        base_labels = np.array(["AC"] * len(d), dtype=object)
        for tol_idx in itertools.combinations(range(len(d)), n_tol):
            pdx = d.copy()
            labels = base_labels.copy()
            labels[list(tol_idx)] = "TOL"
            pdx["Group"] = labels
            _, _, perm_coef = fast_interaction_fit(pdx)
            perm_interaction.append(abs(perm_coef) if np.isfinite(perm_coef) else np.nan)
            pac = pdx[pdx["Group"].eq("AC")]
            ptol = pdx[pdx["Group"].eq("TOL")]
            pr_ac, _ = safe_corr(pac["amyloid_score"], pac["immune_score"], "pearson")
            pr_tol, _ = safe_corr(ptol["amyloid_score"], ptol["immune_score"], "pearson")
            perm_pearson_diff.append(abs(pr_tol - pr_ac) if np.isfinite(pr_tol) and np.isfinite(pr_ac) else np.nan)

        perm_interaction = np.asarray(perm_interaction, dtype=float)
        perm_pearson_diff = np.asarray(perm_pearson_diff, dtype=float)
        interaction_emp_p = (
            (np.sum(perm_interaction[np.isfinite(perm_interaction)] >= obs_abs_interaction) + 1)
            / (np.sum(np.isfinite(perm_interaction)) + 1)
            if np.isfinite(obs_abs_interaction)
            else np.nan
        )
        pearson_diff_emp_p = (
            (np.sum(perm_pearson_diff[np.isfinite(perm_pearson_diff)] >= obs_abs_pearson_diff) + 1)
            / (np.sum(np.isfinite(perm_pearson_diff)) + 1)
            if np.isfinite(obs_abs_pearson_diff)
            else np.nan
        )

        rows.append(
            {
                "analysis_set": label,
                "module": module,
                "n_samples": int(d.shape[0]),
                "n_AC": int(ac.shape[0]),
                "n_TOL": int(tol.shape[0]),
                "mean_score_AC": float(ac["immune_score"].mean()),
                "mean_score_TOL": float(tol["immune_score"].mean()),
                "delta_mean_TOL_minus_AC": float(tol["immune_score"].mean() - ac["immune_score"].mean()),
                "mannwhitney_p_score_group": float(mannwhitneyu(ac["immune_score"], tol["immune_score"], alternative="two-sided").pvalue),
                "welch_p_score_group": float(ttest_ind(ac["immune_score"], tol["immune_score"], equal_var=False).pvalue),
                "cliffs_delta_TOL_vs_AC": float(cliffs_delta(tol["immune_score"], ac["immune_score"])),
                "pearson_r_AC": pearson_ac,
                "pearson_p_AC": pearson_p_ac,
                "pearson_r_TOL": pearson_tol,
                "pearson_p_TOL": pearson_p_tol,
                "pearson_r_diff_TOL_minus_AC": pearson_diff,
                "pearson_r_diff_ci_low": pearson_diff_ci_low,
                "pearson_r_diff_ci_high": pearson_diff_ci_high,
                "pearson_diff_empirical_p_two_sided": interaction_emp_p if pd.isna(pearson_diff_emp_p) else pearson_diff_emp_p,
                "spearman_r_AC": spearman_ac,
                "spearman_p_AC": spearman_p_ac,
                "spearman_r_TOL": spearman_tol,
                "spearman_p_TOL": spearman_p_tol,
                "spearman_r_diff_TOL_minus_AC": spearman_diff,
                "fisher_z_pearson_diff_TOL_minus_AC": fisher_diff,
                "ols_slope_AC": slope_ac,
                "ols_slope_TOL": slope_tol,
                "ols_interaction_coef_TOL_minus_AC": interaction_coef,
                "ols_interaction_p": interaction_p,
                "ols_interaction_empirical_p_two_sided": interaction_emp_p,
                "ols_interaction_ci_low": interaction_ci_low,
                "ols_interaction_ci_high": interaction_ci_high,
                "ols_r2": r2,
            }
        )
    return pd.DataFrame(rows), df


def group_stats(data, label, rng):
    metrics = [
        "strict_amyloid_TPM_fraction",
        "strict_amyloid_TPM",
        "n_strict_amyloid_expressed_TPM_ge_1",
        "strict_weighted_index_sum",
        "continuous_expression_weighted_amyloid_index",
    ] + [f"{m}_score" for m in MODULES]
    rows = []
    for metric in metrics:
        d = data[["Group", metric]].dropna()
        ac = d.loc[d["Group"].eq("AC"), metric].astype(float)
        tol = d.loc[d["Group"].eq("TOL"), metric].astype(float)
        if len(ac) < 2 or len(tol) < 2:
            continue
        ci_low, ci_high = bootstrap_group_delta(data, metric, rng)
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
                "delta_mean_TOL_minus_AC": tol.mean() - ac.mean(),
                "delta_mean_TOL_minus_AC_ci_low": ci_low,
                "delta_mean_TOL_minus_AC_ci_high": ci_high,
                "mannwhitney_p": mannwhitneyu(ac, tol, alternative="two-sided").pvalue,
                "welch_t_p": ttest_ind(ac, tol, equal_var=False, nan_policy="omit").pvalue,
                "cliffs_delta_TOL_vs_AC": cliffs_delta(tol, ac),
            }
        )
    return pd.DataFrame(rows)


def write_report(summary, interactions, group_summary, module_coverage, sample_scores, excluded):
    primary = interactions[
        (interactions["analysis_set"].eq("coverage_QC_protein_coding_annotation_fraction_ge_0_5"))
        & (interactions["module"].eq("myeloid_neutrophil"))
    ]
    primary_row = primary.iloc[0].to_dict() if len(primary) else {}
    burden_qc = group_summary[
        (group_summary["analysis_set"].eq("coverage_QC_protein_coding_annotation_fraction_ge_0_5"))
        & (group_summary["metric"].eq("strict_amyloid_TPM_fraction"))
    ]
    burden_row = burden_qc.iloc[0].to_dict() if len(burden_qc) else {}

    def fmt(x, digits=3):
        if x is None or pd.isna(x):
            return "NA"
        return f"{float(x):.{digits}g}"

    sample_lines = []
    for _, r in sample_scores[["sample", "subject", "raw_Group", "Group", "keep_coverage_qc", "strict_amyloid_TPM_fraction"]].iterrows():
        qc = "QC" if r["keep_coverage_qc"] else "excluded from QC"
        sample_lines.append(
            f"- `{r['sample']}` / `{r['subject']}`: {r['raw_Group']} -> {r['Group']}, {qc}, "
            f"strict fraction={fmt(r['strict_amyloid_TPM_fraction'], 4)}"
        )
    excluded_lines = (
        "\n".join(
            f"- `{r['sample']}` ({r['Group']}), protein_coding_annotation_fraction={fmt(r['protein_coding_annotation_fraction'], 4)}"
            for _, r in excluded.iterrows()
        )
        if len(excluded)
        else "- none"
    )

    report = f"""# Amyloid-Inflammation Coupling Test

## What Was Done

Method 9 tested whether baseline amyloid burden is coupled to inflammatory RNA modules, rather than only asking whether amyloid burden is higher in ChILI/TOL.

Data were recalculated from `project.zarr` for `time == "BC"` only, with `control -> AC` and `chili -> TOL`. Protein variants were collapsed to `sample x transcript` before expression weighting to avoid H1/H2 double-counting. The primary amyloid score was `strict_amyloid_TPM_fraction`, where strict amyloid means `Consensus_collapsed == "Amyloid"`. The continuous score used `Amyloid_Index_max`.

Immune module scores were mean z-scored `log2(gene TPM + 1)` across predefined blood/inflammation panels:

- myeloid/neutrophil
- interferon
- antigen presentation
- proteasome/UPR

The main model was `immune_score ~ amyloid_score_z * Group`; the interaction term asks whether the amyloid/immune slope differs in TOL/chili versus AC/control. Two analysis sets were run: all BC samples and coverage-QC subset excluding `protein_coding_annotation_fraction < 0.5`. Empirical interaction p-values were computed by exhaustive group-label permutation with the same group sizes.

## Samples Used

{chr(10).join(sample_lines)}

## QC Exclusions

{excluded_lines}

## Main Result

Primary burden endpoint in the coverage-QC subset:

- AC/control mean `strict_amyloid_TPM_fraction`: {fmt(burden_row.get('mean_AC'), 4)}
- TOL/chili mean `strict_amyloid_TPM_fraction`: {fmt(burden_row.get('mean_TOL'), 4)}
- Delta TOL - AC: {fmt(burden_row.get('delta_mean_TOL_minus_AC'), 4)}
- Bootstrap 95% CI: {fmt(burden_row.get('delta_mean_TOL_minus_AC_ci_low'), 4)} to {fmt(burden_row.get('delta_mean_TOL_minus_AC_ci_high'), 4)}
- Welch p-value: {fmt(burden_row.get('welch_t_p'), 4)}

Primary coupling test for myeloid/neutrophil activation in the coverage-QC subset:

- Pearson r AC/control: {fmt(primary_row.get('pearson_r_AC'), 4)}
- Pearson r TOL/chili: {fmt(primary_row.get('pearson_r_TOL'), 4)}
- Pearson r difference TOL - AC: {fmt(primary_row.get('pearson_r_diff_TOL_minus_AC'), 4)}
- OLS interaction coefficient TOL - AC: {fmt(primary_row.get('ols_interaction_coef_TOL_minus_AC'), 4)}
- OLS nominal interaction p-value: {fmt(primary_row.get('ols_interaction_p'), 4)}
- Exact empirical interaction p-value: {fmt(primary_row.get('ols_interaction_empirical_p_two_sided'), 4)}
- Bootstrap interaction 95% CI: {fmt(primary_row.get('ols_interaction_ci_low'), 4)} to {fmt(primary_row.get('ols_interaction_ci_high'), 4)}

## Interpretation

The burden result remains directionally consistent with the hypothesis: baseline ChILI/TOL has higher strict amyloid TPM fraction than AC/control after excluding the low-coverage sample. The coupling result is also directionally consistent with the amyloid-inflammation hypothesis: in the coverage-QC subset, amyloid burden is more strongly associated with the myeloid/neutrophil module in TOL/chili than in AC/control.

This remains discovery/exploratory, not validation. The exact label-permutation p-value is small, but the nominal OLS interaction p-value is weak and the QC TOL/chili group has only 4 samples. Across modules, use the interaction table to separate prioritization signals from claims: a positive interaction means amyloid burden is more positively associated with the immune module in TOL/chili than in AC/control; a negative interaction means the opposite.

## Limitations

- Very small sample size: all BC is 6 AC vs 5 TOL; coverage-QC is 6 AC vs 4 TOL.
- `SRR32060276` is biologically and statistically influential because it has very low protein/amyloid annotation coverage and was excluded from the primary QC subset.
- Module scores are compact predefined panels, not externally validated latent factors in this cohort.
- Interaction models with 10-11 samples are unstable; effect size and direction are more informative than p-values.
- The analysis is discovery/exploratory, not validation.

## Next Best Step

Validate the coupling signal in an independent or larger irAE blood RNA-seq cohort, using the same frozen amyloid burden endpoint and predefined immune module panels. If no external cohort is available, the next internal step is leave-one-subject-out stability for the interaction coefficients plus a counts-based module approach that adjusts for batch and treatment where degrees of freedom permit.
"""
    (OUT / "AMYLOID_INFLAMMATION_COUPLING_TEST_REPORT.md").write_text(report, encoding="utf-8")

    compact = {
        "method": "Amyloid-Inflammation Coupling Test",
        "output_dir": str(OUT),
        "primary_endpoint": "strict_amyloid_TPM_fraction",
        "primary_coupling_model": "immune_score ~ amyloid_score_z * Group",
        "summary": summary,
        "primary_burden_qc": burden_row,
        "primary_myeloid_coupling_qc": primary_row,
        "module_coverage": module_coverage.to_dict("records"),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(compact, indent=2, default=str), encoding="utf-8")


def main():
    rng = np.random.default_rng(RNG_SEED)
    sample_meta, sample_burden, amy_tx, tpm_gene, top_contrib = build_bc_data()

    all_samples = sample_burden["sample"].tolist()
    module_scores_all, module_coverage = score_modules(tpm_gene, all_samples)
    sample_scores = sample_burden.merge(module_scores_all, on="sample", how="left")

    sample_scores.to_csv(OUT / "sample_scores_amyloid_inflammation.tsv", sep="\t", index=False)
    sample_meta.to_csv(OUT / "sample_metadata_BC.tsv", sep="\t", index=False)
    sample_burden.to_csv(OUT / "sample_level_amyloid_burden_recomputed.tsv", sep="\t", index=False)
    amy_tx.to_csv(OUT / "collapsed_sample_transcript_amyloid_annotations.tsv", sep="\t", index=False)
    module_coverage.to_csv(OUT / "module_gene_coverage.tsv", sep="\t", index=False)
    top_contrib.sort_values("delta_strict_contribution_TOL_minus_AC", ascending=False).to_csv(
        OUT / "top_strict_amyloid_contributor_panel.tsv", sep="\t", index=False
    )

    analysis_sets = {
        "all_BC_samples": sample_scores.copy(),
        "coverage_QC_protein_coding_annotation_fraction_ge_0_5": sample_scores[sample_scores["keep_coverage_qc"]].copy(),
    }
    interaction_tables = []
    group_tables = []
    scored_tables = []
    for label, df in analysis_sets.items():
        interactions, scored = fit_interactions(df, label, rng)
        gs = group_stats(df, label, rng)
        interaction_tables.append(interactions)
        group_tables.append(gs)
        scored["analysis_set"] = label
        scored_tables.append(scored)
        sub = OUT / label
        sub.mkdir(exist_ok=True)
        df.to_csv(sub / "sample_scores_amyloid_inflammation.tsv", sep="\t", index=False)
        interactions.to_csv(sub / "coupling_interaction_results.tsv", sep="\t", index=False)
        gs.to_csv(sub / "group_statistics.tsv", sep="\t", index=False)

    interactions = pd.concat(interaction_tables, ignore_index=True)
    group_summary = pd.concat(group_tables, ignore_index=True)
    scored_by_set = pd.concat(scored_tables, ignore_index=True)
    interactions.to_csv(OUT / "coupling_interaction_results.tsv", sep="\t", index=False)
    group_summary.to_csv(OUT / "group_statistics.tsv", sep="\t", index=False)
    scored_by_set.to_csv(OUT / "model_input_scores_by_analysis_set.tsv", sep="\t", index=False)

    excluded = sample_scores[~sample_scores["keep_coverage_qc"]][
        ["sample", "subject", "Group", "protein_coding_annotation_fraction"]
    ].copy()
    excluded.to_csv(OUT / "coverage_qc_excluded_samples.tsv", sep="\t", index=False)

    summary = {
        "n_all_BC": int(sample_scores.shape[0]),
        "group_counts_all_BC": sample_scores["Group"].astype(str).value_counts().to_dict(),
        "n_coverage_QC": int(sample_scores["keep_coverage_qc"].sum()),
        "group_counts_coverage_QC": sample_scores.loc[sample_scores["keep_coverage_qc"], "Group"].astype(str).value_counts().to_dict(),
        "excluded_samples": excluded.to_dict("records"),
    }
    write_report(summary, interactions, group_summary, module_coverage, sample_scores, excluded)

    print("DONE")
    print(OUT)
    print(interactions[["analysis_set", "module", "pearson_r_AC", "pearson_r_TOL", "ols_interaction_coef_TOL_minus_AC", "ols_interaction_empirical_p_two_sided"]])


if __name__ == "__main__":
    main()
