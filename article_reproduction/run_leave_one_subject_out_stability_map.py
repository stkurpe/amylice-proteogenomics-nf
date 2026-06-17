import itertools
import json
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.stats import mannwhitneyu, ttest_ind

warnings.filterwarnings("ignore")


ROOT = Path(__file__).resolve().parent
ZARR_PATH = ROOT / "project.zarr"
ANNOTATION_PATH = ROOT / "GSE287540_SraRunTable.csv"
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "leave_one_subject_out_stability_map"
GENCODE_CACHE = (
    ROOT
    / "analysis_outputs"
    / "amyloid_bc_control_vs_chili"
    / "competitive_gene_set_test_gencode_v48_protein_coding"
    / "gencode_transcript_annotation_used.tsv"
)

TIMEPOINT = "BC"
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
GROUP_LABELS = {"AC": "control", "TOL": "chili"}
COVERAGE_MIN = 0.5
PSEUDOCOUNT = 1.0
RNG_SEED = 20260604
N_BOOT = 20000

PRIMARY_METRIC = "strict_amyloid_TPM_fraction"
METRICS = [
    "strict_amyloid_TPM_fraction",
    "strict_amyloid_TPM",
    "n_strict_amyloid_expressed_TPM_ge_1",
    "strict_weighted_index_sum",
    "continuous_expression_weighted_amyloid_index",
]


def clean_enst(x):
    if pd.isna(x):
        return np.nan
    m = re.search(r"ENST\d+(?:\.\d+)?", str(x))
    return m.group(0).split(".")[0] if m else np.nan


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


def zscore(s):
    s = pd.Series(s, dtype=float)
    sd = s.std(skipna=True, ddof=0)
    if sd == 0 or pd.isna(sd):
        return pd.Series(np.nan, index=s.index)
    return (s - s.mean(skipna=True)) / sd


def cliffs_delta_tOL_vs_ac(tol, ac):
    tol = np.asarray(tol, dtype=float)
    ac = np.asarray(ac, dtype=float)
    tol = tol[np.isfinite(tol)]
    ac = ac[np.isfinite(ac)]
    if len(tol) == 0 or len(ac) == 0:
        return np.nan
    gt = sum(np.sum(x > ac) for x in tol)
    lt = sum(np.sum(x < ac) for x in tol)
    return (gt - lt) / (len(tol) * len(ac))


def hedges_g_tOL_minus_ac(tol, ac):
    tol = np.asarray(tol, dtype=float)
    ac = np.asarray(ac, dtype=float)
    tol = tol[np.isfinite(tol)]
    ac = ac[np.isfinite(ac)]
    if len(tol) < 2 or len(ac) < 2:
        return np.nan
    pooled_var = ((len(tol) - 1) * tol.var(ddof=1) + (len(ac) - 1) * ac.var(ddof=1)) / (len(tol) + len(ac) - 2)
    if pooled_var <= 0:
        return np.nan
    d = (tol.mean() - ac.mean()) / math.sqrt(pooled_var)
    correction = 1 - (3 / (4 * (len(tol) + len(ac)) - 9))
    return d * correction


def bootstrap_ci_delta(tol, ac, rng):
    tol = np.asarray(tol, dtype=float)
    ac = np.asarray(ac, dtype=float)
    tol = tol[np.isfinite(tol)]
    ac = ac[np.isfinite(ac)]
    if len(tol) < 2 or len(ac) < 2:
        return np.nan, np.nan, np.nan
    deltas = np.empty(N_BOOT)
    for i in range(N_BOOT):
        deltas[i] = rng.choice(tol, size=len(tol), replace=True).mean() - rng.choice(ac, size=len(ac), replace=True).mean()
    return np.median(deltas), np.quantile(deltas, 0.025), np.quantile(deltas, 0.975)


def exact_permutation_p(values, groups):
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups)
    keep = np.isfinite(values) & pd.notna(groups)
    values = values[keep]
    groups = groups[keep]
    n_tol = int(np.sum(groups == "TOL"))
    n = len(values)
    if n_tol == 0 or n_tol == n:
        return np.nan, np.nan
    observed = values[groups == "TOL"].mean() - values[groups == "AC"].mean()
    deltas = []
    for tol_idx in itertools.combinations(range(n), n_tol):
        mask = np.zeros(n, dtype=bool)
        mask[list(tol_idx)] = True
        deltas.append(values[mask].mean() - values[~mask].mean())
    deltas = np.asarray(deltas)
    one_sided = (np.sum(deltas >= observed) + 1) / (len(deltas) + 1)
    two_sided = (np.sum(np.abs(deltas) >= abs(observed)) + 1) / (len(deltas) + 1)
    return one_sided, two_sided


def summarize_metric(df, metric, analysis_set, rng):
    d = df.dropna(subset=[metric, "Group"]).copy()
    ac = d.loc[d["Group"].eq("AC"), metric].astype(float).to_numpy()
    tol = d.loc[d["Group"].eq("TOL"), metric].astype(float).to_numpy()
    boot_median, ci_low, ci_high = bootstrap_ci_delta(tol, ac, rng)
    perm_one, perm_two = exact_permutation_p(d[metric], d["Group"])
    return {
        "analysis_set": analysis_set,
        "metric": metric,
        "n_AC": len(ac),
        "n_TOL": len(tol),
        "mean_AC": np.mean(ac) if len(ac) else np.nan,
        "mean_TOL": np.mean(tol) if len(tol) else np.nan,
        "median_AC": np.median(ac) if len(ac) else np.nan,
        "median_TOL": np.median(tol) if len(tol) else np.nan,
        "delta_mean_TOL_minus_AC": (np.mean(tol) - np.mean(ac)) if len(ac) and len(tol) else np.nan,
        "bootstrap_delta_mean_TOL_minus_AC_median": boot_median,
        "bootstrap_delta_mean_TOL_minus_AC_ci95_low": ci_low,
        "bootstrap_delta_mean_TOL_minus_AC_ci95_high": ci_high,
        "empirical_permutation_p_one_sided_TOL_gt_AC": perm_one,
        "empirical_permutation_p_two_sided": perm_two,
        "mannwhitney_p_two_sided": mannwhitneyu(tol, ac, alternative="two-sided").pvalue if len(tol) >= 2 and len(ac) >= 2 else np.nan,
        "welch_t_p_two_sided": ttest_ind(tol, ac, equal_var=False, nan_policy="omit").pvalue if len(tol) >= 2 and len(ac) >= 2 else np.nan,
        "cliffs_delta_TOL_vs_AC": cliffs_delta_tOL_vs_ac(tol, ac),
        "hedges_g_TOL_minus_AC": hedges_g_tOL_minus_ac(tol, ac),
        "direction_supports_hypothesis": bool((np.mean(tol) - np.mean(ac)) > 0) if len(ac) and len(tol) else False,
    }


def load_sample_meta(root):
    annotation = pd.read_csv(ANNOTATION_PATH)
    annotation["raw_Group"] = annotation["Group"].astype(str)
    annotation["Group"] = annotation["raw_Group"].map(GROUP_MAP)
    annotation["Run"] = annotation["Run"].astype(str)
    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    meta = annotation[
        annotation["Run"].isin(zarr_samples)
        & annotation["time"].eq(TIMEPOINT)
        & annotation["Group"].isin(["AC", "TOL"])
    ].copy()
    meta = meta.rename(columns={"Run": "sample"})
    selected = [s for s in zarr_samples if s in set(meta["sample"])]
    meta = meta.set_index("sample").loc[selected].reset_index()
    return selected, meta


def build_collapsed_annotations(root, selected_samples):
    amy_dfs = []
    pf_dfs = []
    for sample in selected_samples:
        sg = root[f"samples/{sample}"]
        amy = zarr_group_to_df(sg["amyloid"])
        pf = zarr_group_to_df(sg["protein_features"])
        amy["ID"] = amy["sample"].astype(str) + "|" + amy["Sequence_ID"].astype(str)
        pf["ID"] = pf["sample"].astype(str) + "|" + pf["protein_id"].astype(str)
        amy_dfs.append(amy)
        pf_dfs.append(pf)

    amyloid_all = pd.concat(amy_dfs, ignore_index=True)
    protein_features = pd.concat(pf_dfs, ignore_index=True)
    combined = protein_features.merge(amyloid_all.drop(columns=["sample"], errors="ignore"), on="ID", how="inner")
    combined["sample"] = combined["sample"].astype(str)
    combined["Transcript_ID_clean"] = combined["Sequence_ID"].map(clean_enst)

    combined["abs_charge_density"] = combined["charge_density"].abs()
    combined["pI_distance_from_7"] = (combined["pI"] - 7).abs()
    for col in ["beta_propensity", "KD_max", "aromaticity", "abs_charge_density", "pI_distance_from_7"]:
        combined[f"{col}_z"] = zscore(combined[col])
    combined["Amyloid_Index_max_input"] = (
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
            Amyloid_Index_max=("Amyloid_Index_max_input", "max"),
            AmyloGram_Prob_max=("AmyloGram_Prob", "max"),
            AMYPred_Prob_max=("AMYPred_Prob", "max"),
            n_protein_variants=("protein_id", "nunique"),
            Consensus_collapsed=("Consensus", consensus_priority),
        )
    )

    expr_ann = zarr_group_to_df(root[f"samples/{selected_samples[0]}/expression"])[
        ["transcript_id", "gene_id", "gene_name"]
    ].copy()
    expr_ann["Transcript_ID_clean"] = expr_ann["transcript_id"].map(clean_enst)
    expr_ann = expr_ann[["Transcript_ID_clean", "gene_id", "gene_name"]].drop_duplicates("Transcript_ID_clean")

    if GENCODE_CACHE.exists():
        gencode = pd.read_csv(GENCODE_CACHE, sep="\t")
        gencode = gencode[["Transcript_ID_clean", "biotype"]].drop_duplicates("Transcript_ID_clean")
        amy_tx = amy_tx.merge(gencode, on="Transcript_ID_clean", how="left")
        amy_tx = amy_tx[amy_tx["biotype"].eq("protein_coding")].copy()
    else:
        amy_tx["biotype"] = np.nan

    amy_tx = amy_tx.merge(expr_ann, on="Transcript_ID_clean", how="left")
    amy_tx["is_strict_amyloid"] = amy_tx["Consensus_collapsed"].eq("Amyloid")
    return amy_tx


def build_sample_burden(root, selected_samples, sample_meta, amy_tx):
    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    sample_idx = [zarr_samples.index(s) for s in selected_samples]
    transcript_ids = pd.Series(root["layers/expression/transcript_ids"][:], name="Transcript_ID").astype(str)
    tpm = pd.DataFrame(
        root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=transcript_ids,
        columns=selected_samples,
    )
    expr_long = tpm.reset_index().melt(id_vars="Transcript_ID", var_name="sample", value_name="TPM")
    expr_long["Transcript_ID_clean"] = expr_long["Transcript_ID"].map(clean_enst)
    expr_long = expr_long.merge(sample_meta[["sample", "subject", "time", "raw_Group", "Group"]], on="sample", how="left")
    expr_amy = expr_long.merge(amy_tx, on=["sample", "Transcript_ID_clean"], how="left")
    expr_amy["has_amyloid_annotation"] = expr_amy["Amyloid_Index_max"].notna()
    expr_amy["strict_amyloid_TPM_component"] = np.where(
        expr_amy["Consensus_collapsed"].eq("Amyloid"), expr_amy["TPM"].fillna(0), 0.0
    )
    expr_amy["strict_weighted_index_component"] = np.where(
        expr_amy["Consensus_collapsed"].eq("Amyloid"),
        expr_amy["TPM"].fillna(0) * expr_amy["Amyloid_Index_max"],
        0.0,
    )
    expr_amy["continuous_weighted_index_component"] = expr_amy["TPM"].fillna(0) * expr_amy["Amyloid_Index_max"]

    def sample_summary(g):
        total_tpm = g["TPM"].sum()
        ann = g[g["has_amyloid_annotation"]]
        ann_tpm = ann["TPM"].sum()
        strict_tpm = ann["strict_amyloid_TPM_component"].sum()
        strict_weighted = ann["strict_weighted_index_component"].sum()
        continuous_weighted = ann["continuous_weighted_index_component"].sum()
        return pd.Series(
            {
                "total_TPM": total_tpm,
                "protein_coding_annotated_TPM": ann_tpm,
                "protein_coding_annotation_fraction": ann_tpm / total_tpm if total_tpm else np.nan,
                "strict_amyloid_TPM": strict_tpm,
                "strict_amyloid_TPM_fraction": strict_tpm / ann_tpm if ann_tpm else np.nan,
                "n_strict_amyloid_expressed_TPM_ge_1": (
                    ann["Consensus_collapsed"].eq("Amyloid") & (ann["TPM"] >= 1)
                ).sum(),
                "strict_weighted_index_sum": strict_weighted,
                "continuous_weighted_index_sum": continuous_weighted,
                "continuous_expression_weighted_amyloid_index": continuous_weighted / ann_tpm if ann_tpm else np.nan,
            }
        )

    burden = expr_amy.groupby("sample").apply(sample_summary).reset_index()
    burden = burden.merge(sample_meta, on="sample", how="left")
    burden["keep_coverage_qc"] = burden["protein_coding_annotation_fraction"] >= COVERAGE_MIN

    strict = expr_amy[expr_amy["Consensus_collapsed"].eq("Amyloid")].copy()
    strict["strict_contribution"] = strict["TPM"].fillna(0) * strict["Amyloid_Index_max"]
    strict["strict_tpm_contribution"] = strict["TPM"].fillna(0)
    return burden, strict


def leave_one_subject_out(burden, analysis_set):
    rows = []
    subjects = burden["subject"].dropna().astype(str).unique()
    for subject in subjects:
        subset = burden[burden["subject"].astype(str).ne(subject)]
        left = burden[burden["subject"].astype(str).eq(subject)]
        left_samples = ",".join(left["sample"].astype(str))
        left_groups = ",".join(left["Group"].astype(str).unique())
        for metric in METRICS:
            d = subset.dropna(subset=[metric, "Group"])
            ac = d.loc[d["Group"].eq("AC"), metric].astype(float)
            tol = d.loc[d["Group"].eq("TOL"), metric].astype(float)
            delta = tol.mean() - ac.mean()
            rows.append(
                {
                    "analysis_set": analysis_set,
                    "metric": metric,
                    "left_out_subject": subject,
                    "left_out_sample": left_samples,
                    "left_out_group": left_groups,
                    "n_AC": len(ac),
                    "n_TOL": len(tol),
                    "mean_AC": ac.mean(),
                    "mean_TOL": tol.mean(),
                    "delta_mean_TOL_minus_AC": delta,
                    "direction_chili_gt_control": bool(delta > 0),
                }
            )
    return pd.DataFrame(rows)


def write_report(sample_burden, stats, stability, stability_summary, top_panel):
    primary = stats[(stats["analysis_set"] == "coverage_QC_protein_coding_annotation_fraction_ge_0_5") & (stats["metric"] == PRIMARY_METRIC)].iloc[0]
    primary_all = stats[(stats["analysis_set"] == "all_BC_samples") & (stats["metric"] == PRIMARY_METRIC)].iloc[0]
    qc_excluded = sample_burden.loc[~sample_burden["keep_coverage_qc"], ["sample", "subject", "Group", "raw_Group", "protein_coding_annotation_fraction"]]
    qc_stab = stability_summary[
        (stability_summary["analysis_set"] == "coverage_QC_protein_coding_annotation_fraction_ge_0_5")
        & (stability_summary["metric"] == PRIMARY_METRIC)
    ].iloc[0]
    all_stab = stability_summary[
        (stability_summary["analysis_set"] == "all_BC_samples") & (stability_summary["metric"] == PRIMARY_METRIC)
    ].iloc[0]

    ac_samples = ", ".join(sample_burden.loc[sample_burden["Group"].eq("AC"), "sample"].astype(str))
    tol_samples = ", ".join(sample_burden.loc[sample_burden["Group"].eq("TOL"), "sample"].astype(str))
    qc_ac_samples = ", ".join(
        sample_burden.loc[sample_burden["Group"].eq("AC") & sample_burden["keep_coverage_qc"], "sample"].astype(str)
    )
    qc_tol_samples = ", ".join(
        sample_burden.loc[sample_burden["Group"].eq("TOL") & sample_burden["keep_coverage_qc"], "sample"].astype(str)
    )

    secondary = stats[stats["analysis_set"].eq("coverage_QC_protein_coding_annotation_fraction_ge_0_5")].copy()
    secondary = secondary.merge(
        stability_summary[["analysis_set", "metric", "n_supports_hypothesis", "n_leave_one_subject_out", "percent_chili_gt_control"]],
        on=["analysis_set", "metric"],
        how="left",
    )
    secondary_rows = []
    for _, row in secondary.iterrows():
        secondary_rows.append(
            f"| `{row['metric']}` | {row['delta_mean_TOL_minus_AC']:.6g} | {row['hedges_g_TOL_minus_AC']:.3g} | "
            f"{row['empirical_permutation_p_one_sided_TOL_gt_AC']:.4g} | "
            f"{int(row['n_supports_hypothesis'])}/{int(row['n_leave_one_subject_out'])} ({row['percent_chili_gt_control']:.1f}%) |"
        )
    secondary_table = "\n".join(secondary_rows)

    report = f"""# Leave-One-Subject-Out Stability Map

## What Was Done

Recomputed baseline (`time == "BC"`) expression-weighted amyloidogenic burden directly from `project.zarr` and `GSE287540_SraRunTable.csv`.
Groups were mapped as `control -> AC` and `chili -> TOL`. Protein H1/H2 variants were collapsed to `sample x transcript` before joining transcript-level TPM, preventing expression double-counting. Strict amyloid label was `Consensus_collapsed == "Amyloid"`; continuous score was `Amyloid_Index_max`.

Primary endpoint: `{PRIMARY_METRIC}`.

## Samples

- All BC samples: {int((sample_burden["Group"] == "AC").sum())} AC/control and {int((sample_burden["Group"] == "TOL").sum())} TOL/chili.
- Coverage-QC subset (`protein_coding_annotation_fraction >= {COVERAGE_MIN}`): {int(((sample_burden["Group"] == "AC") & sample_burden["keep_coverage_qc"]).sum())} AC/control and {int(((sample_burden["Group"] == "TOL") & sample_burden["keep_coverage_qc"]).sum())} TOL/chili.
- QC-excluded samples: {", ".join(qc_excluded["sample"].astype(str).tolist()) if len(qc_excluded) else "none"}.

All BC AC/control samples: {ac_samples}.

All BC TOL/chili samples: {tol_samples}.

Coverage-QC AC/control samples: {qc_ac_samples}.

Coverage-QC TOL/chili samples: {qc_tol_samples}.

## Main Result

Primary endpoint, all BC samples:

- AC/control mean: {primary_all["mean_AC"]:.6g}
- TOL/chili mean: {primary_all["mean_TOL"]:.6g}
- Delta mean TOL - AC: {primary_all["delta_mean_TOL_minus_AC"]:.6g}
- Hedges g: {primary_all["hedges_g_TOL_minus_AC"]:.3g}
- Exact empirical one-sided p (`TOL > AC`): {primary_all["empirical_permutation_p_one_sided_TOL_gt_AC"]:.4g}
- Bootstrap 95% CI for delta: [{primary_all["bootstrap_delta_mean_TOL_minus_AC_ci95_low"]:.6g}, {primary_all["bootstrap_delta_mean_TOL_minus_AC_ci95_high"]:.6g}]
- LOSO direction stability: {int(all_stab["n_supports_hypothesis"])}/{int(all_stab["n_leave_one_subject_out"])} ({all_stab["percent_chili_gt_control"]:.1f}%)

Primary endpoint, coverage-QC subset:

- AC/control mean: {primary["mean_AC"]:.6g}
- TOL/chili mean: {primary["mean_TOL"]:.6g}
- Delta mean TOL - AC: {primary["delta_mean_TOL_minus_AC"]:.6g}
- Hedges g: {primary["hedges_g_TOL_minus_AC"]:.3g}
- Exact empirical one-sided p (`TOL > AC`): {primary["empirical_permutation_p_one_sided_TOL_gt_AC"]:.4g}
- Bootstrap 95% CI for delta: [{primary["bootstrap_delta_mean_TOL_minus_AC_ci95_low"]:.6g}, {primary["bootstrap_delta_mean_TOL_minus_AC_ci95_high"]:.6g}]
- LOSO direction stability: {int(qc_stab["n_supports_hypothesis"])}/{int(qc_stab["n_leave_one_subject_out"])} ({qc_stab["percent_chili_gt_control"]:.1f}%)

## Secondary Endpoints, Coverage-QC Subset

| Metric | Delta mean TOL - AC | Hedges g | Empirical p, one-sided TOL > AC | LOSO direction |
|---|---:|---:|---:|---:|
{secondary_table}

The top strict amyloid contributor panel was regenerated as `top_strict_amyloid_contributor_panel.tsv`.

## Interpretation

For the primary endpoint, the direction is higher in ChILI/TOL than control/AC in both all-sample and coverage-QC analyses. However, exact permutation p-values and bootstrap confidence intervals remain weak/wide, consistent with a pilot discovery signal rather than a validated biomarker. LOSO stability shows whether the direction survives removal of each subject and is the preferred robustness readout for this small cohort.

## Limitations

- Small sample size: 6 vs 5 all BC; 6 vs 4 after coverage QC.
- `SRR32060276` has very low protein/amyloid annotation coverage and inflates uncertainty in all-sample analyses.
- Amyloid labels are computational predictions and need orthogonal validation.
- This is discovery/exploratory stability analysis; no independent validation cohort is included here.

## Next Best Step

Validate the primary endpoint and the top strict amyloid contributor panel in an independent baseline cohort or by a prespecified cross-validation/permutation framework, then re-score top contributors with slower specialized aggregation predictors.
"""
    (OUT / "LEAVE_ONE_SUBJECT_OUT_STABILITY_MAP_REPORT.md").write_text(report, encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    selected_samples, sample_meta = load_sample_meta(root)
    amy_tx = build_collapsed_annotations(root, selected_samples)
    sample_burden, strict = build_sample_burden(root, selected_samples, sample_meta, amy_tx)

    sample_meta.to_csv(OUT / "sample_metadata_BC.tsv", sep="\t", index=False)
    amy_tx.to_csv(OUT / "collapsed_sample_transcript_amyloid_annotations.tsv", sep="\t", index=False)
    sample_burden.to_csv(OUT / "sample_level_amyloid_burden.tsv", sep="\t", index=False)

    analysis_sets = {
        "all_BC_samples": sample_burden,
        "coverage_QC_protein_coding_annotation_fraction_ge_0_5": sample_burden[sample_burden["keep_coverage_qc"]].copy(),
    }

    stats_rows = []
    stability_frames = []
    for label, df in analysis_sets.items():
        for metric in METRICS:
            stats_rows.append(summarize_metric(df, metric, label, rng))
        stability_frames.append(leave_one_subject_out(df, label))
    stats = pd.DataFrame(stats_rows)
    stability = pd.concat(stability_frames, ignore_index=True)
    stability_summary = (
        stability.groupby(["analysis_set", "metric"], as_index=False)
        .agg(
            n_leave_one_subject_out=("direction_chili_gt_control", "size"),
            n_supports_hypothesis=("direction_chili_gt_control", "sum"),
            min_delta_mean_TOL_minus_AC=("delta_mean_TOL_minus_AC", "min"),
            median_delta_mean_TOL_minus_AC=("delta_mean_TOL_minus_AC", "median"),
            max_delta_mean_TOL_minus_AC=("delta_mean_TOL_minus_AC", "max"),
        )
    )
    stability_summary["percent_chili_gt_control"] = (
        100 * stability_summary["n_supports_hypothesis"] / stability_summary["n_leave_one_subject_out"]
    )

    stats.to_csv(OUT / "group_statistics.tsv", sep="\t", index=False)
    stability.to_csv(OUT / "leave_one_subject_out_stability.tsv", sep="\t", index=False)
    stability_summary.to_csv(OUT / "leave_one_subject_out_stability_summary.tsv", sep="\t", index=False)

    top_rows = []
    for label, df in analysis_sets.items():
        strict_set = strict[strict["sample"].isin(df["sample"])].copy()
        panel = (
            strict_set.groupby("Transcript_ID_clean", as_index=False)
            .agg(
                gene_name=("gene_name", "first"),
                gene_id=("gene_id", "first"),
                biotype=("biotype", "first"),
                mean_TPM_AC=("TPM", lambda x: x[strict_set.loc[x.index, "Group"].eq("AC")].mean()),
                mean_TPM_TOL=("TPM", lambda x: x[strict_set.loc[x.index, "Group"].eq("TOL")].mean()),
                mean_strict_contribution_AC=("strict_contribution", lambda x: x[strict_set.loc[x.index, "Group"].eq("AC")].mean()),
                mean_strict_contribution_TOL=("strict_contribution", lambda x: x[strict_set.loc[x.index, "Group"].eq("TOL")].mean()),
                Amyloid_Index_max=("Amyloid_Index_max", "max"),
                n_samples_observed=("sample", "nunique"),
            )
        )
        panel["analysis_set"] = label
        panel["delta_mean_strict_contribution_TOL_minus_AC"] = (
            panel["mean_strict_contribution_TOL"] - panel["mean_strict_contribution_AC"]
        )
        panel["delta_mean_TPM_TOL_minus_AC"] = panel["mean_TPM_TOL"] - panel["mean_TPM_AC"]
        top_rows.append(panel.sort_values("delta_mean_strict_contribution_TOL_minus_AC", ascending=False).head(50))
    top_panel = pd.concat(top_rows, ignore_index=True)
    top_panel.to_csv(OUT / "top_strict_amyloid_contributor_panel.tsv", sep="\t", index=False)

    summary = {
        "method": "Leave-One-Subject-Out Stability Map",
        "output_dir": str(OUT),
        "selected_samples": selected_samples,
        "group_counts_all": sample_burden["Group"].value_counts().to_dict(),
        "group_counts_coverage_qc": sample_burden[sample_burden["keep_coverage_qc"]]["Group"].value_counts().to_dict(),
        "coverage_min": COVERAGE_MIN,
        "qc_excluded_samples": sample_burden.loc[
            ~sample_burden["keep_coverage_qc"],
            ["sample", "subject", "Group", "raw_Group", "protein_coding_annotation_fraction"],
        ].to_dict("records"),
        "primary_metric": PRIMARY_METRIC,
        "primary_group_statistics": stats[stats["metric"].eq(PRIMARY_METRIC)].to_dict("records"),
        "primary_loso_stability": stability_summary[stability_summary["metric"].eq(PRIMARY_METRIC)].to_dict("records"),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    write_report(sample_burden, stats, stability, stability_summary, top_panel)

    print("DONE")
    print(OUT)
    print(stats[stats["metric"].eq(PRIMARY_METRIC)].to_string(index=False))
    print(stability_summary[stability_summary["metric"].eq(PRIMARY_METRIC)].to_string(index=False))


if __name__ == "__main__":
    main()
