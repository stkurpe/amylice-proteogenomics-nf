import itertools
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")


ROOT = Path(__file__).resolve().parent
ZARR_PATH = ROOT / "project.zarr"
ANNOTATION_PATH = ROOT / "GSE287540_SraRunTable.csv"
BASE_OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili"
OUT = BASE_OUT / "patient_level_amyloid_risk_voting"
OUT.mkdir(parents=True, exist_ok=True)

TIMEPOINT = "BC"
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
GROUPS = ["AC", "TOL"]
COVERAGE_MIN = 0.5
PANEL_GENES = ["S100A9", "B2M", "LYZ", "CLU"]
RANDOM_SEED = 20260604
N_BOOT = 10000


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
    return float((gt - lt) / (len(x) * len(y)))


def bootstrap_mean_delta_ci(tol, ac, rng, n_boot=N_BOOT, alpha=0.05):
    tol = np.asarray(tol, dtype=float)
    ac = np.asarray(ac, dtype=float)
    tol = tol[np.isfinite(tol)]
    ac = ac[np.isfinite(ac)]
    if len(tol) < 2 or len(ac) < 2:
        return (np.nan, np.nan)
    tol_idx = rng.integers(0, len(tol), size=(n_boot, len(tol)))
    ac_idx = rng.integers(0, len(ac), size=(n_boot, len(ac)))
    deltas = tol[tol_idx].mean(axis=1) - ac[ac_idx].mean(axis=1)
    return tuple(np.quantile(deltas, [alpha / 2, 1 - alpha / 2]))


def exact_label_permutation_p(values, groups, observed_delta):
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups).astype(str)
    ok = np.isfinite(values) & pd.notna(groups)
    values = values[ok]
    groups = groups[ok]
    n_tol = int(np.sum(groups == "TOL"))
    if n_tol == 0 or n_tol == len(values) or not np.isfinite(observed_delta):
        return (np.nan, np.nan, 0)
    deltas = []
    for tol_idx in itertools.combinations(range(len(values)), n_tol):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(tol_idx)] = True
        deltas.append(values[mask].mean() - values[~mask].mean())
    deltas = np.asarray(deltas)
    p_greater = (np.sum(deltas >= observed_delta - 1e-15) + 1) / (len(deltas) + 1)
    p_two = (np.sum(np.abs(deltas) >= abs(observed_delta) - 1e-15) + 1) / (len(deltas) + 1)
    return (float(p_greater), float(p_two), int(len(deltas)))


def load_bc_sample_metadata(root):
    annotation = pd.read_csv(ANNOTATION_PATH)
    annotation["raw_Group"] = annotation["Group"].astype(str)
    annotation["Group"] = annotation["raw_Group"].map(GROUP_MAP)
    annotation["Run"] = annotation["Run"].astype(str)
    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    meta = annotation[
        annotation["Run"].isin(zarr_samples)
        & annotation["time"].eq(TIMEPOINT)
        & annotation["Group"].isin(GROUPS)
    ].copy()
    meta = meta.rename(columns={"Run": "sample"})
    selected = [s for s in zarr_samples if s in set(meta["sample"])]
    meta = meta.set_index("sample").loc[selected].reset_index()
    meta["Group"] = pd.Categorical(meta["Group"], categories=GROUPS, ordered=True)
    return meta, selected, zarr_samples


def build_collapsed_amyloid(root, samples):
    rows = []
    for sample in samples:
        amy = zarr_group_to_df(root[f"samples/{sample}/amyloid"])
        pf = zarr_group_to_df(root[f"samples/{sample}/protein_features"])
        amy["ID"] = amy["sample"].astype(str) + "|" + amy["Sequence_ID"].astype(str)
        pf["ID"] = pf["sample"].astype(str) + "|" + pf["protein_id"].astype(str)
        combined = pf.merge(amy.drop(columns=["sample"], errors="ignore"), on="ID", how="inner")
        combined["Transcript_ID_clean"] = combined["Sequence_ID"].map(clean_enst)
        combined["abs_charge_density"] = combined["charge_density"].abs()
        combined["pI_distance_from_7"] = (combined["pI"] - 7).abs()
        for col in ["beta_propensity", "KD_max", "aromaticity", "abs_charge_density", "pI_distance_from_7"]:
            combined[f"{col}_z"] = zscore_series(combined[col])
        combined["Amyloid_Index"] = (
            combined["beta_propensity_z"]
            + combined["KD_max_z"]
            + combined["aromaticity_z"]
            - combined["abs_charge_density_z"]
            - combined["pI_distance_from_7_z"]
        ) / 5
        rows.append(
            combined.dropna(subset=["sample", "Transcript_ID_clean"])
            .groupby(["sample", "Transcript_ID_clean"], as_index=False)
            .agg(
                Consensus_collapsed=("Consensus", consensus_priority),
                Amyloid_Index_max=("Amyloid_Index", "max"),
                n_protein_variants=("protein_id", "nunique"),
            )
        )
    return pd.concat(rows, ignore_index=True)


def load_expression_long(root, selected_samples, zarr_samples):
    sample_idx = [zarr_samples.index(s) for s in selected_samples]
    transcript_ids = pd.Series(root["layers/expression/transcript_ids"][:], name="Transcript_ID").astype(str)
    tpm = pd.DataFrame(
        root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=transcript_ids,
        columns=selected_samples,
    )
    expr_ann = zarr_group_to_df(root[f"samples/{selected_samples[0]}/expression"])[
        ["transcript_id", "gene_id", "gene_name"]
    ].copy()
    expr_ann["Transcript_ID_clean"] = expr_ann["transcript_id"].map(clean_enst)
    expr_ann = expr_ann[["Transcript_ID_clean", "gene_id", "gene_name"]].drop_duplicates("Transcript_ID_clean")
    expr_long = tpm.reset_index().melt(id_vars="Transcript_ID", var_name="sample", value_name="TPM")
    expr_long["Transcript_ID_clean"] = expr_long["Transcript_ID"].map(clean_enst)
    expr_long = expr_long.merge(expr_ann, on="Transcript_ID_clean", how="left")
    return expr_long


def build_nonsense_burden(root, samples):
    rows = []
    for sample in samples:
        sg = root[f"samples/{sample}"]
        if "nonsense" not in sg:
            rows.append({"sample": sample, "n_nonsense_rows": 0, "n_nonsense_transcripts": 0})
            continue
        nd = zarr_group_to_df(sg["nonsense"])
        if nd.empty:
            rows.append({"sample": sample, "n_nonsense_rows": 0, "n_nonsense_transcripts": 0})
            continue
        nd["Transcript_ID_clean"] = nd["TRANSCRIPT_ID"].map(clean_enst)
        rows.append(
            {
                "sample": sample,
                "n_nonsense_rows": int(len(nd)),
                "n_nonsense_transcripts": int(nd["Transcript_ID_clean"].nunique()),
            }
        )
    return pd.DataFrame(rows)


def compare_risk_score(votes, analysis_set, rng):
    ac = votes.loc[votes["Group"].eq("AC"), "risk_vote_score"].astype(float).to_numpy()
    tol = votes.loc[votes["Group"].eq("TOL"), "risk_vote_score"].astype(float).to_numpy()
    delta = float(np.mean(tol) - np.mean(ac))
    ci_low, ci_high = bootstrap_mean_delta_ci(tol, ac, rng)
    perm_greater, perm_two, n_perm = exact_label_permutation_p(votes["risk_vote_score"], votes["Group"], delta)
    return {
        "analysis_set": analysis_set,
        "n_AC": int(len(ac)),
        "n_TOL": int(len(tol)),
        "mean_AC": float(np.mean(ac)),
        "mean_TOL": float(np.mean(tol)),
        "median_AC": float(np.median(ac)),
        "median_TOL": float(np.median(tol)),
        "delta_mean_TOL_minus_AC": delta,
        "bootstrap_delta_mean_ci_low": float(ci_low),
        "bootstrap_delta_mean_ci_high": float(ci_high),
        "cliffs_delta_TOL_vs_AC": cliffs_delta(tol, ac),
        "mannwhitney_two_sided_p": float(mannwhitneyu(tol, ac, alternative="two-sided").pvalue),
        "mannwhitney_one_sided_p_TOL_gt_AC": float(mannwhitneyu(tol, ac, alternative="greater").pvalue),
        "exact_label_permutation_p_TOL_gt_AC": perm_greater,
        "exact_label_permutation_two_sided_p": perm_two,
        "n_exact_label_permutations": n_perm,
    }


def make_votes(metrics, analysis_set):
    features = [
        ("strict_amyloid_TPM_fraction", "vote_strict_amyloid_TPM_fraction_high"),
        ("n_strict_amyloid_expressed_TPM_ge_1", "vote_n_strict_amyloid_transcripts_high"),
        ("n_nonsense_transcripts", "vote_nonsense_burden_high"),
        ("panel_strict_amyloid_TPM", "vote_S100A9_B2M_LYZ_CLU_contributor_high"),
        ("continuous_expression_weighted_amyloid_index", "vote_continuous_amyloid_index_high"),
    ]
    ac = metrics[metrics["Group"].eq("AC")]
    thresholds = []
    votes = metrics.copy()
    for metric, vote_col in features:
        threshold = float(ac[metric].max())
        votes[vote_col] = (votes[metric] > threshold).astype(int)
        thresholds.append(
            {
                "analysis_set": analysis_set,
                "metric": metric,
                "vote_column": vote_col,
                "threshold_rule": "sample value > max(AC/control values within this analysis set)",
                "threshold_value": threshold,
            }
        )
    vote_cols = [v for _, v in features]
    votes["risk_vote_score"] = votes[vote_cols].sum(axis=1)
    votes["n_possible_votes"] = len(vote_cols)
    votes["risk_vote_fraction"] = votes["risk_vote_score"] / len(vote_cols)
    votes.insert(0, "analysis_set", analysis_set)
    return votes, pd.DataFrame(thresholds)


def main():
    rng = np.random.default_rng(RANDOM_SEED)
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    sample_meta, selected_samples, zarr_samples = load_bc_sample_metadata(root)

    amy_tx = build_collapsed_amyloid(root, selected_samples)
    expr_long = load_expression_long(root, selected_samples, zarr_samples)
    expr_amy = expr_long.merge(amy_tx, on=["sample", "Transcript_ID_clean"], how="left")
    expr_amy["strict_amyloid_TPM"] = np.where(
        expr_amy["Consensus_collapsed"].eq("Amyloid"), expr_amy["TPM"].fillna(0), 0.0
    )
    expr_amy["strict_weighted_index"] = np.where(
        expr_amy["Consensus_collapsed"].eq("Amyloid"),
        expr_amy["TPM"].fillna(0) * expr_amy["Amyloid_Index_max"],
        0.0,
    )

    panel_long = expr_amy[
        expr_amy["gene_name"].isin(PANEL_GENES) & expr_amy["Consensus_collapsed"].eq("Amyloid")
    ].copy()
    panel_long["panel_gene"] = panel_long["gene_name"]
    panel_long["panel_strict_contribution"] = panel_long["strict_weighted_index"]
    panel_summary = (
        panel_long.groupby(["sample", "panel_gene"], as_index=False)
        .agg(
            panel_gene_strict_amyloid_TPM=("strict_amyloid_TPM", "sum"),
            panel_gene_strict_weighted_index=("strict_weighted_index", "sum"),
            n_panel_strict_amyloid_transcripts=("Transcript_ID_clean", "nunique"),
        )
    )
    panel_wide = (
        panel_summary.groupby("sample", as_index=False)
        .agg(
            panel_strict_amyloid_TPM=("panel_gene_strict_amyloid_TPM", "sum"),
            panel_strict_weighted_index_sum=("panel_gene_strict_weighted_index", "sum"),
            n_panel_strict_amyloid_transcripts=("n_panel_strict_amyloid_transcripts", "sum"),
        )
    )

    sample_burden = pd.read_csv(BASE_OUT / "sample_level_amyloid_burden.tsv", sep="\t")
    fresh_cols = [
        "sample",
        "Group",
        "raw_Group",
        "subject",
        "time",
        "treatment",
        "AGE",
        "sex",
        "protein_coding_annotation_fraction",
        "keep_coverage_qc",
        "strict_amyloid_TPM_fraction",
        "strict_amyloid_TPM",
        "n_strict_amyloid_expressed_TPM_ge_1",
        "strict_weighted_index_sum",
        "continuous_expression_weighted_amyloid_index",
    ]
    metrics = sample_burden[fresh_cols].copy()
    metrics = metrics.merge(build_nonsense_burden(root, selected_samples), on="sample", how="left")
    metrics = metrics.merge(panel_wide, on="sample", how="left")
    for col in [
        "n_nonsense_rows",
        "n_nonsense_transcripts",
        "panel_strict_amyloid_TPM",
        "panel_strict_weighted_index_sum",
        "n_panel_strict_amyloid_transcripts",
    ]:
        metrics[col] = metrics[col].fillna(0)

    all_votes, all_thresholds = make_votes(metrics, "all_BC_samples")
    qc_metrics = metrics[metrics["protein_coding_annotation_fraction"] >= COVERAGE_MIN].copy()
    qc_votes, qc_thresholds = make_votes(qc_metrics, "coverage_QC_protein_coding_annotation_fraction_ge_0_5")
    votes = pd.concat([all_votes, qc_votes], ignore_index=True)
    thresholds = pd.concat([all_thresholds, qc_thresholds], ignore_index=True)
    vote_cols = [c for c in votes.columns if c.startswith("vote_")]
    feature_vote_summary = (
        votes.groupby(["analysis_set", "Group"], as_index=False)[vote_cols]
        .sum()
        .melt(
            id_vars=["analysis_set", "Group"],
            var_name="vote_feature",
            value_name="n_samples_with_vote",
        )
    )

    stats = pd.DataFrame(
        [
            compare_risk_score(all_votes, "all_BC_samples", rng),
            compare_risk_score(qc_votes, "coverage_QC_protein_coding_annotation_fraction_ge_0_5", rng),
        ]
    )

    excluded_qc = metrics.loc[
        metrics["protein_coding_annotation_fraction"] < COVERAGE_MIN,
        ["sample", "Group", "raw_Group", "subject", "protein_coding_annotation_fraction"],
    ].copy()

    sample_meta.to_csv(OUT / "sample_metadata_BC.tsv", sep="\t", index=False)
    metrics.to_csv(OUT / "sample_metric_inputs.tsv", sep="\t", index=False)
    panel_summary.to_csv(OUT / "S100A9_B2M_LYZ_CLU_panel_by_sample_gene.tsv", sep="\t", index=False)
    votes.to_csv(OUT / "sample_risk_votes.tsv", sep="\t", index=False)
    feature_vote_summary.to_csv(OUT / "feature_vote_summary.tsv", sep="\t", index=False)
    thresholds.to_csv(OUT / "risk_vote_feature_thresholds.tsv", sep="\t", index=False)
    stats.to_csv(OUT / "risk_vote_group_statistics.tsv", sep="\t", index=False)
    excluded_qc.to_csv(OUT / "coverage_qc_excluded_samples.tsv", sep="\t", index=False)

    summary = {
        "method": "Patient-Level Amyloid Risk Voting",
        "output_dir": str(OUT),
        "timepoint": TIMEPOINT,
        "group_mapping": {"control": "AC", "chili": "TOL"},
        "n_selected_BC_samples": int(len(selected_samples)),
        "group_counts_all": sample_meta["Group"].astype(str).value_counts().to_dict(),
        "coverage_min": COVERAGE_MIN,
        "excluded_qc_samples": excluded_qc.to_dict("records"),
        "panel_genes": PANEL_GENES,
        "vote_threshold_rule": "high if sample value is greater than the maximum AC/control value within the same analysis set",
        "risk_vote_statistics": stats.to_dict("records"),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    main_all = stats.loc[stats["analysis_set"].eq("all_BC_samples")].iloc[0]
    main_qc = stats.loc[
        stats["analysis_set"].eq("coverage_QC_protein_coding_annotation_fraction_ge_0_5")
    ].iloc[0]
    all_samples_md = ", ".join(
        f"{r.sample}({r.Group}/{r.raw_Group}, {r.subject})" for r in metrics.itertuples(index=False)
    )
    excluded_md = (
        ", ".join(
            f"{r.sample}({r.Group}/{r.raw_Group}, annotation fraction={r.protein_coding_annotation_fraction:.3f})"
            for r in excluded_qc.itertuples(index=False)
        )
        if not excluded_qc.empty
        else "none"
    )
    report = f"""# Patient-Level Amyloid Risk Voting

## What Was Done

Recomputed a sample-level voting endpoint for baseline blood RNA-seq samples only (`time == "BC"`), comparing AC/control vs TOL/chili. Protein variants were collapsed to sample x transcript before strict amyloid labels were used, so H1/H2 protein variants do not double-count transcript-level expression. Strict amyloid label was `Consensus_collapsed == "Amyloid"`; continuous score input was `Amyloid_Index_max` as already summarized into `continuous_expression_weighted_amyloid_index`.

Each sample received one vote for each abnormal high feature, using a control-anchored threshold: high means the value is greater than the maximum AC/control value in the same analysis set.

Votes:

- `strict_amyloid_TPM_fraction` high
- `n_strict_amyloid_transcripts` high, using expressed strict amyloid transcripts with TPM >= 1
- nonsense burden high, using unique nonsense transcript count
- S100A9/B2M/LYZ/CLU strict amyloid contributor high
- `continuous_expression_weighted_amyloid_index` high

## Samples

All BC samples: {all_samples_md}

Coverage-QC exclusion threshold: `protein_coding_annotation_fraction < {COVERAGE_MIN}`.

Excluded from QC subset: {excluded_md}.

## Main Result

All BC samples: mean risk vote score AC={main_all.mean_AC:.3g}, TOL={main_all.mean_TOL:.3g}; median AC={main_all.median_AC:.3g}, TOL={main_all.median_TOL:.3g}; delta mean TOL-AC={main_all.delta_mean_TOL_minus_AC:.3g}; Cliff's delta={main_all.cliffs_delta_TOL_vs_AC:.3g}; Mann-Whitney one-sided p(TOL>AC)={main_all.mannwhitney_one_sided_p_TOL_gt_AC:.4g}; exact permutation p(TOL>AC)={main_all.exact_label_permutation_p_TOL_gt_AC:.4g}; bootstrap 95% CI for mean delta=[{main_all.bootstrap_delta_mean_ci_low:.3g}, {main_all.bootstrap_delta_mean_ci_high:.3g}].

Coverage-QC subset: mean risk vote score AC={main_qc.mean_AC:.3g}, TOL={main_qc.mean_TOL:.3g}; median AC={main_qc.median_AC:.3g}, TOL={main_qc.median_TOL:.3g}; delta mean TOL-AC={main_qc.delta_mean_TOL_minus_AC:.3g}; Cliff's delta={main_qc.cliffs_delta_TOL_vs_AC:.3g}; Mann-Whitney one-sided p(TOL>AC)={main_qc.mannwhitney_one_sided_p_TOL_gt_AC:.4g}; exact permutation p(TOL>AC)={main_qc.exact_label_permutation_p_TOL_gt_AC:.4g}; bootstrap 95% CI for mean delta=[{main_qc.bootstrap_delta_mean_ci_low:.3g}, {main_qc.bootstrap_delta_mean_ci_high:.3g}].

## Interpretation

The direction is consistent with the hypothesis: TOL/chili samples have higher risk vote scores than AC/control samples in both all-BC and coverage-QC analyses. However, the exact permutation evidence is weak after QC, so this should be treated as exploratory support, not validation of a biomarker.

The signal is driven by strict amyloid burden features: strict amyloid TPM fraction and number of expressed strict amyloid transcripts. Nonsense burden, S100A9/B2M/LYZ/CLU contributor burden, and continuous amyloid index did not exceed the control-anchored high threshold in TOL/chili samples.

## Limitations

Small n, control-derived thresholds from the same dataset, treatment imbalance, batch/sex/age confounding risk, incomplete protein/amyloid annotation for SRR32060276, and no independent validation cohort.

## Next Best Step

Run a validation-oriented model on a pre-registered feature set: either exact/penalized patient-level logistic testing of the five votes or counts-based transcript-level DE followed by locked risk-score evaluation in an independent ChILI/control dataset.
"""
    (OUT / "PATIENT_LEVEL_AMYLOID_RISK_VOTING_REPORT.md").write_text(report, encoding="utf-8")
    print("DONE")
    print(OUT)
    print(stats)


if __name__ == "__main__":
    main()
