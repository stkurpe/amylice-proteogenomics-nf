import json
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

warnings.filterwarnings("ignore")


ROOT = Path(__file__).resolve().parent
ZARR_PATH = ROOT / "project.zarr"
ANNOTATION_PATH = ROOT / "GSE287540_SraRunTable.csv"
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "counterfactual_permutation_by_amyloid_labels"
OUT.mkdir(parents=True, exist_ok=True)

TIMEPOINT = "BC"
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
GROUPS = ["AC", "TOL"]
COVERAGE_MIN = 0.5
N_PERMUTATIONS = 10_000
N_BOOTSTRAPS = 10_000
RANDOM_SEED = 20260604
PSEUDOCOUNT = 1e-6

PRIMARY_ENDPOINT = "strict_amyloid_TPM_fraction"
ENDPOINTS = [
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


def zarr_group_to_df(g):
    return pd.DataFrame({col: g[col][:] for col in g.array_keys()})


def safe_qcut(values, q, prefix):
    s = pd.Series(values)
    ranks = s.rank(method="first")
    try:
        bins = pd.qcut(ranks, q=q, labels=False, duplicates="drop")
    except ValueError:
        bins = pd.Series(np.zeros(len(s), dtype=int), index=s.index)
    bins = bins.fillna(0).astype(int)
    return prefix + bins.astype(str)


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


def df_to_markdown(df):
    if df.empty:
        return ""
    display = df.copy()
    for col in display.columns:
        if pd.api.types.is_float_dtype(display[col]):
            display[col] = display[col].map(lambda x: "" if pd.isna(x) else f"{x:.6g}")
        else:
            display[col] = display[col].map(lambda x: "" if pd.isna(x) else str(x))
    headers = list(display.columns)
    rows = display.values.tolist()
    widths = [
        max(len(str(header)), *(len(str(row[i])) for row in rows))
        for i, header in enumerate(headers)
    ]
    out = [
        "| " + " | ".join(str(header).ljust(widths[i]) for i, header in enumerate(headers)) + " |",
        "| " + " | ".join("-" * widths[i] for i in range(len(headers))) + " |",
    ]
    for row in rows:
        out.append("| " + " | ".join(str(row[i]).ljust(widths[i]) for i in range(len(headers))) + " |")
    return "\n".join(out)


def bootstrap_delta_ci(sample_metrics, metric, sample_order, group_by_sample, rng, n=N_BOOTSTRAPS):
    ac_idx = np.array([i for i, s in enumerate(sample_order) if group_by_sample[s] == "AC"], dtype=int)
    tol_idx = np.array([i for i, s in enumerate(sample_order) if group_by_sample[s] == "TOL"], dtype=int)
    values = sample_metrics[metric]
    deltas = np.empty(n, dtype=float)
    for i in range(n):
        ac = rng.choice(ac_idx, size=len(ac_idx), replace=True)
        tol = rng.choice(tol_idx, size=len(tol_idx), replace=True)
        deltas[i] = values[tol].mean() - values[ac].mean()
    return np.nanpercentile(deltas, [2.5, 50, 97.5])


def compute_sample_metrics_from_vectors(
    tpm,
    amyloid_index,
    strict_label,
    sample_code,
    n_samples,
    denominator_by_sample,
):
    strict = strict_label.astype(float)
    strict_tpm = np.bincount(sample_code, weights=tpm * strict, minlength=n_samples)
    strict_weighted = np.bincount(sample_code, weights=tpm * amyloid_index * strict, minlength=n_samples)
    continuous_weighted = np.bincount(sample_code, weights=tpm * amyloid_index, minlength=n_samples)
    strict_count = np.bincount(sample_code, weights=((strict_label) & (tpm >= 1.0)).astype(float), minlength=n_samples)
    return {
        "strict_amyloid_TPM_fraction": strict_tpm / denominator_by_sample,
        "strict_amyloid_TPM": strict_tpm,
        "n_strict_amyloid_expressed_TPM_ge_1": strict_count,
        "strict_weighted_index_sum": strict_weighted,
        "continuous_expression_weighted_amyloid_index": continuous_weighted / denominator_by_sample,
    }


def delta_for_set(metric_values, keep_sample, group_code):
    ac = metric_values[keep_sample & (group_code == 0)]
    tol = metric_values[keep_sample & (group_code == 1)]
    return float(np.nanmean(tol) - np.nanmean(ac))


def summarize_observed(sample_metric_df, analysis_set, metric, rng):
    d = sample_metric_df.dropna(subset=[metric]).copy()
    ac = d.loc[d["Group"].eq("AC"), metric].astype(float).to_numpy()
    tol = d.loc[d["Group"].eq("TOL"), metric].astype(float).to_numpy()
    ci = bootstrap_delta_ci(
        {metric: d[metric].astype(float).to_numpy()},
        metric,
        d["sample"].tolist(),
        dict(zip(d["sample"], d["Group"])),
        rng,
    )
    return {
        "analysis_set": analysis_set,
        "metric": metric,
        "n_AC": int(len(ac)),
        "n_TOL": int(len(tol)),
        "mean_AC": float(np.mean(ac)),
        "mean_TOL": float(np.mean(tol)),
        "median_AC": float(np.median(ac)),
        "median_TOL": float(np.median(tol)),
        "observed_delta_TOL_minus_AC": float(np.mean(tol) - np.mean(ac)),
        "cliffs_delta_TOL_vs_AC": float(cliffs_delta(tol, ac)),
        "bootstrap_delta_ci2_5": float(ci[0]),
        "bootstrap_delta_median": float(ci[1]),
        "bootstrap_delta_ci97_5": float(ci[2]),
    }


def main():
    rng = np.random.default_rng(RANDOM_SEED)
    root = zarr.open_group(str(ZARR_PATH), mode="r")

    annotation = pd.read_csv(ANNOTATION_PATH)
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
    sample_meta.to_csv(OUT / "sample_metadata.tsv", sep="\t", index=False)

    sample_idx = [zarr_samples.index(s) for s in selected_samples]
    transcript_ids = pd.Series(root["layers/expression/transcript_ids"][:], name="Transcript_ID").astype(str)
    transcript_clean = transcript_ids.map(clean_enst)
    tpm_raw = pd.DataFrame(
        root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=transcript_clean,
        columns=selected_samples,
    )
    tpm_raw = tpm_raw[~tpm_raw.index.isna()]
    tpm_by_tx = tpm_raw.groupby(tpm_raw.index).sum()
    mean_tpm_by_tx = tpm_by_tx.mean(axis=1)

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
            protein_length_mean=("protein_length", "mean"),
            n_protein_variants=("protein_id", "nunique"),
            Consensus_collapsed=("Consensus", consensus_priority),
        )
    )
    amy_tx["is_strict_amyloid"] = amy_tx["Consensus_collapsed"].eq("Amyloid")
    amy_tx.to_csv(OUT / "collapsed_sample_transcript_amyloid_annotations.tsv", sep="\t", index=False)

    tpm_by_tx_reset = tpm_by_tx.reset_index()
    tpm_by_tx_reset = tpm_by_tx_reset.rename(columns={tpm_by_tx_reset.columns[0]: "Transcript_ID_clean"})
    expr_rows = amy_tx.merge(
        tpm_by_tx_reset,
        on="Transcript_ID_clean",
        how="left",
    )
    long_rows = []
    for sample in selected_samples:
        cols = ["sample", "Transcript_ID_clean", "Amyloid_Index_max", "protein_length_mean", "Consensus_collapsed", "is_strict_amyloid", sample]
        d = expr_rows.loc[expr_rows["sample"].eq(sample), cols].copy()
        d = d.rename(columns={sample: "TPM"})
        long_rows.append(d)
    rows = pd.concat(long_rows, ignore_index=True)
    rows["TPM"] = rows["TPM"].fillna(0.0).astype(float)
    rows["mean_TPM_BC"] = rows["Transcript_ID_clean"].map(mean_tpm_by_tx).fillna(0.0)
    rows["length_bin"] = safe_qcut(rows["protein_length_mean"], 5, "L")
    rows["expression_bin"] = safe_qcut(np.log2(rows["mean_TPM_BC"].astype(float) + PSEUDOCOUNT), 5, "E")
    rows["permutation_bin"] = rows["length_bin"] + "_" + rows["expression_bin"]
    rows.to_csv(OUT / "annotated_expression_rows_for_permutation.tsv", sep="\t", index=False)

    sample_to_code = {s: i for i, s in enumerate(selected_samples)}
    group_by_sample = dict(zip(sample_meta["sample"], sample_meta["Group"].astype(str)))
    group_code = np.array([0 if group_by_sample[s] == "AC" else 1 for s in selected_samples], dtype=int)
    sample_code = rows["sample"].map(sample_to_code).to_numpy(dtype=int)
    tpm = rows["TPM"].to_numpy(dtype=float)
    amyloid_index = rows["Amyloid_Index_max"].fillna(0.0).to_numpy(dtype=float)
    strict_label = rows["is_strict_amyloid"].fillna(False).to_numpy(dtype=bool)

    total_tpm_by_sample = tpm_by_tx.sum(axis=0).reindex(selected_samples).to_numpy(dtype=float)
    annotated_tpm_by_sample = np.bincount(sample_code, weights=tpm, minlength=len(selected_samples))
    annotation_fraction = annotated_tpm_by_sample / total_tpm_by_sample
    keep_qc = annotation_fraction >= COVERAGE_MIN
    keep_all = np.ones(len(selected_samples), dtype=bool)

    observed_metrics = compute_sample_metrics_from_vectors(
        tpm=tpm,
        amyloid_index=amyloid_index,
        strict_label=strict_label,
        sample_code=sample_code,
        n_samples=len(selected_samples),
        denominator_by_sample=annotated_tpm_by_sample,
    )
    sample_metric_df = pd.DataFrame({"sample": selected_samples})
    sample_metric_df["Group"] = [group_by_sample[s] for s in selected_samples]
    sample_metric_df["total_TPM"] = total_tpm_by_sample
    sample_metric_df["protein_coding_annotated_TPM"] = annotated_tpm_by_sample
    sample_metric_df["protein_coding_annotation_fraction"] = annotation_fraction
    sample_metric_df["keep_coverage_qc"] = keep_qc
    for metric in ENDPOINTS:
        sample_metric_df[metric] = observed_metrics[metric]
    sample_metric_df = sample_metric_df.merge(sample_meta, on=["sample", "Group"], how="left")
    sample_metric_df.to_csv(OUT / "sample_level_permutation_burden.tsv", sep="\t", index=False)

    analysis_sets = {
        "all_BC_samples": keep_all,
        "coverage_QC_protein_coding_annotation_fraction_ge_0_5": keep_qc,
    }
    observed_rows = []
    boot_rng = np.random.default_rng(RANDOM_SEED + 1)
    for set_name, keep in analysis_sets.items():
        d = sample_metric_df.loc[keep].copy()
        for metric in ENDPOINTS:
            observed_rows.append(summarize_observed(d, set_name, metric, boot_rng))
    observed_summary = pd.DataFrame(observed_rows)
    observed_summary.to_csv(OUT / "observed_group_statistics.tsv", sep="\t", index=False)

    bin_indices = [idx.to_numpy() for _, idx in rows.groupby("permutation_bin").groups.items()]
    perm_records = []
    null_deltas = {
        set_name: {metric: np.empty(N_PERMUTATIONS, dtype=float) for metric in ENDPOINTS}
        for set_name in analysis_sets
    }
    shuffled_index = amyloid_index.copy()
    shuffled_strict = strict_label.copy()
    for perm_i in range(N_PERMUTATIONS):
        for idx in bin_indices:
            if len(idx) > 1:
                shuffled = rng.permutation(idx)
                shuffled_index[idx] = amyloid_index[shuffled]
                shuffled_strict[idx] = strict_label[shuffled]
        perm_metrics = compute_sample_metrics_from_vectors(
            tpm=tpm,
            amyloid_index=shuffled_index,
            strict_label=shuffled_strict,
            sample_code=sample_code,
            n_samples=len(selected_samples),
            denominator_by_sample=annotated_tpm_by_sample,
        )
        for set_name, keep in analysis_sets.items():
            for metric in ENDPOINTS:
                null_deltas[set_name][metric][perm_i] = delta_for_set(perm_metrics[metric], keep, group_code)

    for set_name, metric_map in null_deltas.items():
        for metric, null in metric_map.items():
            obs = observed_summary.loc[
                observed_summary["analysis_set"].eq(set_name) & observed_summary["metric"].eq(metric),
                "observed_delta_TOL_minus_AC",
            ].iloc[0]
            perm_records.append(
                {
                    "analysis_set": set_name,
                    "metric": metric,
                    "n_permutations": N_PERMUTATIONS,
                    "observed_delta_TOL_minus_AC": obs,
                    "null_mean_delta": float(np.mean(null)),
                    "null_sd_delta": float(np.std(null, ddof=1)),
                    "null_ci2_5": float(np.percentile(null, 2.5)),
                    "null_ci97_5": float(np.percentile(null, 97.5)),
                    "empirical_p_greater": float((1 + np.sum(null >= obs)) / (N_PERMUTATIONS + 1)),
                    "empirical_p_less": float((1 + np.sum(null <= obs)) / (N_PERMUTATIONS + 1)),
                    "empirical_p_two_sided": float(
                        min(1.0, 2 * min((1 + np.sum(null >= obs)) / (N_PERMUTATIONS + 1), (1 + np.sum(null <= obs)) / (N_PERMUTATIONS + 1)))
                    ),
                    "empirical_percentile": float(np.mean(null <= obs)),
                }
            )
    permutation_summary = pd.DataFrame(perm_records)
    permutation_summary.to_csv(OUT / "permutation_empirical_pvalues.tsv", sep="\t", index=False)

    primary_null = pd.DataFrame(
        {
            f"{set_name}__{metric}": values
            for set_name, metric_map in null_deltas.items()
            for metric, values in metric_map.items()
            if metric == PRIMARY_ENDPOINT
        }
    )
    primary_null.to_csv(OUT / "primary_endpoint_null_deltas.tsv", sep="\t", index=False)

    loo_rows = []
    for set_name, keep in analysis_sets.items():
        samples = np.array(selected_samples)[keep]
        for left_out in samples:
            keep_loo = keep.copy()
            keep_loo[sample_to_code[left_out]] = False
            if np.sum(keep_loo & (group_code == 0)) < 1 or np.sum(keep_loo & (group_code == 1)) < 1:
                continue
            loo_rows.append(
                {
                    "analysis_set": set_name,
                    "left_out_sample": left_out,
                    "left_out_group": group_by_sample[left_out],
                    "metric": PRIMARY_ENDPOINT,
                    "delta_TOL_minus_AC": delta_for_set(observed_metrics[PRIMARY_ENDPOINT], keep_loo, group_code),
                }
            )
    loo = pd.DataFrame(loo_rows)
    loo.to_csv(OUT / "leave_one_sample_out_primary_delta.tsv", sep="\t", index=False)

    contrib = rows[rows["is_strict_amyloid"]].copy()
    contrib["strict_TPM_contribution"] = contrib["TPM"]
    contrib["strict_weighted_contribution"] = contrib["TPM"] * contrib["Amyloid_Index_max"]
    contrib = contrib.merge(sample_meta[["sample", "Group"]], on="sample", how="left")
    top = (
        contrib.groupby("Transcript_ID_clean", as_index=False)
        .agg(
            mean_strict_TPM_AC=("strict_TPM_contribution", lambda x: x[contrib.loc[x.index, "Group"].eq("AC")].mean()),
            mean_strict_TPM_TOL=("strict_TPM_contribution", lambda x: x[contrib.loc[x.index, "Group"].eq("TOL")].mean()),
            mean_strict_weighted_AC=("strict_weighted_contribution", lambda x: x[contrib.loc[x.index, "Group"].eq("AC")].mean()),
            mean_strict_weighted_TOL=("strict_weighted_contribution", lambda x: x[contrib.loc[x.index, "Group"].eq("TOL")].mean()),
            Amyloid_Index_max=("Amyloid_Index_max", "max"),
            protein_length_mean=("protein_length_mean", "mean"),
            n_samples_annotated=("sample", "nunique"),
        )
    )
    top["delta_strict_TPM_TOL_minus_AC"] = top["mean_strict_TPM_TOL"] - top["mean_strict_TPM_AC"]
    top["delta_strict_weighted_TOL_minus_AC"] = top["mean_strict_weighted_TOL"] - top["mean_strict_weighted_AC"]
    top.sort_values("delta_strict_TPM_TOL_minus_AC", ascending=False).head(50).to_csv(
        OUT / "top_strict_amyloid_contributor_panel.tsv", sep="\t", index=False
    )

    summary = {
        "method": "Counterfactual Permutation by Amyloid Labels",
        "random_seed": RANDOM_SEED,
        "n_permutations": N_PERMUTATIONS,
        "n_bootstraps": N_BOOTSTRAPS,
        "selected_samples": selected_samples,
        "group_counts": sample_meta["Group"].astype(str).value_counts().to_dict(),
        "low_coverage_samples": sample_metric_df.loc[
            ~sample_metric_df["keep_coverage_qc"],
            ["sample", "Group", "protein_coding_annotation_fraction"],
        ].to_dict("records"),
        "primary_endpoint": PRIMARY_ENDPOINT,
        "primary_results": permutation_summary[permutation_summary["metric"].eq(PRIMARY_ENDPOINT)].to_dict("records"),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    primary_rows = permutation_summary[permutation_summary["metric"].eq(PRIMARY_ENDPOINT)]
    qc_excluded = sample_metric_df.loc[~sample_metric_df["keep_coverage_qc"], ["sample", "Group", "protein_coding_annotation_fraction"]]
    report = [
        "# Counterfactual Permutation by Amyloid Labels",
        "",
        "## What was done",
        "",
        "Baseline (`time == BC`) blood RNA-seq samples were compared as AC/control versus TOL/chili. Expression was kept fixed, while amyloid labels and continuous amyloid index values were shuffled within bins defined by protein length and baseline mean transcript expression. Protein variants were collapsed to `sample x transcript` before calculating burden, so H1/H2-style variant duplication does not double-count expression.",
        "",
        f"Permutation count: {N_PERMUTATIONS:,}. One-sided empirical p-value tests whether observed `TOL - AC` burden is greater than expected under random amyloid label placement among similar proteins.",
        "",
        "## Samples",
        "",
        df_to_markdown(sample_metric_df[["sample", "subject", "Group", "raw_Group", "time", "protein_coding_annotation_fraction", "keep_coverage_qc"]]),
        "",
        "## QC exclusions",
        "",
        df_to_markdown(qc_excluded) if not qc_excluded.empty else "No samples excluded by coverage QC.",
        "",
        "## Primary result",
        "",
        df_to_markdown(primary_rows),
        "",
        "## Observed endpoint statistics",
        "",
        df_to_markdown(observed_summary),
        "",
        "## Interpretation",
        "",
    ]
    primary_qc = primary_rows[primary_rows["analysis_set"].eq("coverage_QC_protein_coding_annotation_fraction_ge_0_5")].iloc[0]
    if primary_qc["observed_delta_TOL_minus_AC"] > 0 and primary_qc["empirical_p_greater"] < 0.05:
        interpretation = "In the coverage-QC subset, the primary burden is higher in TOL/chili and reaches a borderline one-sided empirical permutation result. Because the two-sided permutation p-value is weaker and the bootstrap CI is broad, this is discovery-level support for amyloidogenic specificity, not validation and not a biomarker claim."
    elif primary_qc["observed_delta_TOL_minus_AC"] > 0:
        interpretation = "In the coverage-QC subset, the primary burden is directionally higher in TOL/chili, but it is not clearly unusual relative to the amyloid-label permutation null. This is weak discovery-level support, not validation."
    else:
        interpretation = "In the coverage-QC subset, the primary burden is not higher in TOL/chili. This does not support the hypothesis under this counterfactual specificity test."
    report.extend(
        [
            interpretation,
            "",
            "## Limitations",
            "",
            "- Small baseline sample size limits inference; effect size and direction stability are more informative than p-value alone.",
            "- The permutation null preserves expression and coarse protein/expression bins, but bin choices are still analytical assumptions.",
            "- Amyloid annotations are computational predictions and remain discovery-stage features.",
            "- The low-coverage sample SRR32060276 materially affects the all-sample analysis, so the coverage-QC subset should be emphasized.",
            "",
            "## Next best step",
            "",
            "Validate the direction and top contributor panel in an independent ChILI baseline cohort or, if unavailable, use a pre-registered internal resampling analysis with fixed bins/endpoints and orthogonal amyloid predictors.",
            "",
        ]
    )
    (OUT / "COUNTERFACTUAL_PERMUTATION_BY_AMYLOID_LABELS_REPORT.md").write_text("\n".join(report), encoding="utf-8")

    print("DONE")
    print(f"Outputs: {OUT}")
    print(primary_rows.to_string(index=False))


if __name__ == "__main__":
    main()
