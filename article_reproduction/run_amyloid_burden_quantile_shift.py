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
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "amyloid_burden_quantile_shift"
GENCODE_CANDIDATES = [
    ROOT
    / "analysis_outputs"
    / "amyloid_bc_control_vs_chili"
    / "rank_based_amyloid_burden_signature"
    / "gencode_v48_transcript_annotation.tsv",
    ROOT
    / "analysis_outputs"
    / "amyloid_bc_control_vs_chili"
    / "competitive_gene_set_test_gencode_v48_protein_coding"
    / "gencode_transcript_annotation_used.tsv",
]

TIMEPOINT = "BC"
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
GROUPS = ["AC", "TOL"]
COVERAGE_MIN = 0.5
N_BOOTSTRAP = 20000
BOOTSTRAP_SEED = 20260604


def clean_enst(x):
    if pd.isna(x):
        return np.nan
    m = re.search(r"ENST\d+(?:\.\d+)?", str(x))
    return m.group(0).split(".")[0] if m else np.nan


def zscore(s):
    s = pd.Series(s, dtype=float)
    sd = s.std(skipna=True, ddof=0)
    if pd.isna(sd) or sd == 0:
        return pd.Series(np.nan, index=s.index)
    return (s - s.mean(skipna=True)) / sd


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


def zarr_group_to_df(g):
    return pd.DataFrame({col: g[col][:] for col in g.array_keys()})


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


def exact_permutation_p(ac, tol):
    ac = np.asarray(ac, dtype=float)
    tol = np.asarray(tol, dtype=float)
    observed = tol.mean() - ac.mean()
    pooled = np.concatenate([ac, tol])
    n_tol = len(tol)
    deltas = []
    for tol_idx in itertools.combinations(range(len(pooled)), n_tol):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(tol_idx)] = True
        deltas.append(pooled[mask].mean() - pooled[~mask].mean())
    deltas = np.asarray(deltas)
    return (
        float(observed),
        float(np.mean(deltas >= observed - 1e-15)),
        float(np.mean(np.abs(deltas) >= abs(observed) - 1e-15)),
    )


def bootstrap_delta_ci(ac, tol, n=N_BOOTSTRAP, seed=BOOTSTRAP_SEED):
    rng = np.random.default_rng(seed)
    ac = np.asarray(ac, dtype=float)
    tol = np.asarray(tol, dtype=float)
    boot = np.empty(n)
    for i in range(n):
        boot[i] = rng.choice(tol, len(tol), replace=True).mean() - rng.choice(ac, len(ac), replace=True).mean()
    return np.quantile(boot, [0.025, 0.5, 0.975])


def compare_metric(sample_metrics, metric, label):
    d = sample_metrics.dropna(subset=[metric, "Group"]).copy()
    ac = d.loc[d["Group"].eq("AC"), metric].astype(float).to_numpy()
    tol = d.loc[d["Group"].eq("TOL"), metric].astype(float).to_numpy()
    observed, emp_one, emp_two = exact_permutation_p(ac, tol)
    ci_low, ci_mid, ci_high = bootstrap_delta_ci(ac, tol)
    return {
        "analysis_set": label,
        "metric": metric,
        "n_AC": int(len(ac)),
        "n_TOL": int(len(tol)),
        "mean_AC": float(np.mean(ac)),
        "mean_TOL": float(np.mean(tol)),
        "median_AC": float(np.median(ac)),
        "median_TOL": float(np.median(tol)),
        "delta_mean_TOL_minus_AC": float(observed),
        "bootstrap_delta_mean_TOL_minus_AC_median": float(ci_mid),
        "bootstrap_delta_mean_TOL_minus_AC_ci95_low": float(ci_low),
        "bootstrap_delta_mean_TOL_minus_AC_ci95_high": float(ci_high),
        "empirical_permutation_p_one_sided_TOL_gt_AC": float(emp_one),
        "empirical_permutation_p_two_sided": float(emp_two),
        "mannwhitney_p_two_sided": float(mannwhitneyu(ac, tol, alternative="two-sided").pvalue),
        "welch_t_p_two_sided": float(ttest_ind(ac, tol, equal_var=False, nan_policy="omit").pvalue),
        "cliffs_delta_TOL_vs_AC": float(cliffs_delta(tol, ac)),
        "direction_supports_hypothesis": bool(observed > 0),
    }


def load_gencode_annotation():
    for path in GENCODE_CANDIDATES:
        if path.exists():
            df = pd.read_csv(path, sep="\t")
            if {"Transcript_ID_clean", "biotype"}.issubset(df.columns):
                return df, path
    raise FileNotFoundError("No local GENCODE transcript annotation cache was found.")


def tail_metrics(values):
    values = pd.Series(values, dtype=float).replace([np.inf, -np.inf], np.nan).dropna()
    if values.empty:
        return {
            "q90": 0.0,
            "q95": 0.0,
            "top_1pct_mass": 0.0,
            "top_1pct_positive_mass": 0.0,
            "n_values": 0,
            "n_positive": 0,
        }
    k = max(1, math.ceil(len(values) * 0.01))
    top = values.sort_values(ascending=False).head(k)
    positive_top = values[values > 0].sort_values(ascending=False).head(k)
    return {
        "q90": float(values.quantile(0.90)),
        "q95": float(values.quantile(0.95)),
        "top_1pct_mass": float(top.sum()),
        "top_1pct_positive_mass": float(positive_top.sum()),
        "n_values": int(len(values)),
        "n_positive": int((values > 0).sum()),
    }


def dataframe_to_markdown(df, floatfmt=".4f"):
    out = df.reset_index()
    headers = [str(c) for c in out.columns]
    rows = []
    for _, row in out.iterrows():
        vals = []
        for val in row:
            if isinstance(val, (float, np.floating)):
                vals.append(format(float(val), floatfmt))
            else:
                vals.append(str(val))
        rows.append(vals)
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    lines.extend("| " + " | ".join(r) + " |" for r in rows)
    return "\n".join(lines)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
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
    sample_meta.to_csv(OUT / "sample_metadata_BC.tsv", sep="\t", index=False)

    sample_idx = [zarr_samples.index(s) for s in selected_samples]
    transcript_ids = pd.Series(root["layers/expression/transcript_ids"][:], name="Transcript_ID").astype(str)
    transcript_clean = transcript_ids.map(clean_enst)
    tpm = pd.DataFrame(
        root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=transcript_clean,
        columns=selected_samples,
    )
    tpm.index.name = "Transcript_ID_clean"
    tpm = tpm.groupby(level=0).sum()

    expr_ann = zarr_group_to_df(root[f"samples/{selected_samples[0]}/expression"])[
        ["transcript_id", "gene_id", "gene_name"]
    ].copy()
    expr_ann["Transcript_ID_clean"] = expr_ann["transcript_id"].map(clean_enst)
    expr_ann = expr_ann[["Transcript_ID_clean", "gene_id", "gene_name"]].drop_duplicates("Transcript_ID_clean")

    gencode, gencode_path = load_gencode_annotation()
    gencode_pc = gencode[gencode["biotype"].eq("protein_coding")].drop_duplicates("Transcript_ID_clean").copy()
    protein_coding_transcripts = set(gencode_pc["Transcript_ID_clean"].dropna()).intersection(set(tpm.index))

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
        combined[f"{col}_z"] = zscore(combined[col])
    combined["Amyloid_Index"] = (
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
            Amyloid_Index_max=("Amyloid_Index", "max"),
            Amyloid_Index_mean=("Amyloid_Index", "mean"),
            AmyloGram_Prob_max=("AmyloGram_Prob", "max"),
            AMYPred_Prob_max=("AMYPred_Prob", "max"),
            n_protein_variants=("protein_id", "nunique"),
            Consensus_collapsed=("Consensus", consensus_priority),
        )
    )
    amy_tx = amy_tx.merge(gencode_pc, on="Transcript_ID_clean", how="inner")
    amy_tx = amy_tx.merge(
        expr_ann.rename(columns={"gene_id": "gene_id_expression", "gene_name": "gene_name_expression"}),
        on="Transcript_ID_clean",
        how="left",
    )
    amy_tx["gene_id"] = amy_tx["gene_id"].fillna(amy_tx["gene_id_expression"])
    amy_tx["gene_name"] = amy_tx["gene_name"].fillna(amy_tx["gene_name_expression"])
    amy_tx = amy_tx.drop(columns=["gene_id_expression", "gene_name_expression"], errors="ignore")
    amy_tx["strict_amyloid"] = amy_tx["Consensus_collapsed"].eq("Amyloid")
    amy_tx.to_csv(OUT / "collapsed_sample_transcript_amyloid_annotations.tsv", sep="\t", index=False)

    rows = []
    long_rows = []
    for sample in selected_samples:
        expr = tpm[[sample]].rename(columns={sample: "TPM"}).reset_index()
        expr = expr[expr["Transcript_ID_clean"].isin(protein_coding_transcripts)].copy()
        ann = amy_tx[amy_tx["sample"].eq(sample)].copy()
        expr = expr.merge(ann, on="Transcript_ID_clean", how="left")
        expr["TPM"] = expr["TPM"].fillna(0).astype(float)
        expr["has_amyloid_annotation"] = expr["Amyloid_Index_max"].notna()
        expr["strict_amyloid"] = expr["Consensus_collapsed"].eq("Amyloid")
        expr["continuous_burden"] = expr["TPM"] * expr["Amyloid_Index_max"]
        expr["strict_burden"] = np.where(expr["strict_amyloid"], expr["continuous_burden"], np.nan)

        protein_coding_tpm = float(expr["TPM"].sum())
        annotated_tpm = float(expr.loc[expr["has_amyloid_annotation"], "TPM"].sum())
        strict = expr[expr["strict_amyloid"]].copy()
        annotated = expr[expr["has_amyloid_annotation"]].copy()

        strict_tail = tail_metrics(strict["strict_burden"])
        continuous_tail = tail_metrics(annotated["continuous_burden"])
        row = {
            "sample": sample,
            "n_gencode_protein_coding_transcripts": int(len(expr)),
            "n_amyloid_annotated_protein_coding_transcripts": int(expr["has_amyloid_annotation"].sum()),
            "n_strict_amyloid_transcripts": int(len(strict)),
            "n_strict_amyloid_expressed_TPM_ge_1": int((strict["TPM"] >= 1).sum()),
            "gencode_protein_coding_TPM": protein_coding_tpm,
            "protein_coding_annotated_TPM": annotated_tpm,
            "protein_coding_annotation_fraction": annotated_tpm / protein_coding_tpm if protein_coding_tpm else np.nan,
            "strict_burden_q90": strict_tail["q90"],
            "strict_burden_q95": strict_tail["q95"],
            "strict_top_1pct_burden_mass": strict_tail["top_1pct_mass"],
            "strict_top_1pct_positive_burden_mass": strict_tail["top_1pct_positive_mass"],
            "n_strict_burden_values": strict_tail["n_values"],
            "n_positive_strict_burden_values": strict_tail["n_positive"],
            "continuous_burden_q90": continuous_tail["q90"],
            "continuous_burden_q95": continuous_tail["q95"],
            "continuous_top_1pct_burden_mass": continuous_tail["top_1pct_mass"],
            "continuous_top_1pct_positive_burden_mass": continuous_tail["top_1pct_positive_mass"],
        }
        rows.append(row)

        strict["sample"] = sample
        strict["strict_burden"] = strict["TPM"] * strict["Amyloid_Index_max"]
        strict["strict_burden_rank_desc"] = strict["strict_burden"].rank(method="min", ascending=False)
        long_rows.append(
            strict[
                [
                    "sample",
                    "Transcript_ID_clean",
                    "gene_id",
                    "gene_name",
                    "TPM",
                    "Amyloid_Index_max",
                    "strict_burden",
                    "strict_burden_rank_desc",
                    "Consensus_collapsed",
                    "n_protein_variants",
                ]
            ]
        )

    sample_metrics = pd.DataFrame(rows).merge(sample_meta, on="sample", how="left")
    sample_metrics["keep_coverage_qc"] = sample_metrics["protein_coding_annotation_fraction"] >= COVERAGE_MIN
    sample_metrics.to_csv(OUT / "sample_quantile_shift_metrics.tsv", sep="\t", index=False)

    strict_long = pd.concat(long_rows, ignore_index=True)
    strict_long = strict_long.merge(sample_meta[["sample", "Group", "raw_Group", "subject", "time", "batch"]], on="sample", how="left")
    strict_long.to_csv(OUT / "strict_amyloid_burden_long.tsv", sep="\t", index=False)

    metrics = [
        "strict_burden_q90",
        "strict_burden_q95",
        "strict_top_1pct_burden_mass",
        "strict_top_1pct_positive_burden_mass",
        "continuous_burden_q90",
        "continuous_burden_q95",
        "continuous_top_1pct_burden_mass",
        "continuous_top_1pct_positive_burden_mass",
    ]
    stat_rows = []
    for label, frame in [
        ("all_BC_samples", sample_metrics),
        ("coverage_qc_protein_annotation_fraction_ge_0.5", sample_metrics[sample_metrics["keep_coverage_qc"]]),
    ]:
        for metric in metrics:
            stat_rows.append(compare_metric(frame, metric, label))
    stats = pd.DataFrame(stat_rows)
    stats.to_csv(OUT / "quantile_shift_group_statistics.tsv", sep="\t", index=False)

    stability_rows = []
    primary_metrics = ["strict_burden_q90", "strict_burden_q95", "strict_top_1pct_burden_mass"]
    for label, frame in [
        ("all_BC_samples", sample_metrics),
        ("coverage_qc_protein_annotation_fraction_ge_0.5", sample_metrics[sample_metrics["keep_coverage_qc"]]),
    ]:
        for sample in frame["sample"]:
            sub = frame[~frame["sample"].eq(sample)]
            if sub["Group"].nunique() == 2 and sub.groupby("Group").size().min() >= 2:
                for metric in primary_metrics:
                    comp = compare_metric(sub, metric, f"{label}_leave_one_out")
                    comp["excluded_sample"] = sample
                    stability_rows.append(comp)
    stability = pd.DataFrame(stability_rows)
    stability.to_csv(OUT / "leave_one_out_stability.tsv", sep="\t", index=False)

    top_panel = (
        strict_long.groupby("Transcript_ID_clean", as_index=False)
        .agg(
            gene_name=("gene_name", "first"),
            gene_id=("gene_id", "first"),
            mean_TPM_AC=("TPM", lambda x: x[strict_long.loc[x.index, "Group"].eq("AC")].mean()),
            mean_TPM_TOL=("TPM", lambda x: x[strict_long.loc[x.index, "Group"].eq("TOL")].mean()),
            mean_strict_burden_AC=("strict_burden", lambda x: x[strict_long.loc[x.index, "Group"].eq("AC")].mean()),
            mean_strict_burden_TOL=("strict_burden", lambda x: x[strict_long.loc[x.index, "Group"].eq("TOL")].mean()),
            n_samples_AC=("sample", lambda x: strict_long.loc[x.index, "sample"][strict_long.loc[x.index, "Group"].eq("AC")].nunique()),
            n_samples_TOL=("sample", lambda x: strict_long.loc[x.index, "sample"][strict_long.loc[x.index, "Group"].eq("TOL")].nunique()),
            Amyloid_Index_max=("Amyloid_Index_max", "max"),
        )
    )
    top_panel["delta_mean_strict_burden_TOL_minus_AC"] = (
        top_panel["mean_strict_burden_TOL"] - top_panel["mean_strict_burden_AC"]
    )
    top_panel["delta_mean_TPM_TOL_minus_AC"] = top_panel["mean_TPM_TOL"] - top_panel["mean_TPM_AC"]
    top_panel = top_panel.sort_values("delta_mean_strict_burden_TOL_minus_AC", ascending=False)
    top_panel.to_csv(OUT / "top_strict_amyloid_quantile_contributor_panel.tsv", sep="\t", index=False)

    qc_excluded = sample_metrics.loc[
        ~sample_metrics["keep_coverage_qc"],
        ["sample", "Group", "raw_Group", "protein_coding_annotation_fraction"],
    ].to_dict("records")
    primary_all = stats[(stats["analysis_set"].eq("all_BC_samples")) & (stats["metric"].eq("strict_burden_q95"))].iloc[0]
    primary_qc = stats[
        (stats["analysis_set"].eq("coverage_qc_protein_annotation_fraction_ge_0.5"))
        & (stats["metric"].eq("strict_burden_q95"))
    ].iloc[0]
    summary = {
        "method": "Amyloid Burden Quantile Shift",
        "timepoint": TIMEPOINT,
        "group_mapping": {"control": "AC", "chili": "TOL"},
        "gencode_annotation_path": str(gencode_path),
        "protein_coding_filter": "GENCODE transcript_type/biotype == protein_coding",
        "variant_collapse": "sample x transcript; Amyloid_Index_max; Consensus_collapsed priority Amyloid > Partial > Discordant > Non-Amyloid",
        "burden_definition": "TPM x Amyloid_Index_max",
        "strict_label": 'Consensus_collapsed == "Amyloid"',
        "primary_tail_metric_for_summary": "strict_burden_q95",
        "n_samples": int(sample_metrics.shape[0]),
        "group_counts_all": sample_metrics["Group"].value_counts().to_dict(),
        "group_counts_qc": sample_metrics[sample_metrics["keep_coverage_qc"]]["Group"].value_counts().to_dict(),
        "coverage_qc_threshold": COVERAGE_MIN,
        "qc_excluded_samples": qc_excluded,
        "primary_all": primary_all.to_dict(),
        "primary_qc": primary_qc.to_dict(),
        "statistics": stats.to_dict("records"),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    direction_table = stats.pivot(index="metric", columns="analysis_set", values="delta_mean_TOL_minus_AC")
    direction_table_md = dataframe_to_markdown(direction_table)
    report = f"""# Amyloid Burden Quantile Shift: BC control/AC vs chili/TOL

## What Was Done

This analysis tests whether the upper tail of expression-weighted amyloidogenic burden is shifted upward in baseline ChILI/TOL blood RNA-seq samples compared with control/AC samples.

Only `time == "BC"` samples were used. Groups were mapped from `GSE287540_SraRunTable.csv` as `control -> AC` and `chili -> TOL`.

The transcript universe was restricted to GENCODE protein-coding transcripts (`biotype == "protein_coding"`) present in the expression layer. Protein variants were collapsed to `sample x transcript` before joining to transcript-level TPM, avoiding H1/H2 double-counting. Strict amyloid label was `Consensus_collapsed == "Amyloid"`. Continuous burden was `TPM x Amyloid_Index_max`.

## Samples

- All BC samples: AC/control n = {int(primary_all['n_AC'])}, TOL/chili n = {int(primary_all['n_TOL'])}
- Coverage-QC subset (`protein_coding_annotation_fraction >= {COVERAGE_MIN}`): AC/control n = {int(primary_qc['n_AC'])}, TOL/chili n = {int(primary_qc['n_TOL'])}
- QC-excluded samples: {', '.join([x['sample'] for x in qc_excluded]) if qc_excluded else 'none'}

## Main Result

Primary summary metric: `strict_burden_q95`, the 95th percentile of `TPM x Amyloid_Index_max` among strict amyloid protein-coding transcripts in each sample.

All BC samples:

- AC/control mean = {primary_all['mean_AC']:.4f}
- TOL/chili mean = {primary_all['mean_TOL']:.4f}
- Delta mean TOL - AC = {primary_all['delta_mean_TOL_minus_AC']:.4f}
- Exact empirical one-sided p (TOL > AC) = {primary_all['empirical_permutation_p_one_sided_TOL_gt_AC']:.4f}
- Exact empirical two-sided p = {primary_all['empirical_permutation_p_two_sided']:.4f}
- Bootstrap 95% CI for delta = [{primary_all['bootstrap_delta_mean_TOL_minus_AC_ci95_low']:.4f}, {primary_all['bootstrap_delta_mean_TOL_minus_AC_ci95_high']:.4f}]
- Cliff's delta TOL vs AC = {primary_all['cliffs_delta_TOL_vs_AC']:.3f}

Coverage-QC subset:

- AC/control mean = {primary_qc['mean_AC']:.4f}
- TOL/chili mean = {primary_qc['mean_TOL']:.4f}
- Delta mean TOL - AC = {primary_qc['delta_mean_TOL_minus_AC']:.4f}
- Exact empirical one-sided p (TOL > AC) = {primary_qc['empirical_permutation_p_one_sided_TOL_gt_AC']:.4f}
- Exact empirical two-sided p = {primary_qc['empirical_permutation_p_two_sided']:.4f}
- Bootstrap 95% CI for delta = [{primary_qc['bootstrap_delta_mean_TOL_minus_AC_ci95_low']:.4f}, {primary_qc['bootstrap_delta_mean_TOL_minus_AC_ci95_high']:.4f}]
- Cliff's delta TOL vs AC = {primary_qc['cliffs_delta_TOL_vs_AC']:.3f}

## Tail Metrics Direction

Delta mean TOL - AC:

{direction_table_md}

## Interpretation

The primary strict 95th-percentile burden metric is higher in ChILI/TOL than in control/AC in both the all-sample and coverage-QC analyses. The same direction is seen for the strict 90th percentile and strict top-1% burden mass, including leave-one-out stability for `strict_burden_q95` and `strict_top_1pct_burden_mass`.

This supports the hypothesis at the level of discovery-direction/stability for strict amyloid-labelled protein-coding transcripts. However, empirical p-values are not small and bootstrap CIs cross zero, so this should not be described as a validated biomarker.

Continuous all-annotated burden metrics do not show the same direction and should be treated as exploratory. The most defensible interpretation is therefore narrow: a weak but directionally stable strict-label upper-tail signal, not broad evidence that all amyloid-index-weighted annotated transcripts are shifted upward.

## Limitations

- Small sample size: all BC n = 11, QC n = 10.
- `SRR32060276` has low protein/amyloid annotation coverage and is excluded from the QC subset.
- Amyloid labels and `Amyloid_Index_max` are computational predictions.
- The burden score can be negative because the integrative amyloid index is z-score based.
- This is a sample-level tail analysis and does not adjust for batch, treatment, age, or sex.

## Next Best Step

Validate the tail endpoint in a larger independent baseline cohort, ideally with improved protein/amyloid annotation coverage and an orthogonal amyloid predictor panel for the top contributors.
"""
    (OUT / "AMYLOID_BURDEN_QUANTILE_SHIFT_REPORT.md").write_text(report, encoding="utf-8")

    print("DONE")
    print(f"Outputs: {OUT}")
    print(stats[stats["metric"].isin(["strict_burden_q90", "strict_burden_q95", "strict_top_1pct_burden_mass"])])


if __name__ == "__main__":
    main()
