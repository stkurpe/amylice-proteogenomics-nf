import itertools
import json
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.stats import mannwhitneyu, rankdata, ttest_ind

warnings.filterwarnings("ignore")


ROOT = Path(__file__).resolve().parent
ZARR_PATH = ROOT / "project.zarr"
ANNOTATION_PATH = ROOT / "GSE287540_SraRunTable.csv"
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "rank_based_amyloid_burden_signature"
FIG = OUT / "figures"
GENCODE_CACHE = OUT / "gencode_v48_transcript_annotation.tsv"

TIMEPOINT = "BC"
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
GROUPS = ["AC", "TOL"]
COVERAGE_MIN = 0.5
N_BOOTSTRAP = 20000
BOOTSTRAP_SEED = 1701


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


def zscore(s):
    s = pd.Series(s, dtype=float)
    sd = s.std(skipna=True, ddof=0)
    if pd.isna(sd) or sd == 0:
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


def zarr_group_to_df(g):
    return pd.DataFrame({col: g[col][:] for col in g.array_keys()})


def load_gencode_transcript_annotation():
    if GENCODE_CACHE.exists():
        return pd.read_csv(GENCODE_CACHE, sep="\t")

    # Same logic as amyloid_expression_burden_analysis.ipynb, cell "Optional transcript annotation from GENCODE".
    gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz"
    gtf = pd.read_csv(
        gtf_url,
        sep="\t",
        comment="#",
        header=None,
        names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"],
        compression="gzip",
    )
    transcripts = gtf[gtf["feature"] == "transcript"].copy()

    def extract_attr(attr, key):
        m = re.search(fr'{key} "([^"]+)"', str(attr))
        return m.group(1) if m else np.nan

    transcripts["Transcript_ID"] = transcripts["attribute"].apply(lambda x: extract_attr(x, "transcript_id"))
    transcripts["Transcript_ID_clean"] = transcripts["Transcript_ID"].map(clean_enst)
    transcripts["gene_id"] = transcripts["attribute"].apply(lambda x: extract_attr(x, "gene_id"))
    transcripts["gene_name"] = transcripts["attribute"].apply(lambda x: extract_attr(x, "gene_name"))
    transcripts["biotype"] = transcripts["attribute"].apply(lambda x: extract_attr(x, "transcript_type"))
    out = transcripts[["Transcript_ID_clean", "gene_id", "gene_name", "biotype"]].drop_duplicates("Transcript_ID_clean")
    out.to_csv(GENCODE_CACHE, sep="\t", index=False)
    return out


def exact_permutation_p(ac, tol):
    ac = np.asarray(ac, dtype=float)
    tol = np.asarray(tol, dtype=float)
    observed = tol.mean() - ac.mean()
    pooled = np.concatenate([ac, tol])
    n_tol = len(tol)
    deltas = []
    idx = range(len(pooled))
    for tol_idx in itertools.combinations(idx, n_tol):
        tol_mask = np.zeros(len(pooled), dtype=bool)
        tol_mask[list(tol_idx)] = True
        deltas.append(pooled[tol_mask].mean() - pooled[~tol_mask].mean())
    deltas = np.asarray(deltas)
    one_sided = np.mean(deltas >= observed - 1e-15)
    two_sided = np.mean(np.abs(deltas) >= abs(observed) - 1e-15)
    return observed, one_sided, two_sided


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
        "n_AC": len(ac),
        "n_TOL": len(tol),
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


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

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
    gencode = load_gencode_transcript_annotation()
    gencode_pc = gencode[gencode["biotype"].eq("protein_coding")].drop_duplicates("Transcript_ID_clean").copy()
    protein_coding_transcripts = set(gencode_pc["Transcript_ID_clean"].dropna())

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
    tx_rows = []
    for sample in selected_samples:
        ann = amy_tx[amy_tx["sample"].eq(sample)].copy()
        expr = tpm[[sample]].rename(columns={sample: "TPM"}).reset_index()
        expr = expr[expr["Transcript_ID_clean"].isin(protein_coding_transcripts)].copy()
        expr = expr.merge(ann, on="Transcript_ID_clean", how="left")
        expr["strict_amyloid"] = expr["Consensus_collapsed"].eq("Amyloid")
        expr["TPM"] = expr["TPM"].fillna(0).astype(float)

        # Percentile is oriented so 1.0 means highest expression and 0.0 means lowest.
        n = len(expr)
        ascending_rank = rankdata(expr["TPM"].to_numpy(), method="average")
        expr["expression_rank_percentile"] = (ascending_rank - 1) / (n - 1) if n > 1 else np.nan
        expr["expression_rank_desc"] = rankdata(-expr["TPM"].to_numpy(), method="average")

        strict = expr[expr["strict_amyloid"]].copy()
        total_tpm = float(tpm[sample].sum())
        protein_coding_tpm = float(expr["TPM"].sum())
        annotated_tpm = float(expr.loc[expr["Amyloid_Index_max"].notna(), "TPM"].sum())
        row = {
            "sample": sample,
            "ranked_universe": "GENCODE_v48_protein_coding_expression_transcripts",
            "n_ranked_protein_coding_transcripts": int(n),
            "n_amyloid_annotated_protein_coding_transcripts": int(expr["Amyloid_Index_max"].notna().sum()),
            "n_strict_amyloid_transcripts": int(len(strict)),
            "strict_amyloid_transcripts_absent": bool(len(strict) == 0),
            "n_strict_amyloid_expressed_TPM_gt_0": int((strict["TPM"] > 0).sum()),
            "n_strict_amyloid_expressed_TPM_ge_1": int((strict["TPM"] >= 1).sum()),
            "total_TPM": total_tpm,
            "gencode_protein_coding_TPM": protein_coding_tpm,
            "protein_coding_annotated_TPM": annotated_tpm,
            "protein_coding_annotation_fraction": annotated_tpm / protein_coding_tpm if protein_coding_tpm else np.nan,
            "mean_rank_percentile_strict_amyloid": float(strict["expression_rank_percentile"].mean()),
            "rank_based_strict_amyloid_burden_score": float(strict["expression_rank_percentile"].mean()) if len(strict) else 0.0,
            "median_rank_percentile_strict_amyloid": float(strict["expression_rank_percentile"].median()),
            "strict_amyloid_top_10pct_fraction": float((strict["expression_rank_percentile"] >= 0.9).mean()),
            "strict_amyloid_top_25pct_fraction": float((strict["expression_rank_percentile"] >= 0.75).mean()),
            "mean_TPM_strict_amyloid": float(strict["TPM"].mean()),
            "sum_TPM_strict_amyloid": float(strict["TPM"].sum()),
        }
        rows.append(row)

        strict["sample"] = sample
        tx_rows.append(
            strict[
                [
                    "sample",
                    "Transcript_ID_clean",
                    "gene_id",
                    "gene_name",
                    "TPM",
                    "expression_rank_percentile",
                    "expression_rank_desc",
                    "Amyloid_Index_max",
                    "Consensus_collapsed",
                    "n_protein_variants",
                ]
            ]
        )

    sample_metrics = pd.DataFrame(rows).merge(sample_meta, on="sample", how="left")
    sample_metrics["keep_coverage_qc"] = sample_metrics["protein_coding_annotation_fraction"] >= COVERAGE_MIN
    sample_metrics.to_csv(OUT / "sample_rank_based_amyloid_burden.tsv", sep="\t", index=False)

    strict_long = pd.concat(tx_rows, ignore_index=True)
    strict_long = strict_long.merge(sample_meta[["sample", "Group", "raw_Group", "subject", "time", "batch"]], on="sample", how="left")
    strict_long.to_csv(OUT / "strict_amyloid_transcript_rank_percentiles.tsv", sep="\t", index=False)

    stats_rows = [
        compare_metric(sample_metrics, "rank_based_strict_amyloid_burden_score", "all_BC_samples"),
        compare_metric(
            sample_metrics[sample_metrics["keep_coverage_qc"]],
            "rank_based_strict_amyloid_burden_score",
            "coverage_qc_protein_annotation_fraction_ge_0.5",
        ),
        compare_metric(sample_metrics, "mean_rank_percentile_strict_amyloid", "sensitivity_strict_label_positive_samples"),
        compare_metric(
            sample_metrics[sample_metrics["keep_coverage_qc"]],
            "mean_rank_percentile_strict_amyloid",
            "sensitivity_strict_label_positive_coverage_qc",
        ),
    ]
    stats = pd.DataFrame(stats_rows)
    stats.to_csv(OUT / "rank_based_group_statistics.tsv", sep="\t", index=False)

    stability_rows = []
    for label, frame in [
        ("all_BC_samples", sample_metrics),
        ("coverage_qc_protein_annotation_fraction_ge_0.5", sample_metrics[sample_metrics["keep_coverage_qc"]]),
    ]:
        for sample in frame["sample"]:
            sub = frame[~frame["sample"].eq(sample)]
            if sub["Group"].nunique() == 2 and sub.groupby("Group").size().min() >= 2:
                comp = compare_metric(sub, "rank_based_strict_amyloid_burden_score", f"{label}_leave_one_out")
                comp["excluded_sample"] = sample
                stability_rows.append(comp)
    stability = pd.DataFrame(stability_rows)
    stability.to_csv(OUT / "leave_one_out_stability.tsv", sep="\t", index=False)

    top_panel = (
        strict_long.groupby("Transcript_ID_clean", as_index=False)
        .agg(
            gene_name=("gene_name", "first"),
            gene_id=("gene_id", "first"),
            mean_rank_percentile_AC=(
                "expression_rank_percentile",
                lambda x: x[strict_long.loc[x.index, "Group"].eq("AC")].mean(),
            ),
            mean_rank_percentile_TOL=(
                "expression_rank_percentile",
                lambda x: x[strict_long.loc[x.index, "Group"].eq("TOL")].mean(),
            ),
            mean_TPM_AC=("TPM", lambda x: x[strict_long.loc[x.index, "Group"].eq("AC")].mean()),
            mean_TPM_TOL=("TPM", lambda x: x[strict_long.loc[x.index, "Group"].eq("TOL")].mean()),
            n_samples_AC=("sample", lambda x: strict_long.loc[x.index, "sample"][strict_long.loc[x.index, "Group"].eq("AC")].nunique()),
            n_samples_TOL=("sample", lambda x: strict_long.loc[x.index, "sample"][strict_long.loc[x.index, "Group"].eq("TOL")].nunique()),
            Amyloid_Index_max=("Amyloid_Index_max", "max"),
        )
    )
    top_panel["delta_mean_rank_percentile_TOL_minus_AC"] = (
        top_panel["mean_rank_percentile_TOL"] - top_panel["mean_rank_percentile_AC"]
    )
    top_panel["delta_mean_TPM_TOL_minus_AC"] = top_panel["mean_TPM_TOL"] - top_panel["mean_TPM_AC"]
    top_panel = top_panel.sort_values("delta_mean_rank_percentile_TOL_minus_AC", ascending=False)
    top_panel.to_csv(OUT / "top_strict_amyloid_rank_contributor_panel.tsv", sep="\t", index=False)

    qc_excluded = sample_metrics.loc[
        ~sample_metrics["keep_coverage_qc"],
        ["sample", "Group", "raw_Group", "protein_coding_annotation_fraction"],
    ].to_dict("records")
    summary = {
        "method": "Rank-Based Amyloid Burden Signature",
        "rank_percentile_orientation": "1.0 = highest expressed among sample-level GENCODE v48 protein-coding transcripts; 0.0 = lowest",
        "ranked_universe": "GENCODE v48 transcript_type == protein_coding transcripts present in the expression layer",
        "gencode_cache": str(GENCODE_CACHE),
        "n_gencode_protein_coding_transcripts": int(gencode_pc.shape[0]),
        "n_gencode_protein_coding_transcripts_in_expression_layer": int(len(protein_coding_transcripts.intersection(set(tpm.index)))),
        "primary_metric": "rank_based_strict_amyloid_burden_score",
        "primary_metric_note": "Mean rank percentile of strict amyloid transcripts; set to 0 when a sample has no strict amyloid transcripts after sample x transcript collapse.",
        "timepoint": TIMEPOINT,
        "group_mapping": {"control": "AC", "chili": "TOL"},
        "n_samples": int(sample_metrics.shape[0]),
        "group_counts_all": sample_metrics["Group"].value_counts().to_dict(),
        "group_counts_qc": sample_metrics[sample_metrics["keep_coverage_qc"]]["Group"].value_counts().to_dict(),
        "coverage_qc_threshold": COVERAGE_MIN,
        "qc_excluded_samples": qc_excluded,
        "statistics": stats.to_dict("records"),
        "leave_one_out_direction_support_count": stability.groupby("analysis_set")["direction_supports_hypothesis"].sum().to_dict()
        if not stability.empty
        else {},
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    primary_all = stats.loc[stats["analysis_set"].eq("all_BC_samples")].iloc[0]
    primary_qc = stats.loc[stats["analysis_set"].eq("coverage_qc_protein_annotation_fraction_ge_0.5")].iloc[0]
    sensitivity_all = stats.loc[stats["analysis_set"].eq("sensitivity_strict_label_positive_samples")].iloc[0]
    sensitivity_qc = stats.loc[stats["analysis_set"].eq("sensitivity_strict_label_positive_coverage_qc")].iloc[0]
    report = f"""# Rank-Based Amyloid Burden Signature: BC control/AC vs chili/TOL

## What Was Done

This analysis tests whether strict consensus amyloid transcripts sit higher in the within-sample expression rank distribution at baseline (`time == "BC"`) in ChILI/TOL versus control/AC samples. It does not compare absolute TPM as the primary endpoint.

For each sample, GENCODE v48 `transcript_type == "protein_coding"` transcripts present in the expression layer were ranked by TPM. The primary metric is `rank_based_strict_amyloid_burden_score`: the mean expression rank percentile of transcripts with collapsed strict amyloid label `Consensus_collapsed == "Amyloid"`, with score set to 0 when a sample has no strict amyloid transcripts after collapse. Rank percentile is oriented so 1.0 means highest expression and 0.0 means lowest expression within that sample.

Protein variants were collapsed to sample x transcript before scoring, avoiding H1/H2 double-counting against transcript-level expression.

GENCODE source/cache: `{GENCODE_CACHE}`.

## Samples

- All BC samples: AC/control n = {int(primary_all['n_AC'])}, TOL/chili n = {int(primary_all['n_TOL'])}
- Coverage-QC subset (`protein_coding_annotation_fraction >= {COVERAGE_MIN}`): AC/control n = {int(primary_qc['n_AC'])}, TOL/chili n = {int(primary_qc['n_TOL'])}
- QC-excluded samples: {', '.join([x['sample'] for x in qc_excluded]) if qc_excluded else 'none'}

## Main Result

Primary metric: `rank_based_strict_amyloid_burden_score`.

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

## Interpretation

The all-sample rank-based burden score is higher in ChILI/TOL because two AC/control samples have no strict amyloid transcripts after sample x transcript collapse. The coverage-QC subset shows the same direction. This supports the hypothesis at the level of effect direction for this discovery endpoint, but the empirical p-values and bootstrap intervals do not justify claiming a validated biomarker.

Sensitivity analysis using only samples with at least one strict amyloid transcript gives the opposite direction:

- All strict-label-positive samples: delta mean TOL - AC = {sensitivity_all['delta_mean_TOL_minus_AC']:.4f}, empirical one-sided p (TOL > AC) = {sensitivity_all['empirical_permutation_p_one_sided_TOL_gt_AC']:.4f}
- Coverage-QC strict-label-positive samples: delta mean TOL - AC = {sensitivity_qc['delta_mean_TOL_minus_AC']:.4f}, empirical one-sided p (TOL > AC) = {sensitivity_qc['empirical_permutation_p_one_sided_TOL_gt_AC']:.4f}

Therefore, the main discovery signal depends on treating absence of strict amyloid transcripts as zero burden. Among samples where strict amyloid transcripts exist, their average expression rank is not higher in ChILI/TOL.

This should be treated as a discovery/pilot signal. Validation requires a larger independent cohort and, ideally, orthogonal amyloidogenicity predictors for top contributors.

## Limitations

- Small sample size limits statistical power.
- The ranked universe now uses independently downloaded GENCODE v48 protein-coding transcripts, but amyloid labels remain available only for transcripts with local protein/amyloid predictions.
- SRR32060276 has low protein/amyloid annotation coverage and is therefore excluded from the QC subset.
- Amyloid labels are computational predictions and need external validation.

## Next Best Step

Use this rank-based burden endpoint as the primary sample-level discovery metric in a larger baseline cohort, with SRR32060276-like samples excluded or repaired by improving annotation coverage. Follow up the top strict amyloid rank contributors with WALTZ/TANGO/AGGRESCAN or an equivalent orthogonal amyloid predictor panel.
"""
    (OUT / "RANK_BASED_AMYLOID_BURDEN_REPORT.md").write_text(report, encoding="utf-8")

    print("DONE")
    print(f"Outputs: {OUT}")
    print(stats)


if __name__ == "__main__":
    main()
