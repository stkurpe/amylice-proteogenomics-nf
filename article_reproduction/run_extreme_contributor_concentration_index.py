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
GENCODE_CANDIDATES = [
    ROOT / "gencode_v48_transcript_annotation.tsv",
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "gencode_v48_transcript_annotation.tsv",
    ROOT
    / "analysis_outputs"
    / "amyloid_bc_control_vs_chili"
    / "rank_based_amyloid_burden_signature"
    / "gencode_v48_transcript_annotation.tsv",
]
OUT = (
    ROOT
    / "analysis_outputs"
    / "amyloid_bc_control_vs_chili"
    / "extreme_contributor_concentration_index"
)
OUT.mkdir(parents=True, exist_ok=True)

TIMEPOINT = "BC"
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
COVERAGE_MIN = 0.5
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


def load_gencode_annotation():
    for path in GENCODE_CANDIDATES:
        if path.exists():
            ann = pd.read_csv(path, sep="\t")
            needed = {"Transcript_ID_clean", "gene_id", "gene_name", "biotype"}
            if needed.issubset(ann.columns):
                ann = ann[["Transcript_ID_clean", "gene_id", "gene_name", "biotype"]]
                ann = ann.drop_duplicates("Transcript_ID_clean")
                return ann, path
    raise FileNotFoundError("Could not find local GENCODE transcript annotation TSV.")


def gini(values):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    x = x[x > 0]
    if len(x) == 0:
        return np.nan
    if len(x) == 1:
        return 0.0
    x = np.sort(x)
    n = len(x)
    return float((2 * np.arange(1, n + 1) @ x) / (n * x.sum()) - (n + 1) / n)


def top_share(values, top_n=10):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    x = x[x > 0]
    total = x.sum()
    if total <= 0:
        return np.nan
    return float(np.sort(x)[::-1][:top_n].sum() / total)


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


def bootstrap_delta_ci(tol, ac, rng, n_boot=N_BOOT, alpha=0.05):
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


def permutation_p(values, groups, observed_delta):
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


def compare_metric(sample_metrics, metric, analysis_set, rng):
    d = sample_metrics.dropna(subset=[metric, "Group"]).copy()
    ac = d.loc[d["Group"].eq("AC"), metric].astype(float).to_numpy()
    tol = d.loc[d["Group"].eq("TOL"), metric].astype(float).to_numpy()
    if len(ac) < 2 or len(tol) < 2:
        return {}
    delta = float(np.mean(tol) - np.mean(ac))
    ci_low, ci_high = bootstrap_delta_ci(tol, ac, rng)
    perm_greater, perm_two, n_perm = permutation_p(d[metric], d["Group"], delta)
    return {
        "analysis_set": analysis_set,
        "metric": metric,
        "n_AC": int(len(ac)),
        "n_TOL": int(len(tol)),
        "mean_AC": float(np.mean(ac)),
        "mean_TOL": float(np.mean(tol)),
        "median_AC": float(np.median(ac)),
        "median_TOL": float(np.median(tol)),
        "delta_mean_TOL_minus_AC": delta,
        "bootstrap_delta_ci_low": float(ci_low),
        "bootstrap_delta_ci_high": float(ci_high),
        "cliffs_delta_TOL_vs_AC": cliffs_delta(tol, ac),
        "mannwhitney_two_sided_p": float(mannwhitneyu(tol, ac, alternative="two-sided").pvalue),
        "mannwhitney_one_sided_p_TOL_gt_AC": float(mannwhitneyu(tol, ac, alternative="greater").pvalue),
        "welch_t_two_sided_p": float(ttest_ind(tol, ac, equal_var=False, nan_policy="omit").pvalue),
        "exact_label_permutation_p_TOL_gt_AC": perm_greater,
        "exact_label_permutation_two_sided_p": perm_two,
        "n_exact_label_permutations": n_perm,
    }


def build_collapsed_amyloid(root, samples, gencode_ann):
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
        for col in [
            "beta_propensity",
            "KD_max",
            "aromaticity",
            "abs_charge_density",
            "pI_distance_from_7",
        ]:
            combined[f"{col}_z"] = zscore_series(combined[col])
        combined["Amyloid_Index"] = (
            combined["beta_propensity_z"]
            + combined["KD_max_z"]
            + combined["aromaticity_z"]
            - combined["abs_charge_density_z"]
            - combined["pI_distance_from_7_z"]
        ) / 5
        collapsed = (
            combined.dropna(subset=["sample", "Transcript_ID_clean"])
            .groupby(["sample", "Transcript_ID_clean"], as_index=False)
            .agg(
                Consensus_collapsed=("Consensus", consensus_priority),
                Amyloid_Index_max=("Amyloid_Index", "max"),
                AmyloGram_Prob_max=("AmyloGram_Prob", "max"),
                AMYPred_Prob_max=("AMYPred_Prob", "max"),
                n_protein_variants=("protein_id", "nunique"),
            )
        )
        rows.append(collapsed)
    collapsed = pd.concat(rows, ignore_index=True)
    collapsed = collapsed.merge(gencode_ann, on="Transcript_ID_clean", how="left")
    collapsed = collapsed[collapsed["biotype"].eq("protein_coding")].copy()
    collapsed["strict_consensus_amyloid"] = collapsed["Consensus_collapsed"].eq("Amyloid")
    return collapsed


def load_expression(root, samples):
    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    sample_idx = [zarr_samples.index(s) for s in samples]
    transcript_ids = pd.Series(root["layers/expression/transcript_ids"][:]).astype(str)
    transcript_clean = transcript_ids.map(clean_enst)
    tpm = pd.DataFrame(
        root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=transcript_clean,
        columns=samples,
    )
    tpm = tpm[~pd.isna(tpm.index)]
    tpm = tpm.groupby(level=0, sort=False).sum()
    return tpm


def compute_coverage(tpm, collapsed, sample_meta, protein_coding_transcripts):
    rows = []
    for sample in sample_meta["sample"]:
        vals = tpm[sample]
        pc_tpm = float(vals.loc[vals.index.intersection(protein_coding_transcripts)].sum())
        covered = set(collapsed.loc[collapsed["sample"].eq(sample), "Transcript_ID_clean"])
        ann_tpm = float(vals.loc[vals.index.intersection(covered)].sum())
        rows.append(
            {
                "sample": sample,
                "protein_coding_TPM_gencode": pc_tpm,
                "protein_coding_annotated_TPM": ann_tpm,
                "protein_coding_annotation_fraction": ann_tpm / pc_tpm if pc_tpm else np.nan,
            }
        )
    coverage = pd.DataFrame(rows).merge(sample_meta, on="sample", how="left")
    coverage["keep_coverage_qc"] = coverage["protein_coding_annotation_fraction"] >= COVERAGE_MIN
    return coverage


def build_contributions(tpm, collapsed):
    strict = collapsed[collapsed["strict_consensus_amyloid"]].copy()
    rows = []
    for sample in sorted(strict["sample"].unique()):
        sample_strict = strict[strict["sample"].eq(sample)].copy()
        sample_tpm = tpm[sample]
        sample_strict["TPM"] = sample_strict["Transcript_ID_clean"].map(sample_tpm).fillna(0.0)
        sample_strict["strict_TPM_contribution"] = sample_strict["TPM"].clip(lower=0)
        sample_strict["strict_weighted_positive_contribution"] = (
            sample_strict["TPM"].clip(lower=0) * sample_strict["Amyloid_Index_max"].clip(lower=0)
        )
        rows.append(sample_strict)
    return pd.concat(rows, ignore_index=True)


def summarize_sample_concentration(contrib, sample_meta, coverage):
    rows = []
    by_sample = {sample: g.copy() for sample, g in contrib.groupby("sample", sort=False)}
    for sample in sample_meta["sample"]:
        g = by_sample.get(sample, contrib.iloc[0:0].copy())
        tpm_contrib = g["strict_TPM_contribution"].to_numpy()
        weighted_contrib = g["strict_weighted_positive_contribution"].to_numpy()
        expressed = g[g["strict_TPM_contribution"] > 0].copy()
        rows.append(
            {
                "sample": sample,
                "n_strict_amyloid_transcripts_annotated": int(len(g)),
                "n_strict_amyloid_transcripts_expressed_TPM_gt_0": int((g["strict_TPM_contribution"] > 0).sum()),
                "n_strict_amyloid_transcripts_expressed_TPM_ge_1": int((g["strict_TPM_contribution"] >= 1).sum()),
                "strict_amyloid_TPM": float(np.sum(tpm_contrib)),
                "strict_amyloid_TPM_contribution_gini": gini(tpm_contrib),
                "strict_amyloid_TPM_top10_share": top_share(tpm_contrib, 10),
                "strict_amyloid_TPM_top5_share": top_share(tpm_contrib, 5),
                "strict_amyloid_TPM_top1_share": top_share(tpm_contrib, 1),
                "strict_weighted_positive_sum": float(np.sum(weighted_contrib)),
                "strict_weighted_positive_contribution_gini": gini(weighted_contrib),
                "strict_weighted_positive_top10_share": top_share(weighted_contrib, 10),
                "top_strict_TPM_transcript": expressed.sort_values("strict_TPM_contribution", ascending=False)
                .head(1)["Transcript_ID_clean"]
                .iloc[0]
                if len(expressed)
                else np.nan,
                "top_strict_TPM_gene": expressed.sort_values("strict_TPM_contribution", ascending=False)
                .head(1)["gene_name"]
                .iloc[0]
                if len(expressed)
                else np.nan,
                "top_strict_TPM": float(expressed["strict_TPM_contribution"].max()) if len(expressed) else np.nan,
            }
        )
    out = pd.DataFrame(rows).merge(sample_meta, on="sample", how="left")
    out = out.merge(
        coverage[["sample", "protein_coding_annotation_fraction", "keep_coverage_qc"]],
        on="sample",
        how="left",
    )
    return out


def summarize_top_panel(contrib, sample_meta):
    c = contrib.merge(sample_meta[["sample", "Group"]], on="sample", how="left")
    rows = []
    for group, g in c.groupby("Group"):
        summary = (
            g.groupby(["Transcript_ID_clean", "gene_name", "gene_id"], dropna=False, as_index=False)
            .agg(
                mean_strict_TPM_contribution=("strict_TPM_contribution", "mean"),
                median_strict_TPM_contribution=("strict_TPM_contribution", "median"),
                mean_weighted_positive_contribution=("strict_weighted_positive_contribution", "mean"),
                n_samples_expressed_TPM_ge_1=("strict_TPM_contribution", lambda x: int((x >= 1).sum())),
                Amyloid_Index_max=("Amyloid_Index_max", "max"),
                n_protein_variants_max=("n_protein_variants", "max"),
            )
            .sort_values("mean_strict_TPM_contribution", ascending=False)
            .head(25)
        )
        summary.insert(0, "Group", group)
        rows.append(summary)
    return pd.concat(rows, ignore_index=True)


def main():
    rng = np.random.default_rng(RANDOM_SEED)
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    gencode_ann, gencode_path = load_gencode_annotation()
    protein_coding_transcripts = set(
        gencode_ann.loc[gencode_ann["biotype"].eq("protein_coding"), "Transcript_ID_clean"].dropna().astype(str)
    )

    metadata = pd.read_csv(ANNOTATION_PATH)
    metadata["raw_Group"] = metadata["Group"].astype(str)
    metadata["Group"] = metadata["raw_Group"].map(GROUP_MAP)
    metadata = metadata.rename(columns={"Run": "sample"})
    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    selected = metadata[
        metadata["sample"].isin(zarr_samples)
        & metadata["time"].eq(TIMEPOINT)
        & metadata["Group"].isin(["AC", "TOL"])
    ].copy()
    selected = selected.set_index("sample").loc[[s for s in zarr_samples if s in set(selected["sample"])]].reset_index()

    tpm = load_expression(root, selected["sample"].tolist())
    collapsed = build_collapsed_amyloid(root, selected["sample"].tolist(), gencode_ann)
    coverage = compute_coverage(tpm, collapsed, selected, protein_coding_transcripts)
    contrib = build_contributions(tpm, collapsed)
    sample_metrics = summarize_sample_concentration(contrib, selected, coverage)
    top_panel = summarize_top_panel(contrib, selected)

    metrics = [
        "strict_amyloid_TPM_contribution_gini",
        "strict_amyloid_TPM_top10_share",
        "strict_amyloid_TPM_top5_share",
        "strict_amyloid_TPM_top1_share",
        "strict_weighted_positive_contribution_gini",
        "strict_weighted_positive_top10_share",
        "strict_amyloid_TPM",
        "n_strict_amyloid_transcripts_expressed_TPM_ge_1",
    ]
    stats = []
    sets = [
        ("all_BC_samples", sample_metrics),
        (
            "coverage_QC_protein_coding_annotation_fraction_ge_0_5",
            sample_metrics[sample_metrics["keep_coverage_qc"]].copy(),
        ),
    ]
    for name, df in sets:
        method_dir = OUT / name
        method_dir.mkdir(exist_ok=True)
        df.to_csv(method_dir / "sample_concentration_metrics.tsv", sep="\t", index=False)
        selected[selected["sample"].isin(df["sample"])].to_csv(method_dir / "sample_metadata.tsv", sep="\t", index=False)
        for metric in metrics:
            row = compare_metric(df, metric, name, rng)
            if row:
                stats.append(row)

    stats = pd.DataFrame(stats)
    coverage.to_csv(OUT / "sample_coverage_qc.tsv", sep="\t", index=False)
    collapsed.to_csv(OUT / "collapsed_sample_transcript_amyloid_annotations.tsv", sep="\t", index=False)
    contrib.to_csv(OUT / "strict_amyloid_contributions_long.tsv", sep="\t", index=False)
    sample_metrics.to_csv(OUT / "sample_concentration_metrics.tsv", sep="\t", index=False)
    stats.to_csv(OUT / "concentration_group_statistics.tsv", sep="\t", index=False)
    top_panel.to_csv(OUT / "top_strict_amyloid_contributor_panel.tsv", sep="\t", index=False)
    gencode_ann.to_csv(OUT / "gencode_transcript_annotation_used.tsv", sep="\t", index=False)

    primary_stats = stats[
        stats["metric"].isin(["strict_amyloid_TPM_contribution_gini", "strict_amyloid_TPM_top10_share"])
    ].copy()
    report = {
        "method": "Extreme Contributor Concentration Index: within-sample concentration of strict amyloid TPM contribution across GENCODE protein-coding transcripts. Primary metrics are Gini coefficient and top-10 contributor share. Protein variants were collapsed to sample x transcript before joining transcript-level expression.",
        "annotation_path": str(ANNOTATION_PATH),
        "gencode_annotation_path": str(gencode_path),
        "selected_BC_samples": selected[
            ["sample", "subject", "time", "raw_Group", "Group", "batch", "AGE", "sex", "treatment"]
        ].to_dict("records"),
        "excluded_by_coverage_qc": coverage.loc[
            ~coverage["keep_coverage_qc"], ["sample", "Group", "protein_coding_annotation_fraction"]
        ].to_dict("records"),
        "primary_results": primary_stats.to_dict("records"),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")

    lines = [
        "# Extreme Contributor Concentration Index",
        "",
        "## What was done",
        "Recomputed from `project.zarr` using only baseline (`time == \"BC\"`) samples. Metadata `Group=control` was mapped to AC and `Group=chili` to TOL. The transcript universe was restricted to local GENCODE `biotype == \"protein_coding\"` transcripts. Protein variants were collapsed to `sample x transcript` before joining transcript-level TPM, avoiding H1/H2 double counting.",
        "",
        "Strict amyloid labels were defined as `Consensus_collapsed == \"Amyloid\"`. The primary concentration metrics were the Gini coefficient of strict amyloid TPM contribution and the top-10 contributor share within each sample. I also computed secondary top-5/top-1 and positive weighted-contribution concentration metrics.",
        "",
        "## Samples and QC",
    ]
    for _, row in coverage.iterrows():
        flag = "kept in coverage-QC subset" if row["keep_coverage_qc"] else "excluded from coverage-QC subset"
        lines.append(
            f"- {row['sample']}: {row['Group']}/{row['raw_Group']}, subject {row['subject']}, "
            f"protein-coding annotation fraction={row['protein_coding_annotation_fraction']:.3f}, {flag}"
        )
    lines.extend(["", "## Main results"])
    for _, row in primary_stats.iterrows():
        lines.append(
            f"- {row['analysis_set']} / {row['metric']}: "
            f"mean AC={row['mean_AC']:.4f}, mean TOL={row['mean_TOL']:.4f}, "
            f"delta TOL-AC={row['delta_mean_TOL_minus_AC']:.4f} "
            f"(bootstrap 95% CI {row['bootstrap_delta_ci_low']:.4f} to {row['bootstrap_delta_ci_high']:.4f}); "
            f"Cliff's delta={row['cliffs_delta_TOL_vs_AC']:.3f}; "
            f"MW two-sided p={row['mannwhitney_two_sided_p']:.3g}; "
            f"exact permutation one-sided p(TOL>AC)={row['exact_label_permutation_p_TOL_gt_AC']:.3g}."
        )
    lines.extend(["", "## Interpretation"])
    qc_primary = primary_stats[
        primary_stats["analysis_set"].eq("coverage_QC_protein_coding_annotation_fraction_ge_0_5")
    ]
    support = (
        (qc_primary["delta_mean_TOL_minus_AC"] > 0).all()
        and (qc_primary["cliffs_delta_TOL_vs_AC"] > 0).all()
    )
    if support:
        lines.append(
            "In the coverage-QC subset, both primary concentration metrics are directionally higher in TOL/chili. "
            "Because the cohort is small, this is discovery-level evidence for a concentration/nucleation-style hypothesis, not a validated biomarker."
        )
    else:
        lines.append(
            "In the coverage-QC subset, the primary concentration metrics are not consistently higher in TOL/chili. "
            "This does not strongly support the concentration/nucleation-style hypothesis."
        )
    lines.extend(
        [
            "",
            "## Limitations",
            "Small n, treatment/batch imbalance, and amyloid annotation coverage limit inference. Concentration metrics can be sensitive to one highly expressed transcript, so I report both all-sample and coverage-QC subsets and provide top contributor panels for inspection.",
            "",
            "## Next best step",
            "Validate the concentration signal in an independent baseline cohort, or test whether the top strict amyloid contributors are reproducible under counts-based models with treatment/batch covariates.",
        ]
    )
    (OUT / "EXTREME_CONTRIBUTOR_CONCENTRATION_INDEX_REPORT.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )

    print(primary_stats.to_string(index=False))
    print(f"Outputs: {OUT}")


if __name__ == "__main__":
    main()
