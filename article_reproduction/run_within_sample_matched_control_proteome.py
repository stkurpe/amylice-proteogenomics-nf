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
BASE_OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili"
METHOD = "within_sample_matched_control_proteome"
OUT = BASE_OUT / METHOD
OUT.mkdir(parents=True, exist_ok=True)
GENCODE_CANDIDATES = [
    BASE_OUT / "gencode_v48_transcript_annotation.tsv",
    BASE_OUT / "rank_based_amyloid_burden_signature" / "gencode_v48_transcript_annotation.tsv",
    ROOT / "gencode_v48_transcript_annotation.tsv",
]

TIMEPOINT = "BC"
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
GROUP_LABELS = {"AC": "control", "TOL": "chili"}
GROUP_ORDER = ["AC", "TOL"]
COVERAGE_MIN = 0.5
PSEUDOCOUNT = 1.0
N_BOOT = 20000
RNG_SEED = 20260601


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
        return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - s.mean(skipna=True)) / sd


def zarr_group_to_df(g):
    return pd.DataFrame({col: g[col][:] for col in g.array_keys()})


def load_gencode_annotation():
    for path in GENCODE_CANDIDATES:
        if path.exists():
            gencode = pd.read_csv(path, sep="\t")
            required = {"Transcript_ID_clean", "gene_id", "gene_name", "biotype"}
            if required.issubset(gencode.columns):
                gencode = gencode[list(required)].drop_duplicates("Transcript_ID_clean")
                return gencode, path
    return pd.DataFrame(columns=["Transcript_ID_clean", "gene_id", "gene_name", "biotype"]), None


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


def bootstrap_ci_delta(ac, tol, rng, n_boot=N_BOOT):
    ac = np.asarray(ac, dtype=float)
    tol = np.asarray(tol, dtype=float)
    ac = ac[np.isfinite(ac)]
    tol = tol[np.isfinite(tol)]
    if len(ac) == 0 or len(tol) == 0:
        return np.nan, np.nan
    vals = np.empty(n_boot)
    for i in range(n_boot):
        vals[i] = rng.choice(tol, len(tol), replace=True).mean() - rng.choice(ac, len(ac), replace=True).mean()
    return np.quantile(vals, [0.025, 0.975])


def exact_label_permutation_p(values, labels, alternative="greater"):
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    keep = np.isfinite(values)
    values = values[keep]
    labels = labels[keep]
    n_tol = int(np.sum(labels == "TOL"))
    n = len(values)
    if n_tol == 0 or n_tol == n:
        return np.nan
    obs = values[labels == "TOL"].mean() - values[labels == "AC"].mean()
    diffs = []
    from itertools import combinations

    idx = np.arange(n)
    for tol_idx in combinations(idx, n_tol):
        tol_mask = np.zeros(n, dtype=bool)
        tol_mask[list(tol_idx)] = True
        diffs.append(values[tol_mask].mean() - values[~tol_mask].mean())
    diffs = np.asarray(diffs)
    if alternative == "greater":
        return float(np.mean(diffs >= obs - 1e-12))
    return float(np.mean(np.abs(diffs) >= abs(obs) - 1e-12))


def compare_groups(sample_df, metric, analysis_set, rng):
    d = sample_df.dropna(subset=[metric, "Group"]).copy()
    ac = d.loc[d["Group"].eq("AC"), metric].astype(float)
    tol = d.loc[d["Group"].eq("TOL"), metric].astype(float)
    if len(ac) < 2 or len(tol) < 2:
        return {}
    ci_low, ci_high = bootstrap_ci_delta(ac, tol, rng)
    values = d[metric].to_numpy(float)
    labels = d["Group"].astype(str).to_numpy()
    return {
        "analysis_set": analysis_set,
        "metric": metric,
        "n_AC": int(len(ac)),
        "n_TOL": int(len(tol)),
        "mean_AC": float(ac.mean()),
        "mean_TOL": float(tol.mean()),
        "median_AC": float(ac.median()),
        "median_TOL": float(tol.median()),
        "delta_mean_TOL_minus_AC": float(tol.mean() - ac.mean()),
        "delta_median_TOL_minus_AC": float(tol.median() - ac.median()),
        "bootstrap_delta_mean_95ci_low": float(ci_low),
        "bootstrap_delta_mean_95ci_high": float(ci_high),
        "cliffs_delta_TOL_vs_AC": float(cliffs_delta(tol, ac)),
        "mannwhitney_two_sided_p": float(mannwhitneyu(tol, ac, alternative="two-sided").pvalue),
        "welch_t_two_sided_p": float(ttest_ind(tol, ac, equal_var=False, nan_policy="omit").pvalue),
        "exact_permutation_one_sided_TOL_gt_AC_p": exact_label_permutation_p(values, labels, "greater"),
        "exact_permutation_two_sided_p": exact_label_permutation_p(values, labels, "two-sided"),
    }


def load_sample_metadata(root):
    annotation = pd.read_csv(ANNOTATION_PATH)
    annotation["raw_Group"] = annotation["Group"].astype(str)
    annotation["Group"] = annotation["raw_Group"].map(GROUP_MAP)
    annotation["Run"] = annotation["Run"].astype(str)
    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    meta = annotation[
        annotation["Run"].isin(zarr_samples)
        & annotation["time"].eq(TIMEPOINT)
        & annotation["Group"].isin(GROUP_ORDER)
    ].copy()
    meta = meta.rename(columns={"Run": "sample"})
    selected_samples = [s for s in zarr_samples if s in set(meta["sample"])]
    meta = meta.set_index("sample").loc[selected_samples].reset_index()
    meta["Group"] = pd.Categorical(meta["Group"], categories=GROUP_ORDER, ordered=True)
    return selected_samples, meta


def build_collapsed_annotation(root, selected_samples):
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
    combined["Amyloid_Index_max_component"] = (
        combined["beta_propensity_z"]
        + combined["KD_max_z"]
        + combined["aromaticity_z"]
        - combined["abs_charge_density_z"]
        - combined["pI_distance_from_7_z"]
    ) / 5

    expr_ann = zarr_group_to_df(root[f"samples/{selected_samples[0]}/expression"])[
        ["transcript_id", "gene_id", "gene_name"]
    ].copy()
    expr_ann["Transcript_ID_clean"] = expr_ann["transcript_id"].map(clean_enst)
    expr_ann = expr_ann[["Transcript_ID_clean", "gene_id", "gene_name"]].drop_duplicates("Transcript_ID_clean")
    gencode, gencode_path = load_gencode_annotation()
    expr_ann = expr_ann.merge(gencode, on="Transcript_ID_clean", how="left", suffixes=("_expr", "_gencode"))
    expr_ann["gene_id"] = expr_ann["gene_id_gencode"].fillna(expr_ann["gene_id_expr"])
    expr_ann["gene_name"] = expr_ann["gene_name_gencode"].fillna(expr_ann["gene_name_expr"])
    expr_ann = expr_ann[["Transcript_ID_clean", "gene_id", "gene_name", "biotype"]]
    expr_ann["subcellular_location"] = np.nan
    expr_ann["gc_content"] = np.nan

    amy_tx = (
        combined.dropna(subset=["sample", "Transcript_ID_clean"])
        .groupby(["sample", "Transcript_ID_clean"], as_index=False)
        .agg(
            Amyloid_Index_max=("Amyloid_Index_max_component", "max"),
            Amyloid_Index_mean=("Amyloid_Index_max_component", "mean"),
            AmyloGram_Prob_max=("AmyloGram_Prob", "max"),
            AMYPred_Prob_max=("AMYPred_Prob", "max"),
            protein_length_mean=("protein_length", "mean"),
            protein_length_median=("protein_length", "median"),
            n_protein_variants=("protein_id", "nunique"),
            Consensus_collapsed=("Consensus", consensus_priority),
        )
    )
    amy_tx = amy_tx.merge(expr_ann, on="Transcript_ID_clean", how="left")
    amy_tx["gencode_biotype"] = amy_tx["biotype"]
    amy_tx["gencode_is_protein_coding"] = amy_tx["biotype"].eq("protein_coding")
    amy_tx = amy_tx[amy_tx["gencode_is_protein_coding"]].copy()
    amy_tx["is_amyloid"] = amy_tx["Consensus_collapsed"].eq("Amyloid")
    amy_tx["is_non_amyloid"] = amy_tx["Consensus_collapsed"].eq("Non-Amyloid")
    return amy_tx, gencode_path


def add_expression(root, selected_samples, amy_tx):
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
    expr_amy = expr_long.merge(amy_tx, on=["sample", "Transcript_ID_clean"], how="inner")
    expr_amy["TPM"] = expr_amy["TPM"].astype(float).fillna(0.0)
    expr_amy["log2_TPM_plus_1"] = np.log2(expr_amy["TPM"] + PSEUDOCOUNT)
    expr_amy["log_protein_length"] = np.log(expr_amy["protein_length_mean"].clip(lower=1))
    return expr_amy, tpm


def choose_matches_for_sample(sample_df):
    amy = sample_df[sample_df["is_amyloid"]].copy()
    controls = sample_df[sample_df["is_non_amyloid"]].copy()
    amy = amy.dropna(subset=["Transcript_ID_clean", "log2_TPM_plus_1", "log_protein_length"])
    controls = controls.dropna(subset=["Transcript_ID_clean", "log2_TPM_plus_1", "log_protein_length"])
    if amy.empty or controls.empty:
        return pd.DataFrame()

    match_space = pd.concat(
        [
            amy[["log2_TPM_plus_1", "log_protein_length"]],
            controls[["log2_TPM_plus_1", "log_protein_length"]],
        ],
        ignore_index=True,
    )
    mu = match_space.mean()
    sd = match_space.std(ddof=0).replace(0, 1).fillna(1)
    for frame in (amy, controls):
        frame["z_log2_tpm"] = (frame["log2_TPM_plus_1"] - mu["log2_TPM_plus_1"]) / sd["log2_TPM_plus_1"]
        frame["z_log_length"] = (frame["log_protein_length"] - mu["log_protein_length"]) / sd["log_protein_length"]

    controls = controls.reset_index(drop=True)
    available = np.ones(len(controls), dtype=bool)
    rows = []
    amy = amy.sort_values(["TPM", "protein_length_mean"], ascending=[False, False]).reset_index(drop=True)
    for _, a in amy.iterrows():
        pool = controls[available]
        if pool.empty:
            break
        distance = np.sqrt(
            (pool["z_log2_tpm"] - a["z_log2_tpm"]) ** 2
            + (pool["z_log_length"] - a["z_log_length"]) ** 2
        )
        if a.get("biotype") == a.get("biotype") and pool["biotype"].notna().any():
            distance = distance + np.where(pool["biotype"].eq(a["biotype"]), 0, 2)
        if a.get("subcellular_location") == a.get("subcellular_location") and pool["subcellular_location"].notna().any():
            distance = distance + np.where(pool["subcellular_location"].eq(a["subcellular_location"]), 0, 2)
        c_idx = int(distance.idxmin())
        available[c_idx] = False
        c = controls.loc[c_idx]
        rows.append(
            {
                "sample": a["sample"],
                "amyloid_transcript": a["Transcript_ID_clean"],
                "control_transcript": c["Transcript_ID_clean"],
                "amyloid_gene_name": a.get("gene_name"),
                "control_gene_name": c.get("gene_name"),
                "amyloid_TPM": float(a["TPM"]),
                "control_TPM": float(c["TPM"]),
                "amyloid_log2_TPM_plus_1": float(a["log2_TPM_plus_1"]),
                "control_log2_TPM_plus_1": float(c["log2_TPM_plus_1"]),
                "paired_log2_expression_shift": float(a["log2_TPM_plus_1"] - c["log2_TPM_plus_1"]),
                "paired_TPM_shift": float(a["TPM"] - c["TPM"]),
                "amyloid_protein_length": float(a["protein_length_mean"]),
                "control_protein_length": float(c["protein_length_mean"]),
                "match_distance": float(distance.loc[c_idx]),
                "amyloid_index": float(a["Amyloid_Index_max"]),
                "amyloid_n_protein_variants": int(a["n_protein_variants"]),
                "control_n_protein_variants": int(c["n_protein_variants"]),
                "biotype_matched": bool(pd.isna(a.get("biotype")) or pd.isna(c.get("biotype")) or a.get("biotype") == c.get("biotype")),
                "subcellular_location_matched": bool(
                    pd.isna(a.get("subcellular_location"))
                    or pd.isna(c.get("subcellular_location"))
                    or a.get("subcellular_location") == c.get("subcellular_location")
                ),
            }
        )
    return pd.DataFrame(rows)


def summarize_samples(pairs, expr_amy, meta, total_tpm_by_sample):
    total_tpm = total_tpm_by_sample.rename("total_TPM")
    annotated_tpm = expr_amy.groupby("sample")["TPM"].sum().rename("protein_coding_annotated_TPM")
    strict = expr_amy[expr_amy["is_amyloid"]].groupby("sample")["TPM"].sum().rename("strict_amyloid_TPM")
    paired_summary = pairs.groupby("sample").agg(
        n_matched_pairs=("paired_log2_expression_shift", "size"),
        mean_paired_log2_expression_shift=("paired_log2_expression_shift", "mean"),
        median_paired_log2_expression_shift=("paired_log2_expression_shift", "median"),
        mean_paired_TPM_shift=("paired_TPM_shift", "mean"),
        median_paired_TPM_shift=("paired_TPM_shift", "median"),
        fraction_pairs_amyloid_gt_control=("paired_log2_expression_shift", lambda x: float(np.mean(np.asarray(x) > 0))),
        mean_match_distance=("match_distance", "mean"),
        median_match_distance=("match_distance", "median"),
        mean_amyloid_TPM=("amyloid_TPM", "mean"),
        mean_control_TPM=("control_TPM", "mean"),
        mean_amyloid_protein_length=("amyloid_protein_length", "mean"),
        mean_control_protein_length=("control_protein_length", "mean"),
    )
    out = meta[["sample"]].set_index("sample").join(paired_summary, how="left")
    out = out.join(total_tpm).join(annotated_tpm).join(strict).reset_index()
    out["n_matched_pairs"] = out["n_matched_pairs"].fillna(0).astype(int)
    out["protein_coding_annotated_TPM"] = out["protein_coding_annotated_TPM"].fillna(0.0)
    out["strict_amyloid_TPM"] = out["strict_amyloid_TPM"].fillna(0.0)
    out["protein_coding_annotation_fraction"] = out["protein_coding_annotated_TPM"] / out["total_TPM"]
    out["strict_amyloid_TPM_fraction"] = out["strict_amyloid_TPM"] / out["protein_coding_annotated_TPM"]
    out = out.merge(meta, on="sample", how="left")
    out["keep_coverage_qc"] = out["protein_coding_annotation_fraction"] >= COVERAGE_MIN
    return out


def leave_one_out(sample_summary, metric, analysis_set):
    rows = []
    d0 = sample_summary.dropna(subset=[metric, "Group"]).copy()
    if analysis_set == "coverage_qc":
        d0 = d0[d0["keep_coverage_qc"]].copy()
    for sample in d0["sample"]:
        d = d0[d0["sample"] != sample]
        ac = d.loc[d["Group"].eq("AC"), metric].astype(float)
        tol = d.loc[d["Group"].eq("TOL"), metric].astype(float)
        if len(ac) == 0 or len(tol) == 0:
            continue
        rows.append(
            {
                "analysis_set": analysis_set,
                "metric": metric,
                "left_out_sample": sample,
                "left_out_group": d0.loc[d0["sample"].eq(sample), "Group"].astype(str).iloc[0],
                "n_AC": len(ac),
                "n_TOL": len(tol),
                "delta_mean_TOL_minus_AC": tol.mean() - ac.mean(),
                "direction_TOL_gt_AC": bool(tol.mean() > ac.mean()),
            }
        )
    return rows


def write_report(sample_summary, stats, loo, meta, gencode_path):
    primary = stats[
        (stats["analysis_set"].eq("coverage_qc"))
        & (stats["metric"].eq("mean_paired_log2_expression_shift"))
    ].iloc[0]
    all_primary = stats[
        (stats["analysis_set"].eq("all_BC_samples"))
        & (stats["metric"].eq("mean_paired_log2_expression_shift"))
    ].iloc[0]
    excluded = sample_summary.loc[
        ~sample_summary["keep_coverage_qc"],
        ["sample", "Group", "raw_Group", "protein_coding_annotation_fraction"],
    ]
    no_pair = sample_summary.loc[
        sample_summary["n_matched_pairs"].eq(0),
        ["sample", "Group", "raw_Group", "protein_coding_annotation_fraction", "strict_amyloid_TPM_fraction"],
    ]
    def markdown_table(df):
        df = df.copy().fillna("")
        cols = list(df.columns)
        rows = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
        for _, row in df.iterrows():
            rows.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
        return "\n".join(rows)

    sample_lines = markdown_table(meta[["sample", "subject", "time", "raw_Group", "Group", "treatment", "batch", "AGE", "sex"]])
    excluded_lines = markdown_table(excluded) if not excluded.empty else "None"
    no_pair_lines = markdown_table(no_pair) if not no_pair.empty else "None"
    direction = "supports" if primary["delta_mean_TOL_minus_AC"] > 0 else "does not support"
    report = f"""# Within-Sample Matched Control Proteome Report

## What was done

Recomputed the baseline (`time == "BC"`) control-vs-ChILI analysis directly from `project.zarr` using a within-sample matched-control proteome design.
For each strict amyloidogenic transcript (`Consensus_collapsed == "Amyloid"`), protein variants were first collapsed to `sample × transcript`, then one non-amyloid transcript from the same sample was selected without replacement by nearest-neighbor matching on observed `log2(TPM + 1)` and log protein length. This avoids H1/H2 double-counting.

GENCODE transcript annotation was used to restrict the matched proteome to `biotype == "protein_coding"`, so the amyloid and matched non-amyloid pools are exactly matched on gene biotype at the analysis-design level. GC/content and subcellular localization were not available in the local files used here, so the executable distance covariates were expression and protein length.

GENCODE source: `{gencode_path}`.

## Samples used

{sample_lines}

## Coverage QC

QC threshold: `protein_coding_annotation_fraction >= {COVERAGE_MIN}`.

Excluded from the QC subset:

{excluded_lines}

Samples with zero strict amyloid matched pairs for the primary paired-shift endpoint:

{no_pair_lines}

## Main result

Primary matched endpoint: sample-level mean paired log2 expression shift, defined as:

`mean(log2(TPM_amyloid + 1) - log2(TPM_matched_non_amyloid + 1))`

All BC samples:

- n AC/control = {int(all_primary['n_AC'])}; n TOL/chili = {int(all_primary['n_TOL'])}
- mean AC = {all_primary['mean_AC']:.6g}; mean TOL = {all_primary['mean_TOL']:.6g}
- delta mean TOL - AC = {all_primary['delta_mean_TOL_minus_AC']:.6g}
- Cliff's delta TOL vs AC = {all_primary['cliffs_delta_TOL_vs_AC']:.6g}
- exact one-sided permutation p (TOL > AC) = {all_primary['exact_permutation_one_sided_TOL_gt_AC_p']:.6g}
- exact two-sided permutation p = {all_primary['exact_permutation_two_sided_p']:.6g}
- bootstrap 95% CI for delta mean = [{all_primary['bootstrap_delta_mean_95ci_low']:.6g}, {all_primary['bootstrap_delta_mean_95ci_high']:.6g}]

Coverage-QC subset:

- n AC/control = {int(primary['n_AC'])}; n TOL/chili = {int(primary['n_TOL'])}
- mean AC = {primary['mean_AC']:.6g}; mean TOL = {primary['mean_TOL']:.6g}
- delta mean TOL - AC = {primary['delta_mean_TOL_minus_AC']:.6g}
- Cliff's delta TOL vs AC = {primary['cliffs_delta_TOL_vs_AC']:.6g}
- exact one-sided permutation p (TOL > AC) = {primary['exact_permutation_one_sided_TOL_gt_AC_p']:.6g}
- exact two-sided permutation p = {primary['exact_permutation_two_sided_p']:.6g}
- bootstrap 95% CI for delta mean = [{primary['bootstrap_delta_mean_95ci_low']:.6g}, {primary['bootstrap_delta_mean_95ci_high']:.6g}]

## Interpretation

Under this within-sample matched-control design, the QC-subset direction {direction} the hypothesis that baseline ChILI/TOL samples have a higher amyloidogenic expression shift than control/AC samples. Because the sample size is small and the confidence interval is wide, this should be treated as discovery/exploratory evidence, not validation and not a proven biomarker.

## Limitations

- Only 11 BC samples from the metadata were present in `project.zarr`; the QC subset contains 10 samples after excluding `SRR32060276`.
- Matching used expression and protein length, with GENCODE `protein_coding` biotype enforced before matching. GC/content and subcellular localization were not locally available.
- Matching on observed expression can deliberately attenuate expression-burden differences; this method asks a conservative residual question rather than reproducing the unadjusted burden endpoint.
- Small n makes p-values unstable; effect size, direction, and leave-one-out stability are more informative here.

## Next best step

Validate the direction in an independent cohort or rerun this matched analysis after adding transcript GC and HPA/subcellular annotations, then use a counts-based model or sample-level mixed/permutation framework for inference.
"""
    (OUT / "WITHIN_SAMPLE_MATCHED_CONTROL_PROTEOME_REPORT.md").write_text(report, encoding="utf-8")


def main():
    rng = np.random.default_rng(RNG_SEED)
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    selected_samples, meta = load_sample_metadata(root)
    meta.to_csv(OUT / "sample_metadata.tsv", sep="\t", index=False)

    amy_tx, gencode_path = build_collapsed_annotation(root, selected_samples)
    amy_tx.to_csv(OUT / "collapsed_sample_transcript_amyloid_annotations.tsv", sep="\t", index=False)
    expr_amy, tpm = add_expression(root, selected_samples, amy_tx)
    expr_amy.to_csv(OUT / "annotated_expression_long.tsv", sep="\t", index=False)

    pairs = pd.concat(
        [choose_matches_for_sample(g) for _, g in expr_amy.groupby("sample", sort=False)],
        ignore_index=True,
    )
    pairs = pairs.merge(meta[["sample", "subject", "raw_Group", "Group", "time", "batch", "AGE", "sex"]], on="sample", how="left")
    pairs.to_csv(OUT / "matched_amyloid_nonamyloid_pairs.tsv", sep="\t", index=False)

    sample_summary = summarize_samples(pairs, expr_amy, meta, tpm.sum(axis=0))
    sample_summary.to_csv(OUT / "sample_level_matched_shift.tsv", sep="\t", index=False)

    metrics = [
        "mean_paired_log2_expression_shift",
        "median_paired_log2_expression_shift",
        "fraction_pairs_amyloid_gt_control",
        "mean_paired_TPM_shift",
        "strict_amyloid_TPM_fraction",
    ]
    rows = []
    for analysis_set, d in [
        ("all_BC_samples", sample_summary),
        ("coverage_qc", sample_summary[sample_summary["keep_coverage_qc"]].copy()),
    ]:
        for metric in metrics:
            rows.append(compare_groups(d, metric, analysis_set, rng))
    stats = pd.DataFrame([r for r in rows if r])
    stats.to_csv(OUT / "matched_shift_group_statistics.tsv", sep="\t", index=False)

    loo_rows = []
    for analysis_set in ["all_BC_samples", "coverage_qc"]:
        for metric in ["mean_paired_log2_expression_shift", "fraction_pairs_amyloid_gt_control"]:
            loo_rows.extend(leave_one_out(sample_summary, metric, analysis_set))
    loo = pd.DataFrame(loo_rows)
    loo.to_csv(OUT / "leave_one_out_stability.tsv", sep="\t", index=False)

    contributor_panel = (
        pairs.groupby(["amyloid_transcript", "amyloid_gene_name"], dropna=False)
        .agg(
            n_samples=("sample", "nunique"),
            mean_amyloid_TPM_AC=("amyloid_TPM", lambda x: x[pairs.loc[x.index, "Group"].eq("AC")].mean()),
            mean_amyloid_TPM_TOL=("amyloid_TPM", lambda x: x[pairs.loc[x.index, "Group"].eq("TOL")].mean()),
            mean_paired_shift_AC=("paired_log2_expression_shift", lambda x: x[pairs.loc[x.index, "Group"].eq("AC")].mean()),
            mean_paired_shift_TOL=("paired_log2_expression_shift", lambda x: x[pairs.loc[x.index, "Group"].eq("TOL")].mean()),
            mean_amyloid_index=("amyloid_index", "mean"),
        )
        .reset_index()
    )
    contributor_panel["delta_paired_shift_TOL_minus_AC"] = (
        contributor_panel["mean_paired_shift_TOL"] - contributor_panel["mean_paired_shift_AC"]
    )
    contributor_panel["delta_amyloid_TPM_TOL_minus_AC"] = (
        contributor_panel["mean_amyloid_TPM_TOL"] - contributor_panel["mean_amyloid_TPM_AC"]
    )
    contributor_panel.sort_values("delta_paired_shift_TOL_minus_AC", ascending=False).to_csv(
        OUT / "top_strict_amyloid_matched_contributor_panel.tsv", sep="\t", index=False
    )

    summary = {
        "method": METHOD,
        "n_selected_BC_samples_in_project_zarr": len(selected_samples),
        "group_counts": meta["Group"].astype(str).value_counts().to_dict(),
        "coverage_threshold": COVERAGE_MIN,
        "gencode_annotation_path": str(gencode_path) if gencode_path else None,
        "gencode_filter": "biotype == protein_coding",
        "excluded_coverage_qc_samples": sample_summary.loc[
            ~sample_summary["keep_coverage_qc"],
            ["sample", "Group", "raw_Group", "protein_coding_annotation_fraction"],
        ].to_dict("records"),
        "matching_covariates_used": ["GENCODE protein_coding biotype restriction", "log2_TPM_plus_1", "log_protein_length"],
        "matching_covariates_requested_but_unavailable": ["gc_content", "subcellular_location"],
        "n_matched_pairs_total": int(len(pairs)),
        "group_statistics": stats.to_dict("records"),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    write_report(sample_summary, stats, loo, meta, gencode_path)

    print("DONE")
    print(f"Outputs: {OUT}")
    print(stats[stats["metric"].eq("mean_paired_log2_expression_shift")].to_string(index=False))


if __name__ == "__main__":
    main()
