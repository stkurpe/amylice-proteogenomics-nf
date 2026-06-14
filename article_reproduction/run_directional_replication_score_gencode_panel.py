import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.stats import binomtest, mannwhitneyu, ttest_ind


ROOT = Path(__file__).resolve().parent
ZARR_PATH = ROOT / "project.zarr"
ANNOTATION_PATH = ROOT / "GSE287540_SraRunTable.csv"
GENCODE_CANDIDATES = [
    ROOT / "gencode_v48_transcript_annotation.tsv",
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "gencode_v48_transcript_annotation.tsv",
    ROOT / "analysis_outputs"
    / "amyloid_bc_control_vs_chili"
    / "rank_based_amyloid_burden_signature"
    / "gencode_v48_transcript_annotation.tsv",
    ROOT / "analysis_outputs"
    / "amyloid_bc_control_vs_chili"
    / "competitive_gene_set_test_gencode_v48_protein_coding"
    / "gencode_transcript_annotation_used.tsv",
]
METHOD_SAFE_NAME = "directional_replication_score_gencode_protein_coding_panel"
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / METHOD_SAFE_NAME
OUT.mkdir(parents=True, exist_ok=True)

TIMEPOINT = "BC"
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
COVERAGE_MIN = 0.5
PANEL_GENES = ["ZNF580", "STMP1", "FXYD5", "SEPTIN7", "AP2B1", "YLPM1", "SNX10", "TACC1", "HBP1", "PSMD8"]
PSEUDOCOUNT = 0.1


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
    if sd == 0 or pd.isna(sd):
        return pd.Series(np.nan, index=s.index)
    return (s - s.mean(skipna=True)) / sd


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


def load_gencode():
    needed = {"Transcript_ID_clean", "gene_id", "gene_name", "biotype"}
    for path in GENCODE_CANDIDATES:
        if path.exists():
            ann = pd.read_csv(path, sep="\t")
            if needed.issubset(ann.columns):
                ann = ann[list(needed)].drop_duplicates("Transcript_ID_clean")
                return ann, path
    raise FileNotFoundError("No local GENCODE transcript annotation with Transcript_ID_clean/gene_name/biotype was found.")


def build_selected_samples(root):
    annotation = pd.read_csv(ANNOTATION_PATH)
    annotation["raw_Group"] = annotation["Group"].astype(str)
    annotation["Group"] = annotation["raw_Group"].map(GROUP_MAP)
    annotation["Run"] = annotation["Run"].astype(str)
    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    sample_meta = annotation[
        annotation["Run"].isin(zarr_samples)
        & annotation["time"].eq(TIMEPOINT)
        & annotation["Group"].isin(["AC", "TOL"])
    ].copy()
    sample_meta = sample_meta.rename(columns={"Run": "sample"})
    selected = [s for s in zarr_samples if s in set(sample_meta["sample"])]
    sample_meta = sample_meta.set_index("sample").loc[selected].reset_index()
    sample_meta.to_csv(OUT / "sample_metadata.tsv", sep="\t", index=False)
    return sample_meta, selected, zarr_samples


def build_tpm(root, selected, zarr_samples):
    sample_idx = [zarr_samples.index(s) for s in selected]
    transcript_ids = pd.Series(root["layers/expression/transcript_ids"][:]).astype(str)
    transcript_clean = transcript_ids.map(clean_enst)
    tpm = pd.DataFrame(
        root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=transcript_clean,
        columns=selected,
    )
    tpm.index.name = "Transcript_ID_clean"
    tpm = tpm[~tpm.index.isna()].groupby(level=0).sum()
    return tpm


def build_collapsed_amyloid(root, selected, gencode_pc):
    amyloid_frames = []
    protein_frames = []
    for sample in selected:
        amyloid_frames.append(zarr_group_to_df(root[f"samples/{sample}/amyloid"]))
        protein_frames.append(zarr_group_to_df(root[f"samples/{sample}/protein_features"]))

    amyloid = pd.concat(amyloid_frames, ignore_index=True)
    protein = pd.concat(protein_frames, ignore_index=True)
    amyloid["ID"] = amyloid["sample"].astype(str) + "|" + amyloid["Sequence_ID"].astype(str)
    protein["ID"] = protein["sample"].astype(str) + "|" + protein["protein_id"].astype(str)
    combined = protein.merge(amyloid.drop(columns=["sample"], errors="ignore"), on="ID", how="inner")
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
    amy_tx["strict_amyloid"] = amy_tx["Consensus_collapsed"].eq("Amyloid")
    amy_tx.to_csv(OUT / "collapsed_sample_transcript_amyloid_annotations.tsv", sep="\t", index=False)
    return amy_tx


def build_coverage(tpm, amy_tx, sample_meta, protein_coding_transcripts):
    rows = []
    for sample in sample_meta["sample"]:
        protein_coding_tpm = float(tpm.loc[tpm.index.intersection(protein_coding_transcripts), sample].sum())
        covered = set(amy_tx.loc[amy_tx["sample"].eq(sample), "Transcript_ID_clean"])
        annotated_tpm = float(tpm.loc[tpm.index.intersection(covered), sample].sum())
        rows.append(
            {
                "sample": sample,
                "protein_coding_TPM_gencode": protein_coding_tpm,
                "protein_coding_annotated_TPM": annotated_tpm,
                "protein_coding_annotation_fraction": annotated_tpm / protein_coding_tpm if protein_coding_tpm else np.nan,
            }
        )
    coverage = pd.DataFrame(rows).merge(sample_meta, on="sample", how="left")
    coverage["keep_coverage_qc"] = coverage["protein_coding_annotation_fraction"] >= COVERAGE_MIN
    coverage.to_csv(OUT / "sample_coverage_qc.tsv", sep="\t", index=False)
    return coverage


def build_panel_long(tpm, amy_tx, gencode_pc, sample_meta):
    panel_tx = gencode_pc[gencode_pc["gene_name"].isin(PANEL_GENES)].copy()
    rows = []
    meta_cols = sample_meta.set_index("sample")[["Group", "raw_Group", "subject", "time", "batch", "treatment"]]
    for sample in sample_meta["sample"]:
        expr = tpm[[sample]].rename(columns={sample: "TPM"}).reset_index()
        expr = panel_tx.merge(expr, on="Transcript_ID_clean", how="left")
        ann = amy_tx[amy_tx["sample"].eq(sample)].copy()
        expr = expr.merge(
            ann[
                [
                    "Transcript_ID_clean",
                    "Consensus_collapsed",
                    "strict_amyloid",
                    "Amyloid_Index_max",
                    "Amyloid_Index_mean",
                    "n_protein_variants",
                ]
            ],
            on="Transcript_ID_clean",
            how="left",
        )
        expr["sample"] = sample
        expr["TPM"] = expr["TPM"].fillna(0).astype(float)
        expr["strict_amyloid"] = expr["strict_amyloid"].fillna(False).astype(bool)
        expr["strict_amyloid_TPM"] = np.where(expr["strict_amyloid"], expr["TPM"], 0.0)
        expr["strict_weighted_index"] = np.where(
            expr["strict_amyloid"],
            expr["TPM"] * expr["Amyloid_Index_max"].fillna(0),
            0.0,
        )
        expr["continuous_weighted_index"] = expr["TPM"] * expr["Amyloid_Index_max"].fillna(0)
        for col in meta_cols.columns:
            expr[col] = meta_cols.loc[sample, col]
        rows.append(expr)
    long_df = pd.concat(rows, ignore_index=True)
    long_df.to_csv(OUT / "panel_transcript_expression_long.tsv", sep="\t", index=False)
    return long_df


def summarize_sample_metrics(long_df, coverage):
    gene_sample = (
        long_df.groupby(["sample", "gene_name", "Group", "raw_Group", "subject", "batch", "treatment"], as_index=False)
        .agg(
            panel_gene_TPM=("TPM", "sum"),
            strict_amyloid_TPM=("strict_amyloid_TPM", "sum"),
            strict_weighted_index_sum=("strict_weighted_index", "sum"),
            continuous_weighted_index_sum=("continuous_weighted_index", "sum"),
            n_strict_amyloid_transcripts=("strict_amyloid", "sum"),
            n_strict_amyloid_expressed_TPM_ge_1=("strict_amyloid_TPM", lambda x: int((x >= 1).sum())),
        )
    )
    gene_sample["strict_amyloid_TPM_fraction"] = gene_sample["strict_amyloid_TPM"] / gene_sample["panel_gene_TPM"].replace(0, np.nan)
    gene_sample["continuous_expression_weighted_amyloid_index"] = (
        gene_sample["continuous_weighted_index_sum"] / gene_sample["panel_gene_TPM"].replace(0, np.nan)
    )
    gene_sample.to_csv(OUT / "panel_gene_sample_metrics.tsv", sep="\t", index=False)

    panel_sample = (
        long_df.groupby(["sample", "Group", "raw_Group", "subject", "batch", "treatment"], as_index=False)
        .agg(
            panel_total_TPM=("TPM", "sum"),
            strict_amyloid_TPM=("strict_amyloid_TPM", "sum"),
            strict_weighted_index_sum=("strict_weighted_index", "sum"),
            continuous_weighted_index_sum=("continuous_weighted_index", "sum"),
            n_strict_amyloid_expressed_TPM_ge_1=("strict_amyloid_TPM", lambda x: int((x >= 1).sum())),
        )
    )
    panel_sample["strict_amyloid_TPM_fraction"] = (
        panel_sample["strict_amyloid_TPM"] / panel_sample["panel_total_TPM"].replace(0, np.nan)
    )
    panel_sample["continuous_expression_weighted_amyloid_index"] = (
        panel_sample["continuous_weighted_index_sum"] / panel_sample["panel_total_TPM"].replace(0, np.nan)
    )
    panel_sample = panel_sample.merge(
        coverage[["sample", "protein_coding_annotation_fraction", "keep_coverage_qc"]],
        on="sample",
        how="left",
    )
    panel_sample.to_csv(OUT / "panel_sample_burden_metrics.tsv", sep="\t", index=False)
    return gene_sample, panel_sample


def compare_groups(values, metric):
    ac = values.loc[values["Group"].eq("AC"), metric].astype(float)
    tol = values.loc[values["Group"].eq("TOL"), metric].astype(float)
    row = {
        "metric": metric,
        "n_AC": int(ac.notna().sum()),
        "n_TOL": int(tol.notna().sum()),
        "mean_AC": float(ac.mean()),
        "mean_TOL": float(tol.mean()),
        "median_AC": float(ac.median()),
        "median_TOL": float(tol.median()),
        "delta_mean_TOL_minus_AC": float(tol.mean() - ac.mean()),
        "delta_median_TOL_minus_AC": float(tol.median() - ac.median()),
        "cliffs_delta_TOL_vs_AC": cliffs_delta(tol, ac),
    }
    if row["n_AC"] >= 2 and row["n_TOL"] >= 2:
        row["mannwhitney_greater_p"] = float(mannwhitneyu(tol, ac, alternative="greater").pvalue)
        row["welch_t_two_sided_p"] = float(ttest_ind(tol, ac, equal_var=False, nan_policy="omit").pvalue)
    else:
        row["mannwhitney_greater_p"] = np.nan
        row["welch_t_two_sided_p"] = np.nan
    return row


def markdown_table(df):
    if df.empty:
        return "None."
    show = df.copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "" if pd.isna(x) else f"{x:.6g}")
        else:
            show[col] = show[col].astype(str)
    header = "| " + " | ".join(show.columns) + " |"
    sep = "| " + " | ".join(["---"] * len(show.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in show.to_numpy(dtype=str)]
    return "\n".join([header, sep] + rows)


def directional_replication(gene_sample, panel_sample):
    all_rows = []
    summary_rows = []
    sample_stats_rows = []
    for analysis_set, samples in {
        "all_BC_samples": panel_sample["sample"].tolist(),
        "coverage_QC_protein_coding_annotation_fraction_ge_0_5": panel_sample.loc[
            panel_sample["keep_coverage_qc"], "sample"
        ].tolist(),
    }.items():
        gs = gene_sample[gene_sample["sample"].isin(samples)].copy()
        ps = panel_sample[panel_sample["sample"].isin(samples)].copy()
        for gene in PANEL_GENES:
            d = gs[gs["gene_name"].eq(gene)]
            ac = d.loc[d["Group"].eq("AC"), "strict_amyloid_TPM"].astype(float)
            tol = d.loc[d["Group"].eq("TOL"), "strict_amyloid_TPM"].astype(float)
            delta_mean = float(tol.mean() - ac.mean())
            delta_median = float(tol.median() - ac.median())
            row = {
                "analysis_set": analysis_set,
                "gene_name": gene,
                "n_AC": int(ac.notna().sum()),
                "n_TOL": int(tol.notna().sum()),
                "mean_strict_amyloid_TPM_AC": float(ac.mean()),
                "mean_strict_amyloid_TPM_TOL": float(tol.mean()),
                "median_strict_amyloid_TPM_AC": float(ac.median()),
                "median_strict_amyloid_TPM_TOL": float(tol.median()),
                "delta_mean_TOL_minus_AC": delta_mean,
                "delta_median_TOL_minus_AC": delta_median,
                "log2FC_mean_TOL_vs_AC": float(math.log2((tol.mean() + PSEUDOCOUNT) / (ac.mean() + PSEUDOCOUNT))),
                "cliffs_delta_TOL_vs_AC": cliffs_delta(tol, ac),
                "replicated_expected_direction": bool(delta_median > 0),
                "tie_no_direction": bool(delta_median == 0),
            }
            if len(ac) >= 2 and len(tol) >= 2:
                row["mannwhitney_greater_p"] = float(mannwhitneyu(tol, ac, alternative="greater").pvalue)
            else:
                row["mannwhitney_greater_p"] = np.nan
            all_rows.append(row)

        dr = pd.DataFrame([r for r in all_rows if r["analysis_set"] == analysis_set])
        informative = dr[~dr["tie_no_direction"]].copy()
        k = int(informative["replicated_expected_direction"].sum())
        n = int(len(informative))
        bt = binomtest(k, n, p=0.5, alternative="greater") if n else None
        ci = bt.proportion_ci(confidence_level=0.95, method="exact") if bt else None
        k_fixed = int(dr["replicated_expected_direction"].sum())
        n_fixed = int(len(PANEL_GENES))
        bt_fixed = binomtest(k_fixed, n_fixed, p=0.5, alternative="greater")
        ci_fixed = bt_fixed.proportion_ci(confidence_level=0.95, method="exact")
        summary_rows.append(
            {
                "analysis_set": analysis_set,
                "n_samples": int(ps["sample"].nunique()),
                "n_AC": int(ps.loc[ps["Group"].eq("AC"), "sample"].nunique()),
                "n_TOL": int(ps.loc[ps["Group"].eq("TOL"), "sample"].nunique()),
                "n_panel_genes": len(PANEL_GENES),
                "n_fixed_panel_replicated_expected_direction": k_fixed,
                "fixed_panel_directional_replication_score": k_fixed / n_fixed,
                "fixed_panel_binomial_greater_p_ties_as_not_replicated": float(bt_fixed.pvalue),
                "fixed_panel_score_ci95_low": float(ci_fixed.low),
                "fixed_panel_score_ci95_high": float(ci_fixed.high),
                "n_informative_genes": n,
                "n_replicated_expected_direction": k,
                "directional_replication_score": k / n if n else np.nan,
                "binomial_sign_test_greater_p": float(bt.pvalue) if bt else np.nan,
                "directional_replication_score_ci95_low": float(ci.low) if ci else np.nan,
                "directional_replication_score_ci95_high": float(ci.high) if ci else np.nan,
                "excluded_samples_qc": ",".join(ps.loc[~ps["keep_coverage_qc"], "sample"]) if "coverage_QC" not in analysis_set else "",
            }
        )
        for metric in [
            "strict_amyloid_TPM_fraction",
            "strict_amyloid_TPM",
            "n_strict_amyloid_expressed_TPM_ge_1",
            "strict_weighted_index_sum",
            "continuous_expression_weighted_amyloid_index",
        ]:
            row = compare_groups(ps, metric)
            row["analysis_set"] = analysis_set
            sample_stats_rows.append(row)

    direction = pd.DataFrame(all_rows)
    summary = pd.DataFrame(summary_rows)
    sample_stats = pd.DataFrame(sample_stats_rows)
    direction.to_csv(OUT / "panel_gene_directional_replication.tsv", sep="\t", index=False)
    summary.to_csv(OUT / "directional_replication_summary.tsv", sep="\t", index=False)
    sample_stats.to_csv(OUT / "panel_sample_burden_statistics.tsv", sep="\t", index=False)
    return direction, summary, sample_stats


def write_report(sample_meta, coverage, direction, summary, sample_stats, gencode_path):
    excluded = coverage.loc[~coverage["keep_coverage_qc"], ["sample", "Group", "raw_Group", "protein_coding_annotation_fraction"]]
    primary = sample_stats[sample_stats["metric"].eq("strict_amyloid_TPM_fraction")].copy()
    qc_summary = summary[summary["analysis_set"].str.startswith("coverage_QC")].iloc[0]
    all_summary = summary[summary["analysis_set"].eq("all_BC_samples")].iloc[0]
    qc_primary = primary[primary["analysis_set"].str.startswith("coverage_QC")].iloc[0]
    all_primary = primary[primary["analysis_set"].eq("all_BC_samples")].iloc[0]

    lines = [
        "# Directional Replication Score: GENCODE Protein-Coding Fixed Panel",
        "",
        "## What Was Done",
        "",
        f"Method: `{METHOD_SAFE_NAME}`.",
        "I re-read `project.zarr`, restricted the comparison to baseline `time == \"BC\"`, mapped `control -> AC` and `chili -> TOL`, and used only GENCODE `biotype == \"protein_coding\"` transcripts.",
        "Protein variants were collapsed to `sample x transcript` before joining expression, so H1/H2 protein variants cannot double-count transcript-level expression.",
        "Strict amyloid labels were `Consensus_collapsed == \"Amyloid\"`; continuous burden used `TPM x Amyloid_Index_max`.",
        "",
        "Fixed validation panel: " + ", ".join(PANEL_GENES) + ".",
        f"GENCODE source: `{gencode_path}`.",
        "",
        "## Samples",
        "",
        f"All BC samples: {sample_meta['sample'].nunique()} total; AC/control={int((sample_meta['Group'] == 'AC').sum())}, TOL/chili={int((sample_meta['Group'] == 'TOL').sum())}.",
        f"Coverage-QC subset: {int(coverage['keep_coverage_qc'].sum())} total; AC/control={int(((coverage['Group'] == 'AC') & coverage['keep_coverage_qc']).sum())}, TOL/chili={int(((coverage['Group'] == 'TOL') & coverage['keep_coverage_qc']).sum())}.",
        "",
        "QC-excluded samples:",
        markdown_table(excluded),
        "",
        "## Main Result",
        "",
        f"All BC samples: {int(all_summary['n_fixed_panel_replicated_expected_direction'])}/{int(all_summary['n_panel_genes'])} fixed-panel genes had higher median strict amyloid TPM in TOL than AC; fixed-panel directional replication score={all_summary['fixed_panel_directional_replication_score']:.3f}, exact binomial p with ties treated as not replicated={all_summary['fixed_panel_binomial_greater_p_ties_as_not_replicated']:.4g}. Among non-tie genes only: {int(all_summary['n_replicated_expected_direction'])}/{int(all_summary['n_informative_genes'])}, sign-test p={all_summary['binomial_sign_test_greater_p']:.4g}.",
        f"Coverage-QC subset: {int(qc_summary['n_fixed_panel_replicated_expected_direction'])}/{int(qc_summary['n_panel_genes'])} fixed-panel genes had higher median strict amyloid TPM in TOL than AC; fixed-panel directional replication score={qc_summary['fixed_panel_directional_replication_score']:.3f}, exact binomial p with ties treated as not replicated={qc_summary['fixed_panel_binomial_greater_p_ties_as_not_replicated']:.4g}. Among non-tie genes only: {int(qc_summary['n_replicated_expected_direction'])}/{int(qc_summary['n_informative_genes'])}, sign-test p={qc_summary['binomial_sign_test_greater_p']:.4g}.",
        "",
        "Primary panel burden endpoint (`strict_amyloid_TPM_fraction`):",
        f"All BC: mean TOL−AC={all_primary['delta_mean_TOL_minus_AC']:.6g}, median TOL−AC={all_primary['delta_median_TOL_minus_AC']:.6g}, Cliff's delta={all_primary['cliffs_delta_TOL_vs_AC']:.3f}, Mann-Whitney greater p={all_primary['mannwhitney_greater_p']:.4g}.",
        f"Coverage-QC: mean TOL−AC={qc_primary['delta_mean_TOL_minus_AC']:.6g}, median TOL−AC={qc_primary['delta_median_TOL_minus_AC']:.6g}, Cliff's delta={qc_primary['cliffs_delta_TOL_vs_AC']:.3f}, Mann-Whitney greater p={qc_primary['mannwhitney_greater_p']:.4g}.",
        "",
        "## Gene-Level Direction",
        "",
        markdown_table(
            direction[
                [
                    "analysis_set",
                    "gene_name",
                    "median_strict_amyloid_TPM_AC",
                    "median_strict_amyloid_TPM_TOL",
                    "delta_median_TOL_minus_AC",
                    "replicated_expected_direction",
                ]
            ]
        ),
        "",
        "## Interpretation",
        "",
        "This is a small fixed-panel validation-style analysis, not biomarker proof. The relevant signal is direction/stability across predefined contributors, with p-values treated as descriptive because the sample size is small.",
        "",
        "## Limitations",
        "",
        "- Panel was fixed from prior/discovery contributors, so this should be read as directional validation/exploration rather than independent biomarker discovery.",
        "- Amyloid calls and `Amyloid_Index_max` are computational predictions.",
        "- The cohort is small and includes treatment/batch heterogeneity.",
        "- One TOL sample, `SRR32060276`, has very low protein/amyloid annotation coverage and is excluded in the QC subset.",
        "",
        "## Next Best Step",
        "",
        "Run the same fixed panel on an independent or held-out baseline cohort with the panel frozen, then report only direction score and effect sizes as the primary validation readout.",
    ]
    (OUT / "DIRECTIONAL_REPLICATION_SCORE_REPORT.md").write_text("\n".join(lines) + "\n")


def main():
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    gencode, gencode_path = load_gencode()
    gencode_pc = gencode[gencode["biotype"].eq("protein_coding")].drop_duplicates("Transcript_ID_clean").copy()
    panel_missing = sorted(set(PANEL_GENES) - set(gencode_pc["gene_name"]))
    if panel_missing:
        raise ValueError(f"Panel genes absent from GENCODE protein-coding annotation: {panel_missing}")

    sample_meta, selected, zarr_samples = build_selected_samples(root)
    tpm = build_tpm(root, selected, zarr_samples)
    amy_tx = build_collapsed_amyloid(root, selected, gencode_pc)
    coverage = build_coverage(tpm, amy_tx, sample_meta, set(gencode_pc["Transcript_ID_clean"]))
    long_df = build_panel_long(tpm, amy_tx, gencode_pc, sample_meta)
    gene_sample, panel_sample = summarize_sample_metrics(long_df, coverage)
    direction, summary, sample_stats = directional_replication(gene_sample, panel_sample)

    metadata = {
        "method": METHOD_SAFE_NAME,
        "timepoint": TIMEPOINT,
        "group_map": GROUP_MAP,
        "panel_genes": PANEL_GENES,
        "strict_label": 'Consensus_collapsed == "Amyloid"',
        "continuous_score": "Amyloid_Index_max",
        "coverage_min": COVERAGE_MIN,
        "gencode_path": str(gencode_path),
        "outputs": sorted(p.name for p in OUT.iterdir() if p.is_file()),
        "summary": summary.to_dict(orient="records"),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(metadata, indent=2))
    write_report(sample_meta, coverage, direction, summary, sample_stats, gencode_path)
    print(summary.to_string(index=False))
    print(f"\nWrote: {OUT}")


if __name__ == "__main__":
    main()
