import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parent
ZARR_PATH = ROOT / "validation_project.zarr"
ANNOTATION_PATH = ROOT / "GSE287540_SraRunTable.csv"
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "validation"
GENCODE_CANDIDATES = [
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "competitive_gene_set_test_gencode_v48_protein_coding" / "gencode_transcript_annotation_used.tsv",
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "rank_based_amyloid_burden_signature" / "gencode_v48_transcript_annotation.tsv",
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "gencode_v48_transcript_annotation.tsv",
    ROOT / "gencode_v48_transcript_annotation.tsv",
]

VALIDATION_SAMPLES = [
    "SRR32060218",
    "SRR32060222",
    "SRR32060224",
    "SRR32060226",
    "SRR32060252",
    "SRR32060255",
    "SRR32060299",
    "SRR32060304",
]
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
GROUPS = ["AC", "TOL"]
TIMEPOINT = "BC"
COVERAGE_MIN = 0.5
N_BOOTSTRAPS = 10_000
N_PERMUTATIONS = 10_000
RANDOM_SEED = 20260606
POSITIVE_CONSENSUS_LABEL = "Partial"
PRIMARY_TPM_COL = "partial_consensus_TPM"
PRIMARY_FRACTION_COL = "partial_consensus_TPM_fraction"
PRIMARY_N_COL = "n_partial_consensus_expressed_TPM_ge_1"
PRIMARY_ENDPOINT = PRIMARY_FRACTION_COL
PRIMARY_METRICS = [PRIMARY_FRACTION_COL, PRIMARY_TPM_COL, PRIMARY_N_COL]


def clean_enst(x):
    if pd.isna(x):
        return np.nan
    match = re.search(r"ENST\d+(?:\.\d+)?", str(x))
    return match.group(0).split(".")[0] if match else np.nan


def configure_positive_label(label):
    global POSITIVE_CONSENSUS_LABEL, PRIMARY_TPM_COL, PRIMARY_FRACTION_COL, PRIMARY_N_COL, PRIMARY_ENDPOINT, PRIMARY_METRICS
    POSITIVE_CONSENSUS_LABEL = label
    if label == "Amyloid":
        prefix = "strict_amyloid"
        count_prefix = "strict_amyloid"
    else:
        prefix = f"{label.lower()}_consensus"
        count_prefix = f"{label.lower()}_consensus"
    PRIMARY_TPM_COL = f"{prefix}_TPM"
    PRIMARY_FRACTION_COL = f"{prefix}_TPM_fraction"
    PRIMARY_N_COL = f"n_{count_prefix}_expressed_TPM_ge_1"
    PRIMARY_ENDPOINT = PRIMARY_FRACTION_COL
    PRIMARY_METRICS = [PRIMARY_FRACTION_COL, PRIMARY_TPM_COL, PRIMARY_N_COL]


def zarr_group_to_df(group):
    return pd.DataFrame({col: group[col][:] for col in group.array_keys()})


def load_gencode():
    needed = {"Transcript_ID_clean", "gene_id", "gene_name", "biotype"}
    for path in GENCODE_CANDIDATES:
        if path.exists():
            ann = pd.read_csv(path, sep="\t")
            if needed.issubset(ann.columns):
                ann = ann[list(needed)].drop_duplicates("Transcript_ID_clean")
                return ann, path
    raise FileNotFoundError("No local GENCODE transcript annotation with Transcript_ID_clean/gene_name/biotype.")


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


def zscore(values):
    s = pd.Series(values, dtype=float)
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


def bootstrap_ci(ac, tol, rng, n=N_BOOTSTRAPS):
    ac = np.asarray(ac, dtype=float)
    tol = np.asarray(tol, dtype=float)
    deltas = np.empty(n, dtype=float)
    for i in range(n):
        ac_b = rng.choice(ac, size=len(ac), replace=True)
        tol_b = rng.choice(tol, size=len(tol), replace=True)
        deltas[i] = np.mean(tol_b) - np.mean(ac_b)
    return np.nanpercentile(deltas, [2.5, 50, 97.5])


def empirical_group_permutation_p(ac, tol, rng, n=N_PERMUTATIONS):
    ac = np.asarray(ac, dtype=float)
    tol = np.asarray(tol, dtype=float)
    observed = float(np.mean(tol) - np.mean(ac))
    values = np.concatenate([ac, tol])
    n_ac = len(ac)
    null = np.empty(n, dtype=float)
    for i in range(n):
        perm = rng.permutation(values)
        null[i] = np.mean(perm[n_ac:]) - np.mean(perm[:n_ac])
    return observed, float((np.sum(null >= observed) + 1) / (n + 1))


def safe_qcut(values, q, prefix):
    s = pd.Series(values)
    ranks = s.rank(method="first")
    try:
        bins = pd.qcut(ranks, q=q, labels=False, duplicates="drop")
    except ValueError:
        bins = pd.Series(np.zeros(len(s), dtype=int), index=s.index)
    bins = bins.fillna(0).astype(int)
    return prefix + bins.astype(str)


def write_markdown_table(df, max_rows=20):
    if df.empty:
        return "_No rows._"
    display = df.head(max_rows).copy()
    for col in display.columns:
        if pd.api.types.is_float_dtype(display[col]):
            display[col] = display[col].map(lambda x: "" if pd.isna(x) else f"{x:.6g}")
        else:
            display[col] = display[col].map(lambda x: "" if pd.isna(x) else str(x))
    return display.to_markdown(index=False)


def write_blocked_report(sample_meta, manifest, reason):
    downloaded = manifest.groupby(["kind", "downloaded"]).size().reset_index(name="n") if not manifest.empty else pd.DataFrame()
    report = [
        "# Validation primary amyloidogenic burden",
        "",
        "## Status",
        "",
        f"Analysis status: blocked before statistical testing. {reason}",
        "",
        "## Validation cohort requested",
        "",
        write_markdown_table(sample_meta[["sample", "subject", "time", "raw_Group", "Group", "treatment", "batch"]]),
        "",
        "## Zarr manifest download status",
        "",
        write_markdown_table(downloaded),
        "",
        "## Pre-specified rule",
        "",
        f"Samples require protein_coding_annotation_fraction >= {COVERAGE_MIN} before statistical testing.",
        "No sample reached this step because expression and amyloid tables were unavailable in validation_project.zarr.",
    ]
    (OUT / "VALIDATION_PRIMARY_REPORT.md").write_text("\n".join(report) + "\n")


def build_sample_metadata():
    annotation = pd.read_csv(ANNOTATION_PATH)
    annotation["Run"] = annotation["Run"].astype(str)
    annotation["raw_Group"] = annotation["Group"].astype(str)
    annotation["Group"] = annotation["raw_Group"].map(GROUP_MAP)
    sample_meta = annotation[
        annotation["Run"].isin(VALIDATION_SAMPLES)
        & annotation["time"].eq(TIMEPOINT)
        & annotation["Group"].isin(GROUPS)
    ].copy()
    sample_meta = sample_meta.rename(columns={"Run": "sample"})
    sample_meta = sample_meta.set_index("sample").loc[VALIDATION_SAMPLES].reset_index()
    sample_meta.to_csv(OUT / "validation_sample_metadata.tsv", sep="\t", index=False)
    return sample_meta


def build_manifest(root):
    if "manifest" not in root:
        manifest = pd.DataFrame(columns=["sample", "kind", "s3", "downloaded"])
    else:
        manifest = zarr_group_to_df(root["manifest"])
    manifest.to_csv(OUT / "validation_project_manifest.tsv", sep="\t", index=False)
    return manifest


def collapse_amyloid_annotations(root, selected_samples, gencode_pc):
    amyloid_frames = []
    protein_frames = []
    for sample in selected_samples:
        sg = root[f"samples/{sample}"]
        if "amyloid" not in sg or "protein_features" not in sg:
            continue
        amyloid_frames.append(zarr_group_to_df(sg["amyloid"]))
        protein_frames.append(zarr_group_to_df(sg["protein_features"]))

    if not amyloid_frames or not protein_frames:
        return pd.DataFrame()

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
            protein_length_mean=("protein_length", "mean"),
            n_protein_variants=("protein_id", "nunique"),
            Consensus_collapsed=("Consensus", consensus_priority),
        )
    )
    amy_tx = amy_tx.merge(gencode_pc, on="Transcript_ID_clean", how="inner")
    amy_tx["is_primary_positive"] = amy_tx["Consensus_collapsed"].eq(POSITIVE_CONSENSUS_LABEL)
    amy_tx.to_csv(OUT / "collapsed_sample_transcript_amyloid_annotations_gencode_pc.tsv", sep="\t", index=False)
    return amy_tx


def build_tpm(root, selected_samples):
    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    sample_idx = [zarr_samples.index(s) for s in selected_samples]
    transcript_ids = pd.Series(root["layers/expression/transcript_ids"][:]).astype(str)
    transcript_clean = transcript_ids.map(clean_enst)
    tpm = pd.DataFrame(
        root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=transcript_clean,
        columns=selected_samples,
    )
    tpm.index.name = "Transcript_ID_clean"
    return tpm[~tpm.index.isna()].groupby(level=0).sum()


def compute_sample_burden(tpm, amy_tx, sample_meta, protein_coding_transcripts):
    rows = []
    amy_by_sample = {sample: amy_tx[amy_tx["sample"].eq(sample)] for sample in sample_meta["sample"]}
    for sample in sample_meta["sample"]:
        pc_tpm = tpm.loc[tpm.index.intersection(protein_coding_transcripts), sample]
        ann = amy_by_sample[sample].set_index("Transcript_ID_clean")
        annotated_ids = pc_tpm.index.intersection(ann.index)
        annotated_tpm = pc_tpm.loc[annotated_ids]
        ann = ann.loc[annotated_ids]
        positive_mask = ann["is_primary_positive"].astype(bool)
        positive_tpm = float(annotated_tpm.loc[positive_mask].sum())
        denominator = float(pc_tpm.sum())
        rows.append(
            {
                "sample": sample,
                "protein_coding_TPM_gencode": denominator,
                "protein_coding_annotated_TPM": float(annotated_tpm.sum()),
                "protein_coding_annotation_fraction": float(annotated_tpm.sum() / denominator) if denominator else np.nan,
                PRIMARY_TPM_COL: positive_tpm,
                PRIMARY_FRACTION_COL: positive_tpm / denominator if denominator else np.nan,
                PRIMARY_N_COL: int(((annotated_tpm >= 1.0) & positive_mask).sum()),
            }
        )
    burden = pd.DataFrame(rows).merge(sample_meta, on="sample", how="left")
    burden["keep_coverage_qc"] = burden["protein_coding_annotation_fraction"] >= COVERAGE_MIN
    burden.to_csv(OUT / "sample_level_primary_amyloid_burden.tsv", sep="\t", index=False)
    burden[["sample", "protein_coding_annotation_fraction", "keep_coverage_qc", "Group", "subject"]].to_csv(
        OUT / "sample_coverage_qc.tsv", sep="\t", index=False
    )
    return burden


def observed_statistics(burden, rng):
    rows = []
    kept = burden[burden["keep_coverage_qc"]].copy()
    for metric in PRIMARY_METRICS:
        ac = kept.loc[kept["Group"].eq("AC"), metric].astype(float).to_numpy()
        tol = kept.loc[kept["Group"].eq("TOL"), metric].astype(float).to_numpy()
        if len(ac) < 1 or len(tol) < 1:
            continue
        ci = bootstrap_ci(ac, tol, rng)
        observed, perm_p = empirical_group_permutation_p(ac, tol, rng)
        rows.append(
            {
                "metric": metric,
                "n_AC": len(ac),
                "n_TOL": len(tol),
                "mean_AC": float(np.mean(ac)),
                "mean_TOL": float(np.mean(tol)),
                "median_AC": float(np.median(ac)),
                "median_TOL": float(np.median(tol)),
                "observed_delta_TOL_minus_AC": observed,
                "mannwhitney_one_sided_TOL_gt_AC_p": float(mannwhitneyu(tol, ac, alternative="greater").pvalue),
                "cliffs_delta_TOL_vs_AC": float(cliffs_delta(tol, ac)),
                "rank_biserial_TOL_vs_AC": float(cliffs_delta(tol, ac)),
                "bootstrap_delta_ci2_5": float(ci[0]),
                "bootstrap_delta_median": float(ci[1]),
                "bootstrap_delta_ci97_5": float(ci[2]),
                "sample_label_permutation_one_sided_p": perm_p,
            }
        )
    stats = pd.DataFrame(rows)
    stats.to_csv(OUT / "primary_metric_group_statistics.tsv", sep="\t", index=False)
    return stats


def counterfactual_by_amyloid_labels(tpm, amy_tx, burden, rng):
    kept_samples = burden.loc[burden["keep_coverage_qc"], "sample"].tolist()
    if not kept_samples:
        return pd.DataFrame(), pd.DataFrame()

    mean_tpm = tpm[kept_samples].mean(axis=1)
    base = amy_tx[amy_tx["sample"].isin(kept_samples)].copy()
    base["mean_TPM_validation"] = base["Transcript_ID_clean"].map(mean_tpm)
    base = base.dropna(subset=["protein_length_mean", "mean_TPM_validation"])
    base["length_bin"] = safe_qcut(base["protein_length_mean"], 5, "L")
    base["expression_bin"] = safe_qcut(np.log1p(base["mean_TPM_validation"]), 5, "E")
    base["permutation_bin"] = base["length_bin"] + "_" + base["expression_bin"]

    observed_delta = float(
        burden.loc[burden["Group"].eq("TOL") & burden["keep_coverage_qc"], PRIMARY_ENDPOINT].mean()
        - burden.loc[burden["Group"].eq("AC") & burden["keep_coverage_qc"], PRIMARY_ENDPOINT].mean()
    )
    burden_idx = burden.set_index("sample")
    denominators = burden_idx.loc[kept_samples, "protein_coding_TPM_gencode"].astype(float).to_numpy()
    group_code = burden_idx.loc[kept_samples, "Group"].map({"AC": 0, "TOL": 1}).astype(int).to_numpy()
    sample_to_code = {sample: i for i, sample in enumerate(kept_samples)}
    sample_code = base["sample"].map(sample_to_code).astype(int).to_numpy()
    tx_codes = {tx: i for i, tx in enumerate(tpm.index)}
    transcript_code = base["Transcript_ID_clean"].map(tx_codes).to_numpy()
    valid = ~pd.isna(transcript_code)
    base = base.loc[valid].reset_index(drop=True)
    sample_code = sample_code[valid]
    transcript_code = transcript_code[valid].astype(int)
    tpm_matrix = tpm[kept_samples].to_numpy()
    row_tpm = tpm_matrix[transcript_code, sample_code]
    positive_labels = base["is_primary_positive"].astype(bool).to_numpy()
    bin_groups = [idx.to_numpy(dtype=int) for _, idx in base.groupby("permutation_bin").groups.items()]

    null = np.empty(N_PERMUTATIONS, dtype=float)
    for i in range(N_PERMUTATIONS):
        permuted = np.empty_like(positive_labels)
        for idx in bin_groups:
            permuted[idx] = rng.permutation(positive_labels[idx])
        positive_tpm = np.bincount(sample_code, weights=row_tpm * permuted, minlength=len(kept_samples))
        fractions = positive_tpm / denominators
        null[i] = np.nanmean(fractions[group_code == 1]) - np.nanmean(fractions[group_code == 0])

    empirical_p = float((np.sum(null >= observed_delta) + 1) / (N_PERMUTATIONS + 1))
    summary = pd.DataFrame(
        [
            {
                "metric": PRIMARY_ENDPOINT,
                "n_permutations": N_PERMUTATIONS,
                "observed_delta_TOL_minus_AC": observed_delta,
                "empirical_p_TOL_minus_AC_ge_observed": empirical_p,
                "observed_percentile_in_null": float(np.mean(null <= observed_delta) * 100),
                "null_mean": float(np.mean(null)),
                "null_ci2_5": float(np.percentile(null, 2.5)),
                "null_ci97_5": float(np.percentile(null, 97.5)),
            }
        ]
    )
    null_df = pd.DataFrame({"permutation": np.arange(1, N_PERMUTATIONS + 1), "null_delta_TOL_minus_AC": null})
    summary.to_csv(OUT / "counterfactual_amyloid_label_permutation_summary.tsv", sep="\t", index=False)
    null_df.to_csv(OUT / "counterfactual_amyloid_label_permutation_null.tsv", sep="\t", index=False)
    return summary, null_df


def leave_one_out(burden):
    kept = burden[burden["keep_coverage_qc"]].copy()
    rows = []
    for sample in kept["sample"]:
        d = kept[~kept["sample"].eq(sample)]
        for metric in PRIMARY_METRICS:
            ac = d.loc[d["Group"].eq("AC"), metric].astype(float)
            tol = d.loc[d["Group"].eq("TOL"), metric].astype(float)
            rows.append(
                {
                    "left_out_sample": sample,
                    "left_out_subject": kept.set_index("sample").loc[sample, "subject"],
                    "metric": metric,
                    "n_AC": int(len(ac)),
                    "n_TOL": int(len(tol)),
                    "delta_TOL_minus_AC": float(tol.mean() - ac.mean()) if len(ac) and len(tol) else np.nan,
                    "direction_TOL_gt_AC": bool((tol.mean() - ac.mean()) > 0) if len(ac) and len(tol) else False,
                }
            )
    loo = pd.DataFrame(rows)
    loo.to_csv(OUT / "leave_one_sample_out_stability.tsv", sep="\t", index=False)
    return loo


def write_success_report(sample_meta, coverage, stats, label_perm, loo, label_qc, gencode_path):
    primary = stats[stats["metric"].eq(PRIMARY_ENDPOINT)]
    report = [
        "# Validation primary amyloidogenic burden",
        "",
        "## Cohort",
        "",
        write_markdown_table(sample_meta[["sample", "subject", "time", "raw_Group", "Group", "treatment", "batch"]]),
        "",
        "## GENCODE filter",
        "",
        f"GENCODE annotation used: `{gencode_path}`. Only `biotype == protein_coding` transcripts were retained.",
        f"Primary-positive amyloid consensus label: `{POSITIVE_CONSENSUS_LABEL}`.",
        "",
        "## Coverage QC",
        "",
        write_markdown_table(coverage[["sample", "Group", "subject", "protein_coding_annotation_fraction", "keep_coverage_qc"]]),
        "",
        f"## Primary {POSITIVE_CONSENSUS_LABEL}-consensus burden statistics",
        "",
        write_markdown_table(stats),
        "",
        "## Amyloid label QC",
        "",
        write_markdown_table(label_qc),
        "",
        "## Counterfactual amyloid-label permutation",
        "",
        write_markdown_table(label_perm),
        "",
        "## Leave-one-sample-out stability",
        "",
        write_markdown_table(loo),
    ]
    if not primary.empty:
        p = primary.iloc[0]
        report.extend(
            [
                "",
                "## Primary endpoint interpretation",
                "",
                (
                    f"For {PRIMARY_ENDPOINT}, observed delta TOL minus AC was "
                    f"{p['observed_delta_TOL_minus_AC']:.6g}; one-sided Mann-Whitney p="
                    f"{p['mannwhitney_one_sided_TOL_gt_AC_p']:.6g}; Cliff's delta="
                    f"{p['cliffs_delta_TOL_vs_AC']:.6g}."
                ),
            ]
        )
    (OUT / "VALIDATION_PRIMARY_REPORT.md").write_text("\n".join(report) + "\n")


def parse_args():
    parser = argparse.ArgumentParser(description="Run validation primary amyloidogenic burden analysis.")
    parser.add_argument("--zarr", default=str(ZARR_PATH), help="Validation Zarr path.")
    parser.add_argument("--out", default=str(OUT), help="Output directory.")
    parser.add_argument("--positive-consensus-label", default=POSITIVE_CONSENSUS_LABEL, help="Consensus label treated as positive.")
    return parser.parse_args()


def main():
    global ZARR_PATH, OUT
    args = parse_args()
    ZARR_PATH = Path(args.zarr)
    OUT = Path(args.out)
    configure_positive_label(args.positive_consensus_label)
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    sample_meta = build_sample_metadata()
    manifest = build_manifest(root)

    summary = {
        "zarr_path": str(ZARR_PATH),
        "out_dir": str(OUT),
        "requested_samples": VALIDATION_SAMPLES,
        "coverage_min": COVERAGE_MIN,
        "n_bootstraps": N_BOOTSTRAPS,
        "n_permutations": N_PERMUTATIONS,
        "random_seed": RANDOM_SEED,
    }

    has_expression = "layers" in root and "expression" in root["layers"]
    has_all_sample_tables = all(
        sample in root["samples"] and "amyloid" in root[f"samples/{sample}"] and "protein_features" in root[f"samples/{sample}"]
        for sample in sample_meta["sample"]
    )
    if not has_expression or not has_all_sample_tables:
        reason = "validation_project.zarr lacks expression and/or amyloid/protein feature tables for the requested samples."
        summary.update({"status": "blocked_missing_validation_inputs", "reason": reason})
        write_blocked_report(sample_meta, manifest, reason)
        (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        return

    gencode, gencode_path = load_gencode()
    gencode_pc = gencode[gencode["biotype"].eq("protein_coding")].copy()
    gencode_pc.to_csv(OUT / "gencode_protein_coding_transcripts_used.tsv", sep="\t", index=False)

    selected_samples = sample_meta["sample"].tolist()
    tpm = build_tpm(root, selected_samples)
    amy_tx = collapse_amyloid_annotations(root, selected_samples, gencode_pc)
    if amy_tx.empty:
        reason = "no amyloid annotations remained after GENCODE protein_coding filtering."
        summary.update({"status": "blocked_no_gencode_pc_amyloid_annotations", "reason": reason})
        write_blocked_report(sample_meta, manifest, reason)
        (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        return

    burden = compute_sample_burden(tpm, amy_tx, sample_meta, set(gencode_pc["Transcript_ID_clean"]))
    label_qc = (
        amy_tx.groupby(["sample", "Consensus_collapsed"], dropna=False)
        .agg(
            n_transcripts=("Transcript_ID_clean", "nunique"),
            n_primary_positive=("is_primary_positive", "sum"),
        )
        .reset_index()
        .merge(sample_meta[["sample", "Group", "subject"]], on="sample", how="left")
    )
    label_qc.to_csv(OUT / "amyloid_label_qc.tsv", sep="\t", index=False)
    stats = observed_statistics(burden, rng)
    label_perm, _ = counterfactual_by_amyloid_labels(tpm, amy_tx, burden, rng)
    loo = leave_one_out(burden)
    write_success_report(sample_meta, burden, stats, label_perm, loo, label_qc, gencode_path)
    summary.update(
        {
            "status": "completed",
            "gencode_annotation": str(gencode_path),
            "positive_consensus_label": POSITIVE_CONSENSUS_LABEL,
            "n_gencode_protein_coding_transcripts": int(len(gencode_pc)),
            "n_samples_requested": int(len(VALIDATION_SAMPLES)),
            "n_samples_metadata": int(len(sample_meta)),
            "n_samples_passing_coverage_qc": int(burden["keep_coverage_qc"].sum()),
        }
    )
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
