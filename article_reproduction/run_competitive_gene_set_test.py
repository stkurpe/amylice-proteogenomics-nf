import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parent
ZARR_PATH = ROOT / "project.zarr"
ANNOTATION_PATH = ROOT / "GSE287540_SraRunTable.csv"
GENCODE_CANDIDATES = [
    ROOT / "gencode_v48_transcript_annotation.tsv",
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "gencode_v48_transcript_annotation.tsv",
    ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "rank_based_amyloid_burden_signature" / "gencode_v48_transcript_annotation.tsv",
]
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "competitive_gene_set_test_gencode_v48_protein_coding"
OUT.mkdir(parents=True, exist_ok=True)

TIMEPOINT = "BC"
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
COVERAGE_MIN = 0.5
MIN_COUNT = 10
MIN_SAMPLES = 3
PSEUDOCOUNT = 1.0
N_GSEA_PERMUTATIONS = 10000
N_RANDOM_SET_PERMUTATIONS = 10000
RANDOM_SEED = 20260601


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


def zarr_group_to_df(g):
    return pd.DataFrame({col: g[col][:] for col in g.array_keys()})


def load_gencode_annotation():
    for path in GENCODE_CANDIDATES:
        if path.exists():
            ann = pd.read_csv(path, sep="\t")
            needed = {"Transcript_ID_clean", "gene_id", "gene_name", "biotype"}
            if needed.issubset(ann.columns):
                ann = ann[["Transcript_ID_clean", "gene_id", "gene_name", "biotype"]].drop_duplicates("Transcript_ID_clean")
                return ann, path
    raise FileNotFoundError("Could not find local GENCODE transcript annotation TSV.")


def signed_neglog10_p(p_values, effects):
    p = np.asarray(p_values, dtype=float)
    effects = np.asarray(effects, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    p = np.clip(p, np.nextafter(0, 1), 1.0)
    return np.sign(effects) * -np.log10(p)


def gsea_es(ranked_scores, is_set):
    scores = np.asarray(ranked_scores, dtype=float)
    hits = np.asarray(is_set, dtype=bool)
    n_hits = int(hits.sum())
    n_total = len(scores)
    n_miss = n_total - n_hits
    if n_hits == 0 or n_miss == 0:
        return np.nan
    weights = np.abs(scores)
    hit_weights = weights[hits]
    if hit_weights.sum() == 0:
        hit_step = hits.astype(float) / n_hits
    else:
        hit_step = np.where(hits, weights / hit_weights.sum(), 0.0)
    miss_step = np.where(~hits, 1.0 / n_miss, 0.0)
    running = np.cumsum(hit_step - miss_step)
    max_es = running.max()
    min_es = running.min()
    return max_es if abs(max_es) >= abs(min_es) else min_es


def normalized_es(observed_es, null_es):
    null_es = np.asarray(null_es, dtype=float)
    if observed_es >= 0:
        denom = np.mean(null_es[null_es >= 0]) if np.any(null_es >= 0) else np.nan
    else:
        denom = abs(np.mean(null_es[null_es < 0])) if np.any(null_es < 0) else np.nan
    return observed_es / denom if denom and np.isfinite(denom) else np.nan


def empirical_p(observed, null_values, alternative):
    null_values = np.asarray(null_values, dtype=float)
    null_values = null_values[np.isfinite(null_values)]
    if len(null_values) == 0 or not np.isfinite(observed):
        return np.nan
    if alternative == "greater":
        return (np.sum(null_values >= observed) + 1) / (len(null_values) + 1)
    if alternative == "less":
        return (np.sum(null_values <= observed) + 1) / (len(null_values) + 1)
    return (np.sum(np.abs(null_values) >= abs(observed)) + 1) / (len(null_values) + 1)


def bootstrap_ci(values, rng, n_boot=10000, alpha=0.05):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 2:
        return (np.nan, np.nan)
    idx = rng.integers(0, len(values), size=(n_boot, len(values)))
    means = values[idx].mean(axis=1)
    return tuple(np.quantile(means, [alpha / 2, 1 - alpha / 2]))


def build_collapsed_amyloid(root, samples):
    rows = []
    for sample in samples:
        amy = zarr_group_to_df(root[f"samples/{sample}/amyloid"])
        pf = zarr_group_to_df(root[f"samples/{sample}/protein_features"])
        if "protein_id" not in pf.columns:
            continue
        amy["ID"] = amy["sample"].astype(str) + "|" + amy["Sequence_ID"].astype(str)
        pf["ID"] = pf["sample"].astype(str) + "|" + pf["protein_id"].astype(str)
        amy = pf[["ID", "sample"]].merge(amy.drop(columns=["sample"], errors="ignore"), on="ID", how="inner")
        amy["Transcript_ID_clean"] = amy["Sequence_ID"].map(clean_enst)
        amy = amy.dropna(subset=["Transcript_ID_clean"])
        collapsed = (
            amy.groupby(["sample", "Transcript_ID_clean"], as_index=False)
            .agg(Consensus_collapsed=("Consensus", consensus_priority))
        )
        rows.append(collapsed)
    collapsed = pd.concat(rows, ignore_index=True)
    collapsed["strict_consensus_amyloid"] = collapsed["Consensus_collapsed"].eq("Amyloid")
    return collapsed


def compute_coverage(root, samples, sample_meta, transcript_clean, protein_coding_transcripts):
    tpm_mat = root["layers/expression/tpm"]
    sample_ids = [str(x) for x in root["layers/expression/sample_ids"][:]]
    sample_idx = [sample_ids.index(s) for s in samples]
    tpm = pd.DataFrame(
        tpm_mat.get_orthogonal_selection((sample_idx, slice(None))).T,
        index=transcript_clean,
        columns=samples,
    )
    ann = build_collapsed_amyloid(root, samples)
    covered = {
        s: set(ann.loc[ann["sample"].eq(s), "Transcript_ID_clean"]).intersection(protein_coding_transcripts)
        for s in samples
    }
    rows = []
    for sample in samples:
        vals = tpm[sample].groupby(level=0).sum()
        protein_coding_tpm = float(vals.loc[vals.index.intersection(protein_coding_transcripts)].sum())
        ann_tpm = float(vals.loc[vals.index.intersection(covered[sample])].sum())
        rows.append(
            {
                "sample": sample,
                "protein_coding_TPM_gencode": protein_coding_tpm,
                "protein_coding_annotated_TPM": ann_tpm,
                "protein_coding_annotation_fraction": ann_tpm / protein_coding_tpm if protein_coding_tpm else np.nan,
            }
        )
    cov = pd.DataFrame(rows).merge(sample_meta, on="sample", how="left")
    cov["keep_coverage_qc"] = cov["protein_coding_annotation_fraction"] >= COVERAGE_MIN
    return cov


def differential_table(root, sample_meta, strict_set, gencode_ann, protein_coding_transcripts, analysis_name, rng):
    sample_ids = [str(x) for x in root["layers/expression/sample_ids"][:]]
    samples = sample_meta["sample"].tolist()
    sample_idx = [sample_ids.index(s) for s in samples]
    transcript_ids = pd.Series(root["layers/expression/transcript_ids"][:]).astype(str)
    transcript_clean = transcript_ids.map(clean_enst)
    counts = pd.DataFrame(
        root["layers/expression/est_counts"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=transcript_ids,
        columns=samples,
    )
    counts["Transcript_ID_clean"] = transcript_clean.values
    counts = counts.dropna(subset=["Transcript_ID_clean"])
    counts = counts.groupby("Transcript_ID_clean", sort=False)[samples].sum()
    counts = counts.loc[counts.index.intersection(protein_coding_transcripts)]

    count_filter = (counts >= MIN_COUNT).sum(axis=1) >= MIN_SAMPLES
    counts = counts.loc[count_filter]
    lib_size = counts.sum(axis=0)
    logcpm = np.log2(counts.div(lib_size, axis=1) * 1e6 + PSEUDOCOUNT)

    ac = sample_meta.loc[sample_meta["Group"].eq("AC"), "sample"].tolist()
    tol = sample_meta.loc[sample_meta["Group"].eq("TOL"), "sample"].tolist()
    mean_ac = logcpm[ac].mean(axis=1)
    mean_tol = logcpm[tol].mean(axis=1)
    log2fc = mean_tol - mean_ac
    p_values = np.array(
        [
            ttest_ind(logcpm.loc[tx, tol], logcpm.loc[tx, ac], equal_var=False, nan_policy="omit").pvalue
            for tx in logcpm.index
        ]
    )
    de = pd.DataFrame(
        {
            "Transcript_ID_clean": logcpm.index,
            "baseMean_counts": counts.mean(axis=1).values,
            "mean_logCPM_AC": mean_ac.values,
            "mean_logCPM_TOL": mean_tol.values,
            "log2FC_TOL_vs_AC": log2fc.values,
            "p_value": p_values,
        }
    )
    de["FDR"] = multipletests(np.nan_to_num(de["p_value"].values, nan=1.0), method="fdr_bh")[1]
    de["rank_stat_signed_neglog10p"] = signed_neglog10_p(de["p_value"], de["log2FC_TOL_vs_AC"])
    de = de.merge(gencode_ann, on="Transcript_ID_clean", how="left")
    de["strict_consensus_amyloid"] = de["Transcript_ID_clean"].isin(strict_set)
    de = de.sort_values("log2FC_TOL_vs_AC", ascending=False).reset_index(drop=True)

    set_values = de.loc[de["strict_consensus_amyloid"], "log2FC_TOL_vs_AC"].to_numpy()
    background_values = de.loc[~de["strict_consensus_amyloid"], "log2FC_TOL_vs_AC"].to_numpy()
    mw = mannwhitneyu(set_values, background_values, alternative="greater")
    auc = mw.statistic / (len(set_values) * len(background_values))
    mean_diff = float(np.mean(set_values) - np.mean(background_values))
    random_means = np.array(
        [
            rng.choice(de["log2FC_TOL_vs_AC"].to_numpy(), size=len(set_values), replace=False).mean()
            for _ in range(N_RANDOM_SET_PERMUTATIONS)
        ]
    )
    mean_ci = bootstrap_ci(set_values, rng)

    ranked = de.sort_values("log2FC_TOL_vs_AC", ascending=False).reset_index(drop=True)
    observed_es = gsea_es(ranked["log2FC_TOL_vs_AC"].to_numpy(), ranked["strict_consensus_amyloid"].to_numpy())
    null_es = []
    labels = sample_meta["Group"].to_numpy().copy()
    for _ in range(N_GSEA_PERMUTATIONS):
        perm = rng.permutation(labels)
        perm_meta = sample_meta.copy()
        perm_meta["perm_group"] = perm
        perm_ac = perm_meta.loc[perm_meta["perm_group"].eq("AC"), "sample"].tolist()
        perm_tol = perm_meta.loc[perm_meta["perm_group"].eq("TOL"), "sample"].tolist()
        perm_fc = logcpm[perm_tol].mean(axis=1) - logcpm[perm_ac].mean(axis=1)
        perm_order = np.argsort(-perm_fc.to_numpy())
        null_es.append(gsea_es(perm_fc.to_numpy()[perm_order], de["strict_consensus_amyloid"].to_numpy()[perm_order]))
    null_es = np.asarray(null_es)
    nes = normalized_es(observed_es, null_es)

    summary = {
        "analysis_set": analysis_name,
        "n_samples": int(len(samples)),
        "n_AC": int(len(ac)),
        "n_TOL": int(len(tol)),
        "n_transcripts_tested": int(len(de)),
        "n_strict_amyloid_transcripts_tested": int(de["strict_consensus_amyloid"].sum()),
        "mean_log2FC_strict_amyloid": float(np.mean(set_values)),
        "mean_log2FC_background": float(np.mean(background_values)),
        "mean_log2FC_difference_amyloid_minus_background": mean_diff,
        "strict_amyloid_mean_log2FC_bootstrap_ci_low": float(mean_ci[0]),
        "strict_amyloid_mean_log2FC_bootstrap_ci_high": float(mean_ci[1]),
        "mannwhitney_one_sided_p_amyloid_gt_background": float(mw.pvalue),
        "auc_probability_amyloid_log2FC_gt_background": float(auc),
        "empirical_p_random_sets_mean_log2FC_greater": float(empirical_p(np.mean(set_values), random_means, "greater")),
        "gsea_ES": float(observed_es),
        "gsea_NES": float(nes),
        "gsea_empirical_p_ES_greater": float(empirical_p(observed_es, null_es, "greater")),
        "n_gsea_label_permutations": int(N_GSEA_PERMUTATIONS),
        "n_random_set_permutations": int(N_RANDOM_SET_PERMUTATIONS),
    }

    return de, ranked, pd.DataFrame([summary]), counts


def main():
    rng = np.random.default_rng(RANDOM_SEED)
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    gencode_ann, gencode_path = load_gencode_annotation()
    protein_coding_transcripts = set(
        gencode_ann.loc[gencode_ann["biotype"].eq("protein_coding"), "Transcript_ID_clean"].dropna().astype(str)
    )
    gencode_ann.to_csv(OUT / "gencode_transcript_annotation_used.tsv", sep="\t", index=False)
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

    transcript_ids = pd.Series(root["layers/expression/transcript_ids"][:]).astype(str)
    transcript_clean = transcript_ids.map(clean_enst)
    coverage = compute_coverage(root, selected["sample"].tolist(), selected, transcript_clean, protein_coding_transcripts)
    coverage.to_csv(OUT / "sample_coverage_qc.tsv", sep="\t", index=False)

    collapsed = build_collapsed_amyloid(root, selected["sample"].tolist())
    collapsed = collapsed.merge(gencode_ann, on="Transcript_ID_clean", how="left")
    collapsed = collapsed[collapsed["biotype"].eq("protein_coding")].copy()
    global_set = (
        collapsed.groupby("Transcript_ID_clean", as_index=False)
        .agg(Consensus_any=("Consensus_collapsed", consensus_priority))
    )
    global_set = global_set.merge(gencode_ann, on="Transcript_ID_clean", how="left")
    global_set["strict_consensus_amyloid"] = global_set["Consensus_any"].eq("Amyloid")
    strict_set = set(global_set.loc[global_set["strict_consensus_amyloid"], "Transcript_ID_clean"])
    global_set.to_csv(OUT / "strict_amyloid_gene_set.tsv", sep="\t", index=False)
    collapsed.to_csv(OUT / "collapsed_sample_transcript_amyloid_annotations.tsv", sep="\t", index=False)

    summaries = []
    for name, meta in [
        ("all_BC_samples", selected),
        ("coverage_QC_protein_coding_annotation_fraction_ge_0_5", selected[selected["sample"].isin(coverage.loc[coverage["keep_coverage_qc"], "sample"])]),
    ]:
        de, ranked, summary, counts = differential_table(
            root,
            meta.reset_index(drop=True),
            strict_set,
            gencode_ann,
            protein_coding_transcripts,
            name,
            rng,
        )
        method_dir = OUT / name
        method_dir.mkdir(exist_ok=True)
        de.to_csv(method_dir / "transcript_logCPM_DE_ranked_by_log2FC.tsv", sep="\t", index=False)
        ranked.to_csv(method_dir / "transcript_ranked_for_gsea.tsv", sep="\t", index=False)
        counts.to_csv(method_dir / "counts_matrix_transcript_collapsed.tsv", sep="\t")
        meta.to_csv(method_dir / "sample_metadata.tsv", sep="\t", index=False)
        summary.to_csv(method_dir / "competitive_gene_set_summary.tsv", sep="\t", index=False)
        summaries.append(summary)

    summary = pd.concat(summaries, ignore_index=True)
    summary.to_csv(OUT / "competitive_gene_set_summary.tsv", sep="\t", index=False)

    report = {
        "method": "Competitive gene-set test using strict amyloid transcripts as set. Transcripts ranked by log2FC TOL/chili vs AC/control from counts-derived logCPM; competitive Mann-Whitney/AUC and random-set empirical p-value, plus phenotype-permutation GSEA.",
        "annotation_path": str(ANNOTATION_PATH),
        "gencode_annotation_path": str(gencode_path),
        "universe": "GENCODE v48 biotype == protein_coding transcripts intersecting the quantified transcript count matrix and count filter.",
        "selected_BC_samples": selected[["sample", "subject", "time", "raw_Group", "Group", "batch", "AGE", "sex", "treatment"]].to_dict("records"),
        "coverage_qc": coverage[["sample", "Group", "protein_coding_annotation_fraction", "keep_coverage_qc"]].to_dict("records"),
        "excluded_by_coverage_qc": coverage.loc[~coverage["keep_coverage_qc"], ["sample", "Group", "protein_coding_annotation_fraction"]].to_dict("records"),
        "n_strict_amyloid_transcripts_in_global_set": int(global_set["strict_consensus_amyloid"].sum()),
        "results": summary.to_dict("records"),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")

    lines = ["# Competitive Gene-Set Test: strict amyloid transcripts in baseline ChILI", ""]
    lines.append("## What was done")
    lines.append(
        "Recomputed from `project.zarr` using only baseline (`time == \"BC\"`) blood RNA-seq samples. "
        "Metadata `Group=control` was mapped to AC and `Group=chili` to TOL. "
        f"Transcript universe was restricted using local GENCODE annotation `{gencode_path}` to `biotype == \"protein_coding\"`. "
        "Strict amyloid transcripts were defined by collapsed `sample x transcript` amyloid annotations with "
        "`Consensus_collapsed == \"Amyloid\"`; protein variants mapping to the same transcript were collapsed to avoid H1/H2 double counting."
    )
    lines.append("")
    lines.append("Transcripts were ranked by counts-derived logCPM `log2FC_TOL_vs_AC`. The competitive test asks whether strict amyloid transcripts have systematically higher log2FC than the transcript background. I report Mann-Whitney/AUC, random-set empirical p-value, and phenotype-permutation GSEA.")
    lines.append("")
    lines.append("## Samples")
    for _, row in coverage.iterrows():
        flag = "kept" if row["keep_coverage_qc"] else "excluded from coverage-QC"
        lines.append(f"- {row['sample']}: {row['Group']}, subject {row['subject']}, coverage={row['protein_coding_annotation_fraction']:.3f}, {flag}")
    lines.append("")
    lines.append("## Results")
    for _, row in summary.iterrows():
        lines.append(
            f"- {row['analysis_set']}: n={int(row['n_samples'])} "
            f"(AC={int(row['n_AC'])}, TOL={int(row['n_TOL'])}); "
            f"strict amyloid tested={int(row['n_strict_amyloid_transcripts_tested'])}/{int(row['n_transcripts_tested'])}. "
            f"Mean log2FC amyloid={row['mean_log2FC_strict_amyloid']:.4f} "
            f"(95% bootstrap CI {row['strict_amyloid_mean_log2FC_bootstrap_ci_low']:.4f} to {row['strict_amyloid_mean_log2FC_bootstrap_ci_high']:.4f}); "
            f"background={row['mean_log2FC_background']:.4f}; "
            f"delta={row['mean_log2FC_difference_amyloid_minus_background']:.4f}; "
            f"AUC={row['auc_probability_amyloid_log2FC_gt_background']:.3f}; "
            f"MW one-sided p={row['mannwhitney_one_sided_p_amyloid_gt_background']:.3g}; "
            f"random-set empirical p={row['empirical_p_random_sets_mean_log2FC_greater']:.3g}; "
            f"GSEA ES={row['gsea_ES']:.4f}, NES={row['gsea_NES']:.3f}, empirical p={row['gsea_empirical_p_ES_greater']:.3g}."
        )
    lines.append("")
    lines.append("## Interpretation")
    best = summary.iloc[-1]
    if (
        best["mean_log2FC_difference_amyloid_minus_background"] > 0
        and best["auc_probability_amyloid_log2FC_gt_background"] > 0.5
        and best["empirical_p_random_sets_mean_log2FC_greater"] < 0.1
    ):
        lines.append("The coverage-QC analysis is directionally consistent with enrichment of strict amyloid transcripts among genes upregulated in ChILI/TOL, but this should be treated as discovery evidence rather than validation.")
    else:
        lines.append("The coverage-QC competitive gene-set result does not provide strong support that strict amyloid transcripts are systematically higher among genes upregulated in ChILI/TOL.")
    lines.append("")
    lines.append("## Limitations")
    lines.append("Small sample size, heterogeneous treatment/batch composition, and transcript-level amyloid annotation coverage limit inference. The permutation p-values are empirical and constrained by the number of possible label rearrangements in these small groups.")
    lines.append("")
    lines.append("## Next best step")
    lines.append("Run a counts-based model with design covariates where feasible, then validate the amyloid gene-set direction in an independent baseline cohort or by targeted qPCR/protein-level follow-up for the top contributing strict amyloid transcripts.")
    (OUT / "COMPETITIVE_GENE_SET_TEST_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(summary.to_string(index=False))
    print(f"Outputs: {OUT}")


if __name__ == "__main__":
    main()
