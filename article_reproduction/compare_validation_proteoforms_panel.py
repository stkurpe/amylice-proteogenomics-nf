import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parent
ZARR_PATH = Path("/Users/user841/Projects/BIoinfServer/Results_analysis/bioinfo_zarr/validation_project.zarr")
ANNOTATION_PATH = ROOT / "GSE287540_SraRunTable.csv"
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "validation" / "proteoform_panel"
GENCODE_PATH = (
    ROOT
    / "analysis_outputs"
    / "amyloid_bc_control_vs_chili"
    / "competitive_gene_set_test_gencode_v48_protein_coding"
    / "gencode_transcript_annotation_used.tsv"
)
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
TIMEPOINT = "BC"
RANDOM_SEED = 20260606
N_BOOTSTRAPS = 10_000


PROTEOFORMS = [
    ("RNF5", "H2_ENST00000375094.4"),
    ("LTB", "H2_ENST00000446745.2"),
    ("AP2B1", "H2_ENST00000590432.5"),
    ("AIMP1", "H2_ENST00000510207.5"),
    ("SLAIN2", "H2_ENST00000512093.5"),
    ("CERT1", "H2_ENST00000508809.2"),
    ("SEPTIN7", "H2_ENST00000485569.1"),
    ("MRPS21", "H2_ENST00000614145.5"),
    ("SLAIN2", "H1_ENST00000512093.5"),
    ("TRAF3IP3", "H1_ENST00000367023.5"),
    ("PPP2R5C", "H2_ENST00000561006.5"),
    ("ZNF302", "H2_ENST00000507959.5"),
    ("STK17B", "H2_ENST00000714423.1"),
    ("RNF215", "H2_ENST00000431544.2"),
    ("TRAPPC2L", "H2_ENST00000567895.5"),
    ("DGKD", "H1_ENST00000442524.4"),
    ("ERAL1", "H2_ENST00000583487.1"),
    ("USP31", "H2_ENST00000563525.2"),
    ("RPS20", "H2_ENST00000519807.5"),
    ("PPP4R2", "H1_ENST00000710398.1"),
    ("PPP4R2", "H1_ENST00000356692.10"),
    ("KMT2B", "H2_ENST00000686920.1"),
    ("ZNF672", "H1_ENST00000423362.1"),
    ("ANKRD49", "H2_ENST00000540349.1"),
    ("TATDN3", "H1_ENST00000527693.1"),
    ("YLPM1", "H2_ENST00000554107.2"),
    ("CCDC91", "H2_ENST00000536154.5"),
    ("ZNF580", "H2_ENST00000545125.1"),
    ("TFDP2", "H2_ENST00000467634.1"),
    ("PSMD8", "H2_ENST00000591250.1"),
]


def clean_enst(x):
    if pd.isna(x):
        return np.nan
    match = re.search(r"ENST\d+(?:\.\d+)?", str(x))
    return match.group(0).split(".")[0] if match else np.nan


def zarr_group_to_df(group):
    return pd.DataFrame({col: group[col][:] for col in group.array_keys()})


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


def bootstrap_ci(ac, tol, rng):
    ac = np.asarray(ac, dtype=float)
    tol = np.asarray(tol, dtype=float)
    if len(ac) == 0 or len(tol) == 0:
        return [np.nan, np.nan, np.nan]
    vals = np.empty(N_BOOTSTRAPS, dtype=float)
    for i in range(N_BOOTSTRAPS):
        vals[i] = rng.choice(tol, len(tol), replace=True).mean() - rng.choice(ac, len(ac), replace=True).mean()
    return np.percentile(vals, [2.5, 50, 97.5])


def compare_metric(df, metric, rng):
    rows = []
    for (gene, protein_id), g in df.groupby(["gene_symbol_input", "protein_id"], sort=False):
        ac = g.loc[g["Group"].eq("AC") & g[metric].notna(), metric].astype(float).to_numpy()
        tol = g.loc[g["Group"].eq("TOL") & g[metric].notna(), metric].astype(float).to_numpy()
        ci = bootstrap_ci(ac, tol, rng)
        p = np.nan
        if len(ac) and len(tol):
            try:
                p = mannwhitneyu(tol, ac, alternative="greater").pvalue
            except ValueError:
                p = np.nan
        rows.append(
            {
                "gene_symbol_input": gene,
                "protein_id": protein_id,
                "metric": metric,
                "n_AC": len(ac),
                "n_TOL": len(tol),
                "mean_AC": float(np.mean(ac)) if len(ac) else np.nan,
                "mean_TOL": float(np.mean(tol)) if len(tol) else np.nan,
                "median_AC": float(np.median(ac)) if len(ac) else np.nan,
                "median_TOL": float(np.median(tol)) if len(tol) else np.nan,
                "delta_TOL_minus_AC": float(np.mean(tol) - np.mean(ac)) if len(ac) and len(tol) else np.nan,
                "mannwhitney_one_sided_TOL_gt_AC_p": float(p) if not pd.isna(p) else np.nan,
                "cliffs_delta_TOL_vs_AC": float(cliffs_delta(tol, ac)) if len(ac) and len(tol) else np.nan,
                "bootstrap_delta_ci2_5": float(ci[0]),
                "bootstrap_delta_median": float(ci[1]),
                "bootstrap_delta_ci97_5": float(ci[2]),
            }
        )
    return pd.DataFrame(rows)


def markdown_table(df):
    display = df.copy()
    for col in display.columns:
        if pd.api.types.is_float_dtype(display[col]):
            display[col] = display[col].map(lambda x: "" if pd.isna(x) else f"{x:.6g}")
    return display.to_markdown(index=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)
    root = zarr.open_group(str(ZARR_PATH), mode="r")

    requested = pd.DataFrame(PROTEOFORMS, columns=["gene_symbol_input", "protein_id"])
    requested["Transcript_ID_clean"] = requested["protein_id"].map(clean_enst)
    gencode = pd.read_csv(GENCODE_PATH, sep="\t")
    gencode = gencode[["Transcript_ID_clean", "gene_id", "gene_name", "biotype"]].drop_duplicates("Transcript_ID_clean")
    requested = requested.merge(gencode, on="Transcript_ID_clean", how="left")
    requested["keep_protein_coding"] = requested["biotype"].eq("protein_coding")
    requested.to_csv(OUT / "requested_proteoforms_gencode_filter.tsv", sep="\t", index=False)
    panel = requested[requested["keep_protein_coding"]].copy()

    meta = pd.read_csv(ANNOTATION_PATH)
    meta["Run"] = meta["Run"].astype(str)
    meta["raw_Group"] = meta["Group"].astype(str)
    meta["Group"] = meta["raw_Group"].map(GROUP_MAP)
    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    sample_meta = meta[
        meta["Run"].isin(zarr_samples)
        & meta["time"].eq(TIMEPOINT)
        & meta["Group"].isin(["AC", "TOL"])
    ].rename(columns={"Run": "sample"})
    sample_meta = sample_meta.set_index("sample").loc[zarr_samples].reset_index()

    pf = pd.concat([zarr_group_to_df(root[f"samples/{sample}/protein_features"]) for sample in zarr_samples], ignore_index=True)
    pf["Transcript_ID_clean"] = pf["protein_id"].map(clean_enst)
    pf["abs_charge_density"] = pf["charge_density"].abs()
    pf["pI_distance_from_7"] = (pf["pI"] - 7).abs()
    for col in ["beta_propensity", "KD_max", "aromaticity", "abs_charge_density", "pI_distance_from_7"]:
        pf[f"{col}_z"] = zscore(pf[col])
    pf["Amyloid_Integrative_Index_validation"] = (
        pf["beta_propensity_z"]
        + pf["KD_max_z"]
        + pf["aromaticity_z"]
        - pf["abs_charge_density_z"]
        - pf["pI_distance_from_7_z"]
    ) / 5

    panel_long = sample_meta[["sample", "subject", "Group", "raw_Group", "time", "treatment", "batch"]].merge(
        panel[["gene_symbol_input", "protein_id", "Transcript_ID_clean", "gene_id", "gene_name", "biotype"]],
        how="cross",
    )
    panel_long = panel_long.merge(
        pf[
            [
                "sample",
                "protein_id",
                "protein_length",
                "Amyloid_Integrative_Index_validation",
                "beta_propensity",
                "KD_max",
                "aromaticity",
                "charge_density",
                "pI",
            ]
        ],
        on=["sample", "protein_id"],
        how="left",
    )
    panel_long["proteoform_present"] = panel_long["Amyloid_Integrative_Index_validation"].notna()

    sample_idx = [zarr_samples.index(s) for s in zarr_samples]
    tx_clean = pd.Series(root["layers/expression/transcript_ids"][:]).astype(str).map(clean_enst)
    tpm = pd.DataFrame(
        root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=tx_clean,
        columns=zarr_samples,
    )
    tpm = tpm[~tpm.index.isna()].groupby(level=0).sum()
    tpm_long = tpm.loc[tpm.index.intersection(panel["Transcript_ID_clean"])].reset_index().melt(
        id_vars="index", var_name="sample", value_name="Transcript_TPM"
    )
    tpm_long = tpm_long.rename(columns={"index": "Transcript_ID_clean"})
    panel_long = panel_long.merge(tpm_long, on=["sample", "Transcript_ID_clean"], how="left")
    panel_long["Transcript_TPM"] = panel_long["Transcript_TPM"].fillna(0.0)
    panel_long.to_csv(OUT / "validation_proteoform_panel_long.tsv", sep="\t", index=False)

    presence = (
        panel_long.groupby(["gene_symbol_input", "protein_id", "Group"], as_index=False)
        .agg(
            n_samples=("sample", "nunique"),
            n_present=("proteoform_present", "sum"),
            present_fraction=("proteoform_present", "mean"),
        )
        .pivot(index=["gene_symbol_input", "protein_id"], columns="Group")
    )
    presence.columns = ["_".join(col).strip("_") for col in presence.columns.to_flat_index()]
    presence = presence.reset_index()
    presence["present_fraction_delta_TOL_minus_AC"] = presence.get("present_fraction_TOL", np.nan) - presence.get(
        "present_fraction_AC", np.nan
    )

    index_stats = compare_metric(panel_long, "Amyloid_Integrative_Index_validation", rng)
    tpm_stats = compare_metric(panel_long, "Transcript_TPM", rng)
    stats = pd.concat([index_stats, tpm_stats], ignore_index=True)
    stats = stats.merge(presence, on=["gene_symbol_input", "protein_id"], how="left")
    stats.to_csv(OUT / "validation_proteoform_panel_group_statistics.tsv", sep="\t", index=False)

    summary = {
        "zarr_path": str(ZARR_PATH),
        "out_dir": str(OUT),
        "n_requested_proteoforms": int(len(requested)),
        "n_protein_coding_retained": int(len(panel)),
        "n_excluded_non_protein_coding_or_unannotated": int((~requested["keep_protein_coding"]).sum()),
        "random_seed": RANDOM_SEED,
        "n_bootstraps": N_BOOTSTRAPS,
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    index_view = stats[stats["metric"].eq("Amyloid_Integrative_Index_validation")].copy()
    tpm_view = stats[stats["metric"].eq("Transcript_TPM")].copy()
    report = [
        "# Validation proteoform panel comparison",
        "",
        "Only GENCODE `protein_coding` transcripts were retained before comparison.",
        "",
        "## GENCODE Filter",
        "",
        markdown_table(
            requested[
                [
                    "gene_symbol_input",
                    "protein_id",
                    "Transcript_ID_clean",
                    "gene_name",
                    "biotype",
                    "keep_protein_coding",
                ]
            ]
        ),
        "",
        "## Amyloid Integrative Index",
        "",
        markdown_table(
            index_view[
                [
                    "gene_symbol_input",
                    "protein_id",
                    "n_AC",
                    "n_TOL",
                    "mean_AC",
                    "mean_TOL",
                    "delta_TOL_minus_AC",
                    "mannwhitney_one_sided_TOL_gt_AC_p",
                    "cliffs_delta_TOL_vs_AC",
                    "present_fraction_AC",
                    "present_fraction_TOL",
                ]
            ]
        ),
        "",
        "## Transcript TPM",
        "",
        markdown_table(
            tpm_view[
                [
                    "gene_symbol_input",
                    "protein_id",
                    "n_AC",
                    "n_TOL",
                    "mean_AC",
                    "mean_TOL",
                    "delta_TOL_minus_AC",
                    "mannwhitney_one_sided_TOL_gt_AC_p",
                    "cliffs_delta_TOL_vs_AC",
                ]
            ]
        ),
    ]
    (OUT / "VALIDATION_PROTEOFORM_PANEL_REPORT.md").write_text("\n".join(report) + "\n")


if __name__ == "__main__":
    main()
