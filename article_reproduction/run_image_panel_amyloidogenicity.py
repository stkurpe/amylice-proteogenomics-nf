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
GENCODE_PATH = (
    ROOT
    / "analysis_outputs"
    / "amyloid_bc_control_vs_chili"
    / "rank_based_amyloid_burden_signature"
    / "gencode_v48_transcript_annotation.tsv"
)
OUT = ROOT / "analysis_outputs" / "amyloid_bc_control_vs_chili" / "image_mixedlm_feature_panel_amyloidogenicity"
OUT.mkdir(parents=True, exist_ok=True)

TIMEPOINT = "BC"
GROUP_MAP = {"control": "AC", "chili": "TOL", "chill": "TOL"}
COVERAGE_MIN = 0.5


# Gene/protein labels read from the supplied figure.
IMAGE_PANEL = [
    ("ZNF580", "H2_ENST00000545125.1"),
    ("SEPTIN7", "H2_ENST00000485569.1"),
    ("TRAF3IP3", "H1_ENST00000367023.5"),
    ("CERT1", "H2_ENST00000508809.2"),
    ("PPP2R5C", "H2_ENST00000555237.5"),
    ("MRPS21", "H2_ENST00000614145.5"),
    ("ATP8A1", "H1_ENST00000510289.1"),
    ("GAR1", "H1_ENST00000394631.7"),
    ("GFUS", "H2_ENST00000524719.1"),
    ("MTURN", "H2_ENST00000324489.5"),
    ("CC2D2A", "H1_ENST00000652443.1"),
    ("CIR1", "H2_ENST00000377973.4"),
    ("SPATS2L", "H1_ENST00000366118.2"),
    ("RARB", "H1_ENST00000692640.1"),
    ("NLRC5", "H2_ENST00000399221.3"),
    ("SLAIN2", "H1_ENST00000512093.5"),
    ("SLAIN2", "H2_ENST00000512093.5"),
    ("STK17B", "H2_ENST00000714423.1"),
    ("YLPM1", "H2_ENST00000554107.2"),
    ("ZNF672", "H1_ENST00000423362.1"),
    ("PSMD8", "H2_ENST00000591250.1"),
    ("PTPN2", "H2_ENST00000591901.5"),
    ("S100P", "H1_ENST00000296370.4"),
    ("PNRC2", "H1_ENST00000334351.8"),
    ("NBRF9", "H1_ENST00000621074.5"),
    ("GOLPH3L", "H1_ENST00000271732.8"),
    ("KIAA0040", "H1_ENST00000619513.1"),
    ("STX6", "H1_ENST00000258301.6"),
    ("TACC1", "H2_ENST00000520611.1"),
    ("AP2B1", "H2_ENST00000590432.5"),
    ("ZNF616", "H2_ENST00000600282.1"),
    ("ZNF302", "H2_ENST00000507959.5"),
    ("KIAA1191", "H2_ENST00000504688.5"),
    ("CEP41", "H2_ENST00000489512.5"),
    ("EPC2", "H2_ENST00000409654.5"),
    ("CCDC169", "H2_ENST00000379864.6"),
    ("SHE", "H1_ENST00000555188.5"),
    ("MCUB", "H1_ENST00000394650.7"),
    ("CLNS1A", "H2_ENST00000525428.6"),
    ("SPI1", "H2_ENST00000713542.1"),
    ("MTHD", "H1_ENST00000521933.1"),
    ("ILF2", "H1_ENST00000361891.9"),
    ("SNX10", "H2_ENST00000409838.1"),
    ("TIMM10", "H2_ENST00000257245.9"),
    ("PIK3C2A", "H2_ENST00000532035.1"),
    ("S100A8", "H1_ENST00000368733.4"),
    ("LGMN", "H2_ENST00000554080.5"),
    ("POLD4", "H2_ENST00000529704.5"),
    ("TIA1", "H2_ENST00000415783.6"),
    ("TIA1", "H1_ENST00000415783.6"),
    ("HNRNPD", "H1_ENST00000313899.12"),
    ("RPUSD3", "H2_ENST00000433555.1"),
    ("HBP1", "H2_ENST00000479011.1"),
    ("AIMP1", "H1_ENST00000510207.5"),
    ("TMA16", "H1_ENST00000513272.5"),
    ("RP9", "H1_ENST00000684207.1"),
    ("METTL6", "H1_ENST00000453846.5"),
    ("EIF4H", "H2_ENST00000677681.1"),
    ("BCL7B", "H2_ENST00000455335.2"),
    ("TNRC18", "H2_ENST00000399434.2"),
    ("STMP1", "H2_ENST00000507606.3"),
    ("FXYD5", "H2_ENST00000392217.3"),
    ("CD59", "H2_ENST00000652678.1"),
    ("GTF2A2", "H2_ENST00000396060.7"),
    ("CD59", "H2_ENST00000395850.9"),
    ("VPS53", "H2_ENST00000571456.2"),
    ("PLIN5", "H2_ENST00000592610.1"),
    ("ANAPC11", "H2_ENST00000582222.5"),
    ("NAA38", "H2_ENST00000576384.1"),
    ("CLEC2D", "H2_ENST00000290855.11"),
]


def clean_enst(x):
    if pd.isna(x):
        return np.nan
    m = re.search(r"ENST\d+(?:\.\d+)?", str(x))
    return m.group(0).split(".")[0] if m else np.nan


def zarr_group_to_df(g):
    return pd.DataFrame({col: g[col][:] for col in g.array_keys()})


def zscore_series(s):
    s = pd.Series(s, dtype=float)
    sd = s.std(skipna=True, ddof=0)
    if sd == 0 or pd.isna(sd):
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


def compare(df, metric, label):
    d = df.dropna(subset=[metric, "Group"])
    ac = d.loc[d["Group"].eq("AC"), metric].astype(float)
    tol = d.loc[d["Group"].eq("TOL"), metric].astype(float)
    if len(ac) < 2 or len(tol) < 2:
        return None
    return {
        "analysis_set": label,
        "metric": metric,
        "n_AC": int(len(ac)),
        "n_TOL": int(len(tol)),
        "mean_AC": float(ac.mean()),
        "mean_TOL": float(tol.mean()),
        "median_AC": float(ac.median()),
        "median_TOL": float(tol.median()),
        "delta_TOL_minus_AC": float(tol.mean() - ac.mean()),
        "cliffs_delta_TOL_vs_AC": cliffs_delta(tol, ac),
        "mannwhitney_two_sided_p": float(mannwhitneyu(tol, ac, alternative="two-sided").pvalue),
        "mannwhitney_one_sided_p_TOL_gt_AC": float(mannwhitneyu(tol, ac, alternative="greater").pvalue),
        "welch_two_sided_p": float(ttest_ind(tol, ac, equal_var=False, nan_policy="omit").pvalue),
    }


def load_metadata(root):
    meta = pd.read_csv(ANNOTATION_PATH)
    meta["raw_Group"] = meta["Group"].astype(str)
    meta["Group"] = meta["raw_Group"].map(GROUP_MAP)
    meta = meta.rename(columns={"Run": "sample"})
    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    selected = meta[
        meta["sample"].isin(zarr_samples)
        & meta["time"].eq(TIMEPOINT)
        & meta["Group"].isin(["AC", "TOL"])
    ].copy()
    return selected.set_index("sample").loc[[s for s in zarr_samples if s in set(selected["sample"])]].reset_index()


def load_tpm(root, samples):
    zarr_samples = [str(x) for x in root["layers/expression/sample_ids"][:]]
    sample_idx = [zarr_samples.index(s) for s in samples]
    tids = pd.Series(root["layers/expression/transcript_ids"][:]).astype(str)
    clean = tids.map(clean_enst)
    tpm = pd.DataFrame(
        root["layers/expression/tpm"].get_orthogonal_selection((sample_idx, slice(None))).T,
        index=clean,
        columns=samples,
    )
    tpm = tpm[~pd.isna(tpm.index)]
    return tpm.groupby(level=0, sort=False).sum()


def build_combined(root, selected, panel):
    rows = []
    for sample in selected["sample"]:
        amy = zarr_group_to_df(root[f"samples/{sample}/amyloid"])
        pf = zarr_group_to_df(root[f"samples/{sample}/protein_features"])
        amy["ID"] = amy["sample"].astype(str) + "|" + amy["Sequence_ID"].astype(str)
        pf["ID"] = pf["sample"].astype(str) + "|" + pf["protein_id"].astype(str)
        combined = pf.merge(amy.drop(columns=["sample"], errors="ignore"), on="ID", how="left")
        combined["sample"] = sample
        rows.append(combined)
    combined = pd.concat(rows, ignore_index=True)
    combined["protein_id_no_version"] = combined["protein_id"].astype(str).str.replace(r"\\.\\d+$", "", regex=True)
    combined["Transcript_ID_clean"] = combined["protein_id"].map(clean_enst)
    combined["abs_charge_density"] = combined["charge_density"].abs()
    combined["pI_distance_from_7"] = (combined["pI"] - 7).abs()
    for col in ["beta_propensity", "KD_max", "aromaticity", "abs_charge_density", "pI_distance_from_7"]:
        combined[f"{col}_z"] = zscore_series(combined[col])
    combined["Amyloid_Integrative_Index"] = (
        combined["beta_propensity_z"]
        + combined["KD_max_z"]
        + combined["aromaticity_z"]
        - combined["abs_charge_density_z"]
        - combined["pI_distance_from_7_z"]
    ) / 5
    combined = combined.merge(selected[["sample", "subject", "Group", "raw_Group", "batch", "AGE", "sex"]], on="sample")
    panel_ids = set(panel["protein_id_no_version"])
    return combined[combined["protein_id_no_version"].isin(panel_ids)].copy()


def compute_coverage(tpm, selected, combined, protein_coding):
    rows = []
    for sample in selected["sample"]:
        vals = tpm[sample]
        pc_tpm = float(vals.loc[vals.index.intersection(protein_coding)].sum())
        covered = set(combined.loc[combined["sample"].eq(sample), "Transcript_ID_clean"].dropna())
        ann_tpm = float(vals.loc[vals.index.intersection(covered)].sum())
        rows.append(
            {
                "sample": sample,
                "protein_coding_TPM_gencode": pc_tpm,
                "panel_annotated_TPM": ann_tpm,
                "panel_annotation_fraction_of_protein_coding_TPM": ann_tpm / pc_tpm if pc_tpm else np.nan,
            }
        )
    return pd.DataFrame(rows).merge(selected, on="sample", how="left")


def main():
    root = zarr.open_group(str(ZARR_PATH), mode="r")
    selected = load_metadata(root)
    gencode = pd.read_csv(GENCODE_PATH, sep="\t")
    protein_coding = set(gencode.loc[gencode["biotype"].eq("protein_coding"), "Transcript_ID_clean"].astype(str))
    panel = pd.DataFrame(IMAGE_PANEL, columns=["image_gene", "protein_id"])
    panel["protein_id_no_version"] = panel["protein_id"].str.replace(r"\\.\\d+$", "", regex=True)
    panel["Transcript_ID_clean"] = panel["protein_id"].map(clean_enst)
    panel = panel.drop_duplicates(["protein_id_no_version"])
    panel = panel.merge(gencode, on="Transcript_ID_clean", how="left")
    panel["is_gencode_protein_coding"] = panel["biotype"].eq("protein_coding")

    tpm = load_tpm(root, selected["sample"].tolist())
    combined = build_combined(root, selected, panel)
    combined = combined.merge(
        panel[["protein_id_no_version", "image_gene", "gene_name", "gene_id", "biotype", "is_gencode_protein_coding"]],
        on="protein_id_no_version",
        how="left",
    )
    combined = combined[combined["is_gencode_protein_coding"]].copy()
    combined["TPM"] = [tpm.loc[tx, s] if tx in tpm.index else np.nan for tx, s in zip(combined["Transcript_ID_clean"], combined["sample"])]
    combined["strict_amyloid"] = combined["Consensus"].eq("Amyloid")
    combined["expression_weighted_amyloid_index"] = combined["TPM"].fillna(0) * combined["Amyloid_Integrative_Index"]

    coverage = compute_coverage(tpm, selected, combined, protein_coding)
    previous_cov = pd.read_csv(
        ROOT
        / "analysis_outputs"
        / "amyloid_bc_control_vs_chili"
        / "extreme_contributor_concentration_index"
        / "sample_coverage_qc.tsv",
        sep="\t",
    )[["sample", "protein_coding_annotation_fraction", "keep_coverage_qc"]]
    coverage = coverage.merge(previous_cov, on="sample", how="left")

    sample_summary = (
        combined.groupby("sample", as_index=False)
        .agg(
            n_image_panel_proteins_observed=("protein_id", "nunique"),
            n_image_panel_transcripts_observed=("Transcript_ID_clean", "nunique"),
            mean_Amyloid_Integrative_Index=("Amyloid_Integrative_Index", "mean"),
            median_Amyloid_Integrative_Index=("Amyloid_Integrative_Index", "median"),
            mean_AmyloGram_Prob=("AmyloGram_Prob", "mean"),
            mean_AMYPred_Prob=("AMYPred_Prob", "mean"),
            strict_amyloid_fraction=("strict_amyloid", "mean"),
            expression_weighted_amyloid_index_sum=("expression_weighted_amyloid_index", "sum"),
            panel_TPM=("TPM", "sum"),
            mean_beta_propensity=("beta_propensity", "mean"),
            mean_KD_max=("KD_max", "mean"),
            mean_aromaticity=("aromaticity", "mean"),
            mean_abs_charge_density=("abs_charge_density", "mean"),
            mean_pI_distance_from_7=("pI_distance_from_7", "mean"),
            mean_chameleon_score=("chameleon_score", "mean"),
        )
    ).merge(selected, on="sample", how="left")
    sample_summary = sample_summary.merge(
        previous_cov[["sample", "protein_coding_annotation_fraction", "keep_coverage_qc"]], on="sample", how="left"
    )
    sample_summary["expression_weighted_amyloid_index_per_panel_TPM"] = (
        sample_summary["expression_weighted_amyloid_index_sum"] / sample_summary["panel_TPM"]
    )

    protein_summary = (
        combined.groupby(["protein_id_no_version", "image_gene", "gene_name", "Transcript_ID_clean"], as_index=False)
        .agg(
            n_samples_observed=("sample", "nunique"),
            mean_AC_Amyloid_Index=("Amyloid_Integrative_Index", lambda x: x[combined.loc[x.index, "Group"].eq("AC")].mean()),
            mean_TOL_Amyloid_Index=("Amyloid_Integrative_Index", lambda x: x[combined.loc[x.index, "Group"].eq("TOL")].mean()),
            mean_AC_AmyloGram_Prob=("AmyloGram_Prob", lambda x: x[combined.loc[x.index, "Group"].eq("AC")].mean()),
            mean_TOL_AmyloGram_Prob=("AmyloGram_Prob", lambda x: x[combined.loc[x.index, "Group"].eq("TOL")].mean()),
            mean_AC_AMYPred_Prob=("AMYPred_Prob", lambda x: x[combined.loc[x.index, "Group"].eq("AC")].mean()),
            mean_TOL_AMYPred_Prob=("AMYPred_Prob", lambda x: x[combined.loc[x.index, "Group"].eq("TOL")].mean()),
            strict_amyloid_fraction_AC=("strict_amyloid", lambda x: x[combined.loc[x.index, "Group"].eq("AC")].mean()),
            strict_amyloid_fraction_TOL=("strict_amyloid", lambda x: x[combined.loc[x.index, "Group"].eq("TOL")].mean()),
            mean_AC_TPM=("TPM", lambda x: x[combined.loc[x.index, "Group"].eq("AC")].mean()),
            mean_TOL_TPM=("TPM", lambda x: x[combined.loc[x.index, "Group"].eq("TOL")].mean()),
        )
    )
    for col in ["Amyloid_Index", "AmyloGram_Prob", "AMYPred_Prob", "TPM"]:
        protein_summary[f"delta_TOL_minus_AC_{col}"] = protein_summary[f"mean_TOL_{col}"] - protein_summary[f"mean_AC_{col}"]
    protein_summary["delta_TOL_minus_AC_strict_amyloid_fraction"] = (
        protein_summary["strict_amyloid_fraction_TOL"] - protein_summary["strict_amyloid_fraction_AC"]
    )

    metrics = [
        "mean_Amyloid_Integrative_Index",
        "mean_AmyloGram_Prob",
        "mean_AMYPred_Prob",
        "strict_amyloid_fraction",
        "expression_weighted_amyloid_index_per_panel_TPM",
        "expression_weighted_amyloid_index_sum",
        "mean_beta_propensity",
        "mean_KD_max",
        "mean_aromaticity",
        "mean_abs_charge_density",
        "mean_pI_distance_from_7",
        "mean_chameleon_score",
    ]
    stats = []
    for label, df in [
        ("all_BC_samples", sample_summary),
        (
            "coverage_QC_protein_coding_annotation_fraction_ge_0_5",
            sample_summary[sample_summary["keep_coverage_qc"]].copy(),
        ),
    ]:
        method_dir = OUT / label
        method_dir.mkdir(exist_ok=True)
        df.to_csv(method_dir / "sample_image_panel_amyloidogenicity.tsv", sep="\t", index=False)
        for metric in metrics:
            row = compare(df, metric, label)
            if row:
                stats.append(row)
    stats = pd.DataFrame(stats)

    panel.to_csv(OUT / "image_panel_genes_proteins_from_figure.tsv", sep="\t", index=False)
    combined.to_csv(OUT / "image_panel_sample_protein_long.tsv", sep="\t", index=False)
    sample_summary.to_csv(OUT / "sample_image_panel_amyloidogenicity.tsv", sep="\t", index=False)
    protein_summary.to_csv(OUT / "protein_image_panel_amyloidogenicity_summary.tsv", sep="\t", index=False)
    stats.to_csv(OUT / "image_panel_group_statistics.tsv", sep="\t", index=False)
    coverage.to_csv(OUT / "sample_panel_coverage.tsv", sep="\t", index=False)

    primary = stats[
        stats["metric"].isin(
            [
                "mean_Amyloid_Integrative_Index",
                "mean_AmyloGram_Prob",
                "mean_AMYPred_Prob",
                "strict_amyloid_fraction",
                "expression_weighted_amyloid_index_per_panel_TPM",
            ]
        )
    ].copy()
    report = {
        "method": "Exploratory image-derived protein panel amyloidogenicity analysis on BC AC/control vs TOL/chili. Protein IDs were read from the supplied figure; GENCODE protein_coding only; H1/H2 protein IDs retained for protein-level panel, with transcript IDs also reported.",
        "n_image_labels": len(IMAGE_PANEL),
        "n_unique_panel_proteins": int(panel["protein_id_no_version"].nunique()),
        "n_unique_panel_proteins_protein_coding": int(panel["is_gencode_protein_coding"].sum()),
        "primary_results": primary.to_dict("records"),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")

    lines = [
        "# Image MixedLM Feature Panel: Amyloidogenicity",
        "",
        "## What was done",
        "Used the gene/protein IDs visible in the supplied MixedLM feature plot. The analysis was restricted to BC samples, AC/control vs TOL/chili, and to GENCODE protein-coding transcripts. H1/H2 protein IDs from the image were retained at protein level; transcript IDs are included so H1/H2 duplicates can be inspected.",
        "",
        "Amyloidogenicity was summarized with the existing integrative index: `z(beta_propensity) + z(KD_max) + z(aromaticity) - z(abs_charge_density) - z(pI_distance_from_7)`, divided by 5. I also report AmyloGram probability, AMYPred probability, strict consensus Amyloid fraction, and expression-weighted panel index.",
        "",
        "## Main results",
    ]
    for _, row in primary.iterrows():
        lines.append(
            f"- {row['analysis_set']} / {row['metric']}: mean AC={row['mean_AC']:.4f}, "
            f"mean TOL={row['mean_TOL']:.4f}, delta TOL-AC={row['delta_TOL_minus_AC']:.4f}, "
            f"Cliff's delta={row['cliffs_delta_TOL_vs_AC']:.3f}, MW p={row['mannwhitney_two_sided_p']:.3g}."
        )
    lines.extend(["", "## Top protein-level shifts by integrative index"])
    top = protein_summary.sort_values("delta_TOL_minus_AC_Amyloid_Index", ascending=False).head(12)
    for _, row in top.iterrows():
        lines.append(
            f"- {row['image_gene']} | {row['protein_id_no_version']}: "
            f"delta index={row['delta_TOL_minus_AC_Amyloid_Index']:.3f}, "
            f"delta AmyloGram={row['delta_TOL_minus_AC_AmyloGram_Prob']:.3f}, "
            f"delta AMYPred={row['delta_TOL_minus_AC_AMYPred_Prob']:.3f}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "This is an exploratory, image-selected panel rather than independent validation. Direction should be interpreted mainly as whether the proteins highlighted by the feature plot become more amyloidogenic in TOL/chili after combining the amyloid-related physicochemical features and expression.",
        ]
    )
    (OUT / "IMAGE_PANEL_AMYLOIDOGENICITY_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(primary.to_string(index=False))
    print(f"Outputs: {OUT}")


if __name__ == "__main__":
    main()
