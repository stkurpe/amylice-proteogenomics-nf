from pathlib import Path
import sys
import re

import numpy as np
import pandas as pd
import zarr
from numcodecs import VLenUTF8


ROOT = Path(__file__).resolve().parent
ZARR_PATH = ROOT / "project.zarr"
VALCALC_WORK = Path("/Users/user841/Projects/BIoinfServer/ValCalculation/work")
GENCODE_PATH = (
    ROOT
    / "analysis_outputs"
    / "amyloid_bc_control_vs_chili"
    / "rank_based_amyloid_burden_signature"
    / "gencode_v48_transcript_annotation.tsv"
)

LOCAL_AMYPRED_SAMPLES = [
    "SRR32060216",
    "SRR32060240",
    "SRR32060276",
    "SRR32060306",
    "SRR32060309",
    "SRR32060312",
    "SRR32060326",
    "SRR32060336",
    "SRR32060341",
]


def write_df_zarr(group: zarr.Group, name: str, df: pd.DataFrame, overwrite: bool = True) -> None:
    if overwrite and name in group:
        del group[name]
    g = group.create_group(name)
    g.attrs["type"] = "dataframe"
    g.attrs["n_rows"] = int(len(df))
    g.attrs["columns"] = list(map(str, df.columns))

    for col in df.columns:
        s = df[col]
        arr_name = str(col)
        chunks = (max(1, min(len(s), 100_000)),)
        if pd.api.types.is_numeric_dtype(s):
            values = s.to_numpy()
            g.create_dataset(arr_name, data=values, shape=values.shape, chunks=chunks)
        else:
            values = s.astype("string").fillna("").astype(str).to_numpy(dtype=object)
            g.create_dataset(
                arr_name,
                data=values,
                shape=values.shape,
                chunks=chunks,
                dtype=object,
                object_codec=VLenUTF8(),
            )


def zarr_df(group: zarr.Group) -> pd.DataFrame:
    return pd.DataFrame({col: group[col][:] for col in group.attrs["columns"]})


def write_string_array(group: zarr.Group, name: str, values: list[str]) -> None:
    arr = np.asarray(values, dtype=object)
    if name in group:
        del group[name]
    group.create_dataset(
        name,
        data=arr,
        shape=arr.shape,
        chunks=(max(1, min(len(arr), 100_000)),),
        dtype=object,
        object_codec=VLenUTF8(),
    )


def make_matrix(per_sample: dict[str, pd.DataFrame], id_col: str, value_cols: list[str]):
    samples = list(per_sample.keys())
    features = sorted(
        set().union(*[set(df[id_col].astype(str)) for df in per_sample.values() if id_col in df.columns])
    )
    feature_index = {f: i for i, f in enumerate(features)}
    matrices = {}

    for value_col in value_cols:
        mat = np.full((len(samples), len(features)), np.nan, dtype="float32")
        for i, sample in enumerate(samples):
            df = per_sample[sample]
            if id_col not in df.columns or value_col not in df.columns:
                continue
            sub = df[[id_col, value_col]].dropna(subset=[id_col])
            for fid, val in zip(sub[id_col].astype(str), pd.to_numeric(sub[value_col], errors="coerce")):
                mat[i, feature_index[fid]] = val
        matrices[value_col] = mat
    return samples, features, matrices


def normalize_sample_columns(df: pd.DataFrame, sample: str) -> pd.DataFrame:
    df = df.copy()
    if "sample" not in df.columns:
        df.insert(0, "sample", sample)
    else:
        df["sample"] = sample
    return df


def clean_pred(s: pd.Series) -> pd.Series:
    return s.astype("string").replace("", pd.NA)


def transcript_id_clean(value):
    match = re.search(r"(ENST\d+(?:\.\d+)?)", str(value))
    if not match:
        return pd.NA
    return match.group(1).split(".")[0]


def recompute_consensus(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "AMYPred_Pred" not in df.columns:
        return df

    amypred = clean_pred(df["AMYPred_Pred"])
    if "AmyloGramPy_Pred" in df.columns:
        amylogram = clean_pred(df["AmyloGramPy_Pred"])
    elif "AmyloGram_Pred" in df.columns:
        amylogram = clean_pred(df["AmyloGram_Pred"])
    else:
        return df

    if "AmyloGram_Pred" in df.columns:
        legacy = clean_pred(df["AmyloGram_Pred"])
        amylogram = amylogram.fillna(legacy)

    consensus = pd.Series("Partial", index=df.index, dtype="string")
    both_present = amypred.notna() & amylogram.notna()
    consensus.loc[both_present & amypred.eq(amylogram) & amypred.eq("Amyloid")] = "Amyloid"
    consensus.loc[both_present & amypred.eq(amylogram) & amypred.eq("Non-Amyloid")] = "Non-Amyloid"
    consensus.loc[both_present & ~amypred.eq(amylogram)] = "Discordant"
    df["Consensus"] = consensus.astype(object)
    return df


def classify_probability(s: pd.Series) -> pd.Series:
    prob = pd.to_numeric(s, errors="coerce")
    out = pd.Series("Intermediate", index=s.index, dtype="string")
    out.loc[prob > 0.9] = "Amyloidogenic"
    out.loc[prob < 0.2] = "Non-Amyloidogenic"
    out.loc[prob.isna()] = pd.NA
    return out


def add_stringent_probability_classes(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "AMYPred_Prob" in df.columns:
        df["AMYPred_Stringent_Class"] = classify_probability(df["AMYPred_Prob"]).astype(object)
    if "AmyloGramPy_Prob" in df.columns:
        df["AmyloGramPy_Stringent_Class"] = classify_probability(df["AmyloGramPy_Prob"]).astype(object)

    if {"AMYPred_Prob", "AmyloGramPy_Prob"}.issubset(df.columns):
        amypred = pd.to_numeric(df["AMYPred_Prob"], errors="coerce")
        amylogram = pd.to_numeric(df["AmyloGramPy_Prob"], errors="coerce")
        both_present = amypred.notna() & amylogram.notna()
        combined = pd.Series("Intermediate", index=df.index, dtype="string")
        combined.loc[amypred.isna() | amylogram.isna()] = pd.NA
        combined.loc[both_present & (amypred > 0.9) & (amylogram > 0.9)] = "Amyloidogenic"
        combined.loc[both_present & (amypred < 0.2) & (amylogram < 0.2)] = "Non-Amyloidogenic"
        df["BothPredictors_Stringent_Class"] = combined.astype(object)
    return df


def add_quantile_probability_classes(
    df: pd.DataFrame,
    thresholds: dict[str, tuple[float, float]] | None = None,
) -> pd.DataFrame:
    df = df.copy()
    if thresholds is None:
        thresholds = {}

    for prob_col, class_col in [
        ("AMYPred_Prob", "AMYPred_Q17_Q83_Class"),
        ("AmyloGramPy_Prob", "AmyloGramPy_Q17_Q83_Class"),
    ]:
        if prob_col not in df.columns or prob_col not in thresholds:
            continue
        low, high = thresholds[prob_col]
        prob = pd.to_numeric(df[prob_col], errors="coerce")
        out = pd.Series("Intermediate", index=df.index, dtype="string")
        out.loc[prob <= low] = "Non-Amyloidogenic"
        out.loc[prob >= high] = "Amyloidogenic"
        out.loc[prob.isna()] = pd.NA
        df[class_col] = out.astype(object)

    if {"AMYPred_Prob", "AmyloGramPy_Prob"}.issubset(df.columns) and {
        "AMYPred_Prob",
        "AmyloGramPy_Prob",
    }.issubset(thresholds):
        amypred = pd.to_numeric(df["AMYPred_Prob"], errors="coerce")
        amylogram = pd.to_numeric(df["AmyloGramPy_Prob"], errors="coerce")
        amypred_low, amypred_high = thresholds["AMYPred_Prob"]
        amylogram_low, amylogram_high = thresholds["AmyloGramPy_Prob"]
        both_present = amypred.notna() & amylogram.notna()
        combined = pd.Series("Intermediate", index=df.index, dtype="string")
        combined.loc[~both_present] = pd.NA
        combined.loc[both_present & (amypred >= amypred_high) & (amylogram >= amylogram_high)] = "Amyloidogenic"
        combined.loc[both_present & (amypred <= amypred_low) & (amylogram <= amylogram_low)] = "Non-Amyloidogenic"
        df["BothPredictors_Q17_Q83_Class"] = combined.astype(object)
    return df


def compute_quantile_thresholds(root: zarr.Group) -> dict[str, tuple[float, float]]:
    values = {"AMYPred_Prob": [], "AmyloGramPy_Prob": []}
    for sample in root.attrs["samples"]:
        df = zarr_df(root["samples"][sample]["amyloid"])
        for col in values:
            if col in df.columns:
                s = pd.to_numeric(df[col], errors="coerce").dropna()
                if len(s):
                    values[col].append(s)

    thresholds: dict[str, tuple[float, float]] = {}
    for col, chunks in values.items():
        if not chunks:
            continue
        all_values = pd.concat(chunks, ignore_index=True)
        thresholds[col] = (float(all_values.quantile(0.17)), float(all_values.quantile(0.83)))
    return thresholds


def load_gencode_annotation() -> pd.DataFrame:
    if not GENCODE_PATH.exists():
        raise FileNotFoundError(f"Missing GENCODE annotation: {GENCODE_PATH}")
    ann = pd.read_csv(GENCODE_PATH, sep="\t")
    return ann[["Transcript_ID_clean", "gene_id", "gene_name", "biotype"]].drop_duplicates("Transcript_ID_clean")


def add_protein_coding_consensus_class(df: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["Transcript_ID_clean"] = df["Sequence_ID"].map(transcript_id_clean)
    df = df.drop(columns=["gene_id", "gene_name", "biotype"], errors="ignore")
    df = df.merge(ann, on="Transcript_ID_clean", how="left")
    df["gencode_biotype"] = df["biotype"].fillna("unannotated")
    df["gencode_gene_id"] = df["gene_id"].astype("string").fillna("").astype(object)
    df["gencode_gene_name"] = df["gene_name"].astype("string").fillna("").astype(object)
    df["is_gencode_protein_coding"] = df["gencode_biotype"].eq("protein_coding")

    if "BothPredictors_Q17_Q83_Class" in df.columns:
        consensus = clean_pred(df["BothPredictors_Q17_Q83_Class"])
    else:
        consensus = pd.Series(pd.NA, index=df.index, dtype="string")
    protein_coding_consensus = consensus.where(df["is_gencode_protein_coding"], pd.NA)
    df["ProteinCoding_BothPredictors_Q17_Q83_Class"] = protein_coding_consensus.astype(object)
    df = df.drop(columns=["gene_id", "gene_name", "biotype"], errors="ignore")
    return df


def main() -> int:
    if not ZARR_PATH.exists():
        print(f"Missing zarr: {ZARR_PATH}", file=sys.stderr)
        return 1

    root = zarr.open_group(str(ZARR_PATH), mode="a")
    samples = list(root.attrs["samples"])
    sample_tables = root["samples"]
    gencode_ann = load_gencode_annotation()

    updated = []
    for sample in LOCAL_AMYPRED_SAMPLES:
        local_csv = VALCALC_WORK / sample / "results" / "amyloid_combined_predictions.csv"
        if sample not in samples:
            print(f"SKIP {sample}: not present in project.zarr")
            continue
        if not local_csv.exists():
            print(f"SKIP {sample}: missing {local_csv}")
            continue

        df = pd.read_csv(local_csv)
        df = normalize_sample_columns(df, sample)
        df = recompute_consensus(df)
        df = add_stringent_probability_classes(df)
        write_df_zarr(sample_tables[sample], "amyloid", df)
        updated.append(sample)
        present = int(df["AMYPred_Pred"].notna().sum()) if "AMYPred_Pred" in df.columns else 0
        print(f"UPDATED {sample}: rows={len(df)} AMYPred_Pred_present={present}")

    thresholds = compute_quantile_thresholds(root)
    root.attrs["q17_q83_probability_thresholds"] = {
        col: {"q17": low, "q83": high} for col, (low, high) in thresholds.items()
    }
    root.attrs["protein_coding_consensus_class"] = {
        "column": "ProteinCoding_BothPredictors_Q17_Q83_Class",
        "source_class": "BothPredictors_Q17_Q83_Class",
        "filter": "GENCODE biotype == protein_coding",
        "gencode_annotation": str(GENCODE_PATH.relative_to(ROOT)),
    }

    amyloid_by_sample = {}
    for sample in samples:
        df = zarr_df(sample_tables[sample]["amyloid"])
        df = recompute_consensus(df)
        df = add_stringent_probability_classes(df)
        df = add_quantile_probability_classes(df, thresholds)
        df = add_protein_coding_consensus_class(df, gencode_ann)
        write_df_zarr(sample_tables[sample], "amyloid", df)
        amyloid_by_sample[sample] = df

    layers = root["layers"]
    if "amyloid" in layers:
        del layers["amyloid"]
    amy_g = layers.create_group("amyloid")
    prob_cols = ["AMYPred_Prob", "AMYPred_Pred", "AmyloGram_Prob", "AmyloGramPy_Prob", "AmyloGramPy_Pred"]
    prob_cols = [c for c in prob_cols if any(c in df.columns for df in amyloid_by_sample.values())]
    sample_axis, feature_axis, mats = make_matrix(amyloid_by_sample, id_col="Sequence_ID", value_cols=prob_cols)
    write_string_array(amy_g, "sample_ids", sample_axis)
    write_string_array(amy_g, "sequence_ids", feature_axis)
    for name, mat in mats.items():
        amy_g.create_dataset(
            name,
            data=mat,
            shape=mat.shape,
            chunks=(min(mat.shape[0], 128), min(mat.shape[1], 4096)),
        )
    amy_g.attrs["axis_0"] = "sample_ids"
    amy_g.attrs["axis_1"] = "sequence_ids"

    print(f"Rebuilt layers/amyloid for {len(samples)} samples; updated {len(updated)} sample tables.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
