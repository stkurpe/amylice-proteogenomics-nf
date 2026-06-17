#!/usr/bin/env python3
"""
Build one Zarr store from per-sample bioinformatics outputs stored on S3.

Expected S3 layout:
  s3://codex-test-ngsdata-calculations/prepared-bioinfo-proteome/{SAMPLE}/
    results_expression/abundance.tsv
    results_proteins/nonsense_candidates.txt
    results_amyloid/amyloid_combined_predictions.csv
    results_hla/hla_consensus.tsv


Install:
  pip install boto3 pandas numpy zarr numcodecs tqdm

Example:
  python build_bioinfo_zarr.py \
    --samples SRR32060315 SRRXXXX \
    --bucket codex-test-ngsdata-calculations \
    --prefix prepared-bioinfo-proteome \
    --out project.zarr

Or:
  python build_bioinfo_zarr.py --samples-file samples.txt --out project.zarr
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path
from typing import Iterable, Optional

import boto3
import numpy as np
import pandas as pd
import zarr
from numcodecs import VLenUTF8
from tqdm import tqdm


FILES = {
    "expression": "results_expression/abundance.tsv",
    "nonsense": "results_proteins/nonsense_candidates.txt",
    "amyloid": "results_amyloid/amyloid_combined_predictions.csv",
    "hla": "results_hla/hla_consensus.tsv",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Aggregate per-sample pipeline outputs into one Zarr store")
    p.add_argument("--samples", nargs="*", help="Sample IDs, e.g. SRR32060315 SRR32060316")
    p.add_argument("--samples-file", help="Text file with one sample ID per line")
    p.add_argument("--bucket", default="codex-test-ngsdata-calculations")
    p.add_argument("--prefix", default="prepared-bioinfo-proteome")
    p.add_argument("--out", default="bioinfo_project.zarr", help="Output Zarr directory")
    p.add_argument("--aws-profile", default=None, help="Optional AWS profile name")
    p.add_argument("--region", default=None, help="Optional AWS region")
    p.add_argument("--tmpdir", default=None)
    p.add_argument("--missing-ok", action="store_true", help="Skip missing files instead of failing")
    return p.parse_args()


def get_samples(args: argparse.Namespace) -> list[str]:
    samples: list[str] = []
    if args.samples:
        samples.extend(args.samples)
    if args.samples_file:
        with open(args.samples_file) as f:
            samples.extend([x.strip() for x in f if x.strip() and not x.startswith("#")])
    samples = list(dict.fromkeys(samples))
    if not samples:
        raise SystemExit("No samples provided. Use --samples or --samples-file")
    return samples


def s3_client(profile: Optional[str], region: Optional[str]):
    session = boto3.Session(profile_name=profile, region_name=region) if profile else boto3.Session(region_name=region)
    return session.client("s3")


def download_s3_file(s3, bucket: str, key: str, local: Path, missing_ok: bool) -> bool:
    local.parent.mkdir(parents=True, exist_ok=True)
    try:
        s3.download_file(bucket, key, str(local))
        return True
    except Exception as e:
        if missing_ok:
            print(f"[WARN] missing or inaccessible: s3://{bucket}/{key} ({e})")
            return False
        raise


def read_table(path: Path, kind: str) -> pd.DataFrame:
    if kind == "amyloid":
        return pd.read_csv(path)
    # abundance.tsv, nonsense_candidates.txt, hla_consensus.tsv are tab-separated in your examples
    return pd.read_csv(path, sep="\t")


def normalize_expression(df: pd.DataFrame) -> pd.DataFrame:
    """Extract stable identifiers from kallisto-like abundance.tsv."""
    df = df.copy()
    parts = df["target_id"].astype(str).str.split("|", expand=True)
    df.insert(0, "transcript_id", parts[0])
    if parts.shape[1] > 1:
        df.insert(1, "gene_id", parts[1])
    if parts.shape[1] > 5:
        df.insert(2, "gene_name", parts[5])
    return df


def normalize_sample_columns(df: pd.DataFrame, sample: str) -> pd.DataFrame:
    df = df.copy()
    if "sample" not in df.columns:
        df.insert(0, "sample", sample)
    return df


def write_df_zarr(group: zarr.Group, name: str, df: pd.DataFrame, overwrite: bool = True) -> None:
    """Store a pandas DataFrame as a columnar Zarr table.

    Uses Zarr v2 + VLenUTF8 for string columns, so the output remains readable
    without Arrow/Parquet. Numeric columns are stored as numeric arrays.
    """
    if overwrite and name in group:
        del group[name]
    g = group.create_group(name)
    g.attrs["type"] = "dataframe"
    g.attrs["n_rows"] = int(len(df))
    g.attrs["columns"] = list(map(str, df.columns))

    for col in df.columns:
        s = df[col]
        arr_name = str(col)
        if pd.api.types.is_numeric_dtype(s):
            values = s.to_numpy()
            g.create_dataset(arr_name, data=values, shape=values.shape, chunks=(min(len(values), 100_000),))
        else:
            values = s.astype("string").fillna("").astype(str).to_numpy(dtype=object)
            g.create_dataset(arr_name, data=values, shape=values.shape, chunks=(min(len(values), 100_000),), dtype=object, object_codec=VLenUTF8())


def make_matrix(
    per_sample: dict[str, pd.DataFrame],
    id_col: str,
    value_cols: list[str],
) -> tuple[list[str], list[str], dict[str, np.ndarray]]:
    """Create sample x feature matrices for selected numeric columns."""
    samples = list(per_sample.keys())
    features = sorted(set().union(*[set(df[id_col].astype(str)) for df in per_sample.values() if id_col in df.columns]))
    feature_index = {f: i for i, f in enumerate(features)}
    matrices: dict[str, np.ndarray] = {}

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


def write_string_array(group: zarr.Group, name: str, values: Iterable[str]) -> None:
    arr = np.asarray(list(values), dtype=object)
    if name in group:
        del group[name]
    group.create_dataset(name, data=arr, shape=arr.shape, chunks=(min(len(arr), 100_000),), dtype=object, object_codec=VLenUTF8())


def main() -> None:
    args = parse_args()
    samples = get_samples(args)
    s3 = s3_client(args.aws_profile, args.region)

    out = Path(args.out)
    if out.exists():
        # Safer explicit cleanup than silently appending incompatible data
        import shutil
        shutil.rmtree(out)

    root = zarr.open_group(str(out), mode="w")
    root.attrs["schema"] = "bioinfo-proteome-zarr-v1"
    root.attrs["bucket"] = args.bucket
    root.attrs["prefix"] = args.prefix
    root.attrs["samples"] = samples
    root.attrs["files"] = FILES

    sample_tables = root.create_group("samples")
    tables = root.create_group("tables")
    layers = root.create_group("layers")

    expression_by_sample: dict[str, pd.DataFrame] = {}
    amyloid_by_sample: dict[str, pd.DataFrame] = {}
    all_nonsense = []
    all_hla = []
    manifest = []

    with tempfile.TemporaryDirectory(dir=args.tmpdir) as tmp:
        tmp = Path(tmp)
        for sample in tqdm(samples, desc="samples"):
            sg = sample_tables.create_group(sample)
            for kind, rel in FILES.items():
                key = f"{args.prefix.rstrip('/')}/{sample}/{rel}"
                local = tmp / sample / rel
                ok = download_s3_file(s3, args.bucket, key, local, args.missing_ok)
                manifest.append({"sample": sample, "kind": kind, "s3": f"s3://{args.bucket}/{key}", "downloaded": ok})
                if not ok:
                    continue

                df = read_table(local, kind)
                if kind == "expression":
                    df = normalize_expression(df)
                    df = normalize_sample_columns(df, sample)
                    expression_by_sample[sample] = df
                elif kind == "amyloid":
                    df = normalize_sample_columns(df, sample)
                    amyloid_by_sample[sample] = df
                elif kind == "nonsense":
                    df = normalize_sample_columns(df, sample)
                    all_nonsense.append(df)
                elif kind == "hla":
                    # file already contains sample in your example, but enforce consistency if absent
                    df = normalize_sample_columns(df, sample)
                    all_hla.append(df)

                write_df_zarr(sg, kind, df)

    # Long aggregated tables
    if all_nonsense:
        write_df_zarr(tables, "nonsense_candidates", pd.concat(all_nonsense, ignore_index=True))
    if all_hla:
        write_df_zarr(tables, "hla_consensus", pd.concat(all_hla, ignore_index=True))

    # Matrix layers: samples x transcripts/features
    if expression_by_sample:
        expr_g = layers.create_group("expression")
        sample_axis, feature_axis, mats = make_matrix(
            expression_by_sample,
            id_col="transcript_id",
            value_cols=["tpm", "est_counts", "length", "eff_length"],
        )
        write_string_array(expr_g, "sample_ids", sample_axis)
        write_string_array(expr_g, "transcript_ids", feature_axis)
        for name, mat in mats.items():
            expr_g.create_dataset(name, data=mat, shape=mat.shape, chunks=(min(mat.shape[0], 128), min(mat.shape[1], 4096)))
        expr_g.attrs["axis_0"] = "sample_ids"
        expr_g.attrs["axis_1"] = "transcript_ids"

    if amyloid_by_sample:
        amy_g = layers.create_group("amyloid")
        prob_cols = ["AMYPred_Prob", "AMYPred_Pred", "AmyloGram_Prob", "AmyloGramPy_Prob", "AmyloGramPy_Pred"]
        prob_cols = [c for c in prob_cols if any(c in df.columns for df in amyloid_by_sample.values())]
        sample_axis, feature_axis, mats = make_matrix(amyloid_by_sample, id_col="Sequence_ID", value_cols=prob_cols)
        write_string_array(amy_g, "sample_ids", sample_axis)
        write_string_array(amy_g, "sequence_ids", feature_axis)
        for name, mat in mats.items():
            amy_g.create_dataset(name, data=mat, shape=mat.shape, chunks=(min(mat.shape[0], 128), min(mat.shape[1], 4096)))
        amy_g.attrs["axis_0"] = "sample_ids"
        amy_g.attrs["axis_1"] = "sequence_ids"

    write_df_zarr(root, "manifest", pd.DataFrame(manifest))
    print(f"Done: {out}")
    print("Groups:", list(root.group_keys()))


if __name__ == "__main__":
    main()
