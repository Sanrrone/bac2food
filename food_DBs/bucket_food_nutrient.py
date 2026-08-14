#!/usr/bin/env python3
"""
bucket_food_nutrient.py

Create a bucket-partitioned Parquet dataset from a FoodData Central food_nutrient table.

Why
---
Partitioning by bucket = nutrient_id % B allows fast "read only the nutrients I need"
using Parquet partition pruning.

Input
-----
--food_nutrient_in : CSV or Parquet with at least columns: fdc_id, nutrient_id, amount
  (If percent_daily_value exists, it will be kept only if you request it.)

Output
------
--out_dir : directory like food_nutrient_bucketed/
  Written as Hive-style partitions:
      out_dir/bucket=0/part-*.parquet
      out_dir/bucket=1/part-*.parquet
      ...
with columns: fdc_id, nutrient_id, amount, (optional percent_daily_value)

Notes
-----
- Designed for large files; CSV is processed in chunks.
- For Parquet input, we stream row-groups by reading with pandas then re-bucketing.
  (Still fast for typical FDC parquet sizes; if your parquet is huge, consider exporting to CSV first.)
- Requires: pandas, pyarrow
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Iterable, List, Optional

import pandas as pd


REQUIRED_COLS = ["fdc_id", "nutrient_id", "amount"]


def ensure_cols(df: pd.DataFrame, want_pdv: bool) -> pd.DataFrame:
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise SystemExit(f"Input missing required columns: {missing}")

    keep = REQUIRED_COLS.copy()
    if want_pdv and "percent_daily_value" in df.columns:
        keep.append("percent_daily_value")

    df = df[keep].copy()
    df["fdc_id"] = pd.to_numeric(df["fdc_id"], errors="coerce").astype("Int64")
    df["nutrient_id"] = pd.to_numeric(df["nutrient_id"], errors="coerce").astype("Int64")
    df["amount"] = pd.to_numeric(df["amount"], errors="coerce")

    if "percent_daily_value" in df.columns:
        df["percent_daily_value"] = pd.to_numeric(df["percent_daily_value"], errors="coerce")

    df = df.dropna(subset=["fdc_id", "nutrient_id", "amount"])
    return df


def write_bucket(df: pd.DataFrame, out_dir: Path, bucket: int, compression: str) -> None:
    bdir = out_dir / f"bucket={bucket}"
    bdir.mkdir(parents=True, exist_ok=True)
    # append-style: create unique file per write to avoid file locking issues
    fname = bdir / f"part-{os.getpid()}-{pd.Timestamp.utcnow().value}.parquet"
    df.to_parquet(fname, index=False, engine="pyarrow", compression=compression)


def bucketize_df(df: pd.DataFrame, buckets: int) -> pd.Series:
    # nutrient_id is Int64 nullable; convert safely
    nid = df["nutrient_id"].astype("int64", errors="ignore")
    # If nullable dtype survived, cast via to_numpy
    try:
        arr = nid.to_numpy(dtype="int64", na_value=-1)
    except Exception:
        arr = pd.to_numeric(df["nutrient_id"], errors="coerce").fillna(-1).astype("int64").to_numpy()
    return pd.Series(arr % int(buckets), index=df.index, dtype="int64")


def iter_csv_chunks(path: str, chunksize: int, want_pdv: bool) -> Iterable[pd.DataFrame]:
    usecols = REQUIRED_COLS + (["percent_daily_value"] if want_pdv else [])
    # usecols only if columns exist; for safety, read header first
    header = pd.read_csv(path, nrows=0)
    cols = [c for c in usecols if c in header.columns]

    reader = pd.read_csv(
        path,
        usecols=cols,
        chunksize=chunksize,
        low_memory=False,
        engine="c",
    )
    for chunk in reader:
        yield chunk


def process_csv(in_path: str, out_dir: Path, buckets: int, chunksize: int, want_pdv: bool, compression: str) -> None:
    for chunk in iter_csv_chunks(in_path, chunksize=chunksize, want_pdv=want_pdv):
        chunk = ensure_cols(chunk, want_pdv=want_pdv)
        if chunk.empty:
            continue
        chunk["bucket"] = bucketize_df(chunk, buckets=buckets)

        for b, sub in chunk.groupby("bucket", sort=False):
            sub = sub.drop(columns=["bucket"])
            write_bucket(sub, out_dir=out_dir, bucket=int(b), compression=compression)


def process_parquet(in_path: str, out_dir: Path, buckets: int, want_pdv: bool, compression: str) -> None:
    # Read in reasonably sized batches by row-groups via pyarrow if available; fallback to pandas read.
    # For most FDC parquet sizes, a full read is okay; if yours is massive, convert to CSV first.
    df = pd.read_parquet(in_path, engine="pyarrow")
    df = ensure_cols(df, want_pdv=want_pdv)
    if df.empty:
        return
    df["bucket"] = bucketize_df(df, buckets=buckets)
    for b, sub in df.groupby("bucket", sort=False):
        sub = sub.drop(columns=["bucket"])
        write_bucket(sub, out_dir=out_dir, bucket=int(b), compression=compression)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--food_nutrient_in", required=True, help="CSV or Parquet input file")
    ap.add_argument("--out_dir", required=True, help="Output directory for bucket partitions")
    ap.add_argument("--buckets", type=int, default=256, help="Number of buckets (default 256)")
    ap.add_argument("--chunksize", type=int, default=2_000_000, help="CSV chunksize rows (default 2,000,000)")
    ap.add_argument("--include_percent_daily_value", action="store_true", help="Keep percent_daily_value if present")
    ap.add_argument("--compression", default="snappy", help="Parquet compression (snappy, zstd, gzip, etc.)")
    ap.add_argument("--overwrite", action="store_true", help="Delete out_dir first if exists")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    if out_dir.exists() and args.overwrite:
        # careful recursive delete
        for p in sorted(out_dir.rglob("*"), reverse=True):
            if p.is_file():
                p.unlink()
            else:
                p.rmdir()
        out_dir.rmdir()

    out_dir.mkdir(parents=True, exist_ok=True)

    in_path = args.food_nutrient_in
    if in_path.lower().endswith((".parquet", ".pq")):
        process_parquet(
            in_path=in_path,
            out_dir=out_dir,
            buckets=int(args.buckets),
            want_pdv=bool(args.include_percent_daily_value),
            compression=str(args.compression),
        )
    else:
        process_csv(
            in_path=in_path,
            out_dir=out_dir,
            buckets=int(args.buckets),
            chunksize=int(args.chunksize),
            want_pdv=bool(args.include_percent_daily_value),
            compression=str(args.compression),
        )

    # write a small manifest
    manifest = out_dir / "_BUCKET_MANIFEST.txt"
    manifest.write_text(
        f"bucketed_from={in_path}\n"
        f"buckets={int(args.buckets)}\n"
        f"columns={','.join(REQUIRED_COLS + (['percent_daily_value'] if args.include_percent_daily_value else []))}\n",
        encoding="utf-8",
    )
    print(f"Wrote bucketed dataset to: {out_dir}")
    print(f"Manifest: {manifest}")


if __name__ == "__main__":
    main()
