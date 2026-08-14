#!/usr/bin/env python3
"""ingest_mw_v2.py — Phase 8 re-parse of McCance & Widdowson 2021.

Fixes the v1 parser's two defects:

  1. v1 derived the foods table from `unique_foods = melted_df[melted_df['amount'] > 0]
     .drop_duplicates()` — any food whose ALL measured nutrients were "Tr" (trace),
     "N" (not present), or NaN got silently dropped. That explained the
     1322 / 2888 = 46 % retention.

  2. v1 only mapped ~45 of the ~600 nutrient columns across 8 of the 14 sheets,
     skipping vitamin fractions, fatty-acid fractions per 100 g FA, factors,
     and most organic acids.

v2 approach:

  - Foods spine = '1.3 Proximates' (every McCance food appears there, even if
    its nutrient values are all 'Tr').
  - Per-sheet header detection via openpyxl scan for 'Food Code' literal, not
    iloc[2:].
  - Walk EVERY numeric column in EVERY data sheet (1.2 - 1.14). Map known
    columns to FDC nutrient_ids; mint new IDs in 200001+ block for the rest.
    Emit extra_nutrient_map.tsv recording each new ID with its source column
    and unit.
  - Keep 'Tr' -> 0 but DO NOT filter out zero amounts (zero is a real
    measurement for the EC graph).
  - Build the foods table directly from the spine, independent of the
    nutrient long-format table.

Outputs (overwrites v1's files):
  mccance_food_injection.parquet            — 2888 foods (up from 1322)
  mccance_food_nutrient_bucketed/bucket=*/  — Hive parquet partitions
  mccance_extra_nutrient_map.tsv            — newly-minted nutrient IDs
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

THIS = Path(__file__).resolve().parent
ROOT = THIS.parent
sys.path.insert(0, str(ROOT))
from _common import header_detect  # noqa: E402

# fdc_id allocation is centralised in food_DBs/fdc_blocks.py. Never write a literal offset:
# ids are accessions looked up in fdc_id_map.tsv, so a food keeps its id across releases.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
import fdc_blocks



EXCEL_DEFAULT = "McCance_Widdowsons_Composition_of_Foods_Integrated_Dataset_2021..xlsx"
# RETIRED: block bases live in food_DBs/fdc_blocks.py. This constant is kept out of
# the file entirely so it cannot drift from the real allocation.
EXTRA_ID_START = 200_001        # nutrient IDs ≥ 200_001 are reserved for novel nutrients

# Data sheets (drop "List of tables" + "1.1 Notes" — they are metadata)
DATA_SHEETS = [
    "1.2 Factors",
    "1.3 Proximates",
    "1.4 Inorganics",
    "1.5 Vitamins",
    "1.6 Vitamin Fractions",
    "1.7 (SFA per 100gFA)",
    "1.8 (SFA per 100gFood)",
    "1.9 (MUFA per 100FA)",
    "1.10 (MUFA per 100gFood)",
    "1.11 (PUFA per 100gFA)",
    "1.12 (PUFA per 100gFood)",
    "1.13 Phytosterols",
    "1.14 Organic Acids",
]

SPINE_SHEET = "1.3 Proximates"

# Per-food identity columns that exist on most sheets
META_COLS = {
    "food code", "food name", "description", "group",
    "previous", "main data references", "footnote",
}

# Known McCance column name -> FDC nutrient_id. (Comprehensive but doesn't have
# to cover everything; unmatched columns get a minted ID.)
KNOWN_NUTRIENTS: dict[str, int] = {
    # Proximates
    "water (g)": 1051,
    "total nitrogen (g)": 1002,
    "protein (g)": 1003,
    "fat (g)": 1004,
    "carbohydrate (g)": 1005,
    "energy (kcal)": 1008,
    "energy (kj)": 1062,
    "starch (g)": 1009,
    "oligosaccharide (g)": 1071,
    "total sugars (g)": 2000,
    "glucose (g)": 1011,
    "galactose (g)": 1015,
    "fructose (g)": 1012,
    "sucrose (g)": 1010,
    "maltose (g)": 1014,
    "lactose (g)": 1013,
    "alcohol (g)": 1018,
    "aoac fibre (g)": 1079,
    "englyst fibre (g)": 1085,
    "dietary fibre (g)": 1079,
    # Inorganics
    "sodium (mg)": 1093,
    "potassium (mg)": 1092,
    "calcium (mg)": 1087,
    "magnesium (mg)": 1090,
    "phosphorus (mg)": 1091,
    "iron (mg)": 1089,
    "copper (mg)": 1098,
    "zinc (mg)": 1095,
    "chloride (mg)": 1088,
    "manganese (mg)": 1101,
    "selenium (µg)": 1103,
    "iodine (µg)": 1100,
    # Vitamins
    "retinol (µg)": 1105,
    "carotene (µg)": 1106,
    "vitamin a (µg)": 1104,
    "vitamin d (µg)": 1110,
    "vitamin e (mg)": 1109,
    "vitamin k1 (µg)": 1185,
    "thiamin (mg)": 1165,
    "riboflavin (mg)": 1166,
    "niacin (mg)": 1167,
    "tryptophan/60 (mg)": 1210,
    "niacin equivalent (mg)": 1167,
    "vitamin b6 (mg)": 1175,
    "vitamin b12 (µg)": 1174,
    "folate (µg)": 1177,
    "pantothenate (mg)": 1170,
    "biotin (µg)": 1176,
    "vitamin c (mg)": 1162,
    # Lipids
    "cholesterol (mg)": 1253,
    "saturated fatty acids (g)": 1258,
    "monounsaturated fatty acids (g)": 1259,
    "polyunsaturated fatty acids (g)": 1260,
    # Phytosterols
    "total phytosterols (mg)": 1283,
    "beta-sitosterol (mg)": 1288,
    "campesterol (mg)": 1287,
    "stigmasterol (mg)": 1286,
    # Organic acids
    "citric acid (g)": 1021,
    "malic acid (g)": 1022,
    "lactic acid (g)": 1038,
    "acetic acid (g)": 1023,
    "oxalic acid (g)": 1029,
    "fumaric acid (g)": 1037,
    "tartaric acid (g)": 1030,
    "succinic acid (g)": 1031,
    "propionic acid (g)": 1024,
    "formic acid (g)": 1036,
}

CAT_MAP = {
    "DG": "Vegetables and Vegetable Products",
    "DR": "Fruits and Fruit Juices",
    "C":  "Cereal Grains and Pasta",
    "A":  "Cereal Grains and Pasta",
    "B":  "Cereal Grains and Pasta",
    "N":  "Nut and Seed Products",
    "M":  "Dairy and Egg Products",
    "O":  "Fats and Oils",
    "H":  "Beef Products",
    "J":  "Finfish and Shellfish Products",
    "K":  "Pork Products",
    "L":  "Poultry Products",
    "S":  "Sweets",
    "P":  "Beverages",
    "Q":  "Alcoholic Beverages",
    "DI": "Spices and Herbs",
}


_TRACE_TOKENS = {"tr", "tr.", "trace", "traces"}
_MISSING_TOKENS = {"n", "n.d", "n.d.", "nd", "-", "", "na", "nan", "n/a"}


def coerce_amount(x) -> float:
    """McCance numeric coercion: 'Tr' -> 0, 'N' / '-' / blank -> NaN."""
    if x is None:
        return float("nan")
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x) if np.isfinite(x) else float("nan")
    s = str(x).strip()
    lo = s.lower()
    if lo in _TRACE_TOKENS:
        return 0.0
    if lo in _MISSING_TOKENS:
        return float("nan")
    # "<0.1" -> 0
    if s.startswith("<"):
        return 0.0
    s = s.replace(",", ".")
    try:
        v = float(s)
        return v if np.isfinite(v) else float("nan")
    except ValueError:
        return float("nan")


def normalize_col(c: str) -> str:
    return re.sub(r"\s+", " ", str(c).strip()).lower()


def read_sheet(excel_path: Path, sheet: str) -> Optional[pd.DataFrame]:
    """Read one McCance sheet with programmatic header detection.

    Some sheets (e.g. 1.4 Inorganics) have a blank-space as the first column
    header — the food-code column. We detect the header row by looking for
    *any* of 'Food Code' / 'Food Name' / 'Group', then forcibly rename the
    first column to 'Food Code' regardless of what it actually says.
    """
    hdr = header_detect.find_header_row(str(excel_path), sheet=sheet,
                                        needles=("Food Code", "Food Name", "Group"))
    if hdr is None:
        return None
    df = pd.read_excel(excel_path, sheet_name=sheet, header=hdr, engine="openpyxl")
    # First column is always the food code, even if its header cell is blank/whitespace
    first = df.columns[0]
    if str(first).strip().lower() not in {"food code"}:
        df = df.rename(columns={first: "Food Code"})
    if "Food Code" not in df.columns:
        return None
    df = df[df["Food Code"].notna()].copy()
    df["Food Code"] = df["Food Code"].astype(str).str.strip()
    df = df[df["Food Code"] != ""]
    df = df[~df["Food Code"].str.lower().isin({"nan", "none"})]
    # Drop spacer rows where Food Code is a non-pattern string (e.g. the
    # "Sodium" label that appears in row 2 of Inorganics)
    df = df[df["Food Code"].str.match(r"^\d+-\d+", na=False)].copy()
    return df



def _mccance_codes(spine):
    """The source's own Food Code, used as the accession key.

    The 2021 integrated dataset repeats '13-669' on two different foods. Passing both through
    unchanged would map two foods onto one accession -- exactly the silent merge this scheme
    exists to prevent -- so repeat occurrences are suffixed. The first keeps the clean code.
    """
    seen, out = {}, []
    for c in spine["Food Code"].astype(str).str.strip():
        seen[c] = seen.get(c, 0) + 1
        out.append(c if seen[c] == 1 else f"{c}#{seen[c]}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--excel", default=EXCEL_DEFAULT)
    ap.add_argument("--out_prefix", default="mccance")
    args = ap.parse_args()

    os.chdir(THIS)
    excel_path = Path(args.excel)
    if not excel_path.exists():
        sys.exit(f"Excel not found: {excel_path}")

    print("==================================================")
    print("McCance v2 ingester")
    print("==================================================")

    # --- 1) Build foods spine from 1.3 Proximates -------------------------
    print(f"[1/5] Building foods spine from {SPINE_SHEET!r} ...")
    spine = read_sheet(excel_path, SPINE_SHEET)
    if spine is None:
        sys.exit("Could not read spine sheet")
    spine = spine[["Food Code", "Food Name", "Description", "Group"]].copy()
    # Keyed on the source's Food Code. range() left McCance with no key at all, so any
    # republished edition would have silently renumbered every food.
    spine["fdc_id"] = fdc_blocks.assign("mccance", _mccance_codes(spine))
    n_foods = len(spine)
    print(f"  -> {n_foods} foods on spine")

    # Map Group -> category (use first-letter trick same as v1 + 2-char fallbacks)
    def map_cat(g) -> str:
        g = str(g).strip().upper()
        for k in (g[:2], g[:1]):
            if k in CAT_MAP:
                return CAT_MAP[k]
        return "Meals, Entrees, and Side Dishes"
    spine["food_category"] = spine["Group"].apply(map_cat)

    # --- 2) Walk every data sheet, collect (Food Code, raw_col, amount) ---
    print("[2/5] Walking nutrient sheets...")
    long_rows: list[tuple[str, str, float]] = []   # (Food Code, raw_col_norm, amount)
    raw_col_to_display: dict[str, str] = {}
    raw_col_sheet: dict[str, str] = {}             # which sheet first contributed each col
    sheet_summary: list[tuple[str, int, int]] = []

    for sheet in DATA_SHEETS:
        df = read_sheet(excel_path, sheet)
        if df is None:
            print(f"  [skip] {sheet} — header not detected")
            continue
        # Identify nutrient columns: anything not in META_COLS that is mostly-numeric
        nutrient_cols = []
        for c in df.columns:
            if normalize_col(c) in META_COLS:
                continue
            nutrient_cols.append(c)
        # Coerce
        for c in nutrient_cols:
            df[c] = df[c].apply(coerce_amount)
        # Melt
        for c in nutrient_cols:
            raw_norm = normalize_col(c)
            raw_col_to_display.setdefault(raw_norm, str(c).strip())
            raw_col_sheet.setdefault(raw_norm, sheet)
        sub = df[["Food Code"] + nutrient_cols].copy()
        sub_long = sub.melt(id_vars="Food Code", var_name="raw_col", value_name="amount")
        sub_long["raw_col"] = sub_long["raw_col"].apply(normalize_col)
        sub_long = sub_long.dropna(subset=["amount"])
        for _, r in sub_long.iterrows():
            long_rows.append((r["Food Code"], r["raw_col"], r["amount"]))
        sheet_summary.append((sheet, len(df), len(nutrient_cols)))
    print(f"  -> read {len(sheet_summary)} sheets, collected {len(long_rows)} measurements")
    for sn, nr, nc in sheet_summary:
        print(f"     {sn:35s}  {nr:5d} rows × {nc:4d} cols")

    # --- 3) Map raw_col -> nutrient_id (known FDC ID or minted) ----------
    print("[3/5] Mapping nutrient columns to IDs...")
    raw_to_id: dict[str, int] = {}
    extra_rows: list[dict] = []      # for extra_nutrient_map.tsv
    next_id = EXTRA_ID_START
    for raw_norm in sorted(raw_col_to_display):
        if raw_norm in KNOWN_NUTRIENTS:
            raw_to_id[raw_norm] = KNOWN_NUTRIENTS[raw_norm]
        else:
            raw_to_id[raw_norm] = next_id
            extra_rows.append({
                "nutrient_id": next_id,
                "source_db": "mccance",
                "source_column_raw": raw_col_to_display[raw_norm],
                "source_column_norm": raw_norm,
                "source_sheet": raw_col_sheet[raw_norm],
                "note": "minted from McCance v2 ingester (unmapped column)",
            })
            next_id += 1
    n_known = sum(1 for k in raw_to_id if k in KNOWN_NUTRIENTS)
    print(f"  -> {n_known} known FDC IDs, {len(extra_rows)} new IDs minted "
          f"(range {EXTRA_ID_START}-{next_id - 1})")

    # --- 4) Build the long-form (fdc_id, nutrient_id, amount) table -------
    print("[4/5] Joining onto fdc_id and bucketing...")
    long_df = pd.DataFrame(long_rows, columns=["Food Code", "raw_col", "amount"])
    long_df = long_df.merge(spine[["Food Code", "fdc_id"]], on="Food Code", how="left")
    # Drop rows that didn't match the spine (very few — typically subtitle rows
    # that survived in non-spine sheets)
    n_pre = len(long_df)
    long_df = long_df.dropna(subset=["fdc_id"])
    if len(long_df) < n_pre:
        print(f"  [warn] dropped {n_pre - len(long_df)} measurements with no spine match")
    long_df["fdc_id"] = long_df["fdc_id"].astype(int)
    long_df["nutrient_id"] = long_df["raw_col"].map(raw_to_id).astype(int)
    long_df = long_df[["fdc_id", "nutrient_id", "amount"]]

    # --- 5) Write outputs --------------------------------------------------
    print("[5/5] Writing parquet outputs...")
    foods = pd.DataFrame({
        "fdc_id": spine["fdc_id"].astype(int),
        "data_type": "foundation_food",
        "description": spine["Food Name"].astype(str).str.strip() + " [McCance]",
        "food_category_id": "5555",
        "food_category": spine["food_category"].astype(str),
    })
    out_food = f"{args.out_prefix}_food_injection.parquet"
    foods.to_parquet(out_food, index=False)

    out_bucket = Path(f"{args.out_prefix}_food_nutrient_bucketed")
    # Wipe any old contents to avoid stale-data mixing
    if out_bucket.exists():
        import shutil
        shutil.rmtree(out_bucket)
    out_bucket.mkdir(exist_ok=True)
    long_df["bucket"] = (long_df["nutrient_id"] % 256).astype(int)
    for b, g in long_df.groupby("bucket", sort=True):
        d = out_bucket / f"bucket={int(b)}"
        d.mkdir(exist_ok=True)
        g[["fdc_id", "nutrient_id", "amount"]].sort_values(["fdc_id", "nutrient_id"]) \
            .to_parquet(d / "mccance_data.parquet", index=False)

    if extra_rows:
        pd.DataFrame(extra_rows).to_csv(f"{args.out_prefix}_extra_nutrient_map.tsv",
                                         sep="\t", index=False)

    print()
    print(f"[OK] foods                       : {len(foods)}")
    print(f"[OK] long-form measurements      : {len(long_df)}")
    print(f"[OK] distinct nutrient ids       : {long_df['nutrient_id'].nunique()}")
    print(f"[OK] of which known FDC          : {n_known}")
    print(f"[OK] of which minted (≥{EXTRA_ID_START}) : {len(extra_rows)}")
    print(f"[OK] outputs:")
    print(f"      {out_food}")
    print(f"      {out_bucket}/")
    print(f"      {args.out_prefix}_extra_nutrient_map.tsv")


if __name__ == "__main__":
    main()
