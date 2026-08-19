#!/usr/bin/env python3
"""
ingest_swedish_livsmedelsdb.py

Ingest the Swedish food composition database (Livsmedelsverket Excel export)
into your "FDC-like" parquet artifacts:

Outputs
1) <prefix>_food_injection.parquet
   Columns: fdc_id, data_type, description, food_category_id, food_category

2) <prefix>_food_nutrient_bucketed/
   Hive partitions: bucket=<0..255>/swedish_data.parquet
   Columns: fdc_id, nutrient_id, amount

Notes
- Values are assumed "per 100 g" (as in typical composition tables).
- We mint unique fdc_id by adding an offset (default 40,000,000) to Food Number.
- Nutrient mapping is a curated subset to USDA/FDC nutrient IDs (extend as needed).
- Handles messy numeric cells like commas, "<", "trace", etc. -> coerces to 0.

Example
  python ingest_swedish_livsmedelsdb.py \
    --excel LivsmedelsDB.xlsx \
    --out_prefix sweden \
    --offset 40000000

"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd

# fdc_id allocation: fdc_blocks.py
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
import fdc_blocks



WS = re.compile(r"\s+")
NONNUM = re.compile(r"[^0-9,\.\-\+eE]")


# --- FDC CATEGORY MAPPING (Filter group -> FDC category label) ---
# Adjust freely; anything unmapped goes to "Meals, Entrees, and Side Dishes".
FILTER_GROUP_TO_FDC_CAT: Dict[str, str] = {
    # Plant
    "Vegetables": "Vegetables and Vegetable Products",
    "Vegetables and vegetable products": "Vegetables and Vegetable Products",
    "Fruits": "Fruits and Fruit Juices",
    "Fruit and berry products": "Fruits and Fruit Juices",
    "Berries": "Fruits and Fruit Juices",
    "Cereals and cereal products": "Cereal Grains and Pasta",
    "Bread": "Cereal Grains and Pasta",
    "Pasta": "Cereal Grains and Pasta",
    "Rice": "Cereal Grains and Pasta",
    "Potatoes": "Vegetables and Vegetable Products",
    "Legumes": "Legumes and Legume Products",
    "Nuts and seeds": "Nut and Seed Products",
    "Mushrooms": "Vegetables and Vegetable Products",
    "Herbs and spices": "Spices and Herbs",
    "Spices": "Spices and Herbs",

    # Animal
    "Meat": "Beef Products",
    "Meat products": "Sausages and Luncheon Meats",
    "Poultry": "Poultry Products",
    "Fish": "Finfish and Shellfish Products",
    "Shellfish": "Finfish and Shellfish Products",
    "Eggs": "Dairy and Egg Products",
    "Milk": "Dairy and Egg Products",
    "Milk and milk products": "Dairy and Egg Products",
    "Cheese": "Dairy and Egg Products",

    # Fats / sweets / beverages / prepared
    "Fats and oils": "Fats and Oils",
    "Sugar and sweets": "Sweets",
    "Cakes and pastries": "Baked Products",
    "Baked products": "Baked Products",
    "Soft drinks": "Beverages",
    "Beverages": "Beverages",
    "Alcoholic beverages": "Alcoholic Beverages",
    "Fast food": "Fast Foods",
    "Ready meals": "Meals, Entrees, and Side Dishes",
    "Soups": "Soups, Sauces, and Gravies",
    "Sauces": "Soups, Sauces, and Gravies",
}


# --- NUTRIENT MAPPING (Swedish column name -> FDC nutrient_id, expected unit) ---
# Extend this dict if your Excel contains additional analytes you care about.
# The IDs here are the common USDA nutrient IDs used in your pipeline.
NUTRIENT_MAP: Dict[str, int] = {
    # Energy / macros
    "Energy (kcal)": 1008,
    "Protein (g)": 1003,
    "Fat, total (g)": 1004,
    "Carbohydrates, available (g)": 1005,
    "Fibre (g)": 1079,
    "Sugars, total (g)": 2000,          # may or may not exist in your export
    "Alcohol (g)": 1018,                # may or may not exist
    "Water (g)": 1051,

    # Minerals
    "Calcium (mg)": 1087,
    "Iron (mg)": 1089,
    "Magnesium (mg)": 1090,
    "Phosphorus (mg)": 1091,
    "Potassium (mg)": 1092,
    "Sodium (mg)": 1093,
    "Zinc (mg)": 1095,
    "Copper (mg)": 1098,
    "Selenium (µg)": 1103,

    # Vitamins
    "Vitamin A (RAE) (µg)": 1104,
    "Vitamin D (µg)": 1110,
    "Vitamin C (mg)": 1162,
    "Thiamin (mg)": 1165,
    "Riboflavin (mg)": 1166,
    "Niacin (mg)": 1167,
    "Vitamin B-6 (mg)": 1175,
    "Folate, total (µg)": 1177,
    "Vitamin B-12 (µg)": 1178,

    # Lipids
    "Cholesterol (mg)": 1253,
    "Fatty acids, saturated (g)": 1258,
    "Fatty acids, monounsaturated (g)": 1259,
    "Fatty acids, polyunsaturated (g)": 1260,
}


def _norm_cat(s: str) -> str:
    s = (s or "").strip()
    s = WS.sub(" ", s)
    return s


def _coerce_num(x) -> float:
    """
    Robust numeric coercion:
    - Accepts comma decimals
    - Removes units/footnotes
    - Handles "<0.1", "trace", etc. -> 0.0
    """
    if x is None:
        return 0.0
    if isinstance(x, (int, float, np.integer, np.floating)):
        if np.isfinite(x):
            return float(x)
        return 0.0
    s = str(x).strip()
    if not s or s.lower() in {"na", "nan", "n/a", "trace", "traces", "sp"}:
        return 0.0
    # If it starts with "<" treat as 0 (conservative)
    if s.startswith("<"):
        return 0.0
    # keep numeric-ish characters
    s = NONNUM.sub("", s)
    s = s.replace(",", ".")
    try:
        v = float(s)
        return v if np.isfinite(v) else 0.0
    except Exception:
        return 0.0


def _find_header_row(excel_path: str, sheet=0, needle: str = "Food Name", max_scan: int = 25) -> int:
    """
    Detect header row by scanning first column block for 'Food Name'.
    Returns 0-indexed row number to pass as header=<row>.
    """
    raw = pd.read_excel(excel_path, sheet_name=sheet, header=None, nrows=max_scan)
    for i in range(min(max_scan, len(raw))):
        row = raw.iloc[i].astype(str).tolist()
        if any(needle == str(cell).strip() for cell in row):
            return i
    # fallback: common for this export
    return 2


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--excel", required=True, help="Swedish LivsmedelsDB Excel file")
    ap.add_argument("--sheet", default="0", help="Sheet index or name (default 0)")
    ap.add_argument("--out_prefix", required=True, help="Prefix for outputs (e.g., sweden)")
    # 41M, NOT 40M. McCance mints range(40_000_000, ...) positionally, so a 40M offset here put
    # every Swedish food on top of an unrelated McCance one: fdc_id 40000001 was "Agar, dried"
    # from McCance and "Beef tallow" from this table. 1,788 Swedish food records were overwritten
    # and their nutrient values left attached to the McCance food. Blocks must not be shared.
    ap.add_argument("--food_category_id", default="7777", help="Custom category id string for injected foods")
    ap.add_argument("--drop_zero", action="store_true", help="Drop zero-valued nutrient rows (recommended)")
    ap.add_argument("--limit_rows", type=int, default=0, help="Debug: process only first N foods")
    args = ap.parse_args()

    sheet = int(args.sheet) if str(args.sheet).isdigit() else args.sheet
    header_row = _find_header_row(args.excel, sheet=sheet)

    df = pd.read_excel(args.excel, sheet_name=sheet, header=header_row)
    # Basic required columns
    need_cols = ["Food Name", "Food Number", "Filter group"]
    for c in need_cols:
        if c not in df.columns:
            raise ValueError(f"Missing required column '{c}'. Found columns: {list(df.columns)[:25]} ...")

    if args.limit_rows and args.limit_rows > 0:
        df = df.head(args.limit_rows).copy()

    # Mint fdc_id
    df["Food Number"] = pd.to_numeric(df["Food Number"], errors="coerce")
    df = df.dropna(subset=["Food Number", "Food Name"]).copy()
    df["Food Number"] = df["Food Number"].astype(int)
    df["fdc_id"] = fdc_blocks.assign("swedish", df["Food Number"].astype(int))

    # Category mapping
    df["filter_group_norm"] = df["Filter group"].astype(str).map(_norm_cat)
    df["food_category"] = df["filter_group_norm"].map(FILTER_GROUP_TO_FDC_CAT).fillna("Meals, Entrees, and Side Dishes")

    # Food injection parquet
    food_inj = pd.DataFrame({
        "fdc_id": df["fdc_id"].astype(int),
        "data_type": "foundation_food",
        "description": df["Food Name"].astype(str).str.strip() + " [SWE]",
        "food_category_id": str(args.food_category_id),
        "food_category": df["food_category"].astype(str),
    }).drop_duplicates(subset=["fdc_id"])

    # Nutrient long table
    available = [c for c in NUTRIENT_MAP.keys() if c in df.columns]
    if not available:
        raise ValueError(
            "None of the mapped nutrient columns were found in the Excel. "
            "Inspect your sheet columns and extend NUTRIENT_MAP."
        )

    # Coerce nutrient columns to numeric
    for c in available:
        df[c] = df[c].map(_coerce_num)

    long = df.melt(
        id_vars=["fdc_id"],
        value_vars=available,
        var_name="swedish_nutrient",
        value_name="amount",
    )
    long["nutrient_id"] = long["swedish_nutrient"].map(NUTRIENT_MAP).astype(int)

    if args.drop_zero:
        long = long[long["amount"] > 0].copy()
    else:
        long = long[long["amount"].notna()].copy()

    # Bucketed output
    out_prefix = Path(args.out_prefix)
    food_out = out_prefix.with_name(out_prefix.name + "_food_injection.parquet")
    bucket_dir = out_prefix.with_name(out_prefix.name + "_food_nutrient_bucketed")
    bucket_dir.mkdir(parents=True, exist_ok=True)

    food_inj.to_parquet(food_out, index=False)

    long["bucket"] = (long["nutrient_id"] % 256).astype(int)
    for b, g in long.groupby("bucket", sort=True):
        part_dir = bucket_dir / f"bucket={int(b)}"
        part_dir.mkdir(parents=True, exist_ok=True)
        g = g.sort_values(["fdc_id", "nutrient_id"])
        g[["fdc_id", "nutrient_id", "amount"]].to_parquet(part_dir / "swedish_data.parquet", index=False)

    print("[OK] Wrote:")
    print(f"  - {food_out}  (foods: {len(food_inj):,})")
    print(f"  - {bucket_dir}/bucket=*  (rows: {len(long):,}, nutrients: {long['nutrient_id'].nunique():,})")
    print("[INFO] Nutrient columns ingested:")
    for c in available:
        print(f"  - {c} -> {NUTRIENT_MAP[c]}")


if __name__ == "__main__":
    main()
