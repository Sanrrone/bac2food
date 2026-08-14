#!/usr/bin/env python3
"""injest_phyfood_v2.py — Phase 8 re-parse of FAO/INFOODS PhyFoodComp 1.0.

v1 hardcoded a 5-entry NUTRIENT_MAP (Water, Fe, Zn, Ca + a coalesced phytate)
and silently dropped the other ~35 nutrient columns including every individual
inositol-phosphate measurement (IP3..IP6, IPSUM), every phytate measurement
method (PHYTCA, PHYTCPP, PHYTCPPI, PHYTCPPD, PHYTC-), and the phytate-to-Fe
and phytate-to-Zn molar ratios (bioavailability indicators).

v1 also assigned WATER(g) -> FDC nutrient_id 1005, which is wrong (1005 is
'Carbohydrate, by difference'; water is 1051).

v2 approach:
  - Walk every sheet matching the '<NN> <name>' pattern (18 nutrient sheets).
  - Per sheet: header detection by 'Food item ID' literal (reused from v1).
  - Column-walk every numeric column outside the metadata set.
  - Mint new nutrient_ids in the 210001+ block for unknown columns (separate
    from McCance's 200001+ block to keep collisions impossible).
  - Output bucketed parquet (not CSV) to match the bac2food pipeline.

Outputs:
  phyfoodcomp_food_injection.parquet
  phyfoodcomp_food_nutrient_bucketed/bucket=*/phyfoodcomp_data.parquet
  phyfoodcomp_extra_nutrient_map.tsv
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

THIS = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS.parent))
from _common import header_detect  # noqa: E402

# fdc_id allocation is centralised in food_DBs/fdc_blocks.py. Never write a literal offset:
# ids are accessions looked up in fdc_id_map.tsv, so a food keeps its id across releases.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
import fdc_blocks



EXCEL_DEFAULT = "PhyFoodComp_1.0.xlsx"
# RETIRED: block bases live in food_DBs/fdc_blocks.py. This constant is kept out of
# the file entirely so it cannot drift from the real allocation.
EXTRA_ID_START = 210_001

CATEGORY_MAP = {
    1: "Cereal Grains and Pasta",
    2: "Vegetables and Vegetable Products",
    3: "Legumes and Legume Products",
    4: "Vegetables and Vegetable Products",
    5: "Fruits and Fruit Juices",
    6: "Nut and Seed Products",
    7: "Beef Products",
    8: "American Indian/Alaska Native Foods",
    10: "Finfish and Shellfish Products",
    11: "Dairy and Egg Products",
    12: "Fats and Oils",
    13: "Beverages",
    14: "Sweets",
    15: "Spices and Herbs",
    16: "Baby Foods",
    17: "Spices and Herbs",
    18: "Spices and Herbs",
    19: "Meals, Entrees, and Side Dishes",
}

# Columns to never treat as numeric (metadata / categorical)
META_COLS = {
    "food item id", "old code (as in the original source)", "food group",
    "subgroup", "foodex2 code", "foodex2 name", "missing facet",
    "exact match", "matching comments", "country, region", "type",
    "food name in own language", "food name in english",
    "processing / influencing factors", "species/subspecies",
    "cultivar/variety/accession name", "season", "other", "n",
    "comments on data processing/methods", "publication year",
    "biblioid", "compiler id", "latest revision in version",
    "analytical/biodiversity",
    "comments on why some data is not entered",
}

# Known anchors that exist in FDC's nutrient.csv. The rest get minted.
KNOWN_NUTRIENTS = {
    "water(g)": 1051,
    "fe(mg)": 1089,
    "zn(mg)": 1095,
    "ca(mg)": 1087,
}


def clean_amount(x) -> float:
    if pd.isna(x):
        return float("nan")
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x) if np.isfinite(x) else float("nan")
    s = str(x).strip()
    if not s or s in {"-", "ND", "N.D.", "n.d.", "nd"}:
        return float("nan")
    if s.lower() in {"tr", "trace"}:
        return 0.0
    if s.lower() in {"<lod", "< lod"}:
        return 0.0
    # strip [a,b] -> a
    s = re.sub(r"[\[\]]", "", s).split(",")[0]
    s = re.sub(r"[^\d\.\-]", "", s)
    try:
        return float(s)
    except ValueError:
        return float("nan")


def normalize_col(c: str) -> str:
    return " ".join(str(c).split()).strip().lower()


def read_data_sheet(excel_path: Path, sheet: str) -> pd.DataFrame:
    hdr = header_detect.find_header_row(str(excel_path), sheet=sheet,
                                        needles=("Food item ID", "Food Item ID"))
    if hdr is None:
        return pd.DataFrame()
    df = pd.read_excel(excel_path, sheet_name=sheet, header=hdr, engine="openpyxl")
    df.columns = [" ".join(str(c).split()) for c in df.columns]
    if "Food item ID" not in df.columns:
        return pd.DataFrame()
    df["Food item ID"] = pd.to_numeric(df["Food item ID"], errors="coerce")
    df = df.dropna(subset=["Food item ID"]).copy()
    df["Food item ID"] = df["Food item ID"].astype(int)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--excel", default=EXCEL_DEFAULT)
    ap.add_argument("--out_prefix", default="phyfoodcomp")
    args = ap.parse_args()

    os.chdir(THIS)
    excel_path = Path(args.excel)
    if not excel_path.exists():
        sys.exit(f"missing: {excel_path}")

    print("==================================================")
    print("PhyFoodComp v2 ingester")
    print("==================================================")

    xls = pd.ExcelFile(excel_path)
    data_sheets = [s for s in xls.sheet_names if re.match(r"^\d{2}\s", s)]
    print(f"[1/4] {len(data_sheets)} data sheets:")
    for s in data_sheets:
        print(f"     {s}")

    print("[2/4] Reading + concatenating sheets...")
    frames = []
    for s in data_sheets:
        df = read_data_sheet(excel_path, s)
        if not df.empty:
            frames.append(df)
            print(f"     {s:40s}  -> {len(df)} foods × {len(df.columns)} cols")
    master = pd.concat(frames, ignore_index=True)
    print(f"  -> total raw rows: {len(master)}")

    # Dedupe by Food item ID (some sheets repeat the same food across cuisines)
    master = master.drop_duplicates(subset=["Food item ID"], keep="first").reset_index(drop=True)
    print(f"  -> unique foods after dedup: {len(master)}")

    # --- foods table ---
    master["fdc_id"] = fdc_blocks.assign("phyfoodcomp", master["Food item ID"].astype(int))
    name_en = master.get("Food name in English", pd.Series([None] * len(master))) \
                    .astype(str).replace("nan", "")
    name_local = master.get("Food name in own language", pd.Series([None] * len(master))) \
                       .astype(str).replace("nan", "")
    species = master.get("Species/Subspecies", pd.Series([""] * len(master))) \
                    .astype(str).replace("nan", "")
    desc = name_en.where(name_en.str.strip() != "", name_local)
    desc = desc.str.strip()
    desc = desc.where(species.str.strip() == "", desc + " [" + species.str.strip() + "]")
    desc = desc + " [PhyFoodComp]"
    fg = pd.to_numeric(master.get("Food Group", 0), errors="coerce").fillna(0).astype(int)

    foods = pd.DataFrame({
        "fdc_id": master["fdc_id"].astype(int),
        "data_type": "foundation_food",
        "description": desc,
        "food_category_id": "9999",
        "food_category": fg.map(CATEGORY_MAP).fillna("Vegetables and Vegetable Products"),
    })

    # --- column-walking nutrient extraction ---
    print("[3/4] Walking nutrient columns + minting IDs...")
    nutrient_cols = [c for c in master.columns
                     if normalize_col(c) not in META_COLS
                     and normalize_col(c) != "fdc_id"]
    raw_to_id: dict[str, int] = {}
    extra_rows: list[dict] = []
    next_id = EXTRA_ID_START
    for c in nutrient_cols:
        n = normalize_col(c)
        if n in KNOWN_NUTRIENTS:
            raw_to_id[c] = KNOWN_NUTRIENTS[n]
        else:
            raw_to_id[c] = next_id
            extra_rows.append({
                "nutrient_id": next_id,
                "source_db": "phyfoodcomp",
                "source_column_raw": c,
                "source_column_norm": n,
                "note": "minted from PhyFoodComp v2 ingester",
            })
            next_id += 1
    n_known = sum(1 for c in nutrient_cols if normalize_col(c) in KNOWN_NUTRIENTS)
    print(f"  -> {len(nutrient_cols)} nutrient columns: "
          f"{n_known} known FDC, {len(extra_rows)} minted "
          f"({EXTRA_ID_START}-{next_id - 1})")

    # Long-form
    sub = master[["fdc_id"] + nutrient_cols].copy()
    for c in nutrient_cols:
        sub[c] = sub[c].apply(clean_amount)
    long_df = sub.melt(id_vars="fdc_id", var_name="raw_col", value_name="amount")
    long_df = long_df.dropna(subset=["amount"])
    long_df["nutrient_id"] = long_df["raw_col"].map(raw_to_id).astype(int)
    long_df = long_df[["fdc_id", "nutrient_id", "amount"]]

    print(f"  -> {len(long_df)} measurements after dropping NaN")
    print(f"  -> {long_df['nutrient_id'].nunique()} distinct nutrient ids")

    # --- write outputs ---
    print("[4/4] Writing parquet outputs...")
    out_food = f"{args.out_prefix}_food_injection.parquet"
    foods.drop_duplicates(subset=["fdc_id"]).to_parquet(out_food, index=False)

    out_bucket = Path(f"{args.out_prefix}_food_nutrient_bucketed")
    if out_bucket.exists():
        shutil.rmtree(out_bucket)
    out_bucket.mkdir(exist_ok=True)
    long_df["bucket"] = (long_df["nutrient_id"] % 256).astype(int)
    for b, g in long_df.groupby("bucket", sort=True):
        d = out_bucket / f"bucket={int(b)}"
        d.mkdir(exist_ok=True)
        g[["fdc_id", "nutrient_id", "amount"]].sort_values(["fdc_id", "nutrient_id"]) \
            .to_parquet(d / "phyfoodcomp_data.parquet", index=False)

    if extra_rows:
        pd.DataFrame(extra_rows).to_csv(f"{args.out_prefix}_extra_nutrient_map.tsv",
                                         sep="\t", index=False)

    print()
    print(f"[OK] foods                  : {len(foods)}")
    print(f"[OK] measurements           : {len(long_df)}")
    print(f"[OK] distinct nutrient ids  : {long_df['nutrient_id'].nunique()}")
    print(f"[OK] of which known FDC     : {n_known}")
    print(f"[OK] of which minted        : {len(extra_rows)}")


if __name__ == "__main__":
    main()
