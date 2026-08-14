#!/usr/bin/env python3
"""injest_japan_v2.py — Phase 8 re-parse of STFCJ (Japan FCT 2020).

v1 hits 6573 [STFCJ]-tagged rows in food.parquet (vs ~2478 unique foods) by
outer-joining 4 sub-files on item_no, but:
  - The fatty-acid file is not ingested at all.
  - Only ~50 nutrient names are mapped; ~120 columns dropped silently.
  - The 'Water' column is mapped to FDC id 1005 (= Carbohydrate, by difference).
    Water is 1051.
  - Per-food per-cell Python loop is O(n_foods × n_cols) — slow.
  - amount > 0 filter drops trace measurements.
  - Output is CSV not bucketed parquet, bypassing the bac2food pipeline.

v2:
  - Foods spine = main_*.xlsx 'Table' sheet (every Japanese FCT food is here).
  - Left-join AA / FA / OA Table 1 sub-files (per 100g EP — same basis as main).
  - Column-walking ingester; mint new IDs in 250001+ block.
  - Output bucketed parquet.
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

# fdc_id allocation is centralised in food_DBs/fdc_blocks.py. Never write a literal offset:
# ids are accessions looked up in fdc_id_map.tsv, so a food keeps its id across releases.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
import fdc_blocks


THIS = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS.parent))

EXTRA_ID_START = 250_001
# RETIRED: block bases live in food_DBs/fdc_blocks.py. This constant is kept out of
# the file entirely so it cannot drift from the real allocation.

# Inter-file dependency graph: spine + joins (all 'per 100 g EP' tables)
SPINE = ("main_1374049_1r12_1.xlsx", "Table")
JOINS = [
    ("aminoacid_1374049_2r11_1.xlsx", "Table 1(per 100 g EP)"),
    ("fatty_acid_1374049_3r11_1.xlsx", "Table 1(per 100 g EP)"),
    ("org_acid_1388558_4r12r.xlsx", "Table"),
]
ANNEXES = [
    ("main_1374049_1r12_1.xlsx", "Annex Table"),
    ("org_acid_1388558_4r12r.xlsx", "Annex (organic acids)"),
]

KNOWN_NUTRIENTS = {
    "water": 1051,
    "protein, calculated from reference nitrogen": 1003,
    "protein": 1003,
    "amino acids, total": 1003,
    "lipid": 1004,
    "fatty acids, total": 1004,
    "ash": 1007,
    "carbohydrate, available, total": 1005,
    "carbohydrate, available": 1005,
    "fiber, total dietary": 1079,
    "starch": 1009,
    "sucrose": 1010,
    "glucose": 1011,
    "fructose": 1012,
    "lactose": 1013,
    "maltose": 1014,
    "sodium": 1093,
    "potassium": 1092,
    "calcium": 1087,
    "magnesium": 1090,
    "phosphorus": 1091,
    "iron": 1089,
    "zinc": 1095,
    "copper": 1098,
    "manganese": 1101,
    "iodine": 1100,
    "selenium": 1103,
    "chromium": 1096,
    "molybdenum": 1102,
    "vitamin a": 1104,
    "retinol": 1105,
    "vitamin d": 1110,
    "vitamin e": 1109,
    "vitamin k": 1185,
    "thiamin": 1165,
    "riboflavin": 1166,
    "niacin": 1167,
    "vitamin b-6": 1175,
    "vitamin b6": 1175,
    "vitamin b-12": 1178,
    "vitamin b12": 1178,
    "folate": 1177,
    "pantothenic acid": 1170,
    "biotin": 1176,
    "vitamin c": 1162,
    "isoleucine": 1212, "leucine": 1213, "lysine": 1214, "methionine": 1215,
    "cystine": 1216, "phenylalanine": 1217, "tyrosine": 1218, "threonine": 1211,
    "tryptophan": 1210, "valine": 1219, "histidine": 1221, "arginine": 1220,
    "alanine": 1222, "aspartic acid": 1223, "glutamic acid": 1224, "glycine": 1225,
    "proline": 1226, "serine": 1227,
    "acetic acid": 1026, "lactic acid": 1038, "citric acid": 1033,
    "malic acid": 1043, "succinic acid": 1046, "tartaric acid": 1045,
    "caffeic acid": 1201, "ferulic acid": 1200,
    "cholesterol": 1253,
}

META_COLS = {
    "item no.", "item no", "item_no", "food group", "food group code",
    "food and description", "description", "indexed search", "remarks",
    "scientific name", "sci_name", "fdc_id", "food_category",
}

CAT_MAP = {
    1: "Cereal Grains and Pasta", 2: "Vegetables and Vegetable Products",
    3: "Sweets", 4: "Legumes and Legume Products", 5: "Nut and Seed Products",
    6: "Vegetables and Vegetable Products", 7: "Fruits and Fruit Juices",
    8: "Vegetables and Vegetable Products",  # mushrooms
    9: "Vegetables and Vegetable Products",  # algae
    10: "Finfish and Shellfish Products", 11: "Beef Products",
    12: "Dairy and Egg Products",
    13: "Fats and Oils", 14: "Sweets", 15: "Sweets",
    16: "Beverages", 17: "Spices and Herbs", 18: "Meals, Entrees, and Side Dishes",
}


def normalize(c: str) -> str:
    s = " ".join(str(c).split())
    return s.strip().lower()


def clean_amount(x) -> float:
    if pd.isna(x): return float("nan")
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x) if np.isfinite(x) else float("nan")
    s = str(x).strip()
    if not s or s in {"-", "ND", "N.D.", "nd"}:
        return float("nan")
    if s.lower() in {"tr", "trace"}:
        return 0.0
    # STFCJ wraps estimated values in parens: "(1.2)" -> 1.2
    s = re.sub(r"[\(\)]", "", s)
    s = re.sub(r"[^\d\.\-]", "", s)
    try:
        return float(s)
    except ValueError:
        return float("nan")


def find_header_row(excel_path: str, sheet: str, scan: int = 15) -> int:
    raw = pd.read_excel(excel_path, sheet_name=sheet, header=None, nrows=scan)
    for i in range(len(raw)):
        row = [normalize(x) for x in raw.iloc[i].tolist()]
        if "item no." in row or "item no" in row:
            return i
    return -1


def read_one(excel_path: Path, sheet: str) -> pd.DataFrame:
    if not excel_path.exists():
        print(f"  [skip] missing: {excel_path}")
        return pd.DataFrame()
    hdr = find_header_row(str(excel_path), sheet)
    if hdr < 0:
        print(f"  [skip] no Item No. header in {excel_path.name} :: {sheet}")
        return pd.DataFrame()
    df = pd.read_excel(excel_path, sheet_name=sheet, header=hdr)
    df.columns = [" ".join(str(c).split()) for c in df.columns]
    # rename
    for raw_name, std in (("Item No.", "item_no"), ("Item No", "item_no"),
                          ("Food and Description", "description"),
                          ("Food and  Description", "description")):
        if raw_name in df.columns:
            df = df.rename(columns={raw_name: std})
    if "item_no" not in df.columns:
        return pd.DataFrame()
    # STFCJ uses alphanumeric item_no like "01001" — keep as strings to
    # avoid losing leading zeros, but also coerce to int for the join key
    df["item_no"] = df["item_no"].astype(str).str.strip()
    df = df[df["item_no"].str.match(r"^\d+", na=False)].copy()
    df["item_no_int"] = df["item_no"].str.extract(r"(\d+)")[0].astype(int)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_prefix", default="japan")
    args = ap.parse_args()
    os.chdir(THIS)

    print("==================================================")
    print("STFCJ v2 ingester")
    print("==================================================")

    print(f"[1/5] foods spine = {SPINE[0]} :: {SPINE[1]!r}")
    spine = read_one(Path(SPINE[0]), SPINE[1])
    if spine.empty:
        sys.exit("Could not read spine")
    print(f"  -> spine has {len(spine)} foods, {len(spine.columns)} columns")

    # Append annex foods
    for f, s in ANNEXES:
        ann = read_one(Path(f), s)
        if not ann.empty:
            # New foods only — items not already in spine
            new = ann[~ann["item_no_int"].isin(spine["item_no_int"])]
            if not new.empty:
                spine = pd.concat([spine, new], ignore_index=True)
                print(f"  +{len(new)} foods from {f} :: {s}")

    print(f"[2/5] joining sub-files on item_no ...")
    master = spine.copy()
    for f, s in JOINS:
        df = read_one(Path(f), s)
        if df.empty:
            continue
        # Take only nutrient columns this file uniquely contributes
        new_cols = [c for c in df.columns
                    if c not in master.columns and c not in {"item_no", "item_no_int"}]
        # Left-join (preserves spine foods even when sub-file is missing them)
        master = master.merge(df[["item_no_int"] + new_cols],
                              on="item_no_int", how="left")
        print(f"     {f} :: {s} → +{len(new_cols)} columns "
              f"({len(df)} sub-file foods, {len(set(spine['item_no_int']) & set(df['item_no_int']))} joined)")

    # Scientific names (optional)
    sci_path = Path("scientific_names.xlsx")
    if sci_path.exists():
        try:
            sci = pd.read_excel(sci_path, sheet_name="Scientific name of food source", header=2)
            sci.columns = [" ".join(str(c).split()) for c in sci.columns]
            if "Item No." in sci.columns:
                sci = sci.rename(columns={"Item No.": "item_no",
                                          "Scientific name": "sci_name"})
            sci = sci.dropna(subset=["item_no"]).copy()
            sci["item_no_int"] = sci["item_no"].astype(str) \
                                      .str.extract(r"(\d+)")[0].astype("Int64")
            sci = sci[["item_no_int", "sci_name"]].dropna()
            master = master.merge(sci, on="item_no_int", how="left")
        except Exception as e:
            print(f"  [warn] sci_names: {e}")

    master["fdc_id"] = fdc_blocks.assign("stfcj", master["item_no_int"].astype(int))
    master["food_category"] = (master["item_no_int"] // 1000).map(CAT_MAP) \
                                .fillna("Meals, Entrees, and Side Dishes")

    # --- nutrient column walk + mint ---
    print("[3/5] walking nutrient columns...")
    nutrient_cols = []
    for c in master.columns:
        if normalize(c) in META_COLS: continue
        if c in {"item_no", "item_no_int", "fdc_id", "food_category", "sci_name", "description"}:
            continue
        # Quick numeric-ness check
        vals = master[c].dropna().head(200)
        if len(vals) == 0:
            continue
        nums = vals.apply(clean_amount)
        if nums.notna().sum() / max(1, len(vals)) < 0.4:
            continue
        nutrient_cols.append(c)

    raw_to_id: Dict[str, int] = {}
    extra_rows: List[Dict] = []
    next_id = EXTRA_ID_START
    for c in nutrient_cols:
        n = normalize(c)
        # Strip trailing unit '... (mg)' style
        n_base = re.sub(r"\s*\(.*?\)\s*$", "", n).strip()
        if n in KNOWN_NUTRIENTS:
            raw_to_id[c] = KNOWN_NUTRIENTS[n]
        elif n_base in KNOWN_NUTRIENTS:
            raw_to_id[c] = KNOWN_NUTRIENTS[n_base]
        else:
            raw_to_id[c] = next_id
            extra_rows.append({
                "nutrient_id": next_id,
                "source_db": "stfcj",
                "source_column_raw": c,
                "source_column_norm": n_base,
                "note": "minted from STFCJ v2 ingester",
            })
            next_id += 1
    UNIT_RE = re.compile(r"\s*\(.*?\)\s*$")
    n_known = sum(1 for c in nutrient_cols
                  if normalize(c) in KNOWN_NUTRIENTS
                  or UNIT_RE.sub("", normalize(c)).strip() in KNOWN_NUTRIENTS)
    print(f"  -> {len(nutrient_cols)} nutrient columns: "
          f"{n_known} known FDC, {len(extra_rows)} minted "
          f"({EXTRA_ID_START}-{next_id - 1})")

    # --- long-form ---
    print("[4/5] building long-form...")
    sub = master[["fdc_id"] + nutrient_cols].copy()
    for c in nutrient_cols:
        sub[c] = sub[c].apply(clean_amount)
    long = sub.melt(id_vars="fdc_id", var_name="raw_col", value_name="amount") \
              .dropna(subset=["amount"])
    long["nutrient_id"] = long["raw_col"].map(raw_to_id).astype(int)
    long = long[["fdc_id", "nutrient_id", "amount"]]
    print(f"  -> {len(long)} measurements, {long['nutrient_id'].nunique()} distinct nutrient IDs")

    # --- write outputs ---
    print("[5/5] writing parquet...")
    desc = master["description"].astype(str)
    if "sci_name" in master.columns:
        desc = desc + master["sci_name"].apply(lambda x: f" [{x}]" if pd.notna(x) else "")
    desc = desc + " [STFCJ]"
    foods = pd.DataFrame({
        "fdc_id": master["fdc_id"].astype(int),
        "data_type": "foundation_food",
        "description": desc,
        "food_category_id": "8181",
        "food_category": master["food_category"],
    }).drop_duplicates(subset=["fdc_id"])
    foods.to_parquet(f"{args.out_prefix}_food_injection.parquet", index=False)

    out_bucket = Path(f"{args.out_prefix}_food_nutrient_bucketed")
    if out_bucket.exists():
        shutil.rmtree(out_bucket)
    out_bucket.mkdir(exist_ok=True)
    long["bucket"] = (long["nutrient_id"] % 256).astype(int)
    for b, g in long.groupby("bucket", sort=True):
        d = out_bucket / f"bucket={int(b)}"
        d.mkdir(exist_ok=True)
        g[["fdc_id", "nutrient_id", "amount"]] \
            .sort_values(["fdc_id", "nutrient_id"]) \
            .to_parquet(d / "japan_data.parquet", index=False)

    if extra_rows:
        pd.DataFrame(extra_rows) \
            .to_csv(f"{args.out_prefix}_extra_nutrient_map.tsv", sep="\t", index=False)

    print()
    print(f"[OK] foods       : {len(foods)}")
    print(f"[OK] measurements: {len(long)}")
    print(f"[OK] nutrient ids: {long['nutrient_id'].nunique()}")
    print(f"[OK] minted ids  : {len(extra_rows)}")


if __name__ == "__main__":
    main()
