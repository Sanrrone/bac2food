#!/usr/bin/env python3
"""
afcd_to_fdc_ingest.py

Ingest Australian Food Composition Database (AFCD) "Nutrient profiles" Excel
into your FDC-like files:

Outputs
1) afcd_food_injection.parquet
   Columns: fdc_id, data_type, description, food_category_id, food_category
usage: python afcd_to_fdc_ingest.py --excel "Swiss_food_composition_database.xlsx"   --out_prefix afcd   --fdc_nutrient_csv /data/bac2food/nutrient.csv
2) afcd_food_nutrient_bucketed/
   Hive-partitioned by bucket=nutrient_id%256 with parquet files containing:
   fdc_id, nutrient_id, amount

3) afcd_id_map.tsv
   Public Food Key -> minted fdc_id (deterministic by sort order)

Notes
- Uses /mnt/data/nutrient.csv (or --fdc_nutrient_csv) to map nutrient columns -> FDC nutrient_id.
- Only ingests nutrient columns that can be mapped; others are ignored with a report.
- Category is heuristically inferred from food name (AFCD sheet lacks category text).
- Uses the sheet "All solids & liquids per 100 g" by default (contains liquids too).
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


WS = re.compile(r"\s+")
PUNCT = re.compile(r"[^a-z0-9\s]+")


def norm(s: str) -> str:
    s = ("" if s is None else str(s)).lower()
    s = s.replace("µ", "u")
    s = PUNCT.sub(" ", s)
    s = WS.sub(" ", s).strip()
    return s


def detect_header_row(excel_path: str, sheet: str, key: str = "Public Food Key", scan_rows: int = 30) -> int:
    raw = pd.read_excel(excel_path, sheet_name=sheet, header=None, nrows=scan_rows)
    for i in range(len(raw)):
        row = raw.iloc[i].astype(str).tolist()
        if any(str(x).strip() == key for x in row):
            return i
    raise ValueError(f"Could not detect header row containing '{key}' in sheet '{sheet}'")


def category_from_name(food_name: str) -> str:
    s = norm(food_name)

    # beverages
    if any(w in s for w in ["juice", "coffee", "tea", "soft drink", "soda", "beverage", "drink", "water"]):
        return "Beverages"
    if any(w in s for w in ["beer", "wine", "vodka", "gin", "rum", "whisky", "whiskey", "liqueur"]):
        return "Alcoholic Beverages"

    # dairy / eggs
    if any(w in s for w in ["milk", "yogurt", "yoghurt", "cheese", "cream", "butter", "whey", "kefir"]):
        return "Dairy and Egg Products"
    if any(w in s for w in ["egg", "omelet", "omelette"]):
        return "Dairy and Egg Products"

    # meats / fish
    if any(w in s for w in ["beef", "pork", "ham", "bacon", "sausage", "chicken", "turkey", "lamb", "veal"]):
        return "Sausages and Luncheon Meats"
    if any(w in s for w in ["fish", "salmon", "tuna", "sardine", "prawn", "shrimp", "shellfish", "oyster", "mussel"]):
        return "Finfish and Shellfish Products"

    # plant groups
    if any(w in s for w in ["lentil", "bean", "chickpea", "peas", "legume", "soy", "tofu", "tempeh"]):
        return "Legumes and Legume Products"
    if any(w in s for w in ["apple", "banana", "pear", "orange", "grape", "berry", "mango", "peach", "plum", "fruit"]):
        return "Fruits and Fruit Juices"
    if any(w in s for w in ["spinach", "lettuce", "broccoli", "carrot", "tomato", "onion", "garlic", "vegetable", "cabbage"]):
        return "Vegetables and Vegetable Products"
    if any(w in s for w in ["rice", "oat", "wheat", "rye", "barley", "pasta", "noodle", "bread", "flour", "cereal", "grain"]):
        return "Cereal Grains and Pasta"
    if any(w in s for w in ["almond", "walnut", "hazelnut", "peanut", "cashew", "pistachio", "flax", "chia", "seed", "nuts"]):
        return "Nut and Seed Products"

    # fats / oils
    if any(w in s for w in ["oil", "olive oil", "canola", "sunflower oil", "margarine", "shortening", "lard"]):
        return "Fats and Oils"

    # sweets / baked / snacks
    if any(w in s for w in ["cake", "cookie", "biscuit", "pastry", "pie", "donut", "doughnut", "muffin", "bread roll"]):
        return "Baked Products"
    if any(w in s for w in ["chocolate", "candy", "sweet", "syrup", "honey", "jam", "jelly", "ice cream"]):
        return "Sweets"
    if any(w in s for w in ["chip", "crisps", "snack"]):
        return "Snacks"

    # default
    return "Meals, Entrees, and Side Dishes"


def build_fdc_nutrient_name_map(fdc_nutrient_csv: str) -> Dict[str, int]:
    ndf = pd.read_csv(fdc_nutrient_csv)
    if "id" not in ndf.columns or "name" not in ndf.columns:
        raise ValueError(f"{fdc_nutrient_csv} must have columns: id, name")
    m: Dict[str, int] = {}
    for rid, nm in zip(ndf["id"], ndf["name"]):
        try:
            nid = int(rid)
        except Exception:
            continue
        key = norm(nm)
        if key and key not in m:
            m[key] = nid
    return m


def parse_col_base_and_unit(col: str) -> Tuple[str, str]:
    """
    Examples:
      "Protein (g)" -> ("Protein", "g")
      "Energy, without dietary fibre (kJ)" -> ("Energy, without dietary fibre", "kJ")
    """
    s = str(col).strip()
    if s.endswith(")") and "(" in s:
        base = s[: s.rfind("(")].strip()
        unit = s[s.rfind("(") + 1 : -1].strip()
        return base, unit
    return s, ""


def build_manual_overrides(fdc_name2id: Dict[str, int]) -> Dict[str, int]:
    """
    Map common AFCD base names (normalized) to FDC nutrient IDs.
    Only applied if exact matching fails.
    """
    def gid(name: str) -> int | None:
        return fdc_name2id.get(norm(name))

    overrides: Dict[str, int] = {}

    # Macros
    for k, v in [
        ("protein", "Protein"),
        ("total fat", "Total lipid (fat)"),
        ("fat total", "Total lipid (fat)"),
        ("carbohydrate", "Carbohydrate, by difference"),
        ("dietary fibre", "Fiber, total dietary"),
        ("total dietary fibre", "Fiber, total dietary"),
        ("total sugars", "Sugars, total"),
        ("sugars", "Sugars, total"),
        ("starch", "Starch"),
        ("alcohol", "Alcohol, ethyl"),
    ]:
        nid = gid(v)
        if nid is not None:
            overrides[norm(k)] = nid

    # Energy
    for k, v in [
        ("energy", "Energy"),
        ("energy without dietary fibre", "Energy"),
    ]:
        nid = gid(v)
        if nid is not None:
            overrides[norm(k)] = nid

    # Minerals
    for k, v in [
        ("calcium", "Calcium, Ca"),
        ("iron", "Iron, Fe"),
        ("magnesium", "Magnesium, Mg"),
        ("phosphorus", "Phosphorus, P"),
        ("potassium", "Potassium, K"),
        ("sodium", "Sodium, Na"),
        ("zinc", "Zinc, Zn"),
        ("copper", "Copper, Cu"),
        ("selenium", "Selenium, Se"),
    ]:
        nid = gid(v)
        if nid is not None:
            overrides[norm(k)] = nid

    # Vitamins (some names differ)
    for k, v in [
        ("vitamin c", "Vitamin C, total ascorbic acid"),
        ("thiamin", "Thiamin"),
        ("riboflavin", "Riboflavin"),
        ("niacin", "Niacin"),
        ("vitamin b 6", "Vitamin B-6"),
        ("vitamin b6", "Vitamin B-6"),
        ("folate", "Folate, total"),
        ("vitamin b 12", "Vitamin B-12"),
        ("vitamin b12", "Vitamin B-12"),
        ("vitamin a", "Vitamin A, RAE"),
        ("vitamin d", "Vitamin D (D2 + D3)"),
        ("vitamin e", "Vitamin E (alpha-tocopherol)"),
    ]:
        nid = gid(v)
        if nid is not None:
            overrides[norm(k)] = nid

    # Fats
    for k, v in [
        ("saturated fat", "Fatty acids, total saturated"),
        ("monounsaturated fat", "Fatty acids, total monounsaturated"),
        ("polyunsaturated fat", "Fatty acids, total polyunsaturated"),
        ("cholesterol", "Cholesterol"),
    ]:
        nid = gid(v)
        if nid is not None:
            overrides[norm(k)] = nid

    return overrides


def map_afcd_columns_to_fdc_nutrients(
    cols: List[str],
    fdc_name2id: Dict[str, int],
    overrides: Dict[str, int],
) -> Dict[str, int]:
    """
    Return dict: original_column_name -> nutrient_id
    """
    col2nid: Dict[str, int] = {}

    for c in cols:
        base, _unit = parse_col_base_and_unit(c)
        base_n = norm(base)

        # direct exact match
        nid = fdc_name2id.get(base_n)

        # try small normalization tweaks
        if nid is None:
            base_n2 = base_n.replace("dietary fibre", "fiber").replace("fibre", "fiber")
            nid = fdc_name2id.get(base_n2)

        # manual override
        if nid is None:
            nid = overrides.get(base_n)

        if nid is not None:
            col2nid[c] = int(nid)

    return col2nid


def to_numeric_series(s: pd.Series) -> pd.Series:
    """
    Robust numeric conversion:
    - handles strings with '<', 'tr', etc -> NaN
    - handles commas
    """
    if s.dtype == object:
        x = s.astype(str)
        x = x.str.replace(",", ".", regex=False)
        x = x.str.replace(r"^\s*(<|tr|trace|n/a|na|nd|nan).*$", "", regex=True, flags=re.I)
        return pd.to_numeric(x, errors="coerce")
    return pd.to_numeric(s, errors="coerce")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--excel", required=True, help="AFCD Nutrient profiles Excel file")
    ap.add_argument("--sheet", default="All solids & liquids per 100 g", help="Sheet to ingest")
    ap.add_argument("--fdc_nutrient_csv", default="/mnt/data/nutrient.csv", help="FDC nutrient.csv (id,name,unit_name,...) used to map column->nutrient_id")
    ap.add_argument("--out_prefix", default="afcd", help="Prefix for output files")
    ap.add_argument("--fdc_id_offset", type=int, default=70000000, help="Start of minted FDC IDs for AFCD foods")
    ap.add_argument("--category_id", default="7777", help="Custom food_category_id for AFCD injection foods")
    ap.add_argument("--data_type", default="foundation_food", help="data_type written into food injection parquet")
    ap.add_argument("--min_amount", type=float, default=0.0, help="Drop nutrients with amount <= this")
    ap.add_argument("--limit_rows", type=int, default=0, help="Debug: limit number of foods")
    args = ap.parse_args()

    excel_path = args.excel
    sheet = args.sheet

    header_row = detect_header_row(excel_path, sheet, key="Public Food Key")
    df = pd.read_excel(excel_path, sheet_name=sheet, header=header_row)
    df = df.dropna(how="all")

    # Required AFCD columns
    required = ["Public Food Key", "Food Name"]
    for r in required:
        if r not in df.columns:
            raise ValueError(f"Sheet '{sheet}' missing required column: {r}")

    if args.limit_rows and args.limit_rows > 0:
        df = df.head(args.limit_rows)

    # Mint stable AFCD fdc_id by sorted Public Food Key
    keys = df["Public Food Key"].astype(str).fillna("").tolist()
    uniq_keys = sorted(set([k.strip() for k in keys if k.strip()]))
    raise SystemExit(
        "This v1 ingest mints fdc_ids from a literal offset, which the block scheme "
        "no longer permits: ids are accessions held in fdc_id_map.tsv. Run afcd_to_fdc_ingest_v2.py "
        "instead, which allocates through food_DBs/fdc_blocks.py."
    )

    df["fdc_id"] = df["Public Food Key"].astype(str).map(lambda x: key2fdc.get(str(x).strip(), np.nan))
    df = df.dropna(subset=["fdc_id"])
    df["fdc_id"] = df["fdc_id"].astype(int)

    # Category from name heuristic
    df["food_category"] = df["Food Name"].astype(str).map(category_from_name)

    # Build nutrient mapping
    fdc_name2id = build_fdc_nutrient_name_map(args.fdc_nutrient_csv)
    overrides = build_manual_overrides(fdc_name2id)

    non_meta_cols = [c for c in df.columns if c not in {"Public Food Key", "Food Name", "Classification", "fdc_id", "food_category"}]
    col2nid = map_afcd_columns_to_fdc_nutrients(non_meta_cols, fdc_name2id, overrides)

    mapped_cols = sorted(col2nid.keys())
    unmapped_cols = sorted([c for c in non_meta_cols if c not in col2nid])

    print(f"[AFCD] rows={len(df):,} unique_foods={len(uniq_keys):,}")
    print(f"[AFCD] mapped_nutrient_columns={len(mapped_cols):,} unmapped_columns={len(unmapped_cols):,}")

    # Melt mapped nutrients
    wide = df[["fdc_id", "Food Name", "food_category"] + mapped_cols].copy()
    for c in mapped_cols:
        wide[c] = to_numeric_series(wide[c])

    long = wide.melt(
        id_vars=["fdc_id", "Food Name", "food_category"],
        value_vars=mapped_cols,
        var_name="afcd_col",
        value_name="amount",
    )

    long = long.dropna(subset=["amount"])
    if args.min_amount is not None:
        long = long[long["amount"] > float(args.min_amount)]

    long["nutrient_id"] = long["afcd_col"].map(col2nid).astype(int)
    long = long[["fdc_id", "nutrient_id", "amount"]].sort_values(["fdc_id", "nutrient_id"])

    # Output paths
    out_prefix = Path(args.out_prefix)
    food_out = out_prefix.with_name(out_prefix.name + "_food_injection.parquet")
    idmap_out = out_prefix.with_name(out_prefix.name + "_id_map.tsv")
    bucket_dir = Path(out_prefix.name + "_food_nutrient_bucketed")
    bucket_dir.mkdir(exist_ok=True)

    # Food injection parquet
    foods = df[["fdc_id", "Food Name", "food_category"]].drop_duplicates().copy()
    foods["data_type"] = args.data_type
    foods["description"] = foods["Food Name"].astype(str) + " [AFCD]"
    foods["food_category_id"] = str(args.category_id)

    foods_out = foods[["fdc_id", "data_type", "description", "food_category_id", "food_category"]]
    foods_out.to_parquet(food_out, index=False)

    # ID map
    pd.DataFrame(
        {"Public Food Key": list(key2fdc.keys()), "fdc_id": list(key2fdc.values())}
    ).sort_values("Public Food Key").to_csv(idmap_out, sep="\t", index=False)

    # Bucketed parquet
    long["bucket"] = (long["nutrient_id"] % 256).astype(int)
    for b, g in long.groupby("bucket", sort=True):
        bdir = bucket_dir / f"bucket={int(b)}"
        bdir.mkdir(parents=True, exist_ok=True)
        g = g.sort_values(["fdc_id", "nutrient_id"])
        g[["fdc_id", "nutrient_id", "amount"]].to_parquet(bdir / "afcd_data.parquet", index=False)

    # Report unmapped columns (for future extension)
    if unmapped_cols:
        rep = out_prefix.with_name(out_prefix.name + "_unmapped_columns.tsv")
        pd.DataFrame({"unmapped_column": unmapped_cols}).to_csv(rep, sep="\t", index=False)
        print(f"[AFCD] wrote unmapped column report: {rep}")

    print("[AFCD] DONE")
    print(f"  - food injection: {food_out}")
    print(f"  - nutrient bucketed: {bucket_dir}/bucket=*/afcd_data.parquet")
    print(f"  - id map: {idmap_out}")


if __name__ == "__main__":
    main()
