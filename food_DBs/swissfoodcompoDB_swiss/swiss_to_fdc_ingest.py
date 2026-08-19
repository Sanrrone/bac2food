#!/usr/bin/env python3
"""
swiss_to_fdc_ingest.py

Ingest Swiss Food Composition Database (Excel) into FDC-like files.
#usage: python swiss_to_fdc_ingest.py --excel "Swiss_food_composition_database.xlsx"   --out_prefix swissfcdb   --fdc_nutrient_csv /data/bac2food/nutrient.csv
Outputs
1) swiss_food_injection.parquet
   Columns: fdc_id, data_type, description, food_category_id, food_category

2) swiss_food_nutrient_bucketed/
   Hive-partitioned by bucket=nutrient_id%256 with parquet files containing:
   fdc_id, nutrient_id, amount

3) swiss_id_map.tsv
   Swiss ID -> minted fdc_id (deterministic by sort order)

Requirements:
- pandas, numpy, openpyxl (to read Excel)
- nutrient.csv from FDC with columns id, name, unit_name
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

# fdc_id allocation: fdc_blocks.py
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
import fdc_blocks


WS = re.compile(r"\s+")
PUNCT = re.compile(r"[^a-z0-9\s]+")


def norm(s: str) -> str:
    """Normalize string for comparison: lower case, replace µ with u, remove punctuation, collapse spaces."""
    s = ("" if s is None else str(s)).lower()
    s = s.replace("µ", "u")
    s = PUNCT.sub(" ", s)
    s = WS.sub(" ", s).strip()
    return s


def detect_header_row(excel_path: str, sheet: str, key: str = "ID", scan_rows: int = 30) -> int:
    """Find row index containing the given column name."""
    raw = pd.read_excel(excel_path, sheet_name=sheet, header=None, nrows=scan_rows)
    for i in range(len(raw)):
        row = raw.iloc[i].astype(str).tolist()
        if any(str(x).strip() == key for x in row):
            return i
    raise ValueError(f"Could not detect header row containing '{key}' in sheet '{sheet}'")


def parse_col_name_and_unit(col: str) -> Tuple[str, str]:
    """
    Extract base name and unit from a column like 'Fat, total (g)'.
    Returns (base_name, unit). If no parentheses, unit is empty string.
    """
    s = str(col).strip()
    if s.endswith(")") and "(" in s:
        base = s[: s.rfind("(")].strip()
        unit = s[s.rfind("(") + 1 : -1].strip()
        return base, unit
    return s, ""


def is_nutrient_column(col: str) -> bool:
    """Heuristic: if column ends with a unit in parentheses, treat as nutrient column."""
    base, unit = parse_col_name_and_unit(col)
    # common units: g, mg, µg, kJ, kcal, etc.
    return bool(unit and re.match(r"^[a-zµg]+$", unit, re.I))


def load_nutrient_map(fdc_nutrient_csv: str) -> Dict[Tuple[str, str], int]:
    """
    Load nutrient.csv and build a mapping from (normalized_name, normalized_unit) to nutrient_id.
    Also store a version without unit for fallback (where unit doesn't matter, e.g., sugars).
    """
    df = pd.read_csv(fdc_nutrient_csv)
    required = {"id", "name", "unit_name"}
    if not required.issubset(df.columns):
        raise ValueError(f"{fdc_nutrient_csv} must have columns: {required}")

    mapping = {}
    for _, row in df.iterrows():
        nid = int(row["id"])
        name_norm = norm(row["name"])
        unit_norm = norm(row["unit_name"]) if pd.notna(row["unit_name"]) else ""
        # store with unit
        mapping[(name_norm, unit_norm)] = nid
        # also store without unit (for fallback)
        # but be careful not to overwrite if same name has multiple units
        if (name_norm, "") not in mapping:
            mapping[(name_norm, "")] = nid
    return mapping


def build_manual_overrides() -> Dict[str, Tuple[str, str]]:
    """
    Return dict: normalized_swiss_base_name -> (target_fdc_name, unit) for manual mapping.
    This helps when exact matching fails due to naming differences.
    """
    overrides = {}

    # Macros
    overrides[norm("protein")] = ("Protein", "g")
    overrides[norm("fat total")] = ("Total lipid (fat)", "g")
    overrides[norm("carbohydrates available")] = ("Carbohydrate, by difference", "g")
    overrides[norm("dietary fibres")] = ("Fiber, total dietary", "g")
    overrides[norm("sugars")] = ("Sugars, total", "g")
    overrides[norm("starch")] = ("Starch", "g")
    overrides[norm("alcohol")] = ("Alcohol, ethyl", "g")
    overrides[norm("salt (nacl)")] = ("Salt", "g")  # or NaCl? FDC has "Salt" (nutrient_id 1093, unit g)
    overrides[norm("water")] = ("Water", "g")

    # Energy
    overrides[norm("energy kilojoules")] = ("Energy (kJ)", "kJ")
    overrides[norm("energy kilocalories")] = ("Energy (kcal)", "kcal")

    # Fatty acids
    overrides[norm("fatty acids saturated")] = ("Fatty acids, total saturated", "g")
    overrides[norm("fatty acids monounsaturated")] = ("Fatty acids, total monounsaturated", "g")
    overrides[norm("fatty acids polyunsaturated")] = ("Fatty acids, total polyunsaturated", "g")
    overrides[norm("linoleic acid")] = ("18:2 n-6 c,c (Linoleic)", "g")  # FDC name?
    overrides[norm("alpha linolenic acid")] = ("18:3 n-3 c,c,c (alpha-Linolenic)", "g")
    overrides[norm("eicosapentaenoic acid epa")] = ("20:5 n-3 (EPA)", "g")
    overrides[norm("docosahexaenoic acid dha")] = ("22:6 n-3 (DHA)", "g")
    overrides[norm("cholesterol")] = ("Cholesterol", "mg")

    # Vitamins
    overrides[norm("vitamin a activity re")] = ("Vitamin A, RAE", "µg")  # careful: RE vs RAE
    overrides[norm("vitamin a activity rae")] = ("Vitamin A, RAE", "µg")
    overrides[norm("retinol")] = ("Retinol", "µg")
    overrides[norm("beta carotene")] = ("Beta-carotene", "µg")
    overrides[norm("vitamin b1 thiamine")] = ("Thiamin", "mg")
    overrides[norm("vitamin b2 riboflavin")] = ("Riboflavin", "mg")
    overrides[norm("vitamin b6 pyridoxine")] = ("Vitamin B-6", "mg")
    overrides[norm("vitamin b12 cobalamin")] = ("Vitamin B-12", "µg")
    overrides[norm("niacin")] = ("Niacin", "mg")
    overrides[norm("folate")] = ("Folate, total", "µg")
    overrides[norm("pantothenic acid")] = ("Pantothenic acid", "mg")
    overrides[norm("vitamin c ascorbic acid")] = ("Vitamin C, total ascorbic acid", "mg")
    overrides[norm("vitamin d calciferol")] = ("Vitamin D (D2 + D3)", "µg")
    overrides[norm("vitamin e α tocopherol")] = ("Vitamin E (alpha-tocopherol)", "mg")

    # Minerals
    overrides[norm("potassium k")] = ("Potassium, K", "mg")
    overrides[norm("sodium na")] = ("Sodium, Na", "mg")
    overrides[norm("chloride cl")] = ("Chloride, Cl", "mg")
    overrides[norm("calcium ca")] = ("Calcium, Ca", "mg")
    overrides[norm("magnesium mg")] = ("Magnesium, Mg", "mg")
    overrides[norm("phosphorus p")] = ("Phosphorus, P", "mg")
    overrides[norm("iron fe")] = ("Iron, Fe", "mg")
    overrides[norm("iodide i")] = ("Iodine, I", "µg")  # FDC likely "Iodine"
    overrides[norm("zinc zn")] = ("Zinc, Zn", "mg")
    overrides[norm("selenium se")] = ("Selenium, Se", "µg")

    return overrides


def map_swiss_columns_to_fdc(
    cols: List[str],
    nutrient_map: Dict[Tuple[str, str], int],
    overrides: Dict[str, Tuple[str, str]],
) -> Dict[str, int]:
    """
    Return dict: original column name -> nutrient_id.
    """
    col2nid = {}
    for col in cols:
        if not is_nutrient_column(col):
            continue
        base, unit = parse_col_name_and_unit(col)
        base_norm = norm(base)
        unit_norm = norm(unit)

        # 1. exact match (name, unit)
        nid = nutrient_map.get((base_norm, unit_norm))

        # 2. try match without unit (if unit not critical)
        if nid is None:
            nid = nutrient_map.get((base_norm, ""))

        # 3. try manual override
        if nid is None and base_norm in overrides:
            target_name, target_unit = overrides[base_norm]
            target_name_norm = norm(target_name)
            target_unit_norm = norm(target_unit)
            nid = nutrient_map.get((target_name_norm, target_unit_norm))
            if nid is None:
                nid = nutrient_map.get((target_name_norm, ""))

        if nid is not None:
            col2nid[col] = nid

    return col2nid


def to_numeric_series(s: pd.Series) -> pd.Series:
    """
    Convert a series to numeric, handling special strings:
    - 'n.d.', '<0.5', 'tr', etc. -> NaN
    - replace comma with dot
    """
    if s.dtype == object:
        x = s.astype(str)
        x = x.str.replace(",", ".", regex=False)
        # remove any leading/trailing non-numeric characters like '<', '>', 'n.d.', 'tr'
        x = x.str.replace(r"^\s*(<|>|n\.d\.|nd|tr|trace|not detected|n/a|na).*$", "", regex=True, flags=re.I)
        return pd.to_numeric(x, errors="coerce")
    return pd.to_numeric(s, errors="coerce")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--excel", required=True, help="Swiss food composition database Excel file")
    ap.add_argument("--sheet", default="Generic Foods", help="Sheet name to ingest (e.g., 'Generic Foods' or 'Branded foods')")
    ap.add_argument("--fdc_nutrient_csv", default="/mnt/data/nutrient.csv", help="FDC nutrient.csv (id,name,unit_name)")
    ap.add_argument("--out_prefix", default="swiss", help="Prefix for output files")
    ap.add_argument("--category_id", default="8888", help="Custom food_category_id for all Swiss foods")
    ap.add_argument("--data_type", default="generic_food", help="data_type written into food injection parquet")
    ap.add_argument("--min_amount", type=float, default=0.0, help="Drop nutrients with amount <= this")
    ap.add_argument("--limit_rows", type=int, default=0, help="Debug: limit number of rows")
    args = ap.parse_args()

    excel_path = args.excel
    sheet = args.sheet

    # Detect header row (should contain "ID")
    header_row = detect_header_row(excel_path, sheet, key="ID")
    df = pd.read_excel(excel_path, sheet_name=sheet, header=header_row)
    df = df.dropna(how="all")

    # Required columns
    if "ID" not in df.columns:
        raise ValueError(f"Sheet '{sheet}' missing required column: ID")
    if "Name" not in df.columns:
        raise ValueError(f"Sheet '{sheet}' missing required column: Name")
    if "Category" not in df.columns:
        raise ValueError(f"Sheet '{sheet}' missing required column: Category")

    if args.limit_rows and args.limit_rows > 0:
        df = df.head(args.limit_rows)

    # Mint stable fdc_id based on sorted ID
    ids = df["ID"].astype(str).fillna("").tolist()
    uniq_ids = sorted(set([x.strip() for x in ids if x.strip()]))
    id2fdc = dict(zip(uniq_ids, fdc_blocks.assign("swiss", uniq_ids)))

    df["fdc_id"] = df["ID"].astype(str).map(lambda x: id2fdc.get(str(x).strip(), np.nan))
    df = df.dropna(subset=["fdc_id"])
    df["fdc_id"] = df["fdc_id"].astype(int)

    # Keep category from the sheet
    df["food_category"] = df["Category"].astype(str).fillna("Uncategorized")

    # Load nutrient mapping from FDC nutrient.csv
    nutrient_map = load_nutrient_map(args.fdc_nutrient_csv)
    overrides = build_manual_overrides()

    # Identify all columns that are nutrient columns (skip metadata like ID, Name, Category, Derivation, Source, etc.)
    all_cols = df.columns.tolist()
    # Remove columns we know are not nutrients
    exclude_patterns = re.compile(r"(Derivation of value|Source|Record has changed|ID|Name|Category|Density|Matrix unit|Synonyms|ID V 4.0|ID SwissFIR)", re.I)
    candidate_cols = [c for c in all_cols if not exclude_patterns.search(c)]

    col2nid = map_swiss_columns_to_fdc(candidate_cols, nutrient_map, overrides)

    mapped_cols = sorted(col2nid.keys())
    unmapped_cols = sorted([c for c in candidate_cols if c not in col2nid])

    print(f"[Swiss] rows={len(df):,} unique_foods={len(uniq_ids):,}")
    print(f"[Swiss] mapped_nutrient_columns={len(mapped_cols):,} unmapped_columns={len(unmapped_cols):,}")

    # Melt mapped nutrients
    wide = df[["fdc_id", "Name", "food_category"] + mapped_cols].copy()
    for c in mapped_cols:
        wide[c] = to_numeric_series(wide[c])

    long = wide.melt(
        id_vars=["fdc_id", "Name", "food_category"],
        value_vars=mapped_cols,
        var_name="swiss_col",
        value_name="amount",
    )

    long = long.dropna(subset=["amount"])
    if args.min_amount is not None:
        long = long[long["amount"] > float(args.min_amount)]

    long["nutrient_id"] = long["swiss_col"].map(col2nid).astype(int)
    long = long[["fdc_id", "nutrient_id", "amount"]].sort_values(["fdc_id", "nutrient_id"])

    # Output paths
    out_prefix = Path(args.out_prefix)
    food_out = out_prefix.with_name(out_prefix.name + "_food_injection.parquet")
    idmap_out = out_prefix.with_name(out_prefix.name + "_id_map.tsv")
    bucket_dir = Path(out_prefix.name + "_food_nutrient_bucketed")
    bucket_dir.mkdir(exist_ok=True)

    # Food injection parquet
    foods = df[["fdc_id", "Name", "food_category"]].drop_duplicates().copy()
    foods["data_type"] = args.data_type
    foods["description"] = foods["Name"].astype(str) + " [Swiss " + args.sheet + "]"
    foods["food_category_id"] = str(args.category_id)

    foods_out = foods[["fdc_id", "data_type", "description", "food_category_id", "food_category"]]
    foods_out.to_parquet(food_out, index=False)

    # ID map
    pd.DataFrame(
        {"Swiss_ID": list(id2fdc.keys()), "fdc_id": list(id2fdc.values())}
    ).sort_values("Swiss_ID").to_csv(idmap_out, sep="\t", index=False)

    # Bucketed parquet
    long["bucket"] = (long["nutrient_id"] % 256).astype(int)
    for b, g in long.groupby("bucket", sort=True):
        bdir = bucket_dir / f"bucket={int(b)}"
        bdir.mkdir(parents=True, exist_ok=True)
        g = g.sort_values(["fdc_id", "nutrient_id"])
        g[["fdc_id", "nutrient_id", "amount"]].to_parquet(bdir / "swiss_data.parquet", index=False)

    # Report unmapped columns
    if unmapped_cols:
        rep = out_prefix.with_name(out_prefix.name + "_unmapped_columns.tsv")
        pd.DataFrame({"unmapped_column": unmapped_cols}).to_csv(rep, sep="\t", index=False)
        print(f"[Swiss] wrote unmapped column report: {rep}")

    print("[Swiss] DONE")
    print(f"  - food injection: {food_out}")
    print(f"  - nutrient bucketed: {bucket_dir}/bucket=*/swiss_data.parquet")
    print(f"  - id map: {idmap_out}")


if __name__ == "__main__":
    main()
