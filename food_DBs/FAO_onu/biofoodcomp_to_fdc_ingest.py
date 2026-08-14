#!/usr/bin/env python3
"""
biofoodcomp_to_fdc_ingest.py

Ingest FAO/INFOODS BioFoodComp4.0 Excel into your FDC-like injection artifacts.

Fixes in this version:
- Robust header detection (handles banner rows like: "All food components are presented...")
- Uses descriptor row human labels when present
- Adds BioFoodComp-specific nutrient-name overrides (energy/protein/fat/carbs/sugars/starch/nitrogen/dry matter)
- Safer food-name column detection
- Better meta-column exclusion

Outputs:
1) <out_prefix>_food_injection.parquet
2) <out_prefix>_food_nutrient_bucketed/ (hive bucket partition)
3) <out_prefix>_unmapped_columns.tsv

python biofoodcomp_to_fdc_ingest.py   --excel BioFoodComp4.0.xlsx   --fdc_nutrient_csv /data/bac2food/nutrient.csv   --out_prefix biofoodcomp   --fdc_id_offset 82000000   --category_id 8888
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
UNIT_PAREN = re.compile(r"\s*\(.*?\)\s*$")


def norm(s: str) -> str:
    s = ("" if s is None else str(s)).lower()
    s = s.replace("µ", "u")
    s = PUNCT.sub(" ", s)
    s = WS.sub(" ", s).strip()
    return s


def detect_header_row(excel_path: str, sheet: str, scan_rows: int = 80) -> int | None:
    """
    Find row index where the true header lives.
    BioFoodComp often has banner text in row0, and header later.
    """
    raw = pd.read_excel(excel_path, sheet_name=sheet, header=None, nrows=scan_rows)
    # Acceptable variants for the food id column (normalized)
    key_variants = {
        "food item id",
        "fooditemid",
        "food id",
        "foodid",
        "food itemid",
        "food_item_id",
    }

    for i in range(len(raw)):
        row = [norm(x) for x in raw.iloc[i].tolist()]
        if any(x in key_variants for x in row):
            return i
    return None


def build_fdc_nutrient_name_map(fdc_nutrient_csv: str) -> Dict[str, int]:
    ndf = pd.read_csv(fdc_nutrient_csv)
    if "id" not in ndf.columns or "name" not in ndf.columns:
        raise ValueError(f"{fdc_nutrient_csv} must contain columns: id, name")
    out: Dict[str, int] = {}
    for rid, nm in zip(ndf["id"], ndf["name"]):
        try:
            nid = int(rid)
        except Exception:
            continue
        k = norm(nm)
        if k and k not in out:
            out[k] = nid
    return out


def add_biofoodcomp_overrides(overrides: dict[str, int], fdc_name2id: dict[str, int]) -> dict[str, int]:
    """
    Map common BioFoodComp human labels to FDC nutrient IDs.
    """

    def gid(fdc_name: str):
        return fdc_name2id.get(norm(fdc_name))

    def put(src: str, fdc_name: str):
        nid = gid(fdc_name)
        if nid is not None:
            overrides[norm(src)] = int(nid)

    # Energy
    put(
        "energy, total metabolizable; calculated from the energy-producing food components",
        "Energy",
    )
    put(
        "energy, total metabolizable; calculated from the energy-producing food components original as from source",
        "Energy",
    )
    put("energy, total metabolizable", "Energy")
    # gross energy is not a perfect match; keep optional
    put("Energy, gross", "Energy")

    # Dry matter -> Solids (FDC has Solids)
    put("Dry matter", "Solids")

    # Nitrogen
    put("Nitrogen, total", "Nitrogen")

    # Protein
    put("protein, total; calculated from total nitrogen", "Protein")
    put("Protein, total; calculated from protein nitrogen", "Protein")
    put("protein, total; calculated from protein nitrogen", "Protein")
    put("protein, total", "Protein")

    # Fat
    put("fat, total", "Total lipid (fat)")
    put("fat, total; derived by analysis using continuous extraction", "Total lipid (fat)")

    # Carbs
    put("carbohydrate, available; calculated by difference", "Carbohydrate, by difference")
    put("carbohydrate, total; calculated by difference", "Carbohydrate, by difference")
    put("Carbohydrate, available by weight", "Carbohydrate, by difference")
    put("Carbohydrate, available; expressed in monosaccharide equivalents", "Carbohydrate, by difference")
    put("carbohydrate, available", "Carbohydrate, by difference")

    # Sugars
    put("Sugar", "Sugars, total")
    put("Sugars, total; expressed in monosaccharide equivalents", "Sugars, total")
    put("sugars, total; expression unknown", "Sugars, total")
    put("Sugars, total", "Sugars, total")

    # Starch
    put("starch, available", "Starch")
    put("Starch, available; expressed in monosaccharide equivalents", "Starch")
    put("Starch, available", "Starch")

    return overrides


def sheet_to_category(sheet_name: str) -> str:
    s = norm(sheet_name)
    if "cereal" in s:
        return "Cereal Grains and Pasta"
    if "starchy roots" in s or "tuber" in s:
        return "Vegetables and Vegetable Products"
    if "legume" in s:
        return "Legumes and Legume Products"
    if "nuts" in s or "seeds" in s:
        return "Nut and Seed Products"
    if "vegetable" in s:
        return "Vegetables and Vegetable Products"
    if "fruit" in s:
        return "Fruits and Fruit Juices"
    if "meat" in s:
        return "Sausages and Luncheon Meats"
    if "egg" in s:
        return "Dairy and Egg Products"
    if "fish" in s or "seafood" in s or "shellfish" in s:
        return "Finfish and Shellfish Products"
    if "herbs" in s or "spices" in s:
        return "Spices and Herbs"
    if "milk" in s or "dairy" in s:
        return "Dairy and Egg Products"
    if "oil" in s or "fat" in s:
        return "Fats and Oils"
    if "beverage" in s or "drink" in s:
        return "Beverages"
    if "misc" in s:
        return "Meals, Entrees, and Side Dishes"
    return "Meals, Entrees, and Side Dishes"


def to_numeric_series(s: pd.Series) -> pd.Series:
    if s.dtype == object:
        x = s.astype(str)
        x = x.str.replace(",", ".", regex=False)
        x = x.str.replace(r"^\s*(<|tr|trace|n/a|na|nd|nan).*$", "", regex=True, flags=re.I)
        return pd.to_numeric(x, errors="coerce")
    return pd.to_numeric(s, errors="coerce")


def choose_name_column(cols: List[str]) -> str:
    """
    BioFoodComp often includes both "Foodname in English" and "Foodname in own language".
    Prefer English if present.
    """
    prefs = [
        "Foodname in English",
        "Food name in English",
        "Foodname in own language",
        "Food name",
        "Food",
    ]
    for p in prefs:
        if p in cols:
            return p

    # fuzzy fallback
    for c in cols:
        nc = norm(c)
        if "foodname" in nc and "english" in nc:
            return c
    for c in cols:
        if "foodname" in norm(c) or norm(c) in {"food name", "food"}:
            return c

    raise ValueError(f"Could not detect a food name column. Columns (first 40): {cols[:40]}")


def _coerce_food_id_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rename Food Item ID column to exactly 'Food Item ID' if found with variant casing/spacing.
    """
    if "Food Item ID" in df.columns:
        return df

    for c in df.columns:
        if norm(c) in {"food item id", "fooditemid", "food id", "foodid", "food_item_id"}:
            return df.rename(columns={c: "Food Item ID"})
    return df


def ingest_sheet(
    excel_path: str,
    sheet: str,
    fdc_name2id: Dict[str, int],
    overrides: Dict[str, int],
    fdc_id_offset: int,
    category_id: str,
    data_type: str,
    min_amount: float,
) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
    """
    Returns:
      foods_df, long_df, unmapped_columns (human labels)
    """
    hdr = detect_header_row(excel_path, sheet)
    if hdr is None:
        return pd.DataFrame(), pd.DataFrame(), []

    df = pd.read_excel(excel_path, sheet_name=sheet, header=hdr)
    if df.empty:
        return pd.DataFrame(), pd.DataFrame(), []

    # Sometimes a wrong hdr produces garbage columns; fallback: if we can't find food id column,
    # try scanning deeper and re-reading.
    df = _coerce_food_id_column(df)
    if "Food Item ID" not in df.columns:
        # fallback: scan more rows
        hdr2 = detect_header_row(excel_path, sheet, scan_rows=200)
        if hdr2 is None:
            return pd.DataFrame(), pd.DataFrame(), []
        df = pd.read_excel(excel_path, sheet_name=sheet, header=hdr2)
        df = _coerce_food_id_column(df)
        if "Food Item ID" not in df.columns:
            return pd.DataFrame(), pd.DataFrame(), []

    # Descriptor row: BioFoodComp uses first row with NaN Food Item ID as "human labels" row
    desc_row_idx = None
    if df["Food Item ID"].isna().any():
        desc_row_idx = df.index[df["Food Item ID"].isna()][0]

    col_to_human: Dict[str, str] = {}
    if desc_row_idx is not None:
        desc_row = df.loc[desc_row_idx]
        for c in df.columns:
            hv = str(desc_row.get(c, "")).strip()
            if hv and hv.lower() not in {"nan", "none"}:
                col_to_human[c] = hv

    # Drop rows without numeric Food Item ID (removes descriptor row + stray notes)
    df["Food Item ID_num"] = pd.to_numeric(df["Food Item ID"], errors="coerce")
    df = df.dropna(subset=["Food Item ID_num"]).copy()
    if df.empty:
        return pd.DataFrame(), pd.DataFrame(), []

    df["Food Item ID_num"] = df["Food Item ID_num"].astype(int)

    name_col = choose_name_column(list(df.columns))
    df["food_name"] = df[name_col].astype(str).fillna("").str.strip()
    df = df[df["food_name"].astype(str).str.len() > 0].copy()
    if df.empty:
        return pd.DataFrame(), pd.DataFrame(), []

    # Mint IDs: stable and non-colliding if you use different offsets per DB
    raise SystemExit(
        "This v1 ingest mints fdc_ids from a literal offset, which the block scheme "
        "no longer permits: ids are accessions held in fdc_id_map.tsv. Run biofoodcomp_to_fdc_ingest_v2.py "
        "instead, which allocates through food_DBs/fdc_blocks.py."
    )

    food_category = sheet_to_category(sheet)
    df["food_category"] = food_category

    # Meta columns to exclude from nutrients
    meta_like = {
        "Food Item ID",
        "Food Item ID_num",
        "fdc_id",
        "food_name",
        "food_category",
        name_col,
    }

    # Exclude obvious non-nutrient text/meta columns by normalized name patterns
    META_PATTERNS = [
        r"^foodname\b",
        r"^food name\b",
        r"\bpublication\b",
        r"\bcompiler\b",
        r"\bbiblio\b",
        r"\bcomments?\b",
        r"\bprocessing\b",
        r"\bcountry\b",
        r"\bregion\b",
        r"\bedible\b",
        r"^den\b",
        r"\bmethod\b",
        r"\bconversion factor\b",
        r"\bplant origin\b",
        r"\banimal origin\b",
        r"\bown language\b",
        r"\benglish\b",
        r"^n$",
    ]
    meta_rx = re.compile("|".join(META_PATTERNS), re.I)

    nutrient_cols = []
    for c in df.columns:
        if c in meta_like:
            continue
        if meta_rx.search(str(c)):
            continue
        nutrient_cols.append(c)

    # Map each nutrient column -> FDC nutrient_id
    col2nid: Dict[str, int] = {}
    unmapped: List[str] = []

    for c in nutrient_cols:
        human = col_to_human.get(c, str(c))
        key = norm(human)
        if not key or key == "nan":
            continue

        # normalize a bit
        key = key.replace("dietary fibre", "fiber").replace("fibre", "fiber")
        key = key.replace("oxalo acetic", "oxaloacetic")
        # strip trailing units
        base = UNIT_PAREN.sub("", human).strip()
        key_base = norm(base).replace("dietary fibre", "fiber").replace("fibre", "fiber")

        nid = fdc_name2id.get(key) or fdc_name2id.get(key_base)

        # override mapping (BioFoodComp quirks)
        if nid is None:
            nid = overrides.get(key) or overrides.get(key_base)

        if nid is None:
            unmapped.append(human)
            continue

        col2nid[c] = int(nid)

    mapped_cols = list(col2nid.keys())
    if not mapped_cols:
        return pd.DataFrame(), pd.DataFrame(), unmapped

    # Melt to long
    wide = df[["fdc_id", "food_name", "food_category"] + mapped_cols].copy()
    for c in mapped_cols:
        wide[c] = to_numeric_series(wide[c])

    long = (
        wide.melt(
            id_vars=["fdc_id", "food_name", "food_category"],
            value_vars=mapped_cols,
            var_name="biofoodcomp_col",
            value_name="amount",
        )
        .dropna(subset=["amount"])
        .copy()
    )

    long = long[long["amount"] > float(min_amount)]
    if long.empty:
        return pd.DataFrame(), pd.DataFrame(), unmapped

    long["nutrient_id"] = long["biofoodcomp_col"].map(col2nid).astype(int)
    long = long[["fdc_id", "nutrient_id", "amount"]].sort_values(["fdc_id", "nutrient_id"])

    foods = df[["fdc_id", "food_name", "food_category"]].drop_duplicates().copy()
    foods["data_type"] = data_type
    foods["description"] = foods["food_name"].astype(str) + " [BioFoodComp]"
    foods["food_category_id"] = str(category_id)
    foods = foods[["fdc_id", "data_type", "description", "food_category_id", "food_category"]]

    return foods, long, unmapped


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--excel", required=True)
    ap.add_argument("--fdc_nutrient_csv", required=True)
    ap.add_argument("--out_prefix", default="biofoodcomp")
    ap.add_argument("--fdc_id_offset", type=int, default=82000000)
    ap.add_argument("--category_id", default="8888")
    ap.add_argument("--data_type", default="foundation_food")
    ap.add_argument("--min_amount", type=float, default=0.0)
    ap.add_argument(
        "--include_sheets_regex",
        default=r"^\s*\d{2}\s",  # matches "01 Cereals", "02 Legumes", ...
        help="Regex of sheet names to ingest (default matches '01 ...', '02 ...' etc).",
    )
    ap.add_argument("--limit_sheets", type=int, default=0)
    args = ap.parse_args()

    xls = pd.ExcelFile(args.excel)
    rx = re.compile(args.include_sheets_regex)

    sheets = [s for s in xls.sheet_names if rx.search(s)]
    if args.limit_sheets and args.limit_sheets > 0:
        sheets = sheets[: args.limit_sheets]

    if not sheets:
        raise SystemExit(
            f"[BioFoodComp] No sheets matched regex: {args.include_sheets_regex}. "
            f"Available sheets: {xls.sheet_names}"
        )

    fdc_name2id = build_fdc_nutrient_name_map(args.fdc_nutrient_csv)
    overrides: Dict[str, int] = {}
    overrides = add_biofoodcomp_overrides(overrides, fdc_name2id)

    all_foods = []
    all_long = []
    all_unmapped = []

    for sh in sheets:
        foods, long, unmapped = ingest_sheet(
            excel_path=args.excel,
            sheet=sh,
            fdc_name2id=fdc_name2id,
            overrides=overrides,
            fdc_id_offset=args.fdc_id_offset,
            category_id=args.category_id,
            data_type=args.data_type,
            min_amount=args.min_amount,
        )

        if not foods.empty and not long.empty:
            all_foods.append(foods)
            all_long.append(long)

        all_unmapped.extend([(sh, u) for u in unmapped])

        mapped_n = int(long["nutrient_id"].nunique()) if not long.empty else 0
        print(
            f"[BioFoodComp] sheet='{sh}' foods={len(foods):,} rows={len(long):,} mapped_nutrients={mapped_n:,}",
            flush=True,
        )

    if not all_foods or not all_long:
        raise SystemExit(
            "[BioFoodComp] No foods produced. "
            "Try a broader --include_sheets_regex or inspect sheets."
        )

    out_prefix = Path(args.out_prefix)
    food_out = out_prefix.with_name(out_prefix.name + "_food_injection.parquet")
    bucket_dir = Path(out_prefix.name + "_food_nutrient_bucketed")
    bucket_dir.mkdir(parents=True, exist_ok=True)

    foods_df = pd.concat(all_foods, ignore_index=True).drop_duplicates(subset=["fdc_id"])
    long_df = pd.concat(all_long, ignore_index=True)

    foods_df.to_parquet(food_out, index=False)

    long_df["bucket"] = (long_df["nutrient_id"] % 256).astype(int)
    for b, g in long_df.groupby("bucket", sort=True):
        bdir = bucket_dir / f"bucket={int(b)}"
        bdir.mkdir(parents=True, exist_ok=True)
        g = g.sort_values(["fdc_id", "nutrient_id"])
        g[["fdc_id", "nutrient_id", "amount"]].to_parquet(bdir / "biofoodcomp_data.parquet", index=False)

    if all_unmapped:
        rep = out_prefix.with_name(out_prefix.name + "_unmapped_columns.tsv")
        (
            pd.DataFrame(all_unmapped, columns=["sheet", "unmapped_human_name"])
            .drop_duplicates()
            .to_csv(rep, sep="\t", index=False)
        )
        print(f"[BioFoodComp] wrote unmapped column report: {rep}", flush=True)

    print("[BioFoodComp] DONE", flush=True)
    print(f"  - food injection: {food_out}")
    print(f"  - nutrient bucketed: {bucket_dir}/bucket=*/biofoodcomp_data.parquet")


if __name__ == "__main__":
    main()