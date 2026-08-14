#!/usr/bin/env python3
"""biofoodcomp_to_fdc_ingest_v2.py — Phase 8 re-parse of FAO BioFoodComp 4.0.

v1 already correctly extracts foods (10133 — matches the source docs) and has
working header/descriptor-row detection. v1's defects are nutrient-column:
~1401 of the ~1500 raw columns get dropped because their normalized human
label doesn't exactly match a name in /data/bac2food/nutrient.csv. The v1
parser writes the dropped names to biofoodcomp_unmapped_columns.tsv but never
mints IDs for them.

v2 changes (minimal patch over v1):
  1. Mint nutrient_ids in the 220001+ block for any unmapped human-labeled
     column whose values are mostly-numeric. The minted IDs go into
     biofoodcomp_extra_nutrient_map.tsv along with the source column + sheet.
  2. Drop the `amount > min_amount` filter — preserve zero measurements
     (a legitimate "this food was tested and contains 0 mg sodium" reading is
     biologically meaningful and shouldn't be filtered).
  3. Keep the v1 descriptor-row + sheet-walking + name-column logic intact.

Outputs (overwrite v1's):
  biofoodcomp_food_injection.parquet
  biofoodcomp_food_nutrient_bucketed/bucket=*/biofoodcomp_data.parquet
  biofoodcomp_unmapped_columns.tsv         — still written for transparency
  biofoodcomp_extra_nutrient_map.tsv       — new: the minted nutrient IDs
"""
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

# fdc_id allocation is centralised in food_DBs/fdc_blocks.py. Never write a literal offset:
# ids are accessions looked up in fdc_id_map.tsv, so a food keeps its id across releases.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
import fdc_blocks


WS = re.compile(r"\s+")
PUNCT = re.compile(r"[^a-z0-9\s]+")
UNIT_PAREN = re.compile(r"\s*\(.*?\)\s*$")

EXTRA_ID_START = 220_001


def norm(s: str) -> str:
    s = ("" if s is None else str(s)).lower().replace("µ", "u")
    s = PUNCT.sub(" ", s)
    return WS.sub(" ", s).strip()


def detect_header_row(excel_path: str, sheet: str, scan_rows: int = 80) -> int | None:
    raw = pd.read_excel(excel_path, sheet_name=sheet, header=None, nrows=scan_rows)
    key_variants = {"food item id", "fooditemid", "food id", "foodid",
                    "food itemid", "food_item_id"}
    for i in range(len(raw)):
        row = [norm(x) for x in raw.iloc[i].tolist()]
        if any(x in key_variants for x in row):
            return i
    return None


def build_fdc_nutrient_name_map(fdc_nutrient_csv: str) -> Dict[str, int]:
    ndf = pd.read_csv(fdc_nutrient_csv)
    out: Dict[str, int] = {}
    for rid, nm in zip(ndf["id"], ndf["name"]):
        try:
            nid = int(rid)
        except (ValueError, TypeError):
            continue
        k = norm(nm)
        if k and k not in out:
            out[k] = nid
    return out


def add_biofoodcomp_overrides(overrides: dict, fdc_name2id: dict) -> dict:
    def put(src: str, fdc_name: str):
        nid = fdc_name2id.get(norm(fdc_name))
        if nid is not None:
            overrides[norm(src)] = int(nid)
    put("energy, total metabolizable", "Energy")
    put("energy, total metabolizable; calculated from the energy-producing food components", "Energy")
    put("Energy, gross", "Energy")
    put("Dry matter", "Solids")
    put("Nitrogen, total", "Nitrogen")
    put("protein, total; calculated from total nitrogen", "Protein")
    put("protein, total", "Protein")
    put("fat, total", "Total lipid (fat)")
    put("fat, total; derived by analysis using continuous extraction", "Total lipid (fat)")
    put("carbohydrate, available; calculated by difference", "Carbohydrate, by difference")
    put("carbohydrate, available", "Carbohydrate, by difference")
    put("Sugars, total", "Sugars, total")
    put("Sugars, total; expressed in monosaccharide equivalents", "Sugars, total")
    put("Starch, available", "Starch")
    put("Starch, available; expressed in monosaccharide equivalents", "Starch")
    return overrides


def sheet_to_category(sheet_name: str) -> str:
    s = norm(sheet_name)
    if "cereal" in s: return "Cereal Grains and Pasta"
    if "starchy roots" in s or "tuber" in s: return "Vegetables and Vegetable Products"
    if "legume" in s: return "Legumes and Legume Products"
    if "nuts" in s or "seeds" in s: return "Nut and Seed Products"
    if "vegetable" in s: return "Vegetables and Vegetable Products"
    if "fruit" in s: return "Fruits and Fruit Juices"
    if "meat" in s: return "Sausages and Luncheon Meats"
    if "egg" in s: return "Dairy and Egg Products"
    if "fish" in s or "seafood" in s or "shellfish" in s: return "Finfish and Shellfish Products"
    if "herbs" in s or "spices" in s: return "Spices and Herbs"
    if "milk" in s or "dairy" in s: return "Dairy and Egg Products"
    if "oil" in s or "fat" in s: return "Fats and Oils"
    if "beverage" in s or "drink" in s: return "Beverages"
    return "Meals, Entrees, and Side Dishes"


def to_numeric_series(s: pd.Series) -> pd.Series:
    if s.dtype == object:
        x = s.astype(str).str.replace(",", ".", regex=False)
        x = x.str.replace(r"^\s*(<|tr|trace|n/a|na|nd|nan).*$", "", regex=True, flags=re.I)
        return pd.to_numeric(x, errors="coerce")
    return pd.to_numeric(s, errors="coerce")


def choose_name_column(cols: List[str]) -> str:
    for p in ("Foodname in English", "Food name in English", "Foodname in own language",
              "Food name", "Food"):
        if p in cols:
            return p
    for c in cols:
        nc = norm(c)
        if "foodname" in nc and "english" in nc: return c
    for c in cols:
        if "foodname" in norm(c) or norm(c) in {"food name", "food"}: return c
    raise ValueError(f"No food name col in: {cols[:40]}")


def _coerce_food_id_column(df: pd.DataFrame) -> pd.DataFrame:
    if "Food Item ID" in df.columns:
        return df
    for c in df.columns:
        if norm(c) in {"food item id", "fooditemid", "food id", "foodid", "food_item_id"}:
            return df.rename(columns={c: "Food Item ID"})
    return df


# Meta-column exclusion patterns: these columns are never nutrients regardless
# of whether a descriptor row provided a label.
META_PATTERNS = re.compile("|".join([
    r"^foodname\b", r"^food name\b", r"\bpublication\b", r"\bcompiler\b",
    r"\bbiblio\b", r"\bcomments?\b", r"\bprocessing\b", r"\bcountry\b",
    r"\bregion\b", r"\bedible\b", r"^den\b", r"\bmethod\b",
    r"\bconversion factor\b", r"\bplant origin\b", r"\banimal origin\b",
    r"\bown language\b", r"\benglish\b", r"^n$",
    r"\bcultivar\b", r"\bvariety\b", r"\baccession\b", r"\bseason\b",
    r"\bsubspecies\b", r"\bspecies\b", r"\bsubgroup\b", r"^type$",
    r"\bsamples?\b",
]), re.I)


def is_mostly_numeric(s: pd.Series, min_frac: float = 0.4) -> bool:
    """Treat a column as a nutrient if ≥ min_frac of its non-NaN values parse as float."""
    vals = s.dropna()
    if len(vals) == 0:
        return False
    nums = to_numeric_series(vals.head(200))
    n_ok = nums.notna().sum()
    return n_ok / max(1, min(200, len(vals))) >= min_frac


def ingest_sheet_v2(
    excel_path: str,
    sheet: str,
    fdc_name2id: Dict[str, int],
    overrides: Dict[str, int],
    minted: Dict[str, int],   # raw_key (norm'd human label) -> minted id
    next_id_box: list[int],
    fdc_id_offset: int,
    category_id: str,
    data_type: str,
) -> Tuple[pd.DataFrame, pd.DataFrame, List[str], List[Dict]]:
    hdr = detect_header_row(excel_path, sheet)
    if hdr is None:
        return pd.DataFrame(), pd.DataFrame(), [], []
    df = pd.read_excel(excel_path, sheet_name=sheet, header=hdr)
    if df.empty:
        return pd.DataFrame(), pd.DataFrame(), [], []
    df = _coerce_food_id_column(df)
    if "Food Item ID" not in df.columns:
        return pd.DataFrame(), pd.DataFrame(), [], []

    # Descriptor row = first row where Food Item ID is NaN; values are human labels
    col_to_human: Dict[str, str] = {}
    if df["Food Item ID"].isna().any():
        desc_idx = df.index[df["Food Item ID"].isna()][0]
        for c in df.columns:
            hv = str(df.loc[desc_idx, c]).strip()
            if hv and hv.lower() not in {"nan", "none"}:
                col_to_human[c] = hv

    df["Food Item ID_num"] = pd.to_numeric(df["Food Item ID"], errors="coerce")
    df = df.dropna(subset=["Food Item ID_num"]).copy()
    if df.empty:
        return pd.DataFrame(), pd.DataFrame(), [], []
    df["Food Item ID_num"] = df["Food Item ID_num"].astype(int)
    name_col = choose_name_column(list(df.columns))
    df["food_name"] = df[name_col].astype(str).str.strip()
    df = df[df["food_name"].str.len() > 0].copy()
    if df.empty:
        return pd.DataFrame(), pd.DataFrame(), [], []
    df["fdc_id"] = fdc_blocks.assign("biofoodcomp", df["Food Item ID_num"].astype(int))
    df["food_category"] = sheet_to_category(sheet)

    meta_like = {"Food Item ID", "Food Item ID_num", "fdc_id", "food_name",
                 "food_category", name_col}

    # Identify candidate nutrient columns
    candidate_cols = []
    for c in df.columns:
        if c in meta_like: continue
        if META_PATTERNS.search(str(c)): continue
        if not is_mostly_numeric(df[c]): continue
        candidate_cols.append(c)

    col2nid: Dict[str, int] = {}
    unmapped_in_sheet: List[str] = []
    new_minted_in_sheet: List[Dict] = []

    for c in candidate_cols:
        human = col_to_human.get(c, str(c))
        key = norm(human).replace("dietary fibre", "fiber").replace("fibre", "fiber") \
                         .replace("oxalo acetic", "oxaloacetic")
        if not key or key == "nan":
            continue
        base = UNIT_PAREN.sub("", human).strip()
        key_base = norm(base).replace("dietary fibre", "fiber").replace("fibre", "fiber")
        nid = fdc_name2id.get(key) or fdc_name2id.get(key_base) \
              or overrides.get(key) or overrides.get(key_base)
        if nid is None:
            # Mint or reuse a previously-minted ID for this normalized human label
            stable_key = key_base or key
            if stable_key in minted:
                nid = minted[stable_key]
            else:
                nid = next_id_box[0]
                minted[stable_key] = nid
                next_id_box[0] += 1
                new_minted_in_sheet.append({
                    "nutrient_id": nid,
                    "source_db": "biofoodcomp",
                    "source_human_label": human,
                    "source_column_raw": str(c),
                    "source_sheet": sheet,
                    "note": "minted (no FDC nutrient.csv match)",
                })
            unmapped_in_sheet.append(human)
        col2nid[c] = int(nid)

    mapped_cols = list(col2nid.keys())
    if not mapped_cols:
        return pd.DataFrame(), pd.DataFrame(), unmapped_in_sheet, new_minted_in_sheet

    wide = df[["fdc_id", "food_name", "food_category"] + mapped_cols].copy()
    for c in mapped_cols:
        wide[c] = to_numeric_series(wide[c])
    long = (wide.melt(id_vars=["fdc_id", "food_name", "food_category"],
                       value_vars=mapped_cols,
                       var_name="raw_col", value_name="amount")
                .dropna(subset=["amount"]))    # only drop NaN; keep zeros
    if long.empty:
        return pd.DataFrame(), pd.DataFrame(), unmapped_in_sheet, new_minted_in_sheet
    long["nutrient_id"] = long["raw_col"].map(col2nid).astype(int)
    long = long[["fdc_id", "nutrient_id", "amount"]] \
              .sort_values(["fdc_id", "nutrient_id"])

    foods = df[["fdc_id", "food_name", "food_category"]].drop_duplicates().copy()
    foods["data_type"] = data_type
    foods["description"] = foods["food_name"] + " [BioFoodComp]"
    foods["food_category_id"] = str(category_id)
    foods = foods[["fdc_id", "data_type", "description",
                   "food_category_id", "food_category"]]

    return foods, long, unmapped_in_sheet, new_minted_in_sheet


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--excel", default="BioFoodComp4.0.xlsx")
    ap.add_argument("--fdc_nutrient_csv", default="/data/bac2food/nutrient.csv")
    ap.add_argument("--out_prefix", default="biofoodcomp")
    ap.add_argument("--fdc_id_offset", type=int, default=82_000_000)
    ap.add_argument("--category_id", default="8888")
    ap.add_argument("--data_type", default="foundation_food")
    ap.add_argument("--include_sheets_regex", default=r"^\s*\d{2}[\s_]")
    args = ap.parse_args()

    import os; os.chdir(Path(__file__).parent)
    xls = pd.ExcelFile(args.excel)
    rx = re.compile(args.include_sheets_regex)
    sheets = [s for s in xls.sheet_names if rx.search(s)]
    print(f"[BioFoodComp v2] {len(sheets)} sheets:")
    for s in sheets:
        print(f"  {s}")

    fdc_name2id = build_fdc_nutrient_name_map(args.fdc_nutrient_csv)
    overrides = add_biofoodcomp_overrides({}, fdc_name2id)
    minted: Dict[str, int] = {}
    next_id_box = [EXTRA_ID_START]

    all_foods, all_long, all_unmapped, all_minted = [], [], [], []
    for sh in sheets:
        foods, long, unmapped, minted_in_sheet = ingest_sheet_v2(
            args.excel, sh, fdc_name2id, overrides, minted, next_id_box,
            args.fdc_id_offset, args.category_id, args.data_type,
        )
        if not foods.empty and not long.empty:
            all_foods.append(foods)
            all_long.append(long)
        all_unmapped.extend([(sh, u) for u in unmapped])
        all_minted.extend(minted_in_sheet)
        print(f"  [{sh}] foods={len(foods):,} measurements={len(long):,} "
              f"nutrients={int(long['nutrient_id'].nunique()) if not long.empty else 0}")

    if not all_foods:
        raise SystemExit("No foods produced.")
    foods_df = pd.concat(all_foods, ignore_index=True).drop_duplicates(subset=["fdc_id"])
    long_df = pd.concat(all_long, ignore_index=True)

    food_out = f"{args.out_prefix}_food_injection.parquet"
    foods_df.to_parquet(food_out, index=False)

    out_bucket = Path(f"{args.out_prefix}_food_nutrient_bucketed")
    if out_bucket.exists():
        shutil.rmtree(out_bucket)
    out_bucket.mkdir(exist_ok=True)
    long_df["bucket"] = (long_df["nutrient_id"] % 256).astype(int)
    for b, g in long_df.groupby("bucket", sort=True):
        d = out_bucket / f"bucket={int(b)}"
        d.mkdir(exist_ok=True)
        g[["fdc_id", "nutrient_id", "amount"]] \
            .sort_values(["fdc_id", "nutrient_id"]) \
            .to_parquet(d / "biofoodcomp_data.parquet", index=False)

    if all_unmapped:
        pd.DataFrame(all_unmapped, columns=["sheet", "unmapped_human_name"]) \
            .drop_duplicates() \
            .to_csv(f"{args.out_prefix}_unmapped_columns.tsv", sep="\t", index=False)
    if all_minted:
        pd.DataFrame(all_minted).to_csv(f"{args.out_prefix}_extra_nutrient_map.tsv",
                                          sep="\t", index=False)
    print()
    print(f"[OK] foods       : {len(foods_df)}")
    print(f"[OK] measurements: {len(long_df)}")
    print(f"[OK] nutrient ids: {long_df['nutrient_id'].nunique()}")
    print(f"[OK] minted ids  : {len(all_minted)} "
          f"({EXTRA_ID_START}-{next_id_box[0] - 1})")


if __name__ == "__main__":
    main()
