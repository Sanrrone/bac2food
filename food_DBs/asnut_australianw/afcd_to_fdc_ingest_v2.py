#!/usr/bin/env python3
"""afcd_to_fdc_ingest_v2.py — Phase 8 re-parse of AFCD Release 3.

v1 correctly extracts 1588 foods from the "All solids & liquids per 100 g"
sheet (the "Liquids only per 100 mL" sheet is a 213-row subset, all keys
already present in the larger sheet — verified). v1's defect is purely
nutrient-column coverage: ~234 of the ~270 raw columns get dropped because
their normalized name doesn't match anything in /data/bac2food/nutrient.csv,
including individual carotenoid forms (α-/β-/γ-/δ-carotene, lycopene cis/trans),
tocopherol/tocotrienol forms (α/β/γ/δ), 25-OH-D2 and 25-OH-D3 metabolites,
biotin, aluminum, and many trace minerals.

v2 changes:
  1. Mint nutrient_ids in the 230001+ block (separate from McCance/PhyFood/Bio
     blocks) for any unmapped column whose values are mostly-numeric.
  2. Drop the `amount > min_amount` filter (keeps zero as a legitimate value).

Imports v1's helper functions where possible.
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import os

THIS = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS))
sys.path.insert(0, str(THIS.parent))          # food_DBs/, for _common

from _common.units import (   # noqa: E402
    fdc_units, conversion_factor, report as unit_report,
)
from _common.minted import load_registry   # noqa: E402

NUTRIENT_CSV = os.environ.get("BAC2FOOD_NUTRIENT_CSV", "/data/bac2food/nutrient.csv")

# Reuse v1 helpers
from afcd_to_fdc_ingest import (   # noqa: E402
    norm, detect_header_row, category_from_name, build_fdc_nutrient_name_map,
    build_manual_overrides, parse_col_base_and_unit, to_numeric_series,
)

# fdc_id allocation: fdc_blocks.py
sys.path.insert(0, str(THIS.parent))
import fdc_blocks  # noqa: E402

EXTRA_ID_START = 230_001

# norm() folds an underscore to a space, so the literal "fdc_id" never matched
# and AFCD's own key column was minted as a nutrient (id 230104, 125 rows).
META_LIKE = {"public food key", "food name", "classification", "fdc_id", "fdc id",
             "food_category", "derivation", "biological name",
             "common name", "n", "samples"}


def is_mostly_numeric(s: pd.Series, min_frac: float = 0.4) -> bool:
    vals = s.dropna()
    if len(vals) == 0:
        return False
    nums = to_numeric_series(vals.head(200))
    n_ok = nums.notna().sum()
    return n_ok / max(1, min(200, len(vals))) >= min_frac


PER_GRAM_N_RE = re.compile(r"/\s*g\s*n\b", re.I)

# Source unit per column, filled in by map_columns_with_minting.
col2unit: Dict[str, str] = {}


def map_columns_with_minting(
    df: pd.DataFrame,
    fdc_name2id: Dict[str, int],
    overrides: Dict[str, int],
    minted: Dict[str, int],
    next_id_box: list[int],
) -> tuple[Dict[str, int], List[Dict], List[str]]:
    """Return (col -> nutrient_id), minted records, list of unmapped column names."""
    col2nid: Dict[str, int] = {}
    new_minted: List[Dict] = []
    unmapped: List[str] = []
    for c in df.columns:
        n = norm(c)
        if n in META_LIKE:
            continue
        if not is_mostly_numeric(df[c]):
            continue
        base, col_unit = parse_col_base_and_unit(c)
        base_n = norm(base)
        base_n2 = base_n.replace("dietary fibre", "fiber").replace("fibre", "fiber")
        # AFCD reports each amino acid twice: once per 100 g and once "(mg/gN)",
        # milligrams per gram of NITROGEN. The second is a different
        # measurement, not a different unit - it cannot be scaled onto the
        # per-100 g basis - so it is minted rather than mapped onto the FDC id,
        # where it would have sat beside the real value and won the MAX-union.
        if PER_GRAM_N_RE.search(str(col_unit or "")):
            nid = None
        else:
            nid = fdc_name2id.get(base_n) or fdc_name2id.get(base_n2) \
                  or overrides.get(base_n) or overrides.get(base_n2)
        if nid is None:
            stable = (base_n2 or base_n)
            if PER_GRAM_N_RE.search(str(col_unit or "")):
                stable += " per g n"
            if stable in minted:
                nid = minted[stable]
            else:
                nid = next_id_box[0]
                minted[stable] = nid
                next_id_box[0] += 1
                new_minted.append({
                    "nutrient_id": nid,
                    "source_db": "afcd",
                    "source_column_raw": str(c).replace("\n", " "),
                    "source_column_norm": stable,
                    "note": "minted (no FDC nutrient.csv match)",
                })
            unmapped.append(str(c).replace("\n", " "))
        col2nid[c] = int(nid)
        col2unit[c] = col_unit
    return col2nid, new_minted, unmapped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--excel", default="AFCD Release 3 - Nutrient profiles.xlsx")
    ap.add_argument("--sheet", default="All solids & liquids per 100 g")
    ap.add_argument("--fdc_nutrient_csv", default="/data/bac2food/nutrient.csv")
    ap.add_argument("--out_prefix", default="afcd")
    ap.add_argument("--category_id", default="7777")
    ap.add_argument("--data_type", default="foundation_food")
    args = ap.parse_args()

    import os; os.chdir(THIS)
    excel = args.excel
    sheet = args.sheet
    header_row = detect_header_row(excel, sheet, key="Public Food Key")
    df = pd.read_excel(excel, sheet_name=sheet, header=header_row).dropna(how="all")
    # AFCD sometimes has duplicate column headers (e.g. Energy with/without
    # dietary fibre with whitespace variants that collide after stripping).
    # Disambiguate by appending .1, .2, ... to duplicates so each column can
    # be addressed individually as a Series.
    seen: dict[str, int] = {}
    new_cols = []
    for c in df.columns:
        k = str(c)
        if k in seen:
            seen[k] += 1
            new_cols.append(f"{k}__dup{seen[k]}")
        else:
            seen[k] = 0
            new_cols.append(k)
    df.columns = new_cols
    print(f"[AFCD v2] {sheet}: {len(df)} rows × {len(df.columns)} cols")

    for r in ("Public Food Key", "Food Name"):
        if r not in df.columns:
            sys.exit(f"missing column: {r}")

    keys = df["Public Food Key"].astype(str).fillna("").tolist()
    uniq_keys = sorted({k.strip() for k in keys if k.strip()})
    key2fdc = dict(zip(uniq_keys, fdc_blocks.assign("afcd", uniq_keys)))
    df["fdc_id"] = df["Public Food Key"].astype(str) \
                     .map(lambda x: key2fdc.get(str(x).strip(), np.nan))
    df = df.dropna(subset=["fdc_id"]).copy()
    df["fdc_id"] = df["fdc_id"].astype(int)
    df["food_category"] = df["Food Name"].astype(str).map(category_from_name)

    fdc_name2id = build_fdc_nutrient_name_map(args.fdc_nutrient_csv)
    overrides = build_manual_overrides(fdc_name2id)
    # Seed the mint from the registry this ingester wrote last time, so a label
    # keeps the id it was first given. Positional minting renumbered the whole
    # block whenever a column was added or dropped, and the enzyme digest
    # references those ids BY NUMBER.
    _extra_map = f"{args.out_prefix}_extra_nutrient_map.tsv"
    minted: Dict[str, int] = load_registry(_extra_map, "source_column_norm")
    if minted:
        print(f"  -> reusing {len(minted)} minted ids from {_extra_map}")
    next_id_box = [max([*minted.values(), EXTRA_ID_START - 1]) + 1]

    col2nid, new_minted, unmapped = map_columns_with_minting(
        df, fdc_name2id, overrides, minted, next_id_box,
    )
    mapped_cols = sorted(col2nid.keys())
    n_known = sum(1 for c, nid in col2nid.items() if nid < EXTRA_ID_START)
    n_minted = sum(1 for c, nid in col2nid.items() if nid >= EXTRA_ID_START)
    print(f"  -> {len(mapped_cols)} nutrient columns: "
          f"{n_known} known FDC, {n_minted} minted ({EXTRA_ID_START}-{next_id_box[0] - 1})")

    wide = df[["fdc_id", "Food Name", "food_category"] + mapped_cols].copy()
    # Defensive: drop duplicate column labels just in case
    wide = wide.loc[:, ~wide.columns.duplicated()]
    # AFCD names the unit in the column header and the amount is then published
    # under FDC's unit for whatever id the column maps to. Where the two differ
    # the value is off by the ratio: the amino acids are mg here and G in FDC,
    # so 7,518 rows were deposited a thousand times too high.
    _fdcu = fdc_units(NUTRIENT_CSV)
    _applied = []
    for c in mapped_cols:
        dst = _fdcu.get(col2nid[c])
        if dst is None:
            continue                      # a minted id has no canonical unit
        f = conversion_factor(col2unit.get(c), dst)
        _applied.append((str(c).replace("\n", " "), str(col2unit.get(c)), str(dst), f))
        if f != 1.0 and c in wide.columns:
            wide[c] = to_numeric_series(wide[c]) * f
    unit_report(_applied)
    for c in mapped_cols:
        if c not in wide.columns:
            continue
        col = wide[c]
        if isinstance(col, pd.DataFrame):
            col = col.iloc[:, 0]   # collapse accidental dup
        wide[c] = to_numeric_series(col)
    long = wide.melt(id_vars=["fdc_id", "Food Name", "food_category"],
                     value_vars=mapped_cols,
                     var_name="afcd_col", value_name="amount").dropna(subset=["amount"])
    long["nutrient_id"] = long["afcd_col"].map(col2nid).astype(int)
    long = long[["fdc_id", "nutrient_id", "amount"]].sort_values(["fdc_id", "nutrient_id"])
    print(f"  -> {len(long)} measurements; {long['nutrient_id'].nunique()} distinct nutrient_ids")

    # Write outputs
    food_out = f"{args.out_prefix}_food_injection.parquet"
    foods = df[["fdc_id", "Food Name", "food_category"]].drop_duplicates().copy()
    foods["data_type"] = args.data_type
    foods["description"] = foods["Food Name"].astype(str) + " [AFCD]"
    foods["food_category_id"] = str(args.category_id)
    foods = foods[["fdc_id", "data_type", "description",
                   "food_category_id", "food_category"]]
    foods.to_parquet(food_out, index=False)

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
            .to_parquet(d / "afcd_data.parquet", index=False)

    pd.DataFrame({"Public Food Key": list(key2fdc),
                  "fdc_id": list(key2fdc.values())}) \
        .to_csv(f"{args.out_prefix}_id_map.tsv", sep="\t", index=False)

    if new_minted:
        pd.DataFrame(new_minted) \
            .to_csv(f"{args.out_prefix}_extra_nutrient_map.tsv",
                    sep="\t", index=False)
    if unmapped:
        pd.DataFrame({"unmapped_column": sorted(set(unmapped))}) \
            .to_csv(f"{args.out_prefix}_unmapped_columns.tsv",
                    sep="\t", index=False)

    print()
    print(f"[OK] foods       : {len(foods)}")
    print(f"[OK] measurements: {len(long)}")
    print(f"[OK] nutrient ids: {long['nutrient_id'].nunique()}")
    print(f"[OK] minted ids  : {len(new_minted)}")


if __name__ == "__main__":
    main()
