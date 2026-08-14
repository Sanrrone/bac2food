#!/usr/bin/env python3
"""ingest_phenol_explorer_v2.py — Phase 8 re-parse of Phenol-Explorer 3.6.

The composition-data.tsv has 7486 rows × 17 cols documenting (food, compound,
method, mean ± min/max) measurements for **508 distinct polyphenols** across
458 foods.

The v1 ingester collapsed all of those compounds into 7 super-group nutrient
IDs (Flavonoids/Flavonols/Flavanols/Anthocyanins/Isoflavonoids/Phenolic acids/
Lignans) and summed means within each group. That sum-then-discard is the
single largest information loss in the entire bac2food pipeline — every
distinct bacterial substrate (quercetin / kaempferol / caffeic acid / etc.)
became a single number per food.

v2:
  - One nutrient_id per distinct compound (508 minted in the 240001+ block).
  - Multiple methods for the same (food, compound) collapse via MEAN of means
    (matches Phase 4 MAX→MEAN decision for cultivar aggregation).
  - All units are already normalized to mg/100g (or mg/100mL for liquids) —
    we keep them as-is; the EC graph treats mg the same as mg.
  - Compound metadata persisted in extra_nutrient_map.tsv so 0_building/
    2_nutri2chebi_from_obo.py can wire ChEBI links by IUPAC name later.
"""
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

# fdc_id allocation is centralised in food_DBs/fdc_blocks.py. Never write a literal offset:
# ids are accessions looked up in fdc_id_map.tsv, so a food keeps its id across releases.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
import fdc_blocks


EXTRA_ID_START = 240_001
# RETIRED: block bases live in food_DBs/fdc_blocks.py. This constant is kept out of
# the file entirely so it cannot drift from the real allocation.


CAT_MAP = {
    "Fruits and fruit products": "Fruits and Fruit Juices",
    "Vegetables": "Vegetables and Vegetable Products",
    "Seasonings": "Spices and Herbs",
    "Seeds": "Nut and Seed Products",
    "Non-alcoholic beverages": "Beverages",
    "Alcoholic beverages": "Alcoholic Beverages",
    "Cereals and cereal products": "Cereal Grains and Pasta",
    "Oils": "Fats and Oils",
    "Coffee and cocoa": "Beverages",
    "Cocoa products": "Sweets",
    "Tea": "Beverages",
    "Other foods": "Meals, Entrees, and Side Dishes",
}


def normalize(s: str) -> str:
    return re.sub(r"\s+", " ", str(s)).strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", default="composition-data.tsv")
    ap.add_argument("--out_prefix", default="pe")
    args = ap.parse_args()

    import os; os.chdir(Path(__file__).parent)

    print("==================================================")
    print("Phenol-Explorer v2 ingester")
    print("==================================================")

    df = pd.read_csv(args.tsv, sep="\t", dtype=str)
    df["mean"] = pd.to_numeric(df["mean"], errors="coerce")
    print(f"[1/4] {len(df)} raw rows; {df['food'].nunique()} foods; "
          f"{df['compound'].nunique()} compounds; "
          f"{df['experimental_method_group'].nunique()} methods")

    # --- foods table -------------------------------------------------------
    foods_index = df[["food", "food_group", "food_sub_group"]].drop_duplicates() \
                                                              .reset_index(drop=True)
    # Keyed on the food NAME, which is what Phenol-Explorer publishes as its identifier;
    # range() gave no stable key and renumbered every food whenever the set changed.
    foods_index["fdc_id"] = fdc_blocks.assign("phenol_explorer", foods_index["food"])
    foods_index["food_category"] = foods_index["food_group"].map(CAT_MAP) \
                                                            .fillna("Other")
    food_to_fdc = dict(zip(foods_index["food"], foods_index["fdc_id"]))
    df["fdc_id"] = df["food"].map(food_to_fdc)

    # --- compound table (one nutrient_id per compound) --------------------
    compounds = sorted(df["compound"].dropna().unique())
    print(f"[2/4] minting {len(compounds)} compound nutrient_ids "
          f"(range {EXTRA_ID_START}-{EXTRA_ID_START + len(compounds) - 1})")
    compound_to_id = {c: EXTRA_ID_START + i for i, c in enumerate(compounds)}
    df["nutrient_id"] = df["compound"].map(compound_to_id)

    # Build the extra_nutrient_map with compound metadata for downstream
    # ChEBI linking in 0_building/2_nutri2chebi_from_obo.py
    extra_records = []
    cmeta = df[["compound", "compound_group", "compound_sub_group", "units"]] \
              .drop_duplicates(subset=["compound"]).set_index("compound")
    for c in compounds:
        unit = cmeta.loc[c, "units"] if c in cmeta.index else ""
        extra_records.append({
            "nutrient_id": compound_to_id[c],
            "source_db": "phenol_explorer",
            "compound": c,
            "compound_group": cmeta.loc[c, "compound_group"] if c in cmeta.index else "",
            "compound_sub_group": cmeta.loc[c, "compound_sub_group"]
                                  if c in cmeta.index else "",
            "units": unit,
            "note": "minted (Phenol-Explorer compound, one ID per compound)",
        })

    # --- collapse multiple methods per (food, compound) via mean of means -
    print("[3/4] collapsing methods (mean of means) ...")
    long = df.dropna(subset=["fdc_id", "nutrient_id", "mean"]) \
             .groupby(["fdc_id", "nutrient_id"], as_index=False)["mean"].mean() \
             .rename(columns={"mean": "amount"})
    long["fdc_id"] = long["fdc_id"].astype(int)
    long["nutrient_id"] = long["nutrient_id"].astype(int)
    print(f"  -> {len(long)} measurements")
    print(f"  -> {long['nutrient_id'].nunique()} distinct nutrient ids")

    # --- foods parquet ----------------------------------------------------
    print("[4/4] writing outputs ...")
    foods = pd.DataFrame({
        "fdc_id": foods_index["fdc_id"].astype(int),
        "data_type": "foundation_food",
        "description": foods_index["food"] + " [Phenol-Explorer]",
        "food_category_id": "9999",
        "food_category": foods_index["food_category"],
    })
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
            .to_parquet(d / "phenol_explorer_data.parquet", index=False)

    pd.DataFrame(extra_records) \
        .to_csv(f"{args.out_prefix}_extra_nutrient_map.tsv", sep="\t", index=False)

    print()
    print(f"[OK] foods       : {len(foods)}")
    print(f"[OK] measurements: {len(long)}")
    print(f"[OK] compound IDs: {long['nutrient_id'].nunique()} "
          f"(was 7 in v1)")


if __name__ == "__main__":
    main()
