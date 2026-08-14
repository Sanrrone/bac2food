#!/usr/bin/env python3
"""merge_phase8_v2.py — Integrate every Phase 8 v2 ingester output into
/data/bac2food/food.parquet + /data/bac2food/food_nutrient_bucketed/.

For each affected source (McCance / PhyFoodComp / BioFoodComp / AFCD /
Phenol-Explorer / STFCJ):
  1. Strip the source's existing rows from /data/bac2food/food.parquet
     (identified by the trailing source-tag substring in `description`).
  2. Append the v2 food_injection.parquet rows.
  3. Delete the source's per-bucket parquet files (the named-by-source
     `<src>_data.parquet`) from /data/bac2food/food_nutrient_bucketed/, then
     copy the v2 bucketed parquets in.
  4. Append the source's extra_nutrient_map.tsv rows into nutrient.csv.

Finally, wipes the modeled index cache at /data/bac2food/index_modeled/.

Safe to re-run; creates a single timestamped backup of food.parquet at start.
"""
from __future__ import annotations

import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq

GLOBAL_FOOD = Path("/data/bac2food/food.parquet")
GLOBAL_BUCKETS = Path("/data/bac2food/food_nutrient_bucketed")
NUTRIENT_CSV = Path("/data/bac2food/nutrient.csv")
INDEX_DIR = Path("/data/bac2food/index_modeled")
FOOD_DBS = Path("/home/svalenzuela/Desktop/bac2food/food_DBs")


@dataclass
class V2Source:
    short: str
    folder: str
    tag_lower: str
    food_injection: str       # e.g. mccance_food_injection.parquet
    bucket_dir: str           # e.g. mccance_food_nutrient_bucketed
    data_filename: str        # e.g. mccance_data.parquet — the file under bucket=X/
    extra_map: str            # e.g. mccance_extra_nutrient_map.tsv
    expected_min_foods: int   # safety check — refuse to swap if v2 has fewer foods


SOURCES = [
    V2Source("mccance", "McCance_Widdowsons_uk", "[mccance]",
             "mccance_food_injection.parquet",
             "mccance_food_nutrient_bucketed",
             "mccance_data.parquet",
             "mccance_extra_nutrient_map.tsv",
             expected_min_foods=2800),
    V2Source("phyfoodcomp", "FAO_onu", "[phyfoodcomp]",
             "phyfoodcomp_food_injection.parquet",
             "phyfoodcomp_food_nutrient_bucketed",
             "phyfoodcomp_data.parquet",
             "phyfoodcomp_extra_nutrient_map.tsv",
             expected_min_foods=3000),
    V2Source("biofoodcomp", "FAO_onu", "[biofoodcomp]",
             "biofoodcomp_food_injection.parquet",
             "biofoodcomp_food_nutrient_bucketed",
             "biofoodcomp_data.parquet",
             "biofoodcomp_extra_nutrient_map.tsv",
             expected_min_foods=9500),
    V2Source("afcd", "asnut_australianw", "[afcd]",
             "afcd_food_injection.parquet",
             "afcd_food_nutrient_bucketed",
             "afcd_data.parquet",
             "afcd_extra_nutrient_map.tsv",
             expected_min_foods=1500),
    V2Source("phenol_explorer", "phenol_explorer_france", "[phenol-explorer]",
             "pe_food_injection.parquet",
             "pe_food_nutrient_bucketed",
             "phenol_explorer_data.parquet",
             "pe_extra_nutrient_map.tsv",
             expected_min_foods=400),
    V2Source("stfcj", "stfcj_japan", "[stfcj]",
             "japan_food_injection.parquet",
             "japan_food_nutrient_bucketed",
             "japan_data.parquet",
             "japan_extra_nutrient_map.tsv",
             expected_min_foods=2000),
]


def check_v2_outputs():
    missing = []
    for s in SOURCES:
        folder = FOOD_DBS / s.folder
        for f in (s.food_injection, s.bucket_dir):
            p = folder / f
            if not p.exists():
                missing.append(str(p))
    if missing:
        print("[!] Missing v2 outputs:")
        for m in missing:
            print(f"    {m}")
        sys.exit(1)


def load_v2_food_count(s: V2Source) -> int:
    p = FOOD_DBS / s.folder / s.food_injection
    return pq.read_table(p, columns=["fdc_id"]).num_rows


def strip_source_rows(food_df: pd.DataFrame, tag_lower: str) -> tuple[pd.DataFrame, int]:
    """Drop rows whose lower-cased description contains the tag substring."""
    mask = food_df["description"].astype(str).str.lower() \
                .str.contains(tag_lower, regex=False, na=False)
    n_dropped = int(mask.sum())
    return food_df[~mask].copy(), n_dropped


def append_extra_nutrients() -> int:
    """Append minted nutrient IDs from each source's extra_nutrient_map.tsv
    into /data/bac2food/nutrient.csv. Returns the count appended."""
    nut = pd.read_csv(NUTRIENT_CSV)
    existing_ids = set(nut["id"].astype(int))
    rows_to_add: list[dict] = []
    for s in SOURCES:
        p = FOOD_DBS / s.folder / s.extra_map
        if not p.exists():
            continue
        em = pd.read_csv(p, sep="\t")
        for _, r in em.iterrows():
            nid = int(r["nutrient_id"])
            if nid in existing_ids:
                continue
            existing_ids.add(nid)
            # Synthesize a nutrient.csv row from the available fields
            name = (r.get("source_column_raw") or r.get("source_column_norm")
                    or r.get("compound") or f"nutrient_{nid}")
            unit = r.get("units", "") if "units" in r else ""
            rows_to_add.append({
                "id": nid, "name": str(name)[:300], "unit_name": str(unit)[:30],
                "nutrient_nbr": nid, "rank": 999999,
            })
    if not rows_to_add:
        return 0
    # Align columns with existing nutrient.csv schema
    add_df = pd.DataFrame(rows_to_add)
    for c in nut.columns:
        if c not in add_df.columns:
            add_df[c] = ""
    add_df = add_df[nut.columns]
    out = pd.concat([nut, add_df], ignore_index=True)
    out.to_csv(NUTRIENT_CSV, index=False)
    return len(rows_to_add)


def swap_bucketed_data(s: V2Source) -> tuple[int, int]:
    """For one source: delete per-source files from each global bucket dir,
    then copy v2 bucketed files in. Returns (n_files_deleted, n_files_copied)."""
    src_root = FOOD_DBS / s.folder / s.bucket_dir
    n_deleted = 0
    n_copied = 0
    # Delete existing per-source files anywhere under GLOBAL_BUCKETS
    for hit in GLOBAL_BUCKETS.rglob(s.data_filename):
        hit.unlink()
        n_deleted += 1
    # Copy v2 files in
    for src_file in src_root.rglob(s.data_filename):
        rel = src_file.parent.name  # "bucket=42"
        dst_dir = GLOBAL_BUCKETS / rel
        dst_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src_file, dst_dir / s.data_filename)
        n_copied += 1
    return n_deleted, n_copied


def main():
    print("==================================================")
    print("Phase 8 v2 merge")
    print("==================================================")

    check_v2_outputs()

    # Safety check
    for s in SOURCES:
        n = load_v2_food_count(s)
        if n < s.expected_min_foods:
            sys.exit(f"[abort] {s.short}: v2 has only {n} foods, "
                     f"expected ≥ {s.expected_min_foods}. "
                     "Did the v2 ingester complete?")

    # Backup
    backup = Path(f"/data/bac2food/food_backup_phase8.parquet")
    if not backup.exists():
        print(f"[1/5] backup: cp {GLOBAL_FOOD} -> {backup}")
        shutil.copy2(GLOBAL_FOOD, backup)
    else:
        print(f"[1/5] backup exists at {backup} (not overwriting)")

    # Strip + append in food.parquet
    print("[2/5] rewriting food.parquet ...")
    food = pd.read_parquet(GLOBAL_FOOD)
    print(f"      pre:  {len(food)} rows")
    total_stripped = 0
    total_added = 0
    for s in SOURCES:
        food, n_drop = strip_source_rows(food, s.tag_lower)
        v2_food = pd.read_parquet(FOOD_DBS / s.folder / s.food_injection)
        food = pd.concat([food, v2_food], ignore_index=True)
        total_stripped += n_drop
        total_added += len(v2_food)
        print(f"      {s.short:18s}  strip {n_drop:6d}  add {len(v2_food):6d}")
    # Drop duplicate fdc_ids (keep v2 = last)
    food = food.drop_duplicates(subset=["fdc_id"], keep="last")
    food.to_parquet(GLOBAL_FOOD, index=False)
    print(f"      post: {len(food)} rows  (stripped {total_stripped}, added {total_added})")

    # Swap bucketed nutrient data
    print("[3/5] swapping bucketed nutrient data ...")
    for s in SOURCES:
        n_del, n_cp = swap_bucketed_data(s)
        print(f"      {s.short:18s}  deleted {n_del} stale files, copied {n_cp} v2 files")

    # Append extra nutrient IDs
    print("[4/5] appending new nutrient IDs to nutrient.csv ...")
    n_added = append_extra_nutrients()
    print(f"      +{n_added} new nutrient IDs")

    # Wipe modeled index cache
    print("[5/5] wiping index_modeled cache ...")
    if INDEX_DIR.exists():
        for ext in ("*.parquet", "*.pkl"):
            for f in INDEX_DIR.glob(ext):
                f.unlink()
                print(f"      rm {f.name}")
    print()
    print("DONE.")


if __name__ == "__main__":
    main()
