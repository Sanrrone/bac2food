#!/usr/bin/env python3
"""
export_matrices.py — Generates raw baseline matrices for microbiome differential feeding analysis.
"""
import pandas as pd
import pyarrow.dataset as ds
from pathlib import Path
import numpy as np

# --- PATHS (Adjust if necessary) ---
PATH_NUTRIENT_TO_EC = "../0_building/3_nutrient_to_ec.tsv"
PATH_NUTRIENT_CSV   = "/data/bac2food/nutrient.csv"
PATH_FOOD_CSV       = "/data/bac2food/food.parquet"
PATH_BUCKETED_DIR   = "/data/bac2food/food_nutrient_bucketed"
OUT_DIR             = Path("./boss_matrices")
# -----------------------------------

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("[*] Loading dictionaries...")
    
    # Load dictionaries for names
    nut_df = pd.read_csv(PATH_NUTRIENT_CSV, low_memory=False)
    nut_map = dict(zip(nut_df['id'].astype(str), nut_df['name'].astype(str)))
    
    food_df = pd.read_parquet(PATH_FOOD_CSV) if str(PATH_FOOD_CSV).lower().endswith(('.parquet', '.pq')) else pd.read_csv(PATH_FOOD_CSV, low_memory=False)
    food_map = dict(zip(food_df['fdc_id'].astype(str), food_df['description'].astype(str)))

    # =====================================================================
    # 1. GENERATE NUTRIENT-ENZYME DB (Binary Matrix)
    # =====================================================================
    print("[*] Generating Nutrient-Enzyme Matrix...")
    n2ec = pd.read_csv(PATH_NUTRIENT_TO_EC, sep="\t")
    # =====================================================================
    # 1. GENERATE NUTRIENT-ENZYME DB (Binary Matrix)
    # =====================================================================
    print("[*] Generating Nutrient-Enzyme Matrix...")
    n2ec = pd.read_csv(PATH_NUTRIENT_TO_EC, sep="\t")
    
    # THE CLEANUP FIX: 
    # Drop any enzyme mapping if the nutrient ID is completely missing from nutrient.csv
    valid_dict_ids = set(nut_map.keys())
    n2ec = n2ec[n2ec['nutrient_id'].astype(str).isin(valid_dict_ids)].copy()
    
    # Map nutrient names (We can remove the fillna("Unknown") now because they are all valid!)
    n2ec['nutrient_name'] = n2ec['nutrient_id'].astype(str).map(nut_map)
    n2ec['nutrient_key'] = n2ec['nutrient_id'].astype(str) + "_" + n2ec['nutrient_name']
    
    # Pivot into a 0/1 Matrix (Rows: Nutrients, Cols: EC Numbers)
    n2ec['value'] = 1
    enzyme_matrix = n2ec.pivot_table(
        index='nutrient_key', 
        columns='ec_number', 
        values='value', 
        fill_value=0
    ).reset_index()
    
    out_enz = OUT_DIR / "Nutrient_Enzyme_DB.csv"
    enzyme_matrix.to_csv(out_enz, index=False)
    print(f"  -> Saved to {out_enz} (Shape: {enzyme_matrix.shape})")
    
    # Get the list of ONLY nutrients that bacteria can actually use
    valid_nutrients = set(n2ec['nutrient_id'].astype(int).unique())

    # =====================================================================
    # 2. GENERATE FOOD-NUTRIENT DB (Restricted to Bacterial Nutrients)
    # =====================================================================
    print(f"[*] Generating Food-Nutrient Matrix (Restricted to {len(valid_nutrients)} bacterial nutrients)...")
    
    # Scan the bucketed parquet files, filtering ONLY for valid bacterial nutrients
    scn = ds.dataset(PATH_BUCKETED_DIR, format="parquet", partitioning="hive").scanner(
        columns=["fdc_id", "nutrient_id", "amount"],
        filter=ds.field("nutrient_id").isin(valid_nutrients)
    )
    
    food_nut_data = []
    for batch in scn.to_batches():
        df = batch.to_pandas()
        # Drop zero or negative amounts
        df = df[df['amount'] > 0]
        food_nut_data.append(df)
        
    full_fn_df = pd.concat(food_nut_data, ignore_index=True)
    
    # Map Names and create the exact same bulletproof ID_Name string
    full_fn_df['food_name'] = full_fn_df['fdc_id'].astype(str).map(food_map)
    full_fn_df['nutrient_name'] = full_fn_df['nutrient_id'].astype(str).map(nut_map).fillna("Unknown")
    full_fn_df['nutrient_key'] = full_fn_df['nutrient_id'].astype(str) + "_" + full_fn_df['nutrient_name']
    
    # Drop rows where food name wasn't found (e.g., branded foods you dropped earlier)
    full_fn_df = full_fn_df.dropna(subset=['food_name'])
    
    # Create the pivot table (Rows: Foods, Cols: Nutrients, Values: Amount)
    # Create the pivot table (Rows: Foods, Cols: Nutrients, Values: Amount)
    print("[*] Pivoting Food-Nutrient table (This may take a moment)...")
    food_matrix = full_fn_df.pivot_table(
        index=['fdc_id', 'food_name'],
        columns='nutrient_key',
        values='amount',
        aggfunc='max',  # Take max if duplicates exist
        fill_value=0.0
    ).reset_index()
    
    # =================================================================
    # THE ALIGNMENT FIX (Memory Optimized)
    # Force the Food matrix to have exactly the same nutrients as the Enzyme matrix
    # =================================================================
    print("[*] Aligning dimensions for strict matrix compatibility (Memory Optimized)...")
    
    # Get the exact list of 434 nutrients from the enzyme matrix
    required_nutrients = list(enzyme_matrix['nutrient_key'])
    
    # Define the exact final column structure we want
    final_columns = ['fdc_id', 'food_name'] + required_nutrients
    
    # reindex() adds missing columns and sorts them in a single, highly optimized C-level operation
    # By specifying fill_value=0.0, it automatically pads the missing 265 nutrients with zeros.
    food_matrix = food_matrix.reindex(columns=final_columns, fill_value=0.0)
    # =================================================================
    
    out_food = OUT_DIR / "Food_Nutrient_DB.csv"
    food_matrix.to_csv(out_food, index=False)
    print(f"  -> Saved to {out_food} (Shape: {food_matrix.shape})")
    
    print("\n[*] DONE!")

if __name__ == "__main__":
    main()