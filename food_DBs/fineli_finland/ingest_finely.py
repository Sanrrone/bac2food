#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq

# fdc_id allocation: fdc_blocks.py
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
import fdc_blocks


def main():
    print("==================================================")
    print("FINELI TO FDC INGESTION SCRIPT")
    print("==================================================")
    
    # 1. LOAD FINELI DATA
    print("[1/5] Loading Fineli composition data...")
    # Fineli uses semi-colons as separators
    fineli_df = pd.read_csv("resultset.csv", sep=";")
    
    print("[2/5] Cleaning European numeric formats and trace values...")
    # Fineli uses commas for decimals and '<0.1' for trace amounts. We must clean this for PyArrow.
    for c in fineli_df.columns[2:]: # Skip 'id' and 'name'
        if fineli_df[c].dtype == object:
            # Replace trace markers with 0.0, replace commas with dots
            fineli_df[c] = fineli_df[c].astype(str).str.replace(r'^<.*', '0.0', regex=True)
            fineli_df[c] = fineli_df[c].str.replace(',', '.')
            fineli_df[c] = pd.to_numeric(fineli_df[c], errors='coerce')
            
    # 2. CREATE FDC-COMPATIBLE FOOD IDS
    print("[3/5] Minting custom FDC IDs for Fineli foods...")
    # Offset Fineli IDs by 10,000,000 to prevent collisions with USDA or Phenol-Explorer
    fineli_df['fdc_id'] = fdc_blocks.assign('fineli', fineli_df['id'].astype(int))
    
    # 3. MAP FINELI COMPOUNDS TO FDC TARGETS
    print("[4/5] Mapping Fineli nutrients to FDC IDs...")
    # Translating Fineli columns directly to USDA FDC Nutrient IDs
    target_mapping = {
        'protein, total (g)': 1003,
        'fat, total (g)': 1004,
        'carbohydrate, available (g)': 1005,
        'starch, total (g)': 1009,
        'sucrose (g)': 1010,
        'glucose (g)': 1011,
        'fructose (g)': 1012,
        'lactose (g)': 1013,
        'maltose (g)': 1014,
        'galactose (g)': 1015,
        'sugar alcohols (g)': 1062,
        'fibre, total (g)': 1079,
        'calcium (mg)': 1087,
        'iron, total (mg)': 1089,
        'magnesium (mg)': 1090,
        'phosphorus (mg)': 1091,
        'potassium (mg)': 1092,
        'sodium (mg)': 1093,
        'zinc (mg)': 1095,
        'vitamin A  retinol activity equivalents (µg)': 1104,
        'folate, total (µg)': 1177,
        'vitamin D (µg)': 1110,
        'vitamin C (ascorbic acid) (mg)': 1162,
        'vitamin B-12 (cobalamin) (µg)': 1174,
        'cholesterol (GC) (mg)': 1253,
        'fatty acids, total saturated (g)': 1258,
        'fatty acids, total monounsaturated cis (g)': 1259,
        'fatty acids, total polyunsaturated (g)': 1260,
        'fatty acids, total trans (g)': 1257,
        'sugars, total (g)': 2000
    }
    
    # Melt the dataframe into long format (fdc_id, nutrient, amount)
    available_cols = [c for c in target_mapping.keys() if c in fineli_df.columns]
    melted_df = fineli_df.melt(id_vars=['fdc_id'], value_vars=available_cols, var_name='fineli_nutrient', value_name='amount')
    
    # Drop NAs and zero values to keep the index blazing fast
    melted_df = melted_df.dropna(subset=['amount'])
    melted_df = melted_df[melted_df['amount'] > 0]
    
    # Map the string names to their numeric FDC IDs
    melted_df['nutrient_id'] = melted_df['fineli_nutrient'].map(target_mapping).astype(int)
    
    # 4. APPEND TO FDC FOOD.PARQUET
    print("[5/5] Generating database files...")
    new_foods = pd.DataFrame({
        'fdc_id': fineli_df['fdc_id'],
        'data_type': 'foundation_food', 
        'description': fineli_df['name'] + " [Fineli]",
        'food_category_id': '8888', 
        'food_category': 'Fineli Database' # Distinct category avoids USDA processing penalties naturally
    })
    
    new_foods.to_parquet("fineli_food_injection.parquet", index=False)
    
    # 5. BUILD THE BUCKETED PARQUET FILES
    melted_df['bucket'] = melted_df['nutrient_id'] % 256
    
    out_dir = Path("fineli_food_nutrient_bucketed")
    out_dir.mkdir(exist_ok=True)
    
    for bucket_id, group in melted_df.groupby('bucket'):
        bucket_dir = out_dir / f"bucket={bucket_id}"
        bucket_dir.mkdir(exist_ok=True)
        # Sort for rapid Arrow scanning
        group = group.sort_values(['fdc_id', 'nutrient_id'])
        group[['fdc_id', 'nutrient_id', 'amount']].to_parquet(bucket_dir / "fineli_data.parquet", index=False)

    print("\n[SUCCESS] Fineli integration files generated!")
    print("1. 'fineli_food_injection.parquet' -> Contains the food metadata.")
    print("2. 'fineli_food_nutrient_bucketed/' -> Contains the partitioned FDC measurements.")

if __name__ == "__main__":
    main()
