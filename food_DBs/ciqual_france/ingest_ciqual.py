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
    print("CIQUAL (FRANCE) TO FDC INGESTION SCRIPT")
    print("==================================================")
    
    print("[1/5] Loading CIQUAL composition data from Excel...")
    # Read directly from the Excel file, targeting the specific data sheet
    # (Update the filename below if yours is slightly different, e.g., .xls instead of .xlsx)
    excel_file = "Table Ciqual 2025_ENG_2025_11_03.xlsx"
    ciqual_df = pd.read_excel(excel_file, sheet_name="food composition")
    
    print("[2/5] Cleaning CIQUAL strings and numeric formats...")
    # Clean the column names (remove newlines so they match our map perfectly)
    ciqual_df.columns = ciqual_df.columns.str.replace('\n', ' ')
    
    # Translating CIQUAL columns directly to USDA FDC Nutrient IDs
    target_mapping = {
        'Protein (g 100g)': 1003,
        'Fat (g 100g)': 1004,
        'Carbohydrate (g 100g)': 1005,
        'Starch (g 100g)': 1009,
        'Sucrose (g 100g)': 1010,
        'Glucose (g 100g)': 1011,
        'Fructose (g 100g)': 1012,
        'Lactose (g 100g)': 1013,
        'Maltose (g 100g)': 1014,
        'Galactose (g 100g)': 1015,
        'Polyols (g 100g)': 1062,
        'Fibres (g 100g)': 1079,
        'Calcium (mg 100g)': 1087,
        'Iron (mg 100g)': 1089,
        'Magnesium (mg 100g)': 1090,
        'Phosphorus (mg 100g)': 1091,
        'Potassium (mg 100g)': 1092,
        'Sodium (mg 100g)': 1093,
        'Zinc (mg 100g)': 1095,
        'Vitamin A activity, retinol equivalent (µg 100mg)': 1104,
        'Vitamin D (µg 100g)': 1110,
        'Vitamin C (mg 100g)': 1162,
        'Vitamin B12 (µg 100g)': 1174,
        'Vitamin B9 or total folates (µg 100g)': 1177,
        'Cholesterol (mg 100g)': 1253,
        'FA saturated (g 100g)': 1258,
        'FA mono (g 100g)': 1259,
        'FA poly (g 100g)': 1260,
        'Sugars (g 100g)': 2000
    }
    
    # Clean the numeric data for these target columns
    available_cols = [c for c in target_mapping.keys() if c in ciqual_df.columns]
    
    for c in available_cols:
        if ciqual_df[c].dtype == object:
            # CIQUAL uses "-" for missing, "< 0.1" for trace, and "traces" for trace
            ciqual_df[c] = ciqual_df[c].astype(str)
            ciqual_df[c] = ciqual_df[c].str.replace(r'^<.*', '0.0', regex=True)
            ciqual_df[c] = ciqual_df[c].str.replace('traces', '0.0', regex=False)
            ciqual_df[c] = ciqual_df[c].str.replace('-', 'NaN', regex=False) 
            ciqual_df[c] = ciqual_df[c].str.replace(',', '.', regex=False)
            ciqual_df[c] = pd.to_numeric(ciqual_df[c], errors='coerce')

    print("[3/5] Minting custom FDC IDs for CIQUAL foods...")
    # Offset CIQUAL IDs by 20,000,000 to prevent collisions with USDA, Phenol-Explorer, or Fineli
    ciqual_df['fdc_id'] = fdc_blocks.assign('ciqual', ciqual_df['alim_code'].astype(int))
    
    print("[4/5] Mapping CIQUAL categories and melting...")
    # Clean category string
    ciqual_df['alim_grp_nom_eng'] = ciqual_df['alim_grp_nom_eng'].str.replace('\n', ' ')
    cat_map = {
        'fruits, vegetables, legumes and nuts': 'Vegetables and Vegetable Products',
        'meat, egg and fish': 'Pork Products', 
        'starters and dishes': 'Meals, Entrees, and Side Dishes',
        'sugar and confectionery': 'Sweets',
        'milk and milk products': 'Dairy and Egg Products',
        'beverages': 'Beverages',
        'cereal products': 'Cereal Grains and Pasta',
        'fats and oils': 'Fats and Oils',
        'baby food': 'Baby Foods',
        'ice cream and sorbet': 'Sweets'
    }
    ciqual_df['food_category'] = ciqual_df['alim_grp_nom_eng'].map(cat_map).fillna('Meals, Entrees, and Side Dishes')
    
    # Melt the dataframe
    melted_df = ciqual_df.melt(id_vars=['fdc_id', 'alim_nom_eng', 'food_category'], value_vars=available_cols, var_name='ciqual_nutrient', value_name='amount')
    
    # Drop NAs and zeroes
    melted_df = melted_df.dropna(subset=['amount'])
    melted_df = melted_df[melted_df['amount'] > 0]
    
    # Map names to numeric FDC IDs
    melted_df['nutrient_id'] = melted_df['ciqual_nutrient'].map(target_mapping).astype(int)
    
    print("[5/5] Generating database files...")
    # Create the food dictionary
    unique_foods = melted_df[['fdc_id', 'alim_nom_eng', 'food_category']].drop_duplicates()
    new_foods = pd.DataFrame({
        'fdc_id': unique_foods['fdc_id'],
        'data_type': 'foundation_food',
        'description': unique_foods['alim_nom_eng'] + " [CIQUAL]",
        'food_category_id': '7777',
        'food_category': unique_foods['food_category']
    })
    
    new_foods.to_parquet("ciqual_food_injection.parquet", index=False)
    
    # Build bucketed Parquets
    melted_df['bucket'] = melted_df['nutrient_id'] % 256
    
    out_dir = Path("ciqual_food_nutrient_bucketed")
    out_dir.mkdir(exist_ok=True)
    
    for bucket_id, group in melted_df.groupby('bucket'):
        bucket_dir = out_dir / f"bucket={bucket_id}"
        bucket_dir.mkdir(exist_ok=True)
        group = group.sort_values(['fdc_id', 'nutrient_id'])
        group[['fdc_id', 'nutrient_id', 'amount']].to_parquet(bucket_dir / "ciqual_data.parquet", index=False)

    print("\n[SUCCESS] CIQUAL integration files generated!")
    print("1. 'ciqual_food_injection.parquet' -> Contains the food metadata.")
    print("2. 'ciqual_food_nutrient_bucketed/' -> Contains the partitioned FDC measurements.")

if __name__ == "__main__":
    main()