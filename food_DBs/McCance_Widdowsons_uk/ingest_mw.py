#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq

def main():
    print("==================================================")
    print("McCANCE & WIDDOWSON (UK) TO FDC INGESTION SCRIPT")
    print("==================================================")
    
    excel_file = "McCance_Widdowsons_Composition_of_Foods_Integrated_Dataset_2021..xlsx"
    
    # The sheets we want to merge to build a complete food profile
    sheets = [
        "1.3 Proximates", "1.4 Inorganics", "1.5 Vitamins", 
        "1.8 (SFA per 100gFood)", "1.10 (MUFA per 100gFood)", 
        "1.12 (PUFA per 100gFood)", "1.13 Phytosterols", "1.14 Organic Acids"
    ]
    
    print("[1/5] Loading and merging Excel sheets...")
    merged_df = None
    
    for sheet in sheets:
        try:
            # Header=0 reads the actual nutrient names. We then drop rows 0 and 1 which are metadata.
            df = pd.read_excel(excel_file, sheet_name=sheet, header=0)
            df = df.iloc[2:].copy() 
            df.rename(columns={df.columns[0]: 'Food Code'}, inplace=True)
            
            # Drop redundant metadata columns if we are merging
            meta_to_drop = ['Food Name', 'Description', 'Group', 'Previous', 'Main data references', 'Footnote']
            if merged_df is not None:
                df = df.drop(columns=[c for c in meta_to_drop if c in df.columns], errors='ignore')
                merged_df = pd.merge(merged_df, df, on='Food Code', how='outer')
            else:
                merged_df = df
        except Exception as e:
            print(f"  -> Skipping {sheet} or encountered error: {e}")
            
    print(f"  -> Successfully merged {len(merged_df)} foods!")

    print("[2/5] Minting custom FDC IDs...")
    # McCance uses alphanumeric IDs like "13-145", so we mint numeric IDs starting at 40,000,000
    merged_df = merged_df.reset_index(drop=True)
    raise SystemExit(
        "This v1 ingest mints fdc_ids from a literal offset, which the block scheme "
        "no longer permits: ids are accessions held in fdc_id_map.tsv. Run ingest_mw_v2.py "
        "instead, which allocates through food_DBs/fdc_blocks.py."
    )
    
    print("[3/5] Cleaning British trace formats and mapping categories...")
    # Map McCance's Group column to FDC Categories
    cat_map = {
        'DG': 'Vegetables and Vegetable Products',
        'C': 'Cereal Grains and Pasta',
        'A': 'Cereal Grains and Pasta', # Breads/Flours
        'DR': 'Fruits and Fruit Juices',
        'N': 'Nut and Seed Products',
        'M': 'Dairy and Egg Products',
        'O': 'Fats and Oils',
        'H': 'Beef Products', # Meats
        'J': 'Finfish and Shellfish Products',
        'S': 'Sweets',
        'P': 'Beverages',
        'Q': 'Alcoholic Beverages',
        'DI': 'Spices and Herbs'
    }
    # We use the first letter of the Group code to roughly approximate FDC categories
    merged_df['food_category'] = merged_df['Group'].astype(str).str[0].map(cat_map).fillna('Meals, Entrees, and Side Dishes')

    print("[4/5] Mapping McCance nutrients to FDC IDs...")
    target_mapping = {
        'Protein (g)': 1003,
        'Fat (g)': 1004,
        'Carbohydrate (g)': 1005,
        'Starch (g)': 1009,
        'Total sugars (g)': 2000,
        'Glucose (g)': 1011,
        'Fructose (g)': 1012,
        'Lactose (g)': 1013,
        'Maltose (g)': 1014,
        'Sucrose (g)': 1010,
        'Galactose (g)': 1015,
        'Alcohol (g)': 1018,
        'AOAC fibre (g)': 1079,
        'Sodium (mg)': 1093,
        'Potassium (mg)': 1092,
        'Calcium (mg)': 1087,
        'Magnesium (mg)': 1090,
        'Phosphorus (mg)': 1091,
        'Iron (mg)': 1089,
        'Copper (mg)': 1098,
        'Zinc (mg)': 1095,
        'Chloride (mg)': 1088,
        'Manganese (mg)': 1101,
        'Selenium (µg)': 1103,
        'Iodine (µg)': 1100,
        'Retinol (µg)': 1105,
        'Carotene (µg)': 1106,
        'Vitamin D (µg)': 1110,
        'Vitamin E (mg)': 1109,
        'Vitamin K1 (µg)': 1185,
        'Thiamin (mg)': 1165,
        'Riboflavin (mg)': 1166,
        'Niacin (mg)': 1167,
        'Vitamin B6 (mg)': 1175,
        'Vitamin B12 (µg)': 1174,
        'Folate (µg)': 1177,
        'Pantothenate (mg)': 1170,
        'Biotin (µg)': 1176,
        'Vitamin C (mg)': 1162,
        'Cholesterol (mg)': 1253,
        'Total Phytosterols (mg)': 1283,
        'Beta-sitosterol (mg)': 1288,
        'Campesterol (mg)': 1287,
        'Stigmasterol (mg)': 1286,
        'Citric acid (g)': 1021,
        'Malic acid (g)': 1022
    }
    
    available_cols = [c for c in target_mapping.keys() if c in merged_df.columns]
    
    # McCance uses "Tr" for Trace and "N" for Not Present.
    for c in available_cols:
        if merged_df[c].dtype == object:
            merged_df[c] = merged_df[c].astype(str)
            merged_df[c] = merged_df[c].replace(['Tr', 'N', '-', 'nan'], ['0.0', '0.0', 'NaN', 'NaN'])
            merged_df[c] = pd.to_numeric(merged_df[c], errors='coerce')
            
    # Melt the dataframe into long format
    melted_df = merged_df.melt(id_vars=['fdc_id', 'Food Name', 'food_category'], 
                               value_vars=available_cols, 
                               var_name='mccance_nutrient', 
                               value_name='amount')
    
    # Drop NAs and zeroes
    melted_df = melted_df.dropna(subset=['amount'])
    melted_df = melted_df[melted_df['amount'] > 0]
    
    # Map names to numeric FDC IDs
    melted_df['nutrient_id'] = melted_df['mccance_nutrient'].map(target_mapping).astype(int)
    
    print("[5/5] Generating database files...")
    # Create the food dictionary
    unique_foods = melted_df[['fdc_id', 'Food Name', 'food_category']].drop_duplicates()
    new_foods = pd.DataFrame({
        'fdc_id': unique_foods['fdc_id'],
        'data_type': 'foundation_food',
        'description': unique_foods['Food Name'] + " [McCance]",
        'food_category_id': '5555',
        'food_category': unique_foods['food_category']
    })
    
    new_foods.to_parquet("mccance_food_injection.parquet", index=False)
    
    # Build bucketed Parquets
    melted_df['bucket'] = melted_df['nutrient_id'] % 256
    
    out_dir = Path("mccance_food_nutrient_bucketed")
    out_dir.mkdir(exist_ok=True)
    
    for bucket_id, group in melted_df.groupby('bucket'):
        bucket_dir = out_dir / f"bucket={bucket_id}"
        bucket_dir.mkdir(exist_ok=True)
        group = group.sort_values(['fdc_id', 'nutrient_id'])
        group[['fdc_id', 'nutrient_id', 'amount']].to_parquet(bucket_dir / "mccance_data.parquet", index=False)

    print("\n[SUCCESS] McCance & Widdowson integration files generated!")
    print("1. 'mccance_food_injection.parquet' -> Contains the food metadata.")
    print("2. 'mccance_food_nutrient_bucketed/' -> Contains the partitioned FDC measurements.")

if __name__ == "__main__":
    main()
