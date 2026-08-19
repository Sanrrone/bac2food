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



# NEVO id offset, category map and nutrient map. Module-level and shared:
# reconstruct_nevo.py imports these so a user's local reconstruction of the NEVO partition is
# identical to what the build pipeline produces. Duplicating them is how a deposit and a local
# build silently diverge.
DEFAULT_CATEGORY = 'Meals, Entrees, and Side Dishes'

CAT_MAP = {
    'Vegetables': 'Vegetables and Vegetable Products',
    'Potatoes and tubers': 'Vegetables and Vegetable Products',
    'Legumes': 'Legumes and Legume Products',
    'Meat substitutes and dairy alternatives': 'Legumes and Legume Products',
    'Fruits': 'Fruits and Fruit Juices',
    'Nuts and seeds': 'Nut and Seed Products',
    'Bread': 'Cereal Grains and Pasta',
    'Cereal products': 'Cereal Grains and Pasta',
    'Milk and milk products': 'Dairy and Egg Products',
    'Cheese': 'Dairy and Egg Products',
    'Meat, meat products and poultry': 'Sausages and Luncheon Meats',
    'Fish': 'Finfish and Shellfish Products',
    'Fats and oils': 'Fats and Oils',
    'Non-alcoholic beverages': 'Beverages',
    'Alcoholic beverages': 'Alcoholic Beverages',
    'Sugar, candy, sweet spreads and sauces': 'Sweets',
    'Cakes and pastries': 'Baked Products',
    'Herbs and spices': 'Spices and Herbs',
    'Soups': 'Soups, Sauces, and Gravies',
    'Sauces and savory spreads': 'Soups, Sauces, and Gravies',
    'Savory snacks': 'Snacks',
    'Miscellaneous': 'Meals, Entrees, and Side Dishes'
}

NUTRIENT_MAP = {
    'PROT (g)': 1003,      # Protein
    'FAT (g)': 1004,       # Total lipid
    'CHO (g)': 1005,       # Carbohydrate
    'STARCH (g)': 1009,    # Starch
    'FIBT (g)': 1079,      # Fiber
    'SUGAR (g)': 2000,     # Sugars, total
    'ALC (g)': 1018,       # Alcohol
    'CA (mg)': 1087,       # Calcium
    'FE (mg)': 1089,       # Iron
    'MG (mg)': 1090,       # Magnesium
    'P (mg)': 1091,        # Phosphorus
    'K (mg)': 1092,        # Potassium
    'NA (mg)': 1093,       # Sodium
    'ZN (mg)': 1095,       # Zinc
    'CU (mg)': 1098,       # Copper
    'SE (µg)': 1103,       # Selenium
    'VITA_RAE (µg)': 1104, # Vitamin A
    'VITD (µg)': 1110,     # Vitamin D
    'VITC (mg)': 1162,     # Vitamin C
    'THIA (mg)': 1165,     # Thiamin
    'RIBF (mg)': 1166,     # Riboflavin
    'NIA (mg)': 1167,      # Niacin
    'VITB6 (mg)': 1175,    # Vitamin B-6
    'FOL (µg)': 1177,      # Folate
    'VITB12 (µg)': 1174,   # Vitamin B-12
    'CHORL (mg)': 1253,    # Cholesterol
    'FASAT (g)': 1258,     # Saturated Fat
    'FAMSCIS (g)': 1259,   # Monounsaturated Fat
    'FAPU (g)': 1260       # Polyunsaturated Fat
}


def main():
    print("==================================================")
    print("NEVO (NETHERLANDS) TO FDC INGESTION SCRIPT")
    print("==================================================")
    
    print("[1/5] Loading NEVO composition data from Excel...")
    # Read directly from the first tab of the Excel file
    excel_file = "NEVO2025_v9.0.xlsx"
    nevo_df = pd.read_excel(excel_file, sheet_name=0)
    
    print("[2/5] Minting custom FDC IDs for NEVO foods...")
    # Offset NEVO IDs by 30,000,000 to prevent collisions
    nevo_df['fdc_id'] = fdc_blocks.assign('nevo', nevo_df['NEVO-code'].astype(int))
    
    print("[3/5] Mapping NEVO categories...")
    cat_map = CAT_MAP
    nevo_df['food_category'] = nevo_df['Food group'].map(cat_map).fillna(DEFAULT_CATEGORY)
    
    print("[4/5] Mapping NEVO nutrients to FDC IDs...")
    target_mapping = NUTRIENT_MAP
    # Extract only available columns
    available_cols = [c for c in target_mapping.keys() if c in nevo_df.columns]
    
    # Ensure numeric types (NEVO might use strings like "sp" for traces or use commas)
    for c in available_cols:
        if nevo_df[c].dtype == object:
            # Replace traces or non-numeric markers with 0.0, replace commas with dots
            nevo_df[c] = nevo_df[c].astype(str).str.replace(r'^[a-zA-Z<>].*', '0.0', regex=True)
            nevo_df[c] = nevo_df[c].str.replace(',', '.', regex=False)
            nevo_df[c] = pd.to_numeric(nevo_df[c], errors='coerce')
    
    # Melt the dataframe into long format
    melted_df = nevo_df.melt(id_vars=['fdc_id', 'Engelse naam/Food name', 'food_category'], 
                             value_vars=available_cols, 
                             var_name='nevo_nutrient', 
                             value_name='amount')
    
    # Drop NAs and zeroes
    melted_df = melted_df.dropna(subset=['amount'])
    melted_df = melted_df[melted_df['amount'] > 0]
    
    # Map to numeric FDC IDs
    melted_df['nutrient_id'] = melted_df['nevo_nutrient'].map(target_mapping).astype(int)
    
    print("[5/5] Generating database files...")
    # Create the food dictionary injection file
    unique_foods = melted_df[['fdc_id', 'Engelse naam/Food name', 'food_category']].drop_duplicates()
    new_foods = pd.DataFrame({
        'fdc_id': unique_foods['fdc_id'],
        'data_type': 'foundation_food', 
        'description': unique_foods['Engelse naam/Food name'] + " [NEVO]",
        'food_category_id': '6666', # Custom category ID for NEVO
        'food_category': unique_foods['food_category']
    })
    
    new_foods.to_parquet("nevo_food_injection.parquet", index=False)
    
    # Build bucketed Parquets
    melted_df['bucket'] = melted_df['nutrient_id'] % 256
    
    out_dir = Path("nevo_food_nutrient_bucketed")
    out_dir.mkdir(exist_ok=True)
    
    for bucket_id, group in melted_df.groupby('bucket'):
        bucket_dir = out_dir / f"bucket={bucket_id}"
        bucket_dir.mkdir(exist_ok=True)
        group = group.sort_values(['fdc_id', 'nutrient_id'])
        group[['fdc_id', 'nutrient_id', 'amount']].to_parquet(bucket_dir / "nevo_data.parquet", index=False)

    print("\n[SUCCESS] NEVO integration files generated!")
    print("1. 'nevo_food_injection.parquet' -> Contains the food metadata.")
    print("2. 'nevo_food_nutrient_bucketed/' -> Contains the partitioned FDC measurements.")

if __name__ == "__main__":
    main()