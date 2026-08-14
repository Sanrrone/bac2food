#!/usr/bin/env python3
import pandas as pd
import numpy as np
import re
import os

# fdc_id allocation is centralised in food_DBs/fdc_blocks.py. Never write a literal offset:
# ids are accessions looked up in fdc_id_map.tsv, so a food keeps its id across releases.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
import fdc_blocks


# --- FDC Nutrient Mapping Dictionary ---
# Maps unique substrings of the WAFCT column headers to FDC IDs
NUTRIENT_MAP = {
    "Water": 1005,
    "Protein, total": 1003,
    "Fat, total": 1004,
    "Carbohydrate, available": 1050,
    "Fibre, total dietary": 1079,
    "Ash": 1007,
    "Calcium": 1087,
    "Iron": 1089,
    "Magnesium": 1090,
    "Phosphorus": 1091,
    "Potassium": 1092,
    "Sodium": 1093,
    "Zinc": 1095,
    "Copper": 1098,
    "Retinol\n": 1104,
    "Beta-carotene": 1107,
    "Vitamin D": 1114,
    "Thiamine": 1165,
    "Riboflavin": 1166,
    "Niacin, preformed": 1167,
    "Vitamin B6": 1175,
    "Folate, total": 1177,
    "Vitamin B12": 1178,
    "Vitamin C": 1162,
    "Cholesterol": 1253,
    "Phytate": 99999 # Mapping to our custom Inositol/Phytate ID!
}

def clean_numeric(val):
    """Clean the text values, removing brackets, traces, and handling blanks."""
    if pd.isna(val): return np.nan
    val = str(val).strip()
    
    # Handle missing/trace strings
    if val in ('', '-', 'N', 'ND', 'N.D.'): return np.nan
    if val.lower() in ('tr', 'trace', '< lod', '<lod'): return 0.001
    
    # Remove FAO estimation brackets e.g., "[1.4]" -> "1.4"
    val = re.sub(r'[\[\]]', '', val)
    # Remove any non-numeric chars except decimals and minus signs
    val = re.sub(r'[^\d\.\-]', '', val)
    try: 
        return float(val)
    except ValueError: 
        return np.nan

def wafct_code(code_str):
    """Normalise a WAFCT code like '01_172' to the digits-only form used as its accession key.

    This used to add a 120,000,000 base and return an fdc_id directly, after a 94,000,000 base
    put 27 WAFCT foods on top of PhyFoodComp ones (PhyFoodComp adds a native code reaching
    19,020,060 to its own base, so it sprawled from 93M to 111M). Blocks are no longer any
    single ingest's business -- fdc_blocks owns them, and they are uniform.
    """
    clean_code = re.sub(r'[^\d]', '', str(code_str))
    if not clean_code: return np.nan
    return int(clean_code)

def main():
    filepath = "WAFCT_2019.xlsx"
    sheet_name = "05 NV_sum_57 (per 100g EP)" # The master 57-nutrient tab
    
    if not os.path.exists(filepath):
        print(f"[!] Error: {filepath} not found in the current directory.")
        return

    print("==================================================")
    print("INGESTING FAO WEST AFRICAN DATABASE (WAFCT 2019)")
    print("==================================================")

    print(f"Reading tab: '{sheet_name}'...")
    
    # WAFCT headers are typically on row 1 (index 0)
    df = pd.read_excel(filepath, sheet_name=sheet_name, header=0)
    
    # Drop rows without a valid 'Code' (skipping category header rows like 'Cereals...')
    df = df.dropna(subset=['Code'])
    df = df[df['Code'].astype(str).str.contains(r'\d+_\d+')]
    
    print("\nFormatting Food Meta...")
    
    # 1. Build Food Meta
    food_meta = pd.DataFrame()
    # Assign once, then reuse. The nutrient loop below needs the same accessions; recomputing
    # them there would only work while minting was a pure function of the code.
    codes = df['Code'].apply(wafct_code)
    food_meta['fdc_id'] = fdc_blocks.assign('wafct', codes)
    code2fdc = dict(zip(codes, food_meta['fdc_id']))
    
    # Combine English name and Scientific name
    names = df['Food name in English'].astype(str)
    species = df['Scientific name'].astype(str).replace('nan', '')
    
    desc = names
    # Add scientific name if present
    mask = species != ''
    desc[mask] = desc[mask] + " [" + species[mask] + "]"
    
    food_meta['description'] = desc + " [WAFCT]"
    food_meta['data_type'] = "foundation_food"
    
    # Extract Group Code (First 2 digits of the code string, e.g., '01_172' -> 1)
    df['group_code'] = df['Code'].astype(str).str.split('_').str[0]
    df['group_code'] = pd.to_numeric(df['group_code'], errors='coerce')
    
    # Map FAO WAFCT Groups to USDA Categories
    cat_map = {
        1: "Cereal Grains and Pasta", 2: "Vegetables and Vegetable Products", # Roots/tubers
        3: "Legumes and Legume Products", 4: "Vegetables and Vegetable Products",
        5: "Fruits and Fruit Juices", 6: "Nut and Seed Products",
        7: "Beef Products", 8: "American Indian/Alaska Native Foods", # Insects and grubs
        9: "Finfish and Shellfish Products", 10: "Dairy and Egg Products",
        11: "Fats and Oils", 12: "Beverages", 13: "Spices and Herbs",
        14: "Meals, Entrees, and Side Dishes" # Recipes
    }
    food_meta['food_category'] = df['group_code'].map(cat_map).fillna("Vegetables and Vegetable Products")

    # 2. Build Food Nutrient Matrix
    print("Mapping Nutrients (Unpivoting)...")
    nutrient_records = []
    
    # Map exact column names using our substring dictionary
    col_mapping = {}
    for actual_col in df.columns:
        for substr, fdc_id in NUTRIENT_MAP.items():
            if substr.lower() in str(actual_col).lower():
                col_mapping[actual_col] = fdc_id
                break
    
    for _, row in df.iterrows():
        fid = code2fdc.get(wafct_code(row['Code']))
        if fid is None or pd.isna(fid): continue
        fid = int(fid)
        
        for actual_col, nutrient_id in col_mapping.items():
            val = clean_numeric(row[actual_col])
            if pd.notna(val) and val > 0:
                nutrient_records.append({
                    "fdc_id": fid,
                    "nutrient_id": nutrient_id,
                    "amount": val
                })

    food_nutrients = pd.DataFrame(nutrient_records)
    print(f"\nExtracted {len(food_meta)} foods and {len(food_nutrients)} nutrient records.")
    
    # 3. Save to CSV
    food_meta.drop_duplicates(subset=['fdc_id']).to_csv("wafct_food.csv", index=False)
    food_nutrients.drop_duplicates().to_csv("wafct_food_nutrient.csv", index=False)
    
    print("[*] Successfully generated wafct_food.csv and wafct_food_nutrient.csv")

if __name__ == "__main__":
    main()
