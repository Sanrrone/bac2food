#!/usr/bin/env python3
import pandas as pd
import numpy as np
import re
import os

# FDC Nutrient Mapping
# We are minting a CUSTOM FDC ID (99999) for Phytic Acid / Inositol Hexaphosphate 
# because standard government FDC databases do not track it.
NUTRIENT_MAP = {
    "WATER(g)": 1005,
    "FE(mg)": 1089,
    "ZN(mg)": 1095,
    "CA(mg)": 1087,
    "phytic_acid_best": 99999  # Custom ID for Inositol/Phytic Acid
}

def clean_numeric(val):
    """Clean the text values, converting traces and <LOD to trace numeric amounts."""
    if pd.isna(val): return np.nan
    val = str(val).strip()
    if val in ('', '-', 'ND', 'N.D.'): return np.nan
    if val.lower() in ('tr', '< lod', '<lod'): return 0.001
    
    # Strip brackets e.g., "[71,10]" taking the first value
    val = re.sub(r'[\[\]]', '', val).split(',')[0]
    # Remove any non-numeric chars except decimals and minus signs
    val = re.sub(r'[^\d\.\-]', '', val)
    try: return float(val)
    except ValueError: return np.nan

def load_fao_sheet(filepath, sheet_name):
    """Dynamically finds the header and loads an FAO data sheet."""
    print(f"Reading tab: '{sheet_name}'...")
    
    # Read the first 10 rows to locate the header row containing "Food item ID"
    df_raw = pd.read_excel(filepath, sheet_name=sheet_name, header=None, nrows=10)
    header_idx = -1
    for i, row in df_raw.iterrows():
        row_str = " ".join([str(x) for x in row.values])
        if "Food item ID" in row_str:
            header_idx = i
            break
            
    if header_idx == -1:
        print(f"[!] Warning: Could not find 'Food item ID' header in {sheet_name}. Skipping.")
        return pd.DataFrame()
        
    # Read the actual dataframe using the discovered header
    df = pd.read_excel(filepath, sheet_name=sheet_name, header=header_idx)
    
    # Clean column names (strip whitespace and newlines)
    df.columns = [" ".join(str(c).split()) for c in df.columns]
    
    if 'Food item ID' not in df.columns:
        return pd.DataFrame()
        
    # Drop rows that don't have a valid numeric Food item ID (skipping sub-headers or footers)
    df = df.dropna(subset=['Food item ID'])
    df['Food item ID'] = pd.to_numeric(df['Food item ID'], errors='coerce')
    df = df.dropna(subset=['Food item ID'])
    
    return df

def main():
    filepath = "PhyFoodComp_1.0.xlsx"
    if not os.path.exists(filepath):
        print(f"[!] Error: {filepath} not found in the current directory.")
        return

    print("==================================================")
    print("INGESTING FAO/INFOODS PhyFoodComp 1.0 (Excel Native)")
    print("==================================================")

    # 1. Discover data tabs (They all start with a 2-digit number, e.g., "01 Cereals...")
    xls = pd.ExcelFile(filepath)
    data_sheets = [s for s in xls.sheet_names if re.match(r'^\d{2}\s', s)]
    
    dfs = []
    for sheet in data_sheets:
        df = load_fao_sheet(filepath, sheet)
        if not df.empty:
            dfs.append(df)

    if not dfs:
        print("[!] No valid data extracted.")
        return

    master_df = pd.concat(dfs, ignore_index=True)
    
    # 2. Coalesce Phytic Acid Measurements
    # PhyFoodComp contains multiple analytical methods. We want the most accurate one available per food.
    phytate_cols = ['IP6(mg)', 'IPSUM (mg)', 'PHYTCA (mg)', 'PHYTCPP(mg)', 'PHYTCPPI(mg)', 'PHYTCPPD(mg)', 'PHYTC-(mg)', 'PHYT-']
    
    master_df['phytic_acid_best'] = np.nan
    for col in phytate_cols:
        if col in master_df.columns:
            numeric_col = master_df[col].apply(clean_numeric)
            master_df['phytic_acid_best'] = master_df['phytic_acid_best'].fillna(numeric_col)

    # 3. Build Food Meta
    print("\nFormatting Food Meta...")
    # Offset IDs by 92,000,000 to prevent collisions with FDC (USDA) and Japan (81m)
    raise SystemExit(
        "This v1 ingest mints fdc_ids from a literal offset, which the block scheme "
        "no longer permits: ids are accessions held in fdc_id_map.tsv. Run injest_phyfood_v2.py "
        "instead, which allocates through food_DBs/fdc_blocks.py."
    )
    
    food_meta = pd.DataFrame()
    food_meta['fdc_id'] = master_df['fdc_id']
    
    names = master_df.get('Food name in English', master_df.get('Food name in own language', pd.Series(['Unknown']*len(master_df)))).astype(str).replace('nan', 'Unknown Food')
    species = master_df.get('Species/Subspecies', pd.Series(['']*len(master_df))).astype(str).replace('nan', '')
    
    desc = names + " [" + species + "]"
    desc = desc.str.replace(r' \[\]', '', regex=True) # clean empty brackets
    food_meta['description'] = desc + " [PhyFoodComp]"
    food_meta['data_type'] = "foundation_food"
    
    # Map FAO Group Codes to USDA Categories
    cat_map = {
        1: "Cereal Grains and Pasta", 2: "Vegetables and Vegetable Products",
        3: "Legumes and Legume Products", 4: "Vegetables and Vegetable Products",
        5: "Fruits and Fruit Juices", 6: "Nut and Seed Products",
        7: "Beef Products", 8: "American Indian/Alaska Native Foods", 
        10: "Finfish and Shellfish Products", 11: "Dairy and Egg Products",
        12: "Fats and Oils", 13: "Beverages", 14: "Sweets", 15: "Spices and Herbs",
        19: "Meals, Entrees, and Side Dishes"
    }
    master_df['Food Group'] = pd.to_numeric(master_df.get('Food Group', 0), errors='coerce')
    food_meta['food_category'] = master_df['Food Group'].map(cat_map).fillna("Vegetables and Vegetable Products")

    # 4. Build Food Nutrient Matrix
    print("Extracting Phytic Acid and baseline minerals...")
    nutrient_records = []
    
    for col, nutrient_id in NUTRIENT_MAP.items():
        if col not in master_df.columns: continue
        
        temp = master_df[['fdc_id', col]].copy()
        if col != 'phytic_acid_best':
            temp[col] = temp[col].apply(clean_numeric)
            
        temp = temp.dropna(subset=[col])
        temp = temp[temp[col] > 0]
        
        for _, row in temp.iterrows():
            nutrient_records.append({
                "fdc_id": int(row['fdc_id']),
                "nutrient_id": nutrient_id,
                "amount": row[col]
            })

    food_nutrients = pd.DataFrame(nutrient_records)
    print(f"\nExtracted {len(food_meta)} foods and {len(food_nutrients)} nutrient records.")
    
    food_meta.drop_duplicates(subset=['fdc_id']).to_csv("phyfoodcomp_food.csv", index=False)
    food_nutrients.drop_duplicates().to_csv("phyfoodcomp_food_nutrient.csv", index=False)
    print("[*] Successfully generated phyfoodcomp_food.csv and phyfoodcomp_food_nutrient.csv")

if __name__ == "__main__":
    main()
