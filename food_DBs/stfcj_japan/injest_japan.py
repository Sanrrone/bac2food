#!/usr/bin/env python3
import pandas as pd
import numpy as np
import re
import os

# --- FDC Nutrient Mapping Dictionary ---
NUTRIENT_MAP = {
    # Proximates
    "Water": 1005,
    "Protein, calculated from reference nitrogen": 1003,
    "Lipid": 1004,
    "Ash": 1007,
    "Carbohydrate, available, total": 1050,
    "Fiber, total dietary": 1079,
    
    # Carbohydrates & Polyols
    "Starch": 1009, "Sucrose": 1010, "Glucose": 1011, "Fructose": 1012,
    "Lactose": 1013, "Maltose": 1014, "Galactose": 1070, "Sorbitol": 1061,
    "Mannitol": 1062, "Trehalose": 1060, 
    
    # Organic Acids & Polyphenols
    "Acetic acid": 1026, "Lactic acid": 1038, "Citric acid": 1033,
    "Malic acid": 1043, "Quinic acid": 1044, "Succinic acid": 1046,
    "Tartaric acid": 1045, "Chlorogenic acid": 1030, "Caffeic acid": 1201,
    "Ferulic acid": 1200, "p-Coumaric acid": 1202,
    
    # Amino Acids
    "Isoleucine": 1212, "Leucine": 1213, "Lysine": 1214, "Methionine": 1215,
    "Cystine": 1216, "Phenylalanine": 1217, "Tyrosine": 1218, "Threonine": 1211,
    "Tryptophan": 1210, "Valine": 1219, "Histidine": 1221, "Arginine": 1220,
    "Alanine": 1222, "Aspartic acid": 1223, "Glutamic acid": 1224, "Glycine": 1225,
    "Proline": 1226, "Serine": 1227, "Hydroxy-proline": 1228,
    "Porphyran": 99998,"Formic acid": 95000,
    "Gluconic acid": 95001,
    "Malonic acid": 95002,
    "alpha-Ketoglutaric acid": 95003,
    "Orotic acid": 95004,
    "Ferulic acid": 95005, # Note: if you already mapped Ferulic to 1200 earlier, keep 1200!
    "Tannin": 95006,
    "Polyphenol": 95007,
}

def clean_numeric(val):
    """Converts Japanese STFCJ text values (Traces, Estimates) to clean floats."""
    if pd.isna(val): return np.nan
    val = str(val).strip()
    if val in ('-', '', 'ND', 'N.D.'): return np.nan
    if val.lower() == 'tr': return 0.001 # Trace amounts
    
    # STFCJ wraps estimated values in parentheses e.g., "(1.2)" -> "1.2"
    val = re.sub(r'[\(\)]', '', val)
    # Remove any non-numeric chars except decimals and minus signs
    val = re.sub(r'[^\d\.\-]', '', val)
    try: return float(val)
    except ValueError: return np.nan

def load_stfcj_excel(filepath, sheet_name):
    """Loads an Excel sheet and dynamically finds the correct header row."""
    if not os.path.exists(filepath):
        print(f"[!] Warning: File not found: {filepath}")
        return pd.DataFrame()
        
    print(f"Reading {filepath} | Tab: '{sheet_name}' ...")
    
    # Read the first 15 rows to locate the header row containing "Item No."
    df_raw = pd.read_excel(filepath, sheet_name=sheet_name, header=None, nrows=15)
    header_idx = -1
    for i, row in df_raw.iterrows():
        # CRITICAL FIX: Loop through x in row.values
        row_str = " ".join([" ".join(str(x).split()) for x in row.values])
        if "Item No." in row_str:
            header_idx = i
            break
            
    if header_idx == -1:
        print(f"[!] Error: Could not find 'Item No.' header row in {filepath} ({sheet_name})")
        return pd.DataFrame()
        
    # Read the actual dataframe using the discovered header index
    df = pd.read_excel(filepath, sheet_name=sheet_name, header=header_idx)
    
    # Clean column names by collapsing all whitespace/newlines into a single space
    df.columns = [" ".join(str(c).split()) for c in df.columns]
    
    if "Item No." in df.columns:
        df = df.rename(columns={"Item No.": "item_no", "Food and Description": "description"})
    else:
        print(f"[!] Fatal Error: 'Item No.' missing after cleaning! Found: {df.columns.tolist()}")
        return pd.DataFrame()
        
    # Drop rows that don't have a valid Item No
    df = df.dropna(subset=['item_no'])
    df['item_no'] = pd.to_numeric(df['item_no'], errors='coerce')
    df = df.dropna(subset=['item_no'])
    df['item_no'] = df['item_no'].astype(int)
    
    return df

def load_scientific_names(filepath):
    if not os.path.exists(filepath): 
        print(f"[!] Warning: Scientific names file not found: {filepath}")
        return pd.DataFrame()
        
    print(f"Reading {filepath} | Tab: 'Scientific name of food source' ...")
    df = pd.read_excel(filepath, sheet_name="Scientific name of food source", header=2)
    
    # CRITICAL FIX: Clean headers here too
    df.columns = [" ".join(str(c).split()) for c in df.columns]
    
    if "Item No." in df.columns:
        df = df.rename(columns={"Item No.": "item_no", "Scientific name": "sci_name"})
    
    df = df.dropna(subset=['item_no'])
    # Handle strings like "14002, 14003" -> Extract the first ID
    df['item_no'] = df['item_no'].astype(str).str.extract(r'(\d+)')[0]
    df = df.dropna(subset=['item_no'])
    df['item_no'] = df['item_no'].astype(int)
    
    return df[['item_no', 'sci_name']]

def main():
    print("==================================================")
    print("INGESTING JAPANESE STFCJ DATABASE (.xlsx Native)")
    print("==================================================")

    # 1. Target exactly the tabs that use the "per 100g EP (Edible Portion)" standard
    targets = [
        ("main_1374049_1r12_1.xlsx", "Table"),
        ("aminoacid_1374049_2r11_1.xlsx", "Table 1(per 100 g EP)"),
        ("org_acid_1388558_4r12r.xlsx", "Table"),
        ("org_acid_1388558_4r12r.xlsx", "Annex (organic acids)")
    ]
    
    dfs = []
    for file, sheet in targets:
        df = load_stfcj_excel(file, sheet)
        if not df.empty: dfs.append(df)
    
    if not dfs:
        print("No data loaded. Check files and paths.")
        return

    print("\nMerging tabs on Item No...")
    master_df = dfs[0]
    for df in dfs[1:]:
        cols_to_use = [c for c in df.columns if c not in master_df.columns or c == 'item_no']
        master_df = pd.merge(master_df, df[cols_to_use], on='item_no', how='outer')

    # Add Scientific Names
    sci_df = load_scientific_names("scientific_names.xlsx")
    if not sci_df.empty:
        master_df = pd.merge(master_df, sci_df, on='item_no', how='left')

    # 2. Build Food Meta
    print("Formatting Food Meta...")
    # Add offset to prevent FDC ID collisions (Japan code = 81)
    raise SystemExit(
        "This v1 ingest mints fdc_ids from a literal offset, which the block scheme "
        "no longer permits: ids are accessions held in fdc_id_map.tsv. Run injest_japan_v2.py "
        "instead, which allocates through food_DBs/fdc_blocks.py."
    ) 
    
    food_meta = master_df[['fdc_id', 'description']].copy()
    if 'sci_name' in master_df.columns:
        food_meta['description'] = food_meta['description'] + master_df['sci_name'].apply(lambda x: f" [{x}]" if pd.notna(x) else "")
    
    food_meta['description'] = food_meta['description'] + " [STFCJ]"
    food_meta['data_type'] = "foundation_food"
    
    # 2-Digit logic for Japanese food categories
    cat_map = {
        1: "Cereal Grains and Pasta", 2: "Vegetables and Vegetable Products",
        3: "Sweets", 4: "Legumes and Legume Products", 5: "Nut and Seed Products",
        6: "Vegetables and Vegetable Products", 7: "Fruits and Fruit Juices",
        8: "Vegetables and Vegetable Products", # Mushrooms
        9: "Vegetables and Vegetable Products", # Seaweeds
        10: "Finfish and Shellfish Products", 11: "Beef Products",
        17: "Spices and Herbs"
    }
    master_df['group_code'] = (master_df['item_no'] // 1000).astype(int)
    food_meta['food_category'] = master_df['group_code'].map(cat_map).fillna("Vegetables and Vegetable Products")

    # 3. Build Food Nutrient Matrix
    print("Mapping Nutrients (Unpivoting)...")
    nutrient_records = []
    available_cols = [col for col in NUTRIENT_MAP.keys() if col in master_df.columns]
    
    for _, row in master_df.iterrows():
        fid = row['fdc_id']
        for col in available_cols:
            val = clean_numeric(row[col])
            if pd.notna(val) and val > 0:
                nutrient_records.append({
                    "fdc_id": fid,
                    "nutrient_id": NUTRIENT_MAP[col],
                    "amount": val
                })

    food_nutrients = pd.DataFrame(nutrient_records)
    print(f"\nExtracted {len(food_meta)} foods and {len(food_nutrients)} nutrient records.")
    
    # 4. Save to CSV
    food_meta.to_csv("japan_food.csv", index=False)
    food_nutrients.to_csv("japan_food_nutrient.csv", index=False)
    
    print("[*] Successfully generated japan_food.csv and japan_food_nutrient.csv")

if __name__ == "__main__":
    main()