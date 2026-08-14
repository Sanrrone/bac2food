#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq

def main():
    print("==================================================")
    print("PHENOL-EXPLORER TO FDC INGESTION SCRIPT")
    print("==================================================")
    
    # 1. LOAD PHENOL EXPLORER DATA
    print("[1/5] Loading Phenol-Explorer composition data...")
    pe_df = pd.read_csv("composition-data.tsv", sep="\t")
    
    # 2. CREATE FDC-COMPATIBLE FOOD IDS
    print("[2/5] Minting custom FDC IDs for Phenol-Explorer foods...")
    unique_foods = pe_df[['food', 'food_group']].drop_duplicates().reset_index(drop=True)
    raise SystemExit(
        "This v1 ingest mints fdc_ids from a literal offset, which the block scheme "
        "no longer permits: ids are accessions held in fdc_id_map.tsv. Run ingest_phenol_explorer_v2.py "
        "instead, which allocates through food_DBs/fdc_blocks.py."
    )
    
    # Map Phenol-Explorer categories to FDC Categories
    cat_map = {
        'Fruits and fruit products': 'Fruits and Fruit Juices',
        'Vegetables': 'Vegetables and Vegetable Products',
        'Seasonings': 'Spices and Herbs',
        'Seeds': 'Nut and Seed Products',
        'Non-alcoholic beverages': 'Beverages',
        'Alcoholic beverages': 'Alcoholic Beverages',
        'Cereals and cereal products': 'Cereal Grains and Pasta',
        'Oils': 'Fats and Oils',
        'Coffee and cocoa': 'Beverages'
    }
    unique_foods['food_category'] = unique_foods['food_group'].map(cat_map).fillna('Vegetables and Vegetable Products')
    
    # Merge the new IDs back into the main dataset
    pe_df = pe_df.merge(unique_foods[['food', 'fdc_id']], on='food', how='left')
    
    # 3. MAP PHENOL-EXPLORER COMPOUNDS TO FDC TARGETS
    print("[3/5] Mapping Phenol-Explorer groups to FDC targets...")
    # These are the exact nutrient IDs Faecalibacterium and others were looking for!
    target_mapping = {
        'Flavonoids': 1347,          # Flavonoids, total
        'Flavonols': 1386,           # Flavonols, total
        'Flavanols': 1340,           # Catechins/Flavanols
        'Anthocyanins': 1411,        # Anthocyanidins
        'Isoflavonoids': 2049,       # Daidzin/Genistin proxy
        'Phenolic acids': 1345,      # Phenolic acids (approximate to FDC)
        'Lignans': 1346              # Lignans
    }
    
    # We will map by compound_group first, then compound_sub_group
    pe_df['fdc_nutrient_id'] = pe_df['compound_sub_group'].map(target_mapping)
    pe_df['fdc_nutrient_id'] = pe_df['fdc_nutrient_id'].fillna(pe_df['compound_group'].map(target_mapping))
    
    # Drop compounds we don't have an FDC target for
    mapped_df = pe_df.dropna(subset=['fdc_nutrient_id']).copy()
    mapped_df['fdc_nutrient_id'] = mapped_df['fdc_nutrient_id'].astype(int)
    
    # Aggregate the means and RENAME to 'nutrient_id' to match FDC format
    nutrient_data = mapped_df.groupby(['fdc_id', 'fdc_nutrient_id'])['mean'].sum().reset_index()
    nutrient_data.rename(columns={'mean': 'amount', 'fdc_nutrient_id': 'nutrient_id'}, inplace=True)
    
    # 4. APPEND TO FDC FOOD.PARQUET
    print("[4/5] Appending new foods to food.parquet...")
    # Create the new food rows
    new_foods = pd.DataFrame({
        'fdc_id': unique_foods['fdc_id'],
        'data_type': 'foundation_food', # Treat PE data as Foundation (High Quality)
        'description': unique_foods['food'] + " [Phenol-Explorer]",
        'food_category_id': '9999', # Custom category ID
        'food_category': unique_foods['food_category']
    })
    
    new_foods.to_parquet("pe_food_injection.parquet", index=False)
    
    # 5. BUILD THE BUCKETED PARQUET FILES
    print("[5/5] Generating bucketed nutrient partitions...")
    # Add the bucket hash
    nutrient_data['bucket'] = nutrient_data['nutrient_id'] % 256
    
    # Save the measurements in exactly the format your scanner expects
    out_dir = Path("pe_food_nutrient_bucketed")
    out_dir.mkdir(exist_ok=True)
    
    for bucket_id, group in nutrient_data.groupby('bucket'):
        bucket_dir = out_dir / f"bucket={bucket_id}"
        bucket_dir.mkdir(exist_ok=True)
        # Sort for rapid Arrow scanning
        group = group.sort_values(['fdc_id', 'nutrient_id'])
        group[['fdc_id', 'nutrient_id', 'amount']].to_parquet(bucket_dir / "phenol_explorer_data.parquet", index=False)

    print("\n[SUCCESS] Phenol-Explorer integration files generated!")
    print("1. 'pe_food_injection.parquet' -> Contains the new food dictionaries.")
    print("2. 'pe_food_nutrient_bucketed/' -> Contains the partitioned FDC measurements.")

if __name__ == "__main__":
    main()