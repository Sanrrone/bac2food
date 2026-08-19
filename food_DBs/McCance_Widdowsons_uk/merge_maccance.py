import pandas as pd
import shutil
import os
from pathlib import Path

# --- PATHS ---
original_food_path = "/data/bac2food/food.parquet"
mccance_food_path = "mccance_food_injection.parquet"
backup_path = "/data/bac2food/archive/food_backup_mccance.parquet"

src_buckets = "mccance_food_nutrient_bucketed"
dest_buckets = "/data/bac2food/food_nutrient_bucketed"
index_dir = Path("/data/bac2food/index_modeled")

print("==================================================")
print("INTEGRATING McCANCE DATA INTO MAIN FDC REPOSITORY")
print("==================================================")

print("\n[1/4] Creating a backup of the original food.parquet...")
shutil.copy(original_food_path, backup_path)

print("[2/4] Fusing the food dictionaries...")
df_original = pd.read_parquet(original_food_path)
df_mccance = pd.read_parquet(mccance_food_path)

df_merged = pd.concat([df_original, df_mccance], ignore_index=True)

# Safety check: Drop duplicates and rogue category columns
if 'food_category' in df_merged.columns and 'food_category_id' in df_merged.columns:
    df_merged = df_merged.drop(columns=['food_category'])
df_merged = df_merged.drop_duplicates(subset=["fdc_id"], keep="last")

df_merged.to_parquet(original_food_path, index=False)
print(f"  -> Success! Food dictionary expanded to {len(df_merged)} items.")

print("\n[3/4] Copying McCance nutrient buckets into main FDC buckets...")
shutil.copytree(src_buckets, dest_buckets, dirs_exist_ok=True)
print("  -> Buckets merged successfully!")

print("\n[4/4] Clearing stale cache and index files in /index_modeled/...")
if index_dir.exists():
    for ext in ["*.parquet", "*.pkl"]:
        for file_path in index_dir.glob(ext):
            try:
                file_path.unlink()
                print(f"  -> Deleted old cache file: {file_path.name}")
            except Exception as e:
                print(f"  -> Error deleting {file_path.name}: {e}")

print("\n==================================================")
print("ALL DONE! McCance & Widdowson is fully integrated.")
print("You can now run ec2food_v42.py with the --rebuild-static-meta flag!")
print("==================================================")
