import pandas as pd
import shutil

# Define your paths
original_food_path = "/data/bac2food/food.parquet"
pe_food_path = "pe_food_injection.parquet"
backup_path = "/data/bac2food/archive/food_backup.parquet"

print("1. Creating a backup of the original food.parquet...")
shutil.copy(original_food_path, backup_path)

print("2. Loading the datasets...")
df_original = pd.read_parquet(original_food_path)
df_pe = pd.read_parquet(pe_food_path)

print("3. Fusing the datasets...")
# Combine them together
df_merged = pd.concat([df_original, df_pe], ignore_index=True)

print("4. Saving the upgraded food.parquet...")
df_merged.to_parquet(original_food_path, index=False)

print(f"Success! Your database just grew from {len(df_original)} to {len(df_merged)} foods.")
