import pandas as pd
import os

# --- Paths ---
FOOD_PARQUET = "/data/bac2food/food.parquet"
BUCKET_DIR = "/data/bac2food/food_nutrient_bucketed"

print("Merging WAFCT Food Meta...")
wafct_food = pd.read_csv("wafct_food.csv")
master_food = pd.read_parquet(FOOD_PARQUET)
master_food = master_food[~master_food['fdc_id'].isin(wafct_food['fdc_id'])]
pd.concat([master_food, wafct_food], ignore_index=True).to_parquet(FOOD_PARQUET, index=False)

print("Merging WAFCT Nutrients...")
wafct_nutrients = pd.read_csv("wafct_food_nutrient.csv")
wafct_nutrients['bucket'] = wafct_nutrients['nutrient_id'] % 256

for bucket_id, group in wafct_nutrients.groupby('bucket'):
    bucket_path = os.path.join(BUCKET_DIR, f"bucket={bucket_id}")
    os.makedirs(bucket_path, exist_ok=True)
    new_file = os.path.join(bucket_path, "wafct_data.parquet")
    group.drop(columns=['bucket']).to_parquet(new_file, index=False)

print("WAFCT successfully integrated!")
