import pandas as pd
import shutil

original_food_path = "/data/bac2food/food.parquet"
fineli_food_path = "fineli_food_injection.parquet"

# Save a quick backup just in case
shutil.copy(original_food_path, "/data/bac2food/archive/food_backup2.parquet")

df_original = pd.read_parquet(original_food_path)
df_fineli = pd.read_parquet(fineli_food_path)

df_merged = pd.concat([df_original, df_fineli], ignore_index=True)
df_merged.to_parquet(original_food_path, index=False)
print("Fineli foods successfully merged!")
