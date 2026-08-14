import pandas as pd

# Load existing master food table
master_food = pd.read_parquet("/data/bac2food/food.parquet")

# Load Japan food table
japan_food = pd.read_csv("japan_food.csv")

# Combine and save back
combined_food = pd.concat([master_food, japan_food], ignore_index=True)
combined_food.to_parquet("/data/bac2food/food.parquet", index=False)
print("Added Japan to master food.parquet!")
