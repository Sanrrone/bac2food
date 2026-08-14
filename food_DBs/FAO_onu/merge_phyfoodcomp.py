import pandas as pd
import os

# --- Set these paths to your master database locations ---
FOOD_PARQUET = "/data/bac2food/food.parquet"
BUCKET_DIR = "/data/bac2food/food_nutrient_bucketed"
NUTRIENT_CSV = "/data/bac2food/nutrient.csv"

def main():
    print("==========================================")
    print("INTEGRATING PHYFOODCOMP INTO MASTER GRAPH")
    print("==========================================")

    # 1. Merge Food Metadata
    print("\n[*] Merging food metadata...")
    phy_food = pd.read_csv("phyfoodcomp_food.csv")
    master_food = pd.read_parquet(FOOD_PARQUET)
    
    # Ensure no duplicates if run twice
    master_food = master_food[~master_food['fdc_id'].isin(phy_food['fdc_id'])]
    combined_food = pd.concat([master_food, phy_food], ignore_index=True)
    combined_food.to_parquet(FOOD_PARQUET, index=False)
    print(f" -> Added {len(phy_food)} PhyFoodComp foods to {FOOD_PARQUET}")

    # 2. Distribute Nutrients into Buckets
    print("\n[*] Merging nutrient data into Parquet buckets...")
    phy_nutrients = pd.read_csv("phyfoodcomp_food_nutrient.csv")
    
    # Calculate the Hive partition bucket (modulo 256)
    phy_nutrients['bucket'] = phy_nutrients['nutrient_id'] % 256

    for bucket_id, group in phy_nutrients.groupby('bucket'):
        bucket_path = os.path.join(BUCKET_DIR, f"bucket={bucket_id}")
        os.makedirs(bucket_path, exist_ok=True)
        
        # PyArrow Dataset automatically reads ALL parquet files inside a bucket folder.
        # We just drop a new file next to the existing ones!
        new_file = os.path.join(bucket_path, "phyfoodcomp_data.parquet")
        
        # Drop the 'bucket' column itself because Hive partitioning implies it
        group.drop(columns=['bucket']).to_parquet(new_file, index=False)

    print(f" -> Distributed {len(phy_nutrients)} nutrient records across partitions.")

    # 3. Register Custom Nutrient 99999
    print("\n[*] Registering Phytic Acid to nutrient.csv...")
    nut_df = pd.read_csv(NUTRIENT_CSV)
    
    if 99999 not in nut_df['id'].values:
        # Append to the CSV
        with open(NUTRIENT_CSV, 'a') as f:
            f.write("99999,Phytic Acid / Myo-Inositol,MG,\n")
        print(" -> Successfully registered Custom ID 99999.")
    else:
        print(" -> Custom ID 99999 already exists. Skipping.")

    print("\n[SUCCESS] Integration complete! You are ready to run the pipeline.")

if __name__ == "__main__":
    main()
