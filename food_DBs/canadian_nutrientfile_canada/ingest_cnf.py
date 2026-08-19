#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
from pathlib import Path

# fdc_id allocation: fdc_blocks.py
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
import fdc_blocks


def main():
    print("==================================================")
    print("CNF (CANADA) TO FDC INGESTION SCRIPT")
    print("==================================================")

    print("[1/5] Loading CNF files with Latin-1 encoding...")
    # Loading the relational CSVs
    df_food = pd.read_csv("FOOD NAME.csv", encoding='latin-1', usecols=['FoodID', 'FoodGroupID', 'FoodDescription'])
    df_group = pd.read_csv("FOOD GROUP.csv", encoding='latin-1', usecols=['FoodGroupID', 'FoodGroupName'])
    df_amt = pd.read_csv("NUTRIENT AMOUNT.csv", usecols=['FoodID', 'NutrientID', 'NutrientValue'])

    print("[2/5] Minting custom FDC IDs and mapping categories...")
    # Offset FoodID by 60,000,000 to prevent collisions
    df_food['fdc_id'] = fdc_blocks.assign('cnf', df_food['FoodID'].astype(int))
    
    # Merge food and group to get the english category names
    df_food = df_food.merge(df_group, on='FoodGroupID', how='left')

    # Mapping CNF groups to exact USDA strings
    cat_map = {
        'Dairy and Egg Products': 'Dairy and Egg Products',
        'Spices and Herbs': 'Spices and Herbs',
        'Babyfoods': 'Baby Foods',
        'Fats and Oils': 'Fats and Oils',
        'Poultry Products': 'Poultry Products',
        'Soups, Sauces and Gravies': 'Soups, Sauces, and Gravies',
        'Sausages and Luncheon meats': 'Sausages and Luncheon Meats',
        'Breakfast cereals': 'Breakfast Cereals',
        'Fruits and fruit juices': 'Fruits and Fruit Juices',
        'Pork Products': 'Pork Products',
        'Vegetables and Vegetable Products': 'Vegetables and Vegetable Products',
        'Nuts and Seeds': 'Nut and Seed Products',
        'Beef Products': 'Beef Products',
        'Beverages': 'Beverages',
        'Finfish and Shellfish Products': 'Finfish and Shellfish Products',
        'Legumes and Legume Products': 'Legumes and Legume Products',
        'Lamb, Veal and Game': 'Lamb, Veal, and Game Products',
        'Baked Products': 'Baked Products',
        'Sweets': 'Sweets',
        'Cereals, Grains and Pasta': 'Cereal Grains and Pasta',
        'Fast Foods': 'Fast Foods',
        'Meals, Entrees and Side Dishes': 'Meals, Entrees, and Side Dishes',
        'Snacks': 'Snacks',
        'American Indian/Alaska Native Foods': 'American Indian/Alaska Native Foods',
        'Restaurant Foods': 'Restaurant Foods'
    }
    df_food['food_category'] = df_food['FoodGroupName'].map(cat_map).fillna('Meals, Entrees, and Side Dishes')

    print("[3/5] Mapping CNF nutrients to FDC IDs...")
    # Mapping CNF's numeric NutrientIDs to USDA FDC numeric IDs
    target_mapping = {
        203: 1003, # Protein
        204: 1004, # Fat
        205: 1005, # Carbohydrate
        810: 1009, # Starch
        210: 1010, # Sucrose
        211: 1011, # Glucose
        212: 1012, # Fructose
        213: 1013, # Lactose
        214: 1014, # Maltose
        287: 1015, # Galactose
        291: 1079, # Fiber
        269: 2000, # Sugars
        301: 1087, # Calcium
        303: 1089, # Iron
        304: 1090, # Magnesium
        305: 1091, # Phosphorus
        306: 1092, # Potassium
        307: 1093, # Sodium
        309: 1095, # Zinc
        312: 1098, # Copper
        317: 1103, # Selenium
        319: 1105, # Retinol
        321: 1106, # Beta Carotene
        339: 1110, # Vitamin D
        401: 1162, # Vitamin C
        404: 1165, # Thiamin
        405: 1166, # Riboflavin
        406: 1167, # Niacin
        415: 1175, # Vitamin B6
        418: 1174, # Vitamin B12
        815: 1177, # Folate DFE
        806: 1177, # Folate natural (fallback)
        601: 1253, # Cholesterol
        606: 1258, # FA saturated
        645: 1259, # FA mono
        646: 1260, # FA poly
        221: 1018  # Alcohol
    }

    # Filter the massive CNF dataset down to just the nutrients we mapped
    df_amt_filtered = df_amt[df_amt['NutrientID'].isin(target_mapping.keys())].copy()
    
    # Apply the FDC nutrient ID mapping
    df_amt_filtered['fdc_nutrient_id'] = df_amt_filtered['NutrientID'].map(target_mapping)
    
    # Merge with the food dictionary to get the custom fdc_id
    df_amt_filtered = df_amt_filtered.merge(df_food[['FoodID', 'fdc_id']], on='FoodID', how='inner')
    
    df_final = df_amt_filtered[['fdc_id', 'fdc_nutrient_id', 'NutrientValue']].rename(columns={
        'fdc_nutrient_id': 'nutrient_id',
        'NutrientValue': 'amount'
    })
    
    # Clean up zeroes
    df_final['amount'] = pd.to_numeric(df_final['amount'], errors='coerce')
    df_final = df_final.dropna(subset=['amount'])
    df_final = df_final[df_final['amount'] > 0]
    
    # In case multiple CNF codes map to the same FDC code (e.g. multiple Folate types), keep the highest value
    df_final = df_final.groupby(['fdc_id', 'nutrient_id'])['amount'].max().reset_index()

    print("[4/5] Generating database files...")
    # Create the food dictionary injection file
    new_foods = pd.DataFrame({
        'fdc_id': df_food['fdc_id'],
        'data_type': 'foundation_food',
        'description': df_food['FoodDescription'] + " [CNF]",
        'food_category_id': '3333', # Custom category marker
        'food_category': df_food['food_category']
    })
    
    new_foods.to_parquet("cnf_food_injection.parquet", index=False)
    
    print("[5/5] Building bucketed Parquets...")
    df_final['bucket'] = df_final['nutrient_id'] % 256
    
    out_dir = Path("cnf_food_nutrient_bucketed")
    out_dir.mkdir(exist_ok=True)
    
    for bucket_id, group in df_final.groupby('bucket'):
        bucket_dir = out_dir / f"bucket={bucket_id}"
        bucket_dir.mkdir(exist_ok=True)
        group = group.sort_values(['fdc_id', 'nutrient_id'])
        group[['fdc_id', 'nutrient_id', 'amount']].to_parquet(bucket_dir / "cnf_data.parquet", index=False)

    print("\n[SUCCESS] CNF integration files generated!")
    print("1. 'cnf_food_injection.parquet' -> Contains the food metadata.")
    print("2. 'cnf_food_nutrient_bucketed/' -> Contains the partitioned FDC measurements.")

if __name__ == "__main__":
    main()
