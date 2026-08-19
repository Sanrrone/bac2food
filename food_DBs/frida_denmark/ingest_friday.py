#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq

# fdc_id allocation: fdc_blocks.py
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
import fdc_blocks


def main():
    print("==================================================")
    print("FRIDA (DENMARK) TO FDC INGESTION SCRIPT")
    print("==================================================")
    
    excel_file = "Frida_5.5_Dataset.xlsx"
    
    print("[1/5] Loading Frida relational data from Excel tabs...")
    # Read directly from the named sheets in the Excel file
    df_norm = pd.read_excel(excel_file, sheet_name="Data_Normalised")
    df_food = pd.read_excel(excel_file, sheet_name="Food")
    df_param = pd.read_excel(excel_file, sheet_name="Parameter")
    
    print("[2/5] Minting custom FDC IDs and mapping categories...")
    # Offset Frida FoodIDs by 50,000,000
    df_food['fdc_id'] = fdc_blocks.assign('frida', df_food['FoodID'].astype(int))
    
    # Map Frida's English "FoodGroup" to FDC Categories
    cat_map = {
        'Vegetables': 'Vegetables and Vegetable Products',
        'Potatoes and other tuberous vegetables': 'Vegetables and Vegetable Products',
        'Legumes': 'Legumes and Legume Products',
        'Fruit and fruit products': 'Fruits and Fruit Juices',
        'Nuts and seeds': 'Nut and Seed Products',
        'Cereals and cereal products': 'Cereal Grains and Pasta',
        'Milk and milk products': 'Dairy and Egg Products',
        'Cheese': 'Dairy and Egg Products',
        'Meat and meat products': 'Sausages and Luncheon Meats',
        'Poultry and poultry products': 'Poultry Products',
        'Fish and fish products': 'Finfish and Shellfish Products',
        'Fats and oils': 'Fats and Oils',
        'Beverages': 'Beverages',
        'Sugar and sugar products': 'Sweets',
        'Spices and herbs': 'Spices and Herbs',
        'Snacks': 'Snacks',
        'Eggs and egg products': 'Dairy and Egg Products',
        'Juice': 'Fruits and Fruit Juices'
    }
    df_food['food_category'] = df_food['FoodGroup'].map(cat_map).fillna('Meals, Entrees, and Side Dishes')
    
    print("[3/5] Mapping Frida nutrients to FDC IDs...")
    # Map the English ParameterName from Frida to the USDA FDC numeric ID
    target_mapping = {
        'Protein': 1003,
        'Fat': 1004,
        'Carbohydrate, available': 1005,
        'Starch': 1009,
        'Sugars, total': 2000,
        'Added sugar': 2000, # Fallback for sugars
        'Glucose': 1011,
        'Fructose': 1012,
        'Lactose': 1013,
        'Maltose': 1014,
        'Sucrose': 1010,
        'Dietary fibre': 1079,
        'Alcohol': 1018,
        'Calcium, Ca': 1087,
        'Iron, Fe': 1089,
        'Magnesium, Mg': 1090,
        'Phosphorus, P': 1091,
        'Potassium, K': 1092,
        'Sodium, Na': 1093,
        'Zinc, Zn': 1095,
        'Copper, Cu': 1098,
        'Selenium, Se': 1103,
        'Iodine, I': 1100,
        'Vitamin A': 1104,
        'Retinol': 1105,
        'ß-carotene': 1106,
        'Vitamin D': 1110,
        'Vitamin E': 1109,
        'Vitamin K': 1185,
        'Thiamin': 1165,
        'Riboflavin': 1166,
        'Niacin': 1167,
        'Vitamin B6': 1175,
        'Vitamin B12': 1174,
        'Folate': 1177,
        'Vitamin C': 1162,
        'Cholesterol': 1253,
        'Fatty acids, sum, saturated': 1258,
        'Fatty acids, sum, monounsaturated': 1259,
        'Fatty acids, sum, polyunsaturated': 1260
    }
    
    # Add the FDC target ID to the parameter dataframe
    df_param['nutrient_id'] = df_param['ParameterName'].map(target_mapping)
    
    # Filter only the parameters we successfully mapped
    mapped_params = df_param.dropna(subset=['nutrient_id'])[['ParameterID', 'nutrient_id']].copy()
    mapped_params['nutrient_id'] = mapped_params['nutrient_id'].astype(int)
    
    print("[4/5] Merging and cleaning nutritional values...")
    # Merge the normalized data with our mapped parameters
    df_merged = df_norm.merge(mapped_params, on='ParameterID', how='inner')
    
    # Merge with the food dictionary to get our custom fdc_id
    df_merged = df_merged.merge(df_food[['FoodID', 'fdc_id']], on='FoodID', how='inner')
    
    # Keep only the columns we need: fdc_id, nutrient_id, amount (ResVal in Frida)
    df_final = df_merged[['fdc_id', 'nutrient_id', 'ResVal']].rename(columns={'ResVal': 'amount'})
    
    # Clean up any bad numeric data
    df_final['amount'] = pd.to_numeric(df_final['amount'], errors='coerce')
    df_final = df_final.dropna(subset=['amount'])
    df_final = df_final[df_final['amount'] > 0]
    
    print("[5/5] Generating database files...")
    # Create the food dictionary injection file
    new_foods = pd.DataFrame({
        'fdc_id': df_food['fdc_id'],
        'data_type': 'foundation_food',
        'description': df_food['FoodName'] + " [Frida]",
        'food_category_id': '4444', # Custom category ID for Denmark
        'food_category': df_food['food_category']
    })
    new_foods.to_parquet("frida_food_injection.parquet", index=False)
    
    # Build bucketed Parquets
    df_final['bucket'] = df_final['nutrient_id'] % 256
    
    out_dir = Path("frida_food_nutrient_bucketed")
    out_dir.mkdir(exist_ok=True)
    
    for bucket_id, group in df_final.groupby('bucket'):
        bucket_dir = out_dir / f"bucket={bucket_id}"
        bucket_dir.mkdir(exist_ok=True)
        group = group.sort_values(['fdc_id', 'nutrient_id'])
        group[['fdc_id', 'nutrient_id', 'amount']].to_parquet(bucket_dir / "frida_data.parquet", index=False)

    print("\n[SUCCESS] Frida integration files generated!")
    print("1. 'frida_food_injection.parquet' -> Contains the food metadata.")
    print("2. 'frida_food_nutrient_bucketed/' -> Contains the partitioned FDC measurements.")

if __name__ == "__main__":
    main()
