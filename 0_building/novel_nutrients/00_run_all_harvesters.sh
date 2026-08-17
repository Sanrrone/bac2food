#!/bin/bash

# Array containing the exact paths to your raw database files
RAW_FILES=(
  "../food_DBs/asnut_australianw/AFCD Release 3 - Nutrient profiles.xlsx"
  "../food_DBs/canadian_nutrientfile_canada/NUTRIENT NAME.csv"
  "../food_DBs/ciqual_france/Table Ciqual 2025_ENG_2025_11_03.xlsx"
  "../food_DBs/FAO_onu/BioFoodComp4.0.xlsx"
  "../food_DBs/FAO_onu/PhyFoodComp_1.0.xlsx"
  "../food_DBs/FAO_onu/WAFCT_2019.xlsx"
  "../food_DBs/fineli_finland/resultset.csv"
  "../food_DBs/frida_denmark/Frida_5.5_Dataset.ods"
  "../food_DBs/livsmedels_sweeden/LivsmedelsDB_202603061604.xlsx"
  "../food_DBs/McCance_Widdowsons_uk/McCance_Widdowsons_Composition_of_Foods_Integrated_Dataset_2021..xlsx"
  "../food_DBs/phenol_explorer_france/composition-data.tsv"
  "../food_DBs/stfcj_japan/main_1374049_1r12_1.xlsx"
  "../food_DBs/stfcj_japan/org_acid_1388558_4r12r.xlsx"
  "../food_DBs/swissfoodcompoDB_swiss/Swiss_food_composition_database.xlsx"
)

# Start assigning new IDs from 96000
CURRENT_ID=96000

echo "Starting Global Dark Matter Harvest..."

for file in "${RAW_FILES[@]}"; do
    echo "=================================================="
    echo "Harvesting: $file"
    echo "Starting ID block: $CURRENT_ID"
    echo "=================================================="
    
    # Run the Python harvester
    python 99_omnivorous_harvester.py \
      --raw_db "$file" \
      --start_id $CURRENT_ID \
      --normalized_tsv "/data/bac2food/global_nutrient.normalized.tsv"
    
    # Jump the ID block by 1,000 for the next database to guarantee no overlap
    # (e.g., CIQUAL gets 96000-96999, Canada gets 97000-97999, etc.)
    CURRENT_ID=$((CURRENT_ID + 1000))
    
    echo ""
done

echo "[*] Harvesting Complete. Check your current directory for the 'novel_*.csv' files!"
