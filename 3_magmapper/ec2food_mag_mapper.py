#!/usr/bin/env python3
"""
ec2food_mag_mapper.py — Metagenomic Food Utilization Predictor
Generates matrices mapping MAGs directly to the foods they can utilize using memory-efficient chunking.
"""
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import re

# --- PATHS ---
DIR_MATRICES = Path("../2_differential/boss_matrices")
PATH_FOOD_NUT_DB = DIR_MATRICES / "Food_Nutrient_DB.zip" # Or .csv if you unzipped it
PATH_NUT_ENZ_DB = DIR_MATRICES / "Nutrient_Enzyme_DB.csv"
#!/usr/bin/env python3

# Regex for validating EC numbers
EC_RE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")

def main():
    ap = argparse.ArgumentParser(description="Map MAGs to Food Databases based on Enzyme overlaps.")
    ap.add_argument("--mag_tsv", required=True, help="Path to your TSV mapping species/MAGs to EC numbers")
    ap.add_argument("--matrix_dir", default="./boss_matrices", help="Directory containing your input databases (Food_Nutrient_DB.zip and Nutrient_Enzyme_DB.csv)")
    ap.add_argument("--out_dir", default="./output_matrices", help="Directory to save the generated output matrices")
    ap.add_argument("--prefix", default="MAG_Food_Utilization", help="Prefix for the generated CSV files")
    args = ap.parse_args()

    print("==========================================================")
    print("Building MAG-to-Food Utilization Matrices (Pathway Count Mode)")
    print("==========================================================")

    # --- SETUP PATHS ---
    matrix_dir = Path(args.matrix_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True) # Create output folder if it doesn't exist

    PATH_FOOD_NUT_DB = matrix_dir / "Food_Nutrient_DB.zip"
    PATH_NUT_ENZ_DB = matrix_dir / "Nutrient_Enzyme_DB.csv"
    
    out_count = out_dir / f"{args.prefix}_ECcount.csv"
    out_bin = out_dir / f"{args.prefix}_binary.csv"

    # =====================================================================
    # 1. LOAD ENZYME DATABASE
    # =====================================================================
    print(f"[*] Loading Enzyme Database: {PATH_NUT_ENZ_DB}...")
    ne_df = pd.read_csv(PATH_NUT_ENZ_DB, index_col="nutrient_key")
    nutrient_columns = list(ne_df.index)

    # =====================================================================
    # 2. PROCESS MAG DATA
    # =====================================================================
    print(f"[*] Parsing MAG genomes from: {args.mag_tsv}...")
    mag_raw = pd.read_csv(args.mag_tsv, sep="\t", dtype="string")
    
    cols = {c.lower(): c for c in mag_raw.columns}
    ec_col = cols.get("ec_number", cols.get("ec", mag_raw.columns[-1]))
    sp_col = cols.get("species", mag_raw.columns[0])
    
    mag_raw = mag_raw.rename(columns={ec_col: "ec_number", sp_col: "species"})
    mag_raw["strain"] = mag_raw.get("strain", pd.Series([""]*len(mag_raw))).fillna("").astype("string").str.strip()
    
    mag_raw = mag_raw[mag_raw["ec_number"].str.match(EC_RE, na=False)].dropna(subset=["species", "ec_number"])
    mag_raw["label"] = mag_raw["species"] + mag_raw["strain"].apply(lambda x: f" | {x}" if x else "")
    
    mag2ec = mag_raw.groupby("label")["ec_number"].apply(set).to_dict()
    mags = list(mag2ec.keys())
    print(f"  -> Discovered {len(mags)} unique MAGs.")

    # =====================================================================
    # 3. BUILD MAG-NUTRIENT CAPABILITY MATRIX (M_matrix)
    # =====================================================================
    print("[*] Translating Genomes into Metabolic Capabilities...")
    M_matrix = pd.DataFrame(0, index=mags, columns=nutrient_columns, dtype=np.int8)
    
    valid_db_ecs = set(ne_df.columns)
    for mag, ecs in mag2ec.items():
        usable_ecs = list(ecs.intersection(valid_db_ecs))
        if usable_ecs:
            M_matrix.loc[mag] = ne_df[usable_ecs].max(axis=1).astype(np.int8)
            
    M_matrix_T = M_matrix.T.astype(np.int32)

    # =====================================================================
    # 4. CHUNKED MATRIX MULTIPLICATION
    # =====================================================================
    if out_count.exists(): out_count.unlink()
    if out_bin.exists(): out_bin.unlink()

    print(f"[*] Streaming Food Database & Calculating Pathway Matches...")
    
    chunk_size = 100_000
    first_chunk = True
    compression_type = 'zip' if str(PATH_FOOD_NUT_DB).endswith('.zip') else 'infer'
    
    for i, chunk in enumerate(pd.read_csv(PATH_FOOD_NUT_DB, compression=compression_type, chunksize=chunk_size)):
        print(f"    -> Processing rows {i * chunk_size:,} to {(i + 1) * chunk_size:,}...")
        
        # 1. Isolate the data
        food_meta = chunk[['fdc_id', 'food_name']]
        
        # 2. Binarize the food array
        F_matrix_bin = (chunk[nutrient_columns] > 0).astype(np.int32)
        
        # 3. Linear Algebra
        match_matrix = F_matrix_bin.dot(M_matrix_T)
        
        # 4. Format outputs
        final_count = pd.concat([food_meta, match_matrix], axis=1)
        final_binary = final_count.copy()
        final_binary[mags] = (match_matrix > 0).astype(np.int8)
        
        # 5. Write to disk
        mode = 'w' if first_chunk else 'a'
        header = first_chunk
        
        final_count.to_csv(out_count, mode=mode, header=header, index=False)
        final_binary.to_csv(out_bin, mode=mode, header=header, index=False)
        
        first_chunk = False

    print(f"\n[*] DONE! Saved outputs to {out_dir.absolute()}")
    print(f"    -> {out_count.name}")
    print(f"    -> {out_bin.name}")

if __name__ == "__main__":
    main()