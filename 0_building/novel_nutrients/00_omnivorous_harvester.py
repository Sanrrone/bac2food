import pandas as pd
import re
import os
import argparse
from pathlib import Path

# --- Multi-lingual Metadata Stopwords ---
IGNORE_TERMS = [
    "food", "index", "code", "item", "name", "desc", "remark", "unit", "weight", "factor", 
    "energy", "kcal", "kj", "joule", "refuse", "edible", "portion", "water", "moisture",
    "protein", "fat", "lipid", "carbohydrate", "ash", "total", "sum", "equivalent", "error", 
    "source", "reference", "date", "id", "group", "category", "type", "method", "sample",
    "aliment", "groupe", "energie", "eau", "proteines", "lipides", "glucides", "cendres",
    "voedingsmiddel", "naam", "eiwit", "koolhydraten", "vet", "vand", "kulhydrat", "energia", "vesi",
    "unnamed"
]

def clean_term(c):
    c = str(c).replace('\n', ' ').replace('\r', ' ')
    return " ".join(c.split())

def is_novel_chemical(term, existing_names):
    term_clean = clean_term(term).lower()
    
    if len(term_clean) < 3 or term_clean.isnumeric(): return False
    for stop in IGNORE_TERMS:
        if stop in term_clean: return False
    if term_clean in existing_names: return False
    if len(term_clean.split()) > 6: return False
        
    return True

def extract_terms_from_df(df):
    """Extracts terms from headers (Wide format) AND specific columns (Long format)."""
    terms = set()
    
    # 1. Grab all column headers (Wide format)
    for col in df.columns:
        terms.add(str(col))
        
    # 2. Grab values if it's a Long format database (e.g., CNF, Phenol-Explorer)
    long_format_keywords = ['nutrient', 'component', 'parameter', 'compound']
    for col in df.columns:
        if any(k in str(col).lower() for k in long_format_keywords):
            # Add all unique string values from this column
            unique_vals = df[col].dropna().astype(str).unique()
            terms.update(unique_vals)
            
    return terms

def load_file(filepath):
    """Dynamically loads the file based on its extension."""
    ext = Path(filepath).suffix.lower()
    
    try:
        if ext == '.csv':
            # Try comma, fallback to semicolon
            try: return [pd.read_csv(filepath, nrows=1000, low_memory=False)]
            except: return [pd.read_csv(filepath, sep=';', nrows=1000, low_memory=False)]
        elif ext == '.tsv':
            return [pd.read_csv(filepath, sep='\t', nrows=1000, low_memory=False)]
        elif ext in ['.xlsx', '.xls', '.ods']:
            # Read all sheets into a dict of DataFrames
            engine = 'odf' if ext == '.ods' else None # Pandas auto-detects xls/xlsx usually
            dfs = pd.read_excel(filepath, sheet_name=None, nrows=1000, engine=engine)
            return list(dfs.values())
        else:
            print(f"[!] Unsupported format: {ext}")
            return []
    except Exception as e:
        print(f"[!] Error reading {filepath}: {e}")
        return []

def main():
    parser = argparse.ArgumentParser(description="Omnivorous Dark Matter Harvester")
    parser.add_argument("--raw_db", required=True, help="Path to ANY database file (.csv, .tsv, .xlsx, .xls, .ods)")
    parser.add_argument("--normalized_tsv", default="/data/bac2food/global_nutrient.normalized.tsv")
    parser.add_argument("--start_id", type=int, default=96000)
    args = parser.parse_args()

    print(f"[*] Loading Master Ontology...")
    norm_df = pd.read_csv(args.normalized_tsv, sep='\t')
    known_set = set(norm_df['name'].dropna().str.lower().str.strip())
    if 'normalized_norm' in norm_df.columns:
        known_set.update(norm_df['normalized_norm'].dropna().str.lower().str.strip())

    print(f"[*] Scanning {args.raw_db}...")
    dfs = load_file(args.raw_db)
    
    raw_terms = set()
    for df in dfs:
        raw_terms.update(extract_terms_from_df(df))
        
    novel_nutrients = []
    current_id = args.start_id
    
    for term in sorted(raw_terms):
        if is_novel_chemical(term, known_set):
            novel_nutrients.append({
                "id": current_id,
                "name": clean_term(term),
                "unit_name": "MG",
                "source": os.path.basename(args.raw_db)
            })
            current_id += 1

    if not novel_nutrients:
        print("[!] No novel biochemicals found. Everything is mapped or filtered!")
        return

    out_df = pd.DataFrame(novel_nutrients)
    out_csv = f"novel_{os.path.basename(args.raw_db)}.csv"
    out_df.to_csv(out_csv, index=False)
    
    print(f"\n[SUCCESS] Found {len(out_df)} potentially novel biochemicals! Saved to {out_csv}")
    print("\n--- DICTIONARY FOR YOUR INGESTION SCRIPT ---")
    for _, row in out_df.iterrows():
        print(f'    "{row["name"]}": {row["id"]},')

if __name__ == "__main__":
    main()
