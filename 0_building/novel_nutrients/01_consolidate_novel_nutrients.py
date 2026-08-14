import pandas as pd
import glob
import re

def clean_embedded_units(name):
    """
    Strips embedded units and junk from column names.
    e.g., 'Vitamin C (mg 100g)' -> 'Vitamin C'
    e.g., 'Calcium, Ca' -> 'Calcium'
    """
    name = str(name)
    # Remove parentheses containing typical unit markers
    name = re.sub(r'\((mg|g|µg|mcg|ug|%|kcal|kj|g/100g|mg/100g).*?\)', '', name, flags=re.IGNORECASE)
    # Remove brackets
    name = re.sub(r'\[.*?\]', '', name)
    # Remove specific artifacts
    name = re.sub(r'§|µg|mg| g ', '', name)
    
    # Split by comma or 'or' and take the first part (e.g., "Vitamin B1 or Thiamin" -> "Vitamin B1")
    # Actually, keeping it intact is safer for normalization, but we can strip trailing spaces
    return " ".join(name.split()).strip()

def is_valid_chemical(name):
    """Drops internal codes and pure numbers."""
    if not name or len(name) < 3: return False
    
    # Drop AFCD codes like F000050
    if re.match(r'^F\d{5,}$', name, flags=re.IGNORECASE): return False
    
    # Drop purely numeric strings like "100.0" or "101.0"
    if re.match(r'^[\d\.]+$', name): return False
    
    return True

def main():
    print("==========================================")
    print("GLOBAL NOVEL NUTRIENT CONSOLIDATOR")
    print("==========================================")

    # 1. Load Master Normalized Set
    master_df = pd.read_csv("/data/bac2food/global_nutrient.normalized.tsv", sep='\t')
    known_set = set(master_df['name'].dropna().str.lower().str.strip())
    if 'normalized_norm' in master_df.columns:
        known_set.update(master_df['normalized_norm'].dropna().str.lower().str.strip())

    # 2. Gather all Novel CSVs
    files = glob.glob("novel_*.csv")
    if not files:
        print("[!] No novel_*.csv files found in current directory.")
        return

    raw_novel_items = []
    for f in files:
        df = pd.read_csv(f)
        for _, row in df.iterrows():
            raw_novel_items.append({
                "raw_name": row["name"],
                "source": row["source"]
            })

    print(f"[*] Loaded {len(raw_novel_items)} raw candidate molecules.")

    # 3. Clean, Deduplicate, and Validate
    consolidated = {}
    
    for item in raw_novel_items:
        clean_name = clean_embedded_units(item["raw_name"])
        clean_lower = clean_name.lower()
        
        # Check if it's valid
        if not is_valid_chemical(clean_name): continue
            
        # Check if, after unit stripping, it actually already exists in FDC!
        if clean_lower in known_set: continue
            
        # Add to our consolidated dictionary to merge duplicates across countries
        if clean_lower not in consolidated:
            consolidated[clean_lower] = {
                "name": clean_name, # Keep original casing
                "sources": set([item["source"]])
            }
        else:
            consolidated[clean_lower]["sources"].add(item["source"])

    # 4. Mint Final Global IDs
    final_rows = []
    start_id = 96000 # The true starting block for global dark matter
    
    for clean_lower, data in sorted(consolidated.items()):
        # Join sources nicely (e.g., "CIQUAL | NEVO")
        sources_str = " | ".join(sorted(list(data["sources"])))
        
        final_rows.append({
            "id": start_id,
            "name": data["name"],
            "unit_name": "MG", # Assuming MG is the safest default
            "emptyspace": "",
            "fakerank": "63001"
        })
        start_id += 1

    final_df = pd.DataFrame(final_rows)
    final_df.to_csv("MASTER_novel_nutrients.csv", index=False)
    
    print(f"\n[SUCCESS] Distilled down to {len(final_df)} truly unique novel molecules!")
    print("Saved to: MASTER_novel_nutrients.csv")

if __name__ == "__main__":
    main()
