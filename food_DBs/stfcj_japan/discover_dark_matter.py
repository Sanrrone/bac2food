import pandas as pd
import re

# The columns we ALREADY mapped (so we don't duplicate them)
MAPPED_JAPAN = [
    "Water", "Protein, calculated from reference nitrogen", "Lipid", "Ash", "Carbohydrate, available, total", 
    "Fiber, total dietary", "Starch", "Sucrose", "Glucose", "Fructose", "Lactose", "Maltose", "Galactose", 
    "Sorbitol", "Mannitol", "Trehalose", "Acetic acid", "Lactic acid", "Citric acid", "Malic acid", 
    "Quinic acid", "Succinic acid", "Tartaric acid", "Chlorogenic acid", "Caffeic acid", "Ferulic acid", 
    "p-Coumaric acid", "Isoleucine", "Leucine", "Lysine", "Methionine", "Cystine", "Phenylalanine", 
    "Tyrosine", "Threonine", "Tryptophan", "Valine", "Histidine", "Arginine", "Alanine", "Aspartic acid", 
    "Glutamic acid", "Glycine", "Proline", "Serine", "Hydroxy-proline", "Item No."
]

# Words that mean it's metadata or a macro-summary, not a specific chemical
IGNORE_WORDS = ["Food", "Index", "Remarks", "Energy", "Yield", "Refuse", "Total", "sum", "equivalent", "Item", "Code"]

def clean_col(c):
    return " ".join(str(c).split())

def is_chemical(col_name):
    # Ignore if it's already mapped
    for m in MAPPED_JAPAN:
        if m.lower() in col_name.lower(): return False
    # Ignore metadata
    for w in IGNORE_WORDS:
        if w.lower() in col_name.lower(): return False
    # Ignore empty or unnamed
    if not col_name or "Unnamed" in col_name: return False
    return True

print("Scanning Japanese STFCJ for unmapped molecules...")
df_org = pd.read_excel("org_acid_1388558_4r12r.xlsx", sheet_name="Annex (organic acids)", header=4)
df_main = pd.read_excel("main_1374049_1r12_1.xlsx", sheet_name="Table", header=5)

new_nutrients = []
start_id = 95000

# Scan Organic Acids
for col in df_org.columns:
    c = clean_col(col)
    if is_chemical(c):
        new_nutrients.append({"id": start_id, "name": c, "unit_name": "MG", "source": "STFCJ_OrgAcids"})
        start_id += 1

# Scan Main Table (looking for polyphenols, tannins, etc.)
for col in df_main.columns:
    c = clean_col(col)
    if is_chemical(c):
        new_nutrients.append({"id": start_id, "name": c, "unit_name": "MG", "source": "STFCJ_Main"})
        start_id += 1

dark_matter_df = pd.DataFrame(new_nutrients)
dark_matter_df.to_csv("proposed_new_nutrients.csv", index=False)
print(f"Found {len(dark_matter_df)} unmapped biochemicals! Saved to proposed_new_nutrients.csv")
