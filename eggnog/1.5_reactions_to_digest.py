import argparse
import re

import pandas as pd

def normalize_for_chebi(text):
    """
    Transforms chemical names into a format more likely to match ChEBI records.
    """
    if not isinstance(text, str):
        return ""
    
    # 1. Basic character cleaning (Printable ASCII only)
    text = "".join(char for char in text if 32 <= ord(char) <= 126)
    
    # 2. Remove common stereoisomer prefixes (L-, D-, (S)-, etc.)
    # ChEBI often has separate entries, but normalization helps initial mapping.
    text = re.sub(r'\b([LD]-|\([RS]\)-|\(z\)-|\(e\)-)', '', text, flags=re.IGNORECASE)
    
    # 3. Remove hydration and state annotations
    # ChEBI names rarely include "(ion)" or "(aqueous)" in the primary label.
    text = re.sub(r'\s*\((ion|aqueous|solution|reduced|oxidized|re-face|si-face)\)', '', text, flags=re.IGNORECASE)
    
    # 4. Collapse whitespace and lowercase
    text = re.sub(r'\s+', ' ', text).strip().lower()
    
    return text

def get_ubiquitous_set():
    """Bookkeeping molecules that should be flagged or dropped."""
    return {
        "h2o", "water", "h+", "proton", "hydron", "o2", "oxygen", "co2",
        "carbon dioxide", "nad+","nad(+)","nad(-)", "nadh", "nadp+", "nadph", "atp", "adp",
        "coa", "coenzyme a", "pi", "ppi", "nh3", "ammonia","nadp(+)","nadp(-)",
        "more","fadh2","fadh","fad","utp","gtp","datp","dutp",
        # BRENDA writes the ambiguous-specificity forms too, and they are every bit as
        # ubiquitous as the resolved ones. Absent until 2026-08-06, when the official
        # download replaced the SPARQL scrape and surfaced 481 rows of them.
        "nad(p)+","nad(p)h","nad(p)(+)",
        # Generic placeholders for an unnamed redox partner. Same class as "more".
        "acceptor","reduced acceptor","donor","reduced donor",
    }


# BRENDA/Rhea wrap MACROMOLECULAR entities in square brackets -- "[oxidized
# NADPH-hemoprotein reductase]", "[thioredoxin]". They are proteins, not dietary compounds,
# and nothing in the food tables can ever match one. 186 such rows were in the old digest.
MACROMOLECULE = re.compile(r"^\[.*\]$")

_ap = argparse.ArgumentParser(description="Normalize substrate names and drop ubiquitous cofactors.")
_ap.add_argument("--input", default="2_digest_dict.tsv")
_ap.add_argument("--output", default="2_digest_norm.tsv")
_args = _ap.parse_args()

input_file = _args.input
output_file = _args.output
ubiquitous = get_ubiquitous_set()

print("Starting normalization...")

# Chunked processing to avoid 'Killed' memory error
reader = pd.read_csv(input_file, sep='\t', chunksize=100000, on_bad_lines='skip')
first_chunk = True

for i, chunk in enumerate(reader):
    # Keep original by creating a new column
    chunk['substrate_normalized'] = chunk['substrate'].map(normalize_for_chebi)
    
    # Filter: Remove rows where the normalized name is ubiquitous
    # (keeps your final file relevant for food/metabolism research)
    chunk = chunk[~chunk['substrate_normalized'].isin(ubiquitous)]
    
    # Filter: Remove short junk (like single letters or '+')
    chunk = chunk[chunk['substrate_normalized'].str.len() > 2]

    # Drop bracketed macromolecules (see MACROMOLECULE above).
    chunk = chunk[~chunk['substrate_normalized'].str.match(MACROMOLECULE, na=False)]

    # Save
    mode = 'w' if first_chunk else 'a'
    chunk.to_csv(output_file, sep='\t', index=False, mode=mode, header=first_chunk)
    first_chunk = False
    
    if i % 10 == 0:
        print(f"Processed {(i+1) * 100000} rows...")

print(f"Done. File saved as {output_file}")
