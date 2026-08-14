import pandas as pd

def clean_text(text):
    if not isinstance(text, str):
        return text
    # Keep only printable ASCII characters (32-126)
    cleaned = "".join(char for char in text if 32 <= ord(char) <= 126)
    return cleaned.strip()

input_file = 'ec_species_substrate.tsv'
output_file = 'cleaned_unique_ec_substrates.tsv'
chunk_size = 50000  # Adjust based on your RAM (smaller = safer)

# Using a set to track unique pairs across chunks without keeping the whole DF in memory
unique_pairs_set = set()

print(f"Starting chunked processing of {input_file}...")

# 1. Process in chunks
chunks = pd.read_csv(input_file, sep='\t', on_bad_lines='skip', chunksize=chunk_size, encoding='utf-8')

for i, chunk in enumerate(chunks):
    # Clean only the columns we need (saves RAM)
    chunk['ec'] = chunk['ec'].astype(str).map(clean_text)
    chunk['substrate'] = chunk['substrate'].astype(str).map(clean_text)
    
    # Drop empty values and ensure EC looks like a number
    chunk = chunk.dropna(subset=['ec', 'substrate'])
    chunk = chunk[chunk['ec'].str.contains(r'^\d', na=False)]
    
    # Add pairs to our set (sets automatically handle duplicates)
    for row in chunk[['ec', 'substrate']].itertuples(index=False):
        unique_pairs_set.add(row)
    
    if i % 10 == 0:
        print(f"Processed { (i+1) * chunk_size } rows...")

# 2. Convert the set back to a DataFrame and save
print("Saving unique pairs...")
df_final = pd.DataFrame(list(unique_pairs_set), columns=['ec', 'substrate'])
df_final.sort_values(by=['ec', 'substrate']).to_csv(output_file, sep='\t', index=False)

print(f"Done! Final unique pairs saved to {output_file}")
