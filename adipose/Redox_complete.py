import pandas as pd
import numpy as np
import re

# Load the data
df = pd.read_excel("/content/cys_summary_with_sites.xlsx", sheet_name="Redox Site Summary")

# Step 1: Identify all replicate columns
sample_cols = [col for col in df.columns if col.startswith("James_Twin")]

# Step 2: Group replicates by biological sample number (e.g., James_Twin_83)
bio_sample_map = {}
for col in sample_cols:
    match = re.match(r'(James_Twin_\d+)_S\d+', col)
    if match:
        key = match.group(1)
        bio_sample_map.setdefault(key, []).append(col)

# Step 3: Create a new DataFrame with merged biological replicates (mean of available replicates)
merged_df = df.copy()
for bio_sample, reps in bio_sample_map.items():
    merged_df[bio_sample] = df[reps].mean(axis=1, skipna=True)

# Step 4: Create list of new biological replicate columns
bio_sample_cols = list(bio_sample_map.keys())

# Step 5: Count unique sites and proteins *per biological sample* (presence in at least one replicate)
unique_sites_per_bio_sample = {}
unique_proteins_per_bio_sample = {}

for col in bio_sample_cols:
    detected = merged_df[col].notna()
    unique_sites_per_bio_sample[col] = merged_df.loc[detected, "Site_Key"].nunique()
    unique_proteins_per_bio_sample[col] = merged_df.loc[detected, "Protein.Ids"].nunique()

# Step 6: Filter to rows with no missing values across all 20 merged samples
complete_rows = merged_df[bio_sample_cols].notna().all(axis=1)
fully_quant_sites = merged_df.loc[complete_rows, "Site_Key"].nunique()
fully_quant_proteins = merged_df.loc[complete_rows, "Protein.Ids"].nunique()

# ✅ Final output
print("✅ Unique cysteine sites per biological sample:\n", unique_sites_per_bio_sample)
print("\n✅ Unique proteins per biological sample:\n", unique_proteins_per_bio_sample)
print(f"\n✅ 100% complete cysteine sites across all biological samples: {fully_quant_sites}")
print(f"✅ 100% complete proteins across all biological samples: {fully_quant_proteins}")

# Optional: Save the fully complete data matrix
merged_complete_df = merged_df.loc[complete_rows, ["Site_Key", "Protein.Ids"] + bio_sample_cols]
merged_complete_df.to_excel("/content/fully_complete_redox_matrix.xlsx", index=False)
