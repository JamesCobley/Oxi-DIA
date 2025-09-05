import pandas as pd
import numpy as np

# === Load the fully complete matrix you already created ===
atlas_in = "/content/fully_complete_redox_matrix.xlsx"
df = pd.read_excel(atlas_in)

# Identify the biological sample columns (everything that's not an ID column)
id_cols = ["Site_Key", "Protein.Ids"]
sample_cols = [c for c in df.columns if c not in id_cols]

# Sanity check: ensure numeric
df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors="coerce")

# === Compute per-site reference stats ===
ref = df[id_cols].copy()
ref["Mean_Oxidation"] = df[sample_cols].mean(axis=1)           # mean across 20 samples
ref["SD_Oxidation"]   = df[sample_cols].std(axis=1, ddof=1)    # sample SD

# Optional: quick flags to make the atlas actionable
# (tweak thresholds to taste)
ref["Is_Fully_Reduced_Majority"] = (df[sample_cols] == 0).mean(axis=1) >= 0.5
ref["Is_Stable_SD_lt_2.5pct"]    = ref["SD_Oxidation"] < 2.5
ref["Is_Mean_gt_5pct"]           = ref["Mean_Oxidation"] > 5

# === Save the atlas ===
atlas_out = "/content/reference_redox_atlas.xlsx"
ref.to_excel(atlas_out, index=False)

print(f"✅ Reference redox atlas written to: {atlas_out}")
print(ref.head(5))
