import pandas as pd
import numpy as np
from scipy.stats import ttest_rel

# === Load your filtered dataframe with common sites ===
df = pd.read_csv('/content/redox_sites.tsv', sep='\t')

# === Define correct columns for Air and Tin ===
air_cols = ['Sample_1_%Oxidized', 'Sample_2_%Oxidized', 'Sample_3_%Oxidized']
tin_cols = ['Sample_4_%Oxidized', 'Sample_5_%Oxidized', 'Sample_6_%Oxidized']

# === Shannon entropy function ===
def shannon_entropy(values, bins=50):
    values = values.dropna()
    hist, _ = np.histogram(values, bins=bins, range=(0, 100), density=True)
    hist = hist[hist > 0]
    return -np.sum(hist * np.log2(hist))

# === Compute entropy for each run ===
entropy_dict = {}
air_entropy = []
tin_entropy = []

for col in air_cols + tin_cols:
    entropy = shannon_entropy(df[col])
    entropy_dict[col] = entropy
    if col in air_cols:
        air_entropy.append(entropy)
    else:
        tin_entropy.append(entropy)

# === Paired t-test ===
air_entropy = np.array(air_entropy)
tin_entropy = np.array(tin_entropy)
t_stat, p_val = ttest_rel(air_entropy, tin_entropy)
mean_diff = np.mean(air_entropy - tin_entropy)
cohen_d = mean_diff / np.std(air_entropy - tin_entropy, ddof=1)

# === Create summary dataframe ===
entropy_df = pd.DataFrame({
    'Sample': air_cols + tin_cols,
    'Condition': ['Air'] * 3 + ['Tin'] * 3,
    'Shannon_Entropy': air_entropy.tolist() + tin_entropy.tolist()
})

# === Print results ===
print("=== Shannon Entropy per Sample ===")
print(entropy_df)

print("\n=== Paired t-test (Shannon Entropy: Air vs Tin) ===")
print(f"t = {t_stat:.4f}, p = {p_val:.4f}")
print(f"Mean difference (Air - Tin): {mean_diff:.4f}")
print(f"Cohen's d: {cohen_d:.4f}")
