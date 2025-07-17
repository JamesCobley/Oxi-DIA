import pandas as pd
import numpy as np
from scipy.special import rel_entr
from scipy.stats import ttest_rel

# === Load dataset ===
df = pd.read_csv('/content/redox_sites.tsv', sep='\t')

# === Define matched condition columns ===
air_cols = ['Sample_1_%Oxidized', 'Sample_2_%Oxidized', 'Sample_3_%Oxidized']
tin_cols = ['Sample_4_%Oxidized', 'Sample_5_%Oxidized', 'Sample_6_%Oxidized']

# === Function to compute probability distribution from values ===
def get_distribution(values, bins=50):
    values = values.dropna().values
    hist, _ = np.histogram(values, bins=bins, range=(0, 100), density=False)
    prob = hist / hist.sum()
    return prob

# === Compute KL divergence for each paired replicate ===
kl_values = []
for air_col, tin_col in zip(air_cols, tin_cols):
    P = get_distribution(df[air_col])
    Q = get_distribution(df[tin_col])
    mask = (P > 0) & (Q > 0)
    kl = np.sum(rel_entr(P[mask], Q[mask]))
    kl_values.append(kl)
    print(f"KL({air_col} || {tin_col}): {kl:.4f}")

# === Summary statistics ===
kl_values = np.array(kl_values)
mean_kl = kl_values.mean()
std_kl = kl_values.std(ddof=1)

# === Paired t-test on KL values vs zero change hypothesis ===
t_stat, p_val = ttest_rel(kl_values, np.zeros_like(kl_values))
cohen_d = mean_kl / std_kl

# === Report summary ===
print("\n=== KL Divergence Summary ===")
print(f"Mean KL Divergence (Air || Tin): {mean_kl:.4f}")
print(f"Standard Deviation: {std_kl:.4f}")
print(f"t = {t_stat:.4f}, p = {p_val:.4f}")
print(f"Cohen's d: {cohen_d:.4f}")
