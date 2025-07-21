import pandas as pd
import numpy as np
from sklearn.metrics import mutual_info_score
from scipy.stats import ttest_rel

# === Step 1: Load the data ===
df = pd.read_csv('/content/redox_sites.tsv', sep='\t')

# === Step 2: Define matched sample columns ===
air_cols = ['Sample_1_%Oxidized', 'Sample_2_%Oxidized', 'Sample_3_%Oxidized']
tin_cols = ['Sample_4_%Oxidized', 'Sample_5_%Oxidized', 'Sample_6_%Oxidized']

# === Step 3: Define MI calculation with discretization ===
def compute_mi(x, y, bins=50):
    # Drop missing values at corresponding indices
    mask = ~x.isna() & ~y.isna()
    x_binned = pd.cut(x[mask], bins=bins, labels=False)
    y_binned = pd.cut(y[mask], bins=bins, labels=False)
    return mutual_info_score(x_binned, y_binned)

# === Step 4: Compute MI for each matched pair ===
mi_values = []
for air_col, tin_col in zip(air_cols, tin_cols):
    mi = compute_mi(df[air_col], df[tin_col])
    mi_values.append(mi)
    print(f"MI({air_col} vs {tin_col}): {mi:.4f}")

# === Step 5: Perform one-sample test vs 0 (null: no mutual structure) ===
t_stat, p_val = ttest_rel(mi_values, np.zeros_like(mi_values))
mean_mi = np.mean(mi_values)
std_mi = np.std(mi_values, ddof=1)
cohen_d = mean_mi / std_mi

# === Step 6: Report summary ===
print("\n=== Mutual Information Summary ===")
print(f"Mean MI: {mean_mi:.4f}")
print(f"Standard Deviation: {std_mi:.4f}")
print(f"t = {t_stat:.4f}, p = {p_val:.4f}")
print(f"Cohen's d: {cohen_d:.4f}")
