import pandas as pd
import numpy as np
from scipy.stats import ttest_rel
from numpy import sqrt, sum, arccos, clip

# === Step 1: Load the data ===
df = pd.read_csv('/content/redox_sites.tsv', sep='\t')

# === Step 2: Define matched sample columns ===
air_cols = ['Sample_1_%Oxidized', 'Sample_2_%Oxidized', 'Sample_3_%Oxidized']
tin_cols = ['Sample_4_%Oxidized', 'Sample_5_%Oxidized', 'Sample_6_%Oxidized']

# === Step 3: Fisher–Rao distance via Hellinger inner product ===
def fisher_rao_distance(x, y, bins=50):
    mask = ~x.isna() & ~y.isna()
    x_vals = x[mask].values
    y_vals = y[mask].values

    hist_x, _ = np.histogram(x_vals, bins=bins, range=(0, 100), density=True)
    hist_y, _ = np.histogram(y_vals, bins=bins, range=(0, 100), density=True)

    P = hist_x / hist_x.sum()
    Q = hist_y / hist_y.sum()

    # Compute the Bhattacharyya coefficient
    inner_product = np.sum(np.sqrt(P * Q))

    # Clip to avoid numerical issues
    inner_product = clip(inner_product, 0, 1)

    # Fisher–Rao distance
    return 2 * arccos(inner_product)

# === Step 4: Compute Fisher–Rao distance for each matched pair ===
fr_values = []
for air_col, tin_col in zip(air_cols, tin_cols):
    fr = fisher_rao_distance(df[air_col], df[tin_col])
    fr_values.append(fr)
    print(f"Fisher–Rao({air_col} vs {tin_col}): {fr:.4f}")

# === Step 5: Stats ===
t_stat, p_val = ttest_rel(fr_values, np.zeros_like(fr_values))
mean_fr = np.mean(fr_values)
std_fr = np.std(fr_values, ddof=1)
cohen_d = mean_fr / std_fr

# === Step 6: Report summary ===
print("\n=== Fisher–Rao Distance Summary ===")
print(f"Mean Fisher–Rao: {mean_fr:.4f}")
print(f"Standard Deviation: {std_fr:.4f}")
print(f"t = {t_stat:.4f}, p = {p_val:.4f}")
print(f"Cohen's d: {cohen_d:.4f}")
