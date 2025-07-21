import pandas as pd
import numpy as np
from scipy.stats import ttest_rel

# === Step 1: Load data ===
df = pd.read_csv('/content/redox_sites.tsv', sep='\t')

# === Step 2: Define matched sample columns ===
air_cols = ['Sample_1_%Oxidized', 'Sample_2_%Oxidized', 'Sample_3_%Oxidized']
tin_cols = ['Sample_4_%Oxidized', 'Sample_5_%Oxidized', 'Sample_6_%Oxidized']

# === Step 3: Function to compute FIM between two distributions ===
def compute_fim(x, y, bins=50):
    # Drop NaNs
    mask = ~x.isna() & ~y.isna()
    x_vals = x[mask].values
    y_vals = y[mask].values

    # Bin values to estimate discrete distributions
    hist_x, _ = np.histogram(x_vals, bins=bins, range=(0, 100), density=True)
    hist_y, _ = np.histogram(y_vals, bins=bins, range=(0, 100), density=True)

    # Normalize to probability distributions
    P = hist_x / np.sum(hist_x)
    Q = hist_y / np.sum(hist_y)

    # Avoid division by zero or log of zero
    mask = (P > 0) & (Q > 0)

    # Fisher Information Metric (symmetric discrete form)
    fim = np.sum((np.sqrt(P[mask]) - np.sqrt(Q[mask]))**2)
    return fim

# === Step 4: Compute FIM for each matched pair ===
fim_values = []
for air_col, tin_col in zip(air_cols, tin_cols):
    fim = compute_fim(df[air_col], df[tin_col])
    fim_values.append(fim)
    print(f"FIM({air_col} vs {tin_col}): {fim:.4f}")

# === Step 5: One-sample t-test vs 0 (null: no curvature difference) ===
t_stat, p_val = ttest_rel(fim_values, np.zeros_like(fim_values))
mean_fim = np.mean(fim_values)
std_fim = np.std(fim_values, ddof=1)
cohen_d = mean_fim / std_fim

# === Step 6: Report summary ===
print("\n=== Fisher Information Metric Summary ===")
print(f"Mean FIM: {mean_fim:.4f}")
print(f"Standard Deviation: {std_fim:.4f}")
print(f"t = {t_stat:.4f}, p = {p_val:.4f}")
print(f"Cohen's d: {cohen_d:.4f}")
