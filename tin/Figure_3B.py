import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# === Load the redox data ===
df = pd.read_csv('/content/redox_sites.tsv', sep='\t')

# === Define site ID (if needed) ===
df['Site_ID'] = df['Protein'] + "_" + df['Residue'].astype(str)

# === Define sample columns explicitly ===
air_cols = ['Sample_1_%Oxidized', 'Sample_2_%Oxidized', 'Sample_3_%Oxidized']
tin_cols = ['Sample_4_%Oxidized', 'Sample_5_%Oxidized', 'Sample_6_%Oxidized']

# === Filter rows with complete data ===
df_filtered = df.dropna(subset=air_cols + tin_cols).copy()

# === Step 1: Compute condition-wise means ===
df_filtered['Air_mean'] = df_filtered[air_cols].mean(axis=1)
df_filtered['Tin_mean'] = df_filtered[tin_cols].mean(axis=1)

# === Step 2: Define transformation T = Tin - Air ===
df_filtered['T'] = df_filtered['Tin_mean'] - df_filtered['Air_mean']

# === Step 3: Apply T again to get T(T(rho)) ===
df_filtered['T_T'] = df_filtered['Tin_mean'] + df_filtered['T']

# === Step 4: Check for involution (T(T(rho)) ≈ rho) ===
df_filtered['Delta_commute'] = (df_filtered['T_T'] - df_filtered['Air_mean']).abs()

# === Report mean deviation ===
mean_error = df_filtered['Delta_commute'].mean()
print(f"Mean commutator-like error: {mean_error:.4f}")

# === Optional: plot histogram ===
plt.figure(figsize=(8, 5))
plt.hist(df_filtered['Delta_commute'], bins=50, color='skyblue')
plt.title('Deviation of T(T(ρ)) from Air — Non-commutativity Measure')
plt.xlabel('Deviation Magnitude')
plt.ylabel('Count')
plt.tight_layout()
plt.savefig('/content/commutator_error.png', dpi=300)
plt.show()
