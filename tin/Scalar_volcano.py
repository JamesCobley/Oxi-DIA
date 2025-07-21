import pandas as pd
import numpy as np
from scipy.stats import ttest_rel
import matplotlib.pyplot as plt

# === Load your filtered dataframe with common sites ===
df = pd.read_csv('/content/redox_sites.tsv', sep='\t')

# === Confirm column names ===
print("Columns:", df.columns.tolist())

# === Map site identifier ===
df['Site_ID'] = df['Protein'] + "_" + df['Residue'].astype(str)

# === Define paired columns (based on Sample IDs) ===
air_cols = ['Sample_1_%Oxidized', 'Sample_2_%Oxidized', 'Sample_3_%Oxidized']
tin_cols = ['Sample_4_%Oxidized', 'Sample_5_%Oxidized', 'Sample_6_%Oxidized']

results = []

for idx, row in df.iterrows():
    air_vals = row[air_cols].values.astype(float)
    tin_vals = row[tin_cols].values.astype(float)

    if np.isnan(air_vals).any() or np.isnan(tin_vals).any():
        continue

    t_stat, p_val = ttest_rel(air_vals, tin_vals)
    fc = np.mean(air_vals) - np.mean(tin_vals)
    log2fc = np.log2((np.mean(air_vals) + 1e-6) / (np.mean(tin_vals) + 1e-6))

    results.append({
        'Site_ID': row['Site_ID'],
        'Protein.Names': row.get('Protein.Names', ''),
        'Gene.Names': row.get('Gene.Names', ''),
        'Log2FC': log2fc,
        'p-value': p_val,
        'Mean_Air': np.mean(air_vals),
        'Mean_Tin': np.mean(tin_vals),
        'Delta_Oxidation': fc
    })

# === Create volcano dataframe ===
volcano_df = pd.DataFrame(results)
volcano_df['-log10(p-value)'] = -np.log10(volcano_df['p-value'])

from statsmodels.stats.multitest import multipletests

# Compute adjusted p-values (Benjamini-Hochberg)
volcano_df['adj_pval'] = multipletests(volcano_df['p-value'], method='fdr_bh')[1]

# Define significance threshold with FDR
volcano_df['Significant'] = (volcano_df['adj_pval'] < 0.05) & (abs(volcano_df['Log2FC']) > 0.5)

# === Save the data ===
volcano_df.to_csv('/content/sitewise_paired_ttest_results.tsv', sep='\t', index=False)

plt.figure(figsize=(10,6))
plt.scatter(volcano_df['Log2FC'], volcano_df['-log10(p-value)'], color='gray', alpha=0.5)

# Highlight FDR-significant hits
sig_mask = volcano_df['Significant']
plt.scatter(volcano_df[sig_mask]['Log2FC'], volcano_df[sig_mask]['-log10(p-value)'],
            color='red', label='FDR < 0.05')

plt.axhline(-np.log10(0.05), color='blue', linestyle='--')
plt.axvline(-0.5, color='black', linestyle='--')
plt.axvline(0.5, color='black', linestyle='--')

plt.xlabel('Log2 Fold Change (Air - Tin)')
plt.ylabel('-Log10 p-value')
plt.title('Volcano Plot of Cysteine Redox Differences (FDR corrected)')
plt.legend()
plt.tight_layout()
plt.savefig('/content/volcano_plot_fdr.png', dpi=300)
plt.show()
