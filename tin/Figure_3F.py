import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# === Step 1: Define target GO terms and names ===
target_go_terms = {
    'GO:0045202': 'Synapse',
    'GO:0007416': 'Synapse assembly',
    'GO:0060074': 'Synapse maturation',
    'GO:0099005': 'Postsynapse',
    'GO:0099004': 'Presynapse',
    'GO:0098978': 'Glutamatergic synapse',
    'GO:0098981': 'Cholinergic synapse'
}

# === Step 2: Filter transformations ===
filtered_df = df[df['Transformation'].isin(['Scaling', 'Deformation'])].copy()
filtered_df['Delta'] = filtered_df['Air_Mean'] - filtered_df['Tin_Mean']

# === Step 3: Explode GO terms ===
filtered_df['GO_Terms'] = filtered_df['GO_Terms'].fillna('')
filtered_df = filtered_df.assign(GO_Term=filtered_df['GO_Terms'].str.split(';')).explode('GO_Term')

# === Step 4: Filter to target GO terms ===
manual_report = filtered_df[filtered_df['GO_Term'].isin(target_go_terms.keys())].copy()

# === Step 5: Aggregate mean delta ===
heatmap_data = (
    manual_report
    .groupby('GO_Term')
    .agg(mean_delta=('Delta', 'mean'))
    .rename(index=target_go_terms)  # replace GO IDs with human-readable names
    .sort_values('mean_delta', ascending=True)  # optional: sort by effect size
)

# === Step 6: Plot heatmap ===
plt.figure(figsize=(6, 4))
sns.heatmap(
    heatmap_data,
    cmap='viridis',
    annot=True,
    fmt=".2f",
    linewidths=0.5,
    cbar_kws={'label': 'Mean ΔOxidation (Air - Tin)'}
)
plt.title("Oxidation Shift by Synapse-Related GO Term")
plt.ylabel("GO Term")
plt.xlabel("")
plt.tight_layout()
plt.savefig("go_term_oxidation_heatmap.png", dpi=300)
plt.show()
