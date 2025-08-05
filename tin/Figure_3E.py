import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_rel

# === Step 1: Load the raw redox data ===
df = pd.read_csv('/content/redox_sites (3).tsv', sep='\t')

# === Step 2: Define sample columns ===
air_cols = ['Sample_1_%Oxidized', 'Sample_2_%Oxidized', 'Sample_3_%Oxidized']
tin_cols = ['Sample_4_%Oxidized', 'Sample_5_%Oxidized', 'Sample_6_%Oxidized']

# === Step 3: Compute means ===
df['Air_Mean'] = df[air_cols].mean(axis=1)
df['Tin_Mean'] = df[tin_cols].mean(axis=1)
df['Delta_Oxidation'] = df['Air_Mean'] - df['Tin_Mean']

# === Step 4: Add site labels ===
df['Label'] = df['Protein'] + '_C' + df['Site'].astype(str)

# === Step 5: Define sites of interest ===
target_sites = {
    'Q03157_C181': 'APHP1',
    'Q61292_C969': 'LAMB2',
    'Q99JR5_C53': 'TINAL',
    'P06802_C607': 'ENPP1',
    'Q8K406_C43': 'LGI3'
}

df['Site_Label'] = df['Label'].map(target_sites)
df_filtered = df[df['Site_Label'].notna()].copy()

# === Step 6: Melt to long format for plotting ===
df_melted = df_filtered.melt(
    id_vars=['Site_Label'],
    value_vars=air_cols + tin_cols,
    var_name='Sample',
    value_name='% Oxidation'
)

# === Step 7: Annotate condition from sample name ===
df_melted['Condition'] = df_melted['Sample'].apply(lambda x: 'Air' if '1' in x or '2' in x or '3' in x else 'Tin')

# === Step 8: Plot per site ===
g = sns.catplot(
    data=df_melted,
    x='Condition',
    y='% Oxidation',
    col='Site_Label',
    kind='box',
    palette='Set2',
    height=4,
    aspect=0.8,
    showfliers=False
)

# Overlay swarm for each box
sns.set_style('whitegrid')
for ax, site in zip(g.axes.flat, df_filtered['Site_Label'].unique()):
    site_data = df_melted[df_melted['Site_Label'] == site]
    sns.swarmplot(
        data=site_data,
        x='Condition',
        y='% Oxidation',
        ax=ax,
        color='black',
        size=6
    )
    ax.set_title(site, fontsize=12)


plt.tight_layout()
plt.savefig("oxidation_per_site_boxplot.png", dpi=300)
plt.show()
