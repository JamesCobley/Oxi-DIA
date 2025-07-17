import pandas as pd
import numpy as np

# === Load your filtered dataframe with common sites ===
df = pd.read_csv('/content/redox_sites.tsv', sep='\t')

# === Define paired columns ===
air_cols = ['Sample_1_%Oxidized', 'Sample_2_%Oxidized', 'Sample_3_%Oxidized']
tin_cols = ['Sample_4_%Oxidized', 'Sample_5_%Oxidized', 'Sample_6_%Oxidized']

# === Map unique site identifier (optional) ===
df['Site_ID'] = df['Protein'] + "_" + df['Residue'].astype(str)

# === Keep rows with complete oxidation data ===
df_shared = df.dropna(subset=air_cols + tin_cols).copy()

# === Compute condition means ===
df_shared['Air_Mean'] = df_shared[air_cols].mean(axis=1)
df_shared['Tin_Mean'] = df_shared[tin_cols].mean(axis=1)

# === Classifier function ===
def classify_transformation(air, tin, threshold=5):
    delta = tin - air
    if np.isclose(delta, 0, atol=1):
        return 'Identity'
    elif abs(delta) <= threshold:
        return 'Scaling'
    elif (air < 30) and (tin > 70):
        return 'Bifurcation'
    elif abs(delta) > threshold:
        return 'Deformation'
    else:
        return 'Unclassified'

# === Apply transformation classification ===
df_shared['Transformation'] = df_shared.apply(
    lambda row: classify_transformation(row['Air_Mean'], row['Tin_Mean']),
    axis=1
)

# === Export classification table ===
df_shared[['Protein', 'Residue', 'Site', 'Air_Mean', 'Tin_Mean', 'Transformation']].to_csv(
    '/content/algebraic_redox_transformation_table.tsv', sep='\t', index=False
)

# === Summary count ===
print(df_shared['Transformation'].value_counts())
