import pandas as pd
import numpy as np
from numpy.linalg import norm
from math import acos, degrees

# === Load your filtered dataframe with common sites ===
df = pd.read_csv('/content/redox_sites (2).tsv', sep='\t')

# === Define paired columns ===
air_cols = ['Sample_1_%Oxidized', 'Sample_2_%Oxidized', 'Sample_3_%Oxidized']
tin_cols = ['Sample_4_%Oxidized', 'Sample_5_%Oxidized', 'Sample_6_%Oxidized']

# === Map unique site identifier ===
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

# === Compute redox vectors ===
df_shared['Air_Vector'] = df_shared[air_cols].values.tolist()
df_shared['Tin_Vector'] = df_shared[tin_cols].values.tolist()

# === Cosine angle computation with flat vector handling ===
def cosine_angle(v1, v2):
    v1 = np.array(v1)
    v2 = np.array(v2)
    if np.allclose(v1, 0) or np.allclose(v2, 0):
        return 0.0  # Define flat vectors as 0°
    cos_sim = np.dot(v1, v2) / (norm(v1) * norm(v2))
    cos_sim = np.clip(cos_sim, -1.0, 1.0)
    return acos(cos_sim)

def redox_direction_angle(air, tin):
    v = np.array([air, tin])
    ref = np.array([1, 1])  # Identity (no redox change)
    if np.allclose(v, 0):
        return 0.0  # Define flat or inactive redox sites as 0° deviation
    cos_sim = np.dot(v, ref) / (norm(v) * norm(ref))
    cos_sim = np.clip(cos_sim, -1.0, 1.0)
    return acos(cos_sim)


df_shared['Angle_Radians'] = df_shared.apply(
    lambda row: redox_direction_angle(row['Air_Mean'], row['Tin_Mean']),
    axis=1
)
df_shared['Angle_Degrees'] = df_shared['Angle_Radians'].apply(degrees)

# === Flag flat vectors ===
def is_flat(v1, v2):
    return np.allclose(v1, 0) or np.allclose(v2, 0)

df_shared['Flat_Vector'] = df_shared.apply(
    lambda row: is_flat(row['Air_Vector'], row['Tin_Vector']),
    axis=1
)

# === Summary: transformation counts ===
print("\nTransformation Class Counts:")
print(df_shared['Transformation'].value_counts())

# === Summary: average angle per transformation class ===
angle_summary = df_shared.groupby('Transformation')['Angle_Degrees'].agg(['mean', 'std', 'count']).round(3)
print("\nAverage Angular Deformation per Transformation Class (Degrees):")
print(angle_summary)

# === Compute Total Mean Redox per Site (Air + Tin) ===
df_shared['Mean_Redox_State'] = df_shared[air_cols + tin_cols].mean(axis=1)

# === Summary: Mean redox state per transformation class ===
redox_summary = df_shared.groupby('Transformation')['Mean_Redox_State'].agg(['mean', 'std', 'count']).round(2)
print("\nMean Redox State per Transformation Class (Air + Tin Average %Oxidation):")
print(redox_summary)

# === Export updated transformation sheet with angle data ===
df_shared[['Protein', 'Residue', 'Site', 'Air_Mean', 'Tin_Mean', 'Mean_Redox_State','Transformation',
           'Angle_Degrees', 'Angle_Radians', 'Flat_Vector']].to_csv(
    '/content/algebraic_redox_transformation_table_with_angles.tsv',
    sep='\t', index=False
)

import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import matplotlib.pyplot as plt
import seaborn as sns

# === Load the updated dataframe (or reuse df_shared if already in memory) ===
df = pd.read_csv('/content/algebraic_redox_transformation_table_with_angles.tsv', sep='\t')

# === Filter to valid transformation classes (if needed) ===
valid_classes = ['Identity', 'Scaling', 'Deformation']
df = df[df['Transformation'].isin(valid_classes)]

# === ANOVA: Angular Deformation ===
anova_angle = stats.f_oneway(
    *[df[df['Transformation'] == t]['Angle_Degrees'] for t in valid_classes]
)
print("\n=== ANOVA: Angle_Degrees by Transformation ===")
print(f"F = {anova_angle.statistic:.4f}, p = {anova_angle.pvalue:.4e}")

# === Tukey HSD post-hoc: Angle ===
tukey_angle = pairwise_tukeyhsd(df['Angle_Degrees'], df['Transformation'])
print("\n=== Tukey HSD: Angle_Degrees ===")
print(tukey_angle)

# === ANOVA: Mean Redox State ===
anova_redox = stats.f_oneway(
    *[df[df['Transformation'] == t]['Mean_Redox_State'] for t in valid_classes]
)
print("\n=== ANOVA: Mean_Redox_State by Transformation ===")
print(f"F = {anova_redox.statistic:.4f}, p = {anova_redox.pvalue:.4e}")

# === Tukey HSD post-hoc: Redox ===
tukey_redox = pairwise_tukeyhsd(df['Mean_Redox_State'], df['Transformation'])
print("\n=== Tukey HSD: Mean_Redox_State ===")
print(tukey_redox)

# === Optional: Boxplots ===
sns.set(style="whitegrid")

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
sns.boxplot(x='Transformation', y='Angle_Degrees', data=df, palette='Set2')
plt.title('Angular Deformation by Transformation')
plt.ylabel('Angle (Degrees)')

plt.subplot(1, 2, 2)
sns.boxplot(x='Transformation', y='Mean_Redox_State', data=df, palette='Set3')
plt.title('Mean Redox State by Transformation')
plt.ylabel('% Oxidation')

plt.tight_layout()
plt.show()
