import pandas as pd
import numpy as np
from numpy.linalg import norm
from math import acos, degrees
import matplotlib.pyplot as plt
import seaborn as sns

# === Load fully complete redox matrix ===
df = pd.read_excel('/content/fully_complete_redox_matrix.xlsx')

# === Extract numeric columns (samples) ===
df_numeric = df.select_dtypes(include=['number'])

# === Compute site-wise redox vector and stats ===
df['Redox_Vector'] = df_numeric.values.tolist()
df['Mean_Redox'] = df_numeric.mean(axis=1)
df['Std_Redox'] = df_numeric.std(axis=1)

# === Reference vector: mean per sample ===
ref_vector = df_numeric.mean(axis=0).values

# === Angle from mean vector ===
def angle_from_mean(v, ref):
    v = np.array(v)
    ref = np.array(ref)
    if np.allclose(v, 0) or np.allclose(ref, 0):
        return 0.0
    cos_sim = np.dot(v, ref) / (norm(v) * norm(ref))
    cos_sim = np.clip(cos_sim, -1.0, 1.0)
    return degrees(acos(cos_sim))

df['Angle_Degrees'] = df['Redox_Vector'].apply(lambda vec: angle_from_mean(vec, ref_vector))

# === Transformation classification with bifurcation ===
def classify_population_geometry(vec, std, angle, low_thresh=30, high_thresh=70):
    vec = np.array(vec)
    low = np.sum(vec < low_thresh)
    high = np.sum(vec > high_thresh)
    total = len(vec)

    if std <= 1 and angle < 1:
        return 'Identity'
    elif low >= total * 0.3 and high >= total * 0.3:
        return 'Bifurcation'
    elif std <= 5:
        return 'Scaling'
    else:
        return 'Deformation'

df['Transformation'] = df.apply(
    lambda row: classify_population_geometry(row['Redox_Vector'], row['Std_Redox'], row['Angle_Degrees']),
    axis=1
)

# === Save transformation-classified dataframe ===
df.to_csv('/content/transformation_table_with_vectors.csv', index=False)
print("✅ Saved transformation-classified redox table to: /content/transformation_table_with_vectors.csv")

# === Print summary ===
print("\nTransformation Class Counts:")
print(df['Transformation'].value_counts())

# === Summary statistics ===
print("\nAverage Angular Deformation per Class:")
print(df.groupby('Transformation')['Angle_Degrees'].agg(['mean', 'std', 'count']).round(2))

print("\nAverage Mean Redox State per Class:")
print(df.groupby('Transformation')['Mean_Redox'].agg(['mean', 'std', 'count']).round(2))

# === Combined Violin Plot ===
sns.set(style="whitegrid")
plt.figure(figsize=(14, 6))

# Subplot 1: Angular deviation
plt.subplot(1, 2, 1)
sns.violinplot(data=df, x='Transformation', y='Angle_Degrees', palette='Set2')
plt.title('Angular Deformation by Transformation Class')
plt.ylabel('Angle (Degrees)')
plt.xlabel('')

# Subplot 2: Mean redox state
plt.subplot(1, 2, 2)
sns.violinplot(data=df, x='Transformation', y='Mean_Redox', palette='Set3')
plt.title('Mean Redox State by Transformation Class')
plt.ylabel('% Oxidation')
plt.xlabel('')

plt.tight_layout()
plt.savefig("/content/transformation_violinplots_angle_redox.png", dpi=300)
plt.show()
