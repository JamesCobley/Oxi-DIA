import pandas as pd
import matplotlib.pyplot as plt

# Load your file
df = pd.read_csv('algebraic_redox_transformation_table_with_angles.tsv', sep='\t')

# Rebuild labels
df['Label'] = df['Protein'] + '_C' + df['Site'].astype(str)

# Define your targets
targets = {
    'A2AJA9_C254': 'AJM1 (Identity)',
    'A2ALS5_C146': 'RPGP1 (Scaling)',
    'A2ASQ1_C1012': 'AGRIN (Deformation)'
}

# Filter and map labels
df_targets = df[df['Label'].isin(targets.keys())].copy()
df_targets['Protein_Label'] = df_targets['Label'].map(targets)

# Plot
fig, ax = plt.subplots(figsize=(6, 5))
for _, row in df_targets.iterrows():
    x = row['Mean_Redox_State']
    y = row['Angle_Degrees']
    ax.scatter(x, y, s=80, color='black')
    ax.text(x + 1.5, y, f"{row['Protein_Label']}\n{y:.1f}°", ha='left', va='center', fontsize=9)

ax.set_xlabel('Mean Redox State (% Oxidation)', fontsize=12)
ax.set_ylabel('Angular Deformation (Degrees)', fontsize=12)
ax.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig("fig3e_redox_vs_angle.png", dpi=300)
plt.show()
