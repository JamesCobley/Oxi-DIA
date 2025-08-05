import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal, mannwhitneyu

# === Load data ===
df = pd.read_csv('/content/algebraic_redox_transformation_table_with_angles.tsv', sep='\t')

# === Filter for valid Transformation and Means ===
df_filtered = df.dropna(subset=['Transformation', 'Air_Mean', 'Tin_Mean']).copy()

# === Define site identifier ===
df_filtered['Site_ID'] = df_filtered['Protein'] + "_C" + df_filtered['Residue'].astype(str)

# === Define transformation T(ρ) = 2*Tin - Air ===
def apply_transformation(rho, air, tin):
    return tin + (tin - air)

# === Apply the transformation twice: T(T(ρ)) ===
df_filtered['T_rho'] = apply_transformation(df_filtered['Air_Mean'], df_filtered['Air_Mean'], df_filtered['Tin_Mean'])
df_filtered['T_T_rho'] = apply_transformation(df_filtered['T_rho'], df_filtered['Air_Mean'], df_filtered['Tin_Mean'])

# === Calculate deviation from identity ===
df_filtered['Delta_commute'] = (df_filtered['T_T_rho'] - df_filtered['Air_Mean']).abs()

# === Summary stats per transformation class ===
print("=== Mean Involution Error by Transformation Class ===")
summary = df_filtered.groupby('Transformation')['Delta_commute'].agg(['mean', 'std', 'count']).sort_values('mean', ascending=False)
print(summary.to_string(float_format="%.3f"))

# === Kruskal–Wallis test (non-parametric ANOVA) ===
groups = [group['Delta_commute'].values for name, group in df_filtered.groupby('Transformation')]
stat, p_kw = kruskal(*groups)
print(f"\nKruskal–Wallis H-test: H = {stat:.4f}, p = {p_kw:.4e}")

# === Pairwise Mann–Whitney U tests ===
transformations = df_filtered['Transformation'].unique()
print("\n=== Pairwise Mann–Whitney U-tests ===")
for i in range(len(transformations)):
    for j in range(i+1, len(transformations)):
        group1 = df_filtered[df_filtered['Transformation'] == transformations[i]]['Delta_commute']
        group2 = df_filtered[df_filtered['Transformation'] == transformations[j]]['Delta_commute']
        u_stat, p_val = mannwhitneyu(group1, group2, alternative='two-sided')
        print(f"{transformations[i]} vs {transformations[j]}: U = {u_stat:.2f}, p = {p_val:.4e}")

# === Plot ===
plt.figure(figsize=(8, 5))
sns.violinplot(data=df_filtered, x='Transformation', y='Delta_commute', inner=None, palette='viridis')
sns.stripplot(data=df_filtered, x='Transformation', y='Delta_commute', color='black', size=4, jitter=0.2)
plt.title('Involution Error by Transformation Class')
plt.ylabel('Δ |T(T(ρ)) - ρ|')
plt.xlabel('Transformation Class')
plt.tight_layout()
plt.savefig('involution_error_by_class.png', dpi=300)
plt.show()
