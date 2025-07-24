import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(6, 5))
sns.set(style="whitegrid")

sns.scatterplot(
    data=merged,
    x='Shannon_Entropy',
    y='Fisher_Info',
    hue='Condition',
    palette='Set2',
    s=100,
    edgecolor='black'
)

plt.title("Redox Geometry Space\nMahalanobis D = 6.41 | MANOVA p = 0.015", fontsize=12)
plt.xlabel('Shannon Entropy (Disorder)', fontsize=12)
plt.ylabel('Fisher Information (Local Structure)', fontsize=12)
plt.legend(title='Condition')
plt.tight_layout()
plt.savefig('redox_geometry_space.png', dpi=300)
plt.show()
