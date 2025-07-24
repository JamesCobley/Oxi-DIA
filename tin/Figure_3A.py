import matplotlib.pyplot as plt
import seaborn as sns

# Ensure you have df_shared with columns: Air_Mean, Tin_Mean, Transformation

# Define a color palette
palette = {'Identity': '#66c2a5', 'Scaling': '#fc8d62', 'Deformation': '#8da0cb'}

# Set up 3 plots
fig, axs = plt.subplots(1, 3, figsize=(18, 5), sharex=True, sharey=True)

# Loop over transformation types
for ax, label in zip(axs, ['Identity', 'Scaling', 'Deformation']):
    subset = df_shared[df_shared['Transformation'] == label]
    sns.scatterplot(
        data=subset,
        x='Air_Mean',
        y='Tin_Mean',
        ax=ax,
        color=palette[label],
        alpha=0.5,
        s=10
    )
    ax.plot([0, 100], [0, 100], '--', color='gray', linewidth=1)  # y=x diagonal
    ax.set_title(f"{label} (n={len(subset)})")
    ax.set_xlabel('%Reduced (Air)')
    ax.set_ylabel('%Reduced (Tin)')

plt.suptitle("Redox Transformation Classes: Air vs. Tin", fontsize=16)
plt.tight_layout()
plt.savefig('/content/transformation_classes.png', dpi=300)
plt.show()
