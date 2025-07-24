plt.figure(figsize=(10, 8))
sns.set(style="whitegrid")

sns.heatmap(
    fisher_df,
    cmap="magma",
    square=True,
    linewidths=0.5,
    cbar_kws={"shrink": 0.5},
    annot=False,
    xticklabels=True,
    yticklabels=True
)

plt.title("Fisher–Rao Distance Heatmap (%Oxidized)", fontsize=14)
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.yticks(rotation=0, fontsize=9)
plt.tight_layout()
plt.savefig("fisher_rao_heatmap.png", dpi=300)
plt.show()
