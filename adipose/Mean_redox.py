import pandas as pd

# Load the fully complete matrix
file_path = "/content/fully_complete_redox_matrix.xlsx"
df = pd.read_excel(file_path)

# Optional: check that non-metadata columns are numeric
df_numeric = df.select_dtypes(include=['number'])

# Compute redox state per biological sample (average oxidation %)
redox_state_per_sample = df_numeric.mean(axis=0)

# Report basic statistics
summary_df = pd.DataFrame({
    "Sample": redox_state_per_sample.index,
    "Average_Redox_%": redox_state_per_sample.values
})
summary_df["SD"] = df_numeric.std(axis=0).values
summary_df["CV_%"] = 100 * (summary_df["SD"] / summary_df["Average_Redox_%"])

# Save the summary
summary_df.to_excel("/content/sample_redox_states.xlsx", index=False)

# Print
print("✅ Average redox state per biological sample:")
print(summary_df.round(2))

import matplotlib.pyplot as plt
import seaborn as sns

# === Compute overall redox state statistics ===
overall_mean = summary_df["Average_Redox_%"].mean()
overall_sd = summary_df["Average_Redox_%"].std()
overall_cv = 100 * (overall_sd / overall_mean)

print("\n📊 Overall Redox State Statistics:")
print(f"Mean Redox State: {overall_mean:.2f}%")
print(f"Standard Deviation: {overall_sd:.2f}")
print(f"Coefficient of Variation (CV): {overall_cv:.2f}%")

# === Scatter plot of average redox state per sample ===
plt.figure(figsize=(10, 5))
sns.scatterplot(data=summary_df, x="Sample", y="Average_Redox_%", s=100, color='blue')

# Optional: Add value labels on top of points
for i, row in summary_df.iterrows():
    plt.text(i, row["Average_Redox_%"] + 0.1, f"{row['Average_Redox_%']:.2f}", 
             ha='center', fontsize=8)

plt.xlabel("Sample")
plt.ylabel("Average Cysteine Oxidation (%)")
plt.xticks(rotation=45, ha='right')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig("/content/sample_redox_scatterplot.png", dpi=300)
plt.show()
