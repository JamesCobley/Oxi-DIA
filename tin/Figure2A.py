import matplotlib.pyplot as plt

# Data
air_values = np.array([1.68, 1.65, 1.64])
tin_values = np.array([1.24, 1.25, 1.44])
subjects = np.arange(1, len(air_values)+1)

# Plot
fig, ax = plt.subplots(figsize=(6, 4))

# Plot paired lines
for i in range(len(air_values)):
    ax.plot([0, 1], [air_values[i], tin_values[i]], marker='o', color='gray', linewidth=1)

# Overlay means ± SD
mean_air = np.mean(air_values)
std_air = np.std(air_values, ddof=1)
mean_tin = np.mean(tin_values)
std_tin = np.std(tin_values, ddof=1)

ax.errorbar(0, mean_air, yerr=std_air, fmt='o', color='blue', capsize=5, label='Air (mean ± SD)')
ax.errorbar(1, mean_tin, yerr=std_tin, fmt='o', color='red', capsize=5, label='Tin (mean ± SD)')

# Add p-value and Cohen’s d
p_val = 0.0429  # Replace with your actual value
cohens_d = 2.69  # Replace with your actual value

ax.text(0.5, max(air_values)*1.05, f"p = {p_val:.4f}\nd = {cohens_d:.2f}",
        ha='center', va='bottom', fontsize=12)


# Formatting
ax.set_xticks([0, 1])
ax.set_xticklabels(['Air', 'Tin'])
ax.set_ylabel('Cysteine oxidation (%)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend()

plt.tight_layout()
plt.savefig('/content/paired_differences.png', dpi=300)
plt.show()
