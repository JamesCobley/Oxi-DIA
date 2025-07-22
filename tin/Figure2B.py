import matplotlib.pyplot as plt
import numpy as np

# === Your data ===
air_entropy = np.array([0.908128, 0.888091, 0.896930])
tin_entropy = np.array([0.835329, 0.837833, 0.850444])
samples = np.arange(1, 4)

# === Plot ===
fig, ax = plt.subplots(figsize=(6, 4))

for i in range(len(samples)):
    ax.plot([0, 1], [air_entropy[i], tin_entropy[i]], marker='o', color='gray', linewidth=1)

# Mean ± SD
mean_air = np.mean(air_entropy)
sd_air = np.std(air_entropy, ddof=1)
mean_tin = np.mean(tin_entropy)
sd_tin = np.std(tin_entropy, ddof=1)

ax.errorbar(0, mean_air, yerr=sd_air, fmt='o', color='blue', capsize=5, label='Air (mean ± SD)')
ax.errorbar(1, mean_tin, yerr=sd_tin, fmt='o', color='red', capsize=5, label='Tin (mean ± SD)')

# Axes and annotation
ax.set_xticks([0, 1])
ax.set_xticklabels(['Air', 'Tin'])
ax.set_ylabel('Shannon Entropy (bits)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend()

# p-value and effect size
p_val = 0.0205
cohen_d = 3.9720

ax.text(0.5, max(air_entropy)*1.03, f'p = {p_val:.4f}\nd = {cohen_d:.2f}', 
        ha='center', va='bottom', fontsize=12)

plt.tight_layout()
plt.savefig('/content/entropy_per_sample.png', dpi=300)
plt.show()
