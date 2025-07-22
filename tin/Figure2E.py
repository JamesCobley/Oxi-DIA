import matplotlib.pyplot as plt
import numpy as np

# Sample labels and MI values
pairs = ['Pair 1', 'Pair 2', 'Pair 3']
mi_values = [0.2341, 0.2368, 0.2495]
mean_mi = np.mean(mi_values)
p_val = 0.0004
cohen_d = 29.0953

# Plot
fig, ax = plt.subplots(figsize=(6, 4))

# Bar plot
ax.bar(pairs, mi_values, color='blue', edgecolor='black')
ax.axhline(y=mean_mi, color='red', linestyle='--', label=f'Mean MI = {mean_mi:.4f}')

# Formatting
ax.set_ylabel('Mutual Information (Air vs Tin)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend()

# Annotate p-value and Cohen’s d
ax.text(1.4, max(mi_values)*1.05, f'p = {p_val:.4f}\nd = {cohen_d:.2f}', 
        ha='center', va='bottom', fontsize=12)

plt.tight_layout()
plt.savefig('/content/mutual_information.png', dpi=300)
plt.show()
