import matplotlib.pyplot as plt
import numpy as np

# Sample labels and FIM values
pairs = ['Pair 1', 'Pair 2', 'Pair 3']
fim_values = [0.0038, 0.0026, 0.0025]
mean_fim = np.mean(fim_values)
p_val = 0.0187
cohen_d = 4.1589

# Plot
fig, ax = plt.subplots(figsize=(6, 4))

# Bar plot
ax.bar(pairs, fim_values, color='blue', edgecolor='black')
ax.axhline(y=mean_fim, color='red', linestyle='--', label=f'Mean FIM = {mean_fim:.4f}')

# Formatting
ax.set_ylabel('Fisher Information Metric (Air vs Tin)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend()

# Annotate p-value and Cohen’s d
ax.text(1.4, max(fim_values)*1.05, f'p = {p_val:.4f}\nd = {cohen_d:.2f}', 
        ha='center', va='bottom', fontsize=12)

plt.tight_layout()
plt.savefig('/content/fisher_information_metric.png', dpi=300)
plt.show()
