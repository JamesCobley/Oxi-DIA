import matplotlib.pyplot as plt

# Sample labels
pairs = ['Pair 1', 'Pair 2', 'Pair 3']
kl_values = [0.0083, 0.0057, 0.0047]

# Plot
fig, ax = plt.subplots(figsize=(6, 4))

# Bar plot
ax.bar(pairs, kl_values, color='blue', edgecolor='black')
ax.axhline(y=np.mean(kl_values), color='red', linestyle='--', label=f'Mean KL = {np.mean(kl_values):.4f}')

# Formatting
ax.set_ylabel('KL Divergence (Air || Tin)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend()

# Annotate p-value and Cohen's d
p_val = 0.0269
cohen_d = 3.4483
ax.text(1.4, max(kl_values)*1.05, f'p = {p_val:.4f}\nd = {cohen_d:.2f}', 
        ha='center', va='bottom', fontsize=12)

plt.tight_layout()
plt.savefig('/content/kl_divergence.png', dpi=300)
plt.show()
