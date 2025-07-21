import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

# === Load matrix ===
df = pd.read_excel("/content/fully_complete_redox_matrix.xlsx")
df_numeric = df.select_dtypes(include=["number"])  # Only sample columns

# === Parameters ===
samples = df_numeric.columns
n_samples = len(samples)
n_cols = 5  # adjust grid layout as desired
n_rows = math.ceil(n_samples / n_cols)

# === Reshape function: force 2D for image display ===
def reshape_for_display(values, width=None):
    if width is None:
        width = int(np.ceil(np.sqrt(len(values))))
    pad_len = width**2 - len(values)
    padded = np.append(values, [np.nan]*pad_len)
    return padded.reshape((width, width))

# === Plot image grid ===
fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 3))

for idx, sample in enumerate(samples):
    row = idx // n_cols
    col = idx % n_cols
    ax = axes[row, col] if n_rows > 1 else axes[col]
    values = df_numeric[sample].values
    image = reshape_for_display(values)
    im = ax.imshow(image, cmap='viridis', vmin=0, vmax=100)
    ax.set_title(sample, fontsize=8)
    ax.axis('off')

# Remove empty subplots if any
for i in range(n_samples, n_rows * n_cols):
    fig.delaxes(axes.flatten()[i])

# Add colorbar
cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
fig.colorbar(im, cax=cbar_ax, label='% Oxidation')

fig.suptitle("Redox Proteome Fingerprints (Cysteine Oxidation)", fontsize=12)
plt.tight_layout(rect=[0, 0, 0.9, 0.95])
plt.savefig("/content/redox_fingerprint_grid.png", dpi=300)
plt.show()
