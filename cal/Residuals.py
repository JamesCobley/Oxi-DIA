import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

# === 1. Load your data ===
# Replace with your actual CSV path
df = pd.read_csv("redox_comparison.csv")  # Ensure columns: Expected_Reduced, Observed_Reduced_MS1, Observed_Reduced_MS2, Observed_Reduced_Combined

# === 2. Safe MAPE implementation ===
def safe_mape(y_true, y_pred):
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    mask = y_true != 0
    if np.any(mask):
        return np.mean(np.abs((y_pred[mask] - y_true[mask]) / y_true[mask]))
    else:
        return np.nan  # or 0 if you prefer

# === 3. Define metric computation ===
def compute_metrics(observed, expected, label):
    r2 = r2_score(expected, observed)
    mse = mean_squared_error(expected, observed)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(expected, observed)
    mape = safe_mape(expected, observed)  # <--- using safe version here
    residuals = observed - expected
    bias = np.mean(residuals)
    sd_residuals = np.std(residuals)

    print(f"\nMetrics for {label}:")
    print(f"  R²     = {r2:.4f}")
    print(f"  RMSE   = {rmse:.4f}")
    print(f"  MAE    = {mae:.4f}")
    print(f"  MAPE   = {mape * 100:.2f}%")
    print(f"  Bias   = {bias:.4f}")
    print(f"  SD(residuals) = {sd_residuals:.4f}")

    return {
        "r2": r2,
        "rmse": rmse,
        "mae": mae,
        "mape": mape,
        "bias": bias,
        "sd": sd_residuals,
        "residuals": residuals
    }

# === 4. Compute metrics for each mode ===
results = {}
for col in ['Observed_Reduced_MS1', 'Observed_Reduced_MS2', 'Observed_Reduced_Combined']:
    label = col.split('_')[-1] if 'Combined' in col else col.split('_')[-1]
    results[label] = compute_metrics(df[col], df['Expected_Reduced'], label)

plt.figure(figsize=(10, 6))
for label, color, marker in zip(['MS1', 'MS2', 'Combined'], ['blue', 'green', 'red'], ['o', 's', '^']):
    plt.scatter(df['Expected_Reduced'], results[label]['residuals'], label=f"{label} Residuals", color=color, marker=marker, alpha=0.7)

plt.axhline(0, color='black', linestyle='--')
plt.axhspan(-0.05, 0.05, color='gray', alpha=0.1, label='±5% Band')
plt.xlabel("Expected Reduced (%)")
plt.ylabel("Residual (Observed - Expected)")
plt.title("Residuals of Oxi-DIA Quantification Across Reduction Range")
plt.ylim(-0.1, 0.1)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('residuals_clean.png', dpi=300)
plt.show()


# === 6. Export residuals ===
df['Residual_MS1'] = results['MS1']['residuals']
df['Residual_MS2'] = results['MS2']['residuals']
df['Residual_Combined'] = results['Combined']['residuals']
df.to_csv("oxi_dia_with_residuals.csv", index=False)
