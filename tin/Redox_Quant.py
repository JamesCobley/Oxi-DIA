# STEP 1: Mount Google Drive
from google.colab import drive
drive.mount('/content/drive')
# STEP 2: Set paths to your DIA-NN output files
# Example paths — replace with your actual file locations inside your Google Drive
file1_path = '/content/drive/MyDrive/report.NEM_H_sites_90.tsv'
file2_path = '/content/drive/MyDrive/report.NEM_L_sites_90.tsv'
# STEP 3: Load the files using pandas (TSV format)
import pandas as pd

df_heavy = pd.read_csv(file1_path, sep='\t')
df_light = pd.read_csv(file2_path, sep='\t')
# STEP 4: Inspect column names and a preview of each dataframe
print("Heavy file columns:", df_heavy.columns.tolist())
print("Light file columns:", df_light.columns.tolist())

print("\nHeavy file preview:")
print(df_heavy.head())

print("\nLight file preview:")
print(df_light.head())


import pandas as pd
import os
import numpy as np

# === Step 5: Load light and heavy data ===
file_light = '/content/drive/MyDrive/report.NEM_L_sites_90.tsv'
file_heavy = '/content/drive/MyDrive/report.NEM_H_sites_90.tsv'

df_light = pd.read_csv(file_light, sep='\t')
df_heavy = pd.read_csv(file_heavy, sep='\t')

# === Step 6: Define site columns ===
site_cols = ['Protein', 'Residue', 'Site']
meta_cols = ['Protein.Names', 'Gene.Names', 'Sequence']

# === Step 7: Drop problematic raw file due to pipette contamination ===
excluded_raw = 'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Tin_6_new_S12.raw'  # pipette tip reuse noted

# Remove from both dataframes
df_light.drop(columns=[c for c in df_light.columns if excluded_raw in c], inplace=True, errors='ignore')
df_heavy.drop(columns=[c for c in df_heavy.columns if excluded_raw in c], inplace=True, errors='ignore')

# === Step 8: Rename run columns ===
light_cols_raw = [c for c in df_light.columns if c not in site_cols + meta_cols]
heavy_cols_raw = [c for c in df_heavy.columns if c not in site_cols + meta_cols]

light_renames = {c: f"{os.path.basename(c).replace('.raw', '')}_L" for c in light_cols_raw}
heavy_renames = {c: f"{os.path.basename(c).replace('.raw', '')}_H" for c in heavy_cols_raw}

df_light.rename(columns=light_renames, inplace=True)
df_heavy.rename(columns=heavy_renames, inplace=True)

# === Step 9: Merge L and H ===
df_light_sum = df_light.groupby(site_cols)[list(light_renames.values())].sum().reset_index()
df_heavy_sum = df_heavy.groupby(site_cols)[list(heavy_renames.values())].sum().reset_index()
df_redox = pd.merge(df_light_sum, df_heavy_sum, on=site_cols, how='outer').fillna(0)

# === Step 10: Compute %Oxidized for each run ===
run_names = sorted(set(c.replace('_L', '') for c in df_redox.columns if c.endswith('_L')))
for run in run_names:
    L = f"{run}_L"
    H = f"{run}_H"
    pct_col = f"{run}_%Oxidized"
    df_redox[pct_col] = 100 * df_redox[H] / (df_redox[L] + df_redox[H])
    df_redox[pct_col] = df_redox[pct_col].replace([np.inf, -np.inf], np.nan)

# === Step 11: Map runs to samples ===
def extract_sample_id(run):
    for part in run.split('_'):
        if part.isdigit():
            return int(part)
    return None

sample_map = {run: f"Sample_{extract_sample_id(run)}" for run in run_names}
oxi_cols = [f"{run}_%Oxidized" for run in run_names]
oxi_col_map = {f"{run}_%Oxidized": sample_map[run] for run in run_names if run in sample_map}
sample_names = sorted(set(oxi_col_map.values()))

# === Step 12: Identify cysteines quantified in all remaining technical runs ===
df_oxi = df_redox[site_cols + oxi_cols].copy()
mask_detected_all = df_oxi[oxi_cols].notna().all(axis=1)
df_common = df_oxi[mask_detected_all].copy()

# === Step 13: Average across technical replicates per sample ===
df_common_avg = df_common[site_cols].copy()
for sample in sample_names:
    sample_cols = [col for col, samp in oxi_col_map.items() if samp == sample]
    df_common_avg[f"{sample}_%Oxidized"] = df_common[sample_cols].mean(axis=1)

# === Step 14: Add metadata ===
df_meta = df_light[site_cols + meta_cols].drop_duplicates(subset=site_cols)
df_common_avg = pd.merge(df_common_avg, df_meta, on=site_cols, how='left')

# === Step 15: Compute summary ===
summary = []
for sample in sample_names:
    col_name = f"{sample}_%Oxidized"
    if col_name in df_common_avg.columns:
        oxi_vals = df_common_avg[col_name].dropna()
        summary.append({
            'Sample': sample,
            'Cysteines_Used': len(oxi_vals),
            'Mean_%Oxidized': round(oxi_vals.mean(), 2),
            'Mean_%Reduced': round(100 - oxi_vals.mean(), 2)
        })

summary_df = pd.DataFrame(summary)

# === Step 16: Save outputs ===
df_common_avg.to_csv('/content/redox_sites.tsv', sep='\t', index=False)
summary_df.to_csv('/content/sample_summary.tsv', sep='\t', index=False)

# === Show output ===
print("Excluded due to pipette error:", excluded_raw)
print(summary_df)
