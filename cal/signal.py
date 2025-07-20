import pandas as pd

# Load your DIA-NN file
file_path = '/content/drive/MyDrive/reportnewcal.parquet'  # <- Change this!
df = pd.read_parquet(file_path)

# Ensure clean values
df = df.copy()
df['Ms1.Area'] = pd.to_numeric(df['Ms1.Area'], errors='coerce')
df['Precursor.Quantity'] = pd.to_numeric(df['Precursor.Quantity'], errors='coerce')
df = df.fillna(0)

# Setup
summary = []
all_runs = sorted(df['Run'].unique())

for run in all_runs:
    df_run = df[df['Run'] == run]

    # Unique peptides
    stripped_seqs = df_run['Stripped.Sequence'].dropna().unique()
    num_total_peptides = len(stripped_seqs)
    
    # Unique C-containing peptides
    cys_seqs = [seq for seq in stripped_seqs if 'C' in seq]
    num_c_peptides = len(cys_seqs)

    # Total intensities
    total_ms1 = df_run['Ms1.Area'].sum()
    total_ms2 = df_run['Precursor.Quantity'].sum()

    # Cys-only intensities
    df_cys = df_run[df_run['Stripped.Sequence'].isin(cys_seqs)]
    c_ms1 = df_cys['Ms1.Area'].sum()
    c_ms2 = df_cys['Precursor.Quantity'].sum()

    summary.append({
        'Run': run,
        'Total_Peptides': num_total_peptides,
        'C_Peptides': num_c_peptides,
        'Total_MS1': total_ms1,
        'C_MS1': c_ms1,
        'Total_MS2': total_ms2,
        'C_MS2': c_ms2
    })

# Final DataFrame
summary_df = pd.DataFrame(summary)
summary_df = summary_df.sort_values(by='Run')

# Display all rows
pd.set_option('display.max_rows', None)
print(summary_df)

# Optional: Save to file
summary_df.to_csv('/content/global_peptide_intensity_summary.csv', index=False)

import pandas as pd
import matplotlib.pyplot as plt

# Paste your table into a CSV or load from DataFrame directly
df = pd.read_csv('/content/global_peptide_intensity_summary.csv')  # Replace if using direct DataFrame

# Compute proportions
df['C_Peptide_%'] = (df['C_Peptides'] / df['Total_Peptides']) * 100
df['C_MS1_%'] = (df['C_MS1'] / df['Total_MS1']) * 100
df['C_MS2_%'] = (df['C_MS2'] / df['Total_MS2']) * 100

# Sort by expected oxidation (infer from Run names)
df['Oxidation'] = df['Run'].str.extract(r'James_(\d+)_').astype(int)
df = df.sort_values(by='Oxidation')

# Plot
plt.figure(figsize=(12, 6))
plt.plot(df['Oxidation'], df['C_MS1_%'], marker='o', label='MS1 Cys %')
plt.plot(df['Oxidation'], df['C_MS2_%'], marker='s', label='MS2 Cys %')
plt.plot(df['Oxidation'], df['C_Peptide_%'], marker='^', label='Cys Peptide %')
plt.axhline(20, color='gray', linestyle='--', label='Reference Line (e.g., expected)')
plt.xlabel('% Target Oxidation (inferred from Run name)')
plt.ylabel('Percentage')
plt.title('Cysteine Contribution vs. Target Oxidation')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('cysteine_contribution.png')
plt.show()
