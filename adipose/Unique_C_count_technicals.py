import pandas as pd
import numpy as np

# Load the DIA-NN output
file_path = '/content/new cys test.pr_matrix (2).tsv'
df = pd.read_csv(file_path, sep='\t')

# Detect sample columns dynamically: any column after 'Precursor.Id'
metadata_cols = list(df.columns[:11])
sample_cols = list(df.columns[11:])

# Clean up sample column names
clean_sample_cols = [
    col.split('\\')[-1].split('/')[-1].replace('.raw', '').replace('.mzML', '')
    for col in sample_cols
]

# Extract relevant columns
df['HasCys'] = df['Stripped.Sequence'].str.contains('C')
df['Peptide_Key'] = df['Stripped.Sequence'] + '_z' + df['Precursor.Charge'].astype(str)

# --- STATS ---

# 1. Count of unique cysteine-containing peptide sequences per sample (sequence-level)
unique_cys_seq_counts = {}
for sample in sample_cols:
    subset = df[(df['HasCys']) & (df[sample].notna())]
    unique_seqs = subset['Stripped.Sequence'].unique()
    unique_cys_seq_counts[sample] = len(unique_seqs)
cys_seq_df = pd.DataFrame.from_dict(unique_cys_seq_counts, orient='index', columns=['Unique Cys-Seqs'])

# 2. Count of unique cysteine precursors (seq + charge) per sample
unique_cys_precursor_counts = {}
for sample in sample_cols:
    subset = df[(df['HasCys']) & (df[sample].notna())]
    unique_precursors = subset['Peptide_Key'].unique()
    unique_cys_precursor_counts[sample] = len(unique_precursors)
cys_prec_df = pd.DataFrame.from_dict(unique_cys_precursor_counts, orient='index', columns=['Unique Cys-Precursors'])

# 3. Combine results
summary_df = pd.concat([cys_seq_df, cys_prec_df], axis=1)
summary_df.index.name = 'Sample'
summary_df.reset_index(inplace=True)

# Display result
summary_df.sort_values(by='Sample')
