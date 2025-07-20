import pandas as pd

# Load full DIA-NN report
df = pd.read_parquet('/content/drive/MyDrive/reportnewcal.parquet')

# Only peptides with C
df['Has_C'] = df['Stripped.Sequence'].str.contains("C")
df_cys = df[df['Has_C']].copy()

# Deduplicate per Run and Stripped.Sequence
df_unique = df_cys.drop_duplicates(subset=['Run', 'Stripped.Sequence'])

# Assign label type
def get_label(seq):
    if "_L" in seq:
        return "L"
    elif "_H" in seq:
        return "H"
    else:
        return "None"

df_unique['Label'] = df_unique['Modified.Sequence'].apply(get_label)

# Summarize per run
summary = df_unique.groupby('Run').agg(
    Total_C_Peptides=('Stripped.Sequence', 'count'),
    Light_C=('Label', lambda x: (x == 'L').sum()),
    Heavy_C=('Label', lambda x: (x == 'H').sum())
).reset_index()

summary['Unlabeled_C'] = summary['Total_C_Peptides'] - summary['Light_C'] - summary['Heavy_C']
summary['%Unlabeled_C'] = (summary['Unlabeled_C'] / summary['Total_C_Peptides']) * 100

# Filter to mixing runs (exclude pure 0 and 100)
mixes_only = summary['Run'].str.contains("James_")
mixes_only = summary[mixes_only]

# Display result
print(mixes_only.to_string(index=False))
