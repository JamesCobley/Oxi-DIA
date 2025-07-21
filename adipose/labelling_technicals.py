# Helper to classify modification status
def get_label_type(modseq):
    if pd.isna(modseq):
        return 'Unlabeled'
    elif 'NEM_L' in modseq:
        return 'NEM_L'
    elif 'NEM_H' in modseq:
        return 'NEM_H'
    else:
        return 'Unlabeled'

df['LabelType'] = df['Modified.Sequence'].apply(get_label_type)
df['HasCys'] = df['Stripped.Sequence'].str.contains('C')
df['Peptide_Key'] = df['Stripped.Sequence'] + '_z' + df['Precursor.Charge'].astype(str)

# Only cysteine-containing peptides
cys_df = df[df['HasCys']]

# Initialize results
efficiency_table = pd.DataFrame(columns=['Sample', 'NEM_L', 'NEM_H', 'Unlabeled'])

for sample in sample_cols:
    sample_data = cys_df[cys_df[sample].notna()]
    counts = sample_data['LabelType'].value_counts()

    efficiency_table = pd.concat([efficiency_table, pd.DataFrame([{
        'Sample': sample,
        'NEM_L': counts.get('NEM_L', 0),
        'NEM_H': counts.get('NEM_H', 0),
        'Unlabeled': counts.get('Unlabeled', 0)
    }])], ignore_index=True)

# Calculate % labeled
efficiency_table['% Labeled'] = 100 * (efficiency_table['NEM_L'] + efficiency_table['NEM_H']) / \
                                (efficiency_table['NEM_L'] + efficiency_table['NEM_H'] + efficiency_table['Unlabeled'])

# Sort and display
efficiency_table = efficiency_table.sort_values(by='% Labeled', ascending=False)
efficiency_table.reset_index(drop=True, inplace=True)
efficiency_table
