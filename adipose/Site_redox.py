import pandas as pd
import numpy as np
import gzip
from Bio import SeqIO

# --- CONFIG ---
tsv_path = '/content/new cys test.pr_matrix (2).tsv'
fasta_path = '/content/uniprotkb_human_AND_model_organism_9606_2024_08_16.fasta.gz'
output_path = '/content/cys_summary_with_sites.xlsx'

# --- LOAD DIA-NN TSV ---
df = pd.read_csv(tsv_path, sep='\t')
metadata_cols = list(df.columns[:11])
sample_cols = list(df.columns[11:])
clean_sample_cols = [
    col.split('\\')[-1].split('/')[-1].replace('.raw', '').replace('.mzML', '')
    for col in sample_cols
]
df.rename(columns=dict(zip(sample_cols, clean_sample_cols)), inplace=True)
sample_cols = clean_sample_cols

# --- CYS MARKING ---
df['HasCys'] = df['Stripped.Sequence'].str.contains('C')
df['Peptide_Key'] = df['Stripped.Sequence'] + '_z' + df['Precursor.Charge'].astype(str)

# --- FASTA PARSING ---
protein_sequences = {}
with gzip.open(fasta_path, 'rt') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        prot_id = record.id.split('|')[1] if '|' in record.id else record.id
        protein_sequences[prot_id] = str(record.seq)

# --- MAP CYSTEINE POSITIONS ---
def map_peptide_to_sites(row):
    seq = row['Stripped.Sequence']
    prot_id = row['Protein.Ids']
    prot_seq = protein_sequences.get(prot_id)
    if not prot_seq:
        return []
    pos = prot_seq.find(seq)
    if pos == -1:
        return []
    return [pos + i + 1 for i, aa in enumerate(seq) if aa == 'C']

df['Cys_Positions'] = df.apply(map_peptide_to_sites, axis=1)
df_exploded = df.explode('Cys_Positions')
df_exploded = df_exploded[df_exploded['Cys_Positions'].notna()]
df_exploded['Site_Key'] = df_exploded['Protein.Ids'] + '_C' + df_exploded['Cys_Positions'].astype(int).astype(str)

# --- CYS SUMMARY SHEET ---
unique_cys_seq_counts = {
    sample: df[(df['HasCys']) & (df[sample].notna())]['Stripped.Sequence'].nunique()
    for sample in sample_cols
}
unique_cys_precursor_counts = {
    sample: df[(df['HasCys']) & (df[sample].notna())]['Peptide_Key'].nunique()
    for sample in sample_cols
}
summary_df = pd.DataFrame({
    'Sample': sample_cols,
    'Unique Cys-Seqs': [unique_cys_seq_counts[s] for s in sample_cols],
    'Unique Cys-Precursors': [unique_cys_precursor_counts[s] for s in sample_cols]
})

# --- LABELING EFFICIENCY ---
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
cys_df = df[df['HasCys']]
efficiency_table = pd.DataFrame(columns=['Sample', 'NEM_L', 'NEM_H', 'Unlabeled'])
for sample in sample_cols:
    sample_data = cys_df[cys_df[sample].notna()]
    counts = sample_data['LabelType'].value_counts()
    row = {
        'Sample': sample,
        'NEM_L': counts.get('NEM_L', 0),
        'NEM_H': counts.get('NEM_H', 0),
        'Unlabeled': counts.get('Unlabeled', 0)
    }
    efficiency_table = pd.concat([efficiency_table, pd.DataFrame([row])], ignore_index=True)
efficiency_table['% Labeled'] = 100 * (efficiency_table['NEM_L'] + efficiency_table['NEM_H']) / \
                                (efficiency_table['NEM_L'] + efficiency_table['NEM_H'] + efficiency_table['Unlabeled'])
efficiency_table = efficiency_table.sort_values(by='% Labeled', ascending=False).reset_index(drop=True)

# --- SITE-LEVEL REDOX ---
df_exploded['Label'] = df_exploded['Modified.Sequence'].apply(get_label_type)
df_labeled = df_exploded[df_exploded['Label'].isin(['NEM_L', 'NEM_H'])]

grouped = df_labeled.groupby(['Site_Key', 'Label'])[sample_cols].sum().reset_index()
pivot = grouped.pivot(index='Site_Key', columns='Label', values=sample_cols)
pivot.columns = [f"{sample}_{label}" for sample, label in pivot.columns]
pivot.reset_index(inplace=True)

meta_cols = ['Site_Key', 'Protein.Ids', 'Protein.Group', 'Protein.Names', 'Genes', 'First.Protein.Description']
metadata_df = df_exploded[meta_cols].drop_duplicates(subset='Site_Key')
merged = pd.merge(metadata_df, pivot, on='Site_Key', how='left')

# Calculate % oxidation
percent_oxidized = merged[meta_cols].copy()
for sample in sample_cols:
    L = merged.get(f"{sample}_NEM_L", pd.Series([np.nan] * len(merged)))
    H = merged.get(f"{sample}_NEM_H", pd.Series([np.nan] * len(merged)))
    with np.errstate(divide='ignore', invalid='ignore'):
        pct_ox = (H / (L + H)) * 100
        pct_ox[L.isna() & H.notna()] = 100
        pct_ox[H.isna() & L.notna()] = 0
        pct_ox[(L.isna()) & (H.isna())] = np.nan
    percent_oxidized[sample] = pct_ox

# --- CYS SITE COVERAGE ---
total_cys_per_protein = {pid: seq.count('C') for pid, seq in protein_sequences.items()}
observed_cys_per_protein = (
    df_exploded[df_exploded['Site_Key'].notna()]
    .groupby('Protein.Ids')['Site_Key'].nunique()
    .to_dict()
)

coverage_data = []
for prot_id, total_cys in total_cys_per_protein.items():
    observed = observed_cys_per_protein.get(prot_id, 0)
    percent = 100 * observed / total_cys if total_cys > 0 else 0
    coverage_data.append({
        'Protein.Ids': prot_id,
        'Total Cysteines': total_cys,
        'Observed Cys Sites': observed,
        'Cys Site Coverage (%)': percent
    })
coverage_df = pd.DataFrame(coverage_data)

# --- EXPORT ---
with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
    summary_df.to_excel(writer, sheet_name='Cys Summary', index=False)
    efficiency_table.to_excel(writer, sheet_name='Labeling Efficiency', index=False)
    percent_oxidized.to_excel(writer, sheet_name='Redox Site Summary', index=False)
    coverage_df.to_excel(writer, sheet_name='Cys Site Coverage', index=False)

print(f"✅ All results saved to: {output_path}")
