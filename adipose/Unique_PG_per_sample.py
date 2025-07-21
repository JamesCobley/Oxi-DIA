import pandas as pd
import re

# Load DIA-NN output
file_path = '/content/new cys test.pr_matrix (2).tsv'
df = pd.read_csv(file_path, sep='\t')

# Identify metadata and sample columns
metadata_cols = ['Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes']
sample_cols = [col for col in df.columns if col not in metadata_cols and not col.startswith('Unnamed')]

# Extract biological sample number from column names (e.g., 83 from James_Twin_83_S3.raw)
def extract_sample_number(colname):
    match = re.search(r'Twin_(\d+)_S\d+\.raw$', colname)
    return match.group(1) if match else None

sample_number_map = {col: extract_sample_number(col) for col in sample_cols if extract_sample_number(col)}

# Map biological sample numbers to their corresponding columns (technicals)
from collections import defaultdict
sample_to_columns = defaultdict(list)
for col, num in sample_number_map.items():
    sample_to_columns[num].append(col)

# Count unique protein groups per biological sample
bio_sample_protein_groups = {}
for sample_num, cols in sample_to_columns.items():
    # Keep rows with a detection in any of the technical replicates
    detected = df[df[cols].notna().any(axis=1)]
    count = detected['Protein.Group'].nunique()
    bio_sample_protein_groups[sample_num] = count

# Convert to DataFrame
counts_df = pd.DataFrame.from_dict(bio_sample_protein_groups, orient='index', columns=['Unique_Protein_Groups'])
counts_df.index.name = 'Biological_Sample'
counts_df.reset_index(inplace=True)

# Calculate summary stats
mean_pg = counts_df['Unique_Protein_Groups'].mean()
std_pg = counts_df['Unique_Protein_Groups'].std()
cv_pg = (std_pg / mean_pg) * 100

# Output
print("✅ Unique Protein Groups per Biological Sample:")
print(counts_df.sort_values(by='Biological_Sample'))

print("\n📊 Summary Statistics:")
print(f"Mean: {mean_pg:.0f}")
print(f"SD: {std_pg:.0f}")
print(f"CV: {cv_pg:.1f}%")

# Optional: Save to Excel
counts_df.to_excel('/content/unique_protein_groups_per_sample.xlsx', index=False)
