import pandas as pd

# === Load annotated biology dataframe ===
bio_path = '/content/annotated_cysteine_with_biology (1) (1).csv'
bio_df = pd.read_csv(bio_path)

# === Load Parquet file ===
parquet_path = '/content/drive/MyDrive/report (2).parquet'  # adjust as needed
df_parquet = pd.read_parquet(parquet_path)

# === Sample map (raw paths from TSV) ===
sample_map_raw = {
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Air_1_new_S1': 'Sample_1',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Air_1_new_S13': 'Sample_1',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Air_1_new_S7': 'Sample_1',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Air_2_new_S15': 'Sample_2',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Air_2_new_S3': 'Sample_2',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Air_2_new_S9': 'Sample_2',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Air_3_new_S11': 'Sample_3',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Air_3_new_S17': 'Sample_3',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Air_3_new_S5': 'Sample_3',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Tin_4_new_S14': 'Sample_4',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Tin_4_new_S2': 'Sample_4',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Tin_4_new_S8': 'Sample_4',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Tin_5_new_S10': 'Sample_5',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Tin_5_new_S16': 'Sample_5',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Tin_5_new_S4': 'Sample_5',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Tin_6_new_S18': 'Sample_6',
    'W:\\H_James\\Astral-L\\June 2025\\James_Mouse_Tin_6_new_S6': 'Sample_6'
}

# === Normalize the sample map keys ===
clean_sample_map = {
    path.split('\\')[-1]: sample
    for path, sample in sample_map_raw.items()
}

# === Assign sample names ===
df_parquet['Sample'] = df_parquet['Run'].map(clean_sample_map)

# === Filter and calculate mean LFQ per protein group ===
df_filtered = df_parquet.dropna(subset=['Sample', 'PG.MaxLFQ'])
mean_lfq = (
    df_filtered
    .groupby('Protein.Group')['PG.MaxLFQ']
    .mean()
    .reset_index()
    .rename(columns={'PG.MaxLFQ': 'Mean_PG_MaxLFQ'})
)

# === Merge only where protein matches ===
if 'Protein.Group' in bio_df.columns:
    merge_key = 'Protein.Group'
elif 'Protein' in bio_df.columns:
    merge_key = 'Protein'
    mean_lfq = mean_lfq.rename(columns={'Protein.Group': 'Protein'})

# Merge with left join to retain only matched entries
bio_df_with_lfq = bio_df.merge(mean_lfq, how='left', on=merge_key)

# === Save to new file ===
output_path = '/content/annotated_with_mean_LFQ.csv'
bio_df_with_lfq.to_csv(output_path, index=False)
print(f"✅ Mean LFQ values added and saved to: {output_path}")
