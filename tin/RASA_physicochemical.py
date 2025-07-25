import pandas as pd
import requests

# === Load File ===
df = pd.read_csv("/content/annotated_cysteine_redox_table_sasa_fixed.csv")
print(f"✅ Loaded {len(df)} annotated cysteine entries")

# === Kyte-Doolittle Hydrophobicity Scale ===
kyte_doolittle = {
    'A': 1.8,  'C': 2.5,  'D': -3.5, 'E': -3.5, 'F': 2.8,
    'G': -0.4, 'H': -3.2, 'I': 4.5,  'K': -3.9, 'L': 3.8,
    'M': 1.9,  'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V': 4.2,  'W': -0.9, 'Y': -1.3
}
polar_residues = set('DERKHNQSTY')

# === Fetch UniProt Sequences ===
def fetch_uniprot_sequence(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        lines = response.text.splitlines()
        return ''.join(lines[1:])  # skip header
    else:
        raise ValueError(f"❌ UniProt fetch failed for {uniprot_id}")

# === Feature Functions ===
def compute_window(seq, site, window=5):
    """Returns hydrophobicity & polarity in ±window around site (1-based index)"""
    start = max(0, site - window - 1)
    end = min(len(seq), site + window)
    window_seq = seq[start:end]
    hydro = [kyte_doolittle.get(res, 0.0) for res in window_seq]
    polar = [res in polar_residues for res in window_seq]
    return round(sum(hydro)/len(hydro), 3), round(sum(polar)/len(polar), 3)

# === Compute Features Per Protein ===
hydro_scores = []
polar_scores = []
norm_positions = []
cys_counts = []

for pid in df["Protein"].unique():
    try:
        print(f"\n🧬 Processing {pid}...")
        seq = fetch_uniprot_sequence(pid)
        cys_count = seq.count("C")
        length = len(seq)

        subset = df[df["Protein"] == pid]
        for i, row in subset.iterrows():
            site = int(row["Site"])
            hydro, polar = compute_window(seq, site)
            norm_pos = round(site / length, 3)

            hydro_scores.append(hydro)
            polar_scores.append(polar)
            norm_positions.append(norm_pos)
            cys_counts.append(cys_count)

            if i % 100 == 0:
                print(f"  - {pid} Cys@{site}: hydro={hydro}, polar={polar}, norm_pos={norm_pos}, total_cys={cys_count}")

    except Exception as e:
        print(f"⚠️ {pid} failed: {e}")
        for _ in range(len(df[df["Protein"] == pid])):
            hydro_scores.append(None)
            polar_scores.append(None)
            norm_positions.append(None)
            cys_counts.append(None)

# === Attach Columns ===
df["Hydrophobicity"] = hydro_scores
df["Polarity"] = polar_scores
df["NormPosition"] = norm_positions
df["CysCount"] = cys_counts

# === Preview ===
print("\n📊 Sanity check on first 5 entries:")
print(df[["Protein", "Site", "Hydrophobicity", "Polarity", "NormPosition", "CysCount"]].head())

# === Save ===
df.to_csv("/content/annotated_cysteine_with_physicochem.csv", index=False)
print("💾 Saved to /content/annotated_cysteine_with_physicochem.csv")
