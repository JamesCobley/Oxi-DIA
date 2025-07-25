import os
import pandas as pd
import requests
import freesasa

# === 1. Load existing annotated file ===
df = pd.read_csv("/content/annotated_cysteine_redox_table.csv")
print(f"✅ Loaded {len(df)} entries")

# === 2. Download AlphaFold PDB if needed ===
def fetch_af_structure(uniprot_id):
    base_url = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_{}.pdb"
    for ver in ["v4", "v3", "v2", "v1"]:
        url = base_url.format(uniprot_id, ver)
        response = requests.get(url)
        if response.status_code == 200:
            filename = f"{uniprot_id}_AF_model.pdb"
            with open(filename, "wb") as f:
                f.write(response.content)
            print(f"🔽 Downloaded AlphaFold {ver} model for {uniprot_id}")
            return filename
    raise FileNotFoundError(f"No AlphaFold model found for {uniprot_id}")

# === 3. Compute per-CYS SASA ===
def get_cys_sasa(pdb_path):
    structure = freesasa.Structure(pdb_path)
    result = freesasa.calc(structure)
    res_areas = result.residueAreas()

    sasa = {}
    for chain_id, residues in res_areas.items():
        for resnum, area in residues.items():
            if area.residueType == "CYS":
                sasa_val = round(area.total, 2)
                sasa[(chain_id, str(resnum).rjust(3))] = sasa_val
    return sasa

# === 4. Fill missing/zero SASA values ===
df["SASA (Å^2)"] = df["SASA (Å^2)"].fillna(0.0)
proteins = df["Protein"].unique()

for protein in proteins:
    print(f"\n🔍 Processing {protein}...")
    try:
        pdb_file = fetch_af_structure(protein)
        sasa_dict = get_cys_sasa(pdb_file)
        print(f"  🧬 Found {len(sasa_dict)} cysteines with SASA")

        updates = []
        for i, row in df[df["Protein"] == protein].iterrows():
            key = (row.get("Chain", "A"), str(row["Site"]).rjust(3))
            if key in sasa_dict:
                new_sasa = sasa_dict[key]
                df.at[i, "SASA (Å^2)"] = new_sasa
                updates.append((key, new_sasa))

        # Print a few updated entries for verification
        print(f"  ✅ Updated {len(updates)} residues. Examples:")
        for k, v in updates[:3]:
            print(f"     - Residue {k} → SASA = {v:.2f} Å²")

    except Exception as e:
        print(f"  ⚠️ Error with {protein}: {e}")

# === 5. Save result ===
output_path = "/content/annotated_cysteine_redox_table_sasa_fixed.csv"
df.to_csv(output_path, index=False)
print(f"\n💾 Done! Updated file saved to: {output_path}")
