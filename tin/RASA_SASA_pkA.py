# === 1. SETUP ===
import os
import requests
import freesasa
import subprocess
import pandas as pd
from Bio.PDB import PDBParser

# === 2. FUNCTIONS ===
def fetch_alphafold_structure(uniprot_id):
    file_name = f"{uniprot_id}.pdb"
    if not os.path.exists(file_name):
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
        response = requests.get(url)
        if response.status_code == 200:
            with open(file_name, "w") as f:
                f.write(response.text)
        else:
            raise ValueError(f"Could not download AlphaFold model for {uniprot_id}")
    return file_name

def compute_cysteine_sasa(pdb_path):
    structure = freesasa.Structure(pdb_path)
    result = freesasa.calc(structure)
    pdb_parser = PDBParser(QUIET=True)
    structure_bio = pdb_parser.get_structure("model", pdb_path)
    
    sasa_dict = {}
    for model in structure_bio:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == "CYS":
                    res_id = str(residue.get_id()[1]).rjust(3)
                    sasa = result.residueAreas().get((chain.id, res_id, "CYS"), 0.0)
                    sasa_dict[(chain.id, res_id)] = round(sasa, 2)
    return sasa_dict

def run_propka(pdb_path):
    from propka import run
    try:
        run.single(pdb_path, optargs=["-q"])  # Quiet mode
        return pdb_path.replace(".pdb", ".pka")
    except Exception as e:
        raise RuntimeError(f"PROPKA failed on {pdb_path}: {e}")

def parse_propka_cysteine_accessibility(pka_file):
    pka_dict = {}
    with open(pka_file, "r") as f:
        for line in f:
            if line.strip().startswith("CYS"):
                parts = line.strip().split()
                if len(parts) >= 4:
                    res_name, res_num, chain_id, pred_pKa = parts[0], parts[1], parts[2], parts[3]
                    pKa_value = float(pred_pKa.rstrip('*'))
                    pka_dict[(chain_id, res_num.rjust(3))] = pKa_value
    return pka_dict

def process_cysteines(uniprot_id, target_sites):
    try:
        pdb_path = fetch_alphafold_structure(uniprot_id)
        sasa_map = compute_cysteine_sasa(pdb_path)
        propka_file = run_propka(pdb_path)
        pka_map = parse_propka_cysteine_accessibility(propka_file)

        results = []
        for site in target_sites:
            site_str = str(site).rjust(3)
            for chain in ["A"]:  # AlphaFold always uses Chain A
                key = (chain, site_str)
                results.append({
                    "Protein": uniprot_id,
                    "Chain": chain,
                    "Site": site,
                    "SASA (Å^2)": sasa_map.get(key, "NA"),
                    "Predicted pKa": pka_map.get(key, "NA")
                })
        return results
    except Exception as e:
        print(f"Error processing {uniprot_id}: {e}")
        return []

# === 3. LOAD YOUR TSV FILE ===
df_input = pd.read_csv("/content/algebraic_redox_transformation_table.tsv", sep="\t")
print("Loaded dataset shape:", df_input.shape)

# === 4. GROUP BY PROTEIN AND CYS SITE ===
grouped = df_input.groupby("Protein")["Site"].apply(list).to_dict()

# === 5. PROCESS EACH PROTEIN ===
all_results = []
for uniprot_id, sites in grouped.items():
    print(f"Processing {uniprot_id} ({len(sites)} sites)...")
    result = process_cysteines(uniprot_id, sites)
    all_results.extend(result)

# === 6. MERGE RESULTS WITH ORIGINAL TABLE ===
df_annot = pd.DataFrame(all_results)
df_merged = pd.merge(df_input, df_annot, on=["Protein", "Site"], how="left")

# === 7. DONE: Show & Save ===
df_merged.to_csv("/content/annotated_cysteine_redox_table.csv", index=False)
df_merged.head(10)
