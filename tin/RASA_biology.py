import pandas as pd
import requests
import time

# === Load your dataset ===
df = pd.read_csv("/content/annotated_cysteine_with_physicochem.csv")
protein_ids = df["Protein"].dropna().unique()  # replace with your actual column

# === Function to query UniProt ===
def fetch_uniprot_info(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            location = []
            go_terms = []
            pathways = []

            # Subcellular location
            for comment in data.get("comments", []):
                if comment["commentType"] == "SUBCELLULAR LOCATION":
                    for loc in comment.get("subcellularLocations", []):
                        location.append(loc["location"]["value"])

            # GO terms
            for xref in data.get("uniProtKBCrossReferences", []):
                if xref["database"] == "GO":
                    go_terms.append(xref["id"])

            # Pathways (Reactome, KEGG)
            for xref in data.get("uniProtKBCrossReferences", []):
                if xref["database"] in ["Reactome", "KEGG"]:
                    pathways.append(xref["id"])

            return {
                "UniProtID": uniprot_id,
                "SubcellularLocation": ";".join(location),
                "GO_Terms": ";".join(go_terms),
                "Pathways": ";".join(pathways)
            }
        else:
            return {"UniProtID": uniprot_id, "SubcellularLocation": "", "GO_Terms": "", "Pathways": ""}
    except Exception as e:
        print(f"Error fetching {uniprot_id}: {e}")
        return {"UniProtID": uniprot_id, "SubcellularLocation": "", "GO_Terms": "", "Pathways": ""}

# === Query UniProt for all proteins ===
results = []
for i, uid in enumerate(protein_ids):
    results.append(fetch_uniprot_info(uid))
    time.sleep(0.5)  # be nice to the UniProt API

# === Save annotations ===
bio_df = pd.DataFrame(results)
bio_df.to_csv("/content/uniprot_annotations.csv", index=False)

import pandas as pd

# Load data
df = pd.read_csv("/content/annotated_cysteine_with_physicochem.csv")
bio_df = pd.read_csv("/content/uniprot_annotations.csv")

# Merge using correct columns
df_merged = df.merge(bio_df, left_on="Protein", right_on="UniProtID", how="left")

# Save the result
df_merged.to_csv("/content/annotated_cysteine_with_biology.csv", index=False)
print("✅ Annotation complete. Ready for enrichment analysis.")

import pandas as pd
from scipy.stats import chi2_contingency
import seaborn as sns
import matplotlib.pyplot as plt

# === Load merged dataframe with biological annotations ===
df = pd.read_csv("/content/annotated_cysteine_with_biology.csv")

# === Drop rows with missing compartment info ===
df = df.dropna(subset=["SubcellularLocation", "Transformation"])

# === Create a contingency table: Transformation × Compartment ===
compartment_counts = pd.crosstab(df["Transformation"], df["SubcellularLocation"])

# === Perform Chi-squared test ===
chi2, p, dof, expected = chi2_contingency(compartment_counts)
print(f"Chi² = {chi2:.4f}, p = {p:.4e}, dof = {dof}")

# === Compute standardized residuals (to see over/under-enrichment) ===
residuals = (compartment_counts - expected) / expected**0.5
residuals = residuals.T  # Transpose for heatmap

# === Filter for compartments with sufficient counts across classes ===
residuals = residuals.loc[compartment_counts.sum(axis=0) > 20]

# === Plot heatmap of residuals ===
plt.figure(figsize=(12, 8))
sns.heatmap(residuals, cmap="coolwarm", center=0, annot=True, fmt=".2f")
plt.title("Subcellular Compartment Enrichment by Transformation Class\n(Standardized Residuals from Chi² Test)")
plt.ylabel("Subcellular Compartment")
plt.xlabel("Transformation")
plt.tight_layout()
plt.savefig("/content/heatmap_residuals.png", dpi=300)
plt.show()
