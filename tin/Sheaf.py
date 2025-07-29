import pandas as pd
import numpy as np
import itertools
from scipy.stats import chi2_contingency

# === Load and Prepare Data ===
df = pd.read_csv("/content/annotated_cysteine_with_biology (1).csv")
df = df.dropna(subset=["GO_Terms", "Transformation", "Site"])
df["GO_Terms"] = (
    df["GO_Terms"]
      .str.split(";")
      .apply(lambda lst: [x.strip() for x in lst])
)
df_exp = df.explode("GO_Terms")

# === 1. χ²‑filter to get only significant GO terms ===
contingency = pd.crosstab(df_exp["Transformation"], df_exp["GO_Terms"])
chi2, p, dof, expected = chi2_contingency(contingency)
residuals = (contingency - expected) / np.sqrt(expected)
significant_terms = residuals.columns[residuals.abs().max(axis=0) > 2].tolist()

# Subset to only those significant GO annotations
df_sig = df_exp[df_exp["GO_Terms"].isin(significant_terms)].copy()

# === 2. Build the cover opens U_T ===
transformations = df_sig["Transformation"].unique().tolist()
opens = {
    T: set(df_sig.loc[df_sig["Transformation"] == T, "Site"].unique())
    for T in transformations
}

# === 3. Build the presheaf F on each open ===
# F(U_T) = union of GO‑terms over sites in U_T
site2gos = df_sig.groupby("Site")["GO_Terms"].apply(set).to_dict()

C0 = {}
for T, sites in opens.items():
    gos_union = set()
    for s in sites:
        gos_union |= site2gos.get(s, set())
    C0[T] = gos_union

# === 4. Build degree‑1 cochains on pairwise overlaps ===
pairs = list(itertools.combinations(transformations, 2))
C1 = {}
for T1, T2 in pairs:
    overlap_sites = opens[T1] & opens[T2]
    gos_overlap = set()
    for s in overlap_sites:
        gos_overlap |= site2gos.get(s, set())
    C1[(T1, T2)] = gos_overlap

# === 5. Define the Čech coboundary δ^0: C^0 → C^1 ===
delta0 = {}
for T1, T2 in pairs:
    A1, A2 = C0[T1], C0[T2]
    A12 = C1[(T1, T2)]
    # δ(a)_{T1,T2} = restriction of A2 minus restriction of A1
    delta0[(T1, T2)] = (A2 & A12) - (A1 & A12)

# === 6. Compute H^0 = ⋂_{T} C^0[T] ===
H0 = set.intersection(*C0.values())

# === 7. Compute H^1 = C^1 / im δ^0 ===
all_C1 = set().union(*C1.values())
im_delta0 = set().union(*delta0.values())
H1 = all_C1 - im_delta0

# === 8. Report results ===
print("=== Sheaf Cohomology of Filtered Redox Functor Cover ===")
print(f"Transformations: {transformations}")
print(f"|C^0| (per T): {{ {', '.join(f'{T}: {len(v)}' for T,v in C0.items())} }}")
print(f"|C^1| (per overlap): {{ {', '.join(f'{pair}: {len(v)}' for pair,v in C1.items())} }}")
print(f"H^0 (common GO terms): {len(H0)}")
print(f"H^1 (obstructed GO terms): {len(H1)}")
print("\nSample H^1 generators:", list(H1)[:10])

import pandas as pd
import numpy as np
import itertools
from scipy.stats import chi2_contingency

# === 1. Load & explode the GO terms ===
df = pd.read_csv("annotated_cysteine_with_biology (1).csv", sep=",")
df = df.dropna(subset=["GO_Terms", "Transformation", "Site", "Protein", "Residue"])
df["GO_Terms"] = df["GO_Terms"].str.split(";").apply(lambda L: [x.strip() for x in L])
df = df.explode("GO_Terms")

# === 2. χ²‑filter to keep only significant GO terms ===
cont = pd.crosstab(df["Transformation"], df["GO_Terms"])
chi2, p, dof, exp = chi2_contingency(cont)
resid = (cont.values - exp) / np.sqrt(exp)
resid_df = pd.DataFrame(resid, index=cont.index, columns=cont.columns)
sig_terms = resid_df.columns[resid_df.abs().max(axis=0) > 2]
df_sig = df[df["GO_Terms"].isin(sig_terms)].copy()

# === 3. Build the “cover” opens and site→GO map ===
transformations = df_sig["Transformation"].unique().tolist()
opens = {T: set(df_sig[df_sig["Transformation"]==T]["Site"]) for T in transformations}
site2gos = df_sig.groupby("Site")["GO_Terms"].apply(set).to_dict()

# === 4. Compute H^1 via Čech coboundary ===
# 4a. C^0
C0 = {T: set().union(*(site2gos[s] for s in opens[T])) for T in transformations}
# 4b. C^1 on pairwise overlaps
pairs = list(itertools.combinations(transformations, 2))
C1 = {pair: set().union(*(site2gos[s] for s in opens[pair[0]] & opens[pair[1]]))
      for pair in pairs}
# 4c. δ^0-image
delta0 = {pair: (C0[pair[1]] & C1[pair]) - (C0[pair[0]] & C1[pair])
          for pair in pairs}
# 4d. H^1 = (∪ C1) \ (∪ im δ^0)
H1 = set().union(*C1.values()) - set().union(*delta0.values())

# === 5. Pick top 3 most‐obstructed GO terms ===
freq = {go: sum(go in gos for gos in C1.values()) for go in H1}
top_terms = sorted(freq, key=lambda g: freq[g], reverse=True)[:3]

# === 6. Build example‐site table for each top GO term ===
records = []
for go in top_terms:
    for T1, T2 in pairs:
        overlap_sites = opens[T1] & opens[T2]
        # only those sites carrying this GO term
        hits = [s for s in overlap_sites if go in site2gos.get(s, set())]
        records.append({
            "GO_ID": go,
            "Overlap": f"{T1}∩{T2}",
            "Example_Sites": hits[:5]
        })
df_examples = pd.DataFrame(records)

# === 7. Attach Protein_C<Residue> labels ===
site_info = df_sig[["Site","Protein"]].drop_duplicates().copy()
site_info["Label"] = site_info["Protein"] + "_C" + site_info["Site"].astype(str)

df_labeled = (
    df_examples
    .explode("Example_Sites")
    .merge(site_info, left_on="Example_Sites", right_on="Site", how="left")
    .groupby(["GO_ID","Overlap"])["Label"]
    .apply(list)
    .reset_index()
    .rename(columns={"Label":"Example_Labels"})
)

# === 8. Show / export ===
print(df_labeled.to_string(index=False))
# Optionally save to CSV:
df_labeled.to_csv("H1_top_GO_examples.csv", index=False)
