import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt
from goatools.obo_parser import GODag

# === Load and Prepare Data ===
df = pd.read_csv("/content/annotated_cysteine_with_biology (1).csv")
df = df.dropna(subset=["GO_Terms", "Transformation", "Site"])
df["GO_Terms"] = (
    df["GO_Terms"]
      .str.split(";")
      .apply(lambda lst: [x.strip() for x in lst])
)
df_exp = df.explode("GO_Terms")

# === Filter Significant GO Terms via χ² ===
contingency = pd.crosstab(df_exp["Transformation"], df_exp["GO_Terms"])
chi2, p, dof, expected = chi2_contingency(contingency)
residuals = (contingency - expected) / np.sqrt(expected)
significant_terms = residuals.columns[residuals.abs().max() > 2].tolist()
df_sig = df_exp[df_exp["GO_Terms"].isin(significant_terms)].copy()

# === Parse GO Ontology into a NetworkX DAG ===
# (Download the latest go-basic.obo from http://purl.obolibrary.org/obo/go/go-basic.obo)
godag = GODag("go-basic.obo")

go_dag = nx.DiGraph()
for term_id, term in godag.items():
    for parent in term.parents:
        go_dag.add_edge(parent.id, term_id)

print(f"Loaded GO DAG with {go_dag.number_of_nodes()} terms and {go_dag.number_of_edges()} edges")

# === Define Category‑Theory Functor ===
class RedoxFunctor:
    def __init__(self, df_filtered, transformation, go_dag):
        self.transformation = transformation

        # Only the rows for this transformation
        df_t = df_filtered[df_filtered["Transformation"] == transformation]

        # Objects of 𝒞 (sites for this T) and 𝒟 (GO terms in this image)
        self.sites = df_t["Site"].unique().tolist()
        self.go_terms = sorted({g for g in df_t["GO_Terms"]})

        # Object mapping: F_obj[site] = set of GO terms
        self.F_obj = defaultdict(set)
        for _, row in df_t.iterrows():
            self.F_obj[row["Site"]].add(row["GO_Terms"])

        # Morphisms in 𝒞 are just identities
        self.C_morphisms = {(s, s) for s in self.sites}
        # Morphisms in 𝒟 are the GO‐ontology edges
        self.D_morphisms = set(go_dag.edges())

    def F_objects(self):
        """Return the object‐map dict: Site → set(GO)."""
        return dict(self.F_obj)

    def F_morphisms(self):
        """
        Return the functor’s action on morphisms.
        Here only identities in C → identities in D.
        """
        fmap = {}
        for s in self.sites:
            for g in self.F_obj[s]:
                fmap[(s, s, 'id')] = (g, g, 'id')
        return fmap

    def draw_functor_graph(self):
        """Visualize the bipartite site→GO mapping on top of the GO‐DAG."""
        G = nx.DiGraph()
        G.add_edges_from(self.D_morphisms)
        for s, gos in self.F_obj.items():
            for g in gos:
                G.add_edge(s, g)

        plt.figure(figsize=(12, 12))
        pos = nx.spring_layout(G, k=0.2, iterations=100)
        nx.draw_networkx_nodes(G, pos,
                               nodelist=self.sites,
                               node_color="tomato",
                               node_size=200)
        nx.draw_networkx_nodes(G, pos,
                               nodelist=self.go_terms,
                               node_color="skyblue",
                               node_size=50)
        go_edges  = [e for e in self.D_morphisms
                     if e[0] in self.go_terms and e[1] in self.go_terms]
        map_edges = [(s, g) for s in self.sites for g in self.F_obj[s]]
        nx.draw_networkx_edges(G, pos, edgelist=go_edges,  alpha=0.3)
        nx.draw_networkx_edges(G, pos, edgelist=map_edges,
                               edge_color="black", width=1.2)
        labels = {n: n for n in (self.sites + self.go_terms)}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=6)

        plt.title(f"Functor F₍{self.transformation}₎: 𝒞→𝒟")
        plt.axis("off")

        # save at 300 dpi:
        out_path = f"/content/functor_{self.transformation}.png"
        plt.savefig(out_path, dpi=300, bbox_inches='tight')
        print(f"Saved graph to {out_path} at 300 dpi")

        plt.show()

# === Instantiate and Draw for Each Transformation ===
for T in df_sig["Transformation"].unique():
    F_T = RedoxFunctor(df_sig, T, go_dag)
    print(f"\n--- Transformation: {T} ---")
    print("Object Mapping (Site → GO‐sets):")
    for site, gos in F_T.F_objects().items():
        print(f"  {site}: {sorted(gos)}")
    F_T.draw_functor_graph()
