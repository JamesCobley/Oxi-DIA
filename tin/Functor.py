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
df["GO_Terms"] = df["GO_Terms"].str.split(";").apply(lambda lst: [x.strip() for x in lst])
df_exp = df.explode("GO_Terms")

# === Filter Significant GO Terms via χ² ===
contingency = pd.crosstab(df_exp["Transformation"], df_exp["GO_Terms"])
chi2, p, dof, expected = chi2_contingency(contingency)
residuals = (contingency - expected) / np.sqrt(expected)
all_sig_terms = residuals.columns[residuals.abs().max() > 2].tolist()

# === Parse GO Ontology into a NetworkX DAG ===
godag = GODag("go-basic.obo")
go_dag = nx.DiGraph((parent.id, child_id)
                    for child_id, term in godag.items()
                    for parent in term.parents)

print(f"Loaded full GO DAG with {go_dag.number_of_nodes()} nodes and {go_dag.number_of_edges()} edges")

# === Filter significant terms to those actually in the ontology ===
valid_sig = [t for t in all_sig_terms if t in go_dag.nodes()]
missing = set(all_sig_terms) - set(valid_sig)
if missing:
    print(f"Warning: {len(missing)} significant GO IDs not found in ontology, dropping:\n", missing)
df_sig = df_exp[df_exp["GO_Terms"].isin(valid_sig)].copy()

# === Prune to only valid_sig terms + their ancestors ===
prune_set = set(valid_sig)
for term in valid_sig:
    prune_set |= nx.ancestors(go_dag, term)

go_sub = go_dag.subgraph(prune_set).copy()
print(f"Pruned GO DAG to {go_sub.number_of_nodes()} nodes and {go_sub.number_of_edges()} edges")

# === Map GO IDs to names for labeling ===
go_labels = {tid: godag[tid].name for tid in go_sub.nodes()}

# === Precompute layout on pruned GO subgraph ===
pos_go = nx.spring_layout(go_sub, k=0.2, iterations=50, seed=0)

# === Functor class with efficient plotting ===
class RedoxFunctor:
    def __init__(self, df_filtered, transformation, go_sub, go_labels):
        self.transformation = transformation
        df_t = df_filtered[df_filtered["Transformation"] == transformation]
        self.sites = df_t["Site"].unique().tolist()
        self.go_terms = sorted({g for g in df_t["GO_Terms"]})
        self.F_obj = defaultdict(set)
        for _, row in df_t.iterrows():
            self.F_obj[row["Site"]].add(row["GO_Terms"])
        self.D_morphisms = set(go_sub.edges())

    def draw_functor_graph(self, out_dir="/content"):
        G = nx.DiGraph()
        G.add_edges_from(self.D_morphisms)
        for s, gos in self.F_obj.items():
            for g in gos:
                G.add_edge(s, g)

        # combine positions
        pos = dict(pos_go)
        if self.sites:
            site_pos = nx.circular_layout(self.sites, scale=1.0)
            pos.update(site_pos)

        plt.figure(figsize=(12,12))
        # GO nodes
        nx.draw_networkx_nodes(G, pos,
            nodelist=self.go_terms,
            node_color="skyblue", node_size=50)
        # Site nodes
        nx.draw_networkx_nodes(G, pos,
            nodelist=self.sites,
            node_color="tomato", node_size=150)
        # Edges
        go_edges = [e for e in self.D_morphisms if e[0] in self.go_terms]
        map_edges = [(s,g) for s in self.sites for g in self.F_obj[s]]
        nx.draw_networkx_edges(G, pos, edgelist=go_edges, alpha=0.3)
        nx.draw_networkx_edges(G, pos, edgelist=map_edges, edge_color="black", width=1)

        # Labels: GO by name, sites by ID
        labels = {**{g: go_labels[g] for g in self.go_terms}, **{s: s for s in self.sites}}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=6)

        plt.title(f"Functor F₍{self.transformation}₎: Sites→GO Processes")
        plt.axis("off")
        out_path = f"{out_dir}/functor_{self.transformation}.png"
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"Saved graph: {out_path}")
        plt.close()

# === Run and save each graph ===
for T in df_sig["Transformation"].unique():
    print(f"\n=== {T} ===")
    F = RedoxFunctor(df_sig, T, go_sub, go_labels)
    F.draw_functor_graph()
