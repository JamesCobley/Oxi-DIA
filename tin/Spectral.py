import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csgraph
from scipy.linalg import eigh
import matplotlib.pyplot as plt

# === Load the redox data ===
df = pd.read_csv('/content/redox_sites.tsv', sep='\t')

# === Define site ID (if needed) ===
df['Site_ID'] = df['Protein'] + "_" + df['Residue'].astype(str)

# === Define sample columns explicitly ===
air_cols = ['Sample_1_%Oxidized', 'Sample_2_%Oxidized', 'Sample_3_%Oxidized']
tin_cols = ['Sample_4_%Oxidized', 'Sample_5_%Oxidized', 'Sample_6_%Oxidized']

# Filter valid rows
df = df.dropna(subset=pct_cols).copy()
df["Air_Mean"] = df[air_cols].mean(axis=1)
df["Tin_Mean"] = df[tin_cols].mean(axis=1)

# Sample to keep memory use low
df = df.sample(n=12606, random_state=1).reset_index(drop=True)

# Create 2D coordinates in a circle
N = len(df)
theta = np.linspace(0, 2 * np.pi, N, endpoint=False)
coords = np.column_stack((np.cos(theta), np.sin(theta)))
z_air = df["Air_Mean"].values
z_tin = df["Tin_Mean"].values
coords_air = np.column_stack((coords, z_air))
coords_tin = np.column_stack((coords, z_tin))

# Build kNN graph
def build_knn(coords, k=5):
    nn = NearestNeighbors(n_neighbors=k+1).fit(coords)
    G = nn.kneighbors_graph(coords, mode='connectivity')
    return G.toarray()

G_air = build_knn(coords_air)
G_tin = build_knn(coords_tin)

# Laplacians
L_air = csgraph.laplacian(G_air)
L_tin = csgraph.laplacian(G_tin)

# Dirichlet energies
E_air = z_air @ L_air @ z_air
E_tin = z_tin @ L_tin @ z_tin

# Spectral decomposition
eigvals_air = eigh(L_air, eigvals_only=True)
eigvals_tin = eigh(L_tin, eigvals_only=True)

# Plot spectra
plt.figure(figsize=(10, 5))
plt.plot(eigvals_air, label="Air")
plt.plot(eigvals_tin, label="Tin")
plt.title("Spectral Eigenmodes (Air vs Tin)")
plt.xlabel("Mode Index")
plt.ylabel("Eigenvalue")
plt.legend()
plt.tight_layout()
plt.savefig('/content/spectral_eigenmodes.png', dpi=300)
plt.show()

print("Dirichlet Energy (Air):", E_air)
print("Dirichlet Energy (Tin):", E_tin)

def spectral_entropy(eigvals):
    eigvals = eigvals[eigvals > 0]  # Avoid zero or negative eigenvalues
    eigvals = eigvals / eigvals.sum()
    return -np.sum(eigvals * np.log(eigvals))

entropy_air = spectral_entropy(eigvals_air)
entropy_tin = spectral_entropy(eigvals_tin)
print("Spectral Entropy (Air):", entropy_air)
print("Spectral Entropy (Tin):", entropy_tin)

gap_air = eigvals_air[1] - eigvals_air[0]
gap_tin = eigvals_tin[1] - eigvals_tin[0]
print("Spectral Gap (Air):", gap_air)
print("Spectral Gap (Tin):", gap_tin)

def laplacian_energy(eigvals):
    mean = np.mean(eigvals)
    return np.sum((eigvals - mean)**2)

lap_energy_air = laplacian_energy(eigvals_air)
lap_energy_tin = laplacian_energy(eigvals_tin)
print("Laplacian Energy (Air):", lap_energy_air)
print("Laplacian Energy (Tin):", lap_energy_tin)

def morse_energy(coords, z, k=5):
    nn = NearestNeighbors(n_neighbors=k+1).fit(coords)
    distances, indices = nn.kneighbors(coords)
    
    energy = 0.0
    for i in range(len(coords)):
        for j in indices[i][1:]:  # skip self (index 0)
            energy += (z[i] - z[j])**2
    return energy / 2  # each pair counted twice

morse_air = morse_energy(coords_air, z_air)
morse_tin = morse_energy(coords_tin, z_tin)

print("Morse Energy (Air):", morse_air)
print("Morse Energy (Tin):", morse_tin)
