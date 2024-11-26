from Bio import AlignIO
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceCalculator
import umap
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN


# 1. Read the MSA file
alignment = AlignIO.read("/Users/z5439708/Desktop/Distance_matrix_PCA/length_filtered_hits90_query_added.aln", "fasta")

# 2. Calculate pairwise distances
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)

# Convert the distance matrix to a numpy array
dm_array = np.array(dm)

# 3. Perform UMAP dimensionality reduction on the distance matrix
umap_model = umap.UMAP(n_neighbors=200, metric='precomputed', min_dist=0.3, random_state=42, n_components=2)
umap_coords = umap_model.fit_transform(dm_array)

# 4. Perform DBSCAN clustering on the UMAP projection
dbscan_model = DBSCAN(eps=10, min_samples=5, metric='euclidean')
dbscan_labels = dbscan_model.fit_predict(umap_coords)

# 5. Read the specific IDs to highlight from a text file
highlight_ids = set()
with open("/Users/z5439708/Desktop/Distance_matrix_PCA/fliG-kegg_id.txt", "r") as file:
    for line in file:
        highlight_ids.add(line.strip())


# 6. Find the indices of the specific IDs to highlight
highlight_indices = [i for i, record in enumerate(alignment) if record.id in highlight_ids]

# 7. Get the coordinates of the specific points to highlight
highlight_points = umap_coords[highlight_indices]

# 8. Plot the clustered points with different colors
plt.figure(figsize=(10, 8))

# Assign colors to clusters
unique_labels = set(dbscan_labels)
colors = plt.cm.nipy_spectral(np.linspace(0, 1, len(unique_labels)))

for label, color in zip(unique_labels, colors):
    if label == -1:
        # Noise points
        col = 'black'
    else:
        # Clustered points
        col = color
    class_member_mask = (dbscan_labels == label)
    xy = umap_coords[class_member_mask]
    plt.scatter(xy[:, 0], xy[:, 1], color=col, edgecolor='k', s=50, label=f'Cluster {label}' if label != -1 else 'Noise')

# Highlight specific IDs
plt.scatter(highlight_points[:, 0], highlight_points[:, 1], color='orange', marker='o', edgecolor='k', s=100, label='Highlighted IDs')

plt.xlabel('UMAP Dimension 1')
plt.ylabel('UMAP Dimension 2')
plt.title('DBSCAN Clustering on UMAP Projection with Highlighted Points')
plt.legend()
plt.show()

