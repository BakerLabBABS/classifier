from Bio import AlignIO
import matplotlib.pyplot as plt
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceCalculator
from sklearn.decomposition import PCA


# 1. Read the MSA file
alignment = AlignIO.read("/Users/z5439708/Downloads/Distance_mat/fliC_homologs_0.9.aln1", "fasta")

# 2. Calculate pairwise distances
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)

# 3. Perform PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(dm)

# 4. Plot PCA results
### Plot centre
center = pca_result.mean(axis=0)
## Plot PCA results
#plt.figure(figsize=(8, 6))
#plt.scatter(pca_result[:, 0], pca_result[:, 1])
#plt.scatter(center[0], center[1], color='red', marker='o', label='Center')
#plt.xlabel('PC1')
#plt.ylabel('PC2')
#plt.title('PCA Plot')
#plt.legend()
#plt.show()

### 5. To get distant points

# Calculate the Euclidean distance between each point and the center
distances = np.linalg.norm(pca_result - center, axis=1)

# Find the indices of the points with maximum distance
max_distance_indices = np.argsort(distances)[-10:]  # Select the top 10 most distant points

# Get the coordinates of the most distant points
most_distant_points = pca_result[max_distance_indices]

# Print the coordinates of the most distant points
print("Most distant points from center:")
for i, point in enumerate(most_distant_points):
    print(f"Point {i+1}: {point}")

### 6. Plot the distant points

# Plot PCA results
plt.figure(figsize=(8, 6))
plt.scatter(pca_result[:, 0], pca_result[:, 1], label='All Points')
plt.scatter(center[0], center[1], color='red', marker='o', label='Center')

# Mark the most distant points
plt.scatter(most_distant_points[:, 0], most_distant_points[:, 1], color='green', marker='x', label='Most Distant Points')

plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('PCA Plot')
plt.legend()
#plt.savefig('/Users/z5439708/Downloads/pca_plot.png')
plt.show()

#### 7. Find the IDs of the distant points

# Calculate the Euclidean distance between each point and the center
distances = np.linalg.norm(pca_result - center, axis=1)

# Find the indices of the points with maximum distance
max_distance_indices = np.argsort(distances)[-10:]  # Select the top 10 most distant points

# Get the IDs of the most distant points from the original dataset
most_distant_ids = [alignment[int(idx)].id for idx in max_distance_indices]

# Print the IDs of the most distant points
print("IDs of the most distant points:")
for i, point_id in enumerate(most_distant_ids):
    print(f"Point {i+1} ID: {point_id}")
