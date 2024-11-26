import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances
from umap import UMAP

# 1. Read the TSV file
file_path = '/Users/z5439708/Desktop/Figures_gene_classifier_paper/presence_absence.tsv'
df = pd.read_csv(file_path, sep='\t', header=None, names=['Species', 'Proteins'])

# Print the first few rows to verify
print(df.head())

# 2. Create the Presence-Absence Matrix
# Split the proteins column into individual proteins
df['Proteins'] = df['Proteins'].apply(lambda x: x.split(','))

# Get all unique proteins
all_proteins = sorted(set(protein for sublist in df['Proteins'] for protein in sublist))

# Initialize a DataFrame for the presence-absence matrix
presence_absence_matrix = pd.DataFrame(0, index=df['Species'], columns=all_proteins)

# Fill the presence-absence matrix
for idx, row in df.iterrows():
    presence_absence_matrix.loc[row['Species'], row['Proteins']] = 1

# Print the presence-absence matrix
print(presence_absence_matrix.head())

# 3. Compute the distance matrix
distance_matrix = pairwise_distances(presence_absence_matrix.values, metric='jaccard')

# 4. Compute UMAP for dimensionality reduction using precomputed distance matrix
umap_model = UMAP(n_neighbors=100, n_components=2, metric='precomputed', min_dist=1, random_state=42)
umap_result = umap_model.fit_transform(distance_matrix)

# 5. Compute protein counts
df['Protein_Count'] = df['Proteins'].apply(len)

# 6. Assign colors based on protein counts
def get_color(protein_count):
    if protein_count <= 15:
        return 'red'
    elif protein_count >= 32:
        return 'green'
    else:
        return 'grey'

colors = df['Protein_Count'].apply(get_color)

# 7. Plot using Matplotlib with colors based on protein counts
fig, ax = plt.subplots(figsize=(12, 8))

# Plot with color based on protein counts
scatter = ax.scatter(umap_result[:, 0], umap_result[:, 1], c=colors, s=50, edgecolor='k', alpha=0.7)

# Customize the plot
plt.xlabel('UMAP Component 1')
plt.ylabel('UMAP Component 2')
plt.title('UMAP Plot of Species with Protein Count-Based Coloring')

# Create a custom legend
import matplotlib.lines as mlines

legend_handles = [
    mlines.Line2D([], [], color='red', marker='o', linestyle='None', markersize=10, label='≤ 15 Proteins'),
    mlines.Line2D([], [], color='green', marker='o', linestyle='None', markersize=10, label='≥ 32 Proteins'),
    mlines.Line2D([], [], color='grey', marker='o', linestyle='None', markersize=10, label='16-31 Proteins')
]

plt.legend(handles=legend_handles, loc='best')

# Save and show the plot
plt.savefig('~/Desktop/Presence_absence_protein_count.pdf')
plt.show()
