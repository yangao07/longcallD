import numpy as np
import matplotlib.pyplot as plt
from sklearn_extra.cluster import KMedoids
import sys

# input data: each line consists of 16 values
# vectors = np.random.randint(0, 100, size=(99, 12))  # Replace with actual data
in_fn = sys.argv[1]
vectors = np.loadtxt(in_fn, delimiter='\t')
# initial_medoids = vectors[[9, 15, 24]]  # Replace with actual initial medoid vectors if different
initial_medoids = vectors[[1, 2, 3]]  # Replace with actual initial medoid vectors if different

optimal_k = 2

# kmedoids = KMedoids(n_clusters=optimal_k, metric='manhattan', init=initial_medoids, random_state=42)
kmedoids = KMedoids(n_clusters=optimal_k, metric='manhattan', random_state=42)
kmedoids.fit(vectors)

cluster_labels = kmedoids.labels_ 
medoid_indices = kmedoids.medoid_indices_

# Step 4: Plot all vectors and highlight the representative vectors
plt.figure(figsize=(10, 6))
for i, vector in enumerate(vectors):
    if i in medoid_indices:
        # Highlight representative vectors with distinct color and size
        plt.plot(range(1, 21), vector, marker='o', color='red', label=f'Mediod Vector {i+1}', linewidth=2, markersize=8)
    else:
        # Plot other vectors in a lighter color
        plt.plot(range(1, 21), vector, marker='o', color='gray', alpha=0.7, linewidth=1, markersize=5)

plt.xlabel('Dimension')
plt.ylabel('Value')

plt.title('20-Dimensional Vectors for X Clusters')
plt.legend(loc='upper right')
plt.grid(True)

# Save the figure to a file
plt.savefig("kmedoids_clusters.png", format='png', dpi=300)  # High resolution