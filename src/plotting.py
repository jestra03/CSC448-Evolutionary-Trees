import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
import seaborn as sns
import os

def plot_tree(linkage_matrix, labels, output_path="results/phylogenetic_tree.png"):
    """
    Plots and saves the phylogenetic tree.
    """
    plt.figure(figsize=(15, 10))

    dendrogram(
        linkage_matrix,
        labels=labels,
        orientation='top',
        leaf_font_size=8,
        leaf_rotation=90,
        show_leaf_counts=True,
        distance_sort='descending',
        show_contracted=True
    )

    plt.title("Phylogenetic Tree")
    plt.xlabel("Sequences")
    plt.ylabel("Distance")
    plt.grid(True, axis='y', linestyle='--', alpha=0.3)
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_heatmap(data, title, output_path):
    """
    Plot gene expression data as a heatmap.
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)  # Ensure directory exists

    plt.figure(figsize=(12, 8))
    sns.heatmap(data, cmap='RdBu_r', center=0, robust=True, vmin=-2, vmax=2)
    plt.title(title)
    plt.xlabel('Time Points')
    plt.ylabel('Genes')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()