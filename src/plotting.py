import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
import numpy as np

def plot_tree(linkage_matrix, labels, output_path="results/phylogenetic_tree.png"):
    """
    Plots and saves the phylogenetic tree.

    Args:
        linkage_matrix: Output from scipy.cluster.hierarchy.linkage
        labels: List of sequence labels
        output_path: Where to save the plot
    """
    plt.figure(figsize=(15, 10))

    # calculate optimal leaf rotation based on number of sequences
    n_seqs = len(labels)
    leaf_rotation = min(90, max(0, n_seqs * 0.5))

    # plot dendrogram
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

    # add grid for better distance reference
    plt.grid(True, axis='y', linestyle='--', alpha=0.3)

    # adjust layout and save
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()