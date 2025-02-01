import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
import os
from plotting import plot_heatmap


def load_gene_data(filepath):
    """
    Load and preprocess gene expression data.
    """
    data = pd.read_csv(filepath, sep='\t', index_col=0)
    ratio_columns = ['R1.Ratio', 'R2.Ratio', 'R3.Ratio', 'R4.Ratio',
                     'R5.Ratio', 'R6.Ratio', 'R7.Ratio']
    data = data[ratio_columns]

    if data.min().min() < 0:
        return data  # Already log-transformed
    else:
        epsilon = 1e-10
        return np.log2(data.clip(lower=epsilon))


def cluster_genes(data, n_clusters=10, method='average', metric='correlation'):
    """
    Cluster genes based on expression pattern similarity.
    """
    # Calculate distance matrix
    distances = pdist(data, metric=metric)

    # Perform hierarchical clustering
    linkage_matrix = linkage(distances, method=method)
    clusters = fcluster(linkage_matrix, n_clusters, criterion='maxclust')

    # Calculate silhouette score
    dist_matrix = squareform(distances)
    sil_score = silhouette_score(dist_matrix, clusters)

    # Sort genes by cluster for visualization
    cluster_order = pd.Series(clusters, index=data.index).sort_values()
    sorted_data = data.loc[cluster_order.index]

    plot_heatmap(sorted_data, f'Clustered Gene Expression (k={n_clusters}, {method}, {metric})',
                 "results/clustered_heatmap.png")

    return clusters, sil_score


def top_variable_genes(data, n_genes=230):
    """
    Identify genes with largest expression changes during diauxic shift.
    Considers both magnitude and pattern of change.
    """
    # Split data into pre-shift (R1-R3) and post-shift (R5-R7)
    pre_shift = data[['R1.Ratio', 'R2.Ratio', 'R3.Ratio']]
    post_shift = data[['R5.Ratio', 'R6.Ratio', 'R7.Ratio']]

    # Calculate mean expression levels
    pre_mean = pre_shift.mean(axis=1)
    post_mean = post_shift.mean(axis=1)

    # Calculate absolute log fold change
    log_fold_change = np.abs(post_mean - pre_mean)

    # Calculate variance to capture pattern changes
    total_var = data.var(axis=1)

    # Combine metrics (both magnitude and pattern of change)
    combined_score = log_fold_change * np.sqrt(total_var)

    return combined_score.nlargest(n_genes).index


def compare_gene_sets(gene_set1, gene_set2):
    """
    Compute overlap and Jaccard similarity.
    """
    set1 = set(gene_set1)
    set2 = set(gene_set2.index if isinstance(gene_set2, pd.DataFrame) else gene_set2)
    overlap = len(set1.intersection(set2))
    jaccard = overlap / len(set1.union(set2))
    return overlap, jaccard
