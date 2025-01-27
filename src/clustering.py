import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform


def generate_phylogenetic_tree(similarity_matrix):
    """
    Generates a phylogenetic tree using UPGMA clustering.

    Args:
        similarity_matrix: Matrix of sequence similarities

    Returns:
        numpy.ndarray: Linkage matrix for plotting dendrogram
    """
    # convert similarity to distance
    distance_matrix = 1 - similarity_matrix

    # ensure matrix is symmetric and has zeros on diagonal
    distance_matrix = (distance_matrix + distance_matrix.T) / 2
    np.fill_diagonal(distance_matrix, 0)

    # verify matrix contains valid distances
    if np.any(distance_matrix < 0) or np.any(distance_matrix > 1):
        raise ValueError("Distance matrix contains invalid values")

    # convert to condensed form
    condensed_dist = squareform(distance_matrix)

    # perform UPGMA clustering
    return linkage(condensed_dist, method='average', optimal_ordering=True)


def find_closest_farthest_pairs(similarity_matrix):
    """
    Finds the closest and farthest pairs of sequences.

    Args:
        similarity_matrix: numpy array of pairwise similarities

    Returns:
        tuple: ((i1, j1, closest_score), (i2, j2, farthest_score))
    """
    n = len(similarity_matrix)
    closest_score = float('-inf')
    farthest_score = float('inf')
    closest_pair = None
    farthest_pair = None

    # search upper triangle of matrix
    for i in range(n):
        for j in range(i + 1, n):
            score = similarity_matrix[i, j]
            if score > closest_score:
                closest_score = score
                closest_pair = (i, j)
            if score < farthest_score:
                farthest_score = score
                farthest_pair = (i, j)

    return (closest_pair + (closest_score,), farthest_pair + (farthest_score,))