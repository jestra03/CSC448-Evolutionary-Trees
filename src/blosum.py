import numpy as np

def load_blosum62(filepath="blosum62.txt"):
    """
    Loads BLOSUM62 matrix into a numpy array for faster access.
    Returns matrix and amino acid index mapping.
    """
    with open(filepath) as f:
        lines = [line.strip() for line in f if not line.startswith("#")]

    # get amino acids from header
    aa_list = lines[0].split()
    n_aa = len(aa_list)

    # create amino acid to index mapping
    aa_to_idx = {aa: idx for idx, aa in enumerate(aa_list)}

    # initialize matrix
    blosum = np.zeros((n_aa, n_aa), dtype=np.int8)

    # parse matrix
    for i, line in enumerate(lines[1:]):
        scores = line.split()[1:]  # Skip AA label
        for j, score in enumerate(scores):
            blosum[i, j] = int(score)

    return blosum, aa_to_idx


def get_score(blosum_matrix, aa_to_idx, aa1, aa2):
    """Get BLOSUM62 score for amino acid pair."""
    try:
        i = aa_to_idx[aa1]
        j = aa_to_idx[aa2]
        return blosum_matrix[i, j]
    except KeyError:
        return -4  # Default score for unknown amino acids