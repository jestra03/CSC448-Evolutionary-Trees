import numpy as np
from blosum import load_blosum62, get_score

def smith_waterman(seq1, seq2, gap_penalty=-12):
    """
    Implements Smith-Waterman algorithm for local sequence alignment using BLOSUM62.

    Args:
        seq1, seq2: Protein sequences to align
        gap_penalty: Gap penalty (default -12)

    Returns:
        tuple: (alignment1, alignment2, score)
    """
    # load BLOSUM62 matrix
    blosum_matrix, aa_to_idx = load_blosum62("blosum62.txt")

    n, m = len(seq1), len(seq2)
    score_matrix = np.zeros((n + 1, m + 1), dtype=np.float32)
    traceback = np.zeros((n + 1, m + 1), dtype=np.int8)

    # fill score matrix
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # calculate match score using BLOSUM62
            match_score = get_score(blosum_matrix, aa_to_idx, seq1[i - 1], seq2[j - 1])

            # calculate scores for all possible moves
            diag = score_matrix[i - 1, j - 1] + match_score
            up = score_matrix[i - 1, j] + gap_penalty
            left = score_matrix[i, j - 1] + gap_penalty

            # find best score
            score_matrix[i, j] = max(0, diag, up, left)

            # store traceback information
            if score_matrix[i, j] == 0:
                traceback[i, j] = 0  # Stop
            elif score_matrix[i, j] == diag:
                traceback[i, j] = 1  # Diagonal
            elif score_matrix[i, j] == up:
                traceback[i, j] = 2  # Up
            else:
                traceback[i, j] = 3  # Left

            # update maximum score
            if score_matrix[i, j] > max_score:
                max_score = score_matrix[i, j]
                max_pos = (i, j)

    # traceback to find alignment
    align1, align2 = [], []
    i, j = max_pos

    while i > 0 and j > 0 and score_matrix[i, j] > 0:
        if traceback[i, j] == 1:  # Diagonal
            align1.append(seq1[i - 1])
            align2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif traceback[i, j] == 2:  # Up
            align1.append(seq1[i - 1])
            align2.append('-')
            i -= 1
        elif traceback[i, j] == 3:  # Left
            align1.append('-')
            align2.append(seq2[j - 1])
            j -= 1
        else:  # Stop
            break

    # reverse the alignments
    align1 = ''.join(reversed(align1))
    align2 = ''.join(reversed(align2))

    return align1, align2, max_score