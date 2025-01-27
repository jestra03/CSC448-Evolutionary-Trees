import numpy as np
import os
from multiprocessing import Pool
from tqdm import tqdm
from smith_waterman import smith_waterman

# move worker function outside
def worker_function(pair_data):
    """Worker function for parallel processing"""
    i, j, seq1, seq2 = pair_data
    _, _, score = smith_waterman(seq1, seq2)
    return i, j, max(score, 0)

def compute_similarity_matrix(sequences, save_path="results/similarity_matrix.npy"):
    """
    Computes or loads the similarity matrix for all sequence pairs.
    Uses parallel processing with memory constraints in mind.
    """
    # check if similarity matrix already exists
    if os.path.exists(save_path):
        print(f"Loading pre-computed similarity matrix from {save_path}")
        return np.load(save_path)

    n = len(sequences)
    matrix = np.zeros((n, n))

    # first compute self-alignments
    print("Computing self-alignments...")
    self_scores = []
    for i in tqdm(range(n), desc="Self alignments"):
        _, _, score = smith_waterman(sequences[i], sequences[i])
        self_scores.append(max(score, 0))
        matrix[i, i] = 1.0

    # prepare pairs for parallel processing (excluding self-alignments)
    pairs = [(i, j, sequences[i], sequences[j])
            for i in range(n)
            for j in range(i+1, n)]

    # use fewer cores to avoid memory issues
    num_cores = min(10, os.cpu_count() or 1)  # Use at most 6 cores
    print(f"Using {num_cores} CPU cores for parallel processing")

    # process pairs in parallel with smaller chunks
    chunk_size = 10  # Process smaller chunks at a time
    with Pool(num_cores) as pool:
        results = list(tqdm(
            pool.imap(worker_function, pairs, chunksize=chunk_size),
            total=len(pairs),
            desc="Computing pairwise alignments"
        ))

    # fill matrix with normalized scores
    print("Finalizing similarity matrix...")
    for i, j, score in results:
        denominator = np.sqrt(self_scores[i] * self_scores[j])
        if denominator > 0:
            normalized_score = score / denominator
        else:
            normalized_score = 0.0
        matrix[i, j] = matrix[j, i] = normalized_score

    # save intermediate results periodically
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    np.save(save_path, matrix)
    print(f"Saved similarity matrix to {save_path}")

    return matrix