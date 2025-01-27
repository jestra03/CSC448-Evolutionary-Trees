import os
from similarity import compute_similarity_matrix
from clustering import generate_phylogenetic_tree
from clustering import find_closest_farthest_pairs
from plotting import plot_tree

def main():
    # get the directory containing the script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # construct paths
    data_path = os.path.join(script_dir, "..", "data", "seq.txt")
    results_dir = os.path.join(script_dir, "..", "results")
    os.makedirs(results_dir, exist_ok=True)

    # load sequences
    with open(data_path) as f:
        sequences = [line.strip() for line in f]

    print(f"Loaded {len(sequences)} sequences")

    # compute similarity matrix
    similarity_matrix = compute_similarity_matrix(sequences)

    closest, farthest = find_closest_farthest_pairs(similarity_matrix)
    print(f"Closest sequences: {closest[0] + 1} and {closest[1] + 1} with score {closest[2]:.3f}")
    print(f"Farthest sequences: {farthest[0] + 1} and {farthest[1] + 1} with score {farthest[2]:.3f}")

    # generate phylogenetic tree
    linkage_matrix = generate_phylogenetic_tree(similarity_matrix)

    # plot the tree
    plot_tree(linkage_matrix, labels=[f"Seq {i}" for i in range(len(sequences))])


if __name__ == "__main__":
    main()