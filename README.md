# Phylogenetic Analysis Program Overview

The program performs a phylogenetic analysis on a set of protein sequences. The overall workflow is as follows:

## File Descriptions

1. **blosum.py**:
   - This file contains functions to load the BLOSUM62 substitution matrix and calculate the score for a given amino acid pair.
   - The BLOSUM62 matrix is used during the sequence alignment step.

2. **clustering.py**:
   - This file defines functions to generate a phylogenetic tree using the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) clustering algorithm.
   - The `generate_phylogenetic_tree()` function takes a similarity matrix as input and returns the linkage matrix for plotting the dendrogram.
   - The `find_closest_farthest_pairs()` function identifies the closest and farthest pairs of sequences based on the similarity matrix.

3. **main.py**:
   - This is the main entry point of the program.
   - It loads the protein sequences from a file, computes the similarity matrix using the Smith-Waterman algorithm and the BLOSUM62 matrix, generates the phylogenetic tree, and plots the resulting dendrogram.
   - The closest and farthest sequence pairs are also reported.

4. **plotting.py**:
   - This file contains the `plot_tree()` function, which takes the linkage matrix and sequence labels as input and generates a dendrogram plot of the phylogenetic tree.
   - The plot is saved to a file in the `results/` directory.

5. **similarity.py**:
   - This file defines the `compute_similarity_matrix()` function, which calculates the pairwise similarity scores for all the input sequences using the Smith-Waterman algorithm.
   - The function uses multiprocessing to speed up the computation and saves the resulting similarity matrix to a file for later use.

## Program Flow

1. Load the protein sequences from a file.
2. Compute the pairwise similarity matrix using the Smith-Waterman algorithm and the BLOSUM62 matrix.
3. Generate the phylogenetic tree using the UPGMA clustering algorithm.
4. Plot the dendrogram of the phylogenetic tree and save it to a file.
5. Identify the closest and farthest sequence pairs based on the similarity matrix.

The program leverages the different modules to handle specific tasks, such as sequence alignment, clustering, and visualization, to perform the complete phylogenetic analysis.
