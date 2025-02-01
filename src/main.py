import os
from similarity import compute_similarity_matrix
from clustering import generate_phylogenetic_tree, find_closest_farthest_pairs
from plotting import plot_tree, plot_heatmap
from gene_expression import load_gene_data, cluster_genes, top_variable_genes, compare_gene_sets

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(script_dir, "..", "data")
    results_dir = os.path.join(script_dir, "..", "results")
    os.makedirs(results_dir, exist_ok=True)

    print("Choose an analysis:")
    print("1. Phylogenetic Tree")
    print("2. Gene Expression")
    choice = input("Enter 1 or 2: ")

    if choice == "1":
        with open(os.path.join(data_path, "seq.txt")) as f:
            sequences = [line.strip() for line in f]

        similarity_matrix = compute_similarity_matrix(sequences)
        linkage_matrix = generate_phylogenetic_tree(similarity_matrix)
        plot_tree(linkage_matrix, labels=[f"Seq {i}" for i in range(len(sequences))])

        closest, farthest = find_closest_farthest_pairs(similarity_matrix)
        print(f"Closest sequences: {closest[0] + 1} and {closest[1] + 1} with score {closest[2]:.3f}")
        print(f"Farthest sequences: {farthest[0] + 1} and {farthest[1] + 1} with score {farthest[2]:.3f}")

    elif choice == "2":
        gene_data = load_gene_data(os.path.join(data_path, "diauxic_raw_ratios.txt"))
        plot_heatmap(gene_data, 'Raw Gene Expression Data', os.path.join(results_dir, "raw_heatmap.png"))

        clusters, silhouette_score = cluster_genes(gene_data, n_clusters=8, method="ward", metric="correlation")

        top_genes = top_variable_genes(gene_data)
        authors_genes = load_gene_data(os.path.join(data_path, "230genes_log_expression.txt"))

        overlap, jaccard = compare_gene_sets(top_genes, authors_genes)
        reclustered_data = gene_data.loc[authors_genes.index]
        reclustered_clusters, reclustered_silhouette = cluster_genes(reclustered_data, n_clusters=8, method="ward", metric="correlation")

        print(f"Initial clustering - Silhouette score: {silhouette_score:.3f}")
        print(f"Top {len(top_genes)} variable genes identified")
        print(f"Overlap with authors' 230 genes: {overlap} / 230, Jaccard coefficient: {jaccard:.3f}")
        print(f"Reclustered with 230 genes - Silhouette score: {reclustered_silhouette:.3f}")
    else:
        print("Invalid choice. Exiting.")

if __name__ == "__main__":
    main()
