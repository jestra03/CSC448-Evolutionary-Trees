## Part 1. Sequences

In this first part you will combine what you learned about aligning sequences, working with substitution matrices, using similarity scores and clustering algorithms to ultimately build an evolutionary tree given orthologous proteins sampled from several strains of bacteria.

0. The protein sequences are in the data/ folder in the github page. Each sequences begins on a new line. 

1. Align two sequences using the Smith–Waterman algorithm:
- you can borrow an implementation from online, however,
- recall that you need to penalize each substitution and you need to use the BLOSUM62 matrix rather than the simple -1 since you are working with protein sequences
- BLOSUM62 substitution matrix is available here: https://anaconda.org/bioconda/blosum
- for gap penalty use -12
- note that the Smith–Waterman algorithm score includes all the penalties for mismatches and gaps as well as positive values for correctly matched, hence it reflects how well the two sequences align or in other words are similiar to each other

2. Use your implementation of 1 to align every pair of sequences in the file:
- convert the alignment score from (1) into a similarity measure for a given pair of sequences (think and explain how you did this given the Smith–Waterman algorithm score)
- given similarity scores for any pair of sequences, construct the similarity matrix (this matrix is the input to the clustering)

3. Use a clustering algorithm of your choice to reconstruct the phylogenic tree:
- plot the tree (the dendrogram function of matplotlib may come in handy)
- what do you notice about the tree, how would you interpret it
- compare your tree to that of a classmate (discuss your observations in the report)

4. Think about and suggest a way (metric) to systematically compare the trees reconstructed by the entire class


5. Write a brief pdf report including:
- plot of the phylogenic tree with your interpretation and discussion
- concise description of your proposed metric to assess similarity between reconstructed trees
- report of the two closets and farthest away proteins (use row number as id)

6. Turn in the pdf and your code via canvas by the *deadline: Tuesday 28th, 11:59pm*


## Part 2. Function

In the second part, you will investigate the biological functions of DNA sequences, that is you will try to figure out what are some of the DNA sequences responsible for. Specifically you will look into protein coding genes and their expression levels in a species of yeast. The goal is to determine which genes are active during the degradation of glucose and ethanol.

0. Biological background: 

Saccharomyces cerevisiae is a species of yeast that brews wine by converting the glucose found in fruit into ethanol.

- If the supply of glucose runs out, *S. cerevisiae* must do something to survive
- It will then invert its metabolism, with the ethanol (alcohol) that it just produced becoming its new food supply.
- This metabolic inversion, called the diauxic shift, can only occur in the presence of oxygen.
- Without oxygen, *S. cerevisiae* hibernates until either glucose or oxygen becomes available


In conclusion, if winemakers don’t seal their barrels, then the yeast in the barrel will metabolize the ethanol that it just produced, ruining the wine.

The diauxic shift is a complex process that affects the expression of many genes.


In 1997, Joseph DeRisi conducted the first massive gene expression experiment by sampling an *S. cerevisiae* culture every two hours for the six hours before and after the diauxic shift. Since there are approximately 6,400 genes in *S. cerevisiae*, and there were seven time points, this experiment resulted in a 6,400 × 7 gene expression matrix.

1. Data: in the github data/diauxic\_raw\_ratios.txt

- Note that it is gene by time point format while in class we mostly discussed gene by cell (sample).

2. Plot the row data as a heatmap

3. Cluster the genes and plot the clustered matrix
- describe what algorithm you chose, how many clusters you got (big parameter! to pay attention to), calculate the Silhouette score
- discuss whether you can readily identify a group of genes that changes expression patterns together

4. Remove genes not of interest

While you used all 6,400 genes in your initial analysis, typically, genes that have not changed their expression substantially or are known not play a role in the process being investigated are removed to increase the statistical power of the test. That is why most studies actually focus on a subset of *highly variable genes*. 

- devise and implement your own metric for assessing which genes have changed the most (both up and down direction!)
- describe it in your pdf explaining why it is reasonable and what pitfalls it may have
- report the 230 most variable genes according to your metric


5. The authors of the study did their own post-processing of the data and selected a set of 230 most variable genes. 
- Find their set in data/230genes_log_expression.txt
- what is the overlap with your list, i.e how many of the genes are the same, make sure to report together the Jaccard coefficient as well


6. Using the reduced list of 230 genes (derived by the authors) redo your clustering analysis
- are the groups of genes active at different stages identifiable now
 

7. Turn in the pdf and your code via canvas by the *deadline: Friday 31st, 11:59pm*
