# Define variable genes from RNA-seq normalized read counts.
Microcolony-seq methodology finds co-existing subpopulations in a bacterial culture based on RNA-seq data. RNA-seq is applied separately on each selected microcolony formed on solid media.
To identify subpopulations, the most variable genes among the selected microcolonies are determined.
The input is a table with normalized read counts of all samples (each column represents one microcolony). Normalization using DESeq2 (Love et al.).
The output is a table with True/False values for each gene to indicate if it is a variable gene.
Variable genes are then used in a PCA plot to determine subpopulations.
