# TAGteam-projectF24

This folder contains the documentation for the TAGteam project from Fall 2025 and Spring 2026. 

## Folder Structure

## Folder Structure

/d_stat_pipeline.R  
Contains the main R pipeline used to process binary TAGteam motif data, calculate Fritz & Purvis D-statistics, rank candidate genes, and generate summary visualizations.

/tree_visualizations  
Contains ranked phylogenetic tree visualizations for candidate genes. Trees display binary TAGteam motif presence/absence patterns across Drosophila species and were used to visually assess phylogenetic clustering and conservation.
all_gene_trees_D.pdf is ranked by most negative D statistic
all_gene_trees_n1prop.pdf is ranked by the proportion of 1s in the phylogenetic tree
all_gene_trees_pval.pdf is ranked by the smalles p-value

/D_stat_results.csv  
Contains the final D-statistic output table, including D values, p-values, adjusted p-values, species counts, and motif-presence proportions used for downstream ranking analyses.

/hog_tree.nw  
Contains the phylogenetic tree in Newick format used for all phylogenetic signal analyses and tree visualizations.

/tagteam-binary-df-rbbpass_full.csv  
Contains the original binary TAGteam motif dataset generated after reciprocal best BLAST hit (RBBH) filtering. Rows represent species and columns represent genes/transcripts, with binary values indicating motif presence or absence.

/txid_gnid_dictionary.fasta  
Contains the transcript ID to gene ID mapping dictionary used to aggregate transcript-level binary values into gene-level classifications.

/histogram.png  
Contains histogram visualizations generated during exploratory data analysis to examine distributions of binary summary statistics across genes and species.

/heatmap.pdf  
Contains heatmap visualizations summarizing phylogenetic signal patterns and D-statistic significance across grouped bins of genes/species.

