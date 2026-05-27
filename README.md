# TAGteam_project_documentation
This repository contains code, data, and documentation for a research project analyzing the conservation and distribution of TAGteam motifs in the upstream regions of genes across 300+ Drosophila species.

This project investigates the abundance, distribution, and evolutionary conservation of TAGteam motifs across more than 300 species in the genus Drosophila. TAGteam motifs are short cis-regulatory DNA sequences involved in early embryonic gene regulation during the maternal-to-zygotic transition (MZT), where developmental control shifts from maternally deposited transcripts to activation of the zygotic genome. Their presence upstream of genes is associated with transcriptional activation during early development.

The project combines computational genomics, phylogenetics, and statistical analysis to examine how TAGteam motifs are conserved or lost across evolutionary lineages. Using custom pipelines written in Python, Bash, and R, genomic annotation and sequence files were processed to extract upstream regulatory regions, identify canonical and degenerate TAGteam motifs, and quantify motif abundance and positional patterns across species.

The first phase of the project focused on visualizing TAGteam motif abundance upstream of genes across ~300 Drosophila species to evaluate large-scale conservation patterns. This work was inspired by findings from [Bosch et al. (2006)](https://pubmed.ncbi.nlm.nih.gov/16624855/), which identified the developmental significance of TAGteam motifs in early embryogenesis.

The second phase explored evolutionary signal within binary motif presence/absence data using phylogenetic comparative methods. For each gene, species were encoded based on whether at least one TAGteam motif was present upstream. Phylogenetic clustering and conservation were then assessed using the [Frits and Purvis D statistic (2010)](https://pubmed.ncbi.nlm.nih.gov/20184650/), allowing comparison of observed motif distributions to expectations under random and Brownian evolutionary models.

The project additionally investigated:

- Genes consistently enriched for TAGteam motifs across species
- Lineage-specific motif gains and losses
- Relationships between motif abundance and phylogenetic structure
- Species- and gene-specific regulatory patterns
- Computational strategies for scaling motif analysis across hundreds of genomes

Results were visualized in R using phylogenetic trees, heatmaps, and comparative statistical summaries to highlight evolutionary trends in developmental regulatory networks across the Drosophila phylogeny.

This research was conducted in the Cornell University Barbash Lab from Summer 2024 through Spring 2026 as part of an undergraduate research project in computational and evolutionary genomics.

## Tools and Technologies

**Languages:** Python, Bash, R
**Libraries:** bedops, bedtools, seqkit, ggplot2, gggenomes, phytools, dplyr, ape, tidyr, purrr, caper, stringr, pheatmap
**Data:** GFF/GTF gene annotations, FASTA genome files for >300 Drosophila species
**Environment:** UNIX-based command line, RStudio, VS Code

Features
- motif scanning
- multi-species analysis
- unbiased motif search
- data-visualizations

## Data Sources

This analysis relies on genome annotations obtained from the Obbard Lab Drosophila annotations. These annotation files are not included in this repository because they consist of large genomic datasets and are not publicly redistributed here.

## Repository Structure

```text
TAGteam_project_MegAdams/
│
├── README.md
│   Project overview, workflow, and documentation
│
├── TAGteam_project_F24/
│   Original Fall 2024 TAGteam motif analysis pipeline
│
│   ├── main.sh # Main script
│   ├── vis.rmd # Visualization script
│   ├── gene_name # Working directory
│   ├── visualizations # .png visualizations
│   └── README.md
│
│
├── D_statistic_F25_S26/
│   Fall 2025 and Spring 2026 phylogenetic signal analyses using binary
│   TAGteam motif presence/absence data
│
│   ├── tree_visualizations # Ranked phylogenetic tree visualizations for candidate genes
│   ├── d_stat_pipeline.R/ # Main R script for D-statistic analysis pipeline
│   ├── D_stat_results.csv/ # Output table containing D-statistic results and significance values
│   ├── hog_tree.nw/ # Phylogenetic tree in Newick format
│   ├── tageteam-binary-df-rbbpass_full.csv/ # Original binary motif presence/absence dataset
│   ├── txid_gnid_dictionary.fasta/ # FASTA dictionary mapping transcript IDs to gene IDs
│   ├── histogram.png/ # Histogram visualization
│   ├── heatmap.pdf/ # Heatmap visualization
│   └── README.md
│
└── results/

```