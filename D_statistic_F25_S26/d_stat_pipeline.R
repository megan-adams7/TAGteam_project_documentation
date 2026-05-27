# TAGteam Binary Data Processing and D-Statistic Pipeline
# -------------------------------------------------------
# Purpose:
#   1. Load transcript-level TAGteam binary data.
#   2. Map transcript IDs to gene IDs.
#   3. Convert the data to wide format with species as rows and genes as columns.
#   4. Identify and compress species/gene cells with multiple transcript values.
#   5. Save both the original wide data and compressed wide data as CSV files.
#   6. Run D-statistic analyses to test phylogenetic signal.
#   7. Generate summary plots, heatmaps, tree PDFs, and candidate-gene rankings.
#
# Notes:
#   - File paths may need to be updated if this script is moved to another computer.

# Loading External Libraries

library(dplyr)
library(tidyr)
library(purrr)
library(ape)
library(caper)
library(stringr)
library(ggplot2)
library(phytools)
library(geiger)
library(pheatmap)

# Create output folders before saving plots/files.
# This prevents errors if the folders do not already exist.
dir.create("out.dir", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("out.dir", "trees"), showWarnings = FALSE, recursive = TRUE)

# Helper function used only for saving wide data frames as CSVs.
# pivot_wider() can create list-columns when a species/gene combination has
# more than one transcript value. Base write.csv() cannot save list-columns
# directly, so this converts each cell to a readable character value.
make_csv_safe <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.list(col)) {
      sapply(col, function(x) {
        if (is.null(x)) {
          return(NA_character_)
        }
        paste(unlist(x), collapse = ";")
      })
    } else {
      col
    }
  })
  df
}

# Load in binary data file
# This file contains transcript-level TAGteam binary.
binary_data <- read.csv("/local/workdir/maa367/d_statistic/tagteam-binary-df-rbbpass_full.csv")

# Load in mapping file of transcript ID to gene ID.
# This lets transcript-level binary values be grouped by gene.
txid_gnid <-read.table("/local/workdir/maa367/d_statistic/txid_gnid_dictionary.fasta", header = FALSE, sep = "\t",
                                   col.names = c("txid", "gnid"))

# Map each transcript ID in the binary dataset to its corresponding gene ID.
binary_data_mapped <- binary_data %>%
  left_join(txid_gnid, by = "txid") 

# Change binary data to wide format.
# Rows = species; columns = genes; values = TAGteam binary.
binary_data_wide <- binary_data_mapped %>%
  dplyr::select(gnid, species, tagteam.binary) %>%
  pivot_wider(names_from = gnid, values_from = tagteam.binary)

# Save the original wide-format binary data before any compression.
# This creates a record of the data immediately after transcript IDs are mapped to genes.
write.csv(
  make_csv_safe(binary_data_wide),
  "binary_data_wide_original.csv",
  row.names = FALSE
)

# Identify transcript names
gene_cols <- setdiff(names(binary_data_wide), "species")

# Summary table to track where species/gene cells contain multiple transcript values.
summary_n <- data.frame(
  species = gene_cols,
  n_multivectors = rep(NA, length(gene_cols)),
  n_with_diff_values = rep(NA, length(gene_cols)),
  stringsAsFactors = FALSE
)


# Compression rule for multiple transcript values in one species/gene cell.
# If any transcript has a TAGteam motif coded as 1, the whole gene/species cell becomes 1.
# If all available transcript values are 0, the compressed value becomes 0.
any1 <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0) return(NA_integer_)
  as.integer(any(v == 1))
}

# Keep copies of the wide data.
# binary_data_wide remains the original wide-format object.
# binary_data_wide1 will be compressed so multi-transcript cells become one value.
# binary_data_wide2 is kept as an extra backup copy if needed later.
binary_data_wide1 <- binary_data_wide
binary_data_wide2 <- binary_data_wide

# Loop through each gene column and compress any species/gene cells
# that contain multiple transcript-level binary values.
# Loop through each gene and calculate Fritz and Purvis' D-statistic.
for (i in seq_along(gene_cols)) {
  g <- gene_cols[i]

  cells <- lapply(binary_data_wide1[[g]], function(x) if (is.null(x)) NA else x)
  lens <- lengths(cells)
  is_multi <- lens >= 2
  n_multi <- sum(is_multi)
  
  n_diff <- if (n_multi == 0) {
    0
  } else {
    sum(vapply(cells[is_multi], function(v) length(unique(v)) > 1L, logical(1L)))
  }
  
  summary_n$n_multivectors[summary_n$species == g] <- n_multi
  summary_n$n_with_diff_values[summary_n$species == g] <- n_diff
  
  if (n_multi) {
    new_vals <- vapply(cells[is_multi], any1, integer(1))
    binary_data_wide1[[g]][is_multi] <- as.list(new_vals)
  }
}

# Save the compressed wide-format binary data.
# This is the version used for downstream D-statistic analysis.
write.csv(
  make_csv_safe(binary_data_wide1),
  "binary_data_wide_compressed.csv",
  row.names = FALSE
)

# Save a summary of how often multiple transcript values appeared
# and how often those values disagreed within a species/gene cell.
write.csv(summary_n, "binary_data_compression_summary.csv", row.names = FALSE)


# Helper functions for summarizing list-cells with multiple transcript values
vec_length <- function(x) {
  if (is.null(x)) return(NA_integer_)
  length(x)
}

vec_mean01 <- function(x) {
  if (is.null(x)) return(NA_real_)
  x <- unlist(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

# Objects to collect vector-length and vector-mean information across all genes
all_vector_lengths <- c()
all_vector_means   <- c()

# Per-gene summary table for multi-transcript cells
vector_summary <- data.frame(
  gene = gene_cols,
  n_vectors = NA_integer_,
  mean_vector_length = NA_real_,
  median_vector_length = NA_real_,
  max_vector_length = NA_integer_,
  mean_of_vector_means = NA_real_,
  stringsAsFactors = FALSE
)

vector_means <- binary_data_wide
vector_means[,] <- NA
vector_means$species <- binary_data_wide$species

for (i in seq_along(gene_cols)) {
  g <- gene_cols[i]
  
  cells <- binary_data_wide[[g]]
  
  # identify cells that are actually vectors of length >= 2
  lens <- sapply(cells, function(x) {
    if (is.null(x)) return(0L)
    length(x)
  })
  # Uncomment these lines for debugging if you want to inspect each gene:
  # print(lens)

  is_multi <- lens >= 2
  # print(is_multi)
  
  # lengths of all multi-value vectors in this gene column
  these_lengths <- lens[is_multi]
  
  # means of all multi-value vectors in this gene column
  these_means <- sapply(cells[is_multi], vec_mean01)
  
  # only assign these means to rows that had multi-value vectors
  vector_means[[g]][is_multi] <- these_means
  
  # collect across all genes
  all_vector_lengths <- c(all_vector_lengths, these_lengths)
  all_vector_means   <- c(all_vector_means, these_means)
  
  # fill per-gene summary
  vector_summary$n_vectors[i] <- sum(is_multi)
  
  if (sum(is_multi) > 0) {
    vector_summary$mean_vector_length[i] <- mean(these_lengths)
    vector_summary$median_vector_length[i] <- median(these_lengths)
    vector_summary$max_vector_length[i] <- max(these_lengths)
    vector_summary$mean_of_vector_means[i] <- mean(these_means, na.rm = TRUE)
  }
}

# overall summary of vector lengths
summary(all_vector_lengths)

# overall summary of vector means
summary(unlist(all_vector_means))

# Save multi-transcript summaries for documentation.
write.csv(vector_summary, "multi_transcript_vector_summary_by_gene.csv", row.names = FALSE)
write.csv(make_csv_safe(vector_means), "multi_transcript_vector_means_by_species_gene.csv", row.names = FALSE)


# Histogram ---------------------------------------------------------------
# Create a histogram of the mean values when multiple binary values

library(ggplot2)

# create dataframe
hist_data <- data.frame(
  mean_binary_value = unlist(all_vector_means)
)

# make plot
p <- ggplot(hist_data, aes(x = mean_binary_value)) +
  geom_histogram(
    binwidth = 0.05,
    color = "black",
    fill = "steelblue"
  ) +
  labs(
    title = "Mean Binary Values Across Species",
    x = "Mean Binary Value",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

# display plot
print(p)

# save as PDF
ggsave(
  filename = "mean_binary_value_histogram.pdf",
  plot = p,
  width = 8,
  height = 6
)


# D-statistics ------------------------------------------------------------


# Species names present in the compressed wide-format binary data
species_names <- binary_data_wide1$species
# Read the phylogenetic tree file
tree <- read.tree("/local/workdir/maa367/d_statistic/hog_tree.nw")
# Modify tree tip labels to match the species naming format in the wide data file.
tree$tip.label <- sub("^([A-Z])[A-Z]+_(.*)$", "\\1_\\L\\2", tree$tip.label, perl = TRUE)

# Drop species that are not present in both the tree and the binary dataset.
drop <- symdiff(tree$tip.label, species_names)
tree <- drop.tip(tree, drop)


# Initialize results table for D-statistic output.
results <- data.frame(
  gene       = gene_cols,
  D          = NA_real_,
  P_random   = NA_real_,
  P_brownian = NA_real_,
  n          = NA_integer_,
  n1         = NA_integer_,
  n0         = NA_integer_,
  status     = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(gene_cols)) {
  g <- gene_cols[i]
  missing_species <- binary_data_wide1$species[ binary_data_wide1[[g]] == "NULL" ]
  
  # Identify species with missing data for this gene and prune them from the tree.
  tips_to_drop <- intersect(missing_species, tree$tip.label)
  tree_sub <- if (length(tips_to_drop)) drop.tip(tree, tips_to_drop) else tree
  #pruned_trees[[g]] <- tree_sub
  
  # Build a two-column data frame for this gene: species and binary trait value.
  dat_g <- binary_data_wide1[, c("species", g)]
  names(dat_g) <- c("species", "binary")
  dat_g$binary[dat_g$binary == "NULL"] <- NA
  dat_g <- subset(dat_g, species %in% tree_sub$tip.label & !is.na(binary))
  dat_g$binary <- as.integer(dat_g$binary)
  dat_g <- as.data.frame(dat_g, stringsAsFactors = FALSE)
  
  # Safety checks: skip genes with too few species or no variation in the binary trait.
  n_sp <- nrow(dat_g)
  n1   <- sum(dat_g$binary == 1)
  n0   <- sum(dat_g$binary == 0)
  
  # Skip genes that cannot be used for D-statistic testing.
  if (n_sp < 4 || n1 == 0 || n0 == 0) {
    results[i, c("D","P_random","P_brownian","n","n1","n0","status")] <-
      list(NA_real_, NA_real_, NA_real_, n_sp, n1, n0, "skipped")
    next
  }
  
  # Combine the gene data and pruned tree into a comparative.data object.
  comp <- comparative.data(
    phy = tree_sub,
    data = dat_g,
    names.col = "species",
    vcv = TRUE,
    warn.dropped = TRUE
  )

  
  # Calculate the D-statistic and associated p-values.
  dres <- phylo.d(comp, binvar = binary)
  
  results[i, c("D","P_random","P_brownian","n","n1","n0","status")] <-
    list(as.numeric(dres$DEstimate),
         as.numeric(dres$Pval1),
         as.numeric(dres$Pval0),
         n_sp, n1, n0, "ok")
  
}



# Multiple testing corrections for the D-statistic p-values.
results$p_adjustp <- p.adjust(results$P_random, "fdr")
results$p_adjustb <- p.adjust(results$P_brownian, "fdr")
results$n_prop <- results$n/300

# Add proportion of species with trait 1
results$n1_prop <- results$n1/results$n

# Save results as csv file
write.csv(results, "D_stat_results.csv", row.names = FALSE)

# Heatmap -----------------------------------------------------------------
# Visualize where strong phylogenetic signal tends to occur across species-count
# and motif-presence proportion bins.

# Filter only genes that passed
results_new <- results %>%
  filter(status == "ok")

# Filter criteria
negative_d <- results_new %>%
  filter(n > 30, D < 0, p_adjustp < 0.05, n1_prop > 0.1, n1_prop < 0.9)
negative_d$gene

# output filtered genes to csv
write.table(negative_d$gene,
            file = "negative_d_gene.csv",
            row.names = FALSE,
            col.names = FALSE,
            sep = ",",
            quote = FALSE)

# eps is kept here in case you want to use a smaller p-value floor later.
eps <- 1e-300 

heat_binned <- results_new %>%
  mutate(
    prop_1  = n1 / n,
    p_safe  = ifelse(is.na(p_adjustp) | p_adjustp <= 0, 0.000001, p_adjustp),
    neglogp = -log10(p_safe),
    
    species_bin = cut(
      n,
      breaks = seq(0, max(n, na.rm = TRUE) + 5, by = 2),
      include.lowest = TRUE,
      right = FALSE
    ),
    
    prop_bin = cut(
      prop_1,
      breaks = seq(0, 1, by = 0.02),
      include.lowest = TRUE,
      right = TRUE
    )
  ) %>%
  group_by(species_bin, prop_bin) %>%
  summarize(
    mean_neglogp = mean(neglogp, na.rm = TRUE),
    .groups = "drop"
  )


heat_wide <- heat_binned %>%
  pivot_wider(
    names_from  = prop_bin,
    values_from = mean_neglogp
  )

# Ordering function used to sort binned interval labels numerically.
extract_lower <- function(x) {
  x_num <- sub("^\\(|^\\[", "", x)  # remove leading ( or [
  x_num <- sub(",.*", "", x_num)    # keep everything before first comma
  as.numeric(x_num)
}


heat_mat <- as.data.frame(heat_wide)
rownames(heat_mat) <- heat_mat$species_bin
heat_mat$species_bin <- NULL
heat_mat <- as.matrix(heat_mat)


keep_rows <- rowSums(!is.na(heat_mat)) > 0
keep_cols <- colSums(!is.na(heat_mat)) > 0
heat_mat <- heat_mat[keep_rows, keep_cols, drop = FALSE]

row_ord <- order(extract_lower(rownames(heat_mat)))
col_ord <- order(extract_lower(colnames(heat_mat)))

heat_mat <- heat_mat[row_ord, col_ord, drop = FALSE]




bk <- seq(-1, 1, length.out = 101)
col_fun <- colorRampPalette(c("navy", "white", "firebrick3"))(length(bk) - 1)

pdf("my_heatmap.pdf", width = 10, height = 8)  

pheatmap(
  heat_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  #color = col_fun,
  fontsize_row = 3,
 # breaks = bk,
  main = "Mean -log(p value)",
  border_color = NA
)

dev.off()                                  



# Few Gene pdf ------------------------------------------------------------
# Plot selected candidate genes as individual tree PDFs.


passed_genes <- c("FBgn0003448", "FBgn0004053", "FBgn0004170", "FBgn0001320", "FBgn0000606", "FBgn0004143", "FBgn0010109", 
                  "FBgn0003510", "FBgn0001180", "FBgn0264270", "FBgn0000490", "FBgn0001168", "FBgn0002985")
names_passed <- c("sna", "zen", "scute", "kni", "eve", "nullo", "dpn", "sry-a", "hb", "sxl", "dpp", "hry", "odd") 
i <- 0
for (g in passed_genes[1:13]) {
  i <- i + 1
  # Raw values for this gene (list of NULL / 0 / 1)
  vals_raw <- binary_data_wide1[[g]]
  
  # Identify species with NULL (to drop)
  is_null <- sapply(vals_raw, is.null)
  missing_species <- binary_data_wide1$species[is_null]
  tips_to_drop    <- intersect(missing_species, tree$tip.label)
  
  # Drop those tips from the tree
  tree2 <- if (length(tips_to_drop) > 0) drop.tip(tree, tips_to_drop) else tree
  
  # Turn vals_raw into a numeric vector (0/1/NA)
  vals_vec <- sapply(vals_raw, function(x) if (is.null(x)) NA_real_ else as.numeric(x))
  
  # Reorder values to match the tips of tree2
  vals <- vals_vec[match(tree2$tip.label, binary_data_wide1$species)]
  
  # Build a color vector: red for 0, blue for 1, grey for NA (if any)
  cols <- ifelse(is.na(vals), "grey70",
                 ifelse(vals == 1, "blue", "red"))
  
  # Get D for this gene from results
  D_g <- results$D[results$gene == g]   # change 'gene' to the right col name if needed
  
  n_tips <- length(tree2$tip.label)
  tip_cex <- min(0.3, 0.5 / log10(n_tips)) 
  
  # Open a high-resolution pdf device
  pdf(
    file   = file.path("out.dir/trees", paste0("tree_1", names_passed[i], ".pdf")),
    width  = 7,   # inches
    height = 7
  )
  
  # Plot the tree: species-only labels, colored by trait, D in title
  plot(
    tree2,
    tip.color = cols,
    cex = tip_cex,
    main = paste0("Tree for ", names_passed[i], "   (D = ", round(D_g, 7), ")")
  )
  
  legend(
    "topleft",
    legend = c("Trait = 0", "Trait = 1"),
    col = c("red", "blue"),
    pch = 19,
    pt.cex = 1.2,
    bty = "n"     # no legend box
  )
  
  dev.off()
}


# All gene pdf -----------------------------------------------------------------
# Plot all filtered candidate genes into one multi-page PDF, ordered by D-statistic.

# Read in file containing gene IDs and corresponding gene symbols
gene_names <- read.csv("genes - Sheet1 (3).csv")

# Merge gene symbols into the negative_d dataframe using the gene column
negative_d <- negative_d %>%
  left_join(gene_names, by = "gene")

# Create ranked dataframes based on different metrics:
# 1. Highest proportion of 1s (n1_prop)
ordered_genes_n1prop <- negative_d %>%
  arrange(desc(n1_prop)) %>%
  dplyr::select(gene, symbol.x)

# 2. Lowest adjusted p-value
ordered_genes_pval <- negative_d %>%
  arrange(p_adjustp) %>%
  dplyr::select(gene, symbol.x)

# 3. Most negative D-statistic
ordered_genes_D <- negative_d %>%
  arrange(D) %>%
  dplyr::select(gene, symbol.x)

# Load ape package for phylogenetic tree visualization
library(ape)

# Open a PDF device to save all generated trees
pdf(
  file = file.path("out.dir/trees", "all_gene_trees_D.pdf"),
  width = 7,
  height = 7
)

# Extract ordered list of gene IDs ranked by D-statistic
ordered_gene_ids <- ordered_genes_D$gene

# Loop through every gene and generate a phylogenetic tree plot
for (i in seq_along(ordered_gene_ids)) {
  
  # Current gene ID
  g <- ordered_gene_ids[i]
  
  # Corresponding gene symbol for labeling plots
  symbol_g <- ordered_genes_D$symbol.x[i]
  
  # Extract binary trait values for the current gene across species
  vals_raw <- binary_data_wide1[[g]]
  
  # Identify species with missing values (NULL entries)
  is_null <- sapply(vals_raw, is.null)
  
  # Get names of species missing data
  missing_species <- binary_data_wide1$species[is_null]
  
  # Find which missing species are present in the phylogenetic tree
  tips_to_drop <- intersect(missing_species, tree$tip.label)
  
  # Remove species with missing data from the tree
  # If no species are missing, use the original tree
  tree2 <- if (length(tips_to_drop) > 0) drop.tip(tree, tips_to_drop) else tree
  
  # Convert list-based binary values into a numeric vector
  # NULL values become NA
  vals_vec <- sapply(vals_raw, function(x) {
    if (is.null(x)) NA_real_ else as.numeric(x)
  })
  
  # Reorder binary values so they match the order of species in the tree
  vals <- vals_vec[match(tree2$tip.label, binary_data_wide1$species)]
  
  # Assign colors to trait values:
  # blue = 1
  # red = 0
  # grey = missing/NA
  cols <- ifelse(is.na(vals), "grey70",
                 ifelse(vals == 1, "blue", "red"))
  
  # Extract D-statistic value for this gene
  D_g <- results$D[results$gene == g]
  
  # Extract adjusted p-value for this gene
  P <- results$p_adjustp[results$gene == g]
  
  # Dynamically scale tip label size based on tree size
  # Larger trees get smaller labels
  n_tips <- length(tree2$tip.label)
  tip_cex <- min(0.3, 0.5 / log10(n_tips))
  
  # Plot phylogenetic tree with colored trait labels
  plot(
    tree2,
    tip.color = cols,
    cex = tip_cex,
    
    # Plot title includes:
    # gene symbol
    # D-statistic
    # adjusted p-value
    main = paste0(
      "Tree for ",
      symbol_g,
      "   (D = ",
      round(D_g, 7),
      ") (P = ",
      round(P, 6),
      ")"
    )
  )
  
  # Add legend explaining color coding
  legend(
    "topleft",
    legend = c("Trait = 0", "Trait = 1"),
    col = c("red", "blue"),
    pch = 19,
    pt.cex = 1.2,
    bty = "n"
  )
}

# Close the PDF device and save the file
dev.off()



# Prioritization --------------------------------------------------------------
# Compare rankings across different criteria and create a combined ranking score.

# Define how many top-ranked genes to keep from each ranking method
top_n <- 600

# Select the top 600 genes ranked by highest proportion of 1s (n1_prop)
# Pull only the gene symbols into a vector
top_n1prop <- ordered_genes_n1prop %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::pull(symbol.x)

# Select the top 600 genes ranked by most negative D-statistic
# Pull only the gene symbols into a vector
top_D <- ordered_genes_D %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::pull(symbol.x)

# Find genes shared between both ranking methods
# intersect() identifies overlapping genes
# Reduce() allows intersect to work across a list of vectors
top_in_all_3 <- Reduce(intersect, list(top_n1prop, top_D))

# Convert overlapping gene symbols into a dataframe
top_in_all_3_df <- data.frame(symbol.x = top_in_all_3)

# View dataframe of genes shared between ranking methods
top_in_all_3_df


# Create a combined ranking system using both:
# 1. D-statistic ranking
# 2. Proportion of 1s ranking
ranked_genes <- negative_d %>%
  
  # Add ranking columns
  dplyr::mutate(
    
    # Rank genes by D-statistic
    # Lower (more negative) D values receive better ranks
    # ties.method = "average" assigns tied values the average rank
    rank_D = rank(D, ties.method = "average"),
    
    # Rank genes by proportion of 1s
    # Negative sign is used so larger n1_prop values rank higher
    rank_n1 = rank(-n1_prop, ties.method = "average"),
    
    # Compute an overall combined ranking score
    # Lower combined_rank values indicate genes that rank well in both metrics
    combined_rank = (rank_D + rank_n1)/2
  ) %>%
  
  # Sort genes from best combined rank to worst
  dplyr::arrange(combined_rank)
