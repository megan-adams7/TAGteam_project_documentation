# Loading External Libraries

library(dplyr)
library(tidyr)
library(purrr)
library(ape)
library(caper)
library(stringr)
library(ggplot2)
library(ape)
library(phytools)
library(geiger)

# Load in binary data file
binary_data <- read.csv("/local/workdir/maa367/d_statistic/tagteam-binary-df-rbbpass_full.csv")

#load in mapping file
txid_gnid <-read.table("/local/workdir/maa367/d_statistic/txid_gnid_dictionary.fasta", header = FALSE, sep = "\t",
                                   col.names = c("txid", "gnid"))

# map transcript to gene id
binary_data_mapped <- binary_data %>%
  left_join(txid_gnid, by = "txid") 

# change binary data to wide format
binary_data_wide <- binary_data_mapped %>%
  dplyr::select(gnid, species, tagteam.binary) %>%
  pivot_wider(names_from = gnid, values_from = tagteam.binary)

# identify transcript names
gene_cols <- setdiff(names(binary_data_wide), "species")

summary_n <- data.frame(
  species = gene_cols,
  n_multivectors = rep(NA, length(gene_cols)),
  n_with_diff_values = rep(NA, length(gene_cols)),
  stringsAsFactors = FALSE
)


any1 <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0) return(NA_integer_)
  as.integer(any(v == 1))
}

binary_data_wide1 <- binary_data_wide
binary_data_wide2 <- binary_data_wide

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


# helper functions
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

# objects to collect info across all genes
all_vector_lengths <- c()
all_vector_means   <- c()

# optional: per-gene summary table
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
  print(lens)
  
  is_multi <- lens >= 2
  print(is_multi)
  
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

h <- hist(unlist(all_vector_means))
hist(vector_summary$mean_of_vector_means)


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


#species names
species_names <- binary_data_wide1$species
# read the tree file
tree <- read.tree("/local/workdir/maa367/d_statistic/hog_tree.nw")
# modify tip labels to match wide data file
tree$tip.label <- sub("^([A-Z])[A-Z]+_(.*)$", "\\1_\\L\\2", tree$tip.label, perl = TRUE)

#drop species that are'nt in both trees
drop <- symdiff(tree$tip.label, species_names)
tree <- drop.tip(tree, drop)


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
  
  #print("tree modification")
  tips_to_drop <- intersect(missing_species, tree$tip.label)
  tree_sub <- if (length(tips_to_drop)) drop.tip(tree, tips_to_drop) else tree
  #pruned_trees[[g]] <- tree_sub
  
  #print("file mod")
  dat_g <- binary_data_wide1[, c("species", g)]
  names(dat_g) <- c("species", "binary")
  dat_g$binary[dat_g$binary == "NULL"] <- NA
  dat_g <- subset(dat_g, species %in% tree_sub$tip.label & !is.na(binary))
  dat_g$binary <- as.integer(dat_g$binary)
  dat_g <- as.data.frame(dat_g, stringsAsFactors = FALSE)
  
  #print("safety check")
  n_sp <- nrow(dat_g)
  n1   <- sum(dat_g$binary == 1)
  n0   <- sum(dat_g$binary == 0)
  
  #print("blank if")
  if (n_sp < 4 || n1 == 0 || n0 == 0) {
    results[i, c("D","P_random","P_brownian","n","n1","n0","status")] <-
      list(NA_real_, NA_real_, NA_real_, n_sp, n1, n0, "skipped")
    next
  }
  
  #print("comp")
  comp <- comparative.data(
    phy = tree_sub,
    data = dat_g,
    names.col = "species",
    vcv = TRUE,
    warn.dropped = TRUE
  )

  
  #print("d stat")
  dres <- phylo.d(comp, binvar = binary)
  
  results[i, c("D","P_random","P_brownian","n","n1","n0","status")] <-
    list(as.numeric(dres$DEstimate),
         as.numeric(dres$Pval1),
         as.numeric(dres$Pval0),
         n_sp, n1, n0, "ok")
  
}



# multiple testing corrections
results$p_adjustp <- p.adjust(results$P_random, "fdr")
results$p_adjustb <- p.adjust(results$P_brownian, "fdr")
results$n_prop <- results$n/300

# add proportion of species with trait 1
results$n1_prop <- results$n1/results$n

#sort by -log(p random)
heat_binned_sorted <- heat_binned[order(-heat_binned$mean_neglogp), ]

library(dplyr)

heat_binned_sorted <- heat_binned_sorted %>%
  mutate(
    prop_chr    = as.character(prop_bin),
    species_chr = as.character(species_bin),
    
    # extract all numeric substrings from each interval
    prop_nums    = lapply(prop_chr,    function(s) regmatches(s, gregexpr("[0-9.]+", s))[[1]]),
    species_nums = lapply(species_chr, function(s) regmatches(s, gregexpr("[0-9.]+", s))[[1]]),
    
    # lower and upper bounds = first and second numbers
    prop_lower    = as.numeric(vapply(prop_nums,    `[`, 1L, FUN.VALUE = character(1))),
    prop_upper    = as.numeric(vapply(prop_nums,    `[`, 2L, FUN.VALUE = character(1))),
    species_lower = as.numeric(vapply(species_nums, `[`, 1L, FUN.VALUE = character(1))),
    species_upper = as.numeric(vapply(species_nums, `[`, 2L, FUN.VALUE = character(1)))
  ) %>%
  dplyr::select(-prop_chr, -species_chr, -prop_nums, -species_nums)

#create emtpy list
vis_species <- c()  

for (i in seq_len(nrow(heat_binned_sorted))) {
  results_sorted1 <- results_new %>%
    filter(
      D < 0,
      p_adjustp < 0.05,
      n >= heat_binned_sorted$species_lower[i],
      n <  heat_binned_sorted$species_upper[i],
      n1_prop >  heat_binned_sorted$prop_lower[i], 
      n1_prop <= heat_binned_sorted$prop_upper[i]
    )
  
  vis_species <- c(vis_species, results_sorted1$gene)
}



# Few Gene pdf ------------------------------------------------------------


passed_genes <- c("FBgn0003448", "FBgn0004053", "FBgn0004170", "FBgn0001320", "FBgn0000606", "FBgn0004143", "FBgn0010109", 
                  "FBgn0003510", "FBgn0001180", "FBgn0264270", "FBgn0000490", "FBgn0001168", "FBgn0002985")
names_passed <- c("sna", "zen", "scute", "kni", "eve", "nullo", "dpn", "sry-a", "hb", "sxl", "dpp", "hry", "odd") 
i <- 0
for (g in passed_genes[1:13]) {
  i <- i + 1
  # 1. Raw values for this gene (list of NULL / 0 / 1)
  vals_raw <- binary_data_wide1[[g]]
  
  # 2. Identify species with NULL (to drop)
  is_null <- sapply(vals_raw, is.null)
  missing_species <- binary_data_wide1$species[is_null]
  tips_to_drop    <- intersect(missing_species, tree$tip.label)
  
  # 3. Drop those tips from the tree
  tree2 <- if (length(tips_to_drop) > 0) drop.tip(tree, tips_to_drop) else tree
  
  # 4. Turn vals_raw into a numeric vector (0/1/NA)
  vals_vec <- sapply(vals_raw, function(x) if (is.null(x)) NA_real_ else as.numeric(x))
  
  # 5. Reorder values to match the tips of tree2
  vals <- vals_vec[match(tree2$tip.label, binary_data_wide1$species)]
  
  # 6. Build a color vector: red for 0, blue for 1, grey for NA (if any)
  cols <- ifelse(is.na(vals), "grey70",
                 ifelse(vals == 1, "blue", "red"))
  
  # 7. Get D for this gene from results
  D_g <- results$D[results$gene == g]   # change 'gene' to the right col name if needed
  
  n_tips <- length(tree2$tip.label)
  tip_cex <- min(0.3, 0.5 / log10(n_tips)) 
  
  # 8. Open a high-resolution pdf device
  pdf(
    file   = file.path("out.dir/trees", paste0("tree_1", names_passed[i], ".pdf")),
    width  = 7,   # inches
    height = 7
  )
  
  # 9. Plot the tree: species-only labels, colored by trait, D in title
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

gene_names <- read.csv("genes - Sheet1 (3).csv")

negative_d <- negative_d %>%
  left_join(gene_names, by = "gene")

ordered_genes_n1prop <- negative_d %>%
  arrange(desc(n1_prop)) %>%
  dplyr::select(gene, symbol.x)
ordered_genes_pval <- negative_d %>%
  arrange(p_adjustp) %>%
  dplyr::select(gene, symbol.x)
ordered_genes_D <- negative_d %>%
  arrange(D) %>%
  dplyr::select(gene, symbol.x)

library(ape)

pdf(
  file = file.path("out.dir/trees", "all_gene_trees_D.pdf"),
  width = 7,
  height = 7
)

ordered_gene_ids <- ordered_genes_D$gene


for (i in seq_along(ordered_gene_ids)) {
  
  g <- ordered_gene_ids[i]
  symbol_g <- ordered_genes_D$symbol.x[i]
  
  vals_raw <- binary_data_wide1[[g]]
  
  is_null <- sapply(vals_raw, is.null)
  missing_species <- binary_data_wide1$species[is_null]
  tips_to_drop <- intersect(missing_species, tree$tip.label)
  
  tree2 <- if (length(tips_to_drop) > 0) drop.tip(tree, tips_to_drop) else tree
  
  vals_vec <- sapply(vals_raw, function(x) {
    if (is.null(x)) NA_real_ else as.numeric(x)
  })
  
  vals <- vals_vec[match(tree2$tip.label, binary_data_wide1$species)]
  
  cols <- ifelse(is.na(vals), "grey70",
                 ifelse(vals == 1, "blue", "red"))
  
  D_g <- results$D[results$gene == g]
  P <- results$p_adjustp[results$gene == g]
  
  n_tips <- length(tree2$tip.label)
  tip_cex <- min(0.3, 0.5 / log10(n_tips))
  
  plot(
    tree2,
    tip.color = cols,
    cex = tip_cex,
    main = paste0("Tree for ", symbol_g, "   (D = ", round(D_g, 7), ") (P = ", round(P, 6), ")")
  )
  
  legend(
    "topleft",
    legend = c("Trait = 0", "Trait = 1"),
    col = c("red", "blue"),
    pch = 19,
    pt.cex = 1.2,
    bty = "n"
  )
}

dev.off()



# Prioritization --------------------------------------------------------------

top_n <- 600

top_n1prop <- ordered_genes_n1prop %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::pull(symbol.x)

top_D <- ordered_genes_D %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::pull(symbol.x)

top_in_all_3 <- Reduce(intersect, list(top_n1prop, top_D))

top_in_all_3_df <- data.frame(symbol.x = top_in_all_3)

top_in_all_3_df


ranked_genes <- negative_d %>%
  dplyr::mutate(
    rank_D = rank(D, ties.method = "average"),
    rank_n1 = rank(-n1_prop, ties.method = "average"),
    
    combined_rank = (rank_D + rank_n1)/2
  ) %>%
  dplyr::arrange(combined_rank)

# Heatmap -----------------------------------------------------------------

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

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

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

# ORDERING FUNCTIONS 
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



# Pagel's Lambda -----------------------------------------------------------

library(geiger)

library(geiger)
library(dplyr)
library(purrr)


lambda_results <- data.frame(
  gene   = gene_cols,
  lambda = NA_real_,
  loglik = NA_real_,
  n      = NA_integer_,
  n1     = NA_integer_,
  n0     = NA_integer_,
  status = NA_character_,
  stringsAsFactors = FALSE
)

species_all <- binary_data_wide$species  # tiny micro-optimization

for (i in seq_along(gene_cols)) {   # change to gene_cols when ready
  
  g <- gene_cols[i]
  
  ## 1. Filter out "NULL"s for this gene
  keep <- binary_data_wide[[g]] != "NULL"
  
  species_keep <- species_all[keep]
  trait_raw    <- binary_data_wide[[g]][keep]   # "0"/"1"
  
  ## basic counts
  n_sp     <- length(trait_raw)
  trait_num <- as.integer(trait_raw)           # 0/1
  n1       <- sum(trait_num == 1L)
  n0       <- sum(trait_num == 0L)
  
  ## skip if too few species or no variation
  if (n_sp < 4 || n1 == 0 || n0 == 0) {
    lambda_results[i, c("lambda","loglik","n","n1","n0","status")] <-
      list(NA_real_, NA_real_, n_sp, n1, n0, "skipped")
    next
  }
  
  ## 2. Recode to 1/2 (required by geiger)
  trait_12 <- trait_num + 1L            # 1/2
  names(trait_12) <- species_keep
  
  ## 3. Match tree to species with data
  # we only need to drop species that *aren't* in species_keep
  tips_to_drop <- setdiff(tree$tip.label, species_keep)
  tree_sub <- if (length(tips_to_drop)) drop.tip(tree, tips_to_drop) else tree
  
  ## 4. Fit Pagel's lambda, catch errors
  fit_lambda <- tryCatch(
    fitDiscrete(tree_sub, trait_12,
                model     = "ARD",
                transform = "lambda"),
    error = function(e) e
  )
  
  # defaults
  status_i   <- "ok"
  lambda_val <- NA_real_
  loglik_val <- NA_real_
  
  if (inherits(fit_lambda, "error")) {
    status_i <- "error_fit"
  } else {
    lambda_raw <- fit_lambda$opt$lambda
    loglik_raw <- fit_lambda$opt$lnL
    
    # if value is NULL or empty, mark as error but still keep NA
    if (is.null(lambda_raw) || length(lambda_raw) == 0 ||
        is.null(loglik_raw) || length(loglik_raw) == 0) {
      status_i <- "error_null"
    } else {
      lambda_val <- as.numeric(lambda_raw[1])
      loglik_val <- as.numeric(loglik_raw[1])
    }
  }
  
  ## 5. Single write to results
  lambda_results[i, c("lambda","loglik","n","n1","n0","status")] <-
    list(lambda_val, loglik_val, n_sp, n1, n0, status_i)
}

lambda_results





for (i in seq_along(gene_cols[500:505])) { 
  
  g <- gene_cols[i]
  
  # 1. Filter out NULLs for this gene
  keep <- binary_data_wide[[g]] != "NULL"
  
  species_keep <- binary_data_wide$species[keep]
  trait_raw    <- binary_data_wide[[g]][keep]   # probably "0"/"1" or 0/1
  
  # 2. Recode to 1/2 (required by geiger)
  # if trait_raw is character "0"/"1":
  trait_num <- as.numeric(trait_raw)         # 0/1
  trait_12  <- trait_num + 1L                # 1/2
  
  # 3. Make named vector
  trait <- setNames(trait_12, species_keep)
  
  # 4. Match tree to species with data
  tips_to_drop <- setdiff(tree$tip.label, species_keep)
  tree_sub <- if (length(tips_to_drop)) drop.tip(tree, tips_to_drop) else tree
  
  # 5. Fit Pagel's lambda
  fit_lambda <- fitDiscrete(tree_sub, trait,
                            model     = "ARD",
                            transform = "lambda")
  
  print(fit_lambda$opt$lambda)
  
  print(fit_lambda$opt$lnL)
  
}






lambda_results_list <- vector("list", length(gene_cols))
idx <- 1  # separate index so we can `next` without leaving gaps

for (i in seq_along(gene_cols[1:5])) { 
  
  g <- gene_cols[i]
  
  # 1. Filter out NULLs for this gene
  keep <- binary_data_wide[[g]] != "NULL"
  
  species_keep <- binary_data_wide$species[keep]
  trait_raw    <- binary_data_wide[[g]][keep]   # "0"/"1" or 0/1
  
  n_species    <- length(species_keep)
  trait_num    <- as.numeric(trait_raw)         # 0/1
  prop_present <- if (n_species > 0) mean(trait_num) else NA_real_
  
  # --- optional but good: skip bad cases early ---
  # too few species or no variation → λ not identifiable / unstable
  if (n_species < 10 || length(unique(trait_num)) < 2) {
    next
  }
  
  # 2. Recode to 1/2 (required by geiger)
  trait_12  <- trait_num + 1L                # 1/2
  
  # 3. Make named vector
  trait <- setNames(trait_12, species_keep)
  
  # 4. Match tree to species with data
  tips_to_drop <- setdiff(tree$tip.label, species_keep)
  tree_sub <- if (length(tips_to_drop)) drop.tip(tree, tips_to_drop) else tree
  
  # 5. Fit Pagel's lambda, but *catch errors*
  fit_lambda <- tryCatch(
    fitDiscrete(tree_sub, trait,
                model     = "ARD",
                transform = "lambda"),
    error = function(e) e
  )
  
  # 6. Skip if there was an error
  if (inherits(fit_lambda, "error")) {
    next
  }
  
  # 7. Safely extract lambda and loglik; skip if they’re missing/empty
  lambda_val <- fit_lambda$opt$lambda
  loglik_val <- fit_lambda$loglik
  
  if (is.null(lambda_val) || length(lambda_val) == 0 ||
      is.null(loglik_val) || length(loglik_val) == 0) {
    next
  }
  
  lambda_val <- as.numeric(lambda_val[1])
  loglik_val <- as.numeric(loglik_val[1])
  
  # 8. Store result row
  lambda_results_list[[idx]] <- data.frame(
    gene         = g,
    lambda       = lambda_val,
    loglik       = loglik_val,
    n_species    = n_species,
    prop_present = prop_present,
    stringsAsFactors = FALSE
  )
  
  idx <- idx + 1
}

# 9. Bind everything into a single data frame
lambda_results <- do.call(rbind, lambda_results_list[1:(idx - 1)])


q <- binary_data_wide %>% 
  dplyr::select(species, FBgn0004170)
