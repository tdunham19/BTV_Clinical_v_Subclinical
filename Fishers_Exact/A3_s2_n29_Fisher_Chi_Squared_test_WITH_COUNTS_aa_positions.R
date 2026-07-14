###################################################################
# Aim 3
# Statistical assessment of aa associated with Clinical/Subclinical
# Running Fishers Exact and Chi Squared WITH COUNTS
# TD 03.31.26
###################################################################

# load required libraries
library(Biostrings)
library(dplyr)
library(tidyr)

# install Biostrings
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# 
# BiocManager::install("Biostrings")

#######################################################

# read in protein MSA
fasta_file <- "A3_s2_n29_CDS_Only_protein_aln.fasta"
msa <- readAAStringSet(fasta_file)

# convert alignment to matrix (rows = sequences, cols = positions)
msa_mat <- as.matrix(msa)

# extract labels from sequence names
seq_names <- names(msa)
labels <- sub("_.*$", "", seq_names)
if (!all(labels %in% c("Clinical", "Subclinical"))) {
  stop("Some sequences do not contain 'Clinical' or 'Subclinical' in the name.")
}

labels <- factor(labels, levels = c("Clinical", "Subclinical"))

# build long table: sample × position × amino acid
aa_long <- data.frame(
  sample   = rep(seq_names, times = ncol(msa_mat)),
  label    = rep(labels,    times = ncol(msa_mat)),
  position = rep(seq_len(ncol(msa_mat)), each = nrow(msa_mat)),
  aa       = as.vector(msa_mat),
  stringsAsFactors = FALSE
)

# remove gaps
aa_long <- aa_long[aa_long$aa != "-", ]

# save for downstream inspection
write.csv(aa_long, "A3_s2_n29_results_amino_acids_by_sample_and_position_071326_WITH_COUNTS.csv", row.names = FALSE)

# initialize results storage
results <- list()

# loop over alignment positions
for (pos in seq_len(ncol(msa_mat))) {
  
  aa_vec <- msa_mat[, pos]
  
  # remove gaps
  keep <- aa_vec != "-"
  aa_vec <- aa_vec[keep]
  lab_vec <- labels[keep]
  
  # skip conserved or degenerate positions
  if (length(unique(aa_vec)) < 2 || length(unique(lab_vec)) < 2) {
    next
  }
  
  # build contingency table
  tab <- table(aa_vec, lab_vec)
  
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    next
  }
  
  # Fisher's exact test
  ft <- fisher.test(tab)
  
  # chi-squared test
  ct <- suppressWarnings(chisq.test(tab))
  
  # convert table to data frame
  tab_df <- as.data.frame.matrix(tab)
  tab_df$aa <- rownames(tab_df)
  
  # rename columns to be explicit
  colnames(tab_df) <- gsub("Clinical", "Clinical_count", colnames(tab_df))
  colnames(tab_df) <- gsub("Subclinical", "Subclinical_count", colnames(tab_df))
  
  # pivot to wide format
  tab_wide <- tab_df %>%
    pivot_longer(cols = -aa, names_to = "group", values_to = "count") %>%
    unite(col_name, aa, group) %>%
    pivot_wider(names_from = col_name, values_from = count)
  
  # combine everything into one row
  results[[length(results) + 1]] <- cbind(
    data.frame(
      position     = pos,
      fisher_p     = ft$p.value,
      fisher_or    = if (!is.null(ft$estimate)) unname(ft$estimate) else NA,
      chisq_p      = ct$p.value,
      chisq_stat   = unname(ct$statistic),
      num_aa       = nrow(tab),
      stringsAsFactors = FALSE
    ),
    tab_wide
  )
}

# combine + adjust results
results_df <- bind_rows(results) %>%
  mutate(
    fisher_p_adj = p.adjust(fisher_p, method = "fdr"),
    chisq_p_adj  = p.adjust(chisq_p, method = "fdr")
  ) %>%
  arrange(fisher_p)

# write output to csv
write.csv(results_df, "A3_s2_n29_results_fisher_exact_chi_squared_WITH_COUNTS_by_position_071326.csv", row.names = FALSE)
