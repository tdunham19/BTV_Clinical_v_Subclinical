###################################################################
# Aim 3
# Statistical assessment of aa associated with Clinical/Subclinical
# TD 01.23.26
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
write.csv(aa_long, "A3_s2_n29_results_amino_acids_by_sample_and_position_013026.csv", row.names = FALSE)

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
  
  # store results
  results[[length(results) + 1]] <- data.frame(
    position   = pos,
    p_value    = ft$p.value,
    odds_ratio = if (!is.null(ft$estimate)) unname(ft$estimate) else NA,
    num_aa     = nrow(tab),
    stringsAsFactors = FALSE
  )
}

# combine + adjust results
results_df <- bind_rows(results) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_value)

# write output to csv
write.csv(results_df, "A3_s2_n29_results_fisher_exact_by_position_013026.csv", row.names = FALSE)
