#######################################################
# Nearest Neighbor Statistic Analysis 
# 
# This script reads in sets of sequences and calculates Snn
# Snn = nearest neighbor statistics 
# according to Hudson (2011) 
#
# Hudson RR. A new statistic for detecting genetic differentiation. 
# Genetics. 2000 Aug;155(4):2011-4. 
# doi: 10.1093/genetics/155.4.2011. PMID: 10924493; PMCID: PMC1461195.
#
# Original script written by: Mark Stenglein 8/28/2025
# Script adapted by: Tillie Dunham 04/01/2026 for Aim 3 Analysis

# A3 - Re-Run all segments TD 06.18.26
# Run with updated alignments and sequences submitted to GenBank

# All samples 
# s1 - n29
# s2 - n29
# s3 - n30
# s4 - n31
# s5 - n31
# s6 - n33
# s7 - n29
# s8 - n33
# s9 - n31
# s10 - n31

#######################################################

#Load Required Libraries
library(tidyverse)
library(readxl)
library(ggplot2)
library(patchwork)
library(ape)
library(profvis)

#######################################################

# CLINICAL/SUBCLINICAL STATUS

################################################################################

################################################################################
# Segment 1 Snn

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata <- metadata %>% select(accession = accession, group, segment_1)
metadata <- metadata %>% drop_na()

num_permutations <- 5000

process_alignment_s1 <-  function (
    fasta_msa = "./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta", 
    prefix="A3_s1_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s1_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn_s1(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn_s1(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.0, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s1 <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment1_snn <- process_alignment_s1("./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta", "Segment 1, n=29, Snn=")

Segment1_snn$plot


# END Segment 1 Snn
################################################################################


################################################################################
# Segment 2 Snn

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata <- metadata %>% select(accession = accession, group, segment_2)
metadata <- metadata %>% drop_na()

num_permutations <- 5000

process_alignment_s2 <-  function (
    fasta_msa = "./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta", 
    prefix="A3_s2_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s2_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn_s2(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn_s2(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.0, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s2 <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment2_snn <- process_alignment_s2("./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta", "Segment 2, n=29, Snn=")

Segment2_snn$plot


# END Segment 2 Snn
################################################################################


################################################################################
# Segment 3 Snn

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata <- metadata %>% select(accession = accession, group, segment_3)
metadata <- metadata %>% drop_na()

num_permutations <- 5000

process_alignment_s3 <-  function (
    fasta_msa = "./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s3_n30_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s3_n30_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn_s3(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn_s3(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.0, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s3 <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment3_snn <- process_alignment_s3("./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta", "Segment 3, n=30, Snn=")

Segment3_snn$plot


# END Segment 3 Snn
################################################################################


################################################################################
# Segment 4 Snn

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata <- metadata %>% select(accession = accession, group, segment_4)
metadata <- metadata %>% drop_na()

num_permutations <- 5000

process_alignment_s4 <-  function (
    fasta_msa = "./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s4_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s4_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn_s4(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn_s4(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.0, y=350, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s4 <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment4_snn <- process_alignment_s4("./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta", "Segment 4, n=31, Snn=")

Segment4_snn$plot


# END Segment 4 Snn
################################################################################


################################################################################
# Segment 5 Snn

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata <- metadata %>% select(accession = accession, group, segment_5)
metadata <- metadata %>% drop_na()

num_permutations <- 5000

process_alignment_s5 <-  function (
    fasta_msa = "./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s5_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s5_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn_s5(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn_s5(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.0, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s5 <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment5_snn <- process_alignment_s5("./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta", "Segment 5, n=31, Snn=")

Segment5_snn$plot


# END Segment 5 Snn
################################################################################


################################################################################
# Segment 6 Snn

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata <- metadata %>% select(accession = accession, group, segment_6)
metadata <- metadata %>% drop_na()

num_permutations <- 5000

process_alignment_s6 <-  function (
    fasta_msa = "./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s6_n33_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s6_n33_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn_s6(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn_s6(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.0, y=250, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s6 <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment6_snn <- process_alignment_s6("./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta", "Segment 6, n=33, Snn=")

Segment6_snn$plot

# END Segment 6 Snn
################################################################################


################################################################################
# Segment 7 Snn

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata <- metadata %>% select(accession = accession, group, segment_7)
metadata <- metadata %>% drop_na()

num_permutations <- 5000

process_alignment_s7 <-  function (
    fasta_msa = "./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s7_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s7_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn_s7(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn_s7(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.0, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s7 <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment7_snn <- process_alignment_s7("./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta", "Segment 7, n=29, Snn=")

Segment7_snn$plot

# END Segment 7 Snn
################################################################################


################################################################################
# Segment 8 Snn

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata <- metadata %>% select(accession = accession, group, segment_8)
metadata <- metadata %>% drop_na()

num_permutations <- 5000

process_alignment_s8 <-  function (
    fasta_msa = "./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s8_n33_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s8_n33_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn_s8(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn_s8(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.0, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s8 <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment8_snn <- process_alignment_s8("./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta", "Segment 8, n=33, Snn=")

Segment8_snn$plot


# END Segment 8 Snn
################################################################################


################################################################################
# Segment 9 Snn

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata <- metadata %>% select(accession = accession, group, segment_9)
metadata <- metadata %>% drop_na()

num_permutations <- 5000

process_alignment_s9 <-  function (
    fasta_msa = "./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s9_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s9_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn_s9(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn_s9(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.0, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s9 <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment9_snn <- process_alignment_s9("./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta", "Segment 9, n=31, Snn=")

Segment9_snn$plot


# END Segment 9 Snn
################################################################################


################################################################################
# Segment 10 Snn

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata <- metadata %>% select(accession = accession, group, segment_10)
metadata <- metadata %>% drop_na()

num_permutations <- 5000

process_alignment_s10 <-  function (
    fasta_msa = "./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s10_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s10_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn_s10(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn_s10(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.0, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s10 <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment10_snn <- process_alignment_s10("./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta", "Segment 10, n=31, Snn=")

Segment10_snn$plot


# END Segment 10 Snn
################################################################################


################################################################################
# Putting it all together 

# -------------------------------
# Create table
# -------------------------------
snn_supp_table <- 
  tibble( segment = c("BTV Segment 1", "BTV Segment 2", "BTV Segment 3", 
                      "BTV Segment 4", "BTV Segment 5", "BTV Segment 6", 
                      "BTV Segment 7", "BTV Segment 8", "BTV Segment 9", "BTV Segment 10"),
                n = c("s1 = n29", "s2 = n29", "s3 = n30",
                      "s4 = n31", "s5 = n31", "s6 = n33", 
                      "s7 = n29", "s8 = n33", "s9 = n31", "s10 = n31"),
          snn  = c(Segment1_snn$snn, Segment2_snn$snn, Segment3_snn$snn, 
                   Segment4_snn$snn, Segment5_snn$snn, Segment6_snn$snn, 
                   Segment7_snn$snn, Segment8_snn$snn, Segment9_snn$snn, Segment10_snn$snn),
          pval = c(Segment1_snn$pval, Segment2_snn$pval, Segment3_snn$pval, 
                   Segment4_snn$pval, Segment5_snn$pval, Segment6_snn$pval, 
                   Segment7_snn$pval, Segment8_snn$pval, Segment9_snn$pval, Segment10_snn$pval)) %>%
  mutate(pval_bonf = p.adjust(pval, method = "bonferroni", n = length(pval)))

write.table(snn_supp_table, 
            file="rerun_NO_Clinical_20_GenBank_Seqs_STATUS_supplemental_table_snn_061826.txt",
            quote=F, sep="\t", row.names=F)
################################################################################ END















################################################################################


















################################################################################


# LURKING VARIABLES TO TEST: SPECIES, AGE, SEX, YEAR COLLECTED, etc. 


################################################################################


################################################################################


# SEGMENT 1


################################################################################
# Segment 1 Snn - testing species as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_species <- metadata %>% select(accession = accession, species, segment_1)
metadata_species <- metadata_species %>% drop_na()

num_permutations <- 5000

process_alignment_s1_species <-  function (
    fasta_msa = "./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta", 
    prefix="A3_s1_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s1_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_species)  
  # create a vector of the group matching the order of the distance matrix
  species <- acc_loc %>% pull(species)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s1_species(msa_dist, species)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(species, replace=F)
    perm_snn             <- calculate_snn_s1_species(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              #  pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s1_species <- function(distance_matrix, species_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(species_vector[nn_i] == species_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment1_species_snn <- process_alignment_s1_species("./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta", "Segment 1 - Species, n=29, Snn=")

Segment1_species_snn$plot

# END Segment 1 Snn - testing species as lurking variable
################################################################################


################################################################################
# Segment 1 Snn - testing serotype as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_serotype <- metadata %>% select(accession = accession, serotype, segment_1)
metadata_serotype <- metadata_serotype %>% drop_na()

num_permutations <- 5000

process_alignment_s1_serotype <-  function (
    fasta_msa = "./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta", 
    prefix="A3_s1_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s1_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_serotype)  
  # create a vector of the group matching the order of the distance matrix
  serotypes <- acc_loc %>% pull(serotype)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s1_serotype(msa_dist, serotypes)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(serotypes, replace=F)
    perm_snn             <- calculate_snn_s1_serotype(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.3, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s1_serotype <- function(distance_matrix, serotype_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(serotype_vector[nn_i] == serotype_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment1_serotype_snn <- process_alignment_s1_serotype("./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta", "Segment 1 - Serotype, n=29, Snn=")

Segment1_serotype_snn$plot

# END Segment 1 Snn - testing serotype as lurking variable
################################################################################


################################################################################
# Segment 1 Snn - testing year as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_year <- metadata %>% select(accession = accession, year, segment_1)
metadata_year <- metadata_year %>% drop_na()

num_permutations <- 5000

process_alignment_s1_year <-  function (
    fasta_msa = "./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta", 
    prefix="A3_s1_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s1_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_year)  
  # create a vector of the group matching the order of the distance matrix
  years <- acc_loc %>% pull(year)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s1_year(msa_dist, years)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(years, replace=F)
    perm_snn             <- calculate_snn_s1_year(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.3, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s1_year <- function(distance_matrix, year_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(year_vector[nn_i] == year_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment1_year_snn <- process_alignment_s1_year("./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta", "Segment 1 - Year, n=29, Snn=")

Segment1_year_snn$plot

# END Segment 1 Snn - testing year as lurking variable
################################################################################


################################################################################
# Segment 1 Snn - testing age as lurking variable & no NAs
# There is not enough age representation across groups to reliably calculate Snn (TD 10.03.25)

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_age <- metadata %>% select(accession = accession, age, segment_1)
metadata_age <- metadata_age %>% drop_na()

num_permutations <- 5000

process_alignment_s1_age <-  function (
    fasta_msa = "./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta", 
    prefix="A3_s1_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s1_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_age)  
  # create a vector of the group matching the order of the distance matrix
  ages <- acc_loc %>% pull(age)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s1_age(msa_dist, ages)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(ages, replace=F)
    perm_snn             <- calculate_snn_s1_age(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- 1/num_permutations
    this_snn_p_val <- (this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.25, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s1_age <- function(distance_matrix, age_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(age_vector[nn_i] == age_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment1_age_snn <- process_alignment_s1_age("./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta", "Segment 1 - Animal Age, n=29, Snn=")

Segment1_age_snn$plot

# END Segment 1 Snn - testing age as lurking variable
################################################################################


################################################################################
# Segment 1 Snn - testing sex as lurking variable & no NAs

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_sex <- metadata %>% select(accession = accession, sex, segment_1)
metadata_sex <- metadata_sex %>% drop_na()

num_permutations <- 5000

process_alignment_s1_sex <-  function (
    fasta_msa = "./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta", 
    prefix="A3_s1_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s1_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_sex)  
  # create a vector of the group matching the order of the distance matrix
  sexs <- acc_loc %>% pull(sex)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s1_sex(msa_dist, sexs)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(sexs, replace=F)
    perm_snn             <- calculate_snn_s1_sex(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s1_sex <- function(distance_matrix, sex_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(sex_vector[nn_i] == sex_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment1_sex_snn <- process_alignment_s1_sex("./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta", "Segment 1 - Animal Sex, n=29, Snn=")

Segment1_sex_snn$plot

# END Segment 1 Snn - testing sex as lurking variable
################################################################################


################################################################################
# Segment 1 Snn - testing state as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_state <- metadata %>% select(accession = accession, state, segment_1)
metadata_state <- metadata_state %>% drop_na(segment_1)

num_permutations <- 5000

process_alignment_s1_state <-  function (
    fasta_msa = "./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta", 
    prefix="A3_s1_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s1_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_state)  
  # create a vector of the group matching the order of the distance matrix
  states <- acc_loc %>% pull(state)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s1_state(msa_dist, states)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(states, replace=F)
    perm_snn             <- calculate_snn_s1_state(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s1_state <- function(distance_matrix, state_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(state_vector[nn_i] == state_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment1_state_snn <- process_alignment_s1_state("./Alignments/A3_s1_n29_CDS_Only_aln_repeat_2.fasta", "Segment 1 - State, n=29, Snn=")

Segment1_state_snn$plot

# END Segment 1 Snn - testing state as lurking variable
################################################################################

################################################################################

# END SEGMENT 1


################################################################################













################################################################################


# SEGMENT 2


################################################################################


################################################################################
# Segment 2 Snn - testing species as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_species <- metadata %>% select(accession = accession, species, segment_2)
metadata_species <- metadata_species %>% drop_na()

num_permutations <- 5000

process_alignment_s2_species <-  function (
    fasta_msa = "./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta", 
    prefix="A3_s2_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s2_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_species)  
  # create a vector of the group matching the order of the distance matrix
  species <- acc_loc %>% pull(species)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s2_species(msa_dist, species)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(species, replace=F)
    perm_snn             <- calculate_snn_s2_species(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s2_species <- function(distance_matrix, species_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(species_vector[nn_i] == species_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment2_species_snn <- process_alignment_s2_species("./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta", "Segment 2 - Species, n=29, Snn=")

Segment2_species_snn$plot

# END Segment 2 Snn - testing species as lurking variable
################################################################################


################################################################################
# Segment 2 Snn - testing serotype as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_serotype <- metadata %>% select(accession = accession, serotype, segment_2)
metadata_serotype <- metadata_serotype %>% drop_na()

num_permutations <- 5000

process_alignment_s2_serotype <-  function (
    fasta_msa = "./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta", 
    prefix="A3_s2_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s2_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_serotype)  
  # create a vector of the group matching the order of the distance matrix
  serotypes <- acc_loc %>% pull(serotype)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s2_serotype(msa_dist, serotypes)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(serotypes, replace=F)
    perm_snn             <- calculate_snn_s2_serotype(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.3, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s2_serotype <- function(distance_matrix, serotype_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(serotype_vector[nn_i] == serotype_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment2_serotype_snn <- process_alignment_s2_serotype("./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta", "Segment 2 - Serotype, n=29, Snn=")

Segment2_serotype_snn$plot

# END Segment 2 Snn - testing serotype as lurking variable
################################################################################


################################################################################
# Segment 2 Snn - testing year as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_year <- metadata %>% select(accession = accession, year, segment_2)
metadata_year <- metadata_year %>% drop_na()

num_permutations <- 5000

process_alignment_s2_year <-  function (
    fasta_msa = "./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta", 
    prefix="A3_s2_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s2_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_year)  
  # create a vector of the group matching the order of the distance matrix
  years <- acc_loc %>% pull(year)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s2_year(msa_dist, years)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(years, replace=F)
    perm_snn             <- calculate_snn_s2_year(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.3, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s2_year <- function(distance_matrix, year_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(year_vector[nn_i] == year_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment2_year_snn <- process_alignment_s2_year("./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta", "Segment 2 - Year, n=29, Snn=")

Segment2_year_snn$plot

# END Segment 2 Snn - testing year as lurking variable
################################################################################


################################################################################
# Segment 2 Snn - testing age as lurking variable
# There is not enough age representation across groups to reliably calculate Snn (TD 10.03.25)

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need
metadata_age <- metadata %>% select(accession = accession, age, segment_2)
metadata_age <- metadata_age %>% drop_na()

num_permutations <- 5000

process_alignment_s2_age <-  function (
    fasta_msa = "./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta", 
    prefix="A3_s2_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s2_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_age)  
  # create a vector of the group matching the order of the distance matrix
  ages <- acc_loc %>% pull(age)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s2_age(msa_dist, ages)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(ages, replace=F)
    perm_snn             <- calculate_snn_s2_age(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- 1/num_permutations
    this_snn_p_val <- (this_snn_p_val_sci)
  }
  
  this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val_bonf)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.25, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val_bonf,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s2_age <- function(distance_matrix, age_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(age_vector[nn_i] == age_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment2_age_snn <- process_alignment_s2_age("./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta", "Segment 2 - Animal Age, n=29, Snn=")

Segment2_age_snn$plot

# END Segment 2 Snn - testing age as lurking variable
################################################################################


################################################################################
# Segment 2 Snn - testing sex as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_sex <- metadata %>% select(accession = accession, sex, segment_2)
metadata_sex <- metadata_sex %>% drop_na()

num_permutations <- 5000

process_alignment_s2_sex <-  function (
    fasta_msa = "./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta", 
    prefix="A3_s2_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s2_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_sex)  
  # create a vector of the group matching the order of the distance matrix
  sexs <- acc_loc %>% pull(sex)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s2_sex(msa_dist, sexs)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(sexs, replace=F)
    perm_snn             <- calculate_snn_s2_sex(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s2_sex <- function(distance_matrix, sex_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(sex_vector[nn_i] == sex_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment2_sex_snn <- process_alignment_s2_sex("./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta", "Segment 2 - Animal Sex, n=29, Snn=")

Segment2_sex_snn$plot

# END Segment 2 Snn - testing sex as lurking variable
################################################################################


################################################################################
# Segment 2 Snn - testing state as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_state <- metadata %>% select(accession = accession, state, segment_2)
metadata_state <- metadata_state %>% drop_na(segment_2)

num_permutations <- 5000

process_alignment_s2_state <-  function (
    fasta_msa = "./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta", 
    prefix="A3_s2_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s2_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_state)  
  # create a vector of the group matching the order of the distance matrix
  states <- acc_loc %>% pull(state)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s2_state(msa_dist, states)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(states, replace=F)
    perm_snn             <- calculate_snn_s2_state(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s2_state <- function(distance_matrix, state_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(state_vector[nn_i] == state_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment2_state_snn <- process_alignment_s2_state("./Alignments/A3_s2_n29_CDS_Only_aln_repeat_2.fasta", "Segment 2 - State, n=29, Snn=")

Segment2_state_snn$plot

# END Segment 2 Snn - testing state as lurking variable
################################################################################


################################################################################


# END SEGMENT 2


################################################################################




















################################################################################


# SEGMENT 3


################################################################################


################################################################################
# Segment 3 Snn - testing species as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_species <- metadata %>% select(accession = accession, species, segment_3)
metadata_species <- metadata_species %>% drop_na()

num_permutations <- 5000

process_alignment_s3_species <-  function (
    fasta_msa = "./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s3_n30_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s3_n30_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_species)  
  # create a vector of the group matching the order of the distance matrix
  species <- acc_loc %>% pull(species)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s3_species(msa_dist, species)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(species, replace=F)
    perm_snn             <- calculate_snn_s3_species(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s3_species <- function(distance_matrix, species_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(species_vector[nn_i] == species_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment3_species_snn <- process_alignment_s3_species("./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta", "Segment 3 - Species, n=30, Snn=")

Segment3_species_snn$plot

# END Segment 3 Snn - testing species as lurking variable
################################################################################


################################################################################
# Segment 3 Snn - testing serotype as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_serotype <- metadata %>% select(accession = accession, serotype, segment_3)
metadata_serotype <- metadata_serotype %>% drop_na()

num_permutations <- 5000

process_alignment_s3_serotype <-  function (
    fasta_msa = "./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s3_n30_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s3_n30_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_serotype)  
  # create a vector of the group matching the order of the distance matrix
  serotypes <- acc_loc %>% pull(serotype)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s3_serotype(msa_dist, serotypes)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(serotypes, replace=F)
    perm_snn             <- calculate_snn_s3_serotype(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.3, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s3_serotype <- function(distance_matrix, serotype_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(serotype_vector[nn_i] == serotype_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment3_serotype_snn <- process_alignment_s3_serotype("./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta", "Segment 3 - Serotype, n=30, Snn=")

Segment3_serotype_snn$plot

# END Segment 3 Snn - testing serotype as lurking variable
################################################################################


################################################################################
# Segment 3 Snn - testing year as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_year <- metadata %>% select(accession = accession, year, segment_3)
metadata_year <- metadata_year %>% drop_na()

num_permutations <- 5000

process_alignment_s3_year <-  function (
    fasta_msa = "./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s3_n30_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s3_n30_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_year)  
  # create a vector of the group matching the order of the distance matrix
  years <- acc_loc %>% pull(year)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s3_year(msa_dist, years)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(years, replace=F)
    perm_snn             <- calculate_snn_s3_year(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s3_year <- function(distance_matrix, year_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(year_vector[nn_i] == year_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment3_year_snn <- process_alignment_s3_year("./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta", "Segment 3 - Year, n=30, Snn=")

Segment3_year_snn$plot

# END Segment 3 Snn - testing year as lurking variable
################################################################################


################################################################################
# Segment 3 Snn - testing age as lurking variable
# There is not enough age representation across groups to reliably calculate Snn (TD 10.03.25)

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need
metadata_age <- metadata %>% select(accession = accession, age, segment_3)
metadata_age <- metadata_age %>% drop_na()

num_permutations <- 5000

process_alignment_s3_age <-  function (
    fasta_msa = "./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s3_n30_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s3_n30_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_age)  
  # create a vector of the group matching the order of the distance matrix
  ages <- acc_loc %>% pull(age)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s3_age(msa_dist, ages)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(ages, replace=F)
    perm_snn             <- calculate_snn_s3_age(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- 1/num_permutations
    this_snn_p_val <- (this_snn_p_val_sci)
  }
  
  this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val_bonf)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.25, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val_bonf,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s3_age <- function(distance_matrix, age_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(age_vector[nn_i] == age_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment3_age_snn <- process_alignment_s3_age("./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta", "Segment 3 - Animal Age, n=30, Snn=")

Segment3_age_snn$plot

# END Segment 3 Snn - testing age as lurking variable
################################################################################


################################################################################
# Segment 3 Snn - testing sex as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_sex <- metadata %>% select(accession = accession, sex, segment_3)
metadata_sex <- metadata_sex %>% drop_na()

num_permutations <- 5000

process_alignment_s3_sex <-  function (
    fasta_msa = "./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s3_n30_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s3_n30_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_sex)  
  # create a vector of the group matching the order of the distance matrix
  sexs <- acc_loc %>% pull(sex)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s3_sex(msa_dist, sexs)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(sexs, replace=F)
    perm_snn             <- calculate_snn_s3_sex(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s3_sex <- function(distance_matrix, sex_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(sex_vector[nn_i] == sex_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment3_sex_snn <- process_alignment_s3_sex("./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta", "Segment 3 - Animal Sex, n=30, Snn=")

Segment3_sex_snn$plot

# END Segment 3 Snn - testing sex as lurking variable
################################################################################

################################################################################
# Segment 3 Snn - testing state as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_state <- metadata %>% select(accession = accession, state, segment_3)
metadata_state <- metadata_state %>% drop_na()

num_permutations <- 5000

process_alignment_s3_state <-  function (
    fasta_msa = "./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s3_n30_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s3_n30_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_state)  
  # create a vector of the group matching the order of the distance matrix
  states <- acc_loc %>% pull(state)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s3_state(msa_dist, states)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(states, replace=F)
    perm_snn             <- calculate_snn_s3_state(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s3_state <- function(distance_matrix, state_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(state_vector[nn_i] == state_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment3_state_snn <- process_alignment_s3_state("./Alignments/A3_s3_n30_CDS_Only_aln_repeat_2.fasta", "Segment 3 - State, n=30, Snn=")

Segment3_state_snn$plot

# END Segment 3 Snn - testing state as lurking variable
################################################################################


################################################################################


# END SEGMENT 3


################################################################################



















################################################################################


# SEGMENT 4


################################################################################


################################################################################
# Segment 4 Snn - testing species as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_species <- metadata %>% select(accession = accession, species, segment_4)
metadata_species <- metadata_species %>% drop_na()

num_permutations <- 5000

process_alignment_s4_species <-  function (
    fasta_msa = "./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s4_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s4_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_species)  
  # create a vector of the group matching the order of the distance matrix
  species <- acc_loc %>% pull(species)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s4_species(msa_dist, species)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(species, replace=F)
    perm_snn             <- calculate_snn_s4_species(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s4_species <- function(distance_matrix, species_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(species_vector[nn_i] == species_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment4_species_snn <- process_alignment_s4_species("./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta", "Segment 4 - Species, n=31, Snn=")

Segment4_species_snn$plot

# END Segment 4 Snn - testing species as lurking variable
################################################################################


################################################################################
# Segment 4 Snn - testing serotype as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_serotype <- metadata %>% select(accession = accession, serotype, segment_4)
metadata_serotype <- metadata_serotype %>% drop_na()

num_permutations <- 5000

process_alignment_s4_serotype <-  function (
    fasta_msa = "./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s4_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s4_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_serotype)  
  # create a vector of the group matching the order of the distance matrix
  serotypes <- acc_loc %>% pull(serotype)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s4_serotype(msa_dist, serotypes)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(serotypes, replace=F)
    perm_snn             <- calculate_snn_s4_serotype(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.3, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s4_serotype <- function(distance_matrix, serotype_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(serotype_vector[nn_i] == serotype_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment4_serotype_snn <- process_alignment_s4_serotype("./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta", "Segment 4 - Serotype, n=31, Snn=")

Segment4_serotype_snn$plot

# END Segment 4 Snn - testing serotype as lurking variable
################################################################################


################################################################################
# Segment 4 Snn - testing year as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_year <- metadata %>% select(accession = accession, year, segment_4)
metadata_year <- metadata_year %>% drop_na()

num_permutations <- 5000

process_alignment_s4_year <-  function (
    fasta_msa = "./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s4_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s4_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_year)  
  # create a vector of the group matching the order of the distance matrix
  years <- acc_loc %>% pull(year)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s4_year(msa_dist, years)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(years, replace=F)
    perm_snn             <- calculate_snn_s4_year(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s4_year <- function(distance_matrix, year_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(year_vector[nn_i] == year_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment4_year_snn <- process_alignment_s4_year("./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta", "Segment 4 - Year, n=31, Snn=")

Segment4_year_snn$plot

# END Segment 4 Snn - testing year as lurking variable
################################################################################


################################################################################
# Segment 4 Snn - testing age as lurking variable
# There is not enough age representation across groups to reliably calculate Snn (TD 10.03.25)

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need
metadata_age <- metadata %>% select(accession = accession, age, segment_4)
metadata_age <- metadata_age %>% drop_na()

num_permutations <- 5000

process_alignment_s4_age <-  function (
    fasta_msa = "./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s4_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s4_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_age)  
  # create a vector of the group matching the order of the distance matrix
  ages <- acc_loc %>% pull(age)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s4_age(msa_dist, ages)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(ages, replace=F)
    perm_snn             <- calculate_snn_s4_age(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- 1/num_permutations
    this_snn_p_val <- (this_snn_p_val_sci)
  }
  
  this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val_bonf)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.25, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val_bonf,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s4_age <- function(distance_matrix, age_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(age_vector[nn_i] == age_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment4_age_snn <- process_alignment_s4_age("./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta", "Segment 4 - Animal Age, n=31, Snn=")

Segment4_age_snn$plot

# END Segment 4 Snn - testing age as lurking variable
################################################################################


################################################################################
# Segment 4 Snn - testing sex as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_sex <- metadata %>% select(accession = accession, sex, segment_4)
metadata_sex <- metadata_sex %>% drop_na()

num_permutations <- 5000

process_alignment_s4_sex <-  function (
    fasta_msa = "./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s4_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s4_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_sex)  
  # create a vector of the group matching the order of the distance matrix
  sexs <- acc_loc %>% pull(sex)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s4_sex(msa_dist, sexs)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(sexs, replace=F)
    perm_snn             <- calculate_snn_s4_sex(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s4_sex <- function(distance_matrix, sex_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(sex_vector[nn_i] == sex_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment4_sex_snn <- process_alignment_s4_sex("./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta", "Segment 4 - Animal Sex, n=31, Snn=")

Segment4_sex_snn$plot

# END Segment 4 Snn - testing sex as lurking variable
################################################################################


################################################################################
# Segment 4 Snn - testing state as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_state <- metadata %>% select(accession = accession, state, segment_4)
metadata_state <- metadata_state %>% drop_na()

num_permutations <- 5000

process_alignment_s4_state <-  function (
    fasta_msa = "./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s4_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s4_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_state)  
  # create a vector of the group matching the order of the distance matrix
  states <- acc_loc %>% pull(state)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s4_state(msa_dist, states)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(states, replace=F)
    perm_snn             <- calculate_snn_s4_state(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s4_state <- function(distance_matrix, state_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(state_vector[nn_i] == state_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment4_state_snn <- process_alignment_s4_state("./Alignments/A3_s4_n31_CDS_Only_aln_repeat_2.fasta", "Segment 4 - State, n=31, Snn=")

Segment4_state_snn$plot

# END Segment 4 Snn - testing state as lurking variable
################################################################################


################################################################################


# END SEGMENT 4


################################################################################

















################################################################################


# SEGMENT 5


################################################################################


################################################################################
# Segment 5 Snn - testing species as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_species <- metadata %>% select(accession = accession, species, segment_5)
metadata_species <- metadata_species %>% drop_na()

num_permutations <- 5000

process_alignment_s5_species <-  function (
    fasta_msa = "./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s5_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s5_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_species)  
  # create a vector of the group matching the order of the distance matrix
  species <- acc_loc %>% pull(species)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s5_species(msa_dist, species)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(species, replace=F)
    perm_snn             <- calculate_snn_s5_species(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s5_species <- function(distance_matrix, species_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(species_vector[nn_i] == species_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment5_species_snn <- process_alignment_s5_species("./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta", "Segment 5 - Species, n=31, Snn=")

Segment5_species_snn$plot

# END Segment 5 Snn - testing species as lurking variable
################################################################################


################################################################################
# Segment 5 Snn - testing serotype as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_serotype <- metadata %>% select(accession = accession, serotype, segment_5)
metadata_serotype <- metadata_serotype %>% drop_na()

num_permutations <- 5000

process_alignment_s5_serotype <-  function (
    fasta_msa = "./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s5_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s5_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_serotype)  
  # create a vector of the group matching the order of the distance matrix
  serotypes <- acc_loc %>% pull(serotype)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s5_serotype(msa_dist, serotypes)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(serotypes, replace=F)
    perm_snn             <- calculate_snn_s5_serotype(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.3, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s5_serotype <- function(distance_matrix, serotype_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(serotype_vector[nn_i] == serotype_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment5_serotype_snn <- process_alignment_s5_serotype("./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta", "Segment 5 - Serotype, n=31, Snn=")

Segment5_serotype_snn$plot

# END Segment 5 Snn - testing serotype as lurking variable
################################################################################


################################################################################
# Segment 5 Snn - testing year as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_year <- metadata %>% select(accession = accession, year, segment_5)
metadata_year <- metadata_year %>% drop_na()

num_permutations <- 5000

process_alignment_s5_year <-  function (
    fasta_msa = "./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s5_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s5_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_year)  
  # create a vector of the group matching the order of the distance matrix
  years <- acc_loc %>% pull(year)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s5_year(msa_dist, years)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(years, replace=F)
    perm_snn             <- calculate_snn_s5_year(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s5_year <- function(distance_matrix, year_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(year_vector[nn_i] == year_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment5_year_snn <- process_alignment_s5_year("./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta", "Segment 5 - Year, n=31, Snn=")

Segment5_year_snn$plot

# END Segment 5 Snn - testing year as lurking variable
################################################################################


################################################################################
# Segment 5 Snn - testing age as lurking variable
# There is not enough age representation across groups to reliably calculate Snn (TD 10.03.25)

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need
metadata_age <- metadata %>% select(accession = accession, age, segment_5)
metadata_age <- metadata_age %>% drop_na()

num_permutations <- 5000

process_alignment_s5_age <-  function (
    fasta_msa = "./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s5_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s5_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_age)  
  # create a vector of the group matching the order of the distance matrix
  ages <- acc_loc %>% pull(age)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s5_age(msa_dist, ages)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(ages, replace=F)
    perm_snn             <- calculate_snn_s5_age(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- 1/num_permutations
    this_snn_p_val <- (this_snn_p_val_sci)
  }
  
  this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val_bonf)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.25, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val_bonf,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s5_age <- function(distance_matrix, age_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(age_vector[nn_i] == age_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment5_age_snn <- process_alignment_s5_age("./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta", "Segment 5 - Animal Age, n=31, Snn=")

Segment5_age_snn$plot

# END Segment 5 Snn - testing age as lurking variable
################################################################################


################################################################################
# Segment 5 Snn - testing sex as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_sex <- metadata %>% select(accession = accession, sex, segment_5)
metadata_sex <- metadata_sex %>% drop_na()

num_permutations <- 5000

process_alignment_s5_sex <-  function (
    fasta_msa = "./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s5_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s5_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_sex)  
  # create a vector of the group matching the order of the distance matrix
  sexs <- acc_loc %>% pull(sex)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s5_sex(msa_dist, sexs)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(sexs, replace=F)
    perm_snn             <- calculate_snn_s5_sex(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s5_sex <- function(distance_matrix, sex_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(sex_vector[nn_i] == sex_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment5_sex_snn <- process_alignment_s5_sex("./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta", "Segment 5 - Animal Sex, n=31, Snn=")

Segment5_sex_snn$plot

# END Segment 5 Snn - testing sex as lurking variable
################################################################################


################################################################################
# Segment 5 Snn - testing state as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_state <- metadata %>% select(accession = accession, state, segment_5)
metadata_state <- metadata_state %>% drop_na()

num_permutations <- 5000

process_alignment_s5_state <-  function (
    fasta_msa = "./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s5_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s5_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_state)  
  # create a vector of the group matching the order of the distance matrix
  states <- acc_loc %>% pull(state)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s5_state(msa_dist, states)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(states, replace=F)
    perm_snn             <- calculate_snn_s5_state(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s5_state <- function(distance_matrix, state_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(state_vector[nn_i] == state_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment5_state_snn <- process_alignment_s5_state("./Alignments/A3_s5_n31_CDS_Only_aln_repeat_2.fasta", "Segment 5 - State, n=31, Snn=")

Segment5_state_snn$plot

# END Segment 5 Snn - testing state as lurking variable
################################################################################


################################################################################


# END SEGMENT 5


################################################################################




















################################################################################


# SEGMENT 6


################################################################################


################################################################################
# Segment 6 Snn - testing species as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_species <- metadata %>% select(accession = accession, species, segment_6)
metadata_species <- metadata_species %>% drop_na()

num_permutations <- 5000

process_alignment_s6_species <-  function (
    fasta_msa = "./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s6_n33_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s6_n33_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_species)  
  # create a vector of the group matching the order of the distance matrix
  species <- acc_loc %>% pull(species)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s6_species(msa_dist, species)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(species, replace=F)
    perm_snn             <- calculate_snn_s6_species(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s6_species <- function(distance_matrix, species_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(species_vector[nn_i] == species_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment6_species_snn <- process_alignment_s6_species("./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta", "Segment 6 - Species, n=33, Snn=")


Segment6_species_snn$plot

# END Segment 6 Snn - testing species as lurking variable
################################################################################


################################################################################
# Segment 6 Snn - testing serotype as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_serotype <- metadata %>% select(accession = accession, serotype, segment_6)
metadata_serotype <- metadata_serotype %>% drop_na()

num_permutations <- 5000

process_alignment_s6_serotype <-  function (
    fasta_msa = "./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s6_n33_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s6_n33_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_serotype)  
  # create a vector of the group matching the order of the distance matrix
  serotypes <- acc_loc %>% pull(serotype)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s6_serotype(msa_dist, serotypes)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(serotypes, replace=F)
    perm_snn             <- calculate_snn_s6_serotype(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.3, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s6_serotype <- function(distance_matrix, serotype_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(serotype_vector[nn_i] == serotype_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment6_serotype_snn <- process_alignment_s6_serotype("./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta", "Segment 6 - Serotype, n=33, Snn=")

Segment6_serotype_snn$plot

# END Segment 6 Snn - testing serotype as lurking variable
################################################################################


################################################################################
# Segment 6 Snn - testing year as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_year <- metadata %>% select(accession = accession, year, segment_6)
metadata_year <- metadata_year %>% drop_na()

num_permutations <- 5000

process_alignment_s6_year <-  function (
    fasta_msa = "./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s6_n33_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s6_n33_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_year)  
  # create a vector of the group matching the order of the distance matrix
  years <- acc_loc %>% pull(year)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s6_year(msa_dist, years)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(years, replace=F)
    perm_snn             <- calculate_snn_s6_year(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.3, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s6_year <- function(distance_matrix, year_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(year_vector[nn_i] == year_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment6_year_snn <- process_alignment_s6_year("./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta", "Segment 6 - Year, n=33, Snn=")

Segment6_year_snn$plot

# END Segment 6 Snn - testing year as lurking variable
################################################################################


################################################################################
# Segment 6 Snn - testing age as lurking variable
# There is not enough age representation across groups to reliably calculate Snn (TD 10.03.25)

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need
metadata_age <- metadata %>% select(accession = accession, age, segment_6)
metadata_age <- metadata_age %>% drop_na()

num_permutations <- 5000

process_alignment_s6_age <-  function (
    fasta_msa = "./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s6_n33_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s6_n33_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_age)  
  # create a vector of the group matching the order of the distance matrix
  ages <- acc_loc %>% pull(age)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s6_age(msa_dist, ages)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(ages, replace=F)
    perm_snn             <- calculate_snn_s6_age(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- 1/num_permutations
    this_snn_p_val <- (this_snn_p_val_sci)
  }
  
  this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val_bonf)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.25, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val_bonf,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s6_age <- function(distance_matrix, age_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(age_vector[nn_i] == age_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment6_age_snn <- process_alignment_s6_age("./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta", "Segment 6 - Animal Age, n=33, Snn=")

Segment6_age_snn$plot

# END Segment 6 Snn - testing age as lurking variable
################################################################################


################################################################################
# Segment 6 Snn - testing sex as lurking variable & no NAs

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_sex <- metadata %>% select(accession = accession, sex, segment_6)
metadata_sex <- metadata_sex %>% drop_na()

num_permutations <- 5000

process_alignment_s6_sex <-  function (
    fasta_msa = "./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s6_n33_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s6_n33_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_sex)  
  # create a vector of the group matching the order of the distance matrix
  sexs <- acc_loc %>% pull(sex)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s6_sex(msa_dist, sexs)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(sexs, replace=F)
    perm_snn             <- calculate_snn_s6_sex(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s6_sex <- function(distance_matrix, sex_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(sex_vector[nn_i] == sex_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment6_sex_snn <- process_alignment_s6_sex("./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta", "Segment 6 - Animal Sex, n=33, Snn=")

Segment6_sex_snn$plot

# END Segment 6 Snn - testing sex as lurking variable
################################################################################


################################################################################
# Segment 6 Snn - testing state as lurking variable & no NAs

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_state <- metadata %>% select(accession = accession, state, segment_6)
metadata_state <- metadata_state %>% drop_na()

num_permutations <- 5000

process_alignment_s6_state <-  function (
    fasta_msa = "./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s6_n33_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s6_n33_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_state)  
  # create a vector of the group matching the order of the distance matrix
  states <- acc_loc %>% pull(state)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s6_state(msa_dist, states)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(states, replace=F)
    perm_snn             <- calculate_snn_s6_state(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s6_state <- function(distance_matrix, state_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(state_vector[nn_i] == state_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment6_state_snn <- process_alignment_s6_state("./Alignments/A3_s6_n33_CDS_Only_aln_repeat_2.fasta", "Segment 6 - State, n=33, Snn=")

Segment6_state_snn$plot

# END Segment 6 Snn - testing state as lurking variable
################################################################################


################################################################################


# END SEGMENT 6


################################################################################



















################################################################################


# SEGMENT 7


################################################################################


################################################################################
# Segment 7 Snn - testing species as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_species <- metadata %>% select(accession = accession, species, segment_7)
metadata_species <- metadata_species %>% drop_na()

num_permutations <- 5000

process_alignment_s7_species <-  function (
    fasta_msa = "./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s7_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s7_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_species)  
  # create a vector of the group matching the order of the distance matrix
  species <- acc_loc %>% pull(species)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s7_species(msa_dist, species)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(species, replace=F)
    perm_snn             <- calculate_snn_s7_species(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=100, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s7_species <- function(distance_matrix, species_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(species_vector[nn_i] == species_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment7_species_snn <- process_alignment_s7_species("./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta", "Segment 7 - Species, n=29, Snn=")

Segment7_species_snn$plot

# END Segment 7 Snn - testing species as lurking variable
################################################################################


################################################################################
# Segment 7 Snn - testing serotype as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_serotype <- metadata %>% select(accession = accession, serotype, segment_7)
metadata_serotype <- metadata_serotype %>% drop_na()

num_permutations <- 5000

process_alignment_s7_serotype <-  function (
    fasta_msa = "./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s7_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s7_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_serotype)  
  # create a vector of the group matching the order of the distance matrix
  serotypes <- acc_loc %>% pull(serotype)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s7_serotype(msa_dist, serotypes)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(serotypes, replace=F)
    perm_snn             <- calculate_snn_s7_serotype(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.3, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s7_serotype <- function(distance_matrix, serotype_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(serotype_vector[nn_i] == serotype_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment7_serotype_snn <- process_alignment_s7_serotype("./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta", "Segment 7 - Serotype, n=29, Snn=")

Segment7_serotype_snn$plot

# END Segment 7 Snn - testing serotype as lurking variable
################################################################################


################################################################################
# Segment 7 Snn - testing year as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_year <- metadata %>% select(accession = accession, year, segment_7)
metadata_year <- metadata_year %>% drop_na()

num_permutations <- 5000

process_alignment_s7_year <-  function (
    fasta_msa = "./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s7_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s7_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_year)  
  # create a vector of the group matching the order of the distance matrix
  years <- acc_loc %>% pull(year)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s7_year(msa_dist, years)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(years, replace=F)
    perm_snn             <- calculate_snn_s7_year(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s7_year <- function(distance_matrix, year_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(year_vector[nn_i] == year_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment7_year_snn <- process_alignment_s7_year("./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta", "Segment 7 - Year, n=29, Snn=")

Segment7_year_snn$plot

# END Segment 7 Snn - testing year as lurking variable
################################################################################


################################################################################
# Segment 7 Snn - testing age as lurking variable
# There is not enough age representation across groups to reliably calculate Snn (TD 10.03.25)

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need
metadata_age <- metadata %>% select(accession = accession, age, segment_7)
metadata_age <- metadata_age %>% drop_na()

num_permutations <- 5000

process_alignment_s7_age <-  function (
    fasta_msa = "./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s7_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s7_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_age)  
  # create a vector of the group matching the order of the distance matrix
  ages <- acc_loc %>% pull(age)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s7_age(msa_dist, ages)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(ages, replace=F)
    perm_snn             <- calculate_snn_s7_age(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- 1/num_permutations
    this_snn_p_val <- (this_snn_p_val_sci)
  }
  
  this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val_bonf)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.25, y=100, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val_bonf,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s7_age <- function(distance_matrix, age_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(age_vector[nn_i] == age_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment7_age_snn <- process_alignment_s7_age("./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta", "Segment 7 - Animal Age, n=29, Snn=")

Segment7_age_snn$plot

# END Segment 7 Snn - testing age as lurking variable
################################################################################


################################################################################
# Segment 7 Snn - testing sex as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_sex <- metadata %>% select(accession = accession, sex, segment_7)
metadata_sex <- metadata_sex %>% drop_na()

num_permutations <- 5000

process_alignment_s7_sex <-  function (
    fasta_msa = "./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s7_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s7_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_sex)  
  # create a vector of the group matching the order of the distance matrix
  sexs <- acc_loc %>% pull(sex)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s7_sex(msa_dist, sexs)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(sexs, replace=F)
    perm_snn             <- calculate_snn_s7_sex(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s7_sex <- function(distance_matrix, sex_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(sex_vector[nn_i] == sex_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment7_sex_snn <- process_alignment_s7_sex("./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta", "Segment 7 - Animal Sex, n=29, Snn=")

Segment7_sex_snn$plot

# END Segment 7 Snn - testing sex as lurking variable
################################################################################

################################################################################
# Segment 7 Snn - testing state as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_state <- metadata %>% select(accession = accession, state, segment_7)
metadata_state <- metadata_state %>% drop_na()

num_permutations <- 5000

process_alignment_s7_state <-  function (
    fasta_msa = "./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s7_n29_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s7_n29_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_state)  
  # create a vector of the group matching the order of the distance matrix
  states <- acc_loc %>% pull(state)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s7_state(msa_dist, states)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(states, replace=F)
    perm_snn             <- calculate_snn_s7_state(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s7_state <- function(distance_matrix, state_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(state_vector[nn_i] == state_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment7_state_snn <- process_alignment_s7_state("./Alignments/A3_s7_n29_CDS_Only_aln_repeat_2.fasta", "Segment 7 - State, n=29, Snn=")

Segment7_state_snn$plot

# END Segment 7 Snn - testing state as lurking variable
################################################################################



################################################################################


# END SEGMENT 7


################################################################################



















################################################################################


# SEGMENT 8


################################################################################


################################################################################
# Segment 8 Snn - testing species as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_species <- metadata %>% select(accession = accession, species, segment_8)
metadata_species <- metadata_species %>% drop_na()

num_permutations <- 5000

process_alignment_s8_species <-  function (
    fasta_msa = "./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s8_n33_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s8_n33_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_species)  
  # create a vector of the group matching the order of the distance matrix
  species <- acc_loc %>% pull(species)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s8_species(msa_dist, species)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(species, replace=F)
    perm_snn             <- calculate_snn_s8_species(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s8_species <- function(distance_matrix, species_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(species_vector[nn_i] == species_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment8_species_snn <- process_alignment_s8_species("./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta", "Segment 8 - Species, n=33, Snn=")

Segment8_species_snn$plot

# END Segment 8 Snn - testing species as lurking variable
################################################################################


################################################################################
# Segment 8 Snn - testing serotype as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_serotype <- metadata %>% select(accession = accession, serotype, segment_8)
metadata_serotype <- metadata_serotype %>% drop_na()

num_permutations <- 5000

process_alignment_s8_serotype <-  function (
    fasta_msa = "./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s8_n33_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s8_n33_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_serotype)  
  # create a vector of the group matching the order of the distance matrix
  serotypes <- acc_loc %>% pull(serotype)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s8_serotype(msa_dist, serotypes)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(serotypes, replace=F)
    perm_snn             <- calculate_snn_s8_serotype(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.3, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s8_serotype <- function(distance_matrix, serotype_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(serotype_vector[nn_i] == serotype_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment8_serotype_snn <- process_alignment_s8_serotype("./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta", "Segment 8 - Serotype, n=33, Snn=")

Segment8_serotype_snn$plot

# END Segment 8 Snn - testing serotype as lurking variable
################################################################################


################################################################################
# Segment 8 Snn - testing year as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_year <- metadata %>% select(accession = accession, year, segment_8)
metadata_year <- metadata_year %>% drop_na()

num_permutations <- 5000

process_alignment_s8_year <-  function (
    fasta_msa = "./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s8_n33_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s8_n33_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_year)  
  # create a vector of the group matching the order of the distance matrix
  years <- acc_loc %>% pull(year)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s8_year(msa_dist, years)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(years, replace=F)
    perm_snn             <- calculate_snn_s8_year(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s8_year <- function(distance_matrix, year_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(year_vector[nn_i] == year_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment8_year_snn <- process_alignment_s8_year("./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta", "Segment 8 - Year, n=33, Snn=")

Segment8_year_snn$plot

# END Segment 8 Snn - testing year as lurking variable
################################################################################


################################################################################
# Segment 8 Snn - testing age as lurking variable
# There is not enough age representation across groups to reliably calculate Snn (TD 10.03.25)

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need
metadata_age <- metadata %>% select(accession = accession, age, segment_8)
metadata_age <- metadata_age %>% drop_na()

num_permutations <- 5000

process_alignment_s8_age <-  function (
    fasta_msa = "./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s8_n33_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s8_n33_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_age)  
  # create a vector of the group matching the order of the distance matrix
  ages <- acc_loc %>% pull(age)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s8_age(msa_dist, ages)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(ages, replace=F)
    perm_snn             <- calculate_snn_s8_age(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- 1/num_permutations
    this_snn_p_val <- (this_snn_p_val_sci)
  }
  
  this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val_bonf)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.25, y=100, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val_bonf,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s8_age <- function(distance_matrix, age_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(age_vector[nn_i] == age_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment8_age_snn <- process_alignment_s8_age("./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta", "Segment 8 - Animal Age, n=33, Snn=")

Segment8_age_snn$plot

# END Segment 8 Snn - testing age as lurking variable
################################################################################


################################################################################
# Segment 8 Snn - testing sex as lurking variable & no NAs

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_sex <- metadata %>% select(accession = accession, sex, segment_8)
metadata_sex <- metadata_sex %>% drop_na()

num_permutations <- 5000

process_alignment_s8_sex <-  function (
    fasta_msa = "./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s8_n33_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s8_n33_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_sex)  
  # create a vector of the group matching the order of the distance matrix
  sexs <- acc_loc %>% pull(sex)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s8_sex(msa_dist, sexs)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(sexs, replace=F)
    perm_snn             <- calculate_snn_s8_sex(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s8_sex <- function(distance_matrix, sex_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(sex_vector[nn_i] == sex_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment8_sex_snn <- process_alignment_s8_sex("./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta", "Segment 8 - Animal Sex, n=33, Snn=")

Segment8_sex_snn$plot

# END Segment 8 Snn - testing sex as lurking variable
################################################################################

################################################################################
# Segment 8 Snn - testing state as lurking variable & no NAs

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_state <- metadata %>% select(accession = accession, state, segment_8)
metadata_state <- metadata_state %>% drop_na()

num_permutations <- 5000

process_alignment_s8_state <-  function (
    fasta_msa = "./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s8_n33_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s8_n33_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_state)  
  # create a vector of the group matching the order of the distance matrix
  states <- acc_loc %>% pull(state)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s8_state(msa_dist, states)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(states, replace=F)
    perm_snn             <- calculate_snn_s8_state(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s8_state <- function(distance_matrix, state_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(state_vector[nn_i] == state_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment8_state_snn <- process_alignment_s8_state("./Alignments/A3_s8_n33_CDS_Only_aln_repeat_2.fasta", "Segment 8 - State, n=33, Snn=")

Segment8_state_snn$plot

# END Segment 8 Snn - testing state as lurking variable
################################################################################


################################################################################


# END SEGMENT 8


################################################################################



















################################################################################


# SEGMENT 9


################################################################################


################################################################################
# Segment 9 Snn - testing species as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_species <- metadata %>% select(accession = accession, species, segment_9)
metadata_species <- metadata_species %>% drop_na()

num_permutations <- 5000

process_alignment_s9_species <-  function (
    fasta_msa = "./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s9_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s9_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_species)  
  # create a vector of the group matching the order of the distance matrix
  species <- acc_loc %>% pull(species)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s9_species(msa_dist, species)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(species, replace=F)
    perm_snn             <- calculate_snn_s9_species(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s9_species <- function(distance_matrix, species_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(species_vector[nn_i] == species_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment9_species_snn <- process_alignment_s9_species("./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta", "Segment 9 - Species, n=31, Snn=")

Segment9_species_snn$plot

# END Segment 9 Snn - testing species as lurking variable
################################################################################


################################################################################
# Segment 9 Snn - testing serotype as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_serotype <- metadata %>% select(accession = accession, serotype, segment_9)
metadata_serotype <- metadata_serotype %>% drop_na()

num_permutations <- 5000

process_alignment_s9_serotype <-  function (
    fasta_msa = "./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s9_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s9_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_serotype)  
  # create a vector of the group matching the order of the distance matrix
  serotypes <- acc_loc %>% pull(serotype)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s9_serotype(msa_dist, serotypes)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(serotypes, replace=F)
    perm_snn             <- calculate_snn_s9_serotype(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.3, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s9_serotype <- function(distance_matrix, serotype_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(serotype_vector[nn_i] == serotype_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment9_serotype_snn <- process_alignment_s9_serotype("./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta", "Segment 9 - Serotype, n=31, Snn=")

Segment9_serotype_snn$plot

# END Segment 9 Snn - testing serotype as lurking variable
################################################################################


################################################################################
# Segment 9 Snn - testing year as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_year <- metadata %>% select(accession = accession, year, segment_9)
metadata_year <- metadata_year %>% drop_na()

num_permutations <- 5000

process_alignment_s9_year <-  function (
    fasta_msa = "./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s9_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s9_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_year)  
  # create a vector of the group matching the order of the distance matrix
  years <- acc_loc %>% pull(year)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s9_year(msa_dist, years)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(years, replace=F)
    perm_snn             <- calculate_snn_s9_year(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- 1/num_permutations
    this_snn_p_val <- (this_snn_p_val_sci)
  }
  
  this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val_bonf)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s9_year <- function(distance_matrix, year_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(year_vector[nn_i] == year_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment9_year_snn <- process_alignment_s9_year("./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta", "Segment 9 - Year, n=31, Snn=")

Segment9_year_snn$plot

# END Segment 9 Snn - testing year as lurking variable
################################################################################


################################################################################
# Segment 9 Snn - testing age as lurking variable
# There is not enough age representation across groups to reliably calculate Snn (TD 10.03.25)

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need
metadata_age <- metadata %>% select(accession = accession, age, segment_9)
metadata_age <- metadata_age %>% drop_na()

num_permutations <- 5000

process_alignment_s9_age <-  function (
    fasta_msa = "./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s9_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s9_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_age)  
  # create a vector of the group matching the order of the distance matrix
  ages <- acc_loc %>% pull(age)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s9_age(msa_dist, ages)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(ages, replace=F)
    perm_snn             <- calculate_snn_s9_age(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- 1/num_permutations
    this_snn_p_val <- (this_snn_p_val_sci)
  }
  
  this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val_bonf)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.25, y=100, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val_bonf,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s9_age <- function(distance_matrix, age_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(age_vector[nn_i] == age_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment9_age_snn <- process_alignment_s9_age("./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta", "Segment 9 - Animal Age, n=31, Snn=")

Segment9_age_snn$plot

# END Segment 9 Snn - testing age as lurking variable
################################################################################


################################################################################
# Segment 9 Snn - testing sex as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_sex <- metadata %>% select(accession = accession, sex, segment_9)
metadata_sex <- metadata_sex %>% drop_na()

num_permutations <- 5000

process_alignment_s9_sex <-  function (
    fasta_msa = "./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s9_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s9_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_sex)  
  # create a vector of the group matching the order of the distance matrix
  sexs <- acc_loc %>% pull(sex)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s9_sex(msa_dist, sexs)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(sexs, replace=F)
    perm_snn             <- calculate_snn_s9_sex(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s9_sex <- function(distance_matrix, sex_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(sex_vector[nn_i] == sex_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment9_sex_snn <- process_alignment_s9_sex("./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta", "Segment 9 - Animal Sex, n=31, Snn=")

Segment9_sex_snn$plot

# END Segment 9 Snn - testing sex as lurking variable
################################################################################

################################################################################
# Segment 9 Snn - testing state as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_state <- metadata %>% select(accession = accession, state, segment_9)
metadata_state <- metadata_state %>% drop_na()

num_permutations <- 5000

process_alignment_s9_state <-  function (
    fasta_msa = "./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s9_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s9_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_state)  
  # create a vector of the group matching the order of the distance matrix
  states <- acc_loc %>% pull(state)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s9_state(msa_dist, states)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(states, replace=F)
    perm_snn             <- calculate_snn_s9_state(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s9_state <- function(distance_matrix, state_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(state_vector[nn_i] == state_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment9_state_snn <- process_alignment_s9_state("./Alignments/A3_s9_n31_CDS_Only_aln_repeat_2.fasta", "Segment 9 - State, n=31, Snn=")

Segment9_state_snn$plot

# END Segment 9 Snn - testing state as lurking variable
################################################################################


################################################################################


# END SEGMENT 9


################################################################################



















################################################################################


# SEGMENT 10


################################################################################


################################################################################
# Segment 10 Snn - testing species as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_species <- metadata %>% select(accession = accession, species, segment_10)
metadata_species <- metadata_species %>% drop_na()

num_permutations <- 5000

process_alignment_s10_species <-  function (
    fasta_msa = "./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s10_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s10_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_species)  
  # create a vector of the group matching the order of the distance matrix
  species <- acc_loc %>% pull(species)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s10_species(msa_dist, species)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(species, replace=F)
    perm_snn             <- calculate_snn_s10_species(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=100, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s10_species <- function(distance_matrix, species_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(species_vector[nn_i] == species_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment10_species_snn <- process_alignment_s10_species("./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta", "Segment 10 - Species, n=31, Snn=")

Segment10_species_snn$plot

# END Segment 10 Snn - testing species as lurking variable
################################################################################


################################################################################
# Segment 10 Snn - testing serotype as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_serotype <- metadata %>% select(accession = accession, serotype, segment_10)
metadata_serotype <- metadata_serotype %>% drop_na()

num_permutations <- 5000

process_alignment_s10_serotype <-  function (
    fasta_msa = "./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s10_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s10_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_serotype)  
  # create a vector of the group matching the order of the distance matrix
  serotypes <- acc_loc %>% pull(serotype)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s10_serotype(msa_dist, serotypes)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(serotypes, replace=F)
    perm_snn             <- calculate_snn_s10_serotype(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.3, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s10_serotype <- function(distance_matrix, serotype_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(serotype_vector[nn_i] == serotype_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment10_serotype_snn <- process_alignment_s10_serotype("./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta", "Segment 10 - Serotype, n=31, Snn=")

Segment10_serotype_snn$plot

# END Segment 10 Snn - testing serotype as lurking variable
################################################################################


################################################################################
# Segment 10 Snn - testing year as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_year <- metadata %>% select(accession = accession, year, segment_10)
metadata_year <- metadata_year %>% drop_na()

num_permutations <- 5000

process_alignment_s10_year <-  function (
    fasta_msa = "./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s10_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s10_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_year)  
  # create a vector of the group matching the order of the distance matrix
  years <- acc_loc %>% pull(year)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s10_year(msa_dist, years)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(years, replace=F)
    perm_snn             <- calculate_snn_s10_year(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s10_year <- function(distance_matrix, year_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(year_vector[nn_i] == year_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment10_year_snn <- process_alignment_s10_year("./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta", "Segment 10 - Year, n=31, Snn=")

Segment10_year_snn$plot

# END Segment 10 Snn - testing year as lurking variable
################################################################################


################################################################################
# Segment 10 Snn - testing age as lurking variable
# There is not enough age representation across groups to reliably calculate Snn (TD 10.03.25)

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need
metadata_age <- metadata %>% select(accession = accession, age, segment_10)
metadata_age <- metadata_age %>% drop_na()

num_permutations <- 5000

process_alignment_s10_age <-  function (
    fasta_msa = "./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s10_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s10_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_age)  
  # create a vector of the group matching the order of the distance matrix
  ages <- acc_loc %>% pull(age)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s10_age(msa_dist, ages)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(ages, replace=F)
    perm_snn             <- calculate_snn_s10_age(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- 1/num_permutations
    this_snn_p_val <- (this_snn_p_val_sci)
  }
  
  this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val_bonf)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.25, y=100, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val_bonf,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s10_age <- function(distance_matrix, age_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(age_vector[nn_i] == age_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment10_age_snn <- process_alignment_s10_age("./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta", "Segment 10 - Animal Age, n=31, Snn=")

Segment10_age_snn$plot

# END Segment 10 Snn - testing age as lurking variable
################################################################################


################################################################################
# Segment 10 Snn - testing sex as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_sex <- metadata %>% select(accession = accession, sex, segment_10)
metadata_sex <- metadata_sex %>% drop_na()


num_permutations <- 5000

process_alignment_s10_sex <-  function (
    fasta_msa = "./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s10_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s10_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_sex)  
  # create a vector of the group matching the order of the distance matrix
  sexs <- acc_loc %>% pull(sex)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s10_sex(msa_dist, sexs)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(sexs, replace=F)
    perm_snn             <- calculate_snn_s10_sex(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s10_sex <- function(distance_matrix, sex_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(sex_vector[nn_i] == sex_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment10_sex_snn <- process_alignment_s10_sex("./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta", "Segment 10 - Animal Sex, n=31, Snn=")

Segment10_sex_snn$plot

# END Segment 10 Snn - testing sex as lurking variable
################################################################################


################################################################################
# Segment 10 Snn - testing state as lurking variable

# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
# metadata <- read_excel("A3_Metadata.xlsx")
metadata <- read_excel("../../../A3_Metadata_All.xlsx")

# pull out just the columns we need and remove any NA variables
metadata_state <- metadata %>% select(accession = accession, state, segment_10)
metadata_state <- metadata_state %>% drop_na()


num_permutations <- 5000

process_alignment_s10_state <-  function (
    fasta_msa = "./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta",
    prefix="A3_s10_n31_CDS_only_aln") {
  
  # DEBUG
  # fasta_msa = "./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta"
  # prefix="A3_s10_n31_CDS_only_aln"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata_state)  
  # create a vector of the group matching the order of the distance matrix
  states <- acc_loc %>% pull(state)
  
  # calculate Snn from distance matrix and sample metadata_species
  this_snn <- calculate_snn_s10_state(msa_dist, states)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(states, replace=F)
    perm_snn             <- calculate_snn_s10_state(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  # DO NOT adjust p-value for multiple comparisons (n=10) with Bonferroni correction
  num_tests <- 10
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  # this_snn_p_val_bonf <- min(this_snn_p_val * num_tests, 1)
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.00, y=150, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14)
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              # pvalue = this_snn_p_val_bonf,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s10_state <- function(distance_matrix, state_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(state_vector[nn_i] == state_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment10_state_snn <- process_alignment_s10_state("./Alignments/A3_s10_n31_CDS_Only_aln_repeat_2.fasta", "Segment 10 - State, n=31, Snn=")

Segment10_state_snn$plot

# END Segment 10 Snn - testing state as lurking variable
################################################################################


################################################################################


# END SEGMENT 10


################################################################################



















################################################################################


# ADJUSTING P-VALUES USING p.adjust


################################################################################


# Adjust Clinical and Subclinical p-values
##########################
# status_pval_supp_table <- 
#   tibble( segment = c("BTV Segment 1", "BTV Segment 2", "BTV Segment 3", 
#                       "BTV Segment 4", "BTV Segment 5", "BTV Segment 6", 
#                       "BTV Segment 7", "BTV Segment 8", "BTV Segment 9", "BTV Segment 10"),
#           snn = c(Segment1_snn$snn, Segment2_snn$snn, Segment3_snn$snn, 
#                   Segment4_snn$snn, Segment5_snn$snn, Segment6_snn$snn, 
#                   Segment7_snn$snn, Segment8_snn$snn, Segment9_snn$snn, Segment10_snn$snn)
#           pval = c(Segment1_snn$pval, Segment2_snn$pval, Segment3_snn$pval, 
#                    Segment4_snn$pval, Segment5_snn$pval, Segment6_snn$pval, 
#                    Segment7_snn$pval, Segment8_snn$pval, Segment9_snn$pval, Segment10_snn$pval)
#   ) %>%
#   mutate(pval_bonf = p.adjust(pval, method = "bonferroni", n = length(pval)))
# 
# write.table(status_pval_supp_table, 
#             file="all_segments_all_samples_status_supplemental_table_snn_031026.txt",
#             quote=F, sep="\t", row.names=F)

##########################



# Adjust Species p-values
##########################
species_pval_supp_table <- 
  tibble( segment = c("BTV Segment 1", "BTV Segment 2", "BTV Segment 3", 
                      "BTV Segment 4", "BTV Segment 5", "BTV Segment 6", 
                      "BTV Segment 7", "BTV Segment 8", "BTV Segment 9", "BTV Segment 10"),
          n = c("s1 = n29", "s2 = n29", "s3 = n30",
                "s4 = n31", "s5 = n31", "s6 = n33", 
                "s7 = n29", "s8 = n33", "s9 = n31", "s10 = n31"),
          snn = c(Segment1_species_snn$snn, Segment2_species_snn$snn, Segment3_species_snn$snn, 
                  Segment4_species_snn$snn, Segment5_species_snn$snn, Segment6_species_snn$snn, 
                  Segment7_species_snn$snn, Segment8_species_snn$snn, Segment9_species_snn$snn, 
                  Segment10_species_snn$snn),
          pval = c(Segment1_species_snn$pval, Segment2_species_snn$pval, Segment3_species_snn$pval, 
                   Segment4_species_snn$pval, Segment5_species_snn$pval, Segment6_species_snn$pval, 
                   Segment7_species_snn$pval, Segment8_species_snn$pval, Segment9_species_snn$pval, 
                   Segment10_species_snn$pval)
  ) %>%
  mutate(pval_bonf = p.adjust(pval, method = "bonferroni", n = length(pval)))

write.table(species_pval_supp_table, 
            file="all_segments_all_samples_SPECIES_supplemental_table_snn_061826.txt",
            quote=F, sep="\t", row.names=F)

##########################



# Adjust Serotype p-values
##########################
serotype_pval_supp_table <- 
  tibble( segment = c("BTV Segment 1", "BTV Segment 2", "BTV Segment 3", 
                      "BTV Segment 4", "BTV Segment 5", "BTV Segment 6", 
                      "BTV Segment 7", "BTV Segment 8", "BTV Segment 9", "BTV Segment 10"),
          n = c("s1 = n29", "s2 = n29", "s3 = n30",
                "s4 = n31", "s5 = n31", "s6 = n33", 
                "s7 = n29", "s8 = n33", "s9 = n31", "s10 = n31"),
          snn = c(Segment1_serotype_snn$snn, Segment2_serotype_snn$snn, Segment3_serotype_snn$snn, 
                  Segment4_serotype_snn$snn, Segment5_serotype_snn$snn, Segment6_serotype_snn$snn, 
                  Segment7_serotype_snn$snn, Segment8_serotype_snn$snn, Segment9_serotype_snn$snn, 
                  Segment10_serotype_snn$snn),
          pval = c(Segment1_serotype_snn$pval, Segment2_serotype_snn$pval, Segment3_serotype_snn$pval, 
                   Segment4_serotype_snn$pval, Segment5_serotype_snn$pval, Segment6_serotype_snn$pval, 
                   Segment7_serotype_snn$pval, Segment8_serotype_snn$pval, Segment9_serotype_snn$pval, 
                   Segment10_serotype_snn$pval)
  ) %>%
  mutate(pval_bonf = p.adjust(pval, method = "bonferroni", n = length(pval)))

write.table(serotype_pval_supp_table, 
            file="all_segments_all_samples_SEROTYPE_supplemental_table_snn_061826.txt",
            quote=F, sep="\t", row.names=F)

##########################



# Adjust Year p-values
##########################
year_pval_supp_table <- 
  tibble( segment = c("BTV Segment 1", "BTV Segment 2", "BTV Segment 3", 
                      "BTV Segment 4", "BTV Segment 5", "BTV Segment 6", 
                      "BTV Segment 7", "BTV Segment 8", "BTV Segment 9", "BTV Segment 10"),
          n = c("s1 = n29", "s2 = n29", "s3 = n30",
                "s4 = n31", "s5 = n31", "s6 = n33", 
                "s7 = n29", "s8 = n33", "s9 = n31", "s10 = n31"),
          snn = c(Segment1_year_snn$snn, Segment2_year_snn$snn, Segment3_year_snn$snn, 
                  Segment4_year_snn$snn, Segment5_year_snn$snn, Segment6_year_snn$snn, 
                  Segment7_year_snn$snn, Segment8_year_snn$snn, Segment9_year_snn$snn, 
                  Segment10_year_snn$snn),
          pval = c(Segment1_year_snn$pval, Segment2_year_snn$pval, Segment3_year_snn$pval, 
                   Segment4_year_snn$pval, Segment5_year_snn$pval, Segment6_year_snn$pval, 
                   Segment7_year_snn$pval, Segment8_year_snn$pval, Segment9_year_snn$pval, 
                   Segment10_year_snn$pval)
  ) %>%
  mutate(pval_bonf = p.adjust(pval, method = "bonferroni", n = length(pval)))

write.table(year_pval_supp_table, 
            file="all_segments_all_samples_YEAR_supplemental_table_snn_061826.txt",
            quote=F, sep="\t", row.names=F)

##########################



# Adjust Sex p-values
##########################
sex_pval_supp_table <- 
  tibble( segment = c("BTV Segment 1", "BTV Segment 2", "BTV Segment 3", 
                      "BTV Segment 4", "BTV Segment 5", "BTV Segment 6", 
                      "BTV Segment 7", "BTV Segment 8", "BTV Segment 9", "BTV Segment 10"),
          n = c("s1 = n29", "s2 = n29", "s3 = n30",
                "s4 = n31", "s5 = n31", "s6 = n33", 
                "s7 = n29", "s8 = n33", "s9 = n31", "s10 = n31"),
          snn = c(Segment1_sex_snn$snn, Segment2_sex_snn$snn, Segment3_sex_snn$snn, 
                  Segment4_sex_snn$snn, Segment5_sex_snn$snn, Segment6_sex_snn$snn, 
                  Segment7_sex_snn$snn, Segment8_sex_snn$snn, Segment9_sex_snn$snn, 
                  Segment10_sex_snn$snn),
          pval = c(Segment1_sex_snn$pval, Segment2_sex_snn$pval, Segment3_sex_snn$pval, 
                   Segment4_sex_snn$pval, Segment5_sex_snn$pval, Segment6_sex_snn$pval, 
                   Segment7_sex_snn$pval, Segment8_sex_snn$pval, Segment9_sex_snn$pval, 
                   Segment10_sex_snn$pval)
  ) %>%
  mutate(pval_bonf = p.adjust(pval, method = "bonferroni", n = length(pval)))

write.table(sex_pval_supp_table, 
            file="all_segments_all_samples_SEX_supplemental_table_snn_061826.txt",
            quote=F, sep="\t", row.names=F)
##########################


# Adjust Age p-values
##########################
age_pval_supp_table <- 
  tibble( segment = c("BTV Segment 1", "BTV Segment 2", "BTV Segment 3", 
                      "BTV Segment 4", "BTV Segment 5", "BTV Segment 6", 
                      "BTV Segment 7", "BTV Segment 8", "BTV Segment 9", "BTV Segment 10"),
          n = c("s1 = n29", "s2 = n29", "s3 = n30",
                "s4 = n31", "s5 = n31", "s6 = n33", 
                "s7 = n29", "s8 = n33", "s9 = n31", "s10 = n31"),
          snn = c(Segment1_age_snn$snn, Segment2_age_snn$snn, Segment3_age_snn$snn, 
                  Segment4_age_snn$snn, Segment5_age_snn$snn, Segment6_age_snn$snn, 
                  Segment7_age_snn$snn, Segment8_age_snn$snn, Segment9_age_snn$snn, 
                  Segment10_age_snn$snn),
          pval = c(Segment1_age_snn$pval, Segment2_age_snn$pval, Segment3_age_snn$pval, 
                   Segment4_age_snn$pval, Segment5_age_snn$pval, Segment6_age_snn$pval, 
                   Segment7_age_snn$pval, Segment8_age_snn$pval, Segment9_age_snn$pval, 
                   Segment10_age_snn$pval)
  ) %>%
  mutate(pval_bonf = p.adjust(pval, method = "bonferroni", n = length(pval)))

write.table(age_pval_supp_table, 
            file="all_segments_all_samples_AGE_supplemental_table_snn_061826.txt",
            quote=F, sep="\t", row.names=F)
##########################


# Adjust State p-values
##########################
state_pval_supp_table <- 
  tibble( segment = c("BTV Segment 1", "BTV Segment 2", "BTV Segment 3", 
                      "BTV Segment 4", "BTV Segment 5", "BTV Segment 6", 
                      "BTV Segment 7", "BTV Segment 8", "BTV Segment 9", "BTV Segment 10"),
          n = c("s1 = n29", "s2 = n29", "s3 = n30",
                "s4 = n31", "s5 = n31", "s6 = n33", 
                "s7 = n29", "s8 = n33", "s9 = n31", "s10 = n31"),
          snn = c(Segment1_state_snn$snn, Segment2_state_snn$snn, Segment3_state_snn$snn, 
                  Segment4_state_snn$snn, Segment5_state_snn$snn, Segment6_state_snn$snn, 
                  Segment7_state_snn$snn, Segment8_state_snn$snn, Segment9_state_snn$snn, 
                  Segment10_state_snn$snn),
          pval = c(Segment1_state_snn$pval, Segment2_state_snn$pval, Segment3_state_snn$pval, 
                   Segment4_state_snn$pval, Segment5_state_snn$pval, Segment6_state_snn$pval, 
                   Segment7_state_snn$pval, Segment8_state_snn$pval, Segment9_state_snn$pval, 
                   Segment10_state_snn$pval)
  ) %>%
  mutate(pval_bonf = p.adjust(pval, method = "bonferroni", n = length(pval)))

write.table(state_pval_supp_table, 
            file="all_segments_all_samples_STATE_supplemental_table_snn_061826.txt",
            quote=F, sep="\t", row.names=F)
##########################

################################################################################



