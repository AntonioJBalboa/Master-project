library(Biostrings)
library(pwalign)
library(BiocGenerics)
library(tidyverse)
library(vegan)
library(dplyr)
library(purrr)
library(ggupset)
library(tibble)
library(ggplot2)
library(tidyr)

# k

evenness_k <- function(OTU_table){
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  shannon_index <- diversity(asv, index = "shannon") #Defining shannon index
  S <- specnumber(asv) #Defining Richness
  J <- ifelse(S > 1, shannon_index / log(S), 0) #Calculating J (con corrección para S=1)
  
  return(J)
}


bray_curtis <- function(OTU_table){
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric)) %>%
    filter(rowSums(.) > 0)  # Delete empty rows to remove errors
  bray <- vegdist(asv, method = "bray")
  return(as.dist(bray))
}


calculate_k_w_single <- function(OTU_table, alpha, beta) {
  
  asv_num <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Check if there is enough data to calculate
  if (ncol(asv_num) < 2 || nrow(asv_num) < 2) {
    return(0)
  }
  
  J_mean <- mean(na.omit(evenness_k(asv_num)))
  BC_mean <- mean(na.omit(bray_curtis(asv_num)))
  
  if (is.nan(J_mean) || is.nan(BC_mean)) {
    return(0)
  }
  
  k_w <- (J_mean^alpha) * (BC_mean^beta)
  return(k_w)
}

# Weights as factor
calculate_k_w_single2 <- function(OTU_table, alpha, beta) {
  
  asv_num <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Check if there is enough data to calculate
  if (ncol(asv_num) < 2 || nrow(asv_num) < 2) {
    return(0)
  }
  
  J_mean <- mean(na.omit(evenness_k(asv_num)))
  BC_mean <- mean(na.omit(bray_curtis(asv_num)))
  
  if (is.nan(J_mean) || is.nan(BC_mean)) {
    return(0)
  }
  
  k_w2 <- (J_mean*alpha) * (BC_mean*beta)
  return(k_w2)
}


# Conservative filter standard
conservative_filter <- function(OTU_table, col_max = 10,pct_samples_present = 0.05,
                                abundance_threshold = 0.00005, prevalence = 2) {
  
  num_samples <- nrow(OTU_table)
  total_reads_dataset <- sum(OTU_table)
  
  if (num_samples == 0 || total_reads_dataset == 0) {
   
    # Handling empty tables
    return(list(preserved = OTU_table[0, 0, drop = FALSE], 
                deleted = OTU_table[0, 0, drop = FALSE]))
  }
  
  col_totals <- colSums(OTU_table)
  col_maxs <- apply(OTU_table, 2, max) # The maximum per ASV
  samples_present <- apply(OTU_table, 2, function(c) sum(c > 0)) 
  
  # Rule 1 
  pct_samples_present <- samples_present / num_samples
  C1_delete <- (col_maxs <= col_max) & (pct_samples_present <= pct_samples_present)
  
  # Rule 2
  sample_sums <- rowSums(OTU_table)
  rel_abun_matrix <- sweep(OTU_table, 1, sample_sums, "/")
  max_rel_abun_in_any_sample <- apply(rel_abun_matrix, 2, max)
  C2_delete <- max_rel_abun_in_any_sample < abundance_threshold
  
  # Rule 3
  C3_delete <- samples_present < prevalence
  
  # Combine
  asv_to_delete_bool <-   C1_delete | C2_delete | C3_delete
  
  asv_ids_to_preserve <- colnames(OTU_table)[!asv_to_delete_bool]
  asv_ids_to_delete <- colnames(OTU_table)[asv_to_delete_bool]
  
  # Two lists. One with the preserved ASVs and another with the filtered ones
  return(list(
    preserved = OTU_table[, asv_ids_to_preserve, drop = FALSE],
    deleted = OTU_table[, asv_ids_to_delete, drop = FALSE]
  ))
}

# Conservative filter 2 (0.00007)
conservative_filter2 <- function(OTU_table, col_max = 10,pct_samples_present = 0.05,
                                abundance_threshold = 0.00007, prevalence = 2) {
  
  num_samples <- nrow(OTU_table)
  total_reads_dataset <- sum(OTU_table)
  
  if (num_samples == 0 || total_reads_dataset == 0) {
    
    # Handling empty tables
    return(list(preserved = OTU_table[0, 0, drop = FALSE], 
                deleted = OTU_table[0, 0, drop = FALSE]))
  }
  
  col_totals <- colSums(OTU_table)
  col_maxs <- apply(OTU_table, 2, max) # The maximum per ASV
  samples_present <- apply(OTU_table, 2, function(c) sum(c > 0)) 
  
  # Rule 1 
  pct_samples_present <- samples_present / num_samples
  C1_delete <- (col_maxs <= col_max) & (pct_samples_present <= pct_samples_present)
  
  # Rule 2
  sample_sums <- rowSums(OTU_table)
  rel_abun_matrix <- sweep(OTU_table, 1, sample_sums, "/")
  max_rel_abun_in_any_sample <- apply(rel_abun_matrix, 2, max)
  C2_delete <- max_rel_abun_in_any_sample < abundance_threshold
  
  # Rule 3
  C3_delete <- samples_present < prevalence
  
  # Combine
  asv_to_delete_bool <-   C1_delete | C2_delete | C3_delete
  
  asv_ids_to_preserve <- colnames(OTU_table)[!asv_to_delete_bool]
  asv_ids_to_delete <- colnames(OTU_table)[asv_to_delete_bool]
  
  # Two lists. One with the preserved ASVs and another with the filtered ones
  return(list(
    preserved = OTU_table[, asv_ids_to_preserve, drop = FALSE],
    deleted = OTU_table[, asv_ids_to_delete, drop = FALSE]
  ))
}


# Conservative filter 3 (0.0001)
conservative_filter3 <- function(OTU_table, col_max = 10,pct_samples_present = 0.05,
                                abundance_threshold = 0.0001, prevalence = 2) {
  
  num_samples <- nrow(OTU_table)
  total_reads_dataset <- sum(OTU_table)
  
  if (num_samples == 0 || total_reads_dataset == 0) {
    
    # Handling empty tables
    return(list(preserved = OTU_table[0, 0, drop = FALSE], 
                deleted = OTU_table[0, 0, drop = FALSE]))
  }
  
  col_totals <- colSums(OTU_table)
  col_maxs <- apply(OTU_table, 2, max) # The maximum per ASV
  samples_present <- apply(OTU_table, 2, function(c) sum(c > 0)) 
  
  # Rule 1 
  pct_samples_present <- samples_present / num_samples
  C1_delete <- (col_maxs <= col_max) & (pct_samples_present <= pct_samples_present)
  
  # Rule 2
  sample_sums <- rowSums(OTU_table)
  rel_abun_matrix <- sweep(OTU_table, 1, sample_sums, "/")
  max_rel_abun_in_any_sample <- apply(rel_abun_matrix, 2, max)
  C2_delete <- max_rel_abun_in_any_sample < abundance_threshold
  
  # Rule 3
  C3_delete <- samples_present < prevalence
  
  # Combine
  asv_to_delete_bool <-   C1_delete | C2_delete | C3_delete
  
  asv_ids_to_preserve <- colnames(OTU_table)[!asv_to_delete_bool]
  asv_ids_to_delete <- colnames(OTU_table)[asv_to_delete_bool]
  
  # Two lists. One with the preserved ASVs and another with the filtered ones
  return(list(
    preserved = OTU_table[, asv_ids_to_preserve, drop = FALSE],
    deleted = OTU_table[, asv_ids_to_delete, drop = FALSE]
  ))
}

# Conservative filter 4 (0.0002)
conservative_filter4 <- function(OTU_table, col_max = 10,pct_samples_present = 0.05,
                                 abundance_threshold = 0.0002, prevalence = 2) {
  
  num_samples <- nrow(OTU_table)
  total_reads_dataset <- sum(OTU_table)
  
  if (num_samples == 0 || total_reads_dataset == 0) {
    
    # Handling empty tables
    return(list(preserved = OTU_table[0, 0, drop = FALSE], 
                deleted = OTU_table[0, 0, drop = FALSE]))
  }
  
  col_totals <- colSums(OTU_table)
  col_maxs <- apply(OTU_table, 2, max) # The maximum per ASV
  samples_present <- apply(OTU_table, 2, function(c) sum(c > 0)) 
  
  # Rule 1 
  pct_samples_present <- samples_present / num_samples
  C1_delete <- (col_maxs <= col_max) & (pct_samples_present <= pct_samples_present)
  
  # Rule 2
  sample_sums <- rowSums(OTU_table)
  rel_abun_matrix <- sweep(OTU_table, 1, sample_sums, "/")
  max_rel_abun_in_any_sample <- apply(rel_abun_matrix, 2, max)
  C2_delete <- max_rel_abun_in_any_sample < abundance_threshold
  
  # Rule 3
  C3_delete <- samples_present < prevalence
  
  # Combine
  asv_to_delete_bool <-   C1_delete | C2_delete | C3_delete
  
  asv_ids_to_preserve <- colnames(OTU_table)[!asv_to_delete_bool]
  asv_ids_to_delete <- colnames(OTU_table)[asv_to_delete_bool]
  
  # Two lists. One with the preserved ASVs and another with the filtered ones
  return(list(
    preserved = OTU_table[, asv_ids_to_preserve, drop = FALSE],
    deleted = OTU_table[, asv_ids_to_delete, drop = FALSE]
  ))
}


# Conservative filter 5 (0.0003)
conservative_filter5 <- function(OTU_table, col_max = 10,pct_samples_present = 0.05,
                                 abundance_threshold = 0.0003, prevalence = 2) {
  
  num_samples <- nrow(OTU_table)
  total_reads_dataset <- sum(OTU_table)
  
  if (num_samples == 0 || total_reads_dataset == 0) {
    
    # Handling empty tables
    return(list(preserved = OTU_table[0, 0, drop = FALSE], 
                deleted = OTU_table[0, 0, drop = FALSE]))
  }
  
  col_totals <- colSums(OTU_table)
  col_maxs <- apply(OTU_table, 2, max) # The maximum per ASV
  samples_present <- apply(OTU_table, 2, function(c) sum(c > 0)) 
  
  # Rule 1 
  pct_samples_present <- samples_present / num_samples
  C1_delete <- (col_maxs <= col_max) & (pct_samples_present <= pct_samples_present)
  
  # Rule 2
  sample_sums <- rowSums(OTU_table)
  rel_abun_matrix <- sweep(OTU_table, 1, sample_sums, "/")
  max_rel_abun_in_any_sample <- apply(rel_abun_matrix, 2, max)
  C2_delete <- max_rel_abun_in_any_sample < abundance_threshold
  
  # Rule 3
  C3_delete <- samples_present < prevalence
  
  # Combine
  asv_to_delete_bool <-   C1_delete | C2_delete | C3_delete
  
  asv_ids_to_preserve <- colnames(OTU_table)[!asv_to_delete_bool]
  asv_ids_to_delete <- colnames(OTU_table)[asv_to_delete_bool]
  
  # Two lists. One with the preserved ASVs and another with the filtered ones
  return(list(
    preserved = OTU_table[, asv_ids_to_preserve, drop = FALSE],
    deleted = OTU_table[, asv_ids_to_delete, drop = FALSE]
  ))
}


# Conservative filter 6 (0.0009)
conservative_filter6 <- function(OTU_table, col_max = 10,pct_samples_present = 0.05,
                                 abundance_threshold = 0.0009, prevalence = 2) {
  
  num_samples <- nrow(OTU_table)
  total_reads_dataset <- sum(OTU_table)
  
  if (num_samples == 0 || total_reads_dataset == 0) {
    
    # Handling empty tables
    return(list(preserved = OTU_table[0, 0, drop = FALSE], 
                deleted = OTU_table[0, 0, drop = FALSE]))
  }
  
  col_totals <- colSums(OTU_table)
  col_maxs <- apply(OTU_table, 2, max) # The maximum per ASV
  samples_present <- apply(OTU_table, 2, function(c) sum(c > 0)) 
  
  # Rule 1 
  pct_samples_present <- samples_present / num_samples
  C1_delete <- (col_maxs <= col_max) & (pct_samples_present <= pct_samples_present)
  
  # Rule 2
  sample_sums <- rowSums(OTU_table)
  rel_abun_matrix <- sweep(OTU_table, 1, sample_sums, "/")
  max_rel_abun_in_any_sample <- apply(rel_abun_matrix, 2, max)
  C2_delete <- max_rel_abun_in_any_sample < abundance_threshold
  
  # Rule 3
  C3_delete <- samples_present < prevalence
  
  # Combine
  asv_to_delete_bool <-   C1_delete | C2_delete | C3_delete
  
  asv_ids_to_preserve <- colnames(OTU_table)[!asv_to_delete_bool]
  asv_ids_to_delete <- colnames(OTU_table)[asv_to_delete_bool]
  
  # Two lists. One with the preserved ASVs and another with the filtered ones
  return(list(
    preserved = OTU_table[, asv_ids_to_preserve, drop = FALSE],
    deleted = OTU_table[, asv_ids_to_delete, drop = FALSE]
  ))
}


# Conservative filter 7 (0.00035)
conservative_filter7 <- function(OTU_table, col_max = 10,pct_samples_present = 0.05,
                                 abundance_threshold = 0.00035, prevalence = 2) {
  
  num_samples <- nrow(OTU_table)
  total_reads_dataset <- sum(OTU_table)
  
  if (num_samples == 0 || total_reads_dataset == 0) {
    
    # Handling empty tables
    return(list(preserved = OTU_table[0, 0, drop = FALSE], 
                deleted = OTU_table[0, 0, drop = FALSE]))
  }
  
  col_totals <- colSums(OTU_table)
  col_maxs <- apply(OTU_table, 2, max) # The maximum per ASV
  samples_present <- apply(OTU_table, 2, function(c) sum(c > 0)) 
  
  # Rule 1 
  pct_samples_present <- samples_present / num_samples
  C1_delete <- (col_maxs <= col_max) & (pct_samples_present <= pct_samples_present)
  
  # Rule 2
  sample_sums <- rowSums(OTU_table)
  rel_abun_matrix <- sweep(OTU_table, 1, sample_sums, "/")
  max_rel_abun_in_any_sample <- apply(rel_abun_matrix, 2, max)
  C2_delete <- max_rel_abun_in_any_sample < abundance_threshold
  
  # Rule 3
  C3_delete <- samples_present < prevalence
  
  # Combine
  asv_to_delete_bool <-   C1_delete | C2_delete | C3_delete
  
  asv_ids_to_preserve <- colnames(OTU_table)[!asv_to_delete_bool]
  asv_ids_to_delete <- colnames(OTU_table)[asv_to_delete_bool]
  
  # Two lists. One with the preserved ASVs and another with the filtered ones
  return(list(
    preserved = OTU_table[, asv_ids_to_preserve, drop = FALSE],
    deleted = OTU_table[, asv_ids_to_delete, drop = FALSE]
  ))
}

# Conservative filter 8 (0.00042)
conservative_filter8 <- function(OTU_table, col_max = 10,pct_samples_present = 0.05,
                                 abundance_threshold = 0.00042, prevalence = 2) {
  
  num_samples <- nrow(OTU_table)
  total_reads_dataset <- sum(OTU_table)
  
  if (num_samples == 0 || total_reads_dataset == 0) {
    
    # Handling empty tables
    return(list(preserved = OTU_table[0, 0, drop = FALSE], 
                deleted = OTU_table[0, 0, drop = FALSE]))
  }
  
  col_totals <- colSums(OTU_table)
  col_maxs <- apply(OTU_table, 2, max) # The maximum per ASV
  samples_present <- apply(OTU_table, 2, function(c) sum(c > 0)) 
  
  # Rule 1 
  pct_samples_present <- samples_present / num_samples
  C1_delete <- (col_maxs <= col_max) & (pct_samples_present <= pct_samples_present)
  
  # Rule 2
  sample_sums <- rowSums(OTU_table)
  rel_abun_matrix <- sweep(OTU_table, 1, sample_sums, "/")
  max_rel_abun_in_any_sample <- apply(rel_abun_matrix, 2, max)
  C2_delete <- max_rel_abun_in_any_sample < abundance_threshold
  
  # Rule 3
  C3_delete <- samples_present < prevalence
  
  # Combine
  asv_to_delete_bool <-   C1_delete | C2_delete | C3_delete
  
  asv_ids_to_preserve <- colnames(OTU_table)[!asv_to_delete_bool]
  asv_ids_to_delete <- colnames(OTU_table)[asv_to_delete_bool]
  
  # Two lists. One with the preserved ASVs and another with the filtered ones
  return(list(
    preserved = OTU_table[, asv_ids_to_preserve, drop = FALSE],
    deleted = OTU_table[, asv_ids_to_delete, drop = FALSE]
  ))
}

# Conservative filter 9 (0.00045)
conservative_filter9 <- function(OTU_table, col_max = 10,pct_samples_present = 0.05,
                                 abundance_threshold = 0.00045, prevalence = 2) {
  
  num_samples <- nrow(OTU_table)
  total_reads_dataset <- sum(OTU_table)
  
  if (num_samples == 0 || total_reads_dataset == 0) {
    
    # Handling empty tables
    return(list(preserved = OTU_table[0, 0, drop = FALSE], 
                deleted = OTU_table[0, 0, drop = FALSE]))
  }
  
  col_totals <- colSums(OTU_table)
  col_maxs <- apply(OTU_table, 2, max) # The maximum per ASV
  samples_present <- apply(OTU_table, 2, function(c) sum(c > 0)) 
  
  # Rule 1 
  pct_samples_present <- samples_present / num_samples
  C1_delete <- (col_maxs <= col_max) & (pct_samples_present <= pct_samples_present)
  
  # Rule 2
  sample_sums <- rowSums(OTU_table)
  rel_abun_matrix <- sweep(OTU_table, 1, sample_sums, "/")
  max_rel_abun_in_any_sample <- apply(rel_abun_matrix, 2, max)
  C2_delete <- max_rel_abun_in_any_sample < abundance_threshold
  
  # Rule 3
  C3_delete <- samples_present < prevalence
  
  # Combine
  asv_to_delete_bool <-   C1_delete | C2_delete | C3_delete
  
  asv_ids_to_preserve <- colnames(OTU_table)[!asv_to_delete_bool]
  asv_ids_to_delete <- colnames(OTU_table)[asv_to_delete_bool]
  
  # Two lists. One with the preserved ASVs and another with the filtered ones
  return(list(
    preserved = OTU_table[, asv_ids_to_preserve, drop = FALSE],
    deleted = OTU_table[, asv_ids_to_delete, drop = FALSE]
  ))
}

# Conservative filter 10 (0.0005)
conservative_filter10 <- function(OTU_table, col_max = 10,pct_samples_present = 0.05,
                                 abundance_threshold = 0.0005, prevalence = 2) {
  
  num_samples <- nrow(OTU_table)
  total_reads_dataset <- sum(OTU_table)
  
  if (num_samples == 0 || total_reads_dataset == 0) {
    
    # Handling empty tables
    return(list(preserved = OTU_table[0, 0, drop = FALSE], 
                deleted = OTU_table[0, 0, drop = FALSE]))
  }
  
  col_totals <- colSums(OTU_table)
  col_maxs <- apply(OTU_table, 2, max) # The maximum per ASV
  samples_present <- apply(OTU_table, 2, function(c) sum(c > 0)) 
  
  # Rule 1 
  pct_samples_present <- samples_present / num_samples
  C1_delete <- (col_maxs <= col_max) & (pct_samples_present <= pct_samples_present)
  
  # Rule 2
  sample_sums <- rowSums(OTU_table)
  rel_abun_matrix <- sweep(OTU_table, 1, sample_sums, "/")
  max_rel_abun_in_any_sample <- apply(rel_abun_matrix, 2, max)
  C2_delete <- max_rel_abun_in_any_sample < abundance_threshold
  
  # Rule 3
  C3_delete <- samples_present < prevalence
  
  # Combine
  asv_to_delete_bool <-   C1_delete | C2_delete | C3_delete
  
  asv_ids_to_preserve <- colnames(OTU_table)[!asv_to_delete_bool]
  asv_ids_to_delete <- colnames(OTU_table)[asv_to_delete_bool]
  
  # Two lists. One with the preserved ASVs and another with the filtered ones
  return(list(
    preserved = OTU_table[, asv_ids_to_preserve, drop = FALSE],
    deleted = OTU_table[, asv_ids_to_delete, drop = FALSE]
  ))
}

# Optimized Add-one procedure

# Filter conservative standard pairwaise

k_validation_pairwaise <- function(OTU_table, ground_truth, alpha, beta, mode, sig_level = 1.96){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- calculate_k_w_single(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  # apply conservative filter
  filter_output <- conservative_filter(asv)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- calculate_k_w_single(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- BC_values <- J_values <- c()
  true_positive_levels <- c() 
  
  
  blast_nuc_matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- calculate_k_w_single(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    j_i <- mean(na.omit(evenness_k(asv_table_i)))
    
    # Evaluate with ground_truth depending on the chosen mode
    
    match_level <- "No Match"
    
    if (mode == "taxonomy") {
      asv_species <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
      asv_genus   <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
      
      if (asv_species %in% ground_truth$species || asv_genus %in% ground_truth$genus) {
        match_level <- "Match"
      }
      
    } else if (mode == "sequence") {
      asv_identity <- seq_dict$Sequence[match(asv_name, seq_dict$ASV)]
      
      match_matrix <- sapply(ground_truth$Sequence, function(gt_seq) {
        
        alignment <- pwalign::pairwiseAlignment(as.character(asv_identity), as.character(gt_seq), type = "local",
                                                   substitutionMatrix = blast_nuc_matrix, gapOpening =5,
                                                gapExtension = 2)
        
        if (score(alignment) <= 0) {
          return(c(pid = 0, qcov = 0)) # No alignment
        }
        
        percent_identity <- pwalign::pid(alignment, type = "PID1")
        
        # Query Coverage
        pattern_string <- as.character(alignment@pattern)
        pattern_no_gaps <- gsub("-", "", pattern_string)
        query_coverage <- (nchar(pattern_no_gaps) / nchar(asv_identity)) * 100
        
        return(c(pid = percent_identity, qcov = query_coverage))
      })
      
      
      if (is.vector(match_matrix)) {
        match_matrix <- t(as.matrix(match_matrix))
      }else {
        match_matrix <- t(match_matrix)
      }
      
      if (nrow(match_matrix) == 0) {
        best_pid <- 0; best_qcov <- 0
      } else {
        best_match_row <- match_matrix[order(match_matrix[, "qcov"],
                                             match_matrix[, "pid"],
                                             decreasing = TRUE)[1], ]
        best_pid <- best_match_row["pid"]
        best_qcov <- best_match_row["qcov"]
      }
      
      
      best_pid_rounded  <- round(best_pid, 2)
      best_qcov_rounded <- round(best_qcov, 2)
      
      # Classification
      if (best_pid >= 98 && best_qcov >= 100) {
        match_level <- "Match"
      }
    }
      
      
    added_asv <- c(added_asv, asv_name)
    k_values <- c(k_values, k_i)
    BC_values <- c(BC_values, bc_i)
    J_values <- c(J_values, j_i)
    true_positive_levels <- c(true_positive_levels, match_level)
    
  }
  
  is_tp <- true_positive_levels != "No Match"
  total_positives_discarded <- sum(is_tp)
  delta_k_all <- k_values - k_0
  
  # Calculate thresholds
  delta_sd <- sd(delta_k_all, na.rm = TRUE)
  sig_threshold_upper <- sig_level * delta_sd
  sig_threshold_lower <- -sig_level * delta_sd
  
  
  # Criteria
  crit_1 <- which(k_values >= k_0) 
  crit_2 <- which(delta_k_all >= sig_threshold_upper)
  crit_3 <- which(delta_k_all >= sig_threshold_lower) 
  
  # Approach 1: Metrics if k >= 0
  reintroduced_flags1 <- is_tp[crit_1]
  
  TP1 <- sum(reintroduced_flags1)
  FP1 <- sum(!reintroduced_flags1)
  
  precision1 <- ifelse((TP1 + FP1) > 0, TP1 / (TP1 + FP1), NA)
  
  
  
  recall1 <- ifelse(total_positives_discarded > 0, TP1 / total_positives_discarded, NA)
  F11 <- ifelse(is.na(precision1) | is.na(recall1) | (precision1 + recall1 == 0),
               NA,
               2 * (precision1 * recall1) / (precision1 + recall1))
  FDR1 <- ifelse((TP1 + FP1) > 0, FP1 / (TP1 + FP1), NA)

  
  metrics_1 <- data.frame(criterion = "k >= k_0",
                          TP = TP1, FP = FP1,
                          precision = precision1, recall = recall1,
                          FDR = FDR1, F1 = F11)
  
  
  
  # Approach 2: Metrics >= sig_upper (significant increase)
  reintroduced_flags2 <- is_tp[crit_2]
  TP2 <- sum(reintroduced_flags2)
  FP2 <- sum(!reintroduced_flags2)
  
  precision2 <- ifelse((TP2 + FP2) > 0, TP2 / (TP2 + FP2), NA)
  recall2 <- ifelse(total_positives_discarded > 0, TP2 / total_positives_discarded, NA)
  F12 <- ifelse(is.na(precision2) | is.na(recall2) | (precision2 + recall2 == 0),
                NA, 2 * (precision2 * recall2) / (precision2 + recall2))
  FDR2 <- ifelse((TP2 + FP2) > 0, FP2 / (TP2 + FP2), NA)

  
  
  metrics_2 <- data.frame(criterion = "k >= sig_upper",
                          TP = TP2, FP = FP2,
                          precision = precision2, recall = recall2,
                          FDR = FDR2, F1 = F12)
  
  
  #Approach 3: Metrics if k >= sig_lower
  reintroduced_flags3 <- is_tp[crit_3]
  TP3 <- sum(reintroduced_flags3)
  FP3 <- sum(!reintroduced_flags3)
  
  precision3 <- ifelse((TP3 + FP3) > 0, TP3 / (TP3 + FP3), NA)
  recall3 <- ifelse(total_positives_discarded > 0, TP3 / total_positives_discarded, NA)
  F13 <- ifelse(is.na(precision3) | is.na(recall3) | (precision3 + recall3 == 0),
                NA, 2 * (precision3 * recall3) / (precision3 + recall3))
  FDR3 <- ifelse((TP3 + FP3) > 0, FP3 / (TP3 + FP3), NA)
  
  
  metrics_3 <- data.frame(criterion = "k >= sig_lower",
                          TP = TP3, FP = FP3,
                          precision = precision3, recall = recall3,
                          FDR = FDR3, F1 = F13)
 
 
  metrics <- rbind(metrics_1,metrics_2, metrics_3)

  
  k_max <- max(k_values, na.rm = TRUE)
  delta_k <- k_max - k_0
  
  # Dataframe por ASV
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    delta_k = delta_k_all,
    BC_values = BC_values,
    J_values = J_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = k_max,
    delta_k = delta_k,
    recovered_asv_ap1 = added_asv[crit_1],
    recovered_asv_ap2 = added_asv[crit_2],
    recovered_asv_ap3 = added_asv[crit_3],
    thresholds = data.frame(upper = sig_threshold_upper, 
                            lower = sig_threshold_lower,
                            sd = delta_sd),
    filtered_values = data.frame(k = k_0, BC = BC_0, J = J_0),
    unfiltered_values = data.frame(k = k_unfiltered, BC = BC_unfiltered, J = J_unfiltered),
    global_metrics = metrics,
    time = time.take
  ))
}

# K_validation with conservative filter 2 pairwaise
k_validation_pairwaise2 <- function(OTU_table, ground_truth, alpha, beta, mode, sig_level = 1.96){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- calculate_k_w_single(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  # apply conservative filter
  filter_output <- conservative_filter2(asv)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- calculate_k_w_single(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- BC_values <- J_values <- c()
  true_positive_levels <- c() 
  
  
  blast_nuc_matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- calculate_k_w_single(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    j_i <- mean(na.omit(evenness_k(asv_table_i)))
    
    # Evaluate with ground_truth depending on the chosen mode
    
    match_level <- "No Match"
    
    if (mode == "taxonomy") {
      asv_species <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
      asv_genus   <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
      
      if (asv_species %in% ground_truth$species || asv_genus %in% ground_truth$genus) {
        match_level <- "Match"
      }
      
    } else if (mode == "sequence") {
      asv_identity <- seq_dict$Sequence[match(asv_name, seq_dict$ASV)]
      
      match_matrix <- sapply(ground_truth$Sequence, function(gt_seq) {
        
        alignment <- pwalign::pairwiseAlignment(as.character(asv_identity), as.character(gt_seq), type = "local",
                                                substitutionMatrix = blast_nuc_matrix, gapOpening =5,
                                                gapExtension = 2)
        
        if (score(alignment) <= 0) {
          return(c(pid = 0, qcov = 0)) # No alignment
        }
        
        percent_identity <- pwalign::pid(alignment, type = "PID1")
        
        # Query Coverage
        pattern_string <- as.character(alignment@pattern)
        pattern_no_gaps <- gsub("-", "", pattern_string)
        query_coverage <- (nchar(pattern_no_gaps) / nchar(asv_identity)) * 100
        
        return(c(pid = percent_identity, qcov = query_coverage))
      })
      
      
      if (is.vector(match_matrix)) {
        match_matrix <- t(as.matrix(match_matrix))
      }else {
        match_matrix <- t(match_matrix)
      }
      
      if (nrow(match_matrix) == 0) {
        best_pid <- 0; best_qcov <- 0
      } else {
        best_match_row <- match_matrix[order(match_matrix[, "qcov"],
                                             match_matrix[, "pid"],
                                             decreasing = TRUE)[1], ]
        best_pid <- best_match_row["pid"]
        best_qcov <- best_match_row["qcov"]
      }
      
      
      best_pid_rounded  <- round(best_pid, 2)
      best_qcov_rounded <- round(best_qcov, 2)
      
      # Classification
      if (best_pid >= 98 && best_qcov >= 100) {
        match_level <- "Match"
      }
    }
    
    
    added_asv <- c(added_asv, asv_name)
    k_values <- c(k_values, k_i)
    BC_values <- c(BC_values, bc_i)
    J_values <- c(J_values, j_i)
    true_positive_levels <- c(true_positive_levels, match_level)
    
  }
  
  is_tp <- true_positive_levels != "No Match"
  total_positives_discarded <- sum(is_tp)
  delta_k_all <- k_values - k_0
  
  # Calculate thresholds
  delta_sd <- sd(delta_k_all, na.rm = TRUE)
  sig_threshold_upper <- sig_level * delta_sd
  sig_threshold_lower <- -sig_level * delta_sd
  
  
  # Criteria
  crit_1 <- which(k_values >= k_0) 
  crit_2 <- which(delta_k_all >= sig_threshold_upper)
  crit_3 <- which(delta_k_all >= sig_threshold_lower) 
  
  # Approach 1: Metrics if k >= 0
  reintroduced_flags1 <- is_tp[crit_1]
  
  TP1 <- sum(reintroduced_flags1)
  FP1 <- sum(!reintroduced_flags1)
  
  precision1 <- ifelse((TP1 + FP1) > 0, TP1 / (TP1 + FP1), NA)
  
  
  
  recall1 <- ifelse(total_positives_discarded > 0, TP1 / total_positives_discarded, NA)
  F11 <- ifelse(is.na(precision1) | is.na(recall1) | (precision1 + recall1 == 0),
                NA,
                2 * (precision1 * recall1) / (precision1 + recall1))
  FDR1 <- ifelse((TP1 + FP1) > 0, FP1 / (TP1 + FP1), NA)
  
  
  metrics_1 <- data.frame(criterion = "k >= k_0",
                          TP = TP1, FP = FP1,
                          precision = precision1, recall = recall1,
                          FDR = FDR1, F1 = F11)
  
  
  
  # Approach 2: Metrics >= sig_upper (significant increase)
  reintroduced_flags2 <- is_tp[crit_2]
  TP2 <- sum(reintroduced_flags2)
  FP2 <- sum(!reintroduced_flags2)
  
  precision2 <- ifelse((TP2 + FP2) > 0, TP2 / (TP2 + FP2), NA)
  recall2 <- ifelse(total_positives_discarded > 0, TP2 / total_positives_discarded, NA)
  F12 <- ifelse(is.na(precision2) | is.na(recall2) | (precision2 + recall2 == 0),
                NA, 2 * (precision2 * recall2) / (precision2 + recall2))
  FDR2 <- ifelse((TP2 + FP2) > 0, FP2 / (TP2 + FP2), NA)
  
  
  
  metrics_2 <- data.frame(criterion = "k >= sig_upper",
                          TP = TP2, FP = FP2,
                          precision = precision2, recall = recall2,
                          FDR = FDR2, F1 = F12)
  
  
  #Approach 3: Metrics if k >= sig_lower
  reintroduced_flags3 <- is_tp[crit_3]
  TP3 <- sum(reintroduced_flags3)
  FP3 <- sum(!reintroduced_flags3)
  
  precision3 <- ifelse((TP3 + FP3) > 0, TP3 / (TP3 + FP3), NA)
  recall3 <- ifelse(total_positives_discarded > 0, TP3 / total_positives_discarded, NA)
  F13 <- ifelse(is.na(precision3) | is.na(recall3) | (precision3 + recall3 == 0),
                NA, 2 * (precision3 * recall3) / (precision3 + recall3))
  FDR3 <- ifelse((TP3 + FP3) > 0, FP3 / (TP3 + FP3), NA)
  
  
  metrics_3 <- data.frame(criterion = "k >= sig_lower",
                          TP = TP3, FP = FP3,
                          precision = precision3, recall = recall3,
                          FDR = FDR3, F1 = F13)
  
  
  metrics <- rbind(metrics_1,metrics_2, metrics_3)
  
  
  k_max <- max(k_values, na.rm = TRUE)
  delta_k <- k_max - k_0
  
  # Dataframe por ASV
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    delta_k = delta_k_all,
    BC_values = BC_values,
    J_values = J_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = k_max,
    delta_k = delta_k,
    recovered_asv_ap1 = added_asv[crit_1],
    recovered_asv_ap2 = added_asv[crit_2],
    recovered_asv_ap3 = added_asv[crit_3],
    thresholds = data.frame(upper = sig_threshold_upper, 
                            lower = sig_threshold_lower,
                            sd = delta_sd),
    filtered_values = data.frame(k = k_0, BC = BC_0, J = J_0),
    unfiltered_values = data.frame(k = k_unfiltered, BC = BC_unfiltered, J = J_unfiltered),
    global_metrics = metrics,
    time = time.take
  ))
}

# K_validation with conservative filter 3 pairwaise
k_validation_pairwaise3 <- function(OTU_table, ground_truth, alpha, beta, mode, sig_level = 1.96){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- calculate_k_w_single(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  # apply conservative filter
  filter_output <- conservative_filter3(asv)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- calculate_k_w_single(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- BC_values <- J_values <- c()
  true_positive_levels <- c() 
  
  
  blast_nuc_matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- calculate_k_w_single(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    j_i <- mean(na.omit(evenness_k(asv_table_i)))
    
    # Evaluate with ground_truth depending on the chosen mode
    
    match_level <- "No Match"
    
    if (mode == "taxonomy") {
      asv_species <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
      asv_genus   <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
      
      if (asv_species %in% ground_truth$species || asv_genus %in% ground_truth$genus) {
        match_level <- "Match"
      }
      
    } else if (mode == "sequence") {
      asv_identity <- seq_dict$Sequence[match(asv_name, seq_dict$ASV)]
      
      match_matrix <- sapply(ground_truth$Sequence, function(gt_seq) {
        
        alignment <- pwalign::pairwiseAlignment(as.character(asv_identity), as.character(gt_seq), type = "local",
                                                substitutionMatrix = blast_nuc_matrix, gapOpening =5,
                                                gapExtension = 2)
        
        if (score(alignment) <= 0) {
          return(c(pid = 0, qcov = 0)) # No alignment
        }
        
        percent_identity <- pwalign::pid(alignment, type = "PID1")
        
        # Query Coverage
        pattern_string <- as.character(alignment@pattern)
        pattern_no_gaps <- gsub("-", "", pattern_string)
        query_coverage <- (nchar(pattern_no_gaps) / nchar(asv_identity)) * 100
        
        return(c(pid = percent_identity, qcov = query_coverage))
      })
      
      
      if (is.vector(match_matrix)) {
        match_matrix <- t(as.matrix(match_matrix))
      }else {
        match_matrix <- t(match_matrix)
      }
      
      if (nrow(match_matrix) == 0) {
        best_pid <- 0; best_qcov <- 0
      } else {
        best_match_row <- match_matrix[order(match_matrix[, "qcov"],
                                             match_matrix[, "pid"],
                                             decreasing = TRUE)[1], ]
        best_pid <- best_match_row["pid"]
        best_qcov <- best_match_row["qcov"]
      }
      
      
      best_pid_rounded  <- round(best_pid, 2)
      best_qcov_rounded <- round(best_qcov, 2)
      
      # Classification
      if (best_pid >= 98 && best_qcov >= 100) {
        match_level <- "Match"
      }
    }
    
    
    added_asv <- c(added_asv, asv_name)
    k_values <- c(k_values, k_i)
    BC_values <- c(BC_values, bc_i)
    J_values <- c(J_values, j_i)
    true_positive_levels <- c(true_positive_levels, match_level)
    
  }
  
  is_tp <- true_positive_levels != "No Match"
  total_positives_discarded <- sum(is_tp)
  delta_k_all <- k_values - k_0
  
  # Calculate thresholds
  delta_sd <- sd(delta_k_all, na.rm = TRUE)
  sig_threshold_upper <- sig_level * delta_sd
  sig_threshold_lower <- -sig_level * delta_sd
  
  
  # Criteria
  crit_1 <- which(k_values >= k_0) 
  crit_2 <- which(delta_k_all >= sig_threshold_upper)
  crit_3 <- which(delta_k_all >= sig_threshold_lower) 
  
  # Approach 1: Metrics if k >= 0
  reintroduced_flags1 <- is_tp[crit_1]
  
  TP1 <- sum(reintroduced_flags1)
  FP1 <- sum(!reintroduced_flags1)
  
  precision1 <- ifelse((TP1 + FP1) > 0, TP1 / (TP1 + FP1), NA)
  
  
  
  recall1 <- ifelse(total_positives_discarded > 0, TP1 / total_positives_discarded, NA)
  F11 <- ifelse(is.na(precision1) | is.na(recall1) | (precision1 + recall1 == 0),
                NA,
                2 * (precision1 * recall1) / (precision1 + recall1))
  FDR1 <- ifelse((TP1 + FP1) > 0, FP1 / (TP1 + FP1), NA)
  
  
  metrics_1 <- data.frame(criterion = "k >= k_0",
                          TP = TP1, FP = FP1,
                          precision = precision1, recall = recall1,
                          FDR = FDR1, F1 = F11)
  
  
  
  # Approach 2: Metrics >= sig_upper (significant increase)
  reintroduced_flags2 <- is_tp[crit_2]
  TP2 <- sum(reintroduced_flags2)
  FP2 <- sum(!reintroduced_flags2)
  
  precision2 <- ifelse((TP2 + FP2) > 0, TP2 / (TP2 + FP2), NA)
  recall2 <- ifelse(total_positives_discarded > 0, TP2 / total_positives_discarded, NA)
  F12 <- ifelse(is.na(precision2) | is.na(recall2) | (precision2 + recall2 == 0),
                NA, 2 * (precision2 * recall2) / (precision2 + recall2))
  FDR2 <- ifelse((TP2 + FP2) > 0, FP2 / (TP2 + FP2), NA)
  
  
  
  metrics_2 <- data.frame(criterion = "k >= sig_upper",
                          TP = TP2, FP = FP2,
                          precision = precision2, recall = recall2,
                          FDR = FDR2, F1 = F12)
  
  
  #Approach 3: Metrics if k >= sig_lower
  reintroduced_flags3 <- is_tp[crit_3]
  TP3 <- sum(reintroduced_flags3)
  FP3 <- sum(!reintroduced_flags3)
  
  precision3 <- ifelse((TP3 + FP3) > 0, TP3 / (TP3 + FP3), NA)
  recall3 <- ifelse(total_positives_discarded > 0, TP3 / total_positives_discarded, NA)
  F13 <- ifelse(is.na(precision3) | is.na(recall3) | (precision3 + recall3 == 0),
                NA, 2 * (precision3 * recall3) / (precision3 + recall3))
  FDR3 <- ifelse((TP3 + FP3) > 0, FP3 / (TP3 + FP3), NA)
  
  
  metrics_3 <- data.frame(criterion = "k >= sig_lower",
                          TP = TP3, FP = FP3,
                          precision = precision3, recall = recall3,
                          FDR = FDR3, F1 = F13)
  
  
  metrics <- rbind(metrics_1,metrics_2, metrics_3)
  
  
  k_max <- max(k_values, na.rm = TRUE)
  delta_k <- k_max - k_0
  
  # Dataframe por ASV
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    delta_k = delta_k_all,
    BC_values = BC_values,
    J_values = J_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = k_max,
    delta_k = delta_k,
    recovered_asv_ap1 = added_asv[crit_1],
    recovered_asv_ap2 = added_asv[crit_2],
    recovered_asv_ap3 = added_asv[crit_3],
    thresholds = data.frame(upper = sig_threshold_upper, 
                            lower = sig_threshold_lower,
                            sd = delta_sd),
    filtered_values = data.frame(k = k_0, BC = BC_0, J = J_0),
    unfiltered_values = data.frame(k = k_unfiltered, BC = BC_unfiltered, J = J_unfiltered),
    global_metrics = metrics,
    time = time.take
  ))
}

# K_validation with conservative filter 4 pairwaise
k_validation_pairwaise4 <- function(OTU_table, ground_truth, alpha, beta, mode, sig_level = 1.96){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- calculate_k_w_single(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  # apply conservative filter
  filter_output <- conservative_filter4(asv)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- calculate_k_w_single(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- BC_values <- J_values <- c()
  true_positive_levels <- c() 
  
  
  blast_nuc_matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- calculate_k_w_single(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    j_i <- mean(na.omit(evenness_k(asv_table_i)))
    
    # Evaluate with ground_truth depending on the chosen mode
    
    match_level <- "No Match"
    
    if (mode == "taxonomy") {
      asv_species <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
      asv_genus   <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
      
      if (asv_species %in% ground_truth$species || asv_genus %in% ground_truth$genus) {
        match_level <- "Match"
      }
      
    } else if (mode == "sequence") {
      asv_identity <- seq_dict$Sequence[match(asv_name, seq_dict$ASV)]
      
      match_matrix <- sapply(ground_truth$Sequence, function(gt_seq) {
        
        alignment <- pwalign::pairwiseAlignment(as.character(asv_identity), as.character(gt_seq), type = "local",
                                                substitutionMatrix = blast_nuc_matrix, gapOpening =5,
                                                gapExtension = 2)
        
        if (score(alignment) <= 0) {
          return(c(pid = 0, qcov = 0)) # No alignment
        }
        
        percent_identity <- pwalign::pid(alignment, type = "PID1")
        
        # Query Coverage
        pattern_string <- as.character(alignment@pattern)
        pattern_no_gaps <- gsub("-", "", pattern_string)
        query_coverage <- (nchar(pattern_no_gaps) / nchar(asv_identity)) * 100
        
        return(c(pid = percent_identity, qcov = query_coverage))
      })
      
      
      if (is.vector(match_matrix)) {
        match_matrix <- t(as.matrix(match_matrix))
      }else {
        match_matrix <- t(match_matrix)
      }
      
      if (nrow(match_matrix) == 0) {
        best_pid <- 0; best_qcov <- 0
      } else {
        best_match_row <- match_matrix[order(match_matrix[, "qcov"],
                                             match_matrix[, "pid"],
                                             decreasing = TRUE)[1], ]
        best_pid <- best_match_row["pid"]
        best_qcov <- best_match_row["qcov"]
      }
      
      
      best_pid_rounded  <- round(best_pid, 2)
      best_qcov_rounded <- round(best_qcov, 2)
      
      # Classification
      if (best_pid >= 98 && best_qcov >= 100) {
        match_level <- "Match"
      }
    }
    
    
    added_asv <- c(added_asv, asv_name)
    k_values <- c(k_values, k_i)
    BC_values <- c(BC_values, bc_i)
    J_values <- c(J_values, j_i)
    true_positive_levels <- c(true_positive_levels, match_level)
    
  }
  
  is_tp <- true_positive_levels != "No Match"
  total_positives_discarded <- sum(is_tp)
  delta_k_all <- k_values - k_0
  
  # Calculate thresholds
  delta_sd <- sd(delta_k_all, na.rm = TRUE)
  sig_threshold_upper <- sig_level * delta_sd
  sig_threshold_lower <- -sig_level * delta_sd
  
  
  # Criteria
  crit_1 <- which(k_values >= k_0) 
  crit_2 <- which(delta_k_all >= sig_threshold_upper)
  crit_3 <- which(delta_k_all >= sig_threshold_lower) 
  
  # Approach 1: Metrics if k >= 0
  reintroduced_flags1 <- is_tp[crit_1]
  
  TP1 <- sum(reintroduced_flags1)
  FP1 <- sum(!reintroduced_flags1)
  
  precision1 <- ifelse((TP1 + FP1) > 0, TP1 / (TP1 + FP1), NA)
  
  
  
  recall1 <- ifelse(total_positives_discarded > 0, TP1 / total_positives_discarded, NA)
  F11 <- ifelse(is.na(precision1) | is.na(recall1) | (precision1 + recall1 == 0),
                NA,
                2 * (precision1 * recall1) / (precision1 + recall1))
  FDR1 <- ifelse((TP1 + FP1) > 0, FP1 / (TP1 + FP1), NA)
  
  
  metrics_1 <- data.frame(criterion = "k >= k_0",
                          TP = TP1, FP = FP1,
                          precision = precision1, recall = recall1,
                          FDR = FDR1, F1 = F11)
  
  
  
  # Approach 2: Metrics >= sig_upper (significant increase)
  reintroduced_flags2 <- is_tp[crit_2]
  TP2 <- sum(reintroduced_flags2)
  FP2 <- sum(!reintroduced_flags2)
  
  precision2 <- ifelse((TP2 + FP2) > 0, TP2 / (TP2 + FP2), NA)
  recall2 <- ifelse(total_positives_discarded > 0, TP2 / total_positives_discarded, NA)
  F12 <- ifelse(is.na(precision2) | is.na(recall2) | (precision2 + recall2 == 0),
                NA, 2 * (precision2 * recall2) / (precision2 + recall2))
  FDR2 <- ifelse((TP2 + FP2) > 0, FP2 / (TP2 + FP2), NA)
  
  
  
  metrics_2 <- data.frame(criterion = "k >= sig_upper",
                          TP = TP2, FP = FP2,
                          precision = precision2, recall = recall2,
                          FDR = FDR2, F1 = F12)
  
  
  #Approach 3: Metrics if k >= sig_lower
  reintroduced_flags3 <- is_tp[crit_3]
  TP3 <- sum(reintroduced_flags3)
  FP3 <- sum(!reintroduced_flags3)
  
  precision3 <- ifelse((TP3 + FP3) > 0, TP3 / (TP3 + FP3), NA)
  recall3 <- ifelse(total_positives_discarded > 0, TP3 / total_positives_discarded, NA)
  F13 <- ifelse(is.na(precision3) | is.na(recall3) | (precision3 + recall3 == 0),
                NA, 2 * (precision3 * recall3) / (precision3 + recall3))
  FDR3 <- ifelse((TP3 + FP3) > 0, FP3 / (TP3 + FP3), NA)
  
  
  metrics_3 <- data.frame(criterion = "k >= sig_lower",
                          TP = TP3, FP = FP3,
                          precision = precision3, recall = recall3,
                          FDR = FDR3, F1 = F13)
  
  
  metrics <- rbind(metrics_1,metrics_2, metrics_3)
  
  
  k_max <- max(k_values, na.rm = TRUE)
  delta_k <- k_max - k_0
  
  # Dataframe por ASV
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    delta_k = delta_k_all,
    BC_values = BC_values,
    J_values = J_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = k_max,
    delta_k = delta_k,
    recovered_asv_ap1 = added_asv[crit_1],
    recovered_asv_ap2 = added_asv[crit_2],
    recovered_asv_ap3 = added_asv[crit_3],
    thresholds = data.frame(upper = sig_threshold_upper, 
                            lower = sig_threshold_lower,
                            sd = delta_sd),
    filtered_values = data.frame(k = k_0, BC = BC_0, J = J_0),
    unfiltered_values = data.frame(k = k_unfiltered, BC = BC_unfiltered, J = J_unfiltered),
    global_metrics = metrics,
    time = time.take
  ))
}

# K_validation with conservative filter 5 pairwaise
k_validation_pairwaise5 <- function(OTU_table, ground_truth, alpha, beta, mode, sig_level = 1.96){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- calculate_k_w_single(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  # apply conservative filter
  filter_output <- conservative_filter5(asv)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- calculate_k_w_single(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- BC_values <- J_values <- c()
  true_positive_levels <- c() 
  
  
  blast_nuc_matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- calculate_k_w_single(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    j_i <- mean(na.omit(evenness_k(asv_table_i)))
    
    # Evaluate with ground_truth depending on the chosen mode
    
    match_level <- "No Match"
    
    if (mode == "taxonomy") {
      asv_species <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
      asv_genus   <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
      
      if (asv_species %in% ground_truth$species || asv_genus %in% ground_truth$genus) {
        match_level <- "Match"
      }
      
    } else if (mode == "sequence") {
      asv_identity <- seq_dict$Sequence[match(asv_name, seq_dict$ASV)]
      
      match_matrix <- sapply(ground_truth$Sequence, function(gt_seq) {
        
        alignment <- pwalign::pairwiseAlignment(as.character(asv_identity), as.character(gt_seq), type = "local",
                                                substitutionMatrix = blast_nuc_matrix, gapOpening =5,
                                                gapExtension = 2)
        
        if (score(alignment) <= 0) {
          return(c(pid = 0, qcov = 0)) # No alignment
        }
        
        percent_identity <- pwalign::pid(alignment, type = "PID1")
        
        # Query Coverage
        pattern_string <- as.character(alignment@pattern)
        pattern_no_gaps <- gsub("-", "", pattern_string)
        query_coverage <- (nchar(pattern_no_gaps) / nchar(asv_identity)) * 100
        
        return(c(pid = percent_identity, qcov = query_coverage))
      })
      
      
      if (is.vector(match_matrix)) {
        match_matrix <- t(as.matrix(match_matrix))
      }else {
        match_matrix <- t(match_matrix)
      }
      
      if (nrow(match_matrix) == 0) {
        best_pid <- 0; best_qcov <- 0
      } else {
        best_match_row <- match_matrix[order(match_matrix[, "qcov"],
                                             match_matrix[, "pid"],
                                             decreasing = TRUE)[1], ]
        best_pid <- best_match_row["pid"]
        best_qcov <- best_match_row["qcov"]
      }
      
      
      best_pid_rounded  <- round(best_pid, 2)
      best_qcov_rounded <- round(best_qcov, 2)
      
      # Classification
      if (best_pid >= 98 && best_qcov >= 100) {
        match_level <- "Match"
      }
    }
    
    
    added_asv <- c(added_asv, asv_name)
    k_values <- c(k_values, k_i)
    BC_values <- c(BC_values, bc_i)
    J_values <- c(J_values, j_i)
    true_positive_levels <- c(true_positive_levels, match_level)
    
  }
  
  is_tp <- true_positive_levels != "No Match"
  total_positives_discarded <- sum(is_tp)
  delta_k_all <- k_values - k_0
  
  # Calculate thresholds
  delta_sd <- sd(delta_k_all, na.rm = TRUE)
  sig_threshold_upper <- sig_level * delta_sd
  sig_threshold_lower <- -sig_level * delta_sd
  
  
  # Criteria
  crit_1 <- which(k_values >= k_0) 
  crit_2 <- which(delta_k_all >= sig_threshold_upper)
  crit_3 <- which(delta_k_all >= sig_threshold_lower) 
  
  # Approach 1: Metrics if k >= 0
  reintroduced_flags1 <- is_tp[crit_1]
  
  TP1 <- sum(reintroduced_flags1)
  FP1 <- sum(!reintroduced_flags1)
  
  precision1 <- ifelse((TP1 + FP1) > 0, TP1 / (TP1 + FP1), NA)
  
  
  
  recall1 <- ifelse(total_positives_discarded > 0, TP1 / total_positives_discarded, NA)
  F11 <- ifelse(is.na(precision1) | is.na(recall1) | (precision1 + recall1 == 0),
                NA,
                2 * (precision1 * recall1) / (precision1 + recall1))
  FDR1 <- ifelse((TP1 + FP1) > 0, FP1 / (TP1 + FP1), NA)
  
  
  metrics_1 <- data.frame(criterion = "k >= k_0",
                          TP = TP1, FP = FP1,
                          precision = precision1, recall = recall1,
                          FDR = FDR1, F1 = F11)
  
  
  
  # Approach 2: Metrics >= sig_upper (significant increase)
  reintroduced_flags2 <- is_tp[crit_2]
  TP2 <- sum(reintroduced_flags2)
  FP2 <- sum(!reintroduced_flags2)
  
  precision2 <- ifelse((TP2 + FP2) > 0, TP2 / (TP2 + FP2), NA)
  recall2 <- ifelse(total_positives_discarded > 0, TP2 / total_positives_discarded, NA)
  F12 <- ifelse(is.na(precision2) | is.na(recall2) | (precision2 + recall2 == 0),
                NA, 2 * (precision2 * recall2) / (precision2 + recall2))
  FDR2 <- ifelse((TP2 + FP2) > 0, FP2 / (TP2 + FP2), NA)
  
  
  
  metrics_2 <- data.frame(criterion = "k >= sig_upper",
                          TP = TP2, FP = FP2,
                          precision = precision2, recall = recall2,
                          FDR = FDR2, F1 = F12)
  
  
  #Approach 3: Metrics if k >= sig_lower
  reintroduced_flags3 <- is_tp[crit_3]
  TP3 <- sum(reintroduced_flags3)
  FP3 <- sum(!reintroduced_flags3)
  
  precision3 <- ifelse((TP3 + FP3) > 0, TP3 / (TP3 + FP3), NA)
  recall3 <- ifelse(total_positives_discarded > 0, TP3 / total_positives_discarded, NA)
  F13 <- ifelse(is.na(precision3) | is.na(recall3) | (precision3 + recall3 == 0),
                NA, 2 * (precision3 * recall3) / (precision3 + recall3))
  FDR3 <- ifelse((TP3 + FP3) > 0, FP3 / (TP3 + FP3), NA)
  
  
  metrics_3 <- data.frame(criterion = "k >= sig_lower",
                          TP = TP3, FP = FP3,
                          precision = precision3, recall = recall3,
                          FDR = FDR3, F1 = F13)
  
  
  metrics <- rbind(metrics_1,metrics_2, metrics_3)
  
  
  k_max <- max(k_values, na.rm = TRUE)
  delta_k <- k_max - k_0
  
  # Dataframe por ASV
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    delta_k = delta_k_all,
    BC_values = BC_values,
    J_values = J_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = k_max,
    delta_k = delta_k,
    recovered_asv_ap1 = added_asv[crit_1],
    recovered_asv_ap2 = added_asv[crit_2],
    recovered_asv_ap3 = added_asv[crit_3],
    thresholds = data.frame(upper = sig_threshold_upper, 
                            lower = sig_threshold_lower,
                            sd = delta_sd),
    filtered_values = data.frame(k = k_0, BC = BC_0, J = J_0),
    unfiltered_values = data.frame(k = k_unfiltered, BC = BC_unfiltered, J = J_unfiltered),
    global_metrics = metrics,
    time = time.take
  ))
}

# K_validation with conservative filter 6 pairwaise
k_validation_pairwaise6 <- function(OTU_table, ground_truth, alpha, beta, mode, sig_level = 1.96){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- calculate_k_w_single(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  # apply conservative filter
  filter_output <- conservative_filter6(asv)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- calculate_k_w_single(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- BC_values <- J_values <- c()
  true_positive_levels <- c() 
  
  
  blast_nuc_matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- calculate_k_w_single(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    j_i <- mean(na.omit(evenness_k(asv_table_i)))
    
    # Evaluate with ground_truth depending on the chosen mode
    
    match_level <- "No Match"
    
    if (mode == "taxonomy") {
      asv_species <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
      asv_genus   <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
      
      if (asv_species %in% ground_truth$species || asv_genus %in% ground_truth$genus) {
        match_level <- "Match"
      }
      
    } else if (mode == "sequence") {
      asv_identity <- seq_dict$Sequence[match(asv_name, seq_dict$ASV)]
      
      match_matrix <- sapply(ground_truth$Sequence, function(gt_seq) {
        
        alignment <- pwalign::pairwiseAlignment(as.character(asv_identity), as.character(gt_seq), type = "local",
                                                substitutionMatrix = blast_nuc_matrix, gapOpening =5,
                                                gapExtension = 2)
        
        if (score(alignment) <= 0) {
          return(c(pid = 0, qcov = 0)) # No alignment
        }
        
        percent_identity <- pwalign::pid(alignment, type = "PID1")
        
        # Query Coverage
        pattern_string <- as.character(alignment@pattern)
        pattern_no_gaps <- gsub("-", "", pattern_string)
        query_coverage <- (nchar(pattern_no_gaps) / nchar(asv_identity)) * 100
        
        return(c(pid = percent_identity, qcov = query_coverage))
      })
      
      
      if (is.vector(match_matrix)) {
        match_matrix <- t(as.matrix(match_matrix))
      }else {
        match_matrix <- t(match_matrix)
      }
      
      if (nrow(match_matrix) == 0) {
        best_pid <- 0; best_qcov <- 0
      } else {
        best_match_row <- match_matrix[order(match_matrix[, "qcov"],
                                             match_matrix[, "pid"],
                                             decreasing = TRUE)[1], ]
        best_pid <- best_match_row["pid"]
        best_qcov <- best_match_row["qcov"]
      }
      
      
      best_pid_rounded  <- round(best_pid, 2)
      best_qcov_rounded <- round(best_qcov, 2)
      
      # Classification
      if (best_pid >= 98 && best_qcov >= 100) {
        match_level <- "Match"
      }
    }
    
    
    added_asv <- c(added_asv, asv_name)
    k_values <- c(k_values, k_i)
    BC_values <- c(BC_values, bc_i)
    J_values <- c(J_values, j_i)
    true_positive_levels <- c(true_positive_levels, match_level)
    
  }
  
  is_tp <- true_positive_levels != "No Match"
  total_positives_discarded <- sum(is_tp)
  delta_k_all <- k_values - k_0
  
  # Calculate thresholds
  delta_sd <- sd(delta_k_all, na.rm = TRUE)
  sig_threshold_upper <- sig_level * delta_sd
  sig_threshold_lower <- -sig_level * delta_sd
  
  
  # Criteria
  crit_1 <- which(k_values >= k_0) 
  crit_2 <- which(delta_k_all >= sig_threshold_upper)
  crit_3 <- which(delta_k_all >= sig_threshold_lower) 
  
  # Approach 1: Metrics if k >= 0
  reintroduced_flags1 <- is_tp[crit_1]
  
  TP1 <- sum(reintroduced_flags1)
  FP1 <- sum(!reintroduced_flags1)
  
  precision1 <- ifelse((TP1 + FP1) > 0, TP1 / (TP1 + FP1), NA)
  
  
  
  recall1 <- ifelse(total_positives_discarded > 0, TP1 / total_positives_discarded, NA)
  F11 <- ifelse(is.na(precision1) | is.na(recall1) | (precision1 + recall1 == 0),
                NA,
                2 * (precision1 * recall1) / (precision1 + recall1))
  FDR1 <- ifelse((TP1 + FP1) > 0, FP1 / (TP1 + FP1), NA)
  
  
  metrics_1 <- data.frame(criterion = "k >= k_0",
                          TP = TP1, FP = FP1,
                          precision = precision1, recall = recall1,
                          FDR = FDR1, F1 = F11)
  
  
  
  # Approach 2: Metrics >= sig_upper (significant increase)
  reintroduced_flags2 <- is_tp[crit_2]
  TP2 <- sum(reintroduced_flags2)
  FP2 <- sum(!reintroduced_flags2)
  
  precision2 <- ifelse((TP2 + FP2) > 0, TP2 / (TP2 + FP2), NA)
  recall2 <- ifelse(total_positives_discarded > 0, TP2 / total_positives_discarded, NA)
  F12 <- ifelse(is.na(precision2) | is.na(recall2) | (precision2 + recall2 == 0),
                NA, 2 * (precision2 * recall2) / (precision2 + recall2))
  FDR2 <- ifelse((TP2 + FP2) > 0, FP2 / (TP2 + FP2), NA)
  
  
  
  metrics_2 <- data.frame(criterion = "k >= sig_upper",
                          TP = TP2, FP = FP2,
                          precision = precision2, recall = recall2,
                          FDR = FDR2, F1 = F12)
  
  
  #Approach 3: Metrics if k >= sig_lower
  reintroduced_flags3 <- is_tp[crit_3]
  TP3 <- sum(reintroduced_flags3)
  FP3 <- sum(!reintroduced_flags3)
  
  precision3 <- ifelse((TP3 + FP3) > 0, TP3 / (TP3 + FP3), NA)
  recall3 <- ifelse(total_positives_discarded > 0, TP3 / total_positives_discarded, NA)
  F13 <- ifelse(is.na(precision3) | is.na(recall3) | (precision3 + recall3 == 0),
                NA, 2 * (precision3 * recall3) / (precision3 + recall3))
  FDR3 <- ifelse((TP3 + FP3) > 0, FP3 / (TP3 + FP3), NA)
  
  
  metrics_3 <- data.frame(criterion = "k >= sig_lower",
                          TP = TP3, FP = FP3,
                          precision = precision3, recall = recall3,
                          FDR = FDR3, F1 = F13)
  
  
  metrics <- rbind(metrics_1,metrics_2, metrics_3)
  
  
  k_max <- max(k_values, na.rm = TRUE)
  delta_k <- k_max - k_0
  
  # Dataframe por ASV
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    delta_k = delta_k_all,
    BC_values = BC_values,
    J_values = J_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = k_max,
    delta_k = delta_k,
    recovered_asv_ap1 = added_asv[crit_1],
    recovered_asv_ap2 = added_asv[crit_2],
    recovered_asv_ap3 = added_asv[crit_3],
    thresholds = data.frame(upper = sig_threshold_upper, 
                            lower = sig_threshold_lower,
                            sd = delta_sd),
    filtered_values = data.frame(k = k_0, BC = BC_0, J = J_0),
    unfiltered_values = data.frame(k = k_unfiltered, BC = BC_unfiltered, J = J_unfiltered),
    global_metrics = metrics,
    time = time.take
  ))
}


# K_validation with conservative filter 7 pairwaise
k_validation_pairwaise7 <- function(OTU_table, ground_truth, alpha, beta, mode, sig_level = 1.96){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- calculate_k_w_single(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  # apply conservative filter
  filter_output <- conservative_filter7(asv)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- calculate_k_w_single(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- BC_values <- J_values <- c()
  true_positive_levels <- c() 
  
  
  blast_nuc_matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- calculate_k_w_single(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    j_i <- mean(na.omit(evenness_k(asv_table_i)))
    
    # Evaluate with ground_truth depending on the chosen mode
    
    match_level <- "No Match"
    
    if (mode == "taxonomy") {
      asv_species <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
      asv_genus   <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
      
      if (asv_species %in% ground_truth$species || asv_genus %in% ground_truth$genus) {
        match_level <- "Match"
      }
      
    } else if (mode == "sequence") {
      asv_identity <- seq_dict$Sequence[match(asv_name, seq_dict$ASV)]
      
      match_matrix <- sapply(ground_truth$Sequence, function(gt_seq) {
        
        alignment <- pwalign::pairwiseAlignment(as.character(asv_identity), as.character(gt_seq), type = "local",
                                                substitutionMatrix = blast_nuc_matrix, gapOpening =5,
                                                gapExtension = 2)
        
        if (score(alignment) <= 0) {
          return(c(pid = 0, qcov = 0)) # No alignment
        }
        
        percent_identity <- pwalign::pid(alignment, type = "PID1")
        
        # Query Coverage
        pattern_string <- as.character(alignment@pattern)
        pattern_no_gaps <- gsub("-", "", pattern_string)
        query_coverage <- (nchar(pattern_no_gaps) / nchar(asv_identity)) * 100
        
        return(c(pid = percent_identity, qcov = query_coverage))
      })
      
      
      if (is.vector(match_matrix)) {
        match_matrix <- t(as.matrix(match_matrix))
      }else {
        match_matrix <- t(match_matrix)
      }
      
      if (nrow(match_matrix) == 0) {
        best_pid <- 0; best_qcov <- 0
      } else {
        best_match_row <- match_matrix[order(match_matrix[, "qcov"],
                                             match_matrix[, "pid"],
                                             decreasing = TRUE)[1], ]
        best_pid <- best_match_row["pid"]
        best_qcov <- best_match_row["qcov"]
      }
      
      
      best_pid_rounded  <- round(best_pid, 2)
      best_qcov_rounded <- round(best_qcov, 2)
      
      # Classification
      if (best_pid >= 98 && best_qcov >= 100) {
        match_level <- "Match"
      }
    }
    
    
    added_asv <- c(added_asv, asv_name)
    k_values <- c(k_values, k_i)
    BC_values <- c(BC_values, bc_i)
    J_values <- c(J_values, j_i)
    true_positive_levels <- c(true_positive_levels, match_level)
    
  }
  
  is_tp <- true_positive_levels != "No Match"
  total_positives_discarded <- sum(is_tp)
  delta_k_all <- k_values - k_0
  
  # Calculate thresholds
  delta_sd <- sd(delta_k_all, na.rm = TRUE)
  sig_threshold_upper <- sig_level * delta_sd
  sig_threshold_lower <- -sig_level * delta_sd
  
  
  # Criteria
  crit_1 <- which(k_values >= k_0) 
  crit_2 <- which(delta_k_all >= sig_threshold_upper)
  crit_3 <- which(delta_k_all >= sig_threshold_lower) 
  
  # Approach 1: Metrics if k >= 0
  reintroduced_flags1 <- is_tp[crit_1]
  
  TP1 <- sum(reintroduced_flags1)
  FP1 <- sum(!reintroduced_flags1)
  
  precision1 <- ifelse((TP1 + FP1) > 0, TP1 / (TP1 + FP1), NA)
  
  
  
  recall1 <- ifelse(total_positives_discarded > 0, TP1 / total_positives_discarded, NA)
  F11 <- ifelse(is.na(precision1) | is.na(recall1) | (precision1 + recall1 == 0),
                NA,
                2 * (precision1 * recall1) / (precision1 + recall1))
  FDR1 <- ifelse((TP1 + FP1) > 0, FP1 / (TP1 + FP1), NA)
  
  
  metrics_1 <- data.frame(criterion = "k >= k_0",
                          TP = TP1, FP = FP1,
                          precision = precision1, recall = recall1,
                          FDR = FDR1, F1 = F11)
  
  
  
  # Approach 2: Metrics >= sig_upper (significant increase)
  reintroduced_flags2 <- is_tp[crit_2]
  TP2 <- sum(reintroduced_flags2)
  FP2 <- sum(!reintroduced_flags2)
  
  precision2 <- ifelse((TP2 + FP2) > 0, TP2 / (TP2 + FP2), NA)
  recall2 <- ifelse(total_positives_discarded > 0, TP2 / total_positives_discarded, NA)
  F12 <- ifelse(is.na(precision2) | is.na(recall2) | (precision2 + recall2 == 0),
                NA, 2 * (precision2 * recall2) / (precision2 + recall2))
  FDR2 <- ifelse((TP2 + FP2) > 0, FP2 / (TP2 + FP2), NA)
  
  
  
  metrics_2 <- data.frame(criterion = "k >= sig_upper",
                          TP = TP2, FP = FP2,
                          precision = precision2, recall = recall2,
                          FDR = FDR2, F1 = F12)
  
  
  #Approach 3: Metrics if k >= sig_lower
  reintroduced_flags3 <- is_tp[crit_3]
  TP3 <- sum(reintroduced_flags3)
  FP3 <- sum(!reintroduced_flags3)
  
  precision3 <- ifelse((TP3 + FP3) > 0, TP3 / (TP3 + FP3), NA)
  recall3 <- ifelse(total_positives_discarded > 0, TP3 / total_positives_discarded, NA)
  F13 <- ifelse(is.na(precision3) | is.na(recall3) | (precision3 + recall3 == 0),
                NA, 2 * (precision3 * recall3) / (precision3 + recall3))
  FDR3 <- ifelse((TP3 + FP3) > 0, FP3 / (TP3 + FP3), NA)
  
  
  metrics_3 <- data.frame(criterion = "k >= sig_lower",
                          TP = TP3, FP = FP3,
                          precision = precision3, recall = recall3,
                          FDR = FDR3, F1 = F13)
  
  
  metrics <- rbind(metrics_1,metrics_2, metrics_3)
  
  
  k_max <- max(k_values, na.rm = TRUE)
  delta_k <- k_max - k_0
  
  # Dataframe por ASV
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    delta_k = delta_k_all,
    BC_values = BC_values,
    J_values = J_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = k_max,
    delta_k = delta_k,
    recovered_asv_ap1 = added_asv[crit_1],
    recovered_asv_ap2 = added_asv[crit_2],
    recovered_asv_ap3 = added_asv[crit_3],
    thresholds = data.frame(upper = sig_threshold_upper, 
                            lower = sig_threshold_lower,
                            sd = delta_sd),
    filtered_values = data.frame(k = k_0, BC = BC_0, J = J_0),
    unfiltered_values = data.frame(k = k_unfiltered, BC = BC_unfiltered, J = J_unfiltered),
    global_metrics = metrics,
    time = time.take
  ))
}

k_validation_pairwaise8 <- function(OTU_table, ground_truth, alpha, beta, mode, sig_level = 1.96){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- calculate_k_w_single(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  # apply conservative filter
  filter_output <- conservative_filter8(asv)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- calculate_k_w_single(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- BC_values <- J_values <- c()
  true_positive_levels <- c() 
  
  
  blast_nuc_matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- calculate_k_w_single(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    j_i <- mean(na.omit(evenness_k(asv_table_i)))
    
    # Evaluate with ground_truth depending on the chosen mode
    
    match_level <- "No Match"
    
    if (mode == "taxonomy") {
      asv_species <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
      asv_genus   <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
      
      if (asv_species %in% ground_truth$species || asv_genus %in% ground_truth$genus) {
        match_level <- "Match"
      }
      
    } else if (mode == "sequence") {
      asv_identity <- seq_dict$Sequence[match(asv_name, seq_dict$ASV)]
      
      match_matrix <- sapply(ground_truth$Sequence, function(gt_seq) {
        
        alignment <- pwalign::pairwiseAlignment(as.character(asv_identity), as.character(gt_seq), type = "local",
                                                substitutionMatrix = blast_nuc_matrix, gapOpening =5,
                                                gapExtension = 2)
        
        if (score(alignment) <= 0) {
          return(c(pid = 0, qcov = 0)) # No alignment
        }
        
        percent_identity <- pwalign::pid(alignment, type = "PID1")
        
        # Query Coverage
        pattern_string <- as.character(alignment@pattern)
        pattern_no_gaps <- gsub("-", "", pattern_string)
        query_coverage <- (nchar(pattern_no_gaps) / nchar(asv_identity)) * 100
        
        return(c(pid = percent_identity, qcov = query_coverage))
      })
      
      
      if (is.vector(match_matrix)) {
        match_matrix <- t(as.matrix(match_matrix))
      }else {
        match_matrix <- t(match_matrix)
      }
      
      if (nrow(match_matrix) == 0) {
        best_pid <- 0; best_qcov <- 0
      } else {
        best_match_row <- match_matrix[order(match_matrix[, "qcov"],
                                             match_matrix[, "pid"],
                                             decreasing = TRUE)[1], ]
        best_pid <- best_match_row["pid"]
        best_qcov <- best_match_row["qcov"]
      }
      
      
      best_pid_rounded  <- round(best_pid, 2)
      best_qcov_rounded <- round(best_qcov, 2)
      
      # Classification
      if (best_pid >= 98 && best_qcov >= 100) {
        match_level <- "Match"
      }
    }
    
    
    added_asv <- c(added_asv, asv_name)
    k_values <- c(k_values, k_i)
    BC_values <- c(BC_values, bc_i)
    J_values <- c(J_values, j_i)
    true_positive_levels <- c(true_positive_levels, match_level)
    
  }
  
  is_tp <- true_positive_levels != "No Match"
  total_positives_discarded <- sum(is_tp)
  delta_k_all <- k_values - k_0
  
  # Calculate thresholds
  delta_sd <- sd(delta_k_all, na.rm = TRUE)
  sig_threshold_upper <- sig_level * delta_sd
  sig_threshold_lower <- -sig_level * delta_sd
  
  
  # Criteria
  crit_1 <- which(k_values >= k_0) 
  crit_2 <- which(delta_k_all >= sig_threshold_upper)
  crit_3 <- which(delta_k_all >= sig_threshold_lower) 
  
  # Approach 1: Metrics if k >= 0
  reintroduced_flags1 <- is_tp[crit_1]
  
  TP1 <- sum(reintroduced_flags1)
  FP1 <- sum(!reintroduced_flags1)
  
  precision1 <- ifelse((TP1 + FP1) > 0, TP1 / (TP1 + FP1), NA)
  
  
  
  recall1 <- ifelse(total_positives_discarded > 0, TP1 / total_positives_discarded, NA)
  F11 <- ifelse(is.na(precision1) | is.na(recall1) | (precision1 + recall1 == 0),
                NA,
                2 * (precision1 * recall1) / (precision1 + recall1))
  FDR1 <- ifelse((TP1 + FP1) > 0, FP1 / (TP1 + FP1), NA)
  
  
  metrics_1 <- data.frame(criterion = "k >= k_0",
                          TP = TP1, FP = FP1,
                          precision = precision1, recall = recall1,
                          FDR = FDR1, F1 = F11)
  
  
  
  # Approach 2: Metrics >= sig_upper (significant increase)
  reintroduced_flags2 <- is_tp[crit_2]
  TP2 <- sum(reintroduced_flags2)
  FP2 <- sum(!reintroduced_flags2)
  
  precision2 <- ifelse((TP2 + FP2) > 0, TP2 / (TP2 + FP2), NA)
  recall2 <- ifelse(total_positives_discarded > 0, TP2 / total_positives_discarded, NA)
  F12 <- ifelse(is.na(precision2) | is.na(recall2) | (precision2 + recall2 == 0),
                NA, 2 * (precision2 * recall2) / (precision2 + recall2))
  FDR2 <- ifelse((TP2 + FP2) > 0, FP2 / (TP2 + FP2), NA)
  
  
  
  metrics_2 <- data.frame(criterion = "k >= sig_upper",
                          TP = TP2, FP = FP2,
                          precision = precision2, recall = recall2,
                          FDR = FDR2, F1 = F12)
  
  
  #Approach 3: Metrics if k >= sig_lower
  reintroduced_flags3 <- is_tp[crit_3]
  TP3 <- sum(reintroduced_flags3)
  FP3 <- sum(!reintroduced_flags3)
  
  precision3 <- ifelse((TP3 + FP3) > 0, TP3 / (TP3 + FP3), NA)
  recall3 <- ifelse(total_positives_discarded > 0, TP3 / total_positives_discarded, NA)
  F13 <- ifelse(is.na(precision3) | is.na(recall3) | (precision3 + recall3 == 0),
                NA, 2 * (precision3 * recall3) / (precision3 + recall3))
  FDR3 <- ifelse((TP3 + FP3) > 0, FP3 / (TP3 + FP3), NA)
  
  
  metrics_3 <- data.frame(criterion = "k >= sig_lower",
                          TP = TP3, FP = FP3,
                          precision = precision3, recall = recall3,
                          FDR = FDR3, F1 = F13)
  
  
  metrics <- rbind(metrics_1,metrics_2, metrics_3)
  
  
  k_max <- max(k_values, na.rm = TRUE)
  delta_k <- k_max - k_0
  
  # Dataframe por ASV
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    delta_k = delta_k_all,
    BC_values = BC_values,
    J_values = J_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = k_max,
    delta_k = delta_k,
    recovered_asv_ap1 = added_asv[crit_1],
    recovered_asv_ap2 = added_asv[crit_2],
    recovered_asv_ap3 = added_asv[crit_3],
    thresholds = data.frame(upper = sig_threshold_upper, 
                            lower = sig_threshold_lower,
                            sd = delta_sd),
    filtered_values = data.frame(k = k_0, BC = BC_0, J = J_0),
    unfiltered_values = data.frame(k = k_unfiltered, BC = BC_unfiltered, J = J_unfiltered),
    global_metrics = metrics,
    time = time.take
  ))
}

k_validation_pairwaise9 <- function(OTU_table, ground_truth, alpha, beta, mode, sig_level = 1.96){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- calculate_k_w_single(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  # apply conservative filter
  filter_output <- conservative_filter9(asv)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- calculate_k_w_single(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- BC_values <- J_values <- c()
  true_positive_levels <- c() 
  
  
  blast_nuc_matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- calculate_k_w_single(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    j_i <- mean(na.omit(evenness_k(asv_table_i)))
    
    # Evaluate with ground_truth depending on the chosen mode
    
    match_level <- "No Match"
    
    if (mode == "taxonomy") {
      asv_species <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
      asv_genus   <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
      
      if (asv_species %in% ground_truth$species || asv_genus %in% ground_truth$genus) {
        match_level <- "Match"
      }
      
    } else if (mode == "sequence") {
      asv_identity <- seq_dict$Sequence[match(asv_name, seq_dict$ASV)]
      
      match_matrix <- sapply(ground_truth$Sequence, function(gt_seq) {
        
        alignment <- pwalign::pairwiseAlignment(as.character(asv_identity), as.character(gt_seq), type = "local",
                                                substitutionMatrix = blast_nuc_matrix, gapOpening =5,
                                                gapExtension = 2)
        
        if (score(alignment) <= 0) {
          return(c(pid = 0, qcov = 0)) # No alignment
        }
        
        percent_identity <- pwalign::pid(alignment, type = "PID1")
        
        # Query Coverage
        pattern_string <- as.character(alignment@pattern)
        pattern_no_gaps <- gsub("-", "", pattern_string)
        query_coverage <- (nchar(pattern_no_gaps) / nchar(asv_identity)) * 100
        
        return(c(pid = percent_identity, qcov = query_coverage))
      })
      
      
      if (is.vector(match_matrix)) {
        match_matrix <- t(as.matrix(match_matrix))
      }else {
        match_matrix <- t(match_matrix)
      }
      
      if (nrow(match_matrix) == 0) {
        best_pid <- 0; best_qcov <- 0
      } else {
        best_match_row <- match_matrix[order(match_matrix[, "qcov"],
                                             match_matrix[, "pid"],
                                             decreasing = TRUE)[1], ]
        best_pid <- best_match_row["pid"]
        best_qcov <- best_match_row["qcov"]
      }
      
      
      best_pid_rounded  <- round(best_pid, 2)
      best_qcov_rounded <- round(best_qcov, 2)
      
      # Classification
      if (best_pid >= 98 && best_qcov >= 100) {
        match_level <- "Match"
      }
    }
    
    
    added_asv <- c(added_asv, asv_name)
    k_values <- c(k_values, k_i)
    BC_values <- c(BC_values, bc_i)
    J_values <- c(J_values, j_i)
    true_positive_levels <- c(true_positive_levels, match_level)
    
  }
  
  is_tp <- true_positive_levels != "No Match"
  total_positives_discarded <- sum(is_tp)
  delta_k_all <- k_values - k_0
  
  # Calculate thresholds
  delta_sd <- sd(delta_k_all, na.rm = TRUE)
  sig_threshold_upper <- sig_level * delta_sd
  sig_threshold_lower <- -sig_level * delta_sd
  
  
  # Criteria
  crit_1 <- which(k_values >= k_0) 
  crit_2 <- which(delta_k_all >= sig_threshold_upper)
  crit_3 <- which(delta_k_all >= sig_threshold_lower) 
  
  # Approach 1: Metrics if k >= 0
  reintroduced_flags1 <- is_tp[crit_1]
  
  TP1 <- sum(reintroduced_flags1)
  FP1 <- sum(!reintroduced_flags1)
  
  precision1 <- ifelse((TP1 + FP1) > 0, TP1 / (TP1 + FP1), NA)
  
  
  
  recall1 <- ifelse(total_positives_discarded > 0, TP1 / total_positives_discarded, NA)
  F11 <- ifelse(is.na(precision1) | is.na(recall1) | (precision1 + recall1 == 0),
                NA,
                2 * (precision1 * recall1) / (precision1 + recall1))
  FDR1 <- ifelse((TP1 + FP1) > 0, FP1 / (TP1 + FP1), NA)
  
  
  metrics_1 <- data.frame(criterion = "k >= k_0",
                          TP = TP1, FP = FP1,
                          precision = precision1, recall = recall1,
                          FDR = FDR1, F1 = F11)
  
  
  
  # Approach 2: Metrics >= sig_upper (significant increase)
  reintroduced_flags2 <- is_tp[crit_2]
  TP2 <- sum(reintroduced_flags2)
  FP2 <- sum(!reintroduced_flags2)
  
  precision2 <- ifelse((TP2 + FP2) > 0, TP2 / (TP2 + FP2), NA)
  recall2 <- ifelse(total_positives_discarded > 0, TP2 / total_positives_discarded, NA)
  F12 <- ifelse(is.na(precision2) | is.na(recall2) | (precision2 + recall2 == 0),
                NA, 2 * (precision2 * recall2) / (precision2 + recall2))
  FDR2 <- ifelse((TP2 + FP2) > 0, FP2 / (TP2 + FP2), NA)
  
  
  
  metrics_2 <- data.frame(criterion = "k >= sig_upper",
                          TP = TP2, FP = FP2,
                          precision = precision2, recall = recall2,
                          FDR = FDR2, F1 = F12)
  
  
  #Approach 3: Metrics if k >= sig_lower
  reintroduced_flags3 <- is_tp[crit_3]
  TP3 <- sum(reintroduced_flags3)
  FP3 <- sum(!reintroduced_flags3)
  
  precision3 <- ifelse((TP3 + FP3) > 0, TP3 / (TP3 + FP3), NA)
  recall3 <- ifelse(total_positives_discarded > 0, TP3 / total_positives_discarded, NA)
  F13 <- ifelse(is.na(precision3) | is.na(recall3) | (precision3 + recall3 == 0),
                NA, 2 * (precision3 * recall3) / (precision3 + recall3))
  FDR3 <- ifelse((TP3 + FP3) > 0, FP3 / (TP3 + FP3), NA)
  
  
  metrics_3 <- data.frame(criterion = "k >= sig_lower",
                          TP = TP3, FP = FP3,
                          precision = precision3, recall = recall3,
                          FDR = FDR3, F1 = F13)
  
  
  metrics <- rbind(metrics_1,metrics_2, metrics_3)
  
  
  k_max <- max(k_values, na.rm = TRUE)
  delta_k <- k_max - k_0
  
  # Dataframe por ASV
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    delta_k = delta_k_all,
    BC_values = BC_values,
    J_values = J_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = k_max,
    delta_k = delta_k,
    recovered_asv_ap1 = added_asv[crit_1],
    recovered_asv_ap2 = added_asv[crit_2],
    recovered_asv_ap3 = added_asv[crit_3],
    thresholds = data.frame(upper = sig_threshold_upper, 
                            lower = sig_threshold_lower,
                            sd = delta_sd),
    filtered_values = data.frame(k = k_0, BC = BC_0, J = J_0),
    unfiltered_values = data.frame(k = k_unfiltered, BC = BC_unfiltered, J = J_unfiltered),
    global_metrics = metrics,
    time = time.take
  ))
}

k_validation_pairwaise10 <- function(OTU_table, ground_truth, alpha, beta, mode, sig_level = 1.96){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- calculate_k_w_single(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  # apply conservative filter
  filter_output <- conservative_filter10(asv)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- calculate_k_w_single(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- BC_values <- J_values <- c()
  true_positive_levels <- c() 
  
  
  blast_nuc_matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- calculate_k_w_single(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    j_i <- mean(na.omit(evenness_k(asv_table_i)))
    
    # Evaluate with ground_truth depending on the chosen mode
    
    match_level <- "No Match"
    
    if (mode == "taxonomy") {
      asv_species <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
      asv_genus   <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
      
      if (asv_species %in% ground_truth$species || asv_genus %in% ground_truth$genus) {
        match_level <- "Match"
      }
      
    } else if (mode == "sequence") {
      asv_identity <- seq_dict$Sequence[match(asv_name, seq_dict$ASV)]
      
      match_matrix <- sapply(ground_truth$Sequence, function(gt_seq) {
        
        alignment <- pwalign::pairwiseAlignment(as.character(asv_identity), as.character(gt_seq), type = "local",
                                                substitutionMatrix = blast_nuc_matrix, gapOpening =5,
                                                gapExtension = 2)
        
        if (score(alignment) <= 0) {
          return(c(pid = 0, qcov = 0)) # No alignment
        }
        
        percent_identity <- pwalign::pid(alignment, type = "PID1")
        
        # Query Coverage
        pattern_string <- as.character(alignment@pattern)
        pattern_no_gaps <- gsub("-", "", pattern_string)
        query_coverage <- (nchar(pattern_no_gaps) / nchar(asv_identity)) * 100
        
        return(c(pid = percent_identity, qcov = query_coverage))
      })
      
      
      if (is.vector(match_matrix)) {
        match_matrix <- t(as.matrix(match_matrix))
      }else {
        match_matrix <- t(match_matrix)
      }
      
      if (nrow(match_matrix) == 0) {
        best_pid <- 0; best_qcov <- 0
      } else {
        best_match_row <- match_matrix[order(match_matrix[, "qcov"],
                                             match_matrix[, "pid"],
                                             decreasing = TRUE)[1], ]
        best_pid <- best_match_row["pid"]
        best_qcov <- best_match_row["qcov"]
      }
      
      
      best_pid_rounded  <- round(best_pid, 2)
      best_qcov_rounded <- round(best_qcov, 2)
      
      # Classification
      if (best_pid >= 98 && best_qcov >= 100) {
        match_level <- "Match"
      }
    }
    
    
    added_asv <- c(added_asv, asv_name)
    k_values <- c(k_values, k_i)
    BC_values <- c(BC_values, bc_i)
    J_values <- c(J_values, j_i)
    true_positive_levels <- c(true_positive_levels, match_level)
    
  }
  
  is_tp <- true_positive_levels != "No Match"
  total_positives_discarded <- sum(is_tp)
  delta_k_all <- k_values - k_0
  
  # Calculate thresholds
  delta_sd <- sd(delta_k_all, na.rm = TRUE)
  sig_threshold_upper <- sig_level * delta_sd
  sig_threshold_lower <- -sig_level * delta_sd
  
  
  # Criteria
  crit_1 <- which(k_values >= k_0) 
  crit_2 <- which(delta_k_all >= sig_threshold_upper)
  crit_3 <- which(delta_k_all >= sig_threshold_lower) 
  
  # Approach 1: Metrics if k >= 0
  reintroduced_flags1 <- is_tp[crit_1]
  
  TP1 <- sum(reintroduced_flags1)
  FP1 <- sum(!reintroduced_flags1)
  
  precision1 <- ifelse((TP1 + FP1) > 0, TP1 / (TP1 + FP1), NA)
  
  
  
  recall1 <- ifelse(total_positives_discarded > 0, TP1 / total_positives_discarded, NA)
  F11 <- ifelse(is.na(precision1) | is.na(recall1) | (precision1 + recall1 == 0),
                NA,
                2 * (precision1 * recall1) / (precision1 + recall1))
  FDR1 <- ifelse((TP1 + FP1) > 0, FP1 / (TP1 + FP1), NA)
  
  
  metrics_1 <- data.frame(criterion = "k >= k_0",
                          TP = TP1, FP = FP1,
                          precision = precision1, recall = recall1,
                          FDR = FDR1, F1 = F11)
  
  
  
  # Approach 2: Metrics >= sig_upper (significant increase)
  reintroduced_flags2 <- is_tp[crit_2]
  TP2 <- sum(reintroduced_flags2)
  FP2 <- sum(!reintroduced_flags2)
  
  precision2 <- ifelse((TP2 + FP2) > 0, TP2 / (TP2 + FP2), NA)
  recall2 <- ifelse(total_positives_discarded > 0, TP2 / total_positives_discarded, NA)
  F12 <- ifelse(is.na(precision2) | is.na(recall2) | (precision2 + recall2 == 0),
                NA, 2 * (precision2 * recall2) / (precision2 + recall2))
  FDR2 <- ifelse((TP2 + FP2) > 0, FP2 / (TP2 + FP2), NA)
  
  
  
  metrics_2 <- data.frame(criterion = "k >= sig_upper",
                          TP = TP2, FP = FP2,
                          precision = precision2, recall = recall2,
                          FDR = FDR2, F1 = F12)
  
  
  #Approach 3: Metrics if k >= sig_lower
  reintroduced_flags3 <- is_tp[crit_3]
  TP3 <- sum(reintroduced_flags3)
  FP3 <- sum(!reintroduced_flags3)
  
  precision3 <- ifelse((TP3 + FP3) > 0, TP3 / (TP3 + FP3), NA)
  recall3 <- ifelse(total_positives_discarded > 0, TP3 / total_positives_discarded, NA)
  F13 <- ifelse(is.na(precision3) | is.na(recall3) | (precision3 + recall3 == 0),
                NA, 2 * (precision3 * recall3) / (precision3 + recall3))
  FDR3 <- ifelse((TP3 + FP3) > 0, FP3 / (TP3 + FP3), NA)
  
  
  metrics_3 <- data.frame(criterion = "k >= sig_lower",
                          TP = TP3, FP = FP3,
                          precision = precision3, recall = recall3,
                          FDR = FDR3, F1 = F13)
  
  
  metrics <- rbind(metrics_1,metrics_2, metrics_3)
  
  
  k_max <- max(k_values, na.rm = TRUE)
  delta_k <- k_max - k_0
  
  # Dataframe por ASV
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    delta_k = delta_k_all,
    BC_values = BC_values,
    J_values = J_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = k_max,
    delta_k = delta_k,
    recovered_asv_ap1 = added_asv[crit_1],
    recovered_asv_ap2 = added_asv[crit_2],
    recovered_asv_ap3 = added_asv[crit_3],
    thresholds = data.frame(upper = sig_threshold_upper, 
                            lower = sig_threshold_lower,
                            sd = delta_sd),
    filtered_values = data.frame(k = k_0, BC = BC_0, J = J_0),
    unfiltered_values = data.frame(k = k_unfiltered, BC = BC_unfiltered, J = J_unfiltered),
    global_metrics = metrics,
    time = time.take
  ))
}



#bootstrap
bootstrap_k <- function(OTU_table, asv_base, asv_candidate, 
                        alpha_weight, beta_weight, 
                        n_boot, alpha) {
  
  # Ensure there are no NAs in the data
  OTU_table[is.na(OTU_table)] <- 0
  asv_base[is.na(asv_base)] <- 0
  
  # Calculation of the initial k and with the candidate ASV
  k0 <- calculate_k_w_single(asv_base, alpha_weight, beta_weight)
  asv_new <- cbind(asv_base, OTU_table[, asv_candidate, drop = FALSE])
  k_i <- calculate_k_w_single(asv_new, alpha_weight, beta_weight)
  delta_k_obs <- k_i - k0
  
  # Bootstrap
  delta_boot <- replicate(n_boot, {
    idx <- sample(1:nrow(OTU_table), replace = TRUE)
    base_boot <- asv_base[idx, , drop = FALSE]
    cand_boot <- OTU_table[idx, asv_candidate, drop = FALSE]
    
    # replace possible NAs
    base_boot[is.na(base_boot)] <- 0
    cand_boot[is.na(cand_boot)] <- 0
    
    k0_b <- suppressWarnings(calculate_k_w_single(base_boot, alpha_weight, beta_weight))
    k_i_b <- suppressWarnings(calculate_k_w_single(cbind(base_boot, cand_boot), alpha_weight, beta_weight))
    
    k_i_b - k0_b
  })
  
  # Calculation of the p-value
  
  p_value <- (sum(delta_boot <= 0) + 1) / (length(delta_boot) + 1)
 
  
  signif <-  p_value < alpha
  
  list(delta_k_obs = delta_k_obs, p_value = p_value, signif = signif)
}


                    #      |------------------------------------|

bootstrap_k_dynamic <- function(OTU_table, asv_base, asv_candidate, 
                                alpha_weight, beta_weight, 
                                n_boot, alpha) {
  
  # Ensure there are no NAs in the data
  OTU_table[is.na(OTU_table)] <- 0
  asv_base[is.na(asv_base)] <- 0
  
  # Calculation of the initial k and with the candidate ASV
  k0 <- calculate_k_w_single(asv_base, alpha_weight, beta_weight)
  asv_new <- cbind(asv_base, OTU_table[, asv_candidate, drop = FALSE])
  k_i <- calculate_k_w_single(asv_new, alpha_weight, beta_weight)
  delta_k_obs <- k_i - k0
  
  # Bootstrap
  delta_boot <- replicate(n_boot, {
    # --- CORRECCIÓN: Usar nrow(asv_base) para el remuestreo ---
    idx <- sample(1:nrow(asv_base), replace = TRUE) 
    base_boot <- asv_base[idx, , drop = FALSE]
    cand_boot <- OTU_table[idx, asv_candidate, drop = FALSE]
    
    # replace possible NAs
    base_boot[is.na(base_boot)] <- 0
    cand_boot[is.na(cand_boot)] <- 0
    
    k0_b <- suppressWarnings(calculate_k_w_single(base_boot, alpha_weight, beta_weight))
    k_i_b <- suppressWarnings(calculate_k_w_single(cbind(base_boot, cand_boot), alpha_weight, beta_weight))
    
    k_i_b - k0_b
  })
  
  # --- CORRECCIÓN: Cálculo del p-value dinámico ---
  if (delta_k_obs > 0) {
    # Test if the INCREASE is significant (H0: delta <= 0)
    p_value <- (sum(delta_boot <= 0) + 1) / (n_boot + 1)
  } else if (delta_k_obs < 0) {
    # Test if the DECREASE is significant (H0: delta >= 0)
    p_value <- (sum(delta_boot >= 0) + 1) / (n_boot + 1)
  } else {
    # No change
    p_value <- 1.0
  }
  
  signif <-  p_value < alpha
  
  # Devuelve un data.frame para que sea fácil de combinar luego
  data.frame(
    added_asv = asv_candidate,
    delta_k_obs = delta_k_obs, 
    p_value = p_value, 
    signif = signif
  )
}


# Plots
plot_k2 <- function(results) {
  
  # Extract data
  df_results <- results$df_results
  k_inicial  <- results$k_inicial
  
  # Order from least to greatest
  df_plot_k <- df_results %>%
    arrange(k_values) 
  
  # As Factor
  df_plot_k$added_asv <- factor(df_plot_k$added_asv, levels = df_plot_k$added_asv)
  
  # Graph
  plot <- ggplot(df_plot_k, aes(x = added_asv, y = k_values)) +
    
    # Line to conect the points
    geom_line(aes(group = 1), color = "gray80", linetype = "solid", alpha = 0.5) +
    
    # Puntos coloreados por los 3 niveles de Ground Truth
    geom_point(aes(colour = TP_level), size = 2.5, alpha = 0.9) +
    
    # k-initial umbral line
    geom_hline(yintercept = k_inicial, linetype = "dashed", color = "red", linewidth = 0.7) +
    
    # red line tag
    annotate("text", x = 1, y = k_inicial, label = paste("k inicial:", round(k_inicial, 4)), 
             vjust = -1, hjust = 0, color = "red", size = 2.5, fontface = "bold") +
    
    # Titles and tags
    labs(title = "Evaluation of k for ALL Filtered ASVs",
         subtitle = "Comparison of re-evaluated ASVs against the initial community k",
         x = "ASV Name (Ordered by k value)",
         y = "New k value (k_i)",
         color = "Ground Truth Check") +
    
    # --- MODIFICACIÓN 2: Paleta de 3 colores y etiquetas ---
    scale_color_manual(
      # Asegúrate que los nombres coincidan con los de df_results
      values = c("Perfect Match" = "blue2", 
                 "Close Match"   = "skyblue", 
                 "No Match"      = "red2"),
      
      # Etiquetas para la leyenda
      labels = c("Perfect Match" = "Real Taxa (Perfect Match)", 
                 "Close Match"   = "Posible Real Taxa (Close Match)", 
                 "No Match"      = "Artifact (No Match)")) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom"
    )
  
  return(plot)
}

plot_k <- function(results) {
  
  # 1. Extraer los datos
  df_results <- results$df_results  # Contiene la columna "TP_level"
  k_inicial  <- results$k_inicial
  
  
  # Extraemos los umbrales heurísticos (calculados en k_validation)
  thresholds <- results$thresholds
  
  # Calculamos el valor k absoluto para las líneas del gráfico
  upper_threshold_line <- k_inicial + thresholds$upper
  lower_threshold_line <- k_inicial + thresholds$lower
  # --- FIN DE LA MODIFICACIÓN ---
  
  # 2. Procesar datos para el plot (SIN FILTRAR)
  df_plot_k <- df_results %>%
    arrange(k_values) # Ordenamos de menor a mayor k
  
  # 3. Convertir el nombre a FACTOR con el orden específico
  df_plot_k$added_asv <- factor(df_plot_k$added_asv, levels = df_plot_k$added_asv)
  
  # 4. Gráfico
  plot <- ggplot(df_plot_k, aes(x = added_asv, y = k_values)) +
    
    # Línea que conecta los puntos
    geom_line(aes(group = 1), color = "gray80", linetype = "solid", alpha = 0.5) +
    
    # Puntos coloreados por los 3 niveles de Ground Truth
    geom_point(aes(colour = TP_level), size = 2.5, alpha = 0.9) +
    
    # Línea de umbral k_inicial (Roja)
    geom_hline(yintercept = k_inicial, linetype = "dashed", color = "red", linewidth = 0.7) +
    
    # Líneas de umbral de significancia (Gris Oscuro)
    geom_hline(yintercept = upper_threshold_line, linetype = "dotted", color = "gray30", linewidth = 0.9) +
    geom_hline(yintercept = lower_threshold_line, linetype = "dotted", color = "gray30", linewidth = 0.9) +
    # --- FIN DE LA MODIFICACIÓN ---
    
    # Etiqueta para la línea roja
    annotate("text", x = 1, y = k_inicial, label = paste("k inicial:", round(k_inicial, 4)), 
             vjust = -1, hjust = 0, color = "red", size = 2.5, fontface = "bold") +
    
    
     annotate("text", x = 1, y = upper_threshold_line, label = "Umbral Sup. (+1.96 SD)", 
              vjust = -0.5, hjust = 0, color = "gray30", size = 2.5) +
     annotate("text", x = 1, y = lower_threshold_line, label = "Umbral Inf. (-1.96 SD)", 
              vjust = 1.5, hjust = 0, color = "gray30", size = 2.5) +
    
    # Títulos y Etiquetas actualizados
    labs(title = "Evaluation of k for ALL Filtered ASVs",
         subtitle = "Comparison of re-evaluated ASVs against the initial community k",
         x = "ASV Name (Ordered by k value)",
         y = "New k value (k_i)",
         color = "Ground Truth Check") +
    
    # Paleta de 3 colores y etiquetas
    scale_color_manual(
      values = c("Match"    = "blue2",
                 "No Match" = "red2"),    
      
      labels = c("Match"    = "True ASV",
                 "No Match" = "False ASV")) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom"
    )
  
  return(plot)
}


k_validation_pairwaise22 <- function(OTU_table, ground_truth, alpha, beta, mode){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- calculate_k_w_single(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  # apply conservative filter
  filter_output <- conservative_filter2(asv)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- calculate_k_w_single(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- BC_values <- J_values <- c()
  true_positive_levels <- c() 
  
  
  blast_nuc_matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- calculate_k_w_single(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    j_i <- mean(na.omit(evenness_k(asv_table_i)))
    
    # Evaluate with ground_truth depending on the chosen mode
    
    match_level <- "No Match"
    
    if (mode == "taxonomy") {
      asv_species <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
      asv_genus   <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
      
      if (asv_species %in% ground_truth$species) {
        match_level <- "Perfect Match" # Is our specie
      } else if (asv_genus %in% ground_truth$genus) {
        match_level <- "Close Match" # Genus Match
      }
      
    } else if (mode == "sequence") {
      asv_identity <- seq_dict$Sequence[match(asv_name, seq_dict$ASV)]
      
      match_matrix <- sapply(ground_truth$Sequence, function(gt_seq) {
        
        alignment <- pwalign::pairwiseAlignment(as.character(asv_identity), as.character(gt_seq), type = "local",
                                                substitutionMatrix = blast_nuc_matrix, gapOpening =5,
                                                gapExtension = 2)
        
        if (score(alignment) <= 0) {
          return(c(pid = 0, qcov = 0)) # No alignment
        }
        
        percent_identity <- pwalign::pid(alignment, type = "PID1")
        
        # Query Coverage
        pattern_string <- as.character(alignment@pattern)
        pattern_no_gaps <- gsub("-", "", pattern_string)
        query_coverage <- (nchar(pattern_no_gaps) / nchar(asv_identity)) * 100
        
        return(c(pid = percent_identity, qcov = query_coverage))
      })
      
      
      if (is.vector(match_matrix)) {
        match_matrix <- t(as.matrix(match_matrix))
      }else {
        match_matrix <- t(match_matrix)
      }
      
      if (nrow(match_matrix) == 0) {
        best_pid <- 0; best_qcov <- 0
      } else {
        best_match_row <- match_matrix[order(match_matrix[, "qcov"],
                                             match_matrix[, "pid"],
                                             decreasing = TRUE)[1], ]
        best_pid <- best_match_row["pid"]
        best_qcov <- best_match_row["qcov"]
      }
      
      
      best_pid_rounded  <- round(best_pid, 2)
      best_qcov_rounded <- round(best_qcov, 2)
      
      # Classification
      if (best_pid_rounded == 100 && best_qcov_rounded == 100) {
        match_level <- "Perfect Match"
        
      } else if (best_pid_rounded > 99.5 && best_qcov_rounded == 100) {
        match_level <- "Close Match"
      }
    }
    
    
    added_asv <- c(added_asv, asv_name)
    k_values <- c(k_values, k_i)
    BC_values <- c(BC_values, bc_i)
    J_values <- c(J_values, j_i)
    true_positive_levels <- c(true_positive_levels, match_level)
    
  }
  
  is_tp <- true_positive_levels != "No Match"
  total_positives_discarded <- sum(is_tp)
  
  
  # Approach 1: Metrics if k > 0
  reintroduced_asvs1 <- added_asv[k_values > k_0]
  reintroduced_flags1 <- is_tp[k_values > k_0]
  
  TP1 <- sum(reintroduced_flags1)
  FP1 <- sum(!reintroduced_flags1)
  
  precision1 <- ifelse((TP1 + FP1) > 0, TP1 / (TP1 + FP1), NA)
  
  
  
  recall1 <- ifelse(total_positives_discarded > 0, TP1 / total_positives_discarded, NA)
  F11 <- ifelse(is.na(precision1) | is.na(recall1) | (precision1 + recall1 == 0),
                NA,
                2 * (precision1 * recall1) / (precision1 + recall1))
  FDR1 <- ifelse((TP1 + FP1) > 0, FP1 / (TP1 + FP1), NA)
  
  metrics_1 <- data.frame(criterion = "k > k_0",
                          TP = TP1, FP = FP1,
                          precision = precision1, recall = recall1,
                          FDR = FDR1, F1 = F11)
  
  
  # Approach 2: Metrics if k >= 0
  reintroduced_asvs2 <- added_asv[k_values >= k_0]
  reintroduced_flags2 <- is_tp[k_values >= k_0]
  
  TP2 <- sum(reintroduced_flags2)
  FP2 <- sum(!reintroduced_flags2)
  
  precision2 <- ifelse((TP2 + FP2) > 0, TP2 / (TP2 + FP2), NA)
  
  
  
  recall2 <- ifelse(total_positives_discarded > 0, TP2 / total_positives_discarded, NA)
  F12 <- ifelse(is.na(precision2) | is.na(recall2) | (precision2 + recall2 == 0),
                NA,
                2 * (precision2 * recall2) / (precision2 + recall2))
  FDR2 <- ifelse((TP2 + FP2) > 0, FP2 / (TP2 + FP2), NA)
  
  metrics_2 <- data.frame(criterion = "k >= k_0",
                          TP = TP2, FP = FP2,
                          precision = precision2, recall = recall2,
                          FDR = FDR2, F1 = F12)
  
  metrics <- rbind(metrics_1,metrics_2)
  
  k_max <- max(k_values, na.rm = TRUE)
  
  delta_k <- k_max - k_0
  delta_k_rel <- ifelse(k_0 != 0, delta_k / k_0, NA)
  
  # Dataframe por ASV
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    BC_values = BC_values,
    J_values = J_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = k_max,
    delta_k = delta_k,
    delta_k_rel = delta_k_rel,
    asvs_recuperadas_ap1 = reintroduced_asvs1,
    asvs_recuperadas_ap2 = reintroduced_asvs2,
    filtered_values = data.frame(k = k_0, BC = BC_0, J = J_0),
    unfiltered_values = data.frame(k = k_unfiltered, BC = BC_unfiltered, J = J_unfiltered),
    global_metrics = metrics,
    time = time.take
  ))
}


k_validation_pairwaise_levels <- function(OTU_table, ground_truth, alpha, beta, mode, sig_level = 1.96){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- calculate_k_w_single(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  # apply conservative filter
  filter_output <- conservative_filter(asv)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- calculate_k_w_single(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- BC_values <- J_values <- c()
  true_positive_levels <- c() 
  
  
  blast_nuc_matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- calculate_k_w_single(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    j_i <- mean(na.omit(evenness_k(asv_table_i)))
    
    # Evaluate with ground_truth depending on the chosen mode
    
    match_level <- "No Match"
    
    if (mode == "taxonomy") {
      asv_species <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
      asv_genus   <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
      
      if (asv_species %in% ground_truth$species) {
        match_level <- "Perfect Match" # Is our specie
      } else if (asv_genus %in% ground_truth$genus) {
        match_level <- "Close Match" # Genus Match
      }
      
    } else if (mode == "sequence") {
      asv_identity <- seq_dict$Sequence[match(asv_name, seq_dict$ASV)]
      
      match_matrix <- sapply(ground_truth$Sequence, function(gt_seq) {
        
        alignment <- pwalign::pairwiseAlignment(as.character(asv_identity), as.character(gt_seq), type = "local",
                                                substitutionMatrix = blast_nuc_matrix, gapOpening =5,
                                                gapExtension = 2)
        
        if (score(alignment) <= 0) {
          return(c(pid = 0, qcov = 0)) # No alignment
        }
        
        percent_identity <- pwalign::pid(alignment, type = "PID1")
        
        # Query Coverage
        pattern_string <- as.character(alignment@pattern)
        pattern_no_gaps <- gsub("-", "", pattern_string)
        query_coverage <- (nchar(pattern_no_gaps) / nchar(asv_identity)) * 100
        
        return(c(pid = percent_identity, qcov = query_coverage))
      })
      
      
      if (is.vector(match_matrix)) {
        match_matrix <- t(as.matrix(match_matrix))
      }else {
        match_matrix <- t(match_matrix)
      }
      
      if (nrow(match_matrix) == 0) {
        best_pid <- 0; best_qcov <- 0
      } else {
        best_match_row <- match_matrix[order(match_matrix[, "qcov"],
                                             match_matrix[, "pid"],
                                             decreasing = TRUE)[1], ]
        best_pid <- best_match_row["pid"]
        best_qcov <- best_match_row["qcov"]
      }
      
      
      best_pid_rounded  <- round(best_pid, 2)
      best_qcov_rounded <- round(best_qcov, 2)
      
      # Classification
      if (best_pid_rounded == 100 && best_qcov_rounded == 100) {
        match_level <- "Perfect Match"
        
      } else if (best_pid_rounded > 99.5 && best_qcov_rounded == 100) {
        match_level <- "Close Match"
      }
    }
    
    
    added_asv <- c(added_asv, asv_name)
    k_values <- c(k_values, k_i)
    BC_values <- c(BC_values, bc_i)
    J_values <- c(J_values, j_i)
    true_positive_levels <- c(true_positive_levels, match_level)
    
  }
  
  is_tp <- true_positive_levels != "No Match"
  total_positives_discarded <- sum(is_tp)
  delta_k_all <- k_values - k_0
  
  # Calculate thresholds
  delta_sd <- sd(delta_k_all, na.rm = TRUE)
  sig_threshold_upper <- sig_level * delta_sd
  sig_threshold_lower <- -sig_level * delta_sd
  
  
  # Criteria
  crit_1 <- which(k_values >= k_0) 
  crit_2 <- which(delta_k_all >= sig_threshold_upper)
  crit_3 <- which(delta_k_all >= sig_threshold_lower) 
  
  # Approach 1: Metrics if k >= 0
  reintroduced_flags1 <- is_tp[crit_1]
  
  TP1 <- sum(reintroduced_flags1)
  FP1 <- sum(!reintroduced_flags1)
  
  precision1 <- ifelse((TP1 + FP1) > 0, TP1 / (TP1 + FP1), NA)
  
  
  
  recall1 <- ifelse(total_positives_discarded > 0, TP1 / total_positives_discarded, NA)
  F11 <- ifelse(is.na(precision1) | is.na(recall1) | (precision1 + recall1 == 0),
                NA,
                2 * (precision1 * recall1) / (precision1 + recall1))
  FDR1 <- ifelse((TP1 + FP1) > 0, FP1 / (TP1 + FP1), NA)
  
  
  metrics_1 <- data.frame(criterion = "k >= k_0",
                          TP = TP1, FP = FP1,
                          precision = precision1, recall = recall1,
                          FDR = FDR1, F1 = F11)
  
  
  
  # Approach 2: Metrics >= sig_upper (significant increase)
  reintroduced_flags2 <- is_tp[crit_2]
  TP2 <- sum(reintroduced_flags2)
  FP2 <- sum(!reintroduced_flags2)
  
  precision2 <- ifelse((TP2 + FP2) > 0, TP2 / (TP2 + FP2), NA)
  recall2 <- ifelse(total_positives_discarded > 0, TP2 / total_positives_discarded, NA)
  F12 <- ifelse(is.na(precision2) | is.na(recall2) | (precision2 + recall2 == 0),
                NA, 2 * (precision2 * recall2) / (precision2 + recall2))
  FDR2 <- ifelse((TP2 + FP2) > 0, FP2 / (TP2 + FP2), NA)
  
  
  
  metrics_2 <- data.frame(criterion = "k >= sig_upper",
                          TP = TP2, FP = FP2,
                          precision = precision2, recall = recall2,
                          FDR = FDR2, F1 = F12)
  
  
  #Approach 3: Metrics if k >= sig_lower
  reintroduced_flags3 <- is_tp[crit_3]
  TP3 <- sum(reintroduced_flags3)
  FP3 <- sum(!reintroduced_flags3)
  
  precision3 <- ifelse((TP3 + FP3) > 0, TP3 / (TP3 + FP3), NA)
  recall3 <- ifelse(total_positives_discarded > 0, TP3 / total_positives_discarded, NA)
  F13 <- ifelse(is.na(precision3) | is.na(recall3) | (precision3 + recall3 == 0),
                NA, 2 * (precision3 * recall3) / (precision3 + recall3))
  FDR3 <- ifelse((TP3 + FP3) > 0, FP3 / (TP3 + FP3), NA)
  
  
  metrics_3 <- data.frame(criterion = "k >= sig_lower",
                          TP = TP3, FP = FP3,
                          precision = precision3, recall = recall3,
                          FDR = FDR3, F1 = F13)
  
  
  metrics <- rbind(metrics_1,metrics_2, metrics_3)
  
  
  k_max <- max(k_values, na.rm = TRUE)
  delta_k <- k_max - k_0
  
  # Dataframe por ASV
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    delta_k = delta_k_all,
    BC_values = BC_values,
    J_values = J_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = k_max,
    delta_k = delta_k,
    recovered_asv_ap1 = added_asv[crit_1],
    recovered_asv_ap2 = added_asv[crit_2],
    recovered_asv_ap3 = added_asv[crit_3],
    thresholds = data.frame(upper = sig_threshold_upper, 
                            lower = sig_threshold_lower,
                            sd = delta_sd),
    filtered_values = data.frame(k = k_0, BC = BC_0, J = J_0),
    unfiltered_values = data.frame(k = k_unfiltered, BC = BC_unfiltered, J = J_unfiltered),
    global_metrics = metrics,
    time = time.take
  ))
}

plot_k_levels <- function(results) {
  
  # 1. Extraer los datos
  df_results <- results$df_results  # Contiene la columna "TP_level"
  k_inicial  <- results$k_inicial
  
  # --- INICIO DE LA MODIFICACIÓN ---
  # Extraemos los umbrales heurísticos (calculados en k_validation)
  thresholds <- results$thresholds
  
  # Calculamos el valor k absoluto para las líneas del gráfico
  upper_threshold_line <- k_inicial + thresholds$upper
  lower_threshold_line <- k_inicial + thresholds$lower
  # --- FIN DE LA MODIFICACIÓN ---
  
  # 2. Procesar datos para el plot (SIN FILTRAR)
  df_plot_k <- df_results %>%
    arrange(k_values) # Ordenamos de menor a mayor k
  
  # 3. Convertir el nombre a FACTOR con el orden específico
  df_plot_k$added_asv <- factor(df_plot_k$added_asv, levels = df_plot_k$added_asv)
  
  # 4. Gráfico
  plot <- ggplot(df_plot_k, aes(x = added_asv, y = k_values)) +
    
    # Línea que conecta los puntos
    geom_line(aes(group = 1), color = "gray80", linetype = "solid", alpha = 0.5) +
    
    # Puntos coloreados por los 3 niveles de Ground Truth
    geom_point(aes(colour = TP_level), size = 2.5, alpha = 0.9) +
    
    # Línea de umbral k_inicial (Roja)
    geom_hline(yintercept = k_inicial, linetype = "dashed", color = "red", linewidth = 0.7) +
    
    # Líneas de umbral de significancia (Gris Oscuro)
    geom_hline(yintercept = upper_threshold_line, linetype = "dotted", color = "gray30", linewidth = 0.9) +
    geom_hline(yintercept = lower_threshold_line, linetype = "dotted", color = "gray30", linewidth = 0.9) +
    # --- FIN DE LA MODIFICACIÓN ---
    
    # Etiqueta para la línea roja
    annotate("text", x = 1, y = k_inicial, label = paste("k inicial:", round(k_inicial, 4)), 
             vjust = -1, hjust = 0, color = "red", size = 2.5, fontface = "bold") +
    
    
    annotate("text", x = 1, y = upper_threshold_line, label = "Umbral Sup. (+1.96 SD)", 
             vjust = -0.5, hjust = 0, color = "gray30", size = 2.5) +
    annotate("text", x = 1, y = lower_threshold_line, label = "Umbral Inf. (-1.96 SD)", 
             vjust = 1.5, hjust = 0, color = "gray30", size = 2.5) +
    
    # Títulos y Etiquetas actualizados
    labs(title = "Evaluation of k for ALL Filtered ASVs",
         subtitle = "Comparison of re-evaluated ASVs against the initial community k",
         x = "ASV Name (Ordered by k value)",
         y = "New k value (k_i)",
         color = "Ground Truth Check") +
    
    # Paleta de 3 colores y etiquetas
    scale_color_manual(
      values = c("Perfect Match" = "blue2", 
                 "Close Match"   = "skyblue", 
                 "No Match"      = "red2"),
      
      labels = c("Perfect Match" = "Real Taxa (Perfect Match)", 
                 "Close Match"   = "Posible Real Taxa (Close Match)", 
                 "No Match"      = "Artifact (No Match)")) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom"
    )
  
  return(plot)
}

plot_k2 <- function(results) {
  
  # 1. Extraer los datos
  df_results <- results$df_results  # Contiene la columna "TP_level"
  k_inicial  <- results$k_inicial
  
  # --- INICIO DE LA MODIFICACIÓN ---
  # Extraemos los umbrales heurísticos (calculados en k_validation)
  thresholds <- results$thresholds
  
  # Calculamos el valor k absoluto para las líneas del gráfico
  upper_threshold_line <- k_inicial + thresholds$upper
  lower_threshold_line <- k_inicial + thresholds$lower
  # --- FIN DE LA MODIFICACIÓN ---
  
  # 2. Procesar datos para el plot (SIN FILTRAR)
  df_plot_k <- df_results %>%
    arrange(k_values) # Ordenamos de menor a mayor k
  
  # 3. Convertir el nombre a FACTOR con el orden específico
  df_plot_k$added_asv <- factor(df_plot_k$added_asv, levels = df_plot_k$added_asv)
  
  # 4. Gráfico
  plot <- ggplot(df_plot_k, aes(x = added_asv, y = k_values)) +
    
    # Línea que conecta los puntos
    geom_line(aes(group = 1), color = "gray80", linetype = "solid", alpha = 0.5) +
    
    # Puntos coloreados por los 2 niveles de Ground Truth
    geom_point(aes(colour = TP_level), size = 2.5, alpha = 0.9) +
    
    # Línea de umbral k_inicial (Roja)
    geom_hline(yintercept = k_inicial, linetype = "dashed", color = "red", linewidth = 0.7) +
    
    # Líneas de umbral de significancia (Gris Oscuro)
    #geom_hline(yintercept = upper_threshold_line, linetype = "dotted", color = "gray30", linewidth = 0.9) +
    #geom_hline(yintercept = lower_threshold_line, linetype = "dotted", color = "gray30", linewidth = 0.9) +
    
    # Etiqueta para la línea roja
    annotate("text", x = 1, y = k_inicial, label = paste("k inicial:", round(k_inicial, 4)), 
             vjust = -1, hjust = 0, color = "red", size = 2.5, fontface = "bold") +
    
    
    #annotate("text", x = 1, y = upper_threshold_line, label = "Umbral Sup. (+1.96 SD)", 
             #vjust = -0.5, hjust = 0, color = "gray30", size = 2.5) +
    #annotate("text", x = 1, y = lower_threshold_line, label = "Umbral Inf. (-1.96 SD)", 
             #vjust = 1.5, hjust = 0, color = "gray30", size = 2.5) +
    
    # Títulos y Etiquetas actualizados
    labs(title = "Evaluation of k for ALL Filtered ASVs",
         subtitle = "Comparison of re-evaluated ASVs against the initial community k",
         x = "ASV Name (Ordered by k value)",
         y = "New k value (k_i)",
         color = "Ground Truth Check") +
    
    # Paleta de 3 colores y etiquetas
    scale_color_manual(
      values = c("Match"    = "blue2",
                 "No Match" = "red2"),    
      
      labels = c("Match"    = "True ASV",
                 "No Match" = "False ASV")) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom"
    )
  
  return(plot)
}

