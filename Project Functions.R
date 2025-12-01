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

chao1 <- function(OTU_table){
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Calculate with vegan
  est <- vegan::estimateR(asv)
  
  #Extract only Chao1
  return(est["S.chao1", ])
}

jaccard <- function(OTU_table){
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Calculate with vegan
  jaccard_dist <- vegan::vegdist(asv, method = "jaccard", binary = TRUE)
  
  return(jaccard_dist)
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

k3 <- function(OTU_table, alpha, beta) {
  
  asv_num <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Check if there is enough data to calculate
  if (ncol(asv_num) < 2 || nrow(asv_num) < 2) {
    return(0)
  }
  
  #Calculation of chao1 and normalization
  chao_vals <- chao1(asv_num)
  chao_mean <- mean(chao_vals, na.rm = TRUE)
  total_asvs <- ncol(asv_num)
  chao_norm <- ifelse(total_asvs > 0, chao_mean / total_asvs, 0)
  
  # Normalize using min-max
  #normalize_01 <- function(x){
    #if (length(unique(x)) == 1) return(rep(0, length(x)))
    #(x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  #}
  #chao_norm <- normalize_01(chao_vals)
  #chao_mean_norm <- mean(chao_norm, na.rm = TRUE)
  
  #Bray Curtis
  BC_mean <- mean(na.omit(bray_curtis(asv_num)))
  
  #k3 <- (chao_mean_norm^alpha) * (BC_mean^beta)
  k3 <- (chao_norm^alpha) * (BC_mean^beta)
  return(k3)
}

k4 <- function(OTU_table, alpha, beta) {
  
  asv_num <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Check if there is enough data to calculate
  if (ncol(asv_num) < 2 || nrow(asv_num) < 2) {
    return(0)
  }
  
  J_mean <- mean(na.omit(evenness_k(asv_num)))
  jaccard_mean <- mean(jaccard(asv_num))
  
  if (is.nan(J_mean) || is.nan(jaccard_mean)) {
    return(0)
  }
  
  k4 <- (J_mean^alpha) * (jaccard_mean^beta)
  return(k4)
}

# Conservative filter standard
conservative_filter <- function(OTU_table, col_max = 10,pct_samples_threshold = 0.05,
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
  C1_delete <- (col_maxs <= col_max) & (pct_samples_present <= pct_samples_threshold)
  
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

k_validation_pairwaise_k2 <- function(OTU_table, ground_truth, alpha, beta, mode, filter_parameters){
  
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
  filter_arguments <- c(list(OTU_table=asv), filter_parameters)
  filter_output <- do.call(conservative_filter, filter_arguments)
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
  
  # Only in sequence mode
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
  
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    delta_k = k_values - k_0,
    BC_values = BC_values,
    J_values = J_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = max(k_values),
    time = round(end.time - start.time, 2)
  ))
}

k_validation_pairwaise_k3 <- function(OTU_table, ground_truth, alpha, beta, mode, filter_parameters){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- k3(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  Chao1_unfiltered <- mean(na.omit(chao1(asv)))
  
  # apply conservative filter
  filter_arguments <- c(list(OTU_table=asv), filter_parameters)
  filter_output <- do.call(conservative_filter, filter_arguments)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- k3(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  Chao1_0 <- mean(na.omit(chao1(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- BC_values <- Chao1_values <- c()
  true_positive_levels <- c() 
  
  # Only in sequence mode
  blast_nuc_matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- k3(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    chao1_i <- mean(na.omit(chao1(asv_table_i)))
    
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
    Chao1_values <- c(Chao1_values, chao1_i)
    true_positive_levels <- c(true_positive_levels, match_level)
    
  }
  
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    delta_k = k_values - k_0,
    BC_values = BC_values,
    Chao1_values = Chao1_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = max(k_values),
    time = round(end.time - start.time, 2)
  ))
}

k_validation_pairwaise_k4 <- function(OTU_table, ground_truth, alpha, beta, mode, filter_parameters){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- k4(asv, alpha, beta)
  JC_unfiltered <- mean(na.omit(jaccard(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  # apply conservative filter
  filter_arguments <- c(list(OTU_table=asv), filter_parameters)
  filter_output <- do.call(conservative_filter, filter_arguments)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- k4(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(jaccard(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- JC_values <- J_values <- c()
  true_positive_levels <- c() 
  
  # Only in sequence mode
  blast_nuc_matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- k4(asv_table_i, alpha, beta)
    jc_i <- mean(na.omit(jaccard(asv_table_i)))
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
    JC_values <- c(JC_values, jc_i)
    J_values <- c(J_values, j_i)
    true_positive_levels <- c(true_positive_levels, match_level)
    
  }
  
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    delta_k = k_values - k_0,
    JC_values = JC_values,
    J_values = J_values,
    TP_level = true_positive_levels
  )
  
  end.time <- Sys.time()
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = max(k_values),
    time = round(end.time - start.time, 2)
  ))
}



# Calculate metrics
calculate_threshold_metrics <- function(validation_output) {
  
  # Extract data
  df <- validation_output$df_results
  k_initial <- validation_output$k_initial
  
  # Match is True
  df$is_true <- df$TP_level == "Match"
  
  # Total True ASVs
  total_true <- sum(df$is_true)
  total_false <- sum(!df$is_true)
  
  # Standard desviation
  sd_noise <- sd(df$delta_k, na.rm = TRUE)
  
  # Thresholds
  # Try with all thresholds for each k value
  thresholds <- sort(unique(c(df$k_values, k_initial)))
  
  # Calculate metrics for each threshold
  metrics_list <- lapply(thresholds, function(threshold){
    
    # Criterion: k >= threshold
    accepted_subset <- df[df$k_values >= threshold, ]
    
    # Count
    TP <- sum(accepted_subset$is_true)
    FP <- sum(!accepted_subset$is_true)
    
    # Metrics
    # Precision
    precision <- ifelse((TP + FP) > 0, TP/(TP+FP), 0)
    
    #Recall
    recall <- ifelse(total_true > 0, TP/ total_true, 0)
    
    #False positive rate
    # False rescued from all falses
    fpr<- ifelse(total_false >0, FP / total_false, 0)
    
    #F1- score
    f1 <- ifelse((precision + recall) > 0,
                 2* (precision * recall) / (precision + recall), 0)
    
    #FDR (False Discovery Rate)
    fdr <- ifelse((TP + FP) > 0, FP/ (TP+FP),0)
    
    # Ratio True/False
    ratio_tf <- ifelse(FP>0, TP/FP, Inf)
    
    #z-score above the current threshold
    z_score <- (threshold - k_initial) / sd_noise
    
    
    return(data.frame(
      Threshold_k = threshold,
      TP = TP,
      FP = FP,
      Total_True = total_true,
      Total_False = total_false,
      Ratio_TF =ratio_tf,
      Precision = round(precision,3),
      Recall = round(recall, 3),
      F1 = round(f1, 3),
      FDR = round(fdr, 3),
      FPR = round(fpr,3),
      Z_score = round(z_score, 2)
    ))
  })
  
  metrics <- do.call(rbind, metrics_list)
  
  return(metrics)
}



plot_k <- function(results) {
  
  # Extract data
  df_results <- results$df_results  # Contiene la columna "TP_level"
  k_inicial  <- results$k_inicial
  
  # data
  df_plot_k <- df_results %>%
    arrange(k_values) # Order
  
  # As factor
  df_plot_k$added_asv <- factor(df_plot_k$added_asv, levels = df_plot_k$added_asv)
  
  # Graph
  plot <- ggplot(df_plot_k, aes(x = added_asv, y = k_values)) +
    
    # Line to conect points
    geom_line(aes(group = 1), color = "gray80", linetype = "solid", alpha = 0.5) +
    
    # Colour for three levels
    geom_point(aes(colour = TP_level), size = 2.5, alpha = 0.9) +
    
    # k initial line
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
    
    # Colurs and tags
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

