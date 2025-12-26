library(Biostrings)
library(pwalign)
library(BiocGenerics)
library(tidyverse)
library(vegan)
library(dplyr)
library(purrr)
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
  J <- ifelse(S > 1, shannon_index / log(S), 0) #Calculating J 
  
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

#Function of k2 approach
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

# Function of k4 approach
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

# Conservative filter standard. The filter combine three conventional filters
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

# k2: Pielou - Bray Curtis
k2_optimization <- function(OTU_table, ground_truth, alpha, beta, mode, filter_parameters){
  
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

# Plot for k2
plot_k2 <- function(results) {
  
  # Extract data
  df_results <- results$df_results
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
    labs(title = "Evaluation of k values for all filtered ASVs in k2 algorithm",
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

#k3: Chao1 - Bray Curtis
k3_optimization <- function(OTU_table, ground_truth, alpha, beta, mode, filter_parameters){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessing and initial filter
  asv <- OTU_table %>% mutate_if(is.integer, as.numeric) %>% select(where(is.numeric))
  
  filter_arguments <- c(list(OTU_table=asv), filter_parameters)
  filter_output <- do.call(conservative_filter, filter_arguments)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  discarded_asv <- colnames(asv_deleted)
  
  # Calculate base line
  # Raw values of the initial community (preserved)
  chao_0_raw <- mean(na.omit(chao1(asv_preserved)))
  bc_0_raw   <- mean(na.omit(bray_curtis(asv_preserved)))
  
  # Initialize lists
  # Baseline in the position 1
  raw_chao_list <- c(chao_0_raw)
  raw_bc_list   <- c(bc_0_raw)
  
  # Lists for results
  added_asv_list <- c("BASELINE") 
  tp_levels_list <- c("Base")
  
  
  # Loop
  for (asv_name in discarded_asv) {
    
    # Temporal table
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    # Calculate Chao and Bray 
    ch_i <- mean(na.omit(chao1(asv_table_i)))
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    
    # Ground Truth Check
    match_level <- "No Match"
    if (mode == "taxonomy") {
      if(exists("seq_dict")){
        asv_sp <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
        asv_gn <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
        if (asv_sp %in% ground_truth$species || asv_gn %in% ground_truth$genus) match_level <- "Match"
      }
    } else if (mode == "sequence") {
      # No need in this case. It is not evaluated in this mode
    }
    
    # Save in the lists
    raw_chao_list <- c(raw_chao_list, ch_i)
    raw_bc_list   <- c(raw_bc_list, bc_i)
    added_asv_list <- c(added_asv_list, asv_name)
    tp_levels_list <- c(tp_levels_list, match_level)
  }
  
  
  # Calculate gain relative to the baseline
  # For the first position the result will be 0
  chao_gain_list <- raw_chao_list - raw_chao_list[1]
  
  # Normalize the gains
  min_gain <- min(chao_gain_list, na.rm=TRUE)
  max_gain <- max(chao_gain_list, na.rm=TRUE)
  
  # Division by zero protection
  if (max_gain == min_gain) {
    chao_norm <- rep(0, length(chao_gain_list)) 
  } else {
    chao_norm <- (chao_gain_list - min_gain) / (max_gain - min_gain)
  }
  
  # Bray Curtis
  bc_factor <- raw_bc_list
  
  # Calculate k
  # Apply weights
  k_final_vector <- (chao_norm ^ alpha) * (bc_factor ^ beta)
  
  # Separate results
  
  k_0_final <- k_final_vector[1]
  k_values_loop <- k_final_vector[-1]
  
  # Construct final dataframe
  df_results <- data.frame(
    added_asv = added_asv_list[-1],
    k_values = k_values_loop,
    
    # Delta k with respect to the beginning
    delta_k = k_values_loop - k_0_final, 
    
    Chao1_Base = rep(raw_chao_list[1], length(k_values_loop)),
    # Metrics
    Chao1_Raw = raw_chao_list[-1],    
    Chao1_Gain = chao_gain_list[-1],        
    Chao1_Norm_Gain = chao_norm[-1],        
    
    BC_Raw = raw_bc_list[-1],
    TP_level = tp_levels_list[-1]
  )
  
  return(list(
    df_results = df_results,
    k_inicial = k_0_final,
    k_max = max(k_values_loop, na.rm=TRUE),
    time = round(Sys.time() - start.time, 2)
  ))
}

# Plot for k3
plot_k3 <- function(results) {
  
  # Extract data
  df_results <- results$df_results
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
    labs(title = "Evaluation of k values for all filtered ASVs in k3 algorithm",
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

# k4: Pielou- Jaccard
k4_optimization <- function(OTU_table, ground_truth, alpha, beta, mode, filter_parameters){
  
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

# Plot for k4
plot_k4 <- function(results) {
  
  # Extract data
  df_results <- results$df_results
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
    labs(title = "Evaluation of k values for all filtered ASVs in k4 algorithm",
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

#Chao1 - Jaccard

k5_optimization <- function(OTU_table, ground_truth, alpha, beta, mode, filter_parameters){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Initial filter and preprocessed
  asv <- OTU_table %>% mutate_if(is.integer, as.numeric) %>% select(where(is.numeric))
  
  filter_arguments <- c(list(OTU_table=asv), filter_parameters)
  filter_output <- do.call(conservative_filter, filter_arguments)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  discarded_asv <- colnames(asv_deleted)
  
  # Calculate initial values
  chao_0_raw <- mean(na.omit(chao1(asv_preserved)))
  jc_0_raw   <- mean(na.omit(jaccard(asv_preserved)))
  
  raw_chao_list <- c(chao_0_raw)
  raw_jc_list   <- c(jc_0_raw)
  
  # Result lists
  added_asv_list <- c("BASELINE") 
  tp_levels_list <- c("Base")
  
  # Loop
  for (asv_name in discarded_asv) {
    
    # Temporal table (Preserved + ASVi)
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    # Calculate metrics
    ch_i <- mean(na.omit(chao1(asv_table_i)))
    jc_i <- mean(na.omit(jaccard(asv_table_i)))
    
    # Ground Truth Check
    match_level <- "No Match"
    if (mode == "taxonomy") {
      # Asegúrate de que seq_dict esté disponible en el entorno global o pásalo como argumento
      if(exists("seq_dict")){
        asv_sp <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
        asv_gn <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
        if (asv_sp %in% ground_truth$species || asv_gn %in% ground_truth$genus) match_level <- "Match"
      }
    } else if (mode == "sequence") {
      # No need this approach
    }
    
    # Save in the lists
    raw_chao_list <- c(raw_chao_list, ch_i)
    raw_jc_list   <- c(raw_jc_list, jc_i)
    added_asv_list <- c(added_asv_list, asv_name)
    tp_levels_list <- c(tp_levels_list, match_level)
  }
  
  # Normalization
  
  # A. Calcular GANANCIA (Gain) relativa al Baseline
  #Calculate "gain" to baseline
  chao_gain_list <- raw_chao_list - raw_chao_list[1]
  
  # Normalization of the gains (min - max)
  min_gain <- min(chao_gain_list, na.rm=TRUE)
  max_gain <- max(chao_gain_list, na.rm=TRUE)
  
  # Error control (/0)
  if (max_gain == min_gain) {
    chao_norm <- rep(0, length(chao_gain_list)) 
  } else {
    chao_norm <- (chao_gain_list - min_gain) / (max_gain - min_gain)
  }
  
  # Jaccard no need normalization
  jc_factor <- raw_jc_list
  
  # K calculator
  # Apply weights
  k_final_vector <- (chao_norm ^ alpha) * (jc_factor ^ beta)
  
  # Results
  k_0_final <- k_final_vector[1]
  k_values_loop <- k_final_vector[-1]
  
  # Final dataframe
  df_results <- data.frame(
    added_asv = added_asv_list[-1],
    k_values = k_values_loop,
    
    # Delta k 
    delta_k = k_values_loop - k_0_final, 
    
    Chao1_Base = rep(raw_chao_list[1], length(k_values_loop)),
    Chao1_Raw = raw_chao_list[-1],          
    Chao1_Gain = chao_gain_list[-1],        
    Chao1_Norm_Gain = chao_norm[-1],        
    
    JC_Raw = raw_jc_list[-1],
    TP_level = tp_levels_list[-1]
  )
  
  return(list(
    df_results = df_results,
    k_inicial = k_0_final,
    k_max = max(k_values_loop, na.rm=TRUE),
    time = round(Sys.time() - start.time, 2)
  ))
}

# Plot for k5
plot_k5 <- function(results) {
  
  # Extract data
  df_results <- results$df_results
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
    labs(title = "Evaluation of k values for all filtered ASVs in k5 algorithm",
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

# Pielou - Aitchison

k6_optimization <- function(OTU_table, ground_truth, alpha, beta, mode, filter_parameters){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))

  
  aitchison_dist <- function(vec_base, vec_cand) {
    # Numeric values
    vec_base <- as.numeric(vec_base)
    vec_cand <- as.numeric(vec_cand)
    
    # Vectors need to have the same length
    # Sum pseudocount
    mat <- rbind(vec_base, vec_cand) + 0.5
    
    # Geometric mean and CLR
    gm <- apply(mat, 1, function(x) exp(mean(log(x))))
    clr_mat <- log(mat / gm)
    
    # Euclidean distance
    return(as.numeric(dist(clr_mat, method = "euclidean")))
  }

  
  asv <- as.data.frame(OTU_table) %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  filter_arguments <- c(list(OTU_table=asv), filter_parameters)
  filter_output <- do.call(conservative_filter, filter_arguments)
  
  asv_preserved <- as.data.frame(filter_output$preserved)
  asv_deleted <- filter_output$deleted
  discarded_asv <- colnames(asv_deleted)
  
  #Initial values
  
  val_J0 <- evenness_k(asv_preserved)
  J_0_raw <- mean(as.numeric(unlist(val_J0)), na.rm = TRUE)
  Aitch_0_raw <- 0
  
  raw_J_list     <- c(J_0_raw)
  raw_Aitch_list <- c(Aitch_0_raw)
  added_asv_list <- c("BASELINE")
  tp_levels_list <- c("Base")
  
  # Common vector
  vec_common <- as.numeric(unlist(asv_preserved))
  
  # ADD-ONE LOOP
  
  for (asv_name in discarded_asv) {
    
    # Take the candidate column
    asv_cand_col <- asv_deleted[, asv_name, drop = FALSE] 
    
    # For Pielou
    asv_table_i <- cbind(asv_preserved, asv_cand_col)
    
    # Pielou
    val_j <- evenness_k(asv_table_i)
    j_i <- mean(as.numeric(unlist(val_j)), na.rm=TRUE)
    
    # Aitchison
    val_new_asv <- as.numeric(unlist(asv_cand_col))
    
    # Repeat 0 to equal the lenght of the new column
    zeros_to_add <- rep(0, length(val_new_asv))
    
    vec_base_aligned <- c(vec_common, zeros_to_add)
    vec_cand_aligned <- c(vec_common, val_new_asv)
    
    aitch_i <- aitchison_dist(vec_base_aligned, vec_cand_aligned)
    
    # Taxonomy Check
    match_level <- "No Match"
    if (mode == "taxonomy") {
      if(exists("seq_dict")){ 
        asv_species <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
        asv_genus   <- seq_dict$Genus[match(asv_name, seq_dict$ASV)]
        
        # Error control (NA)
        if (!is.na(asv_species) && !is.na(asv_genus)) {
          if (asv_species %in% ground_truth$species || asv_genus %in% ground_truth$genus) {
            match_level <- "Match"
          }
        }
      }
    }
    
    added_asv_list <- c(added_asv_list, asv_name)
    raw_J_list     <- c(raw_J_list, j_i)
    raw_Aitch_list <- c(raw_Aitch_list, aitch_i)
    tp_levels_list <- c(tp_levels_list, match_level)
  } 

  # Aitchison normalization (only candidates)
  candidates_aitch <- raw_Aitch_list[-1]
  
  if(length(candidates_aitch) > 0){
    min_cand <- min(candidates_aitch, na.rm=TRUE)
    max_cand <- max(candidates_aitch, na.rm=TRUE)
    
    if (max_cand == min_cand) {
      norm_candidates <- rep(0, length(candidates_aitch))
    } else {
      norm_candidates <- (candidates_aitch - min_cand) / (max_cand - min_cand)
    }
    Aitch_norm_ZOOM <- c(0, norm_candidates)
  } else {
    Aitch_norm_ZOOM <- c(0)
  }
  
  # Calculate k final
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)
  
  k_final_vector <- (raw_J_list ^ alpha) * (Aitch_norm_ZOOM ^ beta)
  
  # Results
  
  if(length(k_final_vector) > 1){
    k_0 <- k_final_vector[1]        
    k_values <- k_final_vector[-1] 
    
    df_results <- data.frame(
      added_asv = added_asv_list[-1],
      k_values  = k_values,
      delta_k   = k_values - k_0, 
      
      J_Raw            = raw_J_list[-1],
      Aitchison_Raw    = raw_Aitch_list[-1],
      Aitchison_Norm   = Aitch_norm_ZOOM[-1], 
      
      TP_level = tp_levels_list[-1]
    )
  } else {
    k_0 <- k_final_vector[1]
    df_results <- data.frame() 
  }
  
  end.time <- Sys.time()
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = if(length(df_results) > 0) max(df_results$k_values, na.rm=TRUE) else k_0,
    time = round(end.time - start.time, 2)
  ))
}

# Plot for k6
plot_k6 <- function(results) {
  
  # Extract data
  df_results <- results$df_results
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
    labs(title = "Evaluation of k values for all filtered ASVs in k6 algorithm",
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


# Calculate metrics
calculate_threshold_metrics <- function(validation_output, penalty_weight = 1) {
  
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
    
    FN <- total_true - TP
    TN <- total_false - FP
    
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
    
    #F0.5 Score
    beta_val <- 0.5
    beta_sq <- beta_val^2
    
    f0_5 <- ifelse((beta_sq*precision)+recall > 0,
                   ((1 + beta_sq) * precision * recall) / ((beta_sq*precision)+recall),
                   0)
    
    #F1-FPR
    f1_fpr <- f1 - (penalty_weight * fpr)
    
    #F1/FPR
    f1_over_fpr <- f1/(1 + penalty_weight*fpr)
    
    #f balanced
    f_balanced <- 3/((1/precision) + (1/recall) + (1/(1-fpr)))
    
    # Matthews correlation coefficient (MCC)
    mcc = ((TP*TN) - (FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    
    return(data.frame(
      Threshold_k = threshold,
      TP = TP,
      FP = FP,
      FN =FN,
      TN = TN,
      Total_True = total_true,
      Total_False = total_false,
      Ratio_TF =ratio_tf,
      Precision = round(precision,3),
      Recall = round(recall, 3),
      F1 = round(f1, 3),
      FDR = round(fdr, 3),
      FPR = round(fpr,3),
      Z_score = round(z_score, 2),
      F0.5 = f0_5,
      "F1-FPR" = f1_fpr,
      "F1/FPR" = f1_over_fpr,
      F_Balanced = f_balanced,
      MCC = mcc
    ))
  })
  
  metrics <- do.call(rbind, metrics_list)
  
  return(metrics)
}



identify_true_asvs <- function(asv_list, seq_dict, ground_truth_taxa) {
  
  index <- match(asv_list, seq_dict$ASV)
  
  asv_genera <- seq_dict$Genus[index]
  
  is_true_match <- !is.na(asv_genera) & (asv_genera %in% ground_truth_taxa$genus)
  
  return(asv_list[is_true_match])
}



calculate_metrics_performance <- function(filtered_community, raw_community, seq_dict, ground_truth, penalty_weight = 1, method_name ="Method"){
  
  original_asvs <- colnames(raw_community)
  
  retained_asvs <- colnames(filtered_community)
  
  true_asvs_list <- identify_true_asvs(original_asvs, seq_dict, ground_truth)
  
  total_true<- length(true_asvs_list)
  total_false <- length(original_asvs) - total_true
  
  tp_asvs <- intersect(retained_asvs, true_asvs_list)
  TP <- length(tp_asvs)
  
  fp_asvs <- setdiff(retained_asvs, true_asvs_list)
  FP <- length(fp_asvs)
  
  FN <- total_true - TP
  TN <- total_false - FP
  
  # Metrics
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  recall <- ifelse(total_true > 0, TP / total_true, 0)
  f1 <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0)
  fdr <- ifelse((TP + FP) > 0, FP/ (TP+FP),0)
  fpr <- ifelse(total_false > 0, FP / total_false, 0)
  
  # Ratio TRUE/FALSE
  ratio_tf <- ifelse(FP > 0, TP / FP, Inf)
  
  #F0.5 Score
  beta_val <- 0.5
  beta_sq <- beta_val^2
  
  f0_5 <- ifelse((beta_sq * precision) + recall > 0,
                 ((1 + beta_sq) * precision * recall) / ((beta_sq * precision) + recall),
                 0)
  #F1-FPR
  f1_fpr <- f1 - (penalty_weight * fpr)
  
  #F1/FPR
  f1_over_fpr <- f1/(1 + penalty_weight*fpr)
  
  #f balanced
  f_balanced <- 3/((1/precision) + (1/recall) + (1/(1-fpr)))
  
  # Matthews correlation coefficient (MCC)
  
  mcc = ((TP*TN) - (FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  
  return(data.frame(
    Method_name = method_name,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    Total_True = total_true,
    Total_False = total_false,
    Ratio_TF =ratio_tf,
    Precision = precision,
    Recall = recall,
    F1 = f1,
    FDR = fdr,
    FPR = fpr,
    F0.5 = f0_5,
    "F1-FPR" = f1_fpr,
    "F1/FPR" = f1_over_fpr,
    F_Balanced = f_balanced,
    MCC = mcc
  ))
}

k2_validation<- function(OTU_table, alpha, beta, filter_parameters, z_threshold, k_threshold){
  
  start.time <- Sys.time() 
  
  asv <- as.data.frame(OTU_table) %>%
    mutate(across(where(is.integer), as.numeric)) %>%
    select(where(is.numeric))
  
  k_unfiltered <- calculate_k_w_single(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  filter_arguments <- c(list(OTU_table=asv), filter_parameters)
  filter_output <- do.call(conservative_filter, filter_arguments)
  
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  k_0 <- calculate_k_w_single(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  added_asv <- k_values <- BC_values <- J_values <- c()
  
  if (length(discarded_asv) > 0) {
    
    for (asv_name in discarded_asv) {
      
      # Temporaly reintroducing the deleted ASV
      asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
      
      # Calculate metrics
      k_i <- calculate_k_w_single(asv_table_i, alpha, beta)
      bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
      j_i <- mean(na.omit(evenness_k(asv_table_i)))
      
      # Save results
      added_asv <- c(added_asv, asv_name)
      k_values <- c(k_values, k_i)
      BC_values <- c(BC_values, bc_i)
      J_values <- c(J_values, j_i)
    }
  }
  
  # Construct results
  if (length(added_asv) > 0) {
    df_results <- data.frame(
      added_asv = added_asv,
      k_values = k_values,
      delta_k = k_values - k_0,
      BC_values = BC_values,
      J_values = J_values
    )
    k_max_val <- max(k_values, na.rm = TRUE)
  } else {
    
    df_results <- data.frame(
      added_asv = character(), k_values = numeric(), 
      delta_k = numeric(), BC_values = numeric(), J_values = numeric()
    )
    k_max_val <- k_0 
  }
  
  cutoff_used <- NA
  method_used <- "None"
  
  if (nrow(df_results) > 0) {
    if (!is.null(k_threshold)) {
      # Fixed relative k threshold
      cutoff_used <- k_threshold
      method_used <- "Fixed k"
    } else {
      # Z-score 
      mean_k <- mean(df_results$k_values, na.rm = TRUE)
      sd_k   <- sd(df_results$k_values, na.rm = TRUE)
      cutoff_used <- mean_k + (z_threshold * sd_k)
      method_used <- paste0("Z-score (", z_threshold, "sd)")
    }
    
    # Identify winners
    recovered_asvs <- df_results$added_asv[df_results$k_values >= cutoff_used]
  } else {
    recovered_asvs <- character(0)
  }
  
  # Construct the final table
  if (length(recovered_asvs) > 0) {
    # # recovered the original data from the winners
    data_recovered <- asv[, recovered_asvs, drop = FALSE]
    # Join the preserved ones
    final_asv_table <- cbind(asv_preserved, data_recovered)
  } else {
    # If nothing is recovered, the final is the same as the preserved one
    final_asv_table <- asv_preserved
  }
  
  end.time <- Sys.time()
  
  return(list(
    
    df_results = df_results,
    
    final_asv_table = final_asv_table,
    recovered_count = length(recovered_asvs),
    threshold_used = cutoff_used,
    
    k_inicial = as.numeric(k_0),
    stats_post_filter = list(k=k_0, BC=BC_0, J=J_0),
    stats_raw = list(k=k_unfiltered, BC=BC_unfiltered, J=J_unfiltered),
    k_max = k_max_val,
    time = round(end.time - start.time, 2)
  ))
}


k4_validation<- function(OTU_table, alpha, beta, filter_parameters, z_threshold, k_threshold){
  
  start.time <- Sys.time() 
  
  asv <- as.data.frame(OTU_table) %>%
    mutate(across(where(is.integer), as.numeric)) %>%
    select(where(is.numeric))
  
  k_unfiltered <- k4(asv, alpha, beta)
  JC_unfiltered <- mean(na.omit(jaccard(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  filter_arguments <- c(list(OTU_table=asv), filter_parameters)
  filter_output <- do.call(conservative_filter, filter_arguments)
  
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  k_0 <- k4(asv_preserved, alpha, beta)
  JC_0 <- mean(na.omit(jaccard(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  added_asv <- k_values <- JC_values <- J_values <- c()
  
  if (length(discarded_asv) > 0) {
    
    for (asv_name in discarded_asv) {
      
      # Temporaly reintroducing the deleted ASV
      asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
      
      # Calculate metrics
      k_i <- k4(asv_table_i, alpha, beta)
      jc_i <- mean(na.omit(jaccard(asv_table_i)))
      j_i <- mean(na.omit(evenness_k(asv_table_i)))
      
      # Save results
      added_asv <- c(added_asv, asv_name)
      k_values <- c(k_values, k_i)
      JC_values <- c(JC_values, jc_i)
      J_values <- c(J_values, j_i)
    }
  }
  
  # Results contruction
  if (length(added_asv) > 0) {
    df_results <- data.frame(
      added_asv = added_asv,
      k_values = k_values,
      delta_k = k_values - k_0,
      JC_values = JC_values,
      J_values = J_values
    )
    k_max_val <- max(k_values, na.rm = TRUE)
  } else {
    
    df_results <- data.frame(
      added_asv = character(), k_values = numeric(), 
      delta_k = numeric(), JC_values = numeric(), J_values = numeric()
    )
    k_max_val <- k_0 
  }
  
  cutoff_used <- NA
  method_used <- "None"
  
  if (nrow(df_results) > 0) {
    if (!is.null(k_threshold)) {
      # Fixed relative k threshold
      cutoff_used <- k_threshold
      method_used <- "Fixed k"
    } else {
      # Z-score 
      mean_k <- mean(df_results$k_values, na.rm = TRUE)
      sd_k   <- sd(df_results$k_values, na.rm = TRUE)
      cutoff_used <- mean_k + (z_threshold * sd_k)
      method_used <- paste0("Z-score (", z_threshold, "sd)")
    }
    
    # Identify winners
    recovered_asvs <- df_results$added_asv[df_results$k_values >= cutoff_used]
  } else {
    recovered_asvs <- character(0)
  }
  
  #  Contruct the final table
  if (length(recovered_asvs) > 0) {
    # recovered the original data from the winners
    data_recovered <- asv[, recovered_asvs, drop = FALSE]
    # Join the preserved ones
    final_asv_table <- cbind(asv_preserved, data_recovered)
  } else {
    # If nothing is recovered, the final is the same as the preserved one
    final_asv_table <- asv_preserved
  }
  
  end.time <- Sys.time()
  
  return(list(
    
    df_results = df_results,
    final_asv_table = final_asv_table,
    recovered_count = length(recovered_asvs),
    threshold_used = cutoff_used,
    k_inicial = as.numeric(k_0),
    stats_post_filter = list(k=k_0, JC=JC_0, J=J_0),
    stats_raw = list(k=k_unfiltered, JC=JC_unfiltered, J=J_unfiltered),
    k_max = k_max_val,
    time = round(end.time - start.time, 2)
  ))
}

# PLot for look the recovery for a z-score distance of initial k
plot_z_threshold <- function(results, z_threshold, top_n = 100) {
  
  df_results <- results$df_results
  k_inicial  <- results$k_inicial
  
  if(nrow(df_results) == 0) {
    warning("No data to graph.")
    return(NULL)
  }
  
  mean_k <- mean(df_results$k_values, na.rm = TRUE)
  sd_k   <- sd(df_results$k_values, na.rm = TRUE)
  
  
  k_cutoff <- mean_k + (z_threshold * sd_k)
  
  
  df_full <- df_results %>%
    mutate(Status = ifelse(k_values >= k_cutoff, "Recovered", "Noise"))
  
  
  df_plot_k <- df_full %>%
    arrange(desc(k_values)) %>%    
    slice_head(n = top_n) %>%      
    arrange(k_values)              
  
  df_plot_k$added_asv <- factor(df_plot_k$added_asv, levels = df_plot_k$added_asv)
  
  
  plot <- ggplot(df_plot_k, aes(x = added_asv, y = k_values)) +
    
   
    geom_line(aes(group = 1), color = "gray80", linetype = "solid", alpha = 0.5) +
    
    # Points
    geom_point(aes(colour = Status), size = 3, alpha = 0.9) + 
    
    #  K inicial
    geom_hline(yintercept = k_inicial, linetype = "dashed", color = "red", linewidth = 0.7) +
    
    # Z-score
    geom_hline(yintercept = k_cutoff, linetype = "dotted", color = "blue", linewidth = 0.8) +
    
    
    annotate("text", x = 1, y = k_inicial, label = paste("k initial:", round(k_inicial, 4)), 
             vjust = -1, hjust = 0, color = "red", size = 3, fontface = "bold") +
    
    annotate("text", x = 1, y = k_cutoff, label = paste0("Z-score threshold (", z_threshold, "σ)"), 
             vjust = 1.5, hjust = 0, color = "blue", size = 3, fontface = "italic") +
    
    
    labs(title = paste("Top", nrow(df_plot_k), "ASV Candidates Evaluation"),
         subtitle = paste("Showing top performing ASVs. Thresholds calculated on full dataset (n=", nrow(df_full), ")"),
         x = "ASV Candidate (Ordered by k)",
         y = "New k value",
         color = "Classification") +
    
    scale_color_manual(
      values = c("Recovered" = "forestgreen", 
                 "Noise"     = "gray50"),       
      labels = c("Recovered" = "Candidate Selected",
                 "Noise"     = "Background Noise")
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "bottom"
    )
  
  return(plot)
}

# PLot for look the recovery for a fixed k thershold
plot_k_threshold <- function(results, k_threshold, top_n = 100) {
  
  df_results <- results$df_results
  k_inicial  <- results$k_inicial
  
  # Direct clasification based on fixed threshold
  df_full <- df_results %>%
    mutate(Status = ifelse(k_values >= k_threshold, "Recovered", "Noise"))
  
  # PData preparaiton
  df_plot_k <- df_full %>%
    arrange(desc(k_values)) %>%      
    slice_head(n = top_n) %>%        
    arrange(k_values)               
  
  df_plot_k$added_asv <- factor(df_plot_k$added_asv, levels = df_plot_k$added_asv)
  
  
  plot <- ggplot(df_plot_k, aes(x = added_asv, y = k_values)) +
    
    # Connection line
    geom_line(aes(group = 1), color = "gray80", linetype = "solid", alpha = 0.5) +
    
    # Points coloured by Status
    geom_point(aes(colour = Status), size = 2, alpha = 0.7) + 
    
    # Line for initial k
    geom_hline(yintercept = k_inicial, linetype = "dashed", color = "red", linewidth = 0.7) +
    
    # Line for fixed threshold 
    geom_hline(yintercept = k_threshold, linetype = "dotted", color = "blue", linewidth = 0.8) +
    
    # Tag for initial k
    annotate("text", x = 1, y = k_inicial, label = paste("k initial:", round(k_inicial, 6)), 
             vjust = -1, hjust = 0, color = "red", size = 3, fontface = "bold") +
    
    # Tag fixed threshold
    annotate("text", x = 1, y = k_threshold, label = paste("k Threshold:", round(k_threshold, 6)), 
             vjust = 1.5, hjust = 0, color = "blue", size = 3, fontface = "italic") +
    
    
    labs(title = paste("Top", nrow(df_plot_k), "ASV Candidates Evaluation"),
         subtitle = paste("Classification based on fixed k threshold of", round(k_threshold, 6)),
         x = "ASV Candidate (Ordered by k)",
         y = "New k value",
         color = "Classification") +
    
    scale_color_manual(
      values = c("Recovered" = "forestgreen", 
                 "Noise"     = "gray50"),        
      labels = c("Recovered" = "Candidate Selected",
                 "Noise"     = "Background Noise")
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "bottom"
    )
  
  return(plot)
}

# To see the recovered ASVs
get_candidates <- function(results, z_threshold, k_threshold) {
  
  df_results <- results$df_results
  
  if(nrow(df_results) == 0) return(NULL)
  
  if (!is.null(k_threshold)) {
    # With k value proporcionated
    k_cutoff <- k_threshold
    txt_method <- "Relative threshold of k"
    
  } else {
    # With Z-score
    mean_k <- mean(df_results$k_values, na.rm = TRUE)
    sd_k   <- sd(df_results$k_values, na.rm = TRUE)
    k_cutoff <- mean_k + (z_threshold * sd_k)
    txt_method <- paste("Z-score (", z_threshold, "sigma )")
  }
  
  # Filtering(some for both methods)
  candidates <- df_results %>%
    filter(k_values >= k_cutoff) %>%
    pull(added_asv)
  
  return(as.character(candidates))
}

# Function to calculate metrics of alpha diversity
calculate_alpha <- function(otu_table, method_name) {
  
  richness <- specnumber(otu_table)
  shannon  <- diversity(otu_table, index = "shannon")
  simpson  <- diversity(otu_table, index = "simpson")
  
  # Chao1 (If Chao1 cannot be calculated, observed richness is used)
  chao1 <- tryCatch({
    t(estimateR(otu_table))[, "S.chao1"]
  }, error = function(e) richness)
  
  # Pielou (avoid division by zero if there is only 1 ASV)
  pielou <- ifelse(richness > 1, shannon / log(richness), 0)
  
  data.frame(
    SampleID = rownames(otu_table),
    Filter_Method = method_name,
    Chao1 = chao1,
    Shannon = shannon,
    Simpson = simpson,
    Pielou = pielou
  )
}

# Function to calculate coords for PCoA
calculate_pcoa_coords <- function(otu_table, method_name, distance) {
  
  # Previous Normalization
  if(distance == "bray") {
    # Relative abundance for bray
    mat <- decostand(otu_table, method = "total")
  } else {
    # Presence/Absence for Jaccard
    mat <- decostand(otu_table, method  = "pa")
  }
  
  # Distance and PCoA
  dist_mat <- vegdist(mat, method = distance)
  pcoa <- cmdscale(dist_mat, k = 2, eig = TRUE)
  
  # Explained variance
  var_expl <- round(pcoa$eig / sum(pcoa$eig) * 100, 1)
  
  # Dataframe of coordinates
  data.frame(
    PC1 = pcoa$points[,1],
    PC2 = pcoa$points[,2],
    Filter_Method = method_name,
    Var1 = var_expl[1], 
    Var2 = var_expl[2]
  )
}

