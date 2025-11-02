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


# Conservative filter
conservative_filter <- function(OTU_table) {
  
  num_samples <- nrow(OTU_table)
  total_reads_dataset <- sum(OTU_table)
  
  if (num_samples == 0 || total_reads_dataset == 0) {
   
    # Handling empty tables
    return(list(preserved = OTU_table[0, 0, drop = FALSE], 
                deleted = OTU_table[0, 0, drop = FALSE]))
  }
  
  col_totals <- colSums(OTU_table)
  col_maxs <- apply(OTU_table, 2, max) # The maximum per ASV
  samples_present <- apply(OTU_table, 2, function(c) sum(c > 0)) # Prevalence
  
  # Rule 1: (Max <= 10 reads) and (present in <= 5% of the samples)
  pct_samples_present <- samples_present / num_samples
  C1_delete <- (col_maxs <= 10) & (pct_samples_present <= 0.05)
  
  # Rule 2: < 0.005% total relative abundance (0.00005)
  relative_abundance <- col_totals / total_reads_dataset
  C2_delete <- relative_abundance < 0.00001
  
  # Rule 3: Appears in < 2 samples (only 1)
  C3_delete <- samples_present < 2
  
  # Combine
  asv_to_delete_bool <- C1_delete | C2_delete | C3_delete
  
  asv_ids_to_preserve <- colnames(OTU_table)[!asv_to_delete_bool]
  asv_ids_to_delete <- colnames(OTU_table)[asv_to_delete_bool]
  
  # Two lists. One with the preserved ASVs and another with the filtered ones
  return(list(
    preserved = OTU_table[, asv_ids_to_preserve, drop = FALSE],
    deleted = OTU_table[, asv_ids_to_delete, drop = FALSE]
  ))
}



# Add-one procedure
k_validation <- function(OTU_table, ground_truth, alpha, beta, mode){
  
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
  true_positive_flags <- c() 
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- calculate_k_w_single(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    j_i <- mean(na.omit(evenness_k(asv_table_i)))
    
    # Evaluate with ground_truth depending on the chosen mode
    
    if (mode == "taxonomy") {
      asv_identity <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
      match_truth <- asv_identity %in% ground_truth_species
    } else if (mode == "sequence") {
      asv_identity <- seq_dict$Sequence[match(asv_name, seq_dict$ASV)]
      match_truth <- any(sapply(ground_truth_sequences$Sequence, function(gt) {
        grepl(asv_identity, gt, fixed = TRUE) || grepl(gt, asv_identity, fixed = TRUE)
      }))
    }
    
    added_asv <- c(added_asv, asv_name)
    k_values <- c(k_values, k_i)
    BC_values <- c(BC_values, bc_i)
    J_values <- c(J_values, j_i)
    true_positive_flags <- c(true_positive_flags, match_truth)
    
  }
  
  # Select ASVs that actually went up k
  reintroduced_asvs <- added_asv[k_values > k_0]
  reintroduced_flags <- true_positive_flags[k_values > k_0]
  
  TP <- sum(reintroduced_flags)
  FP <- sum(!reintroduced_flags)
  
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
  recall <- ifelse(length(reintroduced_asvs) > 0, TP / length(reintroduced_asvs), NA)
  F1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall == 0),
               NA,
               2 * (precision * recall) / (precision + recall))
  FDR <- ifelse((TP + FP) > 0, FP / (TP + FP), NA)
  
  k_max <- max(k_values, na.rm = TRUE)
  
  delta_k <- k_max - k_0
  delta_k_rel <- ifelse(k_0 != 0, delta_k / k_0, NA)
  
  # Dataframe por ASV
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    BC_values = BC_values,
    J_values = J_values,
    TP_flag = true_positive_flags
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = k_max,
    delta_k = delta_k,
    delta_k_rel = delta_k_rel,
    asvs_recuperadas = reintroduced_asvs,
    filtered_values = data.frame(k = k_0, BC = BC_0, J = J_0),
    unfiltered_values = data.frame(k = k_unfiltered, BC = BC_unfiltered, J = J_unfiltered),
    global_metrics = data.frame(TP = TP, FP = FP, precision = precision,
                                recall = recall, FDR = FDR, F1 = F1),
    time = time.take
  ))
}

# Add-one procedure
k_validation2 <- function(OTU_table, ground_truth, alpha, beta, mode){
  
  start.time <- Sys.time() 
  mode <- match.arg(mode, choices = c("taxonomy", "sequence"))
  
  # Preprocessed to ensure all columns are numeric
  asv <- OTU_table %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered OTU Table Values
  k_unfiltered <- calculate_k_w_single2(asv, alpha, beta)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_k(asv)))
  
  # apply conservative filter
  filter_output <- conservative_filter(asv)
  asv_preserved <- filter_output$preserved
  asv_deleted <- filter_output$deleted
  
  discarded_asv <- colnames(asv_deleted)
  keep_asv <- colnames(asv_preserved)
  
  # Initial values after conservative filtering
  k_0 <- calculate_k_w_single2(asv_preserved, alpha, beta)
  BC_0 <- mean(na.omit(bray_curtis(asv_preserved)))
  J_0 <- mean(na.omit(evenness_k(asv_preserved)))
  
  # Create vectors to record the results
  added_asv <- k_values <- BC_values <- J_values <- c()
  true_positive_flags <- c() 
  
  # Pipeline add-one
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_preserved, asv[, asv_name, drop = FALSE])
    
    k_i <- calculate_k_w_single2(asv_table_i, alpha, beta)
    bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
    j_i <- mean(na.omit(evenness_k(asv_table_i)))
    
    # Evaluate with ground_truth depending on the chosen mode
    
    if (mode == "taxonomy") {
      asv_identity <- seq_dict$taxa[match(asv_name, seq_dict$ASV)]
      match_truth <- asv_identity %in% ground_truth_species
    } else if (mode == "sequence") {
      asv_identity <- seq_dict$Sequence[match(asv_name, seq_dict$ASV)]
      match_truth <- any(sapply(ground_truth_sequences$Sequence, function(gt) {
        grepl(asv_identity, gt, fixed = TRUE) || grepl(gt, asv_identity, fixed = TRUE)
      }))
    }
    
    added_asv <- c(added_asv, asv_name)
    k_values <- c(k_values, k_i)
    BC_values <- c(BC_values, bc_i)
    J_values <- c(J_values, j_i)
    true_positive_flags <- c(true_positive_flags, match_truth)
    
  }
  
  # Select ASVs that actually went up k
  reintroduced_asvs <- added_asv[k_values > k_0]
  reintroduced_flags <- true_positive_flags[k_values > k_0]
  
  TP <- sum(reintroduced_flags)
  FP <- sum(!reintroduced_flags)
  
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
  recall <- ifelse(length(reintroduced_asvs) > 0, TP / length(reintroduced_asvs), NA)
  F1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall == 0),
               NA,
               2 * (precision * recall) / (precision + recall))
  FDR <- ifelse((TP + FP) > 0, FP / (TP + FP), NA)
  
  k_max <- max(k_values, na.rm = TRUE)
  
  delta_k <- k_max - k_0
  delta_k_rel <- ifelse(k_0 != 0, delta_k / k_0, NA)
  
  # Dataframe por ASV
  df_results <- data.frame(
    added_asv = added_asv,
    k_values = k_values,
    BC_values = BC_values,
    J_values = J_values,
    TP_flag = true_positive_flags
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    df_results = df_results,
    k_inicial = k_0,
    k_max = k_max,
    delta_k = delta_k,
    delta_k_rel = delta_k_rel,
    asvs_recuperadas = reintroduced_asvs,
    filtered_values = data.frame(k = k_0, BC = BC_0, J = J_0),
    unfiltered_values = data.frame(k = k_unfiltered, BC = BC_unfiltered, J = J_unfiltered),
    global_metrics = data.frame(TP = TP, FP = FP, precision = precision,
                                recall = recall, FDR = FDR, F1 = F1),
    time = time.take
  ))
}



bootstrap_k <- function(OTU_table, asv_base, asv_candidate, 
                        approach = 1, alpha_weight = 0.5, beta_weight = 0.5, 
                        n_boot = 1500, alpha = 0.05) {
  
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
  if (approach == 1) {
    p_value <- (sum(delta_boot >= delta_k_obs) + 1) / (length(delta_boot) + 1)
  } else if (enfoque == 2) {
    p_value <- (sum(delta_boot >= delta_k_obs) + 1) / (length(delta_boot) + 1)
  }
  
  signif <- if (enfoque == 1) p_value < alpha else p_value > (1 - alpha)
  
  list(delta_k_obs = delta_k_obs, p_value = p_value, signif = signif)
}







                          |------------------------------------|
# Plots
plot_metrics <- function(df_results) {
  
  # Asegurarse de que las ASVs aparezcan en orden de reintroducción
  df_long <- df_results %>%
    mutate(order = seq_along(added_asv)) %>%
    pivot_longer(
      cols = c(k_values, precision, recall, F1),
      names_to = "metric",
      values_to = "value"
    )
  
  # Etiquetas más limpias
  metric_labels <- c(
    k_values = "Índice k",
    precision = "Precisión",
    recall = "Recall (Sensibilidad)",
    F1 = "F1 Score"
  )
  
  # Plot multipanel con facetado
  ggplot(df_long, aes(x = order, y = value, color = metric, group = metric)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    facet_wrap(~metric, scales = "free_y", ncol = 1,
               labeller = as_labeller(metric_labels)) +
    scale_x_continuous(
      breaks = seq_along(df_results$added_asv),
      labels = df_results$added_asv
    ) +
    labs(
      title = "Evolución del índice k y métricas de clasificación durante la reintroducción secuencial",
      x = "ASV reintroducida",
      y = "Valor de la métrica"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, size = 8),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )
}

