library(tidyverse)
library(progress)
library(stats)
library(kSamples)

# --------------------------------- HKS/AD TEST ---------------------------------

hierarchical_permutation_test <- function(df,
                                          value_col,
                                          subject_col,
                                          cell_col,
                                          group_col,
                                          metric = "ks",
                                          n_cells_per_subject = NULL,
                                          n_obs_per_cell = NULL,
                                          n_resamples = 1000,
                                          n_permutations = 1000,
                                          replace = FALSE,
                                          seed = NULL) {
  
  
  if (!is.null(seed)) set.seed(seed)
  df <- as_tibble(df)
  
  # --------------------------------- 1. Group Validation ---------------------------------
  
  groups <- unique(df[[group_col]])
  if (length(groups) != 2) {
    stop(sprintf("The column '%s' must contain exactly 2 groups.", group_col))
  }
  
  # --------------------------------- 2. Setting Number of Observations ---------------------------------
  
  actual_min_obs <- df %>%
    group_by(!!sym(subject_col), !!sym(cell_col)) %>%
    count() %>%
    pull(n) %>%
    min()
  
  if (is.null(n_obs_per_cell)) {
    n_obs <- actual_min_obs
    message(sprintf("INFO: Auto-detected n_obs_per_cell = %d", n_obs))
  } else {
    if (!replace && n_obs_per_cell > actual_min_obs) {
      stop(sprintf("ERROR: Requested %d observations, but min available is %d. Set 'replace = TRUE'.", 
                   n_obs_per_cell, actual_min_obs))
    } else if (replace && n_obs_per_cell > actual_min_obs) {
      warning(sprintf("WARNING: Requested %d observations (available: %d). Oversampling used.", 
                      n_obs_per_cell, actual_min_obs))
    }
    
    n_obs <- n_obs_per_cell
  }
  
  # --------------------------------- 3. Setting Number of Cells ---------------------------------
  
  actual_min_cells <- df %>%
    group_by(!!sym(subject_col)) %>%
    summarise(n_unique = n_distinct(!!sym(cell_col)), .groups = "drop") %>%
    pull(n_unique) %>%
    min()
  
  if (is.null(n_cells_per_subject)) {
    n_cells <- actual_min_cells
    message(sprintf("INFO: Auto-detected n_cells_per_subject = %d", n_cells))
  } else {
    if (!replace && n_cells_per_subject > actual_min_cells) {
      stop(sprintf("ERROR: Requested %d cells, but min available is %d. Set 'replace = TRUE'.", 
                   n_cells_per_subject, actual_min_cells))
    } else if (replace && n_cells_per_subject > actual_min_cells) {
      warning(sprintf("WARNING: Requested %d cells (available: %d). Oversampling used.", 
                      n_cells_per_subject, actual_min_cells))
    }
    
    n_cells <- n_cells_per_subject
  }
  
  # --------------------------------- 4. Resampling and Test ---------------------------------
  
  hierarchical_resample <- function(data) {
    
    # Observations
    balanced_obs <- data %>%
      group_by(!!sym(subject_col), !!sym(cell_col)) %>%
      slice_sample(n = n_obs, replace = replace) %>%
      ungroup()
    
    # Cells
    selected_cells <- balanced_obs %>%
      select(!!sym(subject_col), !!sym(cell_col)) %>%
      distinct() %>%
      group_by(!!sym(subject_col)) %>%
      slice_sample(n = n_cells, replace = replace) %>%
      ungroup
  
    # Merge
    inner_join(balanced_obs, 
               selected_cells, 
               by = c(subject_col, cell_col),
               relationship = "many-to-many")
  }
  
  compute_stat <- function(data) {
    g1 <- data[[value_col]][data[[group_col]] == groups[1]]
    g2 <- data[[value_col]][data[[group_col]] == groups[2]]
    if (metric == "ks") {
      return(suppressWarnings(ks.test(g1, g2)$statistic))
    } else if (metric == "ad") {
      return(suppressWarnings(kSamples::ad.test(g1, g2)$ad[1, 1]))
    }
  }
  
  # --------------------------------- 5. Bootstrap --------------------------------- 
  
  message("Calculating observed distribution...")
  pb1 <- progress_bar$new(total = n_resamples)
  obs_stats <- replicate(n_resamples, {
    pb1$tick()
    resampled <- hierarchical_resample(df)
    compute_stat(resampled)
  })
  
  obs_median <- median(obs_stats)
  
  message("Calculating null distribution...")
  subjects <- unique(df[[subject_col]])
  subject_labels <- df %>%
    select(all_of(c(subject_col, group_col))) %>%
    distinct()
  
  pb2 <- progress_bar$new(total = n_permutations)
  null_stats <- replicate(n_permutations, {
    pb2$tick()
    
    shuffled_labels <- sample(subject_labels[[group_col]])
    mapping <- setNames(shuffled_labels, subject_labels[[subject_col]])
    
    permuted_df <- df 
    permuted_df[[group_col]] <- mapping[as.character(permuted_df[[subject_col]])]
    
    resampled_perm <- hierarchical_resample(permuted_df)
    compute_stat(resampled_perm)
  })
  
  p_val <- mean(null_stats >= obs_median)
  
  return(list(
    p_value = p_val,
    obs_median = obs_median,
    obs_dist = obs_stats,
    null_dist = null_stats,
    metric = metric
  ))
}

# --------------------------------- PLOTTING ---------------------------------

plot_hks_results <- function(result_list, title = "Hierarchical Permutation Test") {
  
  df_plot <- data.frame(
    stat = c(result_list$null_dist, result_list$obs_dist),
    type = factor(c(rep("Null Distribution (Permuted)", length(result_list$null_dist)), 
                    rep("Observed Distribution (Resampled)", length(result_list$obs_dist))),
                  levels = c("Null Distribution (Permuted)", "Observed Distribution (Resampled)"))
  )
  
  metric_label <- if(is.null(result_list$metric)) "KS Statistic" else paste(toupper(result_list$metric), "Statistic")
  
  ggplot(df_plot, aes(x = stat, fill = type)) +

    geom_histogram(aes(y = after_stat(density)), 
                   alpha = 0.3, bins = 60, position = "identity", color = "black") +

    geom_vline(aes(xintercept = result_list$obs_median, color = "Observed Median"), 
               linetype = "dashed", linewidth = 1, show.legend = TRUE) +
    
    scale_fill_manual(values = c("Null Distribution (Permuted)" = "black", 
                                 "Observed Distribution (Resampled)" = "royalblue")) +
    scale_color_manual(name = "", values = c("Observed Median" = "red")) +
    
    labs(title = title, 
         subtitle = sprintf("p-value: %.5f", result_list$p_value),
         x = metric_label,
         y = "Density",
         fill = "Distributions") +
    
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),           
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
      axis.line = element_line(colour = "black"), 
      legend.position = "bottom"
    )
}

# --------------------------------- SIMULATED DATA ---------------------------------

generate_hierarchical_data <- function(n_basal,
                                       n_treatment,
                                       params_basal,
                                       params_treatment) {
  
  
  conditions <- list(
    list(name = "B", n = n_basal, params = params_basal),
    list(name = "T", n = n_treatment, params = params_treatment)
  )
  
  simulated_df <- map_df(conditions, function(cond) {
    map_df(seq_len(cond$n), function(s_idx) {
      subject_id <- paste0("S", cond$name, s_idx - 1)
      n_cells <- sample(10:60, 1)
      
      map_df(seq_len(n_cells), function(c_idx) {
        cell_id <- paste0("CELL_", subject_id, "_", c_idx - 1)
        n_samples <- sample(25:100, 1)
        samples <- rlnorm(n = n_samples,
                          meanlog = cond$params$mean,
                          sdlog = cond$params$sigma)
        
        tibble(
          subject = subject_id,
          subunit = cell_id,
          group = cond$name,
          variable = samples
        )
      })
    }) 
  })
  
  return(simulated_df)
}

# --------------------------------- TESTING ---------------------------------

# Generating the simulated dataframe
p_basal <- list(mean = 0.5, sigma = 0.3)
p_treatment <- list(mean = 0.5, sigma = 0.35)

df1 <- generate_hierarchical_data(
  n_basal = 12,
  n_treatment = 8,
  params_basal = p_basal,
  params_treatment = p_treatment
)

# Performing the HKS test
hks_result <- hierarchical_permutation_test(df = df1,
                                            value_col = 'variable',
                                            subject_col = 'subject',
                                            cell_col = 'subunit',
                                            group_col = 'group',
                                            metric = 'ad',
                                            n_resamples = 2000,
                                            n_permutations = 2000,
                                            replace = TRUE,
                                            seed = 123)

plot_hks_results(result_list = hks_result)

# Same test but with the AD metric
had_result <- hierarchical_permutation_test(df = df1,
                                            value_col = 'variable',
                                            subject_col = 'subject',
                                            cell_col = 'subunit',
                                            group_col = 'group',
                                            metric = 'ad',
                                            n_resamples = 2000,
                                            n_permutations = 2000,
                                            replace = TRUE,
                                            seed = 123)

plot_hks_results(result_list = had_result)
