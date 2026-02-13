library(R6)
library(data.table)
library(parallel)
library(progress)
library(ggplot2) 
library(tidyverse)
library(pbapply)

HierarchicalPermutationTest <- R6Class("HierarchicalPermutationTest",
                                       public = list(
                                         
                                         # Class properties
                                         data_structure = NULL, 
                                         original_labels = NULL,
                                         n_subjects = NULL,
                                         n_cells = NULL,
                                         n_obs = NULL,
                                         metric = NULL,
                                         replace = FALSE,
                                         observed_distribution = NULL,
                                         null_distribution = NULL,
                                         observed_stat_median = NULL,
                                         p_value = NULL,
                                         
                                         # ---------------------- 1. Initialization ----------------------
                                         
                                         initialize = function(df, 
                                                               value_col,
                                                               subject_col,
                                                               cell_col,
                                                               group_col,
                                                               metric = "ks",
                                                               n_cells_per_subject = NULL,
                                                               n_obs_per_cell = NULL,
                                                               replace = FALSE) {
                                           
                                           self$metric <- metric
                                           self$replace <- replace 
                                           
                                           # Convert to data.table for fast preprocessing
                                           dt <- as.data.table(df)
                                           
                                           groups <- unique(dt[[group_col]])
                                           if (length(groups) != 2) stop("The group column must contain exactly 2 groups.")
                                           
                                           # ---------------------- 2. Setting number of observations ---------------------
                                           
                                           actual_min_obs <- min(dt[, .N, by = c(subject_col, cell_col)]$N)
                                           
                                           if (is.null(n_obs_per_cell)) {
                                             self$n_obs <- actual_min_obs
                                             message(sprintf("INFO: Auto-detected n_obs_per_cell = %d", self$n_obs))
                                           } else {
                                             if (!replace && n_obs_per_cell > actual_min_obs) {
                                               stop(sprintf("ERROR: Requested %d obs, but min is %d. Set replace = TRUE.", n_obs_per_cell, actual_min_obs))
                                             }
                                             self$n_obs <- n_obs_per_cell
                                           }
                                           
                                           # ---------------------- 3. Setting number of cells ---------------------
                                           
                                           # Robustly get the count column (handles whether it defaults to V1 or N)
                                           cell_counts <- dt[, uniqueN(get(cell_col)), by = subject_col]
                                           actual_min_cells <- min(cell_counts[[2]]) # Grab the 2nd column (the counts)
                                           
                                           if (is.null(n_cells_per_subject)) {
                                             self$n_cells <- actual_min_cells
                                             message(sprintf("INFO: Auto-detected n_cells_per_subject = %d", self$n_cells))
                                           } else {
                                             if (!replace && n_cells_per_subject > actual_min_cells) {
                                               stop(sprintf("ERROR: Requested %d cells, but min is %d. Set replace = TRUE.", n_cells_per_subject, actual_min_cells))
                                             }
                                             self$n_cells <- n_cells_per_subject
                                           }
                                           
                                           # ---------------------- 4. Data optimization ---------------------
                                           
                                           message("INFO: Optimizing data structure for fast resampling...")
                                           
                                           subjects <- unique(dt[[subject_col]])
                                           group_map <- setNames(c(0, 1), groups)
                                           
                                           self$data_structure <- lapply(subjects, function(sub) {
                                             sub_dt <- dt[get(subject_col) == sub]
                                             grp_lbl <- group_map[[ sub_dt[[group_col]][1] ]]
                                             cells_data <- split(sub_dt[[value_col]], sub_dt[[cell_col]])
                                             
                                             list(grp_lbl, cells_data)
                                           })
                                           
                                           self$original_labels <- sapply(self$data_structure, `[[`, 1)
                                           self$n_subjects <- length(self$data_structure)
                                         },
                                         
                                         # ---------------------- 5. Fast KS computation ----------------------
                                         
                                         compute_ks_statistic = function(x, y) {
                                           z <- sort(c(x, y))
                                           cdf_x <- ecdf(x)(z)
                                           cdf_y <- ecdf(y)(z)
                                           
                                           max(abs(cdf_x - cdf_y))
                                         },
                                         
                                         # ---------------------- 6. Resampling procedure ----------------------
                                         
                                         fast_resample_and_compute = function(current_labels, seed) {
                                           set.seed(seed)
                                           
                                           # Buffer for pre-allocation
                                           est_size <- (self$n_subjects / 2) * self$n_cells * self$n_obs * 1.5
                                           g0_samples <- vector("numeric", est_size)
                                           g1_samples <- vector("numeric", est_size)
                                           g0_idx = 0
                                           g1_idx = 0
                                           
                                           for (i in seq_len(self$n_subjects)) {
                                             subj_data <- self$data_structure[[i]]
                                             cells_list <- subj_data[[2]]
                                             assigned_group <- current_labels[i]
                                             
                                             n_available_cells <- length(cells_list)
                                             
                                             # Select cells
                                             chosen_cell_indices <- sample.int(n_available_cells, self$n_cells, replace = self$replace)
                                             chosen_cells <- cells_list[chosen_cell_indices]
                                             
                                             subj_samples <- unlist(lapply(chosen_cells, function(x){
                                               if (length(x) == 1) return(rep(x, self$n_obs))
                                               sample(x, self$n_obs, replace = self$replace)
                                             }), use.names = FALSE)
                                             
                                             # Append to correct group
                                             len_s <- length(subj_samples)
                                             if (assigned_group == 0) {
                                               if ((g0_idx + len_s) > length(g0_samples)) length(g0_samples) <- length(g0_samples) * 2 
                                               g0_samples[(g0_idx + 1):(g0_idx + len_s)] <- subj_samples
                                               g0_idx <- g0_idx + len_s
                                             } else {
                                               if ((g1_idx + len_s) > length(g1_samples)) length(g1_samples) <- length(g1_samples) * 2
                                               g1_samples[(g1_idx + 1):(g1_idx + len_s)] <- subj_samples
                                               g1_idx <- g1_idx + len_s
                                             }
                                           } 
                                           
                                           # Trim results
                                           g0_final <- g0_samples[1:g0_idx]
                                           g1_final <- g1_samples[1:g1_idx]
                                           
                                           if (length(g0_final) == 0 || length(g1_final) == 0) return(0)
                                           
                                           if (self$metric == "ks") {
                                             return(self$compute_ks_statistic(g0_final, g1_final))
                                           } else {
                                             return(suppressWarnings(kSamples::ad.test(g0_final, g1_final)$ad[1, 1]))
                                           }
                                         },
                                         
                                         # ------------------ 7. Run (Parallelized) ------------------
                                         
                                         run = function(n_resamples = 1000, n_permutations = 1000, n_cores = NULL, seed = NULL) {
                                           
                                           if (!is.null(seed)) set.seed(seed)
                                           
                                           # Detect CPU cores
                                           if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
                                           if (n_cores < 1) n_cores <- 1
                                           
                                           # Export data for workers
                                           worker_data <- self$data_structure
                                           worker_n_cells <- self$n_cells
                                           worker_n_obs <- self$n_obs
                                           worker_replace <- self$replace
                                           worker_metric <- self$metric
                                           worker_n_sub <- self$n_subjects
                                           
                                           # Define worker function
                                           worker_fun <- function(arg_list) {
                                             current_labels <- arg_list$labels
                                             seed_val <- arg_list$seed
                                             set.seed(seed_val)
                                             
                                             g0_vals <- numeric(0)
                                             g1_vals <- numeric(0)
                                             total_obs_est <- (worker_n_sub / 2) * worker_n_cells * worker_n_obs * 1.5
                                             g0_vals <- vector("numeric", total_obs_est)
                                             g1_vals <- vector("numeric", total_obs_est)
                                             c0 <- 0
                                             c1 <- 0
                                             
                                             for (i in seq_len(worker_n_sub)) {
                                               subj_obj <- worker_data[[i]]
                                               cells_list <- subj_obj[[2]]
                                               grp <- current_labels[i] 
                                               
                                               n_avail <- length(cells_list)
                                               chosen_idx <- sample.int(n_avail, worker_n_cells, replace = worker_replace)
                                               chosen_cells <- cells_list[chosen_idx]
                                               
                                               obs_vec <- unlist(lapply(chosen_cells, function(x) {
                                                 if (length(x) == 1) return(rep(x, worker_n_obs))
                                                 sample(x, worker_n_obs, replace = worker_replace)
                                               }), use.names = FALSE)
                                               
                                               len <- length(obs_vec)
                                               if (grp == 0) {
                                                 if (c0 + len > length(g0_vals)) length(g0_vals) <- length(g0_vals) * 2
                                                 g0_vals[(c0 + 1):(c0 + len)] <- obs_vec
                                                 c0 <- c0 + len
                                               } else {
                                                 if (c1 + len > length(g1_vals)) length(g1_vals) <- length(g1_vals) * 2
                                                 g1_vals[(c1 + 1):(c1 + len)] <- obs_vec
                                                 c1 <- c1 + len
                                               }
                                             }
                                             
                                             g0_fin <- g0_vals[1:c0]
                                             g1_fin <- g1_vals[1:c1]
                                             
                                             if (worker_metric == "ks") {
                                               w <- c(g0_fin, g1_fin)
                                               z <- sort(w)
                                               max(abs(ecdf(g0_fin)(z) - ecdf(g1_fin)(z)))
                                             } else {
                                               suppressWarnings(kSamples::ad.test(g0_fin, g1_fin)$ad[1, 1])
                                             }
                                           }
                                           
                                           # ----------------------- 7.1 Observed distribution -----------------------
                                           
                                           cl_backend <- NULL
                                           
                                           if (.Platform$OS.type == "unix") {
                                             cl_backend <- n_cores
                                           } else {
                                             cl_backend <- makeCluster(n_cores)
                                             clusterExport(cl_backend, c("worker_data", "worker_n_cells", "worker_n_obs", 
                                                                         "worker_replace", "worker_metric", "worker_n_sub"), 
                                                           envir = environment())
                                             clusterEvalQ(cl_backend, library(kSamples)) 
                                           }
                                          
                                           message(sprintf("Calculating Observed Distribution (%d resamples)...", n_resamples))
                                           obs_seeds <- sample.int(1e9, n_resamples)
                                           obs_args <- lapply(obs_seeds, function(s) list(labels = self$original_labels, seed = s))
                                           
                                           self$observed_distribution <- unlist(pblapply(obs_args, worker_fun, cl = cl_backend))
                                           self$observed_stat_median <- median(self$observed_distribution)
                                           
                                           # ----------------------- 7.2 Null distribution -----------------------
                                           
                                           message(sprintf("Calculating Null Distribution (%d permutations)...", n_permutations))
                                           null_seeds <- sample.int(1e9, n_permutations)
                                           null_args <- lapply(null_seeds, function(s) {
                                             set.seed(s)
                                             list(labels = sample(self$original_labels), seed = s)
                                           })
                                           
                                           self$null_distribution <- unlist(pblapply(null_args, worker_fun, cl = cl_backend))
                                           
                                           # Cleanup Windows cluster
                                           if (.Platform$OS.type != "unix") {
                                             stopCluster(cl_backend)
                                           }
                                           
                                           # ----------------------- 7.3 P-value -----------------------
                                           
                                           self$p_value <- mean(self$null_distribution >= self$observed_stat_median)
                                           return(self$p_value)
                                         },
                                          
                                         # ---------------------------- 8. Plotting ----------------------------
                                         
                                         plot_results = function(title = NULL) {
                                           if (is.null(self$p_value)) stop("Run $run() first.")
                                           
                                           df_plot <- data.frame(
                                             stat = c(self$null_distribution, self$observed_distribution),
                                             type = factor(c(rep("Null Distribution (Permuted)", length(self$null_distribution)), 
                                                             rep("Observed Distribution (Resampled)", length(self$observed_distribution))),
                                                           levels = c("Null Distribution (Permuted)", "Observed Distribution (Resampled)"))
                                           )
                                           
                                           metric_name <- ifelse(self$metric == "ks", "KS Statistic", "Anderson-Darling Statistic")
                                           main_title <- ifelse(is.null(title), 
                                                                sprintf("Nested Permutation Test (%s)", toupper(self$metric)),
                                                                title)
                                           
                                           sub_title <- sprintf("p-value: %.5f", self$p_value)
                                           
                                           p <- ggplot(df_plot, aes(x = stat, fill = type)) +
                                             geom_density(alpha = 0.4, color = NA) +
                                             geom_vline(aes(xintercept = self$observed_stat_median, color = "Observed Median"), 
                                                        linetype = "dashed", linewidth = 1, show.legend = TRUE) +
                                             scale_fill_manual(values = c("Null Distribution (Permuted)" = "grey", 
                                                                          "Observed Distribution (Resampled)" = "blue")) +
                                             scale_color_manual(name = "", values = c("Observed Median" = "red")) +
                                        
                                             labs(title = main_title, 
                                                  subtitle = sub_title, 
                                                  x = metric_name, 
                                                  y = "Density", 
                                                  fill = "") +
                                             
                                             theme_minimal() +
                                             theme(
                                               legend.position = "bottom",
                                               
                                               plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                                               plot.subtitle = element_text(hjust = 0.5, size = 12),
                                               panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
                                               axis.line = element_blank()
                                             )
                                           
                                           print(p)
                                         }
                                       )
)

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
# We will use the same data (fMRI) obtained from seaborn (see notebook folder)

# Load the same .csv used in the .ipynb
dfmri_frontal <- read.csv("C:/Users/marce/Escritorio/Proyectos/hks/R/dfmri_frontal.csv")
dfmri_parietal <- read.csv("C:/Users/marce/Escritorio/Proyectos/hks/R/dfmri_parietal.csv")

tester_frontal <- HierarchicalPermutationTest$new(
  df = dfmri_frontal,
  value_col = "signal",
  subject_col = "unique_subject_id",
  cell_col = "recording_id",
  group_col = "event",
  metric = "ks",
  replace = TRUE
)

tester_parietal <- HierarchicalPermutationTest$new(
  df = dfmri_parietal,
  value_col = "signal",
  subject_col = "unique_subject_id",
  cell_col = "recording_id",
  group_col = "event",
  metric = "ks",
  replace = TRUE
)

p_val_frontal <- tester_frontal$run(n_resamples = 100000, n_permutations = 100000, seed = 123)
p_val_parietal <- tester_parietal$run(n_resamples = 100000, n_permutations = 100000, seed = 123)

print(paste("P-value (Frontal):", p_val_frontal))
print(paste("P-value (Parietal):", p_val_parietal))

tester_frontal$plot_results()
tester_parietal$plot_results()
