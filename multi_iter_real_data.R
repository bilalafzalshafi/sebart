# Multi-iteration model comparison for real datasets

source("sebart.R")
source("bart_bc.R")
source("real_data_helpers.R")

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

n_iterations <- 100  # Use 10 iterations for publication-quality results
                   # Set to 2-3 for testing/debugging

# Calculate Continuous Ranked Probability Score
# CRPS = E|Y' - y| - 0.5 * E|Y' - Y''|
# where Y' and Y'' are independent samples from predictive distribution
# calculate_crps <- function(post_samples, y_true) {
#   n_samples <- nrow(post_samples)
#   n_obs <- length(y_true)
#   
#   # First term: E|Y' - y|
#   first_term <- mean(abs(post_samples - matrix(y_true, nrow = n_samples, 
#                                                 ncol = n_obs, byrow = TRUE)))
#   
#   # Second term: 0.5 * E|Y' - Y''|
#   # Approximate by comparing random pairs of posterior samples
#   n_pairs <- min(1000, n_samples)
#   idx1 <- sample(n_samples, n_pairs, replace = TRUE)
#   idx2 <- sample(n_samples, n_pairs, replace = TRUE)
#   
#   second_term <- 0
#   for (i in 1:n_obs) {
#     second_term <- second_term + mean(abs(post_samples[idx1, i] - post_samples[idx2, i]))
#   }
#   second_term <- 0.5 * second_term / n_obs
#   
#   return(first_term - second_term)
# }

run_real_data_iteration <- function(dataset_name, target_feature, 
                                   models_to_test = c("sebart", "bart", "bart_bc", "sblm"),
                                   iteration_seed = 123,
                                   alpha_level = 0.10) {
  
  real_data <- load_real_dataset(dataset_name, target_feature, seed = iteration_seed)
  
  results <- list()
  
  for (model_name in models_to_test) {
    
    # start_time <- Sys.time()

    if (model_name == "sebart") {
      fit <- sebart(y = real_data$y_train, X = real_data$X_train, 
                    X_test = real_data$X_test, ntree = 200, nsave = 1000, nburn = 1000,
                    verbose = FALSE)
      y_pred <- apply(fit$post_ypred, 2, mean)
      pred_intervals <- apply(fit$post_ypred, 2, quantile, c(alpha_level/2, 1-alpha_level/2))
      # post_samples <- fit$post_ypred
      
    } else if (model_name == "bart_bc") {
      fit <- bart_bc(y = real_data$y_train, X = real_data$X_train, 
                      X_test = real_data$X_test, ntree = 200, nsave = 1000, nburn = 1000,
                      verbose = FALSE)
      y_pred <- apply(fit$post_ypred, 2, mean)
      pred_intervals <- apply(fit$post_ypred, 2, quantile, c(alpha_level/2, 1-alpha_level/2))
      # post_samples <- fit$post_ypred
      
    } else if (model_name == "bart") {
      library(dbarts)
      fit <- bart(x.train = real_data$X_train, y.train = real_data$y_train,
                  x.test = real_data$X_test, ndpost = 1000, nskip = 1000,
                  ntree = 200, verbose = FALSE)
      y_pred <- apply(fit$yhat.test, 2, mean)
      pred_intervals <- apply(fit$yhat.test, 2, quantile, c(alpha_level/2, 1-alpha_level/2))
      # post_samples <- fit$yhat.test
      
    } else if (model_name == "sblm") {
      fit <- SeBR::sblm(y = real_data$y_train, X = real_data$X_train, 
                  X_test = real_data$X_test)
      y_pred <- apply(fit$post_ypred, 2, mean)
      pred_intervals <- apply(fit$post_ypred, 2, quantile, c(alpha_level/2, 1-alpha_level/2))
      # post_samples <- fit$post_ypred
    }
    
    # end_time <- Sys.time()
    
    # Calculate metrics
    rmse <- sqrt(mean((y_pred - real_data$y_test)^2))
    
    coverage_90 <- mean(real_data$y_test >= pred_intervals[1, ] & 
                        real_data$y_test <= pred_intervals[2, ])
    
    width_90 <- mean(pred_intervals[2, ] - pred_intervals[1, ])
    
    # Calculate CRPS
    # crps <- calculate_crps(post_samples, real_data$y_test)
    
    results[[model_name]] <- list(
      RMSE = rmse,
      Coverage_90 = coverage_90, 
      Width_90 = width_90,
      # CRPS = crps,
      # Time_sec = as.numeric(difftime(end_time, start_time, units = "secs")),
      Model = model_name,
      Dataset = dataset_name,
      Target = target_feature,
      Iteration = iteration_seed,
      Predictions = y_pred
    )
  }
  
  return(results)
}

run_real_data_study <- function(dataset_name, target_feature, n_iterations = 10) {
  study_results <- list()
  
  for (i in 1:n_iterations) {    
    iteration_results <- run_real_data_iteration(dataset_name, target_feature, 
                                                iteration_seed = i)
    study_results[[i]] <- iteration_results
  }
  
  all_metrics <- do.call(rbind, lapply(study_results, function(iter) {
    do.call(rbind, lapply(iter, function(model_result) {
      data.frame(
        Model = model_result$Model,
        Dataset = model_result$Dataset,
        Target = model_result$Target,
        Iteration = model_result$Iteration,
        RMSE = model_result$RMSE,
        Coverage_90 = model_result$Coverage_90,
        Width_90 = model_result$Width_90,
        # CRPS = model_result$CRPS,
        # Time_sec = model_result$Time_sec,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  return(list(
    iterations = study_results,
    metrics = all_metrics,
    dataset_name = dataset_name,
    target_feature = target_feature
  ))
}

create_real_data_performance_plots <- function(all_results) {
  
  combined_metrics <- do.call(rbind, lapply(all_results, function(x) x$metrics))
  combined_metrics <- combined_metrics[!is.na(combined_metrics$RMSE), ]
  
  if (nrow(combined_metrics) == 0) {
    warning("No valid results to plot")
    return(NULL)
  }
  
  # Ensure correct model order
  combined_metrics$Model <- factor(combined_metrics$Model, 
                                  levels = c("sebart", "bart", "bart_bc", "sblm"))
  
  # Define professional theme for all plots
  publication_theme <- theme_bw() +
    theme(
      panel.grid.major = element_line(color = "gray90", size = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 0.5),
      strip.background = element_rect(fill = "white", color = "black", size = 0.5),
      strip.text = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
      legend.position = "none"
    )
  
  # Single horizontal box plot for interval widths (averaged across all datasets)
  width_coverage_data <- combined_metrics %>%
    group_by(Model) %>%
    summarise(
      Mean_Coverage = mean(Coverage_90),
      .groups = 'drop'
    )
  
  # Calculate position for coverage labels
  max_width <- max(combined_metrics$Width_90) * 1.15
  
  width_plot <- ggplot(combined_metrics, aes(x = Width_90, y = Model)) +
    geom_boxplot(fill = "white", color = "black", outlier.size = 1, 
                 outlier.shape = 1, size = 0.5) +
    geom_text(data = width_coverage_data, 
              aes(x = max_width, y = Model, 
                  label = paste0(round(Mean_Coverage * 100), "%")),
              size = 3.5, color = "blue", fontface = "bold", hjust = 0) +
    labs(title = "90% Prediction Interval Widths",
         x = "Interval Width") +
    publication_theme +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.20)))
  
  # Calculate RMSE summary table
  rmse_summary <- combined_metrics %>%
    group_by(Model, Dataset, Target) %>%
    summarise(
      Mean_RMSE = mean(RMSE, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(Dataset_Target = paste(Dataset, Target, sep = "_")) %>%
    select(-Dataset, -Target) %>%
    pivot_wider(names_from = Dataset_Target, values_from = Mean_RMSE) %>%
    mutate(Average = round(rowMeans(select(., -Model), na.rm = TRUE), 3))
  
  # Round all numeric columns
  rmse_summary <- rmse_summary %>%
    mutate(across(where(is.numeric), ~round(., 3)))
  
  return(list(
    width_plot = width_plot,
    rmse_table = rmse_summary
  ))
}

main_real_data_study <- function() {
  all_results <- list()
  
  datasets_to_test <- list(
    "forestfires" = c("area"),
    "bikesharing" = c("hum", "windspeed"),  # removed "temp"
    "parkinsons" = c("NHR", "DFA")
  )
  
  for (dataset_name in names(datasets_to_test)) {
    available_targets <- get_target_features(dataset_name)
    
    if (is.null(available_targets)) next
    
    targets_to_test <- intersect(datasets_to_test[[dataset_name]], available_targets)
    
    for (target_feature in targets_to_test) {
      study_name <- paste(dataset_name, target_feature, sep = "_")
      
      all_results[[study_name]] <- run_real_data_study(dataset_name, target_feature, n_iterations)
    }
  }
    
  results <- create_real_data_performance_plots(all_results)
  
  if (!is.null(results)) {
    # Save width plot
    ggsave("real_data_width_90.png", results$width_plot, 
           width = 8, height = 5, dpi = 300, bg = "white")
    
    # Print RMSE table
    cat("\nMean RMSE by Model and Dataset\n")
    cat("================================\n")
    print(as.data.frame(results$rmse_table), row.names = FALSE)
    
    # Save RMSE table
    write.csv(results$rmse_table, "rmse_summary_table.csv", row.names = FALSE)
  }
  
  # Summary statistics
  combined_metrics <- do.call(rbind, lapply(all_results, function(x) x$metrics))
  combined_metrics <- combined_metrics[!is.na(combined_metrics$RMSE), ]
  
  if (nrow(combined_metrics) > 0) {
    # Simple summary by model
    summary_by_model <- combined_metrics %>%
      group_by(Model) %>%
      summarise(
        RMSE = round(mean(RMSE, na.rm = TRUE), 3),
        Coverage = round(mean(Coverage_90, na.rm = TRUE) * 100, 1),
        Width = round(mean(Width_90, na.rm = TRUE), 3),
        .groups = 'drop'
      ) %>%
      arrange(match(Model, c("sebart", "bart", "bart_bc", "sblm")))
    
    cat("\n\nOverall Model Performance\n")
    cat("-------------------------\n")
    print(as.data.frame(summary_by_model), row.names = FALSE)
    
    # Save results
    write.csv(combined_metrics, "results_raw.csv", row.names = FALSE)
    write.csv(summary_by_model, "results_summary.csv", row.names = FALSE)
  }
  
  saveRDS(all_results, "real_data_results.rds")
  return(all_results)
}