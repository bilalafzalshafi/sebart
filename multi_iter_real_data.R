# Multi-iteration model comparison for real datasets

source("sbart.R")
source("bbart_bc.R")
source("real_data_helpers.R")  # Loads real datasets

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Set parameters
n_iterations <- 100

#' Run single iteration on real data
run_real_data_iteration <- function(dataset_name, target_feature, 
                                   models_to_test = c("sbart", "bbart_bc", "bart", "sblm"),
                                   iteration_seed = 123,
                                   alpha_level = 0.10) {
  
  # Load data for this iteration
  real_data <- load_real_dataset(dataset_name, target_feature, seed = iteration_seed)
  
  results <- list()
  
  for (model_name in models_to_test) {
    cat("  Testing", model_name, "...\n")
    
    start_time <- Sys.time()
    
    tryCatch({
      if (model_name == "sbart") {
        fit <- sbart(y = real_data$y_train, X = real_data$X_train, 
                     X_test = real_data$X_test, ntree = 200, nsave = 1000, nburn = 1000,
                     verbose = FALSE)
        y_pred <- apply(fit$post_ypred, 2, mean)
        pred_intervals_80 <- apply(fit$post_ypred, 2, quantile, c(0.10, 0.90))
        pred_intervals <- apply(fit$post_ypred, 2, quantile, c(alpha_level/2, 1-alpha_level/2))
        pred_intervals_95 <- apply(fit$post_ypred, 2, quantile, c(0.025, 0.975))
        
      } else if (model_name == "bbart_bc") {
        fit <- bbart_bc(y = real_data$y_train, X = real_data$X_train, 
                        X_test = real_data$X_test, ntree = 200, nsave = 1000, nburn = 1000,
                        verbose = FALSE)
        y_pred <- apply(fit$post_ypred, 2, mean)
        pred_intervals_80 <- apply(fit$post_ypred, 2, quantile, c(0.10, 0.90))
        pred_intervals <- apply(fit$post_ypred, 2, quantile, c(alpha_level/2, 1-alpha_level/2))
        pred_intervals_95 <- apply(fit$post_ypred, 2, quantile, c(0.025, 0.975))
        
      } else if (model_name == "bart") {
        library(dbarts)
        fit <- bart(x.train = real_data$X_train, y.train = real_data$y_train,
                    x.test = real_data$X_test, ndpost = 1000, nskip = 1000,
                    ntree = 200, verbose = FALSE)
        y_pred <- apply(fit$yhat.test, 2, mean)
        pred_intervals_80 <- apply(fit$yhat.test, 2, quantile, c(0.10, 0.90))
        pred_intervals <- apply(fit$yhat.test, 2, quantile, c(alpha_level/2, 1-alpha_level/2))
        pred_intervals_95 <- apply(fit$yhat.test, 2, quantile, c(0.025, 0.975))
        
      } else if (model_name == "sblm") {
        fit <- sblm(y = real_data$y_train, X = real_data$X_train, 
                    X_test = real_data$X_test)
        y_pred <- apply(fit$post_ypred, 2, mean)
        pred_intervals_80 <- apply(fit$post_ypred, 2, quantile, c(0.10, 0.90))
        pred_intervals <- apply(fit$post_ypred, 2, quantile, c(alpha_level/2, 1-alpha_level/2))
        pred_intervals_95 <- apply(fit$post_ypred, 2, quantile, c(0.025, 0.975))
      }
      
      end_time <- Sys.time()
      
      # Calculate performance metrics
      rmse <- sqrt(mean((y_pred - real_data$y_test)^2))
      
      # Coverage rates for multiple interval levels
      coverage_80 <- mean(real_data$y_test >= pred_intervals_80[1, ] & 
                         real_data$y_test <= pred_intervals_80[2, ])
      coverage_90 <- mean(real_data$y_test >= pred_intervals[1, ] & 
                         real_data$y_test <= pred_intervals[2, ])
      coverage_95 <- mean(real_data$y_test >= pred_intervals_95[1, ] & 
                         real_data$y_test <= pred_intervals_95[2, ])
      
      # Average interval widths
      width_80 <- mean(pred_intervals_80[2, ] - pred_intervals_80[1, ])
      width_90 <- mean(pred_intervals[2, ] - pred_intervals[1, ])
      width_95 <- mean(pred_intervals_95[2, ] - pred_intervals_95[1, ])
      
      results[[model_name]] <- list(
        RMSE = rmse,
        Coverage_80 = coverage_80,
        Coverage_90 = coverage_90, 
        Coverage_95 = coverage_95,
        Width_80 = width_80,
        Width_90 = width_90,
        Width_95 = width_95,
        Time_sec = as.numeric(difftime(end_time, start_time, units = "secs")),
        Model = model_name,
        Dataset = dataset_name,
        Target = target_feature,
        Iteration = iteration_seed,
        Predictions = y_pred  # Store for plotting
      )
      
    }, error = function(e) {
      cat("    Error in", model_name, ":", e$message, "\n")
      results[[model_name]] <- list(
        RMSE = NA, Coverage_80 = NA, Coverage_90 = NA, Coverage_95 = NA,
        Width_80 = NA, Width_90 = NA, Width_95 = NA, Time_sec = NA,
        Model = model_name, Dataset = dataset_name, Target = target_feature, 
        Iteration = iteration_seed, Predictions = NULL
      )
    })
  }
  
  return(results)
}

#' Run study for single dataset-target combination
run_real_data_study <- function(dataset_name, target_feature, n_iterations = 10) {
  
  cat("Running study for", dataset_name, "-", target_feature, "\n")
  cat("Iterations:", n_iterations, "\n")
  
  study_results <- list()
  
  for (i in 1:n_iterations) {
    cat("Iteration", i, "of", n_iterations, "\n")
    
    iteration_results <- run_real_data_iteration(dataset_name, target_feature, 
                                                iteration_seed = i)
    study_results[[i]] <- iteration_results
  }
  
  # Combine results into data frame for analysis
  all_metrics <- do.call(rbind, lapply(study_results, function(iter) {
    do.call(rbind, lapply(iter, function(model_result) {
      data.frame(
        Model = model_result$Model,
        Dataset = model_result$Dataset,
        Target = model_result$Target,
        Iteration = model_result$Iteration,
        RMSE = model_result$RMSE,
        Coverage_80 = model_result$Coverage_80,
        Coverage_90 = model_result$Coverage_90,
        Coverage_95 = model_result$Coverage_95,
        Width_80 = model_result$Width_80,
        Width_90 = model_result$Width_90,
        Width_95 = model_result$Width_95,
        Time_sec = model_result$Time_sec,
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

#' Create performance plots for real data results
create_real_data_performance_plots <- function(all_results) {
  
  # Combine all metrics
  combined_metrics <- do.call(rbind, lapply(all_results, function(x) x$metrics))
  combined_metrics <- combined_metrics[!is.na(combined_metrics$RMSE), ]
  
  if (nrow(combined_metrics) == 0) {
    warning("No valid results to plot")
    return(NULL)
  }
  
  # Create study labels combining dataset and target
  combined_metrics$Study <- paste(combined_metrics$Dataset, combined_metrics$Target, sep = "_")
  
  # Coverage plots for multiple levels
  coverage_90_plot <- ggplot(combined_metrics, aes(x = Model, y = Coverage_90, fill = Model)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
    facet_wrap(~ Study, scales = "free_x") +
    labs(title = "90% Coverage Comparison", 
         subtitle = "Red line = 90% target",
         y = "Coverage Rate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  coverage_80_plot <- ggplot(combined_metrics, aes(x = Model, y = Coverage_80, fill = Model)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    facet_wrap(~ Study, scales = "free_x") +
    labs(title = "80% Coverage Comparison", 
         subtitle = "Red line = 80% target",
         y = "Coverage Rate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  coverage_95_plot <- ggplot(combined_metrics, aes(x = Model, y = Coverage_95, fill = Model)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    facet_wrap(~ Study, scales = "free_x") +
    labs(title = "95% Coverage Comparison", 
         subtitle = "Red line = 95% target",
         y = "Coverage Rate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # RMSE plot  
  rmse_plot <- ggplot(combined_metrics, aes(x = Model, y = RMSE, fill = Model)) +
    geom_boxplot() +
    facet_wrap(~ Study, scales = "free") +
    labs(title = "RMSE Comparison", y = "RMSE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Interval width plots for multiple levels
  width_90_plot <- ggplot(combined_metrics, aes(x = Model, y = Width_90, fill = Model)) +
    geom_boxplot() +
    facet_wrap(~ Study, scales = "free") +
    labs(title = "90% Interval Width Comparison", y = "Mean Interval Width") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  width_comparison_plot <- combined_metrics %>%
    select(Model, Study, Width_80, Width_90, Width_95) %>%
    gather(key = "Interval_Level", value = "Width", -Model, -Study) %>%
    mutate(Interval_Level = factor(Interval_Level, levels = c("Width_80", "Width_90", "Width_95"),
                                  labels = c("80%", "90%", "95%"))) %>%
    ggplot(aes(x = Model, y = Width, fill = Interval_Level)) +
    geom_boxplot() +
    facet_wrap(~ Study, scales = "free") +
    labs(title = "Interval Width by Coverage Level", 
         y = "Mean Interval Width", fill = "Interval Level") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Time plot
  time_plot <- ggplot(combined_metrics, aes(x = Model, y = Time_sec, fill = Model)) +
    geom_boxplot() +
    facet_wrap(~ Study, scales = "free") +
    labs(title = "Computation Time Comparison", y = "Time (seconds)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(
    coverage_80 = coverage_80_plot,
    coverage_90 = coverage_90_plot,
    coverage_95 = coverage_95_plot,
    rmse = rmse_plot,
    width_90 = width_90_plot,
    width_comparison = width_comparison_plot,
    time = time_plot
  ))
}

#' Create prediction scatter plots for real data
create_real_data_prediction_plots <- function(all_results, study_name) {
  
  study_results <- all_results[[study_name]]
  if (is.null(study_results)) {
    warning("Study not found:", study_name)
    return(NULL)
  }
  
  # Use first iteration for plotting
  first_iter <- study_results$iterations[[1]]
  
  # Load the same data to get true y values
  parts <- strsplit(study_name, "_")[[1]]
  dataset_name <- parts[1]
  target_feature <- paste(parts[-1], collapse = "_")
  
  real_data <- load_real_dataset(dataset_name, target_feature, seed = 1)
  
  plots <- list()
  
  for (model_name in names(first_iter)) {
    if (!is.null(first_iter[[model_name]]$Predictions)) {
      
      plot_data <- data.frame(
        y_true = real_data$y_test,
        y_pred = first_iter[[model_name]]$Predictions
      )
      
      p <- ggplot(plot_data, aes(x = y_true, y = y_pred)) +
        geom_point(alpha = 0.5) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        labs(title = paste(model_name, "-", study_name),
             x = "True values", y = "Predicted values") +
        theme_minimal()
      
      plots[[model_name]] <- p
    }
  }
  
  return(plots)
}

#' Main real data study function
main_real_data_study <- function() {
  cat("Starting real data model comparison study...\n\n")
  
  all_results <- list()
  
  # Define datasets and targets to test
  datasets_to_test <- list(
    "forestfires" = c("area"),
    "bikesharing" = c("temp", "hum", "windspeed"),
    "concrete" = c("strength"),
    "glass" = c("Na", "Al", "Ca"),
    "studentperformance" = c("G3_norm"),
    "parkinsons" = c("NHR", "DFA"),
    "yeast" = c("mcg", "alm")
  )
  
  for (dataset_name in names(datasets_to_test)) {
    available_targets <- tryCatch({
      get_target_features(dataset_name)
    }, error = function(e) {
      cat("Skipping", dataset_name, "- not available\n")
      return(NULL)
    })
    
    if (is.null(available_targets)) next
    
    targets_to_test <- intersect(datasets_to_test[[dataset_name]], available_targets)
    
    for (target_feature in targets_to_test) {
      study_name <- paste(dataset_name, target_feature, sep = "_")
      
      tryCatch({
        all_results[[study_name]] <- run_real_data_study(dataset_name, target_feature, n_iterations)
      }, error = function(e) {
        cat("Error in study", study_name, ":", e$message, "\n")
      })
    }
  }
  
  if (length(all_results) == 0) {
    cat("No successful studies completed.\n")
    return(NULL)
  }
  
  cat("\nGenerating plots...\n")
  
  # Performance plots
  perf_plots <- create_real_data_performance_plots(all_results)
  
  if (!is.null(perf_plots)) {
    ggsave("real_data_coverage_80.png", perf_plots$coverage_80, width = 12, height = 8)
    ggsave("real_data_coverage_90.png", perf_plots$coverage_90, width = 12, height = 8)
    ggsave("real_data_coverage_95.png", perf_plots$coverage_95, width = 12, height = 8)
    ggsave("real_data_rmse.png", perf_plots$rmse, width = 12, height = 8)
    ggsave("real_data_width_90.png", perf_plots$width_90, width = 12, height = 8)
    ggsave("real_data_width_comparison.png", perf_plots$width_comparison, width = 12, height = 8)
    ggsave("real_data_time.png", perf_plots$time, width = 12, height = 8)
  }
  
  # # Prediction plots for each study
  # for (study_name in names(all_results)) {
  #   pred_plots <- create_real_data_prediction_plots(all_results, study_name)
    
  #   if (length(pred_plots) > 0) {
  #     combined_plot <- do.call(grid.arrange, c(pred_plots, ncol = 2))
  #     ggsave(paste0("real_data_predictions_", study_name, ".png"), 
  #            combined_plot, width = 12, height = 8)
  #   }
  # }
  
  cat("\n=== REAL DATA SUMMARY STATISTICS ===\n")
  combined_metrics <- do.call(rbind, lapply(all_results, function(x) x$metrics))
  combined_metrics <- combined_metrics[!is.na(combined_metrics$RMSE), ]
  
  if (nrow(combined_metrics) > 0) {
    summary_stats <- combined_metrics %>%
      group_by(Model, Dataset, Target) %>%
      summarise(
        Mean_RMSE = mean(RMSE, na.rm = TRUE),
        Coverage_80 = mean(Coverage_80, na.rm = TRUE),
        Coverage_90 = mean(Coverage_90, na.rm = TRUE), 
        Coverage_95 = mean(Coverage_95, na.rm = TRUE),
        Width_80 = mean(Width_80, na.rm = TRUE),
        Width_90 = mean(Width_90, na.rm = TRUE),
        Width_95 = mean(Width_95, na.rm = TRUE),
        Mean_Time = mean(Time_sec, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Create formatted summary table
    formatted_summary <- combined_metrics %>%
      group_by(Model, Dataset, Target) %>%
      summarise(
        RMSE = mean(RMSE, na.rm = TRUE),
        "80%" = paste0(round(mean(Width_80, na.rm = TRUE), 3), " (", 
                      round(mean(Coverage_80, na.rm = TRUE), 2), ")"),
        "90%" = paste0(round(mean(Width_90, na.rm = TRUE), 3), " (", 
                      round(mean(Coverage_90, na.rm = TRUE), 2), ")"),
        "95%" = paste0(round(mean(Width_95, na.rm = TRUE), 3), " (", 
                      round(mean(Coverage_95, na.rm = TRUE), 2), ")"),
        Time_sec = mean(Time_sec, na.rm = TRUE),
        .groups = 'drop'
      )
    
    print("=== DETAILED SUMMARY ===")
    print("Format: Width (Coverage) for each interval level")
    print(formatted_summary)
    
    print("\n=== STANDARD SUMMARY ===")
    print(summary_stats)
    
    write.csv(combined_metrics, "real_data_detailed_metrics.csv", row.names = FALSE)
    write.csv(summary_stats, "real_data_summary_statistics.csv", row.names = FALSE)
    write.csv(formatted_summary, "real_data_formatted_summary.csv", row.names = FALSE)
  }
  
  saveRDS(all_results, "real_data_results.rds")
  
  cat("\nReal data study completed. Results saved to files.\n")
  return(all_results)
}

# Interactive usage
if (interactive()) {
  cat("Starting real data study...\n")
  cat("This may take some time depending on n_iterations =", n_iterations, "\n")
  cat("Make sure you have downloaded the required datasets.\n\n")
  
  # You can also run individual studies:
  # single_result <- run_real_data_study("wine_quality", "alcohol", 5)
  all_results <- main_real_data_study()
} else {
  cat("Source this file and run main_real_data_study() to start the study.\n")
}