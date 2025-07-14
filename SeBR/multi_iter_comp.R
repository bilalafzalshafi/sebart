# Multi-iteration comparison script for semiparametric regression models

library(dbarts)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

source("helper_funs.R")
source("source_sba.R") 
source("sbart.R")
source("bbart_bc.R")
source("slice.R")
source("simulation_helpers.R")

# Set up simulation parameters
set.seed(123)
n_iterations <- 100  # Number of simulation iterations
n_train <- 200
n_test <- 1000
p <- 10

# Function to generate data for a given scenario
generate_data <- function(scenario_name, n_train, n_test, p) {
  simulate_sbart_data(n_train = n_train, n_test = n_test, p = p, 
                     scenario = scenario_name)
}

# Function to fit all models
fit_all_models <- function(data) {
  results <- list()
  timing <- list()
  
  X_train <- data$X_train
  y_train <- data$y_train
  X_test <- data$X_test
  y_test <- data$y_test
  
  # 1. Regular BART (no transformation)
  cat("Fitting BART...")
  timing$bart <- system.time({
    tryCatch({
      fit_bart <- bart(x.train = X_train, y.train = y_train, x.test = X_test, 
                       ntree = 200, ndpost = 1000, nskip = 1000, verbose = FALSE)
      results$bart <- list(
        fitted.values = fit_bart$yhat.test.mean,
        post_ypred = fit_bart$yhat.test,
        model = 'bart'
      )
    }, error = function(e) {
      cat("BART failed:", e$message, "\n")
      results$bart <<- NULL
    })
  })
  
  # 2. Semiparametric BART (sbart)
  cat("Fitting SBART...")
  timing$sbart <- system.time({
    tryCatch({
      fit_sbart <- sbart(y = y_train, X = X_train, X_test = X_test, 
                         ntree = 200, nsave = 1000, nburn = 1000, verbose = FALSE)
      results$sbart <- fit_sbart
    }, error = function(e) {
      cat("SBART failed:", e$message, "\n")
      results$sbart <<- NULL
    })
  })
  
  # 3. Semiparametric Bayesian Linear Model (sblm)
  cat("Fitting SBLM...")
  timing$sblm <- system.time({
    tryCatch({
      fit_sblm <- sblm(y = y_train, X = X_train, X_test = X_test)
      results$sblm <- fit_sblm
    }, error = function(e) {
      cat("SBLM failed:", e$message, "\n")
      results$sblm <<- NULL
    })
  })
  
  # 4. Bayesian BART with Box-Cox (bbart_bc)
  cat("Fitting BBART_BC...")
  timing$bbart_bc <- system.time({
    tryCatch({
      fit_bbart_bc <- bbart_bc(y = y_train, X = X_train, X_test = X_test,
                               ntree = 200, nsave = 1000, nburn = 1000, verbose = FALSE)
      results$bbart_bc <- fit_bbart_bc
    }, error = function(e) {
      cat("BBART_BC failed:", e$message, "\n")
      results$bbart_bc <<- NULL
    })
  })
  
  cat("Models fitted.\n")
  return(list(results = results, timing = timing))
}

# Function to calculate performance metrics
calculate_metrics <- function(results, y_test, alpha = 0.1) {
  metrics <- data.frame(
    Model = character(),
    RMSE = numeric(),
    Coverage = numeric(),
    Mean_Width = numeric(),
    Time_sec = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (model_name in names(results$results)) {
    model_result <- results$results[[model_name]]
    if (is.null(model_result)) next
    
    # RMSE
    rmse <- sqrt(mean((model_result$fitted.values - y_test)^2))
    
    # Coverage and interval width
    if (!is.null(model_result$post_ypred)) {
      lower <- apply(model_result$post_ypred, 2, quantile, alpha/2)
      upper <- apply(model_result$post_ypred, 2, quantile, 1 - alpha/2)
      coverage <- mean(y_test >= lower & y_test <= upper)
      mean_width <- mean(upper - lower)
    } else {
      coverage <- NA
      mean_width <- NA
    }
    
    # Timing
    time_sec <- results$timing[[model_name]][3]
    
    metrics <- rbind(metrics, data.frame(
      Model = model_name,
      RMSE = rmse,
      Coverage = coverage,
      Mean_Width = mean_width,
      Time_sec = time_sec,
      stringsAsFactors = FALSE
    ))
  }
  
  return(metrics)
}

# Main simulation loop
run_simulation_study <- function(scenario_name, n_iterations) {
  cat("Running simulation study for scenario:", scenario_name, "\n")
  cat("Number of iterations:", n_iterations, "\n\n")
  
  all_metrics <- data.frame()
  all_predictions <- list()
  
  for (i in 1:n_iterations) {
    if (i %% 10 == 0) cat("Iteration", i, "of", n_iterations, "\n")
    
    # Generate data
    data <- generate_data(scenario_name, n_train, n_test, p)
    
    # Fit models
    model_fits <- fit_all_models(data)
    
    # Calculate metrics
    metrics <- calculate_metrics(model_fits, data$y_test)
    metrics$Iteration <- i
    metrics$Scenario <- scenario_name
    
    all_metrics <- rbind(all_metrics, metrics)
    
    # Store predictions for plotting
    for (model_name in names(model_fits$results)) {
      if (!is.null(model_fits$results[[model_name]])) {
        if (is.null(all_predictions[[model_name]])) {
          all_predictions[[model_name]] <- list()
        }
        all_predictions[[model_name]][[i]] <- list(
          fitted = model_fits$results[[model_name]]$fitted.values,
          true = data$y_test,
          post_pred = model_fits$results[[model_name]]$post_ypred
        )
      }
    }
  }
  
  return(list(metrics = all_metrics, predictions = all_predictions))
}

# Function to create performance plots
create_performance_plots <- function(all_results) {
  # Combine all metrics
  combined_metrics <- do.call(rbind, lapply(all_results, function(x) x$metrics))
  
  # Remove failed runs
  combined_metrics <- combined_metrics[!is.na(combined_metrics$RMSE), ]
  
  # Coverage plot
  p_coverage <- ggplot(combined_metrics, aes(x = Model, y = Coverage, fill = Model)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
    facet_wrap(~Scenario, scales = "free") +
    labs(title = "90% Prediction Interval Coverage", 
         y = "Empirical Coverage") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Interval width plot
  p_width <- ggplot(combined_metrics, aes(x = Model, y = Mean_Width, fill = Model)) +
    geom_boxplot() +
    facet_wrap(~Scenario, scales = "free") +
    labs(title = "90% Prediction Interval Width", 
         y = "Mean Interval Width") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # RMSE plot
  p_rmse <- ggplot(combined_metrics, aes(x = Model, y = RMSE, fill = Model)) +
    geom_boxplot() +
    facet_wrap(~Scenario, scales = "free") +
    labs(title = "Root Mean Square Error", 
         y = "RMSE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Timing plot
  p_time <- ggplot(combined_metrics, aes(x = Model, y = Time_sec, fill = Model)) +
    geom_boxplot() +
    facet_wrap(~Scenario, scales = "free") +
    labs(title = "Computation Time", 
         y = "Time (seconds)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(coverage = p_coverage, width = p_width, rmse = p_rmse, time = p_time))
}

# Function to create prediction vs truth plots
create_prediction_plots <- function(all_results, scenario_name) {
  # Use first iteration for visualization
  predictions <- all_results[[scenario_name]]$predictions
  
  plots <- list()
  for (model_name in names(predictions)) {
    if (length(predictions[[model_name]]) > 0) {
      # Use first successful iteration
      pred_data <- predictions[[model_name]][[1]]
      
      plot_data <- data.frame(
        True = pred_data$true,
        Predicted = pred_data$fitted
      )
      
      p <- ggplot(plot_data, aes(x = True, y = Predicted)) +
        geom_point(alpha = 0.5) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        labs(title = paste(model_name, "-", scenario_name),
             x = "True Values", y = "Predicted Values") +
        theme_minimal()
      
      plots[[model_name]] <- p
    }
  }
  
  return(plots)
}

# Run the full simulation study
main_simulation <- function() {
  cat("Starting model comparison study...\n\n")
  
  # Run simulation for each scenario
  all_results <- list()
  
  # Start with a smaller set for initial testing
  scenarios_to_test <- c("identity", "box_cox", "step")
  
  for (scenario in scenarios_to_test) {
    all_results[[scenario]] <- run_simulation_study(scenario, n_iterations)
  }
  
  cat("\nGenerating plots...\n")
  
  # Create performance plots
  perf_plots <- create_performance_plots(all_results)
  
  # Save plots
  ggsave("coverage_comparison.png", perf_plots$coverage, width = 12, height = 8)
  ggsave("width_comparison.png", perf_plots$width, width = 12, height = 8)
  ggsave("rmse_comparison.png", perf_plots$rmse, width = 12, height = 8)
  ggsave("time_comparison.png", perf_plots$time, width = 12, height = 8)
  
  # Create prediction plots for each scenario
  for (scenario in scenarios_to_test) {
    pred_plots <- create_prediction_plots(all_results, scenario)
    
    if (length(pred_plots) > 0) {
      combined_plot <- do.call(grid.arrange, c(pred_plots, ncol = 2))
      ggsave(paste0("predictions_", scenario, ".png"), combined_plot, width = 12, height = 8)
    }
  }
  
  # Print summary statistics
  cat("\n=== SUMMARY STATISTICS ===\n")
  combined_metrics <- do.call(rbind, lapply(all_results, function(x) x$metrics))
  combined_metrics <- combined_metrics[!is.na(combined_metrics$RMSE), ]
  
  summary_stats <- combined_metrics %>%
    group_by(Model, Scenario) %>%
    summarise(
      Mean_RMSE = mean(RMSE, na.rm = TRUE),
      Median_Coverage = median(Coverage, na.rm = TRUE),
      Mean_Width = mean(Mean_Width, na.rm = TRUE),
      Mean_Time = mean(Time_sec, na.rm = TRUE),
      .groups = 'drop'
    )
  
  print(summary_stats)
  
  # Save results
  saveRDS(all_results, "simulation_results.rds")
  write.csv(combined_metrics, "detailed_metrics.csv", row.names = FALSE)
  write.csv(summary_stats, "summary_statistics.csv", row.names = FALSE)
  
  cat("\nSimulation completed. Results saved to files.\n")
  return(all_results)
}

# Run the simulation
if (interactive()) {
  cat("Starting simulation study...\n")
  cat("This may take some time depending on n_iterations =", n_iterations, "\n")
  results <- main_simulation()
} else {
  cat("Source this file and run main_simulation() to start the study.\n")
}