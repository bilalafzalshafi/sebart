# Multi-iteration comparison script

library(SeBR)
library(dbarts)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

source("sebart.R")
source("bbart_bc.R")
source("simulation_helpers.R")

set.seed(123)
n_iterations <- 100
n_train <- 200
n_test <- 1000
p <- 10

generate_data <- function(scenario_name, n_train, n_test, p) {
  simulate_sebart_data(n_train = n_train, n_test = n_test, p = p, 
                     scenario = scenario_name)
}

fit_all_models <- function(data) {
  results <- list()
  timing <- list()
  
  X_train <- data$X_train
  y_train <- data$y_train
  X_test <- data$X_test
  y_test <- data$y_test
  
  timing$bart <- system.time({
      fit_bart <- dbarts::bart(x.train = X_train, y.train = y_train, x.test = X_test, 
                       ntree = 200, ndpost = 1000, nskip = 1000, verbose = FALSE)
      results$bart <- list(
        fitted.values = fit_bart$yhat.test.mean,
        post_ypred = fit_bart$yhat.test,
        model = 'bart'
      )
  })
  
  timing$sebart <- system.time({
    fit_sebart <- sebart(y = y_train, X = X_train, X_test = X_test, 
                        ntree = 200, nsave = 1000, nburn = 1000, verbose = FALSE)
    results$sebart <- fit_sebart
  })
  
  timing$sblm <- system.time({
    fit_sblm <- SeBR::sblm(y = y_train, X = X_train, X_test = X_test)
    results$sblm <- fit_sblm
  })
  
  timing$bbart_bc <- system.time({
      fit_bbart_bc <- bbart_bc(y = y_train, X = X_train, X_test = X_test,
                               ntree = 200, nsave = 1000, nburn = 1000, verbose = FALSE)
      results$bbart_bc <- fit_bbart_bc
  })
  
  return(list(results = results, timing = timing))
}

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
    
    rmse <- sqrt(mean((model_result$fitted.values - y_test)^2))
    
    if (!is.null(model_result$post_ypred)) {
      lower <- apply(model_result$post_ypred, 2, quantile, alpha/2)
      upper <- apply(model_result$post_ypred, 2, quantile, 1 - alpha/2)
      coverage <- mean(y_test >= lower & y_test <= upper)
      mean_width <- mean(upper - lower)
    } else {
      coverage <- NA
      mean_width <- NA
    }
    
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

run_simulation_study <- function(scenario_name, n_iterations) {
  all_metrics <- data.frame()
  all_predictions <- list()
  
  for (i in 1:n_iterations) {
    data <- generate_data(scenario_name, n_train, n_test, p)
    
    model_fits <- fit_all_models(data)
    
    metrics <- calculate_metrics(model_fits, data$y_test)
    metrics$Iteration <- i
    metrics$Scenario <- scenario_name
    
    all_metrics <- rbind(all_metrics, metrics)
    
    # for (model_name in names(model_fits$results)) {
    #   if (!is.null(model_fits$results[[model_name]])) {
    #     if (is.null(all_predictions[[model_name]])) {
    #       all_predictions[[model_name]] <- list()
    #     }
    #     all_predictions[[model_name]][[i]] <- list(
    #       fitted = model_fits$results[[model_name]]$fitted.values,
    #       true = data$y_test,
    #       post_pred = model_fits$results[[model_name]]$post_ypred
    #     )
    #   }
    # }
  }
  
  return(list(metrics = all_metrics, predictions = all_predictions))
}

create_width_plot <- function(all_results) {
  combined_metrics <- do.call(rbind, lapply(all_results, function(x) x$metrics))
  
  combined_metrics <- combined_metrics[!is.na(combined_metrics$RMSE), ]
  
  # Filter for sigmoid scenario only
  combined_metrics <- combined_metrics[combined_metrics$Scenario == "sigmoid", ]
  
  # Ensure correct model order
  combined_metrics$Model <- factor(combined_metrics$Model, 
                                  levels = c("sebart", "bart", "bbart_bc", "sblm"))
  
  # Define professional theme (matching real data version)
  publication_theme <- theme_bw() +
    theme(
      panel.grid.major = element_line(color = "gray90", size = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 0.5),
      strip.background = element_rect(fill = "white", color = "black", size = 0.5),
      strip.text = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
      legend.position = "none"
    )
  
  # Calculate coverage percentages for annotation
  coverage_stats <- combined_metrics %>%
    group_by(Model) %>%
    summarise(
      mean_coverage = mean(Coverage, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Calculate max width for positioning labels
  max_width <- max(combined_metrics$Mean_Width) * 1.15
  
  p_width <- ggplot(combined_metrics, aes(x = Mean_Width, y = Model)) +
    geom_boxplot(fill = "white", color = "black", outlier.size = 1, 
                 outlier.shape = 1, size = 0.5) +
    geom_text(data = coverage_stats, 
              aes(x = max_width, y = Model,
                  label = paste0(round(mean_coverage * 100), "%")),
              size = 3.5, color = "blue", fontface = "bold", hjust = 0) +
    labs(title = "90% Prediction Interval Widths (Sigmoid Scenario)",
         x = "Interval Width") +
    publication_theme +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.20)))
  
  return(p_width)
}

# create_performance_plots <- function(all_results) {
#   combined_metrics <- do.call(rbind, lapply(all_results, function(x) x$metrics))
#   
#   combined_metrics <- combined_metrics[!is.na(combined_metrics$RMSE), ]
#   
#   p_coverage <- ggplot(combined_metrics, aes(y = Model, 
#                                             x = Coverage)) +
#     geom_boxplot(fill = "lightgray", color = "black", width = 0.4, 
#                 outlier.shape = 16, outlier.size = 1.5) +
#     geom_vline(xintercept = 0.9, linetype = "dashed", color = "red", size = 0.8) +
#     stat_summary(fun = median, geom = "text", 
#                 aes(label = paste0(round(after_stat(x) * 100), "%")),
#                 hjust = -0.2, color = "blue", size = 4, fontface = "bold") +
#     facet_wrap(~Scenario, scales = "free_x", ncol = 3) +
#     scale_x_continuous(expand = c(0.1, 0)) +
#     labs(title = "90% prediction interval coverage", 
#         x = "Empirical coverage", y = NULL) +
#     theme_minimal() +
#     theme(
#       plot.background = element_rect(fill = "white", color = NA),
#       panel.background = element_rect(fill = "white", color = NA),
#       legend.position = "none",
#       panel.grid.major.x = element_line(color = "gray90", size = 0.3),
#       panel.grid.major.y = element_blank(),
#       panel.grid.minor = element_blank(),
#       strip.background = element_rect(fill = "gray95", color = "black", size = 0.3),
#       strip.text = element_text(face = "bold", size = 10),
#       axis.text.y = element_text(size = 10),
#       axis.text.x = element_text(size = 9),
#       plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#       panel.spacing.x = unit(0.5, "cm")
#     )
#   
#   p_rmse <- ggplot(combined_metrics, aes(y = Model, 
#                                         x = RMSE)) +
#     geom_boxplot(fill = "lightgray", color = "black", width = 0.4, 
#                 outlier.shape = 16, outlier.size = 1.5) +
#     facet_wrap(~Scenario, scales = "free_x", ncol = 3) +
#     scale_x_continuous(expand = c(0.1, 0)) +
#     labs(title = "Root mean square error", 
#         x = "RMSE", y = NULL) +
#     theme_minimal() +
#     theme(
#       plot.background = element_rect(fill = "white", color = NA),
#       panel.background = element_rect(fill = "white", color = NA),
#       legend.position = "none",
#       panel.grid.major.x = element_line(color = "gray90", size = 0.3),
#       panel.grid.major.y = element_blank(),
#       panel.grid.minor = element_blank(),
#       strip.background = element_rect(fill = "gray95", color = "black", size = 0.3),
#       strip.text = element_text(face = "bold", size = 10),
#       axis.text.y = element_text(size = 10),
#       axis.text.x = element_text(size = 9),
#       plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#       panel.spacing.x = unit(0.5, "cm")
#     )
#   
#   p_time <- ggplot(combined_metrics, aes(y = Model, 
#                                         x = Time_sec)) +
#     geom_boxplot(fill = "lightgray", color = "black", width = 0.4, 
#                 outlier.shape = 16, outlier.size = 1.5) +
#     facet_wrap(~Scenario, scales = "free_x", ncol = 3) +
#     scale_x_continuous(expand = c(0.1, 0)) +
#     labs(title = "Computation time", 
#         x = "Time (seconds)", y = NULL) +
#     theme_minimal() +
#     theme(
#       plot.background = element_rect(fill = "white", color = NA),
#       panel.background = element_rect(fill = "white", color = NA),
#       legend.position = "none",
#       panel.grid.major.x = element_line(color = "gray90", size = 0.3),
#       panel.grid.major.y = element_blank(),
#       panel.grid.minor = element_blank(),
#       strip.background = element_rect(fill = "gray95", color = "black", size = 0.3),
#       strip.text = element_text(face = "bold", size = 10),
#       axis.text.y = element_text(size = 10), 
#       axis.text.x = element_text(size = 9),
#       plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#       panel.spacing.x = unit(0.5, "cm")
#     )
#   
#   return(list(coverage = p_coverage, width = p_width, rmse = p_rmse, time = p_time))
# }

# create_prediction_plots <- function(all_results, scenario_name) {
#   # Use first iteration for visualization
#   predictions <- all_results[[scenario_name]]$predictions
#   
#   plots <- list()
#   for (model_name in names(predictions)) {
#     if (length(predictions[[model_name]]) > 0) {
#       # Use first successful iteration
#       pred_data <- predictions[[model_name]][[1]]
#       
#       plot_data <- data.frame(
#         True = pred_data$true,
#         Predicted = pred_data$fitted
#       )
#       
#       p <- ggplot(plot_data, aes(x = True, y = Predicted)) +
#         geom_point(alpha = 0.5) +
#         geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
#         labs(title = paste(model_name, "-", scenario_name),
#              x = "True values", y = "Predicted values") +
#         theme_minimal()
#       
#       plots[[model_name]] <- p
#     }
#   }
#   
#   return(plots)
# }

main_simulation <- function() {
  
  all_results <- list()
  
  # Only test 3 specific scenarios
  scenarios_to_test <- c("beta", "arctangent", "sigmoid")
  
  for (scenario in scenarios_to_test) {
    all_results[[scenario]] <- run_simulation_study(scenario, n_iterations)
  }
  
  # Create the width plot with coverage annotations
  width_plot <- create_width_plot(all_results)
  ggsave("width_comparison.png", width_plot, width = 8, height = 5, dpi = 300, bg = "white")
  
  # perf_plots <- create_performance_plots(all_results)
  # 
  # ggsave("coverage_comparison.png", perf_plots$coverage, width = 12, height = 8)
  # ggsave("rmse_comparison.png", perf_plots$rmse, width = 12, height = 8)
  # ggsave("time_comparison.png", perf_plots$time, width = 12, height = 8)
  # 
  # for (scenario in scenarios_to_test) {
  #   pred_plots <- create_prediction_plots(all_results, scenario)
  #   
  #   if (length(pred_plots) > 0) {
  #     combined_plot <- do.call(grid.arrange, c(pred_plots, ncol = 2))
  #     ggsave(paste0("predictions_", scenario, ".png"), combined_plot, width = 12, height = 8)
  #   }
  # }
  
  cat("\n=== MEAN RMSE TABLE ===\n")
  combined_metrics <- do.call(rbind, lapply(all_results, function(x) x$metrics))
  combined_metrics <- combined_metrics[!is.na(combined_metrics$RMSE), ]
  
  # Create RMSE table for the 4 models over 3 transformations
  rmse_table <- combined_metrics %>%
    group_by(Model, Scenario) %>%
    summarise(Mean_RMSE = mean(RMSE, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Scenario, values_from = Mean_RMSE) %>%
    arrange(match(Model, c("sebart", "bart", "bbart_bc", "sblm")))
  
  # Round values for cleaner display
  rmse_table <- rmse_table %>%
    mutate(across(where(is.numeric), ~round(., 3)))
  
  print(rmse_table)
  write.csv(rmse_table, "mean_rmse_table.csv", row.names = FALSE)
  
  # summary_stats <- combined_metrics %>%
  #   group_by(Model, Scenario) %>%
  #   summarise(
  #     Mean_RMSE = mean(RMSE, na.rm = TRUE),
  #     Median_Coverage = median(Coverage, na.rm = TRUE),
  #     Mean_Width = mean(Mean_Width, na.rm = TRUE),
  #     Mean_Time = mean(Time_sec, na.rm = TRUE),
  #     .groups = 'drop'
  #   )
  # 
  # print(summary_stats)
  # 
  # saveRDS(all_results, "simulation_results.rds")
  # write.csv(combined_metrics, "detailed_metrics.csv", row.names = FALSE)
  # write.csv(summary_stats, "summary_statistics.csv", row.names = FALSE)
  
  return(all_results)
}