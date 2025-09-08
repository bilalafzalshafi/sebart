# sebart diagnostic functions

library(SeBR)
source("sebart.R")
source("simulation_helpers.R")

#' Point estimate plot for sebart
#' 
#' @param sebart_fit sebart fit object or fitted values
#' @param y_test True test values
#' @param main_title Plot title
#' @return RMSE value
plot_sebart_point_estimates <- function(sebart_fit, y_test, main_title = "Point estimates: testing data") {
  
  if (is.list(sebart_fit) && "fitted.values" %in% names(sebart_fit)) {
    y_hat <- sebart_fit$fitted.values
  } else if (is.list(sebart_fit) && "post_ypred" %in% names(sebart_fit)) {
    y_hat <- apply(sebart_fit$post_ypred, 2, mean)
  } else if (is.numeric(sebart_fit)) {
    y_hat <- sebart_fit
  } else {
    stop("sebart_fit must contain fitted.values, post_ypred, or be numeric vector")
  }
  
  rmse_val <- sqrt(mean((y_hat - y_test)^2))
  
  plot_range <- range(c(y_test, y_hat))
  
  plot(y_test, y_hat, 
       xlab = 'y_test', ylab = 'y_hat', 
       main = main_title,
       pch = 16, cex = 0.6,
       xlim = plot_range, ylim = plot_range,
       cex.main = 1.2, cex.lab = 1.1)
  
  abline(0, 1, col = 'red', lwd = 2, lty = 2)
  
  text(min(plot_range) + 0.1 * diff(plot_range), 
       max(plot_range) - 0.1 * diff(plot_range),
       paste("RMSE:", round(rmse_val, 3)), 
       cex = 0.9, col = "blue")
  
  grid(col = "lightgray", lty = "dotted")
  
  return(rmse_val)
}

#' Prediction interval plot for sebart
#' 
#' @param post_ypred Posterior predictive samples (nsave x n_test matrix)
#' @param y_test True test values
#' @param alpha_level Alpha level for prediction intervals (default 0.10 for 90%)
#' @return Coverage rate
plot_pptest <- function(post_ypred, y_test, alpha_level = 0.10) {
  
  lower_quantile <- alpha_level / 2
  upper_quantile <- 1 - alpha_level / 2
  
  y_hat <- apply(post_ypred, 2, mean)  # Posterior predictive mean
  y_lower <- apply(post_ypred, 2, quantile, probs = lower_quantile)
  y_upper <- apply(post_ypred, 2, quantile, probs = upper_quantile)
  
  coverage <- mean(y_test >= y_lower & y_test <= y_upper)
  
  all_values <- c(y_test, y_hat, y_lower, y_upper)
  plot_range <- range(all_values, na.rm = TRUE)
  
  plot(plot_range, plot_range, type = 'n',
       xlab = 'y_test', ylab = 'y_hat', 
       main = "Prediction intervals: testing data",
       cex.main = 1.2, cex.lab = 1.1)
  
  abline(0, 1, col = 'gray', lwd = 1)
  
  for (i in 1:length(y_test)) {
    lines(c(y_test[i], y_test[i]), c(y_lower[i], y_upper[i]), 
          col = 'gray', lwd = 1)
  }
  
  points(y_test, y_hat, pch = 2, cex = 0.6, col = 'black')
  
  grid(col = "lightgray", lty = "dotted")
  
  return(coverage)
}

#' Posterior predictive ECDF plot for sebart
#' 
#' @param sebart_fit Output from sebart() function
#' @param y_test True test data values  
#' @param main_title Title for the plot
#' @param n_samples Number of posterior predictive samples to plot (for speed)
#' @return Creates a plot showing posterior predictive adequacy
plot_sebart_ppd <- function(sebart_fit, y_test, main_title = "Posterior predictive ECDF: testing data", n_samples = 100) {
  
  post_ypred <- sebart_fit$post_ypred
  
  if (is.null(post_ypred)) {
    stop("sebart_fit must contain post_ypred component")
  }
  
  n_mcmc_samples <- nrow(post_ypred)
  if (n_mcmc_samples > n_samples) {
    sample_indices <- sample(1:n_mcmc_samples, n_samples)
    post_ypred <- post_ypred[sample_indices, ]
  }
  
  all_values <- c(y_test, as.vector(post_ypred))
  y_range <- range(all_values, na.rm = TRUE)
  y_grid <- seq(y_range[1], y_range[2], length.out = 200)
  
  plot(y_grid, y_grid, type = 'n', ylim = c(0, 1),
       xlab = 'y', ylab = 'F_y', main = main_title,
       cex.main = 1.2, cex.lab = 1.1)
  
  # Add ECDF for each posterior predictive sample (gray lines)
  for (i in 1:nrow(post_ypred)) {
    pred_sample <- post_ypred[i, ]
    pred_sample <- pred_sample[!is.na(pred_sample)]  # Remove any NAs
    
    if (length(pred_sample) > 0) {
      ecdf_pred <- ecdf(pred_sample)
      lines(y_grid, ecdf_pred(y_grid), col = 'gray', type = 's', lwd = 0.5)
    }
  }
  
  # Add ECDF for observed test data (black line, thicker)
  ecdf_obs <- ecdf(y_test)
  lines(y_grid, ecdf_obs(y_grid), col = 'black', type = 's', lwd = 3)
  
  grid(col = "lightgray", lty = "dotted")
}

#' Posterior transformation plot for sebart
#' 
#' @param sebart_fit Output from sebart() function
#' @param y_values Values at which transformation is evaluated (default: unique y values)
#' @param true_g_values True transformation values (optional, for comparison)
#' @param main_title Title for the plot
#' @param n_draws Number of posterior draws to show (for visual clarity)
#' @param standardize If TRUE, standardize both posterior and true transformations for comparison
#' @return Creates a plot showing posterior draws of g
plot_sebart_transformation <- function(sebart_fit, y_values = NULL, true_g_values = NULL, 
                                    main_title = "Posterior draws of transformation", 
                                    n_draws = 50, standardize = FALSE) {
  
  post_g <- sebart_fit$post_g
  
  if (is.null(post_g)) {
    stop("sebart_fit must contain post_g component")
  }
  
  if (is.null(y_values)) {
    y_values <- sort(unique(sebart_fit$y))
  }
  
  if (ncol(post_g) != length(y_values)) {
    warning("Dimension mismatch between post_g and y_values. Using first ncol(post_g) y values.")
    y_values <- y_values[1:ncol(post_g)]
  }
  
  g_mean <- colMeans(post_g)

  cat("Range of g_mean:", round(range(g_mean), 4), "\n")
  
  if (standardize && !is.null(true_g_values)) {
    # Ensure true_g_values has same length as g_mean
    if (length(true_g_values) != length(g_mean)) {
      warning("Length mismatch between true_g_values and posterior samples. Truncating to shorter length.")
      min_len <- min(length(true_g_values), length(g_mean))
      true_g_values <- true_g_values[1:min_len]
      g_mean <- g_mean[1:min_len]
      post_g <- post_g[, 1:min_len]
      y_values <- y_values[1:min_len]
    }
    
    # Standardize posterior mean: center at zero, scale to unit variance
    g_mean_std <- (g_mean - mean(g_mean)) / sd(g_mean)
    
    # Standardize true transformation: center at zero, scale to unit variance
    true_g_std <- (true_g_values - mean(true_g_values)) / sd(true_g_values)
    
    # Standardize all posterior draws using the same centering/scaling as posterior mean
    post_g_std <- t(apply(post_g, 1, function(g_draw) {
      (g_draw - mean(g_mean)) / sd(g_mean)
    }))
    
    g_mean_plot <- g_mean_std
    true_g_plot <- true_g_std
    post_g_plot <- post_g_std
    
  } else {
    g_mean_plot <- g_mean
    true_g_plot <- true_g_values
    post_g_plot <- post_g
  }
  
  # Subset
  n_total_draws <- nrow(post_g_plot)
  if (n_total_draws > n_draws) {
    draw_indices <- sample(1:n_total_draws, n_draws)
    post_g_subset <- post_g_plot[draw_indices, ]
  } else {
    post_g_subset <- post_g_plot
  }
  
  all_g_values <- c(as.vector(post_g_subset), g_mean_plot)
  if (!is.null(true_g_plot)) {
    all_g_values <- c(all_g_values, true_g_plot)
  }
  
  y_range <- range(y_values)
  g_range <- range(all_g_values, na.rm = TRUE)
  
  plot(y_range, g_range, type = 'n',
       xlab = 'y', ylab = ifelse(standardize && !is.null(true_g_values), 'g(y) [standardized]', 'g(y)'), 
       main = main_title,
       cex.main = 1.2, cex.lab = 1.1)
  
  for (i in 1:nrow(post_g_subset)) {
    lines(y_values, post_g_subset[i, ], col = 'gray', lwd = 0.5)
  }
  
  lines(y_values, g_mean_plot, col = 'black', lwd = 3)
  
  if (!is.null(true_g_plot)) {
    points(y_values, true_g_plot, pch = 2, cex = 0.8, col = 'black')
    
    legend_text <- if (standardize) {
      c('Posterior mean (std)', 'Truth (std)')
    } else {
      c('Posterior mean', 'Truth')
    }
    
    legend('bottomright', legend_text, 
           lty = c(1, NA), pch = c(NA, 2), lwd = c(3, NA), 
           col = c('black', 'black'), cex = 0.8)
  }
  
  grid(col = "lightgray", lty = "dotted")
  
  if (!is.null(true_g_plot)) {
    correlation <- cor(g_mean_plot, true_g_plot)
    if (standardize) {
      cat("Correlation between standardized posterior mean and true transformation:", round(correlation, 3), "\n")
    } else {
      cat("Correlation between posterior mean and true transformation:", round(correlation, 3), "\n")
    }
  }
}

# #' Plot raw (non-standardized) transformation comparison for sebart
# #' 
# #' @param sebart_fit Output from sebart() function
# #' @param y_values Values at which transformation is evaluated (default: unique y values)
# #' @param true_g_values True transformation values (optional, for comparison)
# #' @param main_title Title for the plot
# #' @param n_draws Number of posterior draws to show (for visual clarity)
# #' @return Creates a plot showing posterior draws of g without standardization
# plot_sebart_transformation_raw <- function(sebart_fit, y_values = NULL, true_g_values = NULL, 
#                                         main_title = "Posterior draws of transformation (raw scale)", 
#                                         n_draws = 50) {
  
#   post_g <- sebart_fit$post_g
  
#   if (is.null(post_g)) {
#     stop("sebart_fit must contain post_g component")
#   }
  
#   if (is.null(y_values)) {
#     y_values <- sort(unique(sebart_fit$y))
#   }
  
#   if (ncol(post_g) != length(y_values)) {
#     warning("Dimension mismatch between post_g and y_values. Using first ncol(post_g) y values.")
#     y_values <- y_values[1:ncol(post_g)]
#   }
  
#   g_mean <- colMeans(post_g)
  
#   # NO STANDARDIZATION - use raw values
#   if (!is.null(true_g_values)) {
#     # Ensure true_g_values has same length as g_mean
#     if (length(true_g_values) != length(g_mean)) {
#       warning("Length mismatch between true_g_values and posterior samples. Truncating to shorter length.")
#       min_len <- min(length(true_g_values), length(g_mean))
#       true_g_values <- true_g_values[1:min_len]
#       g_mean <- g_mean[1:min_len]
#       post_g <- post_g[, 1:min_len]
#       y_values <- y_values[1:min_len]
#     }
    
#     # Use raw values - no standardization
#     g_mean_plot <- g_mean
#     true_g_plot <- true_g_values
#     post_g_plot <- post_g
    
#   } else {
#     g_mean_plot <- g_mean
#     true_g_plot <- NULL
#     post_g_plot <- post_g
#   }
  
#   # Subset posterior draws
#   n_total_draws <- nrow(post_g_plot)
#   if (n_total_draws > n_draws) {
#     draw_indices <- sample(1:n_total_draws, n_draws)
#     post_g_subset <- post_g_plot[draw_indices, ]
#   } else {
#     post_g_subset <- post_g_plot
#   }
  
#   all_g_values <- c(as.vector(post_g_subset), g_mean_plot)
#   if (!is.null(true_g_plot)) {
#     all_g_values <- c(all_g_values, true_g_plot)
#   }
  
#   y_range <- range(y_values)
#   g_range <- range(all_g_values, na.rm = TRUE)
  
#   plot(y_range, g_range, type = 'n',
#        xlab = 'y', ylab = 'g(y) [raw scale]', 
#        main = main_title,
#        cex.main = 1.2, cex.lab = 1.1)
  
#   # Plot posterior draws in gray
#   for (i in 1:nrow(post_g_subset)) {
#     lines(y_values, post_g_subset[i, ], col = 'gray', lwd = 0.5)
#   }
  
#   # Plot posterior mean in black
#   lines(y_values, g_mean_plot, col = 'black', lwd = 3)
  
#   # Plot true transformation if available
#   if (!is.null(true_g_plot)) {
#     points(y_values, true_g_plot, pch = 2, cex = 0.8, col = 'black')
    
#     legend('bottomright', legend = c('Posterior mean (raw)', 'Truth (raw)'), 
#            lty = c(1, NA), pch = c(NA, 2), lwd = c(3, NA), 
#            col = c('black', 'black'), cex = 0.8)
#   }
  
#   grid(col = "lightgray", lty = "dotted")
  
#   # Compute and report correlation if both are available
#   if (!is.null(true_g_plot)) {
#     correlation <- cor(g_mean_plot, true_g_plot)
#     cat("Correlation between raw posterior mean and true transformation:", round(correlation, 3), "\n")
    
#     # Also report scale comparison metrics
#     scale_ratio <- sd(g_mean_plot) / sd(true_g_plot)
#     location_diff <- mean(g_mean_plot) - mean(true_g_plot)
#     cat("Scale ratio (learned/true):", round(scale_ratio, 3), "\n")
#     cat("Location difference (learned - true):", round(location_diff, 3), "\n")
#   }
# }

#' Create all plots
#' 
#' @param sebart_fit Output from sebart() function
#' @param y_test True test data values
#' @param true_g_values True transformation values (optional)
#' @param alpha_level Alpha level for prediction intervals
#' @param layout Plot layout (default c(2,2) for 4 plots)
create_all_sebart_diagnostics <- function(sebart_fit, y_test, true_g_values = NULL, 
                                       alpha_level = 0.10, layout = c(2, 2)) {
  
  # Adjust layout if variable selection is enabled
  if (sebart_fit$var_select) {
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))
  } else {
    par(mfrow = layout, mar = c(4, 4, 3, 2))
  }
  
  # Plot 1: Point estimates
  rmse <- plot_sebart_point_estimates(sebart_fit, y_test)
  
  # Plot 2: Prediction intervals
  coverage <- plot_pptest(sebart_fit$post_ypred, y_test, alpha_level)
  
  # Plot 3: Posterior predictive ECDF
  plot_sebart_ppd(sebart_fit, y_test)
  
  # Plot 4: Transformation
  y_values <- sort(unique(sebart_fit$y))
  plot_sebart_transformation(sebart_fit, y_values, true_g_values, standardize = FALSE)

  # Plot 5: Variable selection (if enabled)
  if (sebart_fit$var_select) {
    # Simple bar plot of variable importance
    p <- length(sebart_fit$var_importance_scores)
    barplot(sebart_fit$var_importance_scores, 
            names.arg = 1:p,
            xlab = 'Variable Index', ylab = 'Importance Score',
            main = 'Variable Importance Scores',
            col = ifelse(1:p %in% sebart_fit$selected_vars, 'lightblue', 'lightgray'),
            cex.main = 1.2, cex.lab = 1.1)
    grid(col = "lightgray", lty = "dotted")
  }

  # # Plot 6: Transformation (raw scale) - commented out
  # plot_sebart_transformation_raw(sebart_fit, y_values, true_g_values)
  
  par(mfrow = c(1, 1))
  
  cat("RMSE:", round(rmse, 3), "\n")
  cat("Coverage:", round(coverage, 3), "\n")
  
  # Add variable selection summary if enabled
  if (sebart_fit$var_select) {
    print_var_selection_summary(sebart_fit)
  }
  
  return(list(rmse = rmse, coverage = coverage))
}

#' Variable selection sensitivity to number of trees
#' 
#' @param sebart_fit Output from sebart() function with variable selection enabled
#' @param ntree_range Vector of tree numbers to test (default: c(50, 100, 200))
#' @param main_title Title for the plot
#' @return Creates a plot showing how variable importance changes with number of trees
plot_var_selection_trees <- function(sebart_fit, ntree_range = c(50, 100, 200), 
                                     main_title = "Variable Importance vs Number of Trees") {
  
  if (!sebart_fit$var_select) {
    stop("sebart_fit must have variable selection enabled")
  }
  
  # Get original data
  y <- sebart_fit$y
  X <- sebart_fit$X
  X_test <- sebart_fit$X_test
  
  # Fit models with different numbers of trees
  importance_results <- list()
  
  for (ntree in ntree_range) {
    fit <- sebart(y = y, X = X, X_test = X_test,
                  ntree = ntree, nsave = 200, nburn = 200,
                  var_select = TRUE, var_select_threshold = 0.01,
                  verbose = FALSE)
    
    importance_results[[as.character(ntree)]] <- fit$var_importance_scores
  }
  
  # Create plot
  p <- length(importance_results[[1]])
  ntree_labels <- names(importance_results)
  
  plot(1, 1, type = 'n', xlim = c(1, p), ylim = c(0, max(unlist(importance_results))),
       xlab = 'Variable Index', ylab = 'Importance Score',
       main = main_title, cex.main = 1.2, cex.lab = 1.1)
  
  colors <- rainbow(length(ntree_range))
  
  for (i in 1:length(ntree_range)) {
    lines(1:p, importance_results[[i]], col = colors[i], lwd = 2, type = 'b', pch = 16)
  }
  
  legend('topright', legend = paste(ntree_labels, 'trees'), 
         col = colors, lwd = 2, pch = 16, cex = 0.8)
  
  grid(col = "lightgray", lty = "dotted")
  
  return(importance_results)
}

#' Display variable selection results
#' 
#' @param sebart_fit Output from sebart() function with variable selection enabled
#' @return Prints variable selection summary
print_var_selection_summary <- function(sebart_fit) {
  
  if (!sebart_fit$var_select) {
    cat("Variable selection was not enabled for this model.\n")
    return(invisible(NULL))
  }
  
  cat("\n=== Variable Selection Summary ===\n")
  cat("Threshold:", sebart_fit$var_select_threshold, "\n")
  cat("Selected variables:", length(sebart_fit$selected_vars), "of", length(sebart_fit$var_importance_scores), "\n")
  cat("Variable indices:", paste(sebart_fit$selected_vars, collapse = ", "), "\n")
  
  cat("\nImportance scores:\n")
  for (i in 1:length(sebart_fit$var_importance_scores)) {
    status <- if (i %in% sebart_fit$selected_vars) " [SELECTED]" else ""
    cat(sprintf("  Variable %d: %.4f%s\n", i, sebart_fit$var_importance_scores[i], status))
  }
  
  cat("\n")
}

#' Demo function with sebart diagnostic workflow
#' 
#' @param scenario Transformation scenario to test
#' @param n_train Number of training observations
#' @param n_test Number of test observations
#' @param p Number of predictors
demo_all_sebart_diagnostics <- function(scenario = "box_cox", n_train = 200, n_test = 500, p = 10, var_select = TRUE) {
  
  cat("Generating data for scenario:", scenario, "\n")
  
  set.seed(123)
  
  sim_data <- simulate_sebart_data(n_train = n_train, n_test = n_test, p = p, 
                                scenario = scenario, seed = 123)
  
  if (!scenario %in% get_available_scenarios()) {
    stop("Unknown scenario. Available: ", paste(get_available_scenarios(), collapse = ", "))
  }
  
  X_train <- sim_data$X_train
  X_test <- sim_data$X_test
  y_train <- sim_data$y_train
  y_test <- sim_data$y_test
  f_true_train <- sim_data$f_true_train
  f_true_test <- sim_data$f_true_test
  z_train <- sim_data$z_train
  z_test <- sim_data$z_test
  
  cat("Fitting sebart model...\n")

  ntree_demo <- if (var_select) 10 else 200
  nsave_demo <- if (var_select) 100 else 500
  nburn_demo <- if (var_select) 100 else 500
  
  sebart_fit <- sebart(y = y_train, X = X_train, X_test = X_test, 
                     ntree = ntree_demo, nsave = nsave_demo, nburn = nburn_demo, 
                     var_select = TRUE, var_select_threshold = 0.2,
                     verbose = FALSE)
  
  cat("Creating all diagnostic plots...\n")
  
  true_g_values <- sim_data$g_true
  y_unique <- sim_data$y_unique
  
  metrics <- create_all_sebart_diagnostics(sebart_fit, y_test, true_g_values)
  
  return(list(
    sebart_fit = sebart_fit,
    y_test = y_test,
    metrics = metrics,
    true_g_values = true_g_values
  ))
}

# Example usage
if (interactive()) {
  results <- demo_all_sebart_diagnostics('arctangent')
}