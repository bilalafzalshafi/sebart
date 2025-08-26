# Load required libraries and functions
library(dbarts)
library(SeBR)
source("sebart.R")
source("simulation_helpers.R")

# Create poster visualization for transformation effect
create_transformation_poster_plot <- function(scenario = "sigmoid", 
                                            n_train = 200, 
                                            n_test = 500,
                                            p = 10,
                                            seed = 123) {
  
  set.seed(seed)
  
  sim_data <- simulate_sebart_data(n_train = 200, 
                                  n_test = 500, 
                                  p = p,
                                  scenario = scenario, 
                                  seed = seed)
  
  # Extract data
  X_train <- sim_data$X_train
  y_train <- sim_data$y_train
  X_test <- sim_data$X_test
  y_test <- sim_data$y_test
  z_train <- sim_data$z_train
  z_test <- sim_data$z_test

  sebart_fit <- sebart(y = y_train, 
                     X = X_train, 
                     X_test = X_test,
                     pilot_method = "bart",
                     ntree = 200, 
                     nsave = 1000, 
                     nburn = 1000, 
                     verbose = FALSE)
  
  # Apply learned transformation to get g(y)
  g_posterior_mean <- colMeans(sebart_fit$post_g)
  y_unique_sorted <- sort(unique(y_train))
  
  cat("Range of g_posterior_mean:", round(range(g_posterior_mean), 4), "\n")
  
  # Interpolate
  g_function <- approxfun(y_unique_sorted, g_posterior_mean, rule = 2)
  g_y_train <- g_function(y_train)
  g_y_test <- g_function(y_test)
  
  # Create visualization using first predictor
  x1_train <- X_train[,1]
  x1_test <- X_test[,1]
  
  par(mfrow = c(1, 2), 
      mar = c(4.5, 4.5, 3, 2), 
      cex.main = 1.4,
      cex.lab = 1.2,
      cex.axis = 1.1)
  
  plot(x1_train, y_train, 
       pch = 16, 
       col = rgb(0, 0, 0, 0.5),
       xlab = expression(X[1]), 
       ylab = "Y",
       main = "Before Transformation",
       xlim = c(0, 1),
       ylim = range(y_train))
  
  plot(x1_train, g_y_train, 
       pch = 16, 
       col = rgb(0, 0, 0, 0.5),
       xlab = expression(X[1]), 
       ylab = "g(Y)",
       main = "After SeBART Transformation",
       xlim = c(0, 1))
  

  
  par(mfrow = c(1, 1))
  
  invisible(list(
    sebart_fit = sebart_fit,
    g_function = g_function,
    correlation_original = cor(x1_train, y_train),
    correlation_transformed = cor(x1_train, g_y_train),
    r2_transformed = r2
  ))
}

if(interactive()) {
  results <- create_single_poster_figure("beta")
  
  # Or create both examples
  # all_results <- create_poster_examples()
  
  # Optionally visualize the learned transformation
  # visualize_learned_transformation(results$sebart_fit, "sigmoid")
}