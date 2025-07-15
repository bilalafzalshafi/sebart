# Test whether BART's setResponse() triggers internal re-standardization
# that could interfere with transformation learning in sbart

library(dbarts)

# Test 1: Check if setResponse re-standardizes data
test_bart_standardization <- function() {
  set.seed(123)
  
  # Generate test data with different scales
  n <- 100
  X <- matrix(rnorm(n * 2), n, 2)
  
  # Original response (large scale)
  y1 <- 1000 + 500 * rnorm(n)  
  
  # Create BART sampler with original response
  sampler1 <- dbarts(y1 ~ ., data = data.frame(y1 = y1, X))
  
  # Get initial fitted values (these will be standardized internally)
  initial_run <- sampler1$run(numBurnIn = 10, numSamples = 1)
  initial_fits <- initial_run$train[1, ]
  
  # New response (different scale)
  y2 <- 10 + 5 * rnorm(n)
  
  # Update response and check if re-standardization occurs
  sampler1$setResponse(y2)
  updated_run <- sampler1$run(numBurnIn = 0, numSamples = 1)
  updated_fits <- updated_run$train[1, ]
  
  # If BART re-standardizes, the fitted values should be on similar scales
  # If not, they should reflect the different input scales
  
  cat("Original y range:", range(y1), "\n")
  cat("Updated y range:", range(y2), "\n")
  cat("Initial fitted values range:", range(initial_fits), "\n")
  cat("Updated fitted values range:", range(updated_fits), "\n")
  
  # Check if the scale difference is preserved
  y_ratio <- diff(range(y1)) / diff(range(y2))
  fit_ratio <- diff(range(initial_fits)) / diff(range(updated_fits))
  
  cat("Input scale ratio:", y_ratio, "\n")
  cat("Output scale ratio:", fit_ratio, "\n")
  cat("Ratios similar (no re-standardization)?", abs(log(y_ratio/fit_ratio)) < 0.1, "\n\n")
  
  return(abs(log(y_ratio/fit_ratio)) < 0.1)
}

# Test 2: Verify transformation learning happens on original scale
test_transformation_scale <- function() {
  set.seed(456)
  
  # Simulate data similar to sbart setup
  n <- 50
  X <- matrix(rnorm(n * 2), n, 2)
  
  # Create data with known transformation (log transform)
  z_true <- X[,1] + 0.5 * X[,2] + rnorm(n)
  y <- exp(z_true)  # y = exp(z), so g(y) = log(y) should recover z
  
  # Manual transformation learning (mimicking sbart approach)
  y0 <- sort(unique(y))
  
  # Bayesian bootstrap for F_Y (original scale)
  alpha_y <- rgamma(n, 1)
  alpha_y <- alpha_y / sum(alpha_y)
  F_Y_bb <- sapply(y0, function(t) sum(alpha_y[y <= t]))
  
  # Check that BB is working on original y scale, not standardized
  cat("Original y range:", range(y), "\n")
  cat("F_Y values range:", range(F_Y_bb), "\n")
  cat("F_Y is proper CDF (0 to 1)?", min(F_Y_bb) >= 0 && max(F_Y_bb) <= 1, "\n")
  
  # Verify transformation recovers known relationship
  # For log transform: g(y) = log(y), so g should map y to z scale
  z_recovered <- qnorm(F_Y_bb)  # Simple normal quantile approximation
  y_grid <- y0
  
  # Check if recovered transformation is approximately log
  log_y <- log(y_grid)
  correlation <- cor(z_recovered, log_y, use = "complete.obs")
  
  cat("Correlation between recovered g(y) and log(y):", correlation, "\n")
  cat("High correlation suggests correct scale?", correlation > 0.9, "\n\n")
  
  return(correlation > 0.9)
}

# Test 3: Full integration test mimicking sbart workflow
test_sbart_integration <- function() {
  set.seed(789)
  
  n <- 100
  X <- matrix(rnorm(n * 2), n, 2)
  
  # True model: z = X %*% beta + epsilon, y = g_inv(z)
  beta_true <- c(1, -0.5)
  z_true <- X %*% beta_true + rnorm(n, sd = 0.5)
  y <- z_true^2  # Square transformation, so g(y) = sqrt(y)
  
  # Initialize BART sampler (mimicking sbart setup)
  z_init <- qnorm((rank(y) - 0.5) / n)  # Initial quantile transform
  z_scaled <- (z_init - mean(z_init)) / sd(z_init)
  
  X_df <- data.frame(X1 = X[,1], X2 = X[,2])
  bart_data <- cbind(data.frame(z_scaled = z_scaled), X_df)
  
  sampler <- dbarts(z_scaled ~ ., data = bart_data)
  
  # Simulate one iteration of transformation update
  # 1. Learn transformation from original y
  y0 <- sort(unique(y))
  alpha_y <- rgamma(n, 1); alpha_y <- alpha_y / sum(alpha_y)
  F_Y <- sapply(y0, function(t) sum(alpha_y[y <= t]))
  
  # 2. Apply transformation
  g_y <- qnorm(pmax(pmin(F_Y, 0.999), 0.001))  # Avoid boundary issues
  g_interp <- approxfun(y0, g_y, rule = 2)
  z_new <- g_interp(y)
  
  # 3. Update BART with transformed data
  mu_new <- mean(z_new)
  sigma_new <- sd(z_new)
  z_new_scaled <- (z_new - mu_new) / sigma_new
  
  # Key test: setResponse with scaled transformed data
  sampler$setResponse(z_new_scaled)
  result <- sampler$run(numBurnIn = 5, numSamples = 1)
  
  # Check that BART produces reasonable fits
  fits <- result$train[1, ] * sigma_new + mu_new  # Unscale
  
  cat("True z range:", range(z_true), "\n")
  cat("Recovered z range:", range(fits), "\n")
  cat("Correlation with true z:", cor(fits, z_true), "\n")
  
  # Verify scales are reasonable
  reasonable_correlation <- cor(fits, z_true) > 0.5
  reasonable_scale <- abs(log(sd(fits) / sd(z_true))) < 1
  
  cat("Reasonable correlation?", reasonable_correlation, "\n")
  cat("Reasonable scale?", reasonable_scale, "\n\n")
  
  return(reasonable_correlation && reasonable_scale)
}

# Run all tests
cat("=== Testing BART Standardization Behavior ===\n")
test1_result <- test_bart_standardization()

cat("=== Testing Transformation Scale ===\n") 
test2_result <- test_transformation_scale()

cat("=== Testing SBART Integration ===\n")
test3_result <- test_sbart_integration()

cat("=== Summary ===\n")
cat("All tests passed:", all(test1_result, test2_result, test3_result), "\n")
if (!all(test1_result, test2_result, test3_result)) {
  cat("POTENTIAL ISSUE: BART standardization may interfere with transformation learning\n")
} else {
  cat("CONFIRMED: BART standardization does not interfere with sbart transformation learning\n")
}