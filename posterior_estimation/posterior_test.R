library(dbarts)

# ============================================================================
# TEST: Quick BART Posterior Estimation
# ============================================================================
# Tests the accuracy of quick BART posterior approximation against a 
# gold standard long BART run.
# ============================================================================

source("posterior_estimation.R")

#' Generate synthetic test data
generate_test_data <- function(n = 100, p = 3, seed = 123) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  z <- X[,1] + 0.5 * X[,2]^2 - 0.3 * X[,3] + 0.2 * X[,1] * X[,2] + rnorm(n, 0, 0.5)
  list(X = X, z = z)
}

#' Create gold standard using long BART run
create_gold_standard <- function(X, z) {
  cat("Creating gold standard...\n")
  suppressMessages({
    gold_bart <- bart(x.train = X, y.train = z, ndpost = 1000, nskip = 500,
                      verbose = FALSE, printevery = 9999)
  })
  gold_bart$yhat.train.mean
}

#' Test quick BART accuracy vs gold standard
test_quick_bart_accuracy <- function() {
  cat("=== Testing Quick BART Accuracy ===\n")
  
  # Generate test data
  data <- generate_test_data()
  X <- data$X
  z <- data$z
  
  # Get gold standard
  gold_predictions <- create_gold_standard(X, z)
  
  # Get quick BART approximation
  cat("Running quick BART...\n")
  quick_result <- estimate_bart_posterior_quick(X, z)
  
  # Compute accuracy
  rmse_vs_gold <- sqrt(mean((quick_result$mean_predictions - gold_predictions)^2))
  
  cat("RMSE vs gold standard:", round(rmse_vs_gold, 4), "\n")
  
  # Test passes if RMSE is reasonable (< 0.5)
  success <- rmse_vs_gold < 0.5
  cat("Test result:", if(success) "PASS" else "FAIL", "\n")
  
  return(list(
    success = success,
    rmse_vs_gold = rmse_vs_gold
  ))
}

#' Run test
run_posterior_test <- function() {
  cat("TESTING QUICK BART POSTERIOR ESTIMATION\n")
  cat("======================================\n")
  
  result <- test_quick_bart_accuracy()
  
  cat("\n=== SUMMARY ===\n")
  cat("Quick BART accuracy:", if(result$success) "PASS" else "FAIL", "\n")
  
  if (result$success) {
    cat("Quick BART approximation is adequate for transformation learning.\n")
  } else {
    cat("Quick BART approximation may be insufficient.\n")
  }
  
  return(result)
}

if (interactive()) {
  test_result <- run_posterior_test()
} else {
  cat("Test loaded. Run: test_result <- run_posterior_test()\n")
}