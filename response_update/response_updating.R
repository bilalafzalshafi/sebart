library(dbarts)

# ============================================================================
# TEST: dbarts Response Data Updating Mechanism
# ============================================================================
# This file tests whether dbarts' setResponse() method works correctly for
# updating response data during ongoing MCMC. This is the core technical
# requirement for our semiparametric BART approach.
#
# Test strategy:
# 1. Fit BART to initial data
# 2. Systematically modify response data 
# 3. Update using setResponse()
# 4. Verify predictions change appropriately
# 5. Confirm sampler state is preserved
# ============================================================================

#' Test basic setResponse functionality
test_basic_response_update <- function() {
  cat("=== Test 1: Basic setResponse functionality ===\n")
  
  set.seed(123)
  n <- 100
  X <- matrix(rnorm(n * 2), n, 2)
  y1 <- X[,1] + X[,2] + rnorm(n, 0, 0.1)  # Initial response
  
  sampler <- dbarts(X, y1, 
                   verbose = FALSE,
                   n.samples = 10,
                   control = dbartsControl(keepTrainingFits = TRUE,
                                          n.burn = 50,
                                          n.trees = 20))
  
  cat("Running initial MCMC...\n")
  result1 <- sampler$run(numBurnIn = 50, numSamples = 10)
  initial_predictions <- colMeans(result1$train)
  initial_sigma <- mean(result1$sigma)
  
  cat("Initial predictions (first 5):", round(initial_predictions[1:5], 3), "\n")
  cat("Initial sigma:", round(initial_sigma, 3), "\n")
  
  # Modify response data (add constant shift)
  shift <- 2.0
  y2 <- y1 + shift
  
  cat("Updating response data (adding shift of", shift, ")...\n")
  
  # Update response using setResponse
  sampler$setResponse(y2, updateState = TRUE)
  
  # Continue MCMC
  cat("Running MCMC after update...\n")
  result2 <- sampler$run(numBurnIn = 0, numSamples = 10)
  updated_predictions <- colMeans(result2$train)
  updated_sigma <- mean(result2$sigma)
  
  cat("Updated predictions (first 5):", round(updated_predictions[1:5], 3), "\n")
  cat("Updated sigma:", round(updated_sigma, 3), "\n")
  
  prediction_shift <- mean(updated_predictions - initial_predictions)
  cat("Average prediction shift:", round(prediction_shift, 3), "\n")
  cat("Expected shift:", shift, "\n")
  
  success <- abs(prediction_shift - shift) < 0.5
  cat("Test result:", if(success) "PASS" else "FAIL", "\n")
  
  return(list(
    success = success,
    prediction_shift = prediction_shift,
    expected_shift = shift,
    initial_sigma = initial_sigma,
    updated_sigma = updated_sigma
  ))
}

#' Test multiple sequential updates
test_multiple_updates <- function() {
  cat("\n=== Test 2: Multiple sequential updates ===\n")
  
  set.seed(456)
  n <- 80
  X <- matrix(rnorm(n * 2), n, 2)
  y_base <- X[,1] + 0.5 * X[,2] + rnorm(n, 0, 0.2)
  
  # Create sampler
  sampler <- dbarts(X, y_base,
                   verbose = FALSE, 
                   n.samples = 5,
                   control = dbartsControl(keepTrainingFits = TRUE,
                                          n.burn = 30,
                                          n.trees = 15))
  
  # Initial burn-in
  sampler$run(numBurnIn = 30, numSamples = 5)
  
  # Track predictions through multiple updates
  transformations <- c(1.0, 1.5, 0.8, 1.2)  # Scale factors
  predictions_history <- matrix(NA, length(transformations), n)
  
  for (i in seq_along(transformations)) {
    scale <- transformations[i]
    y_transformed <- y_base * scale
    
    cat("Update", i, ": scaling response by", scale, "\n")
    
    # Update response
    sampler$setResponse(y_transformed, updateState = TRUE)
    
    # Run a few samples
    result <- sampler$run(numBurnIn = 0, numSamples = 5)
    predictions_history[i, ] <- colMeans(result$train)
    
    cat("  Mean prediction:", round(mean(predictions_history[i, ]), 3), "\n")
  }
  
  # Check if predictions scale appropriately
  base_mean <- mean(predictions_history[1, ])
  scaling_errors <- numeric(length(transformations) - 1)
  
  for (i in 2:length(transformations)) {
    expected_mean <- base_mean * transformations[i] / transformations[1]
    actual_mean <- mean(predictions_history[i, ])
    scaling_errors[i-1] <- abs(actual_mean - expected_mean) / abs(expected_mean)
    
    cat("Transformation", i, ": expected mean =", round(expected_mean, 3),
        ", actual mean =", round(actual_mean, 3),
        ", relative error =", round(scaling_errors[i-1], 3), "\n")
  }
  
  # Test passes if relative errors are reasonable
  max_error <- max(scaling_errors)
  success <- max_error < 0.3  # Allow 30% relative error
  
  cat("Maximum relative error:", round(max_error, 3), "\n")
  cat("Test result:", if(success) "PASS" else "FAIL", "\n")
  
  return(list(
    success = success,
    max_error = max_error,
    predictions_history = predictions_history,
    transformations = transformations
  ))
}

#' Run all tests
run_all_tests <- function() {
  cat("TESTING dbarts RESPONSE DATA UPDATING\n")
  cat("=====================================\n")
  
  results <- list()
  
  # Run tests
  results$basic <- test_basic_response_update()
  results$multiple <- test_multiple_updates()
  
  # Summary
  cat("\n=== SUMMARY ===\n")
  all_passed <- all(sapply(results, function(x) x$success))
  
  for (test_name in names(results)) {
    status <- if(results[[test_name]]$success) "PASS" else "FAIL"
    cat(test_name, ":", status, "\n")
  }
  
  cat("\nOverall result:", if(all_passed) "ALL TESTS PASSED" else "SOME TESTS FAILED", "\n")
  
  if (all_passed) {
    cat("\ndbarts setResponse() mechanism works correctly!\n")
  } else {
    cat("\nSome issues detected. Check test output above.\n")
  }
  
  return(results)
}

# Run the tests
if (interactive()) {
  test_results <- run_all_tests()
} else {
  cat("Test suite loaded. Run: test_results <- run_all_tests()\n")
}