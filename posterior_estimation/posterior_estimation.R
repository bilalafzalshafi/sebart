library(dbarts)

# ============================================================================
# QUICK BART POSTERIOR ESTIMATION
# ============================================================================
# Simplified implementation of posterior estimation using a short BART run.
# This provides a reasonable approximation of p̂(θ|data) for transformation
# learning in semiparametric BART.
# ============================================================================

#' Estimate BART posterior using quick BART run
#' 
#' Runs a short BART chain to get a reasonable approximation of the posterior.
#' This captures data information while keeping computational cost manageable.
#' 
#' @param X Covariate matrix
#' @param z Transformed response vector  
#' @param n_samples Number of posterior samples to draw
#' @param n_burn Number of burn-in samples
#' @return List with posterior approximation
estimate_bart_posterior_quick <- function(X, z, n_samples = 50, n_burn = 25) {
  
  # Suppress BART output for cleaner integration
  suppressMessages({
    bart_fit <- bart(x.train = X, 
                     y.train = z,
                     ndpost = n_samples,
                     nskip = n_burn, 
                     verbose = FALSE,
                     printevery = 9999)
  })
  
  # Extract posterior mean predictions
  mean_predictions <- bart_fit$yhat.train.mean
  
  # Extract posterior mean of sigma
  sigma_mean <- mean(bart_fit$sigma)
  
  # For uncertainty, use posterior standard deviation of predictions
  pred_sd <- apply(bart_fit$yhat.train, 2, sd)
  
  list(
    mean_predictions = mean_predictions,
    prediction_sd = pred_sd,
    sigma = sigma_mean,
    posterior_samples = bart_fit$yhat.train,  # Full posterior samples
    sigma_samples = bart_fit$sigma,
    n_samples = n_samples
  )
}

cat("Quick BART posterior estimation loaded\n")