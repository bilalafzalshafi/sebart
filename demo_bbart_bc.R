# End-to-end demonstration of Bayesian BART with Box-Cox transformation (bbart_bc)

# Load required libraries
library(dbarts)
library(SeBR)
source("bbart_bc.R")

# Set seed for reproducibility
set.seed(456)

# Simulate data with Box-Cox transformation
n <- 200
p <- 5
X <- matrix(rnorm(n*p), n, p)

# Create nonlinear function on latent scale
f_true <- sin(2*X[,1]) + X[,2]^2 - X[,3]*X[,4] + rnorm(n, sd=0.1)

# Apply Box-Cox transformation (square-root family)
true_lambda <- 0.5
y <- g_inv_bc(f_true, lambda = true_lambda)

# Create test data
n_test <- 50
X_test <- matrix(rnorm(n_test*p), n_test, p)
f_test_true <- sin(2*X_test[,1]) + X_test[,2]^2 - X_test[,3]*X_test[,4]
y_test_true <- g_inv_bc(f_test_true + rnorm(n_test, sd=0.1), lambda = true_lambda)

cat("Data simulation complete.\n")
cat("Training data: n =", n, ", p =", p, "\n")
cat("Test data: n_test =", n_test, "\n")
cat("True lambda:", true_lambda, "\n")
cat("Response range: [", round(min(y), 2), ",", round(max(y), 2), "]\n")
cat("Latent scale range: [", round(min(f_true), 2), ",", round(max(f_true), 2), "]\n\n")

# Fit model with lambda estimation
cat("Fitting bbart_bc model with lambda estimation...\n")
fit_estimated <- bbart_bc(y = y, X = X, X_test = X_test,
                         ntree = 100,
                         lambda = NULL,  # Estimate lambda
                         sample_lambda = TRUE,
                         nsave = 500,
                         nburn = 200,
                         nskip = 0,
                         verbose = TRUE)

cat("\nModel fitting complete!\n")
print(names(fit_estimated))

# Summary statistics
cat("\nModel summary:\n")
cat("True lambda:", true_lambda, "\n")
cat("Estimated lambda mean:", round(mean(fit_estimated$post_lambda), 3), "\n")
cat("Estimated lambda 95% CI: [", 
    round(quantile(fit_estimated$post_lambda, 0.025), 3), ",",
    round(quantile(fit_estimated$post_lambda, 0.975), 3), "]\n")
cat("Posterior predictive mean range: [", 
    round(min(fit_estimated$fitted.values), 2), ",", 
    round(max(fit_estimated$fitted.values), 2), "]\n")

# Fit model with fixed lambda for comparison
cat("\nFitting model with fixed lambda...\n")
fit_fixed <- bbart_bc(y = y, X = X, X_test = X_test,
                     ntree = 100,
                     lambda = true_lambda,
                     sample_lambda = FALSE,
                     nsave = 300,
                     nburn = 100,
                     verbose = FALSE)

# Plot 1: Lambda posterior distribution
cat("\nCreating diagnostic plots...\n")
png("bbart_lambda_posterior.png", width=800, height=600)
hist(fit_estimated$post_lambda, breaks=30, col='lightblue',
     main='Posterior Distribution of Box-Cox Parameter Lambda',
     xlab='Lambda', ylab='Frequency')
abline(v=true_lambda, col='red', lwd=3, lty=2)
abline(v=mean(fit_estimated$post_lambda), col='blue', lwd=2)
legend('topright', c('True lambda', 'Posterior mean'), 
       col=c('red', 'blue'), lwd=c(3,2), lty=c(2,1))
dev.off()

# Plot 2: Lambda trace plot
png("bbart_lambda_trace.png", width=800, height=600)
plot(fit_estimated$post_lambda, type='l', 
     xlab='MCMC Iteration', ylab='Lambda',
     main='Trace Plot of Box-Cox Parameter Lambda')
abline(h=true_lambda, col='red', lwd=2, lty=2)
abline(h=mean(fit_estimated$post_lambda), col='blue', lwd=2)
dev.off()

# Plot 3: Prediction comparison (estimated vs fixed lambda)
pi_y_est <- t(apply(fit_estimated$post_ypred, 2, quantile, c(0.05, 0.95)))
pi_y_fix <- t(apply(fit_fixed$post_ypred, 2, quantile, c(0.05, 0.95)))

png("bbart_prediction_comparison.png", width=1000, height=600)
par(mfrow=c(1,2))

# Estimated lambda predictions
plot(1:nrow(pi_y_est), pi_y_est[,1], type='n', ylim=range(pi_y_est),
     xlab='Test observation', ylab='y', 
     main='Estimated Lambda: 90% Prediction Intervals')
segments(1:nrow(pi_y_est), pi_y_est[,1], 1:nrow(pi_y_est), pi_y_est[,2], 
         col='lightblue', lwd=2)
points(1:nrow(pi_y_est), fit_estimated$fitted.values, pch=16, col='darkblue')

# Fixed lambda predictions  
plot(1:nrow(pi_y_fix), pi_y_fix[,1], type='n', ylim=range(pi_y_fix),
     xlab='Test observation', ylab='y', 
     main='Fixed Lambda: 90% Prediction Intervals')
segments(1:nrow(pi_y_fix), pi_y_fix[,1], 1:nrow(pi_y_fix), pi_y_fix[,2], 
         col='lightcoral', lwd=2)
points(1:nrow(pi_y_fix), fit_fixed$fitted.values, pch=16, col='darkred')

par(mfrow=c(1,1))
dev.off()

# Plot 4: Box-Cox transformation
y0 <- sort(unique(y))
png("bbart_transformation.png", width=800, height=600)
plot(y0, fit_estimated$post_g[1,], type='n', ylim=range(fit_estimated$post_g),
     xlab='y', ylab='g(y)', main="Box-Cox Transformation: Posterior Draws")

# Plot several transformation draws
for(s in sample(nrow(fit_estimated$post_g), 30)) {
  lines(y0, fit_estimated$post_g[s,], col='gray')
}

# Plot posterior mean transformation
lines(y0, colMeans(fit_estimated$post_g), lwd=3, col='blue')

# Plot true transformation
true_g <- g_bc(y0, lambda = true_lambda)
lines(y0, true_g, lwd=3, col='red', lty=2)

legend('bottomright', c('Posterior draws', 'Posterior mean', 'True g(y)'), 
       col=c('gray', 'blue', 'red'), lwd=c(1,3,3), lty=c(1,1,2))
dev.off()

# Plot 5: Prediction accuracy comparison
png("bbart_prediction_accuracy.png", width=800, height=600)
plot(fit_estimated$fitted.values, fit_fixed$fitted.values,
     xlab='Estimated Lambda Predictions', ylab='Fixed Lambda Predictions',
     main='Prediction Comparison: Estimated vs Fixed Lambda',
     pch=16, col='darkblue')
abline(0, 1, col='red', lwd=2)
correlation <- cor(fit_estimated$fitted.values, fit_fixed$fitted.values)
text(min(fit_estimated$fitted.values), max(fit_fixed$fitted.values),
     paste("Correlation:", round(correlation, 3)), pos=4)
dev.off()

# Calculate interval widths and coverage (if we had true test values)
interval_widths_est <- pi_y_est[,2] - pi_y_est[,1]
interval_widths_fix <- pi_y_fix[,2] - pi_y_fix[,1]

cat("\nPrediction interval width comparison:\n")
cat("Estimated lambda - Mean width:", round(mean(interval_widths_est), 3), "\n")
cat("Fixed lambda - Mean width:", round(mean(interval_widths_fix), 3), "\n")

# Lambda estimation accuracy
lambda_bias <- mean(fit_estimated$post_lambda) - true_lambda
lambda_mse <- mean((fit_estimated$post_lambda - true_lambda)^2)

cat("\nLambda estimation performance:\n")
cat("Bias:", round(lambda_bias, 4), "\n")
cat("MSE:", round(lambda_mse, 4), "\n")
cat("95% CI contains true value:", 
    quantile(fit_estimated$post_lambda, 0.025) <= true_lambda & 
    true_lambda <= quantile(fit_estimated$post_lambda, 0.975), "\n")

# Demonstrate different fixed lambda values
cat("\nTesting different fixed lambda values...\n")
lambda_values <- c(0.25, 0.5, 1.0)
lambda_results <- list()

for(i in seq_along(lambda_values)) {
  cat("Fitting with lambda =", lambda_values[i], "...\n")
  fit_lambda <- bbart_bc(y = y[1:50], X = X[1:50,], X_test = X_test[1:10,],
                        lambda = lambda_values[i], sample_lambda = FALSE,
                        nsave = 100, nburn = 50, verbose = FALSE)
  lambda_results[[i]] <- fit_lambda
}

# Summary of results
cat("\nFixed lambda comparison (first 10 test predictions):\n")
for(i in seq_along(lambda_values)) {
  cat("Lambda =", lambda_values[i], "- Mean prediction:", 
      round(mean(lambda_results[[i]]$fitted.values), 3), "\n")
}

cat("\nBBART with Box-Cox demonstration complete!\n")
cat("Plots saved: bbart_lambda_posterior.png, bbart_lambda_trace.png,\n")
cat("             bbart_prediction_comparison.png, bbart_transformation.png,\n")
cat("             bbart_prediction_accuracy.png\n")