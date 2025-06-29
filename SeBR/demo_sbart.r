# End-to-end demonstration of semiparametric Bayesian BART (sbart)

# Load required libraries
library(dbarts)

# Source helper functions (adjust paths as needed)
source("helper_funs.R")
source("sbart.R")
source("source_sba.R")  # for g_fun, g_inv_approx, bb, etc.

# Set seed for reproducibility
set.seed(123)

# Simulate data from a transformed nonlinear model
n <- 200
p <- 5
X <- matrix(rnorm(n*p), n, p)

# Create nonlinear function
f_true <- sin(2*X[,1]) + X[,2]^2 - X[,3]*X[,4] + rnorm(n, sd=0.1)

# Apply transformation (log-normal with nonlinear mean)
y <- exp(f_true + 0.5*rnorm(n))

# Create test data
n_test <- 50
X_test <- matrix(rnorm(n_test*p), n_test, p)
f_test_true <- sin(2*X_test[,1]) + X_test[,2]^2 - X_test[,3]*X_test[,4]
y_test_true <- exp(f_test_true + 0.5*rnorm(n_test))

cat("Data simulation complete.\n")
cat("Training data: n =", n, ", p =", p, "\n")
cat("Test data: n_test =", n_test, "\n")
cat("Response range: [", round(min(y), 2), ",", round(max(y), 2), "]\n\n")

# Fit the semiparametric Bayesian BART model
cat("Fitting sbart model...\n")
fit <- sbart(y = y, X = X, X_test = X_test,
             fixedX = FALSE,
             approx_g = FALSE,
             pilot_ndraws = 200,
             pilot_nburn = 50,
             ntree = 100,
             nsave = 500,
             nburn = 200,
             verbose = TRUE)

cat("\nModel fitting complete!\n")
print(names(fit))

# Summary statistics
cat("\nModel summary:\n")
cat("Posterior predictive mean range: [", 
    round(min(fit$fitted.values), 2), ",", 
    round(max(fit$fitted.values), 2), "]\n")
cat("Posterior sigma mean:", round(mean(fit$post_sigma), 3), "\n")
cat("Posterior sigma range: [", 
    round(min(fit$post_sigma), 3), ",", 
    round(max(fit$post_sigma), 3), "]\n")

# Plot 1: Prediction intervals on test data
cat("\nCreating prediction plots...\n")
pi_y <- t(apply(fit$post_ypred, 2, quantile, c(0.05, 0.95)))

png("sbart_predictions.png", width=800, height=600)
plot(1:nrow(pi_y), pi_y[,1], type='n', ylim=range(pi_y),
     xlab='Test observation', ylab='y', 
     main='SBART: 90% Prediction Intervals on Test Data')
segments(1:nrow(pi_y), pi_y[,1], 1:nrow(pi_y), pi_y[,2], col='lightblue', lwd=2)
points(1:nrow(pi_y), fit$fitted.values, pch=16, col='darkblue')
legend('topright', c('90% PI', 'Posterior Mean'), 
       col=c('lightblue', 'darkblue'), lwd=c(2,1), pch=c(NA,16))
dev.off()

# Plot 2: Transformation evolution
y0 <- sort(unique(y))
png("sbart_transformation.png", width=800, height=600)
plot(y0, fit$post_g[1,], type='n', ylim=range(fit$post_g),
     xlab='y', ylab='g(y)', main="SBART: Posterior Draws of Transformation")
for(s in sample(nrow(fit$post_g), 50)) {
  lines(y0, fit$post_g[s,], col='gray', alpha=0.5)
}
lines(y0, colMeans(fit$post_g), lwd=3, col='red')
legend('bottomright', c('Posterior draws', 'Posterior mean'), 
       col=c('gray', 'red'), lwd=c(1,3))
dev.off()

# Plot 3: Posterior predictive checks
png("sbart_pp_check.png", width=800, height=600)
y_grid <- seq(min(y), max(y), length.out=100)
plot(y_grid, y_grid, type='n', ylim=c(0,1),
     xlab='y', ylab='F(y)', main='SBART: Posterior Predictive ECDF Check')

# Plot several posterior predictive ECDFs
for(s in sample(nrow(fit$post_ypred), 20)) {
  lines(y_grid, ecdf(fit$post_ypred[s,])(y_grid), 
        col='lightgray', type='s')
}

# Plot observed ECDF
lines(y_grid, ecdf(y)(y_grid), col='black', type='s', lwd=3)
legend('bottomright', c('Posterior predictive', 'Observed'), 
       col=c('lightgray', 'black'), lwd=c(1,3))
dev.off()

# Plot 4: Sigma trace plot
png("sbart_sigma_trace.png", width=800, height=600)
plot(fit$post_sigma, type='l', xlab='MCMC Iteration', ylab='Sigma',
     main='SBART: Posterior Trace of Error Standard Deviation')
abline(h=mean(fit$post_sigma), col='red', lwd=2)
dev.off()

# Calculate prediction performance (if we had true test values)
# Coverage probability for 90% intervals
# coverage <- mean((pi_y[,1] <= y_test_true) & (pi_y[,2] >= y_test_true))
# cat("90% Prediction interval coverage:", round(coverage, 3), "\n")

# Interval width summary
interval_widths <- pi_y[,2] - pi_y[,1]
cat("Prediction interval width summary:\n")
print(summary(interval_widths))

cat("\nSBART demonstration complete!\n")
cat("Plots saved: sbart_predictions.png, sbart_transformation.png,\n")
cat("             sbart_pp_check.png, sbart_sigma_trace.png\n")