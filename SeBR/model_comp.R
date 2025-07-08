# Comparison script for semiparametric regression models
library(dbarts)
source("helper_funs.R")
source("source_sba.R") 
source("sbart.R")
source("bbart_bc.R")
source("slice.R")

set.seed(123)

# Simulate data from transformed linear model
# dat = simulate_tlm(n = 200, p = 10, g_type = 'step')
# y = dat$y
# X = dat$X 
# y_test = dat$y_test
# X_test = dat$X_test
# best rmse: sblm, best coverage: sbart

n = 200; n_test = 1000; p = 10
X = matrix(runif(n * p), n, p)
X_test = matrix(runif(n_test * p), n_test, p)

# Square root transformation
# f_true = 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
# f_test = 10 * sin(pi * X_test[,1] * X_test[,2]) + 20 * (X_test[,3] - 0.5)^2 + 10 * X_test[,4] + 5 * X_test[,5]
# y = (f_true + rnorm(n, sd = 0.5))^2  # Square transformation
# y_test = (f_test + rnorm(n_test, sd = 0.5))^2
# best rmse: bbart_bc, best coverage: bbart_bc
# seem to recover square root function well

# Log transformation
# f_true = 2 + 0.1 * (10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5])
# f_test = 2 + 0.1 * (10 * sin(pi * X_test[,1] * X_test[,2]) + 20 * (X_test[,3] - 0.5)^2 + 10 * X_test[,4] + 5 * X_test[,5])
# y = exp(f_true + rnorm(n, sd = 0.1))
# y_test = exp(f_test + rnorm(n_test, sd = 0.1))
# best rmse: bbart_bc, best coverage: bbart_bc

# No transformation, just regression comparison
# f_true = 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
# f_test = 10 * sin(pi * X_test[,1] * X_test[,2]) + 20 * (X_test[,3] - 0.5)^2 + 10 * X_test[,4] + 5 * X_test[,5]
# y = f_true + rnorm(n, sd = 1)  # Identity transformation
# y_test = f_test + rnorm(n_test, sd = 1)
# best rmse: bart, best coverage: sbart

# Complex nonlinear function + step transformation
# f_true = 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
# f_test = 10 * sin(pi * X_test[,1] * X_test[,2]) + 20 * (X_test[,3] - 0.5)^2 + 10 * X_test[,4] + 5 * X_test[,5]
# z = f_true + rnorm(n, sd = 0.5)
# z_test = f_test + rnorm(n_test, sd = 0.5)
# y = ifelse(z < 10, z^2, ifelse(z < 20, 2*z + 50, exp((z-20)*0.1) + 90))
# y_test = ifelse(z_test < 10, z_test^2, ifelse(z_test < 20, 2*z_test + 50, exp((z_test-20)*0.1) + 90))
# best rmse: sbart, best coverage: sbart, fastest: bart

# Sigmoid transformation (bounded, non-Box-Cox)
f_true = 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
f_test = 10 * sin(pi * X_test[,1] * X_test[,2]) + 20 * (X_test[,3] - 0.5)^2 + 10 * X_test[,4] + 5 * X_test[,5]
z = f_true + rnorm(n, sd = 0.5)
z_test = f_test + rnorm(n_test, sd = 0.5)
y = 20 / (1 + exp(-z/2))  # Sigmoid, range [0, 20]
y_test = 20 / (1 + exp(-z_test/2))
# best rmse: sbart, best coverage: bbart_bc, fastest: bart


# Storage for results
results = list()
timing = list()

cat("Fitting models...\n")

# 1. Regular BART (no transformation)
cat("Fitting regular BART...\n")
timing$bart = system.time({
  fit_bart = bart(x.train = X, y.train = y, x.test = X_test, 
                  ntree = 200, ndpost = 1000, nskip = 1000)
})
results$bart = list(
  fitted.values = fit_bart$yhat.test.mean,
  post_ypred = fit_bart$yhat.test,
  model = 'bart'
)

# 2. Semiparametric BART (sbart)
cat("Fitting sbart...\n")
timing$sbart = system.time({
  fit_sbart = sbart(y = y, X = X, X_test = X_test, 
                    ntree = 200, nsave = 1000, nburn = 1000)
})
results$sbart = fit_sbart

# 3. Semiparametric Bayesian Linear Model (sblm)
cat("Fitting sblm...\n")
timing$sblm = system.time({
  fit_sblm = sblm(y = y, X = X, X_test = X_test)
})
results$sblm = fit_sblm

# 4. Bayesian BART with Box-Cox (bbart_bc)
cat("Fitting bbart_bc...\n")
timing$bbart_bc = system.time({
  fit_bbart_bc = bbart_bc(y = y, X = X, X_test = X_test,
                          ntree = 200, nsave = 1000, nburn = 1000)
})
results$bbart_bc = fit_bbart_bc

# Performance evaluation
cat("\nEvaluating performance...\n")

# Calculate RMSE for each model
rmse = function(pred, true) sqrt(mean((pred - true)^2))

performance = data.frame(
  Model = c("BART", "SBART", "SBLM", "BBART_BC"),
  RMSE = c(
    rmse(results$bart$fitted.values, y_test),
    rmse(results$sbart$fitted.values, y_test), 
    rmse(results$sblm$fitted.values, y_test),
    rmse(results$bbart_bc$fitted.values, y_test)
  ),
  Time_sec = c(
    timing$bart[3], timing$sbart[3], 
    timing$sblm[3], timing$bbart_bc[3]
  )
)

# Calculate coverage rates (90% prediction intervals)
calc_coverage = function(post_pred, y_true, alpha = 0.1) {
  lower = apply(post_pred, 2, quantile, alpha/2)
  upper = apply(post_pred, 2, quantile, 1 - alpha/2)
  mean(y_true >= lower & y_true <= upper)
}

performance$Coverage_90 = c(
  calc_coverage(results$bart$post_ypred, y_test),
  calc_coverage(results$sbart$post_ypred, y_test),
  calc_coverage(results$sblm$post_ypred, y_test), 
  calc_coverage(results$bbart_bc$post_ypred, y_test)
)

# Add interval width (measure of uncertainty)
calc_mean_width = function(post_pred, alpha = 0.1) {
  lower = apply(post_pred, 2, quantile, alpha/2)
  upper = apply(post_pred, 2, quantile, 1 - alpha/2)
  mean(upper - lower)
}

performance$Mean_Width = c(
  calc_mean_width(results$bart$post_ypred),
  calc_mean_width(results$sbart$post_ypred),
  calc_mean_width(results$sblm$post_ypred),
  calc_mean_width(results$bbart_bc$post_ypred)
)

print(performance)

# Visualization
par(mfrow = c(2, 2))

# Plot 1: Prediction accuracy
plot(y_test, results$bart$fitted.values, main = "BART", 
     xlab = "True", ylab = "Predicted", pch = 16)
abline(0, 1, col = "red")

plot(y_test, results$sbart$fitted.values, main = "SBART", 
     xlab = "True", ylab = "Predicted", pch = 16)
abline(0, 1, col = "red")

plot(y_test, results$sblm$fitted.values, main = "SBLM", 
     xlab = "True", ylab = "Predicted", pch = 16)
abline(0, 1, col = "red")

plot(y_test, results$bbart_bc$fitted.values, main = "BBART_BC", 
     xlab = "True", ylab = "Predicted", pch = 16)
abline(0, 1, col = "red")

dev.new()

# Plot transformation estimates (for models that estimate transformations)
par(mfrow = c(1, 3))

# Plot SBART transformation
if(!is.null(results$sbart$post_g)) {
  y_unique = sort(unique(y))
  g_mean = colMeans(results$sbart$post_g)
  
  if(length(y_unique) == length(g_mean)) {
    plot(y_unique, g_mean, main = "SBART Transformation", 
         xlab = "y", ylab = "g(y)", type = "l", lwd = 2)
    
    if(exists("dat") && "g_true" %in% names(dat) && "y_grid" %in% names(dat)) {
      points(dat$y_grid, dat$g_true, col = "red", pch = 16, cex = 0.5)
      legend("topleft", c("Estimated", "True"), col = c("black", "red"), 
             lty = c(1, NA), pch = c(NA, 16))
    }
  }
}

# Plot SBLM transformation  
if(!is.null(results$sblm$post_g)) {
  y_unique = sort(unique(y))
  g_mean = colMeans(results$sblm$post_g)
  
  if(length(y_unique) == length(g_mean)) {
    plot(y_unique, g_mean, main = "SBLM Transformation", 
         xlab = "y", ylab = "g(y)", type = "l", lwd = 2)
    
    if(exists("dat") && "g_true" %in% names(dat) && "y_grid" %in% names(dat)) {
      points(dat$y_grid, dat$g_true, col = "red", pch = 16, cex = 0.5)
    }
  }
}

# Plot BBART_BC transformation (Box-Cox)
if(!is.null(results$bbart_bc$post_lambda)) {
  lambda_mean = mean(results$bbart_bc$post_lambda)
  y_seq = seq(min(y), max(y), length.out = 100)
  
  g_bc_vals = g_bc(y_seq, lambda = lambda_mean)
  plot(y_seq, g_bc_vals, main = paste("BBART_BC (λ =", round(lambda_mean, 2), ")"), 
       xlab = "y", ylab = "g(y)", type = "l", lwd = 2)
    
  if(exists("dat") && "g_true" %in% names(dat) && "y_grid" %in% names(dat)) {
    points(dat$y_grid, dat$g_true, col = "red", pch = 16, cex = 0.5)
  }
}

# Summary
cat("\n=== SUMMARY ===\n")
cat("Best RMSE:", performance$Model[which.min(performance$RMSE)], 
    "(", round(min(performance$RMSE), 4), ")\n")
cat("Best Coverage:", performance$Model[which.max(performance$Coverage_90)], 
    "(", round(max(performance$Coverage_90), 3), ")\n")
cat("Fastest:", performance$Model[which.min(performance$Time_sec)], 
    "(", round(min(performance$Time_sec), 1), "sec)\n")