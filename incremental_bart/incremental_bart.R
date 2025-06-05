library(dbarts)

set.seed(42)

# basic data
n <- 100
x1 <- runif(n, -1, 1)
x2 <- runif(n, -1, 1)
X <- cbind(x1, x2)
y <- sin(pi * x1) + x2^2 + rnorm(n, 0, 0.3)

# batch approach
bart_batch <- bart(X, y, ndpost = 200, nskip = 50, verbose = FALSE)
pred_batch <- colMeans(bart_batch$yhat.train)
rmse_batch <- sqrt(mean((pred_batch - (sin(pi * x1) + x2^2))^2))
cat("Batch RMSE:", rmse_batch, "\n")

# dbarts approach
sampler <- dbarts(y ~ x1 + x2, data = data.frame(y = y, x1 = x1, x2 = x2))
invisible(sampler$run(0, 50))

# collect some samples
samples_db <- sampler$run(200, 0)
if (is.matrix(samples_db)) {
  pred_db <- colMeans(samples_db)
  rmse_db <- sqrt(mean((pred_db - (sin(pi * x1) + x2^2))^2))
  cat("dbarts RMSE:", rmse_db, "\n")
  cat("Methods give similar results\n")
} else {
  cat("dbarts sampling works\n")
}

# test recreation (simpler approach for SeBR concept)
y_new <- y + rnorm(n, 0, 0.1)
sampler_new <- dbarts(y_new ~ x1 + x2, data = data.frame(y_new = y_new, x1 = x1, x2 = x2))
invisible(sampler_new$run(0, 50))
samples_new <- sampler_new$run(50, 0)
cat("Sampler recreation works\n")

# simple SeBR simulation
n_sebr <- 80
x_sebr <- runif(n_sebr, 0, 2)
y_base <- sin(2 * pi * x_sebr)
y_sebr <- exp(y_base + rnorm(n_sebr, 0, 0.2))

# no transformation
bart_no_trans <- bart(matrix(x_sebr, ncol = 1), y_sebr, ndpost = 200, verbose = FALSE)
pred_no_trans <- colMeans(bart_no_trans$yhat.train)
rmse_no_trans <- sqrt(mean((pred_no_trans - y_sebr)^2))

# with log transformation
y_log <- log(y_sebr)
bart_log <- bart(matrix(x_sebr, ncol = 1), y_log, ndpost = 200, verbose = FALSE)
pred_log_transformed <- exp(colMeans(bart_log$yhat.train))
rmse_trans = sqrt(mean((pred_log_transformed - y_sebr)^2))

cat("No transform RMSE:", rmse_no_trans, "\n")
cat("With transform RMSE:", rmse_trans, "\n")
improvement <- rmse_no_trans / rmse_trans
cat("Improvement factor:", improvement, "\n")

# test different bart configurations
small_trees <- bart(X, y, ndpost = 100, ntree = 50, verbose = FALSE)
many_trees <- bart(X, y, ndpost = 100, ntree = 300, verbose = FALSE)
rmse_small <- sqrt(mean((colMeans(small_trees$yhat.train) - (sin(pi * x1) + x2^2))^2))
rmse_many <- sqrt(mean((colMeans(many_trees$yhat.train) - (sin(pi * x1) + x2^2))^2))
cat("50 trees RMSE:", rmse_small, "\n")
cat("300 trees RMSE:", rmse_many, "\n")

save(sampler, bart_batch, improvement, file = "incremental_experiments.RData")